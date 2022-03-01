import numpy as np

"""
Assumptions:

- provided data contains no errors (e.g. negative span) and is in the correct order (e.g. sweep before span);
- asymmetric surface coordinates run from - span to 0 instead of - span / 2 to span / 2;
- if a vertical surface is used the inputs are pre-transformed by 90 degrees along the roll axis;
- induced drag constant along span (elliptical distribution), C_L constant along wing span;
- Aerodynamic center function of (local) X and Y, but not of Z (dihedral height effect neglected);
- Z-component of aerodynamic center due to dihedral neglected;
- sin(sweep) * sin(dihedral) = 0;
- only the turn rate aerodynamic effects are taken into account;
"""

"""
Limitations:

- Scipy could not be used for integration as most methods return tuples, which Scipy cannot handle;
"""


# TODO: implement control surfaces
# TODO: LE translated in both X, Y and Z


class AerodynamicSurface:

    def __init__(self, path, symmetric=True, vertical=False, main=None):

        # 1. load data provided by user.
        with open(f'{path}/curves.txt') as file:
            self.data_aerodynamic = file.readlines()

        with open(f'{path}/values.txt') as file:
            self.data_geometry = np.genfromtxt(file, dtype=str)

        # 2. attributes loaded from the values file.
        self.chord_root = float(self.data_geometry[0, 1])
        self.span = float(self.data_geometry[1, 1])
        self.sweep = float(self.data_geometry[2, 1]) * (np.pi / 180)
        self.taper = float(self.data_geometry[3, 1])
        self.dihedral = float(self.data_geometry[4, 1]) * (np.pi / 180)
        self.ac = float(self.data_geometry[11, 1])

        self.X_le = float(self.data_geometry[5, 1])
        self.Y_le = float(self.data_geometry[6, 1])
        self.Z_le = float(self.data_geometry[7, 1])

        self.C_M_X_ac = float(self.data_geometry[8, 1])
        self.C_M_Y_ac = float(self.data_geometry[9, 1])
        self.C_M_Z_ac = float(self.data_geometry[10, 1])

        # 3. attributes from input.
        self.symmetric = symmetric
        self.vertical = vertical
        self.main = main

        # 4. attributes calculated, loaded or modified by the initialization function.
        self.position = np.array([[self.X_le], [self.Y_le], [self.Z_le]])

        self.alpha = None
        self.downwash = None
        self.C_L = None
        self.C_D = None

        self.area = None
        self.mac = None

        self.X_ac = None
        self.Y_ac = None
        self.Z_ac = None

        self.initialize()

        # 5. attributes in the body frame calculated by the transformation function.
        self.angles = None
        self.rates = None
        self.velocity = None
        self.acceleration = None

        # 6. attributes related to dynamic pressure.
        self.magnitude = None
        self.rho = None

        # 7. attributes calculated and transformed by the forces and moments functions.
        self.output_body = np.zeros((3, 3))
        self.output_earth = None

        # Code below only used for verification tests.
        self.calculate_local_velocity = None
        self.calculate_local_forces = None

    def initialize(self):

        def curves():
            data, start = [], False

            # 1. detect header in curves file and append labels.
            for i in range(len(self.data_aerodynamic)):
                if 'alpha' in self.data_aerodynamic[i]:
                    start = True
                    data.append(self.data_aerodynamic[i].strip().split())

                # 2. detect complete rows and append data.
                elif start:
                    row = [float(j) for j in self.data_aerodynamic[i].strip().split()]

                    if len(row) == len(data[0]):
                        data.append(row)

            # 3. exclude header.
            data = np.array(data[1:])

            # 4. extract alpha, C_L, C_D lists and assign to the class attributes.
            self.alpha = data[:, 0] * (np.pi / 180)
            self.C_L = data[:, 2]
            self.C_D = data[:, 5]

        def area():
            self.area = (self.span / 2) * self.chord_root * (1 + self.taper)

        def position():
            if self.vertical:
                self.position = self.transformation(self.position, transform_vertical=True)

        def mean_aerodynamic_chord():
            if self.symmetric:
                self.Y_ac = - (self.span * (1 + (2 * self.taper))) / (6 * (1 + self.taper))
            else:
                self.Y_ac = - (self.span * (1 + (2 * self.taper))) / (3 * (1 + self.taper))

            self.mac = self.chord(self.Y_ac)

        def mean_aerodynamic_center():
            self.X_ac, self.Y_ac, self.Z_ac = self.aerodynamic_center(self.Y_ac)

        def downwash():
            self.main.C_L_alpha = np.mean(np.gradient(self.main.C_L))
            self.main.ar = self.main.span / self.main.mac

            nominator = (7 * self.main.C_L_alpha)
            denominator = (4 * np.pi * self.main.ar * (1 + abs(self.Z_ac - self.main.Z_ac)) *
                           (self.main.taper * abs(self.X_ac - self.main.X_ac)) ** (1 / 4))

            self.downwash = nominator / denominator

        curves()
        area()
        position()

        mean_aerodynamic_chord()
        mean_aerodynamic_center()

        if self.main is not None:
            downwash()

    def chord(self, coordinate):
        if self.symmetric:
            return self.chord_root - self.chord_root * (1 - self.taper) * ((2 * abs(coordinate)) / self.span)
        else:
            return self.chord_root + self.chord_root * (1 - self.taper) * (coordinate / self.span)

    def aerodynamic_center(self, coordinate):
        X_initial_position = self.position[0, 0] - (self.chord_root / 2)
        Y_initial_position = self.position[1, 0] + coordinate
        Z_initial_position = self.position[2, 0]

        X_shift_due_ac = (1 - (2 * self.ac)) * (self.chord(coordinate) / 2)
        X_shift_due_sweep = - abs(coordinate) * np.sin(self.sweep)
        Z_shift_due_dihedral = - abs(coordinate) * self.dihedral

        X_final_position = X_initial_position + X_shift_due_ac + X_shift_due_sweep
        Y_final_position = Y_initial_position
        Z_final_position = Z_initial_position + Z_shift_due_dihedral

        return X_final_position, Y_final_position, Z_final_position

    def transformation(self, inputs, transform_vertical=False, earth_to_body=True):

        def transformation_vertical(rotation=1):
            return np.array([[1, 0, 0],
                             [0, 0, rotation],
                             [0, - rotation, 0]])

        def transformation_roll(roll):
            return np.array([[1, 0, 0],
                             [0, np.cos(roll), np.sin(roll)],
                             [0, - np.sin(roll), np.cos(roll)]])

        def transformation_pitch(pitch):
            return np.array([[np.cos(pitch), 0, - np.sin(pitch)],
                             [0, 1, 0],
                             [np.sin(pitch), 0, np.cos(pitch)]])

        def transformation_yaw(yaw):
            return np.array([[np.cos(yaw), np.sin(yaw), 0],
                             [- np.sin(yaw), np.cos(yaw), 0],
                             [0, 0, 1]])

        # 1. check if vertical transformation applies due to surface rotation.
        if transform_vertical:
            if earth_to_body:
                transformation_matrix = transformation_vertical()
            else:
                transformation_matrix = transformation_vertical(rotation=-1)

        # 2. check if transformation from the Earth to the body reference frame applies (order: roll, pitch, yaw).
        elif earth_to_body:
            transformation_matrix = np.dot(transformation_roll(self.angles[0, 0]),
                                           transformation_pitch(self.angles[1, 0]))
            transformation_matrix = np.dot(transformation_matrix, transformation_yaw(self.angles[2, 0]))

        # 3. check if transformation from the body to the Earth reference frame applies (order: yaw, pitch, roll).
        else:
            transformation_matrix = np.dot(transformation_yaw(- self.angles[2, 0]),
                                           transformation_pitch(- self.angles[1, 0]))
            transformation_matrix = np.dot(transformation_matrix, transformation_roll(- self.angles[0, 0]))

        return np.dot(transformation_matrix, inputs)

    def calculate_forces(self, rho):

        def roll_dot_effect(coordinate):
            return coordinate * self.rates[0, 0]

        def pitch_dot_effect(coordinate):
            return - self.aerodynamic_center(coordinate)[0] * self.rates[1, 0]

        def yaw_dot_effect(coordinate):
            return - coordinate * self.rates[2, 0]

        def calculate_local_velocity(coordinate):
            # 1. local effects due to turn rates.
            X_dot_due_yaw_dot = yaw_dot_effect(coordinate)
            Z_dot_due_roll_dot = roll_dot_effect(coordinate)
            Z_dot_due_pitch_dot = pitch_dot_effect(coordinate)

            # 2. local effects due to dihedral.
            Y_dot_due_dihedral = self.velocity[2, 0] * self.dihedral * - (coordinate / abs(coordinate))
            Z_dot_due_dihedral = self.velocity[1, 0] * self.dihedral * (coordinate / abs(coordinate))

            X_dot = self.velocity[0, 0] + X_dot_due_yaw_dot
            Y_dot = self.velocity[1, 0] + Y_dot_due_dihedral
            Z_dot = self.velocity[2, 0] + Z_dot_due_roll_dot + Z_dot_due_pitch_dot + Z_dot_due_dihedral

            # 3. local velocity and local angles.
            local_velocity = np.array([[X_dot], [Y_dot], [Z_dot]])
            local_magnitude = np.linalg.norm(local_velocity)

            local_alpha = np.arctan(local_velocity[2, 0] / local_velocity[0, 0])
            local_alpha_dot = (np.arctan((local_velocity[2, 0] + self.acceleration[2, 0]) /
                                         (local_velocity[0, 0] + self.acceleration[0, 0])) - local_alpha)

            local_beta = np.arctan(local_velocity[1, 0] / local_velocity[0, 0])
            local_beta_dot = (np.arctan((local_velocity[1, 0] + self.acceleration[1, 0]) /
                                        (local_velocity[0, 0] + self.acceleration[0, 0])) - local_beta)

            return local_velocity, local_magnitude, local_alpha, local_alpha_dot, local_beta, local_beta_dot

        def calculate_local_forces(coordinate):
            # 1. local velocities.
            local_velocity, local_magnitude, local_alpha, local_alpha_dot, local_beta, local_beta_dot = \
                calculate_local_velocity(coordinate)

            # 2. local geometry.
            local_chord = self.chord(coordinate)

            # 3. local downwash.
            if self.main is not None:
                local_delay = local_alpha_dot * (abs(self.X_ac - self.main.X_ac) / local_velocity[0, 0])
                local_downwash = (local_alpha - local_delay) * self.downwash
            else:
                local_downwash = 0

            # 4. local lift and drag coefficients.
            local_C_L = np.interp(local_alpha - local_downwash, self.alpha, self.C_L)
            local_C_D = np.interp(local_alpha - local_downwash, self.alpha, self.C_D)

            # 5. local lift and drag forces.
            local_L = (rho / 2) * local_magnitude ** 2 * local_C_L * local_chord * \
                      (1 + ((coordinate / abs(coordinate)) * local_beta * np.sin(2 * self.sweep)))
            local_D = (rho / 2) * local_magnitude ** 2 * local_C_D * local_chord * \
                      (1 + ((coordinate / abs(coordinate)) * local_beta * np.sin(2 * self.sweep)))

            # 6. local normal and tangential forces.
            local_F_T = - (np.sin(local_alpha) * local_L + np.cos(local_alpha) * local_D)
            local_F_N = - (np.cos(local_alpha) * local_L + np.sin(local_alpha) * local_D)

            # 7. local forces in the body reference frame.
            F_X = local_F_T
            F_Y = local_F_N * self.dihedral * - (coordinate / abs(coordinate))
            F_Z = local_F_N

            return F_X, F_Y, F_Z

        step_size = 0.1

        if self.symmetric:
            steps = int((self.span - self.span % step_size) / (2 * step_size)) + 1

            for i in range(-steps, 0):
                self.output_body[:, 0] += calculate_local_forces((i * step_size) + (step_size / 2))

            for i in range(0, steps):
                self.output_body[:, 1] += calculate_local_forces((i * step_size) + (step_size / 2))

        else:
            steps = int((self.span - self.span % step_size) / step_size) + 1

            for i in range(-steps, 0):
                self.output_body[:, 0] += calculate_local_forces((i * step_size) + (step_size / 2))

        # Code below only used for verification tests.
        self.calculate_local_velocity = calculate_local_velocity(coordinate=self.span / 2)
        self.calculate_local_forces = calculate_local_forces(coordinate=self.span / 2)

    def calculate_moments(self):
        M_Y_due_F_X = self.Z_ac * (self.output_body[0, 0] + self.output_body[0, 1])
        M_Z_due_F_X = self.Y_ac * (self.output_body[0, 1] - self.output_body[0, 0])

        M_X_due_F_Y = self.Z_ac * (self.output_body[1, 0] + self.output_body[1, 1]) * - 1
        M_Z_due_F_Y = self.X_ac * (self.output_body[1, 0] + self.output_body[1, 1])

        M_X_due_F_Z = self.Y_ac * (self.output_body[2, 0] - self.output_body[2, 1])
        M_Y_due_F_Z = self.X_ac * (self.output_body[2, 0] + self.output_body[2, 1]) * - 1

        M_Y_due_ac = (self.C_M_Y_ac / 2) * self.rho * self.magnitude ** 2 * self.area * self.mac

        self.output_body[0, 2] = M_X_due_F_Y + M_X_due_F_Z
        self.output_body[1, 2] = M_Y_due_F_X + M_Y_due_F_Z + M_Y_due_ac
        self.output_body[2, 2] = M_Z_due_F_X + M_Z_due_F_Y

    def calculate_outputs(self, inputs, magnitude, rho):
        self.magnitude = magnitude
        self.rho = rho

        if self.vertical:
            inputs = self.transformation(inputs, transform_vertical=True)

        self.angles = inputs[:, 0].reshape(-1, 1)
        self.rates = inputs[:, 1].reshape(-1, 1)
        self.velocity = self.transformation(inputs[:, 2].reshape(-1, 1))
        self.acceleration = self.transformation(inputs[:, 3].reshape(-1, 1))

        self.calculate_forces(rho)
        self.calculate_moments()

        self.output_earth = self.transformation(self.output_body, earth_to_body=False)

        if self.vertical:
            self.output_earth = self.transformation(self.output_body, transform_vertical=True, earth_to_body=False)

        return self.output_body, self.output_earth
