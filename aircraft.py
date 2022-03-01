import numpy as np

from surface import AerodynamicSurface


class Aircraft:

    def __init__(self, horizontal=True, vertical=True):
        self.wing = AerodynamicSurface('Data/Wing')
        self.aircraft = [self.wing]

        if horizontal:
            self.horizontal = AerodynamicSurface('Data/Horizontal')
            self.aircraft.append(self.horizontal)

        if vertical:
            self.vertical = AerodynamicSurface('Data/Vertical', symmetric=False)
            self.aircraft.append(self.vertical)

        self.mass = 1
        self.g = 9.81
        self.weight = self.mass * self.g

        self.I_XX = None
        self.I_YY = None
        self.I_ZZ = None

        self.X_dot = 10
        self.Y_dot = 0
        self.Z_dot = -0.5

        self.X_dot_dot = 0
        self.Y_dot_dot = 0
        self.Z_dot_dot = 0

        self.roll = 0
        self.pitch = 0.05
        self.yaw = 0

        self.roll_dot = 0
        self.pitch_dot = 0
        self.yaw_dot = 0

        self.ailerons = 0
        self.elevator = 0
        self.rudder = 0

        self.angles = None
        self.velocities = None
        self.rho = None

        self.magnitude_velocity = np.sqrt(self.X_dot ** 2 + self.Y_dot ** 2 + self.Z_dot ** 2)
        self.magnitude_reference = self.calculate_magnitudes()

        self.derivatives_symmetric = np.zeros((3, 6))
        self.derivatives_asymmetric = np.zeros((3, 5))

        self.calculate_symmetric_derivatives()
        self.calculate_asymmetric_derivatives()

    def setup(self, rho=1.225):
        self.angles = np.array([[self.roll, self.roll_dot],
                                [self.pitch, self.pitch_dot],
                                [self.yaw, self.yaw_dot]])

        self.velocities = np.array([[self.X_dot, self.X_dot_dot],
                                    [self.Y_dot, self.Y_dot_dot],
                                    [self.Z_dot, self.Z_dot_dot]])

        self.rho = rho
        self.magnitude_velocity = np.sqrt(self.X_dot ** 2 + self.Y_dot ** 2 + self.Z_dot ** 2)

    def calculate_magnitudes(self):
        self.setup()

        magnitude_array = np.zeros((3, 2))

        for surface in self.aircraft:
            magnitude = surface.calculate_outputs(self.angles, self.velocities, self.rho)[0]

            magnitude_array[0, 0] += magnitude[0, 0] + magnitude[0, 1]  # C_X
            magnitude_array[1, 0] += magnitude[1, 0] + magnitude[1, 1]  # C_Y
            magnitude_array[2, 0] += magnitude[2, 0] + magnitude[2, 1]  # C_Z

            magnitude_array[0, 1] += magnitude[0, 2]  # C_M_X (C_L in reader)
            magnitude_array[1, 1] += magnitude[1, 2]  # C_M_Y (C_M in reader)
            magnitude_array[2, 1] += magnitude[2, 2]  # C_M_Z (C_N in reader)

        return magnitude_array

    def calculate_symmetric_derivatives(self):

        def calculate_derivative(entry):
            magnitude_new = self.calculate_magnitudes()
            magnitude_difference = magnitude_new - self.magnitude_reference

            self.derivatives_symmetric[entry, 0] = (2 * magnitude_difference[0, 0]) / (
                    self.rho * (self.magnitude_velocity ** 2) * self.wing.area)
            self.derivatives_symmetric[entry, 1] = (2 * magnitude_difference[2, 0]) / (
                    self.rho * (self.magnitude_velocity ** 2) * self.wing.area)
            self.derivatives_symmetric[entry, 2] = (2 * magnitude_difference[1, 1]) / (
                    self.rho * (self.magnitude_velocity ** 2) * self.wing.area * self.wing.mac)

        def derivative_initial():
            self.derivatives_symmetric[0, 0] = 0
            self.derivatives_symmetric[0, 1] = - (2 * self.weight) / (
                    self.rho * (self.magnitude_velocity ** 2) * self.wing.area)
            self.derivatives_symmetric[0, 2] = 0

        def derivative_x_dot(change=1, entry=1):
            self.X_dot += change
            calculate_derivative(entry)
            self.X_dot -= change

        def derivative_alpha(change=1, entry=2):
            self.Z_dot += change
            calculate_derivative(entry)
            self.Z_dot -= change

        def derivative_alpha_dot(change=1, entry=3):
            self.Z_dot_dot += change
            calculate_derivative(entry)
            self.X_dot -= change

        def derivative_pitch_dot(change=1, entry=4):
            self.pitch_dot += change
            calculate_derivative(entry)
            self.X_dot -= change

        # def derivative_elevator(change=1, entry=5):
        #     self.elevator += change
        #     calculate_derivative(entry)
        #     self.elevator -= change

        derivative_initial()
        derivative_x_dot()
        derivative_alpha()
        derivative_alpha_dot()
        derivative_pitch_dot()
        # derivative_elevator()

    def calculate_asymmetric_derivatives(self):
        pass
