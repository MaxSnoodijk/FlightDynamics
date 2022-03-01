import matplotlib.pyplot as plt
import numpy as np
import control as ctrl

from coefficients import coefficients


class Solver:

    def __init__(self):
        # Obtaining coefficients.
        (self.symmetric_coefficients, self.symmetric_misc_coefficients,
         self.asymmetric_coefficients, self.asymmetric_misc_coefficients) = coefficients(example=True)

        # Equations of motion in state-space format, calculated by state-space functions.
        self.A_symmetric = None
        self.B_symmetric = None
        self.C_symmetric = None
        self.D_symmetric = None

        self.A_asymmetric = None
        self.B_asymmetric = None
        self.C_asymmetric = None
        self.D_asymmetric = None

        # Eigenvalues of the equations of motion, calculated by state-space functions.
        self.eigenvalues_symmetric = None
        self.eigenvalues_asymmetric = None

        # Initialization of the class.
        self.state_space_symmetric_eom()
        self.state_space_asymmetric_eom()

    def state_space_symmetric_eom(self):
        # Forces and moments coefficients.
        C_X_0, C_Z_0, C_M_0 = self.symmetric_coefficients[0, :]
        C_X_u, C_Z_u, C_M_u = self.symmetric_coefficients[1, :]
        C_X_alpha, C_Z_alpha, C_M_alpha = self.symmetric_coefficients[2, :]
        C_X_alpha_dot, C_Z_alpha_dot, C_M_alpha_dot = self.symmetric_coefficients[3, :]
        C_X_q, C_Z_q, C_M_q = self.symmetric_coefficients[4, :]
        C_X_elevator, C_Z_elevator, C_M_elevator = self.symmetric_coefficients[5, :]
        C_X_trim, C_Z_trim, C_M_trim = self.symmetric_coefficients[6, :]

        # Aerodynamic coefficients.
        velocity = self.symmetric_misc_coefficients[0]

        # Geometric coefficients.
        m_c, mac, K_YY = self.symmetric_misc_coefficients[1:]

        # Setting up the matrices to transform to the state-space format.
        P = - np.array([[2 * m_c * (mac / velocity), 0, 0, 0],
                        [0, ((2 * m_c) - C_Z_alpha_dot) * (mac / velocity), 0, 0],
                        [0, 0, mac / velocity, 0],
                        [0, - C_M_alpha_dot * (mac / velocity), 0, 2 * m_c * K_YY * (mac / velocity)]])

        Q = - np.array([[C_X_u, C_X_alpha, C_Z_0, 0],
                        [C_Z_u, C_Z_alpha, - C_X_0, C_Z_q + (2 * m_c)],
                        [0, 0, 0, 1],
                        [C_M_u, C_M_alpha, 0, C_M_q]])

        R = - np.array([[C_X_elevator, C_X_trim],
                        [C_Z_elevator, C_Z_trim],
                        [0, 0],
                        [C_M_elevator, C_M_trim]])

        # Assigning state-space matrices to class attributes.
        self.A_symmetric = np.dot(np.linalg.inv(P), Q)
        self.B_symmetric = np.dot(np.linalg.inv(P), R)

        self.C_symmetric = np.array([[1, 0, 0, 0],
                                     [0, 1, 0, 0],
                                     [0, 0, 1, 0],
                                     [0, 0, 0, 1]])

        self.D_symmetric = np.array([[0, 0],
                                     [0, 0],
                                     [0, 0],
                                     [0, 0]])

        # Calculating and assigning eigenvalues to class attribute.
        self.eigenvalues_symmetric = np.linalg.eigvals(self.A_symmetric)

    def state_space_asymmetric_eom(self):
        # Asymmetric attributes.
        C_Y_beta, C_l_beta, C_n_beta = self.asymmetric_coefficients[0, :]
        C_Y_beta_dot, C_l_beta_dot, C_n_beta_dot = self.asymmetric_coefficients[1, :]
        C_Y_p, C_l_p, C_n_p = self.asymmetric_coefficients[2, :]
        C_Y_r, C_l_r, C_n_r = self.asymmetric_coefficients[3, :]
        C_Y_ailerons, C_l_ailerons, C_n_ailerons = self.asymmetric_coefficients[4, :]
        C_Y_rudder, C_l_rudder, C_n_rudder = self.asymmetric_coefficients[5, :]

        # Aerodynamic coefficients.
        velocity, C_L = self.asymmetric_misc_coefficients[:2]

        # Geometric coefficients.
        m_b, span, K_XX, K_ZZ, K_XZ = self.asymmetric_misc_coefficients[2:]

        # Setting up the matrices to transform to the state-space format.
        P = - np.array([[(2 * m_b - C_Y_beta_dot) * (span / velocity), 0, 0, 0],
                        [0, span / (2 * velocity), 0, 0],
                        [0, 0, 4 * m_b * K_XX * (span / velocity), -4 * m_b * K_XZ * (span / velocity)],
                        [- C_n_beta_dot * (span / velocity), 0,
                         -4 * m_b * K_XZ * (span / velocity), 4 * m_b * K_XX * (span / velocity)]])

        Q = - np.array([[C_Y_beta, C_L, C_Y_p, C_Y_r - (4 * m_b)],
                        [0, 0, 1, 0],
                        [C_l_beta, 0, C_l_p, C_l_r],
                        [C_n_beta, 0, C_n_p, C_n_r]])

        R = - np.array([[C_Y_ailerons, C_Y_rudder],
                        [0, 0],
                        [C_l_ailerons, C_l_rudder],
                        [C_n_ailerons, C_n_rudder]])

        # Assigning state-space matrices to class attributes.
        self.A_asymmetric = np.dot(np.linalg.inv(P), Q)
        self.B_asymmetric = np.dot(np.linalg.inv(P), R)

        self.C_asymmetric = np.array([[1, 0, 0, 0],
                                      [0, 1, 0, 0],
                                      [0, 0, 1, 0],
                                      [0, 0, 0, 1]])

        self.D_asymmetric = np.array([[0, 0],
                                      [0, 0],
                                      [0, 0],
                                      [0, 0]])

        # Calculating and assigning eigenvalues to class attribute.
        self.eigenvalues_asymmetric = np.linalg.eigvals(self.A_asymmetric)

    def control_symmetric_eom(self, input_type, input_deflection=0, start=0, stop=5, step=0.01):
        # Creating the control system.
        system = ctrl.ss(self.A_symmetric, self.B_symmetric, self.C_symmetric, self.D_symmetric)
        time = np.arange(start, stop + step, step)

        # Calculating the responses.
        if input_type == 'impulse':
            time, response = ctrl.impulse_response(system, time)
        elif input_type == 'initial':
            time, response = ctrl.initial_response(system, time)
        elif input_type == 'step':
            time, response = ctrl.step_response(system, time)
        else:
            time, response = ctrl.forced_response(system, time, input_deflection)

        return time, response

    def control_asymmetric_eom(self, input_type, input_deflection=0, start=0, stop=10, step=0.01):
        # Creating the control system.
        system = ctrl.ss(self.A_asymmetric, self.B_asymmetric, self.C_asymmetric, self.D_asymmetric)
        time = np.arange(start, stop + step, step)

        # Calculating the responses.
        if input_type == 'impulse':
            time, response = ctrl.impulse_response(system, time)
        elif input_type == 'initial':
            time, response = ctrl.initial_response(system, time)
        elif input_type == 'step':
            time, response = ctrl.step_response(system, time)
        else:
            time, response = ctrl.forced_response(system, time, input_deflection)

        return time, response

    def plot_eigenvalues(self):
        figure, axis = plt.subplots(1, 2)

        symmetric_real = [i.real for i in self.eigenvalues_symmetric]
        symmetric_imaginary = [i.imag for i in self.eigenvalues_symmetric]

        asymmetric_real = [i.real for i in self.eigenvalues_asymmetric]
        asymmetric_imaginary = [i.imag for i in self.eigenvalues_asymmetric]

        axis[0].scatter(symmetric_real, symmetric_imaginary)
        axis[1].scatter(asymmetric_real, asymmetric_imaginary)

        axis[0].set(xlabel='Real [-]', ylabel='Imaginary [-]', title='Symmetric eigenvalues')
        axis[1].set(xlabel='Real [-]', ylabel='Imaginary [-]', title='Asymmetric eigenvalues')

        plt.show()

    def plot_responses(self, input_type, input_deflection=0, symmetric=True):
        figure, axis = plt.subplots(2, 2)

        if symmetric:
            time, response = self.control_symmetric_eom(input_type, input_deflection=input_deflection)

            axis[0, 0].set(xlabel='Time [s]', ylabel='Velocity [-]')
            axis[0, 1].set(xlabel='Time [s]', ylabel='Alpha [rad]')
            axis[1, 0].set(xlabel='Time [s]', ylabel='Theta [rad]')
            axis[1, 1].set(xlabel='Time [s]', ylabel='Pitch rate [-]')

        else:
            time, response = self.control_asymmetric_eom(input_type, input_deflection=input_deflection)

            axis[0, 0].set(xlabel='Time [s]', ylabel='Velocity [-]')
            axis[0, 1].set(xlabel='Time [s]', ylabel='Alpha [rad]')
            axis[1, 0].set(xlabel='Time [s]', ylabel='Theta [rad]')
            axis[1, 1].set(xlabel='Time [s]', ylabel='Pitch rate [-]')

        axis[0, 0].plot(time, response[0])
        axis[0, 1].plot(time, response[1])
        axis[1, 0].plot(time, response[2])
        axis[1, 1].plot(time, response[3])

        plt.show()
