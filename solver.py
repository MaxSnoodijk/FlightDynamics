import numpy as np

from coefficients import Coefficients


class Solver:

    def __init__(self):
        symmetric_matrix, asymmetric_matrix, misc_matrix = Coefficients()

        # 1. symmetric attributes.
        self.C_X_0, self.C_Z_0, self.C_M_0 = symmetric_matrix[0, :]
        self.C_X_u, self.C_Z_u, self.C_M_u = symmetric_matrix[1, :]
        self.C_X_alpha, self.C_Z_alpha, self.C_M_alpha = symmetric_matrix[2, :]
        self.C_X_alpha_dot, self.C_Z_alpha_dot, self.C_M_alpha_dot = symmetric_matrix[3, :]
        self.C_X_q, self.C_Z_q, self.C_M_q = symmetric_matrix[4, :]
        self.C_X_elevator, self.C_Z_elevator, self.C_M_elevator = symmetric_matrix[5, :]

        # 2. asymmetric attributes.
        self.C_Y_beta, self.C_l_beta, self.C_n_beta = asymmetric_matrix[0, :]
        self.C_Y_p, self.C_l_p, self.C_n_p = asymmetric_matrix[1, :]
        self.C_Y_r, self.C_l_r, self.C_n_r = asymmetric_matrix[2, :]
        self.C_Y_ailerons, self.C_l_ailerons, self.C_n_ailerons = asymmetric_matrix[3, :]
        self.C_Y_rudder, self.C_l_rudder, self.C_n_rudder = asymmetric_matrix[4, :]

        # 3. aerodynamic attributes.
        self.velocity, self.density = misc_matrix[:2]

        # 4. aircraft attributes.
        self.K_XX, self.K_YY, self.K_ZZ, self.K_XZ = misc_matrix[2:]

    def solve_symmetric_eom(self):
        pass

    def solve_asymmetric_eom(self):
        pass
