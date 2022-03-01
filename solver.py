import numpy as np
import sympy as sym

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
        self.K_XX, self.K_YY, self.K_ZZ, self.K_XZ, self.m = misc_matrix[2:]

        # 5. eigenvalues to be calculated.
        self.eigenvalues_symmetric = None
        self.eigenvalues_asymmetric = None

    def solve_symmetric_eom(self):
        eigen = sym.Symbol['eigen']

        u_row = [self.C_X_u - (2 * self.m * eigen),
                 self.C_X_alpha,
                 self.C_X_0,
                 0]
        alpha_row = [self.C_Z_u,
                     self.C_Z_alpha + ((self.C_Z_alpha_dot - (2 * self.m)) * eigen),
                     - self.C_X_0,
                     self.C_Z_q + (2 * self.m)]
        theta_row = [0,
                     0,
                     -eigen,
                     1]
        q_row = [self.C_M_u,
                 self.C_M_alpha + (self.C_M_alpha_dot * eigen),
                 0,
                 self.C_M_q - (2 * self.m * self.K_YY * eigen)]

        symmetric_eom_np = np.array([u_row, alpha_row, theta_row, q_row])
        symmetric_eom_sym = sym.Matrix(symmetric_eom_np)
        symmetric_det = symmetric_eom_sym.det()

        self.eigenvalues_symmetric = sym.solve(symmetric_det)

    def solve_asymmetric_eom(self):
        pass


k = sym.Symbol('k')

u = np.array([2, -1, k])
v = np.array([1, 0, 2])
w = np.array([k, 3, k])

z = [u, v, w] #axis with variables

#A = Matrix([u , v, w])
A = sym.Matrix(z)

equation = A.det()
print(sym.solve(equation, k)[0])
