import unittest
import numpy as np

from solver import Solver


class TestSolver(unittest.TestCase):

    def solve_symmetric_eom(self):
        # eigen = sym.Symbol['eigen']
        #
        # u_row = [self.C_X_u - (2 * self.m_c * eigen),
        #          self.C_X_alpha,
        #          self.C_Z_0,
        #          0]
        # alpha_row = [self.C_Z_u,
        #              self.C_Z_alpha + ((self.C_Z_alpha_dot - (2 * self.m_c)) * eigen),
        #              - self.C_X_0,
        #              self.C_Z_q + (2 * self.m_c)]
        # theta_row = [0,
        #              0,
        #              -eigen,
        #              1]
        # q_row = [self.C_M_u,
        #          self.C_M_alpha + (self.C_M_alpha_dot * eigen),
        #          0,
        #          self.C_M_q - (2 * self.m_c * self.K_YY * eigen)]
        #
        # symmetric_eom_np = np.array([u_row, alpha_row, theta_row, q_row])
        # symmetric_eom_sym = sym.Matrix(symmetric_eom_np)
        # symmetric_det = symmetric_eom_sym.det()
        #
        # self.eigenvalues_symmetric = sym.solve(symmetric_det)
        pass

    def solve_asymmetric_eom(self):
        pass


if __name__ == '__main__':
    unittest.main()
