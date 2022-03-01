import unittest
import numpy as np

from surface import AerodynamicSurface


class TestAerodynamicSurface(unittest.TestCase):

    def load_surface(self, symmetric=True, vertical=False, translated=False):
        inputs = np.array([[0.349, 0.070, 15.215, 0.1],
                           [0.175, 0.035, 13.108, 0.2],
                           [0.698, 0.035, -1.300, 0.3]])

        input_rho = 1.225
        input_magnitude = np.linalg.norm(inputs[:, 2])

        if symmetric:
            if translated:
                surface = AerodynamicSurface('Data/Testing/Translated', symmetric=symmetric, vertical=vertical)
            else:
                surface = AerodynamicSurface('Data/Testing/Symmetric', symmetric=symmetric, vertical=vertical)
        else:
            surface = AerodynamicSurface('Data/Testing/Asymmetric', symmetric=symmetric, vertical=vertical)

        surface.calculate_outputs(inputs, input_magnitude, input_rho)

        return surface

    def test_initialization(self, margin=0.01):

        def test_symmetric():
            surface = self.load_surface()

            # 1. attributes independent of symmetry condition.
            alpha_true = -0.087  # first alpha entry in curves text file.
            alpha_calc = round(surface.alpha[0], 3)

            C_L_true = -0.418  # first C_L entry in curves text file.
            C_L_calc = round(surface.C_L[0], 3)

            C_D_true = 0.019  # first C_D entry in curves text file.
            C_D_calc = round(surface.C_D[0], 3)

            C_M_Y_ac_true = 0.056  # C_M_Y_ac entry in values text file.
            C_M_Y_ac_calc = round(surface.C_M_Y_ac, 3)

            # 2. tests.
            self.assertEqual(alpha_true, alpha_calc)
            self.assertEqual(C_L_true, C_L_calc)
            self.assertEqual(C_D_true, C_D_calc)
            self.assertEqual(C_M_Y_ac_true, C_M_Y_ac_calc)

            # 3. attributes dependent on symmetry condition.
            area_true = 0.180  # calculated by hand using formula [1].
            area_calc = round(surface.area, 3)

            mac_true = 0.156  # calculated by hand using formula [2].
            mac_calc = round(surface.mac, 3)

            X_ac_true = -0.044  # calculated by hand using formula [3.1].
            X_ac_calc = round(surface.X_ac, 3)

            Y_ac_true = -0.267  # calculated by hand using formula [3.2].
            Y_ac_calc = round(surface.Y_ac, 3)

            Z_ac_true = -0.023  # calculated by hand using formula [3.3].
            Z_ac_calc = round(surface.Z_ac, 3)

            # 4. tests.
            self.assertAlmostEqual(area_true, area_calc, delta=margin * area_true)
            self.assertAlmostEqual(mac_true, mac_calc, delta=margin * mac_true)

            self.assertAlmostEqual(X_ac_true, X_ac_calc, delta=margin * X_ac_true)
            self.assertAlmostEqual(Y_ac_true, Y_ac_calc, delta=margin * Y_ac_true)
            self.assertAlmostEqual(Z_ac_true, Z_ac_calc, delta=margin * Z_ac_true)

        def test_asymmetric():
            surface = self.load_surface(symmetric=False)

            # 1. attributes dependent on symmetry condition.
            area_true = 0.090  # calculated by hand using formula [1].
            area_calc = round(surface.area, 3)

            mac_true = 0.156  # calculated by hand using formula [2].
            mac_calc = round(surface.mac, 3)

            X_ac_true = -0.044  # calculated by hand using formula [3.1].
            X_ac_calc = round(surface.X_ac, 3)

            Y_ac_true = -0.267  # calculated by hand using formula [3.2].
            Y_ac_calc = round(surface.Y_ac, 3)

            Z_ac_true = -0.023  # calculated by hand using formula [3.3].
            Z_ac_calc = round(surface.Z_ac, 3)

            # 2. tests.
            self.assertAlmostEqual(area_true, area_calc, delta=margin * area_true)
            self.assertAlmostEqual(mac_true, mac_calc, delta=margin * mac_true)

            self.assertAlmostEqual(X_ac_true, X_ac_calc, delta=margin * X_ac_true)
            self.assertAlmostEqual(Y_ac_true, Y_ac_calc, delta=margin * Y_ac_true)
            self.assertAlmostEqual(Z_ac_true, Z_ac_calc, delta=margin * Z_ac_true)

        test_symmetric()
        test_asymmetric()

    def test_transformation_input(self, margin=0.01):

        def test_horizontal():
            surface = self.load_surface()

            # 1. attributes dependent on translation and rotation of surface.
            position_true = np.array([[0.040], [0], [0]])  # unmodified leading edge entries in values text file.
            position_calc = np.round(surface.position, 3)

            # 2. attributes dependent on rotation of surface.
            angles_true = np.array([[0.349], [0.175], [0.698]])  # unmodified input angles.
            angles_calc = np.round(surface.angles, 3)

            rates_true = np.array([[0.070], [0.035], [0.035]])  # unmodified input rates.
            rates_calc = np.round(surface.rates, 3)

            velocity_true = np.array([[20.001], [1.006], [1.992]])  # calculated by hand using formula [4].
            velocity_calc = np.round(surface.velocity, 3)

            acceleration_true = np.array([[0.150], [0.197], [0.281]])  # calculated by hand using formula [4].
            acceleration_calc = np.round(surface.acceleration, 3)

            # 3. tests.
            self.assertTrue(np.allclose(position_true, position_calc, rtol=margin))
            self.assertTrue(np.allclose(angles_true, angles_calc, rtol=margin))
            self.assertTrue(np.allclose(rates_true, rates_calc, rtol=margin))
            self.assertTrue(np.allclose(velocity_true, velocity_calc, rtol=margin))
            self.assertTrue(np.allclose(acceleration_true, acceleration_calc, rtol=margin))

        def test_vertical():

            def test_untranslated():
                surface = self.load_surface(vertical=True)

                # 1. attributes dependent on translation and rotation of surface.
                position_true = np.array([[0.040], [0], [0]])  # calculated by hand using formula [4.1].
                position_calc = np.round(surface.position, 3)

                # 2. attributes dependent on rotation of surface.
                angles_true = np.array([[0.349], [0.698], [-0.175]])  # calculated by hand using formula [4.1].
                angles_calc = np.round(surface.angles, 3)

                rates_true = np.array([[0.070], [0.035], [-0.035]])  # calculated by hand using formula [4.1].
                rates_calc = np.round(surface.rates, 3)

                velocity_true = np.array([[20.076], [1.195], [-0.720]])  # calculated by hand using formula [4].
                velocity_calc = np.round(surface.velocity, 3)

                acceleration_true = np.array([[0.164], [0.252], [-0.223]])  # calculated by hand using formula [4].
                acceleration_calc = np.round(surface.acceleration, 3)

                # 3. tests.
                self.assertTrue(np.allclose(position_true, position_calc, rtol=margin))
                self.assertTrue(np.allclose(angles_true, angles_calc, rtol=margin))
                self.assertTrue(np.allclose(rates_true, rates_calc, rtol=margin))
                self.assertTrue(np.allclose(velocity_true, velocity_calc, rtol=margin))
                self.assertTrue(np.allclose(acceleration_true, acceleration_calc, rtol=margin))

            def test_translated():
                surface = self.load_surface(vertical=True, translated=True)

                # 1. attributes dependent on translation and rotation of surface.
                position_true = np.array([[0.040], [-0.050], [-0.500]])  # calculated by hand using formula [4.1].
                position_calc = np.round(surface.position, 3)

                # 2. tests.
                self.assertTrue(np.allclose(position_true, position_calc, rtol=margin))

            test_untranslated()
            test_translated()

        test_horizontal()
        test_vertical()

    def test_calculate_local_velocity(self, margin=0.01):
        surface = self.load_surface()
        surface_variables = surface.calculate_local_velocity

        # 1. attributes dependant on coordinate and turn rates.
        local_velocity_true = np.array([[19.980], [0.832], [2.124]])
        local_velocity_calc = surface_variables[0]

        local_magnitude_true = 20.110
        local_magnitude_calc = surface_variables[1]

        local_alpha_true = 0.106
        local_alpha_calc = surface_variables[2]

        local_beta_true = 0.042
        local_beta_calc = surface_variables[3]

        # 2. tests.
        self.assertTrue(np.allclose(local_velocity_true, local_velocity_calc, rtol=margin))

        self.assertAlmostEqual(local_magnitude_true, local_magnitude_calc, delta=margin * local_magnitude_true)
        self.assertAlmostEqual(local_alpha_true, local_alpha_calc, delta=margin * local_alpha_true)
        self.assertAlmostEqual(local_beta_true, local_beta_calc, delta=margin * local_beta_true)

    def test_calculate_local_forces(self):
        # surface = self.load_surface()
        # surface_variables = surface.calculate_local_forces
        pass

    def test_calculate_forces(self):

        def test_symmetric():
            surface_symmetric = None

        def test_asymmetric():
            surface_asymmetric = None

        test_symmetric()
        test_asymmetric()

    def test_calculate_moments(self):

        def test_untranslated():
            surface = None

        def test_translated():
            surface = None

        test_untranslated()
        test_translated()

    def test_output(self):
        pass

    def test_transformation_output(self):
        pass


if __name__ == '__main__':
    unittest.main()
