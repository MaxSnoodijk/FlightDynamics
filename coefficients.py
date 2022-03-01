import numpy as np


def coefficients(example=False):

    def example_coefficients():
        return 1

    def symmetric_coefficients():
        return 1

    def asymmetric_coefficients():
        return 1

    def misc_coefficients():
        return 1

    if example:
        return example_coefficients()
    else:
        symmetric = symmetric_coefficients()
        asymmetric = asymmetric_coefficients()
        misc = misc_coefficients()
        return symmetric, asymmetric, misc
