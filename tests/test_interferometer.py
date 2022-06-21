import numpy as np
from unittest import TestCase

from interferometer import random_unitary, triangle_decomposition, square_decomposition

class TestInterferometer(TestCase):

    def test_triangle_interferometer(self):
        U = random_unitary(5)
        I = triangle_decomposition(U)
        self.assertTrue(abs(np.max(I.calculate_transformation() - U)) < 1e-14)

    def test_square_interferometer(self):
        U = random_unitary(5)
        I = square_decomposition(U)
        self.assertTrue(abs(np.max(I.calculate_transformation() - U)) < 1e-14)
    