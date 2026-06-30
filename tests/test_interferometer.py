import numpy as np
from unittest import TestCase
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from interferometer.main import Interferometer, Beamsplitter, triangle_decomposition, square_decomposition, random_unitary

class TestInterferometer(TestCase):

    def test_triangle_interferometer(self):
        U = random_unitary(5)
        I = triangle_decomposition(U)
        self.assertTrue(abs(np.max(I.calculate_transformation() - U)) < 1e-14)

    def test_square_interferometer(self):
        U = random_unitary(5)
        I = square_decomposition(U)
        self.assertTrue(abs(np.max(I.calculate_transformation() - U)) < 1e-14)

    def test_loss(self):
        U = random_unitary(10)
        I = square_decomposition(U)
        U_pure = I.calculate_transformation()
        for bs in I.BS_list:
            bs.alpha = 3.0

        U_loss = I.calculate_transformation()
        input = np.ones(10) * (1/np.sqrt(10))

        output_loss = U_loss @ input
        output_ideal = U_pure @ input
        self.assertTrue(np.sum(np.square(np.abs(output_loss))) < np.sum(np.square(np.abs(output_ideal))))

        

def main():
    test = TestInterferometer()
    test.test_square_interferometer()
    test.test_triangle_interferometer()
    test.test_loss()


if __name__ == "__main__":
    main()