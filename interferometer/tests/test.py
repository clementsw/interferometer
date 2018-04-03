# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 18:04:27 2018

@author: clementsw
"""

import interferometer
import numpy as np
from unittest import TestCase

class TestInterferometer(TestCase):
    def test_triangle_interferometer(self):
        U = interferometer.random_unitary(5)
        I = interferometer.triangle_decomposition(U)
        self.assertTrue(abs(np.max(I.unitary_transformation() - U)) < 1e-14)
        
    def test_square_interferometer(self):
        U = interferometer.random_unitary(5)
        I = interferometer.square_decomposition(U)
        self.assertTrue(abs(np.max(I.unitary_transformation() - U)) < 1e-14)
    
