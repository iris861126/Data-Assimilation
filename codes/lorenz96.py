"""
Definition of the Lorenz 40-variable (Lorenz 96) model
"""
import numpy as np

def f(t, y, F):
    return (np.roll(y, -1) - np.roll(y, 2)) * np.roll(y, 1) - y + F
