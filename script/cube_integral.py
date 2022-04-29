from mpmath import *

"""
Evaluates the analytical solution of the gravity potential
for a cube centered around the origin with edge size 2
for the coordinates (X, Y, Z)
"""

# Constants (not used, since only the actual result of the integral is of interest)
GRAVITATIONAL_CONSTANT = 6.67430e-11
DEFAULT_CONSTANT_DENSITY = 2670.0

# Precision of the calculation
mp.dps = 20

# Coordinates of the point P to evaluate
X = 1
Y = 1
Z = 1

print("Cube with edge length = 2, centered around the origin.")
print("Evaluating gravity potential at ({}, {}, {})".format(X, Y, Z))
print("Precision of the evaluation: {}".format(mp.dps))

# The integral of the gravity potential of a cube
cube_integral = lambda x, y, z: 1.0 / sqrt((X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z))

# The calculation
result = mp.quad(cube_integral, [-1, 1], [-1, 1], [-1, 1])

print()
print("Gravity Potential Result: {}".format(result))
