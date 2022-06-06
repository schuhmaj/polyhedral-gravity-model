from sympy import *

init_printing()

GRAVITATIONAL_CONSTANT = 6.67430e-11
DEFAULT_CONSTANT_DENSITY = 1.0

X = 0
Y = 0
Z = 0

u = symbols('u')
x1, x2, x3 = symbols('x y z')
v = x1 * x2 * x3
r = sqrt(x1*x1 + x2*x2 + x3*x3)

inner_sum = ((v / u) * log(u + r) - (u * u / 2) * atan(v / (u * u * r)))

total_inner = inner_sum.subs({u: x1}) + inner_sum.subs({u: x2}) + inner_sum.subs({u: x3})

outer_sum = total_inner.evalf(subs={x1: 1-X}) - total_inner.evalf(subs={x1: -1-X})

outer_sum = outer_sum.evalf(subs={x2: 1-Y}) - outer_sum.evalf(subs={x2: -1-Y})

outer_sum = outer_sum.evalf(subs={x3: 1-Z}) - outer_sum.evalf(subs={x3: -1-Z})

print(outer_sum * GRAVITATIONAL_CONSTANT * DEFAULT_CONSTANT_DENSITY)
