# Import causes Segfault?!!!
import numpy as np
import polyhedral_gravity as pg

cube_vertices = [
    [-1, -1, -1],
    [1, -1, -1],
    [1, 1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [1, 1, 1],
    [-1, 1, 1]
]

cube_faces = [
    [1, 3, 2],
    [0, 3, 1],
    [0, 1, 5],
    [0, 5, 4],
    [0, 7, 3],
    [0, 4, 7],
    [1, 2, 6],
    [1, 6, 5],
    [2, 3, 6],
    [3, 7, 6],
    [4, 5, 6],
    [4, 6, 7]
]

tsoulis_vertices = [
    np.array([-20.0, 0.0, 25.0]),
    np.array([0.0, 0.0, 25.0]),
    np.array([0.0, 10.0, 25.0]),
    np.array([-20.0, 10.0, 25.0]),
    np.array([-20.0, 0.0, 15.0]),
    np.array([0.0, 0.0, 15.0]),
    np.array([0.0, 10.0, 15.0]),
    np.array([-20.0, 10.0, 15.0])
]

tsoulis_faces = [
    [0, 1, 3],
    [1, 2, 3],
    [0, 4, 5],
    [0, 5, 1],
    [0, 7, 4],
    [0, 3, 7],
    [1, 5, 6],
    [1, 6, 2],
    [3, 6, 7],
    [2, 6, 3],
    [4, 6, 5],
    [4, 7, 6]
]

computation_points = []

for x in range(-50, 50):
    for y in range(-50, 50):
        computation_points.append((x, y, 0.0))

density = 2670.0
polyhedron = pg.Polyhedron(tsoulis_vertices, tsoulis_faces)

result = pg.evaluate(polyhedron, density, computation_points)

u = 0
for x in range(-50, 50):
    for y in (-50, 50):
        print("Potential at ({}, {}, {}) = {}".format(x, y, 1.0, result[u].potential))
        u += 1
