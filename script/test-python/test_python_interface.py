import numpy as np
import datetime
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

computation_points = np.zeros((1000, 3))

density = 2670.0

eros_vertices = '../../example-config/data/Eros.node'
eros_face = '../../example-config/data/Eros.face'

start = datetime.datetime.now()
result = pg.evaluate([eros_vertices, eros_face], density, computation_points)
end = datetime.datetime.now()

duration = end - start
average = duration.microseconds / 1000.0

print("Average: {} microsecond".format(average))

# u = 0
# for x in range(-50, 50):
#     for y in (-50, 50):
#         potential, acceleration, tensor = result[u]
#         print("Potential at ({}, {}, {}) = {}".format(x, y, 1.0, potential))
#         u += 1
