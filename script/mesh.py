# core stuff
import pickle as pk
import numpy as np

# meshing
import tetgen


def read_pk_file(filename):
    """
    Reads in a .pk file and returns the vertices and triangles (faces)
    :param filename: String
    :return: tuple of vertices and faces
    """
    with open(filename, "rb") as f:
        mesh_points, mesh_triangles = pk.load(f)
    mesh_points = np.array(mesh_points)
    mesh_triangles = np.array(mesh_triangles)
    # Characteristic dimension
    L = max(mesh_points[:, 0]) - min(mesh_points[:, 0])

    # Non dimensional units
    mesh_points = mesh_points / L * 2 * 0.8
    return mesh_points, mesh_triangles


def write_to_node_file(filename, nodes, elems):
    with open(filename, "a") as f:
        f.write("# Node count, 3 dimensions, no attribute, no boundary marker\n")
        # TODO Very ugly division! Can cause problems, fix this
        f.write("{} {} {} {}\n".format(int(nodes.size / 3), 3, 0, 0))
        f.write("# Node index, node coordinates\n")
        index = 1
        for n in nodes:
            f.write("{} {} {} {}\n".format(index, n[0], n[1], n[2]))
            index += 1


def main():
    print("Reading in file")
    mesh_points, mesh_triangles = read_pk_file("../mesh/Eros.pk")
    # print("Tetrahralize")
    # tgen = tetgen.TetGen(mesh_points, mesh_triangles)
    # nodes, elems = tgen.tetrahedralize()
    print("Writing to node file")
    write_to_node_file("../mesh/Eros.node", mesh_points, mesh_triangles)
    # tgen.grid.plot(show_edges=True)
    print("Finished")


if __name__ == "__main__":
    main()
