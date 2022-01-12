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


def write_to_node_faces_ele_file(filename, nodes, faces, ele):
    with open(filename + ".node", "a") as f:
        f.write("# Node count, 3 dimensions, no attribute, no boundary marker\n")
        # TODO Very ugly division! Can cause problems, fix this
        f.write("{} {} {} {}\n".format(int(nodes.size / 3), 3, 0, 0))
        f.write("# Node index, node coordinates\n")
        index = 0
        for n in nodes:
            f.write("{} {} {} {}\n".format(index, n[0], n[1], n[2]))
            index += 1
    with open(filename + ".face", "a") as f:
        f.write("# Number of faces, boundary marker off\n")
        # TODO Very ugly division! Can cause problems, fix this
        f.write("{} {}\n".format(int(faces.size / 3), 0))
        f.write("# Face index, nodes of face\n")
        index = 0
        for fac in faces:
            f.write("{} {} {} {}\n".format(index, fac[0], fac[1], fac[2]))
            index += 1
    with open(filename + ".ele", "a") as f:
        f.write("# number of tetrahedra, number of nodes per tet, no region attribute\n")
        # TODO Very ugly division! Can cause problems, fix this
        f.write("{} {} {}\n".format(int(ele.size / 4), 4, 0))
        f.write("# tetrahedra index, nodes\n")
        index = 0
        for tet in ele:
            f.write("{} {} {} {} {}\n".format(index, tet[0], tet[1], tet[2], tet[3]))
            index += 1


def write_to_ele_file(filename, ele):
    with open(filename + ".ele", "a") as f:
        f.write("# number of tetrahedra, number of nodes per tet, no region attribute\n")
        # TODO Very ugly division! Can cause problems, fix this
        f.write("{} {} {}\n".format(int(ele.size / 4), 4, 0))
        f.write("# tetrahedra index, nodes\n")
        index = 0
        for tet in ele:
            f.write("{} {} {} {} {}\n".format(index, tet[0], tet[1], tet[2], tet[3]))
            index += 1


def main():
    print("Reading file...")
    mesh_points, mesh_triangles = read_pk_file("../mesh/Eros.pk")
    print("Tetrahralize...")
    tgen = tetgen.TetGen(mesh_points, mesh_triangles)
    nodes, elems = tgen.tetrahedralize()
    print("Writing to files..")
    tgen.write("../mesh/Eros_python.vtk")
    write_to_node_faces_ele_file("../mesh/Eros", nodes, mesh_triangles, elems)
    # tgen.grid.plot(show_edges=True)
    print("Finished.")


if __name__ == "__main__":
    main()
