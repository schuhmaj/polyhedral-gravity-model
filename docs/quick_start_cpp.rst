Quick Start C++
===============

As Executable
-------------

General
~~~~~~~

After the build, the gravity model can be run by executing:

.. code-block::

    ./polyhedralGravity <YAML-Configuration-File>

where the YAML-Configuration-File contains the required parameters.
Examples for Configuration Files and Polyhedral Source Files can be
found in this repository in the folder `/example-config/`.

Config File
~~~~~~~~~~~

The configuration should look similar to the given example below.
It is required to specify the source-files of the polyhedron's mesh (more info
about the supported files below), the density
of the polyhedron, and the wished computation points where the
gravity tensor shall be computed.

.. code-block:: yaml

    ---
    gravityModel:
      input:
        polyhedron:                                 #polyhedron source-file(s)
          - "../example-config/data/tsoulis.node"   #.node contains the vertices
          - "../example-config/data/tsoulis.face"   #.face contains the triangular faces
        density: 2670.0                             #constant density in [kg/m^3]
        points:                                     #Location of the computation point(s) P
          - [0, 0, 0]                               #Here it is situated at the origin


Polyhedron Source Files
~~~~~~~~~~~~~~~~~~~~~~~

The implementation supports multiple common mesh formats for
the polyhedral source. These include:

====================== ==================================================== ==================================================================================================================================================
File Suffix            Name                                                 Comment
====================== ==================================================== ==================================================================================================================================================
  `.node` and `.face`                     TetGen's files                     These two files need to be given as a pair to the input. [Documentation of TetGen's files](https://wias-berlin.de/software/tetgen/fformats.html)
        `.mesh`                         Medit's mesh files                   Single file containing every needed mesh information.
        `.ply`          The Polygon File format/ Stanfoard Triangle format   Single file containing every needed mesh information. Blender File Format.
        `.off`                          Object File Format                   Single file containing every needed mesh information.
        `.stl`                       Stereolithography format                Single file containing every needed mesh information. Blender File Format.
====================== ==================================================== ==================================================================================================================================================

**Notice!** Only the ASCII versions of those respective files are supported! This is especially
important for e.g. the `.ply` files which also can be in binary format.

Good tools to convert your Polyhedron to a supported format (also for interchanging
ASCII and binary format) are e.g.:

- `Meshio <https://github.com/nschloe/meshio>`__ for Python
- `OpenMesh <https://openmesh-python.readthedocs.io/en/latest/readwrite.html>`__ for Python


As Library
----------

TODO