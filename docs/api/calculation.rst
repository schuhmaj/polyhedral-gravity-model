Calculation
===========

Overview
--------

The calculation is the heart of the Polyhedral Gravity Model
since it contains the important function :code:`evaluate(..)`
used evaluation of the polyhedral gravity model at computation
point(s) P.

Further, it implements the `Möller–Trumbore intersection algorithm <https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm>`__
capable of checking if the plane unit normals of the polyhedron are pointing outwards.


GravityModel
------------

.. doxygennamespace:: polyhedralGravity::GravityModel


MeshChecking
-----------

.. doxygennamespace:: polyhedralGravity::MeshChecking