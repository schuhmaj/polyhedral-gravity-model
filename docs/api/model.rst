Model
=====

Overview
--------

The Model contains two types of data containers. First, there
is the class :code:`Polyhedron` which consists of a vector
vertices and a vector of triangular faces. Second, there
are Utility Container whose single purpose is improving
code readability in providing named access to their
members.

Polyhedron
----------

.. doxygenclass:: polyhedralGravity::Polyhedron

Utility Container
-----------------

.. doxygenstruct:: polyhedralGravity::Distance

.. doxygenstruct:: polyhedralGravity::TranscendentalExpression

.. doxygenstruct:: polyhedralGravity::HessianPlane

.. doxygenclass:: polyhedralGravity::GravityModelResult