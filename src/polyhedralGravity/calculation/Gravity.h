#pragma once

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityResult.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "spdlog/spdlog.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"

namespace polyhedralGravity {

    /**
     * Alias for x, y, z coordinates.
     */
    using CartesianArray = std::array<double, 3>;
    /**
     * Alias for vectors of cartesian coordinates.
     * @example PlanesVector[i] returns the i-th plane
     */
    using PlanesVector = std::vector<CartesianArray>;
    /**
     * Alias for the segments (always three) of one plane
     * @example SegmentsOfPlaneArray[j] returns the j-th segment
     */
    using SegmentsOfPlaneArray = std::array<CartesianArray, 3>;
    /**
     * Alias for two-dimensional structure of cartesian vectors.
     * The second dimension is fixed to size three.
     * @example SegmentsVector[i][j] returns the j-th segment of the i-th plane
     */
    using SegmentsVector = std::vector<SegmentsOfPlaneArray>;

    /**
     * TODO? Make the whole thing to a namespace, only stamp coupling between methods more practical?
     */
    class Gravity {

        /**
         * The Polyhedron for which we evaluated the gravity at point P
         */
        const Polyhedron _polyhedron;

        /**
         * The constant density of the polyhedron in [kg/m^3].
         * The density is initialized with the default constant density 2670.0 from Tsoulis Paper (above (4)).
         */
        const double _density{util::defaultConstantDensity};

        /**
         * The result of the evaluation of the gravity model
         */
        GravityResult _gravityResult{};

    public:

        /**
         * Construct as new Gravity Calculation with an polyhedron as input
         * @param polyhedron - Poylhedron
         */
        explicit Gravity(const Polyhedron &polyhedron)
                : _polyhedron{polyhedron} {}

        /**
         * Construct as new Gravity Calculation with an polyhedron as input
         * @param polyhedron - Poylhedron
         * @param density - the constant density of the Polyhedron
         */
        Gravity(const Polyhedron &polyhedron, double density)
                : _polyhedron{polyhedron},
                  _density{density} {}

        /**
         * Executes the whole calculation.
         * TODO: return values or via getResult()?
         */
        void calculate();

        /**
         * Calculates the G_ij vectors according to Tsoulis equation (18). Each of these vectors G_ij represents one
         * line segment of the polyhedron.
         *
         * Subscript i stands for the corresponding plane/ face of the polyhedron. The subscript j stands for
         * the specific vector whose endpoints are two of the three vertices making up the plane. The vertices are used
         * in a modulo j fashion, meaning G_i0 is formed as the subtraction of vertice_1 - vertice_0, G_i1 as the
         * subtraction of vertice_2 - vertice_1, and G_i2 as vertice_0 - vertice_2.
         *
         * The dimension of i will be equal to the number of faces, whereas the dimension j will be equal to 3 as the
         * given polyhedral's faces always consist of three segments/ nodes (triangles).
         * @return G_ij vectors
         */
        SegmentsVector calculateGij();

        /**
         * Calculate the plane unit normals N_i vectors according to Tsoulis equation (19).
         * Geometrically, N_i stands perpendicular on plane i.
         *
         * The dimension of i will be equal to the number of faces.
         * @param g - the G_ij vectors
         * @return plane unit normals
         */
        PlanesVector calculatePlaneUnitNormals(const SegmentsVector &g);

        /**
         * Calculates the segment unit normals n_ij according to Tsoulis equation (20).
         * Geometrically, n_ij stands perpendicular on the segment ij of plane i and segment j. These segments j
         * are build up from the G_ij vector and therefor identically numbered.
         *
         * The dimension of i will be equal to the number of faces, whereas the dimension j mirrors the number of
         * segments forming one face. Since we always use triangles, j will be 3.
         * @param g - the G_ij vectors
         * @param planeUnitNormals - the plane unit normals
         * @return segment unit normals
         */
        SegmentsVector calculateSegmentUnitNormals(const SegmentsVector &g, const PlanesVector &planeUnitNormals);

        /**
         * TODO? Maybe do this just in time instead of calculating everything at once and storing in a vector
         * Calculates the plane normal orientations, sigma_p, according to equation (21).
         * The sigma_p values represents the relative position of computation point P with respect to the
         * pointing direction of N_p. E. g. if N_p points to the half-space containing P, the inner product of
         * N_p and -G_i1 will be positive, leading to a negative sigma_p.
         *
         * In equation (21), the used -G_i1 corresponds to opposite position vector of the first vertices building
         * the plane i.
         * @param planeUnitNormals - the plane unit normals
         * @return sigma_p
         */
        std::vector<double> calculatePlaneNormalOrientations(const PlanesVector &planeUnitNormals);


        /**
         * Transforms the edges of the polyhedron to the Hessian Plane form.
         * @param p - the reference point for which the transformation should be executed (default origin {0, 0, 0})
         * @return vector of Hessian Normal Planes
         */
        std::vector<HessianPlane> calculateFacesToHessianPlanes(const CartesianArray &p = {0, 0, 0});

        /**
         * TODO Inline?
         * Calculates the Hessian Plane form spanned by three given points p, q, and r.
         * @param p - first point on the plane
         * @param q - second point on the plane
         * @param r - third point on the plane
         * @param origin - default {0, 0, 0}, but reference for the Hessian Plane form can be adapted
         * @return HessianPlane
         * @related Cross-Product method https://tutorial.math.lamar.edu/classes/calciii/eqnsofplanes.aspx
         */
        HessianPlane computeHessianPlane(const CartesianArray &p, const CartesianArray &q,
                                         const CartesianArray &r, const CartesianArray &origin = {0.0, 0.0, 0.0});


        /**
         * Calculates the plane distances h_p of computation point P from each plane S_p
         * according to the following equation:
         * h_p = D / sqrt(A^2+B^2+C^2)
         *
         * Geometrically, h_p is the distance from P' (the orthogonal projection of computation point P
         * onto the plane S_p) to computation point P.
         * @return plane distances h_p
         */
        std::vector<double> calculatePlaneDistances(const std::vector<HessianPlane> &plane);

        /**
         * Calculates the origins P' for each plane S_p according to equation (22) of Tsoulis paper.
         * P' is the orthogonal projection of the computation point P onto the plane S_p. S_p is the p-th
         * plane, i.e the p-th face of the polyhedron.
         * @param hessianPlanes - the Hessian Plane Form for every plane
         * @param planeUnitNormals - the plane unit normals N_i for every plane
         * @param planeDistances - the plane distance h_p for every plane
         * @return P' for each plane S_p in a vector
         */
        PlanesVector calculateOrthogonalProjectionPointsOnPlane(const std::vector<HessianPlane> &hessianPlanes,
                                                                const PlanesVector &planeUnitNormals,
                                                                const std::vector<double> &planeDistances);


        /**
         * Calculates the segment normal orientations, sigma_pq, according to equation (23).
         * These values represent the orientations of the segment normals n_ij.
         * E.g. if sigma_pq is -1 then n_ij points to the half-plane containing the orthogonal projection Point P'_i.
         * If sigma_pq is 1 then P'_i resides in the other half-space and in case of 0, P'_i lies on the line of the
         * segment G_ij.
         *
         * (G_ij is the notation for a segment)
         * (One can exchange i and p, as well as j and q)
         * @param segmentUnitNormals - the segment unit normal n_ij
         * @param orthogonalProjectionPointsOnPlane - the orthogonal projection points P'_i of P on each plane i
         * @return sigma_pq
         */
        std::vector<std::array<double, 3>> calculateSegmentNormalOrientations(
                const SegmentsVector &segmentUnitNormals, const PlanesVector &orthogonalProjectionPointsOnPlane);

        /**
         * Calculates the origins P'' for each line segment G_pq according to equation (24), (25) and (26) of Tsoulis
         * paper. P'' is the orthogonal projection of the point P' onto the straight line defined by the line
         * segment G_pq.
         * @param orthogonalProjectionPointsOnPlane - the P' for every plane
         * @return the P'' for every line segment of the polyhedron
         */
        SegmentsVector calculateOrthogonalProjectionPointsOnSegments(
                const PlanesVector &orthogonalProjectionPointsOnPlane,
                const std::vector<std::array<double, 3>> &segmentNormalOrientation);

        /**
         * Calculates the point P'' for a given Segment consisting of vertices v1 and v2 and the orthogonal projection
         * point P' for the plane consisting of those vertices. Solves the three equations given in (24), (25) and (26).
         * @param v1 - first endpoint of segment
         * @param v2 - second endpoint of segment
         * @param pPrime - the orthogonal projection P' of P on this plane
         * @return P'' for this segment
         * @note If sigma_pq is zero then P'' == P', this is not checked by this method, but has to be assured first
         */
        CartesianArray calculateOrthogonalProjectionOnSegment(const CartesianArray &v1, const CartesianArray &v2,
                                                              const CartesianArray &pPrime);
        /**
         * Calculates the distance h_pg between the orthogonal projection P' of the computation point P
         * for a given plane and the orthogonal projection P'' of P' for a line segment.
         * @param orthogonalProjectionPointsOnPlane - the P' for every plane
         * @param orthogonalProjectionPointsOnSegment - the P'' for every segment
         * @return a two-dimensional vector of the distances h_pq
         */
        std::vector<std::array<double, 3>> calculateSegmentDistances(
                const PlanesVector &orthogonalProjectionPointsOnPlane,
                const SegmentsVector &orthogonalProjectionPointsOnSegment);

    };

}