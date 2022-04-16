#pragma once

#include <utility>
#include <array>
#include <vector>
#include <algorithm>
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "spdlog/spdlog.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"
#include "polyhedralGravity/util/UtilityThrust.h"

namespace polyhedralGravity {

    /**
     * Alias for an array of size 2
     * @example for pair's of the same type
     */
    using Array2 = std::array<double, 2>;
    /**
     * Alias for an array of size 3
     * @example for x, y, z coordinates.
     */
    using Array3 = std::array<double, 3>;
    /**
     * Alias for a triplet of arrays of size 3
     * @example for the segment of a triangular face
     */
    using Array3Triplet = std::array<Array3, 3>;

    /**
     * Namespace containing the methods used to evaluate the polyhedrale Gravity Model
     * TODO: insert ::detail namespace for 'private' methods
     */
    namespace GravityModel {

        /**
         * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
         * point P.
         * @param polyhedron - the polyhedron consisting of vertices and triangular faces
         * @param density - the constant density in [kg/m^3]
         * @param computationPoint - the computation Point P (default: {0,0,0})
         * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration
         * at computation Point P
         */
        GravityModelResult evaluate(
                const Polyhedron &polyhedron,
                double density,
                const Array3 &computationPoint = {0.0, 0.0, 0.0});

        /**
         * Calculates the segment vectors G_ij according to Tsoulis equation (18).
         * Each of these vectors G_ij represents one line segment of the polyhedron.
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
        std::vector<Array3Triplet> calculateSegmentVectors(const Polyhedron &polyhedron);

        /**
         * Calculate the plane unit normals N_i vectors according to Tsoulis equation (19).
         * Geometrically, N_i stands perpendicular on plane i.
         *
         * The dimension of i will be equal to the number of faces.
         * @param segmentVectors - the G_ij vectors of each segment
         * @return plane unit normals
         */
        std::vector<Array3> calculatePlaneUnitNormals(const std::vector<Array3Triplet> &segmentVectors);

        /**
         * Calculates the segment unit normals n_ij according to Tsoulis equation (20).
         * Geometrically, n_ij stands perpendicular on the segment ij of plane i and segment j. These segments j
         * are build up from the G_ij vector and therefor identically numbered.
         *
         * The dimension of i will be equal to the number of faces, whereas the dimension j mirrors the number of
         * segments forming one face. Since we always use triangles, j will be 3.
         * @param segmentVectors - the G_ij vectors of each segment
         * @param planeUnitNormals - the plane unit normals
         * @return segment unit normals
         */
        std::vector<Array3Triplet> calculateSegmentUnitNormals(const std::vector<Array3Triplet> &segmentVectors,
                                                               const std::vector<Array3> &planeUnitNormals);

        /**
         * Calculates the plane normal orientations, sigma_p, according to equation (21).
         * The sigma_p values represents the relative position of computation point P with respect to the
         * pointing direction of N_p. E. g. if N_p points to the half-space containing P, the inner product of
         * N_p and -G_i1 will be positive, leading to a negative sigma_p.
         *
         * In equation (21), the used -G_i1 corresponds to opposite position vector of the first vertices building
         * the plane i.
         * @param face - the plane unit normals
         * @return sigma_p
         */
        std::vector<double> calculatePlaneNormalOrientations(const Polyhedron &planeUnitNormal,
                                                             const std::vector<Array3> &face);


        /**
         * Transforms the edges of the polyhedron to the Hessian Plane form.
         * @param p - the reference point for which the transformation should be executed (default origin {0, 0, 0})
         * @return vector of Hessian Normal Planes
         */
        std::vector<HessianPlane> calculateFacesToHessianPlanes(
                const Polyhedron &polyhedron,
                const Array3 &p = {0.0, 0.0, 0.0});

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
        std::vector<Array3>
        calculateOrthogonalProjectionPointsOnPlane(const std::vector<HessianPlane> &hessianPlanes,
                                                   const std::vector<Array3> &planeUnitNormals,
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
        std::vector<Array3> calculateSegmentNormalOrientations(const Polyhedron &polyhedron,
                                                               const std::vector<Array3Triplet> &segmentUnitNormals,
                                                               const std::vector<Array3> &orthogonalProjectionPointsOnPlane);


        /**
         * Calculates the origins P'' for each line segment G_pq according to equation (24), (25) and (26) of Tsoulis
         * paper. P'' is the orthogonal projection of the point P' onto the straight line defined by the line
         * segment G_pq.
         * @param orthogonalProjectionPointsOnPlane - the P' for every plane
         * @return the P'' for every line segment of the polyhedron
         */
        std::vector<Array3Triplet> calculateOrthogonalProjectionPointsOnSegments(const Polyhedron &polyhedron,
                                                                                 const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
                                                                                 const std::vector<Array3> &segmentNormalOrientation);

        /**
         * Calculates the point P'' for a given Segment consisting of vertices v1 and v2 and the orthogonal projection
         * point P' for the plane consisting of those vertices. Solves the three equations given in (24), (25) and (26).
         * @param v1 - first endpoint of segment
         * @param v2 - second endpoint of segment
         * @param pPrime - the orthogonal projection P' of P on this plane
         * @return P'' for this segment
         * @note If sigma_pq is zero then P'' == P', this is not checked by this method, but has to be assured first
         */
        Array3 calculateOrthogonalProjectionOnSegment(const Array3 &v1, const Array3 &v2,
                                                      const Array3 &pPrime);

        /**
         * Calculates the distance h_pg between the orthogonal projection P' of the computation point P
         * for a given plane and the orthogonal projection P'' of P' for a line segment.
         * @param orthogonalProjectionPointsOnPlane - the P' for every plane
         * @param orthogonalProjectionPointsOnSegment - the P'' for every segment
         * @return a two-dimensional vector of the distances h_pq
         */
        std::vector<Array3> calculateSegmentDistances(
                const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
                const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment);

        /**
         * Calculates the 3D distances between l1_pq and l2_pq between the computation point P and the line
         * segment endpoints of each polyhedral segment.
         * Calculates the 1D distances s1_pq and s2_pq between orthogonal projection of P on the line
         * segment P''_pq and the line segment endpoints for each polyhedral segment.
         * The results are stored in the Distance struct for each segment.
         * @param segmentVectors - the segment vectors
         * @param orthogonalProjectionPointsOnSegment - the P'' for every segment
         * @return Distance struct containing l1, l2, s1, s2
         */
        std::vector<std::array<Distance, 3>>
        calculateDistances(const Polyhedron &polyhedron, const std::vector<Array3Triplet> &segmentVectors,
                           const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment);

        /**
         * TODO Contains enorm duplicate!
         * Calculates the Transcendental Expressions LN_pq and AN_pq for every line segment of the polyhedron.
         * LN_pq is calculated according to (14) using the natural logarithm and AN_pq is calculated according
         * to (15) using the arctan.
         * @param distances - the 3D and 1D distances l1, l2, s1, s2 for every segment
         * @param planeDistances - the plane distances h_p for every plane
         * @param segmentDistances - the segment distances h_pq for every segment
         * @param segmentNormalOrientation - the segment normal orientation sigma_pq for every segment
         * @param orthogonalProjectionPointsOnPlane - the orthogonal Projection Points P' for every plane
         * @return the Transcendental Expressions LN and AN for every segment
         */
        std::vector<std::array<TranscendentalExpression, 3>>
        calculateTranscendentalExpressions(const Polyhedron &polyhedron,
                                           const std::vector<std::array<Distance, 3>> &distances,
                                           const std::vector<double> &planeDistances,
                                           const std::vector<Array3> &segmentDistances,
                                           const std::vector<Array3> &segmentNormalOrientation,
                                           const std::vector<Array3> &orthogonalProjectionPointsOnPlane);

        /**
         * TODO Contains enorm duplicate!
         * Calculates the singularities (correction) terms according to the Flow text.
         * @param segmentVectors - the segment vectors
         * @param segmentNormalOrientation - the segment normal orientations sigma_pq
         * @param orthogonalProjectionPointsOnPlane - the orthogonal projection points per Plane P'
         * @param planeDistances - the plane distances h_p
         * @return the singularities terms
         */
        std::vector<std::pair<double, Array3>>
        calculateSingularityTerms(const Polyhedron &polyhedron, const std::vector<Array3Triplet> &segmentVectors,
                                  const std::vector<Array3> &segmentNormalOrientation,
                                  const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
                                  const std::vector<double> &planeDistances,
                                  const std::vector<double> &planeNormalOrientation,
                                  const std::vector<Array3> &planeUnitNormals);

        /**
         * Computes the segment vectors G_ij for one plane of the polyhedron according to Tsoulis (18).
         * The segment vectors G_ij represent the vector from one vertex of the face to the neighboring vertex and
         * depict every line segment of the triangular face (A-B-C)
         * @param vertex0 - the first vertex A
         * @param vertex1 - the second vertex B
         * @param vertex2 - the third vertex C
         * @return the segment vectors for a plane
         */
        Array3Triplet
        computeSegmentVectorsForPlane(const Array3 &vertex0, const Array3 &vertex1, const Array3 &vertex2);


        /**
         * Computes the plane unit normal N_p for one plane p of the polyhedron according to Tsoulis (19).
         * The plane unit normal is the outward pointing normal of the face from the polyhedron.
         * @param segmentVector1 - first edge
         * @param segmentVector2 - second edge
         * @return plane unit normal
         */
        Array3 computePlaneUnitNormalForPlane(const Array3 &segmentVector1, const Array3 &segmentVector2);

        /**
         * Computes the segment unit normals n_pq for one plane p of the polyhedron according to Tsoulis (20).
         * The segment unit normal n_pq represent the normal of one line segment of a polyhedrale face.
         * @param segmentVectors - the segment vectors of the face G_p(0-2)
         * @param planeUnitNormal - the plane unit normal N_p
         * @return segment unit normals n_pq for plane p with q = {0, 1, 2}
         */
        Array3Triplet computeSegmentUnitNormalForPlane(const Array3Triplet &segmentVectors, const Array3 &planeUnitNormal);

        /**
         * Computes the plane unit normal orientation sigma_p for one plane p of the polyhedron
         * according to Tsoulis (21).
         * The plane unit normal orientation values represents the relative position of computation point P
         * with respect to the pointing direction of N_p. E. g. if N_p points to the half-space containing P, the
         * inner product of N_p and -G_i1 will be positive, leading to a negative sigma_p.
         * @param planeUnitNormal - the plane unit normal N_p
         * @param vertex0 - the first vertex of the plane
         * @return plane normal orientation
         */
        double computePlaneNormalOrientationForPlane(const Array3 &planeUnitNormal, const Array3 &vertex0);

        /**
         * Calculates the Hessian Plane form spanned by three given points p, q, and r.
         * @param p - first point on the plane
         * @param q - second point on the plane
         * @param r - third point on the plane
         * @param origin - default {0, 0, 0}, but reference for the Hessian Plane form can be adapted
         * @return HessianPlane
         * @related Cross-Product method https://tutorial.math.lamar.edu/classes/calciii/eqnsofplanes.aspx
         */
        HessianPlane computeHessianPlane(const Array3 &p, const Array3 &q,
                                         const Array3 &r, const Array3 &origin = {0.0, 0.0, 0.0});

        /**
         * Calculates the plane distances h_p of computation point P to the plane S_p given in Hessian Form
         * according to the following equation:
         * h_p = D / sqrt(A^2+B^2+C^2)
         * @param hessianPlane - Hessian Plane Form of S_p
         * @return plane distance h_p
         */
        double computePlaneDistanceForPlane(const HessianPlane &hessianPlane);

    };

}