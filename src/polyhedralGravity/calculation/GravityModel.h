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
#include "thrust/transform_reduce.h"
#include "thrust/execution_policy.h"
#include "polyhedralGravity/util/UtilityThrust.h"
#include "polyhedralGravity/output/Logging.h"

namespace polyhedralGravity {

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
     * @note Naming scheme corresponds to the following:
     * evaluate()           --> main Method for evaluating the gravity model
     * compute*ForPlane()   --> methods which calculate a property of ONE plane
     * compute*ForSegment() --> methods which calculate a property of ONE segment
     *
     * Typically a calculate() functions calls repeatedly foreach plane p the compute*ForPlane() method which then calls
     * the compute*ForSegment() foreach segment q (triangle faces --> three times)
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
        Array3Triplet
        computeSegmentUnitNormalForPlane(const Array3Triplet &segmentVectors, const Array3 &planeUnitNormal);

        /**
         * Computes the plane unit normal orientation sigma_p for one plane p of the polyhedron
         * according to Tsoulis (21).
         * The plane unit normal orientation values represents the relative position of computation point P
         * with respect to the pointing direction of N_p. E. g. if N_p points to the half-space containing P, the
         * inner product of N_p and -G_i1 will be positive, leading to a negative sigma_p.
         * If sigma_p is zero than P and P' lie geometrically in the same plane --> P == P'.
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
         * @return HessianPlane
         * @related Cross-Product method https://tutorial.math.lamar.edu/classes/calciii/eqnsofplanes.aspx
         */
        HessianPlane computeHessianPlane(const Array3 &p, const Array3 &q, const Array3 &r);

        /**
         * Calculates the plane distances h_p of computation point P to the plane S_p given in Hessian Form
         * according to the following equation:
         * h_p = D / sqrt(A^2+B^2+C^2)
         * @param hessianPlane - Hessian Plane Form of S_p
         * @return plane distance h_p
         */
        double computePlaneDistanceForPlane(const HessianPlane &hessianPlane);

        /**
         * Computes P' for a given plane p according to equation (22) of Tsoulis paper.
         * P' is the orthogonal projection of the computation point P onto the plane S_p.
         * @param planeUnitNormal - the plane unit normal N_p
         * @param planeDistance - the distance from P to the plane h_p
         * @param hessianPlane - the Hessian Plane Form
         * @return P' for this plane
         */
        Array3 computeOrthogonalProjectionPointsOnPlaneForPlane(
                const Array3 &planeUnitNormal,
                double planeDistance,
                const HessianPlane &hessianPlane);

        /**
         * Computes the segment normal orientations sigma_pq for a given plane p.
         * @param vertices - the vertices of this plane
         * @param projectionPointOnPlane - the projection point P' for this plane
         * @param segmentUnitNormalsForPlane - the segment unit normals sigma_pq for this plane
         * @return the segment normal orientations for the plane p
         */
        Array3 computeSegmentNormalOrientationsForPlane(const Array3Triplet &vertices,
                                                        const Array3 &projectionPointOnPlane,
                                                        const Array3Triplet &segmentUnitNormalsForPlane);

        /**
         * Computes the orthogonal projection Points P'' foreach segment q of a given plane p.
         * @param projectionPointOnPlane - the projection Point P'
         * @param segmentNormalOrientations - the segment normal orientations sigma_pq for this plane p
         * @param face - the vertices of the plane p
         * @return the orthogonal projection points of P on the segment P'' foreach segment q of p
         */
        Array3Triplet computeOrthogonalProjectionPointsOnSegmentsForPlane(
                const Array3 &projectionPointOnPlane,
                const Array3 &segmentNormalOrientations,
                const Array3Triplet &face);

        /**
         * Calculates the point P'' for a given Segment consisting of vertices v1 and v2 and the orthogonal projection
         * point P' for the plane consisting of those vertices. Solves the three equations given in (24), (25) and (26).
         * @param vertex1 - first endpoint of segment
         * @param vertex2 - second endpoint of segment
         * @param orthogonalProjectionPointOnPlane - the orthogonal projection P' of P on this plane
         * @return P'' for this segment
         * @note If sigma_pq is zero then P'' == P', this is not checked by this method, but has to be assured first
         */
        Array3 computeOrthogonalProjectionOnSegmentForSegment(const Array3 &vertex1, const Array3 &vertex2,
                                                              const Array3 &orthogonalProjectionPointOnPlane);

        /**
         * Computes the segment distances h_pq between P' for a given plane p and P'' for a given segment q of plane p.
         * @param orthogonalProjectionPointOnPlane - the orthogonal projection point P' for p
         * @param orthogonalProjectionPointOnSegments - the orthogonal projection points P'' for each segment q of p
         * @return distances h_pq for plane p
         */
        Array3 computeSegmentDistancesForPlane(
                const Array3 &orthogonalProjectionPointOnPlane,
                const Array3Triplet &orthogonalProjectionPointOnSegments);

        /**
         * Computes the 3D distances between l1_pq and l2_pq between the computation point P and the line
         * segment endpoints of each polyhedral segment for one plane.
         * Computes the 1D distances s1_pq and s2_pq between orthogonal projection of P on the line
         * segment P''_pq and the line segment endpoints for each polyhedral segment for one plane
         * @param segmentVectorsForPlane - the segment vectors G_pq for plane p
         * @param orthogonalProjectionPointsOnSegmentForPlane - the orthogonal projection Points P'' for plane p
         * @param face - the vertices of plane p
         * @return distances l1_pq and l2_pq and s1_pq and s2_pq foreach segment q of plane p
         */
        std::array<Distance, 3> computeDistancesForPlane(
                const Array3Triplet &segmentVectorsForPlane,
                const Array3Triplet &orthogonalProjectionPointsOnSegmentForPlane,
                const Array3Triplet &face);

        /**
         * Calculates the Transcendental Expressions LN_pq and AN_pq for every line segment of the polyhedron for
         * a given plane p.
         * LN_pq is calculated according to (14) using the natural logarithm and AN_pq is calculated according
         * to (15) using the arctan.
         * @param distancesForPlane - the distances l1, l2, s1, s2 foreach segment q of plane p
         * @param planeDistance - the plane distance h_p for plane p
         * @param segmentDistancesForPlane - the segment distance h_pq for segment q of plane p
         * @param segmentNormalOrientationsForPlane - the segment normal orientations n_pq for a plane p
         * @param orthogonalProjectionPointOnPlane - the orthogonal projection point P' for plane p
         * @param face - the vertices of plane p
         * @return LN_pq and AN_pq foreach segment q of plane p
         */
        std::array<TranscendentalExpression, 3>
        computeTranscendentalExpressionsForPlane(const std::array<Distance, 3> &distancesForPlane,
                                                 double planeDistance, const Array3 &segmentDistancesForPlane,
                                                 const Array3 &segmentNormalOrientationsForPlane,
                                                 const Array3 &projectionPointVertexNorms);

        /**
         * Calculates the singularities (correction) terms according to the Flow text for a given plane p.
         * @param segmentVectorsForPlane - the segment vectors for a given plane
         * @param segmentNormalOrientationForPlane - the segment orientation sigma_pq
         * @param projectionPointVertexNorms - the projection point P'
         * @param planeUnitNormal - the plane unit normal N_p
         * @param planeDistance - the plane distance h_p
         * @param planeNormalOrientation - the plane normal orientation sigma_p
         * @param face - the vertices of plane p
         * @return the singularities for a plane p
         */
        std::pair<double, Array3> computeSingularityTermsForPlane(const Array3Triplet &segmentVectorsForPlane,
                                                                  const Array3 &segmentNormalOrientationForPlane,
                                                                  const Array3 &projectionPointVertexNorms,
                                                                  const Array3 &planeUnitNormal,
                                                                  double planeDistance,
                                                                  double planeNormalOrientation);

        /**
         * Computes the norms of the orthogonal projection point P' on a plane p with each vertex of that plane p.
         * The values are later used to determine if P' is situated at a vertex.
         * @param orthogonalProjectionPointOnPlane - the orthogonal projection point P'
         * @param face - the vertices of plane p
         * @return the norms of p and each vertex
         */
        Array3 computeOrthogonalProjectionPointVertexNormForPlane(
                const Array3 &orthogonalProjectionPointOnPlane,
                const Array3Triplet &face);


        /**
         * An iterator transforming the polyhedron's coordinates on demand by a given offset.
         * This function returns a pair of transform iterators (first --> begin(), second --> end()).
         * @param polyhedron - reference to the polyhedron
         * @param offset - the offset to apply
         * @return pair of transform iterators
         */
        inline auto transformPolyhedron(const Polyhedron &polyhedron, const Array3 &offset) {
            //The offset is captured by value to ensure its lifetime!
            const auto lambdaOffsetApplication = [&polyhedron, offset](
                    const std::array<size_t, 3> &face) -> Array3Triplet {
                using namespace util;
                return {polyhedron.getVertex(face[0]) - offset,
                        polyhedron.getVertex(face[1]) - offset,
                        polyhedron.getVertex(face[2]) - offset};
            };
            auto first = thrust::make_transform_iterator(polyhedron.getFaces().begin(), lambdaOffsetApplication);
            auto last = thrust::make_transform_iterator(polyhedron.getFaces().end(), lambdaOffsetApplication);
            return std::make_pair(first, last);
        }

    };

}