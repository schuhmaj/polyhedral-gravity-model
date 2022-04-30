#include "GravityModel.h"

namespace polyhedralGravity {

    GravityModelResult GravityModel::evaluate(
            const Polyhedron &polyhedron, double density, const Array3 &computationPoint) {
        using namespace util;
        /*
         * Calculate V and Vx, Vy, Vz and Vxx, Vyy, Vzz, Vxy, Vxz, Vyz
         */
        auto polyhedronIterator = transformPolyhedron(polyhedron, computationPoint);

        GravityModelResult result{};
        result = thrust::transform_reduce(//TODO thrust::device,
                polyhedronIterator.first, polyhedronIterator.second, [](const Array3Triplet &face) {
                    using namespace util;
                    Array3Triplet segmentVectors = computeSegmentVectorsForPlane(face[0],
                                                                                 face[1],
                                                                                 face[2]);
                    Array3 planeUnitNormal = computePlaneUnitNormalForPlane(segmentVectors[0],
                                                                            segmentVectors[1]);
                    Array3Triplet segmentUnitNormals = computeSegmentUnitNormalForPlane(
                            segmentVectors,
                            planeUnitNormal);
                    double planeNormalOrientation = computePlaneNormalOrientationForPlane(
                            planeUnitNormal, face[0]);
                    HessianPlane hessianPlane = computeHessianPlane(face[0], face[1],
                                                                    face[2]);
                    double planeDistance = computePlaneDistanceForPlane(hessianPlane);
                    Array3 orthogonalProjectionPointOnPlane =
                            computeOrthogonalProjectionPointsOnPlaneForPlane(
                                    planeUnitNormal, planeDistance, hessianPlane);
                    Array3 segmentNormalOrientations = computeSegmentNormalOrientationsForPlane(
                            face, orthogonalProjectionPointOnPlane, segmentUnitNormals);
                    Array3Triplet orthogonalProjectionPointsOnSegmentsForPlane =
                            computeOrthogonalProjectionPointsOnSegmentsForPlane(
                                    orthogonalProjectionPointOnPlane,
                                    segmentNormalOrientations, face);
                    Array3 segmentDistances = computeSegmentDistancesForPlane(
                            orthogonalProjectionPointOnPlane,
                            orthogonalProjectionPointsOnSegmentsForPlane);
                    std::array<Distance, 3> distances = computeDistancesForPlane(
                            segmentVectors, orthogonalProjectionPointsOnSegmentsForPlane,
                            face);
                    std::array<TranscendentalExpression, 3> transcendentalExpressions =
                            computeTranscendentalExpressionsForPlane(distances, planeDistance,
                                                                     segmentDistances,
                                                                     segmentNormalOrientations,
                                                                     orthogonalProjectionPointOnPlane,
                                                                     face);
                    std::pair<double, Array3> singularities =
                            computeSingularityTermsForPlane(segmentVectors,
                                                            segmentNormalOrientations,
                                                            orthogonalProjectionPointOnPlane,
                                                            planeUnitNormal, planeDistance,
                                                            planeNormalOrientation, face);
                    //Sum 1 - potential and acceleration
                    auto zipIteratorSum1PA = util::zipPair(segmentNormalOrientations,
                                                           segmentDistances,
                                                           transcendentalExpressions);
                    const double sum1PA = std::accumulate(zipIteratorSum1PA.first,
                                                          zipIteratorSum1PA.second,
                                                          0.0, [](double acc,
                                                                  const auto &tuple) {
                                const double &sigma_pq = thrust::get<0>(tuple);
                                const double &h_pq = thrust::get<1>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<2>(
                                        tuple);
                                return acc + sigma_pq * h_pq *
                                             transcendentalExpressions.ln;
                            });

                    //Sum 1 - tensor
                    auto zipIteratorSum1T = util::zipPair(segmentUnitNormals,
                                                          transcendentalExpressions);
                    const Array3 sum1T = std::accumulate(
                            zipIteratorSum1T.first, zipIteratorSum1T.second,
                            Array3{0.0, 0.0, 0.0},
                            [](const Array3 &acc, const auto &tuple) {
                                const Array3 &npq = thrust::get<0>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(
                                        tuple);
                                return acc + (npq * transcendentalExpressions.ln);
                            });

                    //Sum 2 - for both the same
                    auto zipIteratorSum2 = util::zipPair(segmentNormalOrientations,
                                                         transcendentalExpressions);
                    const double sum2 = std::accumulate(zipIteratorSum2.first,
                                                        zipIteratorSum2.second,
                                                        0.0,
                                                        [](double acc, const auto &tuple) {
                                                            const double &sigma_pq = thrust::get<0>(
                                                                    tuple);
                                                            const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(
                                                                    tuple);
                                                            return acc + sigma_pq *
                                                                         transcendentalExpressions.an;
                                                        });

                    //Sum up everything (potential and acceleration)
                    const double planeSumPA =
                            sum1PA + planeDistance * sum2 + singularities.first;

                    //Sum up everything (tensor)
                    const Array3 subSum =
                            (sum1T + (planeUnitNormal * (planeNormalOrientation * sum2))) +
                            singularities.second;
                    const Array3 first = planeUnitNormal * subSum;
                    const Array3 reorderedNp = {planeUnitNormal[0], planeUnitNormal[0],
                                                planeUnitNormal[1]};
                    const Array3 reorderedSubSum = {subSum[1], subSum[2], subSum[2]};
                    const Array3 second = reorderedNp * reorderedSubSum;

                    //Accumulate and return
                    return GravityModelResult{
                            planeNormalOrientation * planeDistance * planeSumPA,
                            planeUnitNormal * planeSumPA,
                            concat(first, second)
                    };
                }, result, [](const GravityModelResult &a, const GravityModelResult &b) {
                    return GravityModelResult {
                        a.gravitationalPotential + b.gravitationalPotential,
                        a.gravitationalPotentialDerivative + b.gravitationalPotentialDerivative,
                        a.gradiometricTensor + b.gradiometricTensor
                    };
                });

        const double prefix = util::GRAVITATIONAL_CONSTANT * density;

        result.gravitationalPotential = (result.gravitationalPotential * prefix) / 2.0;
        result.gravitationalPotentialDerivative = abs(result.gravitationalPotentialDerivative * prefix);
        result.gradiometricTensor = result.gradiometricTensor * prefix;
        SPDLOG_INFO("V= {}", result.gravitationalPotential);
        SPDLOG_INFO("Vx= {}", result.gravitationalPotentialDerivative[0]);
        SPDLOG_INFO("Vy= {}", result.gravitationalPotentialDerivative[1]);
        SPDLOG_INFO("Vz= {}", result.gravitationalPotentialDerivative[2]);
        SPDLOG_INFO("Vxx= {}", result.gradiometricTensor[0]);
        SPDLOG_INFO("Vyy= {}", result.gradiometricTensor[1]);
        SPDLOG_INFO("Vzz= {}", result.gradiometricTensor[2]);
        SPDLOG_INFO("Vxy= {}", result.gradiometricTensor[3]);
        SPDLOG_INFO("Vxz= {}", result.gradiometricTensor[4]);
        SPDLOG_INFO("Vyz= {}", result.gradiometricTensor[5]);

        result.p = computationPoint;

        return result;
    }

    Array3Triplet GravityModel::computeSegmentVectorsForPlane(
            const Array3 &vertex0, const Array3 &vertex1, const Array3 &vertex2) {
        using util::operator-;
        //Calculate G_ij
        return {vertex1 - vertex0, vertex2 - vertex1, vertex0 - vertex2};
    }

    Array3 GravityModel::computePlaneUnitNormalForPlane(const Array3 &segmentVector1, const Array3 &segmentVector2) {
        using namespace util;
        //Calculate N_i as (G_i1 * G_i2) / |G_i1 * G_i2| with * being the cross product
        const Array3 crossProduct = cross(segmentVector1, segmentVector2);
        const double norm = euclideanNorm(crossProduct);
        return crossProduct / norm;
    }

    Array3Triplet GravityModel::computeSegmentUnitNormalForPlane(
            const Array3Triplet &segmentVectors, const Array3 &planeUnitNormal) {
        Array3Triplet segmentUnitNormal{};
        //Calculate n_ij as (G_ij * N_i) / |G_ig * N_i| with * being the cross product
        std::transform(segmentVectors.cbegin(), segmentVectors.end(), segmentUnitNormal.begin(),
                       [&planeUnitNormal](const Array3 &segmentVector) -> Array3 {
                           using namespace util;
                           const Array3 crossProduct = cross(segmentVector, planeUnitNormal);
                           const double norm = euclideanNorm(crossProduct);
                           return crossProduct / norm;
                       });
        return segmentUnitNormal;
    }

    double GravityModel::computePlaneNormalOrientationForPlane(const Array3 &planeUnitNormal, const Array3 &vertex0) {
        using namespace util;
        //Calculate N_i * -G_i1 where * is the dot product and then use the inverted sgn
        //We abstain on the double multiplication with -1 in the line above and beyond since two
        //times multiplying with -1 equals no change
        return sgn(dot(planeUnitNormal, vertex0), util::EPSILON);
    }

    HessianPlane GravityModel::computeHessianPlane(const Array3 &p, const Array3 &q, const Array3 &r) {
        using namespace util;
        constexpr Array3 origin{0.0, 0.0, 0.0};
        const auto crossProduct = cross(p - q, p - r);
        const auto res = (origin - p) * crossProduct;
        const double d = res[0] + res[1] + res[2];

        return {crossProduct[0], crossProduct[1], crossProduct[2], d};
    }

    double GravityModel::computePlaneDistanceForPlane(const HessianPlane &hessianPlane) {
        //Compute h_p as D/sqrt(A^2 + B^2 + C^2)
        return std::abs(hessianPlane.d / std::sqrt(
                hessianPlane.a * hessianPlane.a +
                hessianPlane.b * hessianPlane.b +
                hessianPlane.c * hessianPlane.c));
    }

    Array3 GravityModel::computeOrthogonalProjectionPointsOnPlaneForPlane(
            const Array3 &planeUnitNormal,
            double planeDistance,
            const HessianPlane &hessianPlane) {
        using namespace util;
        //Calculate the projection point by (22) P'_ = N_i / norm(N_i) * h_i
        // norm(N_i) is always 1 since N_i is a "normed" vector --> we do not need this division
        Array3 orthogonalProjectionPoint = planeUnitNormal * planeDistance;

        //Calculate alpha, beta and gamma as D/A, D/B and D/C (Notice that we "forget" the minus before those
        // divisions. In consequence, the conditions for signs are reversed below!!!)
        // These values represent the intersections of each polygonal plane with the axes
        Array3 intersections = {hessianPlane.a == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.a,
                                hessianPlane.b == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.b,
                                hessianPlane.c == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.c};

        //Determine the signs of the coordinates of P' according to the intersection values alpha, beta, gamma
        // denoted as __ below, i.e. -alpha, -beta, -gamma denoted -__
        for (unsigned int index = 0; index < 3; ++index) {
            if (intersections[index] < 0) {
                //If -__ >= 0 --> __ < 0 then coordinates are positive, we calculate abs(orthogonalProjectionPoint[..])
                orthogonalProjectionPoint[index] = std::abs(orthogonalProjectionPoint[index]);
            } else {
                //The coordinates need to be controlled
                if (planeUnitNormal[index] > 0) {
                    //If -__ < 0 --> __ >= 0 then the coordinate is negative -orthogonalProjectionPoint[..]
                    orthogonalProjectionPoint[index] = -1.0 * orthogonalProjectionPoint[index];
                } else {
                    //Else the coordinate is positive orthogonalProjectionPoint[..]
                    orthogonalProjectionPoint[index] = orthogonalProjectionPoint[index];
                }
            }
        }
        return orthogonalProjectionPoint;
    }

    Array3 GravityModel::computeSegmentNormalOrientationsForPlane(
            const Array3Triplet &vertices,
            const Array3 &projectionPointOnPlane,
            const Array3Triplet &segmentUnitNormalsForPlane) {
        using namespace util;
        std::array<double, 3> segmentNormalOrientations{};
        //Equation (23)
        //Calculate x_P' - x_ij^1 (x_P' is the projectionPoint and x_ij^1 is the first vertices of one segment,
        //i.e. the coordinates of the training-planes' nodes --> projectionPointOnPlane - vertex
        //Calculate n_ij * x_ij with * being the dot product and use the inverted sgn to determine the value of sigma_pq
        //--> sgn((dot(segmentNormalOrientation, projectionPointOnPlane - vertex)), util::EPSILON) * -1.0
        std::transform(segmentUnitNormalsForPlane.cbegin(), segmentUnitNormalsForPlane.cend(),
                       vertices.cbegin(), segmentNormalOrientations.begin(),
                       [&projectionPointOnPlane](const Array3 &segmentUnitNormal, const Array3 &vertex) {
                           using namespace util;
                           return sgn((dot(segmentUnitNormal, projectionPointOnPlane - vertex)), util::EPSILON) *
                                  -1.0;
                       });
        return segmentNormalOrientations;
    }

    Array3Triplet GravityModel::computeOrthogonalProjectionPointsOnSegmentsForPlane(
            const Array3 &projectionPointOnPlane,
            const Array3 &segmentNormalOrientations,
            const Array3Triplet &face) {
        auto counterJ = thrust::counting_iterator<unsigned int>(0);
        Array3Triplet orthogonalProjectionPointOnSegmentPerPlane{};

        //Running over the segments of this plane
        thrust::transform(segmentNormalOrientations.begin(), segmentNormalOrientations.end(),
                          counterJ, orthogonalProjectionPointOnSegmentPerPlane.begin(),
                          [&](const double &segmentNormalOrientation, const unsigned int j) {
                              //We actually only accept +0.0 or -0.0, so the equal comparison is ok
                              if (segmentNormalOrientation == 0.0) {
                                  //Geometrically trivial case, in neither of the half space --> already on segment
                                  return projectionPointOnPlane;
                              } else {
                                  //In one of the half space, evaluate the projection point P'' for the segment
                                  //with the endpoints v1 and v2
                                  const auto &vertex1 = face[j];
                                  const auto &vertex2 = face[(j + 1) % 3];
                                  return computeOrthogonalProjectionOnSegmentForSegment(vertex1, vertex2,
                                                                                        projectionPointOnPlane);
                              }
                          });
        return orthogonalProjectionPointOnSegmentPerPlane;
    }

    Array3 GravityModel::computeOrthogonalProjectionOnSegmentForSegment(const Array3 &vertex1, const Array3 &vertex2,
                                                                        const Array3 &orthogonalProjectionPointOnPlane) {
        using namespace util;
        //Preparing our the planes/ equations in matrix form
        const Array3 matrixRow1 = vertex2 - vertex1;
        const Array3 matrixRow2 = cross(vertex1 - orthogonalProjectionPointOnPlane, matrixRow1);
        const Array3 matrixRow3 = cross(matrixRow2, matrixRow1);
        const Array3 d = {dot(matrixRow1, orthogonalProjectionPointOnPlane),
                          dot(matrixRow2, orthogonalProjectionPointOnPlane),
                          dot(matrixRow3, vertex1)};
        Matrix<double, 3, 3> columnMatrix = transpose(Matrix<double, 3, 3>{matrixRow1, matrixRow2, matrixRow3});
        //Calculation and solving the equations of above
        const double determinant = det(columnMatrix);
        return Array3{
                det(Matrix<double, 3, 3>{d, columnMatrix[1], columnMatrix[2]}),
                det(Matrix<double, 3, 3>{columnMatrix[0], d, columnMatrix[2]}),
                det(Matrix<double, 3, 3>{columnMatrix[0], columnMatrix[1], d})
        } / determinant;
    }

    Array3 GravityModel::computeSegmentDistancesForPlane(const Array3 &orthogonalProjectionPointOnPlane,
                                                         const Array3Triplet &orthogonalProjectionPointOnSegments) {
        std::array<double, 3> segmentDistances{};
        //The inner loop with the running j --> iterating over the segments
        //Using the values P'_i and P''_ij for the calculation of the distance
        std::transform(orthogonalProjectionPointOnSegments.cbegin(), orthogonalProjectionPointOnSegments.cend(),
                       segmentDistances.begin(), [&](const Array3 &orthogonalProjectionPointOnSegment) {
                    using namespace util;
                    return euclideanNorm(orthogonalProjectionPointOnSegment - orthogonalProjectionPointOnPlane);
                });
        return segmentDistances;
    }

    std::array<Distance, 3> GravityModel::computeDistancesForPlane(
            const Array3Triplet &segmentVectorsForPlane,
            const Array3Triplet &orthogonalProjectionPointsOnSegmentForPlane,
            const Array3Triplet &face) {
        std::array<Distance, 3> distancesForPlane{};
        auto counter = thrust::counting_iterator<unsigned int>(0);
        auto zip = util::zipPair(segmentVectorsForPlane, orthogonalProjectionPointsOnSegmentForPlane);

        thrust::transform(zip.first, zip.second, counter,
                          distancesForPlane.begin(), [&face](const auto &tuple, unsigned int j) {
                    using namespace util;
                    Distance distance{};
                    //segment vector G_pq
                    const Array3 &segmentVector = thrust::get<0>(tuple);
                    //orthogonal projection point on segment P'' for plane p and segment q
                    const Array3 &orthogonalProjectionPointsOnSegment = thrust::get<1>(tuple);

                    //Calculate the 3D distances between P (0, 0, 0) and
                    // the segment endpoints face[j] and face[(j + 1) % 3])
                    distance.l1 = euclideanNorm(face[j]);
                    distance.l2 = euclideanNorm(face[(j + 1) % 3]);
                    //Calculate the 1D distances between P'' (every segment has its own) and
                    // the segment endpoints face[j] and face[(j + 1) % 3])
                    distance.s1 = euclideanNorm(orthogonalProjectionPointsOnSegment - face[j]);
                    distance.s2 = euclideanNorm(orthogonalProjectionPointsOnSegment - face[(j + 1) % 3]);

                    /*
                     * Additional remark:
                     * Details on these conditions are in the second paper referenced in the README.md (Tsoulis, 2021)
                     * The numbering of these conditions is equal to the numbering scheme of the paper
                     * Assign a sign to those magnitudes depending on the relative position of P'' to the two
                     * segment endpoints
                     */

                    //4. Option: |s1 - l1| == 0 && |s2 - l2| == 0 Computation point P is located from the beginning on
                    // the direction of a specific segment (P coincides with P' and P'')
                    if (std::abs(distance.s1 - distance.l1) < EPSILON &&
                        std::abs(distance.s2 - distance.l2) < EPSILON) {
                        //4. Option - Case 2: P is located on the segment from its right side
                        // s1 = -|s1|, s2 = -|s2|, l1 = -|l1|, l2 = -|l2|
                        if (distance.s2 < distance.s1) {
                            distance.s1 *= -1.0;
                            distance.s2 *= -1.0;
                            distance.l1 *= -1.0;
                            distance.l2 *= -1.0;
                            return distance;
                        } else if (std::abs(distance.s2 - distance.s1) < EPSILON) {
                            //4. Option - Case 1: P is located inside the segment (s2 == s1)
                            // s1 = -|s1|, s2 = |s2|, l1 = -|l1|, l2 = |l2|
                            distance.s1 *= -1.0;
                            distance.l1 *= -1.0;
                            return distance;
                        }
                        //4. Option - Case 3: P is located on the segment from its left side
                        // s1 = |s1|, s2 = |s2|, l1 = |l1|, l2 = |l2| --> Nothing to do!
                    } else {
                        const double norm = euclideanNorm(segmentVector);
                        if (distance.s1 < norm && distance.s2 < norm) {
                            //1. Option: |s1| < |G_ij| && |s2| < |G_ij| Point P'' is situated inside the segment
                            // s1 = -|s1|, s2 = |s2|, l1 = |l1|, l2 = |l2|
                            distance.s1 *= -1.0;
                            return distance;
                        } else if (distance.s2 < distance.s1) {
                            //2. Option: |s2| < |s1| Point P'' is on the right side of the segment
                            // s1 = -|s1|, s2 = -|s2|, l1 = |l1|, l2 = |l2|
                            distance.s1 *= -1.0;
                            distance.s2 *= -1.0;
                            return distance;
                        }
                        //3. Option: |s1| < |s2| Point P'' is on the left side of the segment
                        // s1 = |s1|, s2 = |s2|, l1 = |l1|, l2 = |l2| --> Nothing to do!
                    }
                    return distance;
                });
        return distancesForPlane;
    }

    std::array<TranscendentalExpression, 3> GravityModel::computeTranscendentalExpressionsForPlane(
            const std::array<Distance, 3> &distancesForPlane,
            double planeDistance,
            const Array3 &segmentDistancesForPlane,
            const Array3 &segmentNormalOrientationsForPlane,
            const Array3 &orthogonalProjectionPointOnPlane,
            const Array3Triplet &face) {
        std::array<TranscendentalExpression, 3> transcendentalExpressionsForPlane{};

        //Zip iterator consisting of 3D and 1D distances l1/l2 and s1/2 for this plane | h_pq | sigma_pq for this plane
        auto zip = util::zipPair(distancesForPlane, segmentDistancesForPlane, segmentNormalOrientationsForPlane);
        auto counter = thrust::counting_iterator<unsigned int>(0);

        thrust::transform(zip.first, zip.second, counter,
                          transcendentalExpressionsForPlane.begin(), [&](const auto &tuple, const unsigned int j) {
                    using namespace util;
                    //distances l1, l2, s1, s1 for this segment q of plane p
                    const Distance &distance = thrust::get<0>(tuple);
                    //segment distance h_pq for this segment q of plane p
                    const double segmentDistance = thrust::get<1>(tuple);
                    //segment normal orientation sigma_pq for this segment q of plane p
                    const double segmentNormalOrientation = thrust::get<2>(tuple);

                    //Result for this segment
                    TranscendentalExpression transcendentalExpressionPerSegment{};

                    //Aliases for the Vertices (endpoints) of this segment
                    const Array3 &v1 = face[(j + 1) % 3];
                    const Array3 &v2 = face[j];

                    //Compute LN_pq according to (14)
                    //If either sigmaPQ has no sign AND either of the distances of P' to the two
                    //segment endpoints is zero OR the 1D and 3D distances are below some threshold
                    //then LN_pq is zero, too
                    //TODO (distance.s1 + distance.s2 < 1e-10 && distance.l1 + distance.l2 < 1e-10)?!
                    if ((segmentNormalOrientation == 0.0 &&
                         (euclideanNorm(orthogonalProjectionPointOnPlane - v1) == 0.0
                          || euclideanNorm(orthogonalProjectionPointOnPlane - v2) == 0.0)) ||
                        (distance.s1 + distance.s2 < 1e-10 && distance.l1 + distance.l2 < 1e-10)) {
                        transcendentalExpressionPerSegment.ln = 0.0;
                    } else {
                        transcendentalExpressionPerSegment.ln =
                                std::log((distance.s2 + distance.l2) / (distance.s1 + distance.l1));
                    }

                    //Compute AN_pq according to (15)
                    //If h_p or h_pq is zero then AN_pq is zero, too
                    //TODO Comparison might be wrong!!!
                    if (planeDistance == 0 || segmentDistance == 0) {
                        transcendentalExpressionPerSegment.an = 0.0;
                    } else {
                        transcendentalExpressionPerSegment.an =
                                std::atan((planeDistance * distance.s2) / (segmentDistance * distance.l2)) -
                                std::atan((planeDistance * distance.s1) / (segmentDistance * distance.l1));
                    }

                    return transcendentalExpressionPerSegment;
                });
        return transcendentalExpressionsForPlane;
    }

    std::pair<double, Array3> GravityModel::computeSingularityTermsForPlane(
            const Array3Triplet &segmentVectorsForPlane,
            const Array3 &segmentNormalOrientationForPlane,
            const Array3 &orthogonalProjectionPointOnPlane,
            const Array3 &planeUnitNormal,
            double planeDistance,
            double planeNormalOrientation,
            const Array3Triplet &face) {
        //1. case: If all sigma_pq for a given plane p are 1.0 then P' lies inside the plane S_p
        if (std::all_of(segmentNormalOrientationForPlane.cbegin(), segmentNormalOrientationForPlane.cend(),
                        [](const double sigma) { return sigma == 1.0; })) {
            using namespace util;
            return std::make_pair(
                    -1.0 * util::PI2 * planeDistance,
                    planeUnitNormal * (-1.0 * util::PI2 * planeNormalOrientation));
        }
        //2. case: If sigma_pq == 0 AND norm(P' - v1) < norm(G_ij) && norm(P' - v2) < norm(G_ij) with G_ij
        // as the vector of v1 and v2
        // then P' is located on one line segment G_p of plane p, but not on any of its vertices
        auto counter2 = thrust::counting_iterator<unsigned int>(0);
        auto secondCaseBegin = util::zip(segmentVectorsForPlane.begin(), segmentNormalOrientationForPlane.begin(),
                                         counter2);
        auto secondCaseEnd = util::zip(segmentVectorsForPlane.end(), segmentNormalOrientationForPlane.end(),
                                       counter2 + 3);
        if (std::any_of(secondCaseBegin, secondCaseEnd, [&](const auto &tuple) {
            using namespace util;
            const Array3 &segmentVector = thrust::get<0>(tuple);
            const double segmentNormalOrientation = thrust::get<1>(tuple);
            const unsigned int j = thrust::get<2>(tuple);

            if (segmentNormalOrientation != 0.0) {
                return false;
            }

            const Array3 &v1 = face[(j + 1) % 3];
            const Array3 &v2 = face[j];
            const double segmentVectorNorm = euclideanNorm(segmentVector);
            return euclideanNorm(orthogonalProjectionPointOnPlane - v1) < segmentVectorNorm &&
                   euclideanNorm(orthogonalProjectionPointOnPlane - v2) < segmentVectorNorm;
        })) {
            using namespace util;
            return std::make_pair(
                    -1.0 * util::PI * planeDistance,
                    planeUnitNormal * (-1.0 * util::PI * planeNormalOrientation));
        }
        //3. case If sigma_pq == 0 AND norm(P' - v1) < 0 || norm(P' - v2) < 0
        // then P' is located at one of G_p's vertices
        auto counter3 = thrust::counting_iterator<int>(0);
        auto thirdCaseBegin = util::zip(segmentNormalOrientationForPlane.begin(), counter3);
        auto thirdCaseEnd = util::zip(segmentNormalOrientationForPlane.end(), counter3 + 3);
        double e1;
        double e2;
        unsigned int j;
        if (std::any_of(thirdCaseBegin, thirdCaseEnd, [&](const auto &tuple) {
            using namespace util;
            const double segmentNormalOrientation = thrust::get<0>(tuple);
            j = thrust::get<1>(tuple);

            if (segmentNormalOrientation != 0.0) {
                return false;
            }

            const Array3 &v1 = face[(j + 1) % 3];
            const Array3 &v2 = face[j];
            e1 = euclideanNorm(orthogonalProjectionPointOnPlane - v1);
            e2 = euclideanNorm(orthogonalProjectionPointOnPlane - v2);
            return e1 == 0.0 || e2 == 0.0;
        })) {
            using namespace util;
            const Array3 &g1 = e1 == 0.0 ? segmentVectorsForPlane[j] : segmentVectorsForPlane[(j - 1 + 3) % 3];
            const Array3 &g2 = e1 == 0.0 ? segmentVectorsForPlane[(j + 1) % 3] : segmentVectorsForPlane[j];
            const double gdot = dot(g1 * -1.0, g2);
            const double theta =
                    gdot == 0.0 ? util::PI_2 : std::acos(gdot / (euclideanNorm(g1) * euclideanNorm(g2)));
            return std::make_pair(
                    -1.0 * theta * planeDistance,
                    planeUnitNormal * (-1.0 * theta * planeNormalOrientation));
        }
        //4. case Otherwise P' is located outside the plane S_p and then the singularity equals zero
        return std::make_pair(0.0, Array3{0, 0, 0});
    }


}
