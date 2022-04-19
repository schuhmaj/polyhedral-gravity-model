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
        //TODO thrust::reduce or std::accumulate?
        result = thrust::reduce(polyhedronIterator.first, polyhedronIterator.second,
                                result, [](const GravityModelResult &acc, const Array3Triplet &face) {
                    using namespace util;
                    Array3Triplet segmentVectors = computeSegmentVectorsForPlane(face[0], face[1], face[2]);
                    Array3 planeUnitNormal = computePlaneUnitNormalForPlane(segmentVectors[0], segmentVectors[1]);
                    Array3Triplet segmentUnitNormals = computeSegmentUnitNormalForPlane(segmentVectors,
                                                                                        planeUnitNormal);
                    double planeNormalOrientation = computePlaneNormalOrientationForPlane(planeUnitNormal, face[0]);
                    HessianPlane hessianPlane = computeHessianPlane(face[0], face[1], face[2]);
                    double planeDistance = computePlaneDistanceForPlane(hessianPlane);
                    Array3 orthogonalProjectionPointOnPlane =
                            computeOrthogonalProjectionPointsOnPlaneForPlane(
                                    planeUnitNormal, planeDistance, hessianPlane);
                    Array3 segmentNormalOrientations = computeSegmentNormalOrientationsForPlane(
                            face, orthogonalProjectionPointOnPlane, segmentUnitNormals);
                    Array3Triplet orthogonalProjectionPointsOnSegmentsForPlane =
                            computeOrthogonalProjectionPointsOnSegmentsForPlane(
                                    orthogonalProjectionPointOnPlane, segmentNormalOrientations, face);
                    Array3 segmentDistances = computeSegmentDistancesForPlane(
                            orthogonalProjectionPointOnPlane, orthogonalProjectionPointsOnSegmentsForPlane);
                    std::array<Distance, 3> distances = computeDistancesForPlane(
                            segmentVectors, orthogonalProjectionPointsOnSegmentsForPlane, face);
                    std::array<TranscendentalExpression, 3> transcendentalExpressions =
                            computeTranscendentalExpressionsForPlane(distances, planeDistance, segmentDistances,
                                                                     segmentNormalOrientations,
                                                                     orthogonalProjectionPointOnPlane, face);
                    std::pair<double, Array3> singularities =
                            computeSingularityTermsForPlane(segmentVectors, segmentNormalOrientations,
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
                            zipIteratorSum1T.first, zipIteratorSum1T.second, Array3{0.0, 0.0, 0.0},
                            [](const Array3 &acc, const auto &tuple) {
                                const Array3 &npq = thrust::get<0>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(tuple);
                                return acc + (npq * transcendentalExpressions.ln);
                            });

                    //Sum 2 - for both the same
                    auto zipIteratorSum2 = util::zipPair(segmentNormalOrientations,
                                                         transcendentalExpressions);
                    const double sum2 = std::accumulate(zipIteratorSum2.first,
                                                        zipIteratorSum2.second,
                                                        0.0, [](double acc, const auto &tuple) {
                                const double &sigma_pq = thrust::get<0>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(
                                        tuple);
                                return acc + sigma_pq * transcendentalExpressions.an;
                            });

                    //Sum up everything (potential and acceleration)
                    const double planeSumPA = sum1PA + planeDistance * sum2 + singularities.first;

                    //Sum up everything (tensor)
                    const Array3 subSum =
                            (sum1T + (planeUnitNormal * (planeNormalOrientation * sum2))) + singularities.second;
                    const Array3 first = planeUnitNormal * subSum;
                    const Array3 reorderedNp = {planeUnitNormal[0], planeUnitNormal[0], planeUnitNormal[1]};
                    const Array3 reorderedSubSum = {subSum[1], subSum[2], subSum[2]};
                    const Array3 second = reorderedNp * reorderedSubSum;

                    //Accumulate and return
                    return GravityModelResult{
                            acc.gravitationalPotential + planeNormalOrientation * planeDistance * planeSumPA,
                            acc.gravitationalPotentialDerivative + planeUnitNormal * planeSumPA,
                            acc.gradiometricTensor + concat(first, second)
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

    std::vector<Array3Triplet> GravityModel::calculateSegmentVectors(const Polyhedron &polyhedron) {
        std::vector<Array3Triplet> segmentVectors{polyhedron.countFaces()};
        //Calculate G_ij for every plane given as input the three vertices of every face
        std::transform(polyhedron.getFaces().cbegin(), polyhedron.getFaces().cend(), segmentVectors.begin(),
                       [&polyhedron](const std::array<size_t, 3> &face) -> Array3Triplet {
                           const Array3 &vertex0 = polyhedron.getVertex(face[0]);
                           const Array3 &vertex1 = polyhedron.getVertex(face[1]);
                           const Array3 &vertex2 = polyhedron.getVertex(face[2]);
                           return computeSegmentVectorsForPlane(vertex0, vertex1, vertex2);
                       });
        return segmentVectors;
    }

    std::vector<Array3> GravityModel::calculatePlaneUnitNormals(const std::vector<Array3Triplet> &segmentVectors) {
        std::vector<Array3> planeUnitNormals{segmentVectors.size()};
        //Calculate N_p for every plane given as input: G_i0 and G_i1 of every plane
        std::transform(segmentVectors.cbegin(), segmentVectors.cend(), planeUnitNormals.begin(),
                       [](const Array3Triplet &segmentVectorsForPlane) -> Array3 {
                           return computePlaneUnitNormalForPlane(segmentVectorsForPlane[0], segmentVectorsForPlane[1]);
                       });
        return planeUnitNormals;
    }

    std::vector<Array3Triplet> GravityModel::calculateSegmentUnitNormals(
            const std::vector<Array3Triplet> &segmentVectors,
            const std::vector<Array3> &planeUnitNormals) {
        std::vector<Array3Triplet> segmentUnitNormals{segmentVectors.size()};
        //Loop" over G_i (running i=p) and N_p calculating n_p for every plane
        std::transform(segmentVectors.cbegin(), segmentVectors.cend(), planeUnitNormals.cbegin(),
                       segmentUnitNormals.begin(),
                       [](const Array3Triplet &segmentVectorsForPlane, const Array3 &planeUnitNormal) {
                           return computeSegmentUnitNormalForPlane(segmentVectorsForPlane, planeUnitNormal);
                       });
        return segmentUnitNormals;
    }

    std::vector<double>
    GravityModel::calculatePlaneNormalOrientations(const Array3 &computationPoint, const Polyhedron &polyhedron,
                                                   const std::vector<Array3> &planeUnitNormals) {
        std::vector<double> planeNormalOrientations(planeUnitNormals.size(), 0.0);
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        //Calculate sigma_p for every plane given as input: N_p and vertex0 of every face
        std::transform(planeUnitNormals.cbegin(), planeUnitNormals.cend(), transformedPolyhedronIt.first,
                       planeNormalOrientations.begin(),
                       [](const Array3 &planeUnitNormal, const Array3Triplet &face) {
                           //The first vertices' coordinates of the given face consisting of G_i's
                           return computePlaneNormalOrientationForPlane(planeUnitNormal, face[0]);
                       });
        return planeNormalOrientations;
    }

    std::vector<HessianPlane>
    GravityModel::calculateFacesToHessianPlanes(const Array3 &computationPoint, const Polyhedron &polyhedron) {
        std::vector<HessianPlane> hessianPlanes{polyhedron.countFaces()};
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        //Calculate for each face/ plane/ triangle (here) the Hessian Plane
        std::transform(transformedPolyhedronIt.first, transformedPolyhedronIt.second, hessianPlanes.begin(),
                       [](const Array3Triplet &face) -> HessianPlane {
                           using namespace util;
                           //The three vertices put up the plane, p is the origin of the reference system default 0,0,0
                           return computeHessianPlane(face[0], face[1], face[2]);
                       });
        return hessianPlanes;
    }

    std::vector<double> GravityModel::calculatePlaneDistances(const std::vector<HessianPlane> &plane) {
        std::vector<double> planeDistances(plane.size(), 0.0);
        //For each plane compute h_p
        std::transform(plane.cbegin(), plane.cend(), planeDistances.begin(),
                       [](const HessianPlane &plane) -> double {
                           return computePlaneDistanceForPlane(plane);
                       });
        return planeDistances;
    }

    std::vector<Array3> GravityModel::calculateOrthogonalProjectionPointsOnPlane(
            const std::vector<HessianPlane> &hessianPlanes,
            const std::vector<Array3> &planeUnitNormals,
            const std::vector<double> &planeDistances) {
        std::vector<Array3> orthogonalProjectionPointsOfP{planeUnitNormals.size()};

        //Zip the three required arguments together: Plane normal N_i, Plane Distance h_i and the Hessian Form
        auto zip = util::zipPair(planeUnitNormals, planeDistances, hessianPlanes);

        //Calculates the Projection Point P' for every plane p
        thrust::transform(zip.first, zip.second, orthogonalProjectionPointsOfP.begin(), [](const auto &tuple) {
            using namespace util;
            const Array3 &planeUnitNormal = thrust::get<0>(tuple);
            const double planeDistance = thrust::get<1>(tuple);
            const HessianPlane &hessianPlane = thrust::get<2>(tuple);
            return computeOrthogonalProjectionPointsOnPlaneForPlane(planeUnitNormal, planeDistance, hessianPlane);
        });
        return orthogonalProjectionPointsOfP;
    }

    std::vector<Array3> GravityModel::calculateSegmentNormalOrientations(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3Triplet> &segmentUnitNormals,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane) {
        std::vector<Array3> segmentNormalOrientations{segmentUnitNormals.size()};

        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(transformedPolyhedronIt.first,
                               orthogonalProjectionPointsOnPlane.begin(), segmentUnitNormals.begin());
        auto last = util::zip(transformedPolyhedronIt.second,
                              orthogonalProjectionPointsOnPlane.end(), segmentUnitNormals.end());

        //Calculates the segment normal orientation sigma_pq for every plane p
        thrust::transform(first, last, segmentNormalOrientations.begin(), [](const auto &tuple) {
            const Array3Triplet &face = thrust::get<0>(tuple);
            const Array3 &projectionPointOnPlaneForPlane = thrust::get<1>(tuple);
            const Array3Triplet &segmentUnitNormalsForPlane = thrust::get<2>(tuple);

            return computeSegmentNormalOrientationsForPlane(face, projectionPointOnPlaneForPlane,
                                                            segmentUnitNormalsForPlane);
        });
        return segmentNormalOrientations;
    }

    std::vector<Array3Triplet> GravityModel::calculateOrthogonalProjectionPointsOnSegments(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<Array3> &segmentNormalOrientation) {
        std::vector<Array3Triplet> orthogonalProjectionPointsOnSegments{orthogonalProjectionPointsOnPlane.size()};

        //Zip the three required arguments together: P' for every plane, sigma_pq for every segment, the faces
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(orthogonalProjectionPointsOnPlane.begin(), segmentNormalOrientation.begin(),
                               transformedPolyhedronIt.first);
        auto last = util::zip(orthogonalProjectionPointsOnPlane.end(), segmentNormalOrientation.end(),
                              transformedPolyhedronIt.second);

        //The outer loop with the running i --> the planes
        thrust::transform(first, last, orthogonalProjectionPointsOnSegments.begin(), [](const auto &tuple) {
            //P' for plane i, sigma_pq[i] with fixed i, the nodes making up plane i
            const Array3 &orthogonalProjectionPointOnPlane = thrust::get<0>(tuple);
            const Array3 &segmentNormalOrientationsForPlane = thrust::get<1>(tuple);
            const Array3Triplet &face = thrust::get<2>(tuple);
            return computeOrthogonalProjectionPointsOnSegmentsForPlane(
                    orthogonalProjectionPointOnPlane, segmentNormalOrientationsForPlane, face);
        });

        return orthogonalProjectionPointsOnSegments;
    }

    std::vector<Array3> GravityModel::calculateSegmentDistances(
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment) {
        std::vector<Array3> segmentDistances{orthogonalProjectionPointsOnPlane.size()};
        //Iterating over planes (P'_i and P''_i are the parameters of the lambda)
        std::transform(orthogonalProjectionPointsOnPlane.cbegin(), orthogonalProjectionPointsOnPlane.cend(),
                       orthogonalProjectionPointsOnSegment.cbegin(), segmentDistances.begin(),
                       [](const Array3 projectionPointOnPlane, const Array3Triplet &projectionPointOnSegments) {
                           return computeSegmentDistancesForPlane(projectionPointOnPlane, projectionPointOnSegments);
                       });
        return segmentDistances;
    }

    std::vector<std::array<Distance, 3>> GravityModel::calculateDistances(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3Triplet> &segmentVectors,
            const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment) {
        std::vector<std::array<Distance, 3>> distances{segmentVectors.size()};

        //Zip the three required arguments together: G_ij for every segment, P'' for every segment
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(segmentVectors.begin(), orthogonalProjectionPointsOnSegment.begin(),
                               transformedPolyhedronIt.first);
        auto last = util::zip(segmentVectors.end(), orthogonalProjectionPointsOnSegment.end(),
                              transformedPolyhedronIt.second);

        thrust::transform(first, last, distances.begin(), [](const auto &tuple) {
            const Array3Triplet &segmentVectorsForPlane = thrust::get<0>(tuple);
            const Array3Triplet &orthogonalProjectionPointsOnSegmentForPlane = thrust::get<1>(tuple);
            const Array3Triplet &face = thrust::get<2>(tuple);
            return computeDistancesForPlane(segmentVectorsForPlane, orthogonalProjectionPointsOnSegmentForPlane, face);
        });
        return distances;
    }

    std::vector<std::array<TranscendentalExpression, 3>> GravityModel::calculateTranscendentalExpressions(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<std::array<Distance, 3>> &distances,
            const std::vector<double> &planeDistances,
            const std::vector<Array3> &segmentDistances,
            const std::vector<Array3> &segmentNormalOrientation,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane) {
        std::vector<std::array<TranscendentalExpression, 3>> transcendentalExpressions{distances.size()};

        //Zip iterator consisting of  3D and 1D distances l1/l2 and s1/2 | h_p | h_pq | sigma_pq | P'_p | faces
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(distances.begin(), planeDistances.begin(), segmentDistances.begin(),
                               segmentNormalOrientation.begin(), orthogonalProjectionPointsOnPlane.begin(),
                               transformedPolyhedronIt.first);
        auto last = util::zip(distances.end(), planeDistances.end(), segmentDistances.end(),
                              segmentNormalOrientation.end(), orthogonalProjectionPointsOnPlane.end(),
                              transformedPolyhedronIt.second);

        thrust::transform(first, last, transcendentalExpressions.begin(), [](const auto &tuple) {
            const std::array<Distance, 3> &distancesForPlane = thrust::get<0>(tuple);
            const double planeDistance = thrust::get<1>(tuple);
            const Array3 &segmentDistancesForPlane = thrust::get<2>(tuple);
            const Array3 &segmentNormalOrientationsForPlane = thrust::get<3>(tuple);
            const Array3 &projectionPointOnPlane = thrust::get<4>(tuple);
            const Array3Triplet &face = thrust::get<5>(tuple);

            return computeTranscendentalExpressionsForPlane(distancesForPlane, planeDistance, segmentDistancesForPlane,
                                                            segmentNormalOrientationsForPlane, projectionPointOnPlane,
                                                            face);
        });
        return transcendentalExpressions;
    }

    std::vector<std::pair<double, Array3>> GravityModel::calculateSingularityTerms(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3Triplet> &segmentVectors,
            const std::vector<Array3> &segmentNormalOrientation,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<double> &planeDistances,
            const std::vector<double> &planeNormalOrientations,
            const std::vector<Array3> &planeUnitNormals) {
        //The result
        std::vector<std::pair<double, Array3>> singularities{planeDistances.size()};

        //Zip iterator consisting of G_ij vectors | sigma_pq | faces | P' | h_p | sigma_p | N_i
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(segmentVectors.begin(), segmentNormalOrientation.begin(),
                               orthogonalProjectionPointsOnPlane.begin(), planeUnitNormals.begin(),
                               planeDistances.begin(), planeNormalOrientations.begin(), transformedPolyhedronIt.first);
        auto last = util::zip(segmentVectors.end(), segmentNormalOrientation.end(),
                              orthogonalProjectionPointsOnPlane.end(), planeUnitNormals.end(),
                              planeDistances.end(), planeNormalOrientations.end(), transformedPolyhedronIt.second);

        thrust::transform(first, last, singularities.begin(), [&](const auto &tuple) {
            const Array3Triplet &segmentVectorsForPlane = thrust::get<0>(tuple);
            const Array3 segmentNormalOrientationForPlane = thrust::get<1>(tuple);
            const Array3 &orthogonalProjectionPointOnPlane = thrust::get<2>(tuple);
            const Array3 &planeUnitNormal = thrust::get<3>(tuple);
            const double planeDistance = thrust::get<4>(tuple);
            const double planeNormalOrientation = thrust::get<5>(tuple);
            const Array3Triplet &face = thrust::get<6>(tuple);

            return computeSingularityTermsForPlane(segmentVectorsForPlane, segmentNormalOrientationForPlane,
                                                   orthogonalProjectionPointOnPlane, planeUnitNormal, planeDistance,
                                                   planeNormalOrientation, face);
        });
        return singularities;
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

        //Calculate the sign of the projections points x, y, z coordinates and apply it
        //if -D/A > 0 --> D/A < 0 --> everything is fine, no change
        //if -D/A < 0 --> D/A > 0 --> change sign if Ni is positive, else no change
        Array3 intersections = {hessianPlane.a == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.a,
                                hessianPlane.b == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.b,
                                hessianPlane.c == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.c};

        for (unsigned int index = 0; index < 3; ++index) {
            if (intersections[index] < 0) {
                orthogonalProjectionPoint[index] = std::abs(orthogonalProjectionPoint[index]);
            } else {
                if (planeUnitNormal[index] > 0) {
                    orthogonalProjectionPoint[index] = -1.0 * orthogonalProjectionPoint[index];
                }
                orthogonalProjectionPoint[index] = orthogonalProjectionPoint[index];
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
                       [&projectionPointOnPlane](const Array3 &segmentNormalOrientation, const Array3 &vertex) {
                           using namespace util;
                           return sgn((dot(segmentNormalOrientation, projectionPointOnPlane - vertex)), util::EPSILON) *
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

                    //Calculate the 3D distances between P (0, 0, 0) and the segment endpoints face[j] and face[(j + 1) % 3])
                    distance.l1 = euclideanNorm(face[j]);
                    distance.l2 = euclideanNorm(face[(j + 1) % 3]);
                    //Calculate the 1D distances between P'' (every segment has its own) and the segment endpoints
                    // face[j] and face[(j + 1) % 3])
                    distance.s1 = euclideanNorm(orthogonalProjectionPointsOnSegment - face[j]);
                    distance.s2 = euclideanNorm(orthogonalProjectionPointsOnSegment - face[(j + 1) % 3]);

                    //Change the sign depending on certain conditions
                    //1D and 3D distance are small
                    if (std::abs(distance.s1 - distance.l1) < 1e-10 &&
                        std::abs(distance.s2 - distance.l2) < 1e-10) {
                        if (distance.s2 < distance.s1) {
                            distance.s1 *= -1.0;
                            distance.s2 *= -1.0;
                            distance.l1 *= -1.0;
                            distance.l2 *= -1.0;
                            return distance;
                        } else if (distance.s2 == distance.s1) {
                            distance.s1 *= -1.0;
                            distance.l1 *= -1.0;
                            return distance;
                        }
                    } else {
                        //condition: P'' lies on the segment described by G_ij
                        const double norm = euclideanNorm(segmentVector);
                        if (distance.s1 < norm && distance.s2 < norm) {
                            distance.s1 *= -1.0;
                            return distance;
                        } else if (distance.s2 < distance.s1) {
                            distance.s1 *= -1.0;
                            distance.s2 *= -1.0;
                            return distance;
                        }
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
