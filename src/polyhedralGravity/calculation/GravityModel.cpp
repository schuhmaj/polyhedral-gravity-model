#include "GravityModel.h"

namespace polyhedralGravity {

    void GravityModel::calculate() {
        using namespace util;
        auto gijVectors = calculateGij();
        auto planeUnitNormals = calculatePlaneUnitNormals(gijVectors);
        auto segmentUnitNormals = calculateSegmentUnitNormals(gijVectors, planeUnitNormals);

        auto planeNormalOrientation = calculatePlaneNormalOrientations(planeUnitNormals);

        auto hessianPlanes = calculateFacesToHessianPlanes();

        auto planeDistances = calculatePlaneDistances(hessianPlanes);

        auto orthogonalProjectionOnPlane =
                calculateOrthogonalProjectionPointsOnPlane(hessianPlanes, planeUnitNormals, planeDistances);

        auto segmentNormalOrientation =
                calculateSegmentNormalOrientations(segmentUnitNormals, orthogonalProjectionOnPlane);

        auto orthogonalProjectionOnSegment =
                calculateOrthogonalProjectionPointsOnSegments(orthogonalProjectionOnPlane, segmentNormalOrientation);

        auto segmentDistances = calculateSegmentDistances(orthogonalProjectionOnPlane, orthogonalProjectionOnSegment);

        auto distances = calculateDistances(gijVectors, orthogonalProjectionOnSegment);

        auto transcendentalExpressions =
                calculateTranscendentalExpressions(distances, planeDistances, segmentDistances,
                                                   segmentNormalOrientation, orthogonalProjectionOnPlane);

        auto singularities =
                calculateSingularityTerms(gijVectors, segmentNormalOrientation, orthogonalProjectionOnPlane,
                                          planeDistances, planeNormalOrientation, planeUnitNormals);

        /*
         * Calculate V
         */

        auto firstV = thrust::make_zip_iterator(thrust::make_tuple(planeNormalOrientation.begin(),
                                                                   planeDistances.begin(),
                                                                   segmentNormalOrientation.begin(),
                                                                   segmentDistances.begin(),
                                                                   transcendentalExpressions.begin(),
                                                                   singularities.begin()));

        auto lastV = thrust::make_zip_iterator(thrust::make_tuple(planeNormalOrientation.end(),
                                                                  planeDistances.end(),
                                                                  segmentNormalOrientation.end(),
                                                                  segmentDistances.end(),
                                                                  transcendentalExpressions.end(),
                                                                  singularities.end()));

        double V = std::accumulate(firstV, lastV, 0.0, [](double acc, const auto &tuple) {
            const double sigma_p = thrust::get<0>(tuple);
            const double h_p = thrust::get<1>(tuple);
            const auto &sigmaPQPerPlane = thrust::get<2>(tuple);
            const auto &segmentDistancePerPlane = thrust::get<3>(tuple);
            const auto &transcendentalExpressionsPerPlane = thrust::get<4>(tuple);
            const std::pair<double, std::array<double, 3>> &singularitiesPerPlane = thrust::get<5>(tuple);

            auto sum1Start = thrust::make_zip_iterator(thrust::make_tuple(sigmaPQPerPlane.begin(),
                                                                          segmentDistancePerPlane.begin(),
                                                                          transcendentalExpressionsPerPlane.begin()));

            auto sum1End = thrust::make_zip_iterator(thrust::make_tuple(sigmaPQPerPlane.end(),
                                                                        segmentDistancePerPlane.end(),
                                                                        transcendentalExpressionsPerPlane.end()));

            const double sum1 = std::accumulate(sum1Start, sum1End, 0.0, [](double acc, const auto &tuple) {
                const double &sigma_pq = thrust::get<0>(tuple);
                const double &h_pq = thrust::get<1>(tuple);
                const TranscendentalExpression &transcendentalExpressions = thrust::get<2>(tuple);
                return acc + sigma_pq * h_pq * transcendentalExpressions.ln;
            });

            auto sum2Start = thrust::make_zip_iterator(thrust::make_tuple(sigmaPQPerPlane.begin(),
                                                                          transcendentalExpressionsPerPlane.begin()));

            auto sum2End = thrust::make_zip_iterator(thrust::make_tuple(sigmaPQPerPlane.end(),
                                                                        transcendentalExpressionsPerPlane.end()));

            const double sum2 = std::accumulate(sum2Start, sum2End, 0.0, [](double acc, const auto &tuple) {
                const double &sigma_pq = thrust::get<0>(tuple);
                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(tuple);
                return acc + sigma_pq * transcendentalExpressions.an;
            });


            return acc + sigma_p * h_p * (sum1 + h_p * sum2 + singularitiesPerPlane.first);
        });

        V = (V * util::GRAVITATIONAL_CONSTANT * _density) / 2.0;
        SPDLOG_INFO("V= {}", V);


        /*
         * Calculate Vx, Vy, Vz
         * TODO This code cries because it is a code clone and very unhappy about this circumstance
         */

        auto firstV2 = thrust::make_zip_iterator(thrust::make_tuple(planeNormalOrientation.begin(),
                                                                    planeDistances.begin(),
                                                                    segmentNormalOrientation.begin(),
                                                                    segmentDistances.begin(),
                                                                    transcendentalExpressions.begin(),
                                                                    singularities.begin(),
                                                                    planeUnitNormals.begin()));

        auto lastV2 = thrust::make_zip_iterator(thrust::make_tuple(planeNormalOrientation.end(),
                                                                   planeDistances.end(),
                                                                   segmentNormalOrientation.end(),
                                                                   segmentDistances.end(),
                                                                   transcendentalExpressions.end(),
                                                                   singularities.end(),
                                                                   planeUnitNormals.end()));

        Array3 V2 = std::accumulate(
                firstV2, lastV2, Array3{0.0, 0.0, 0.0},
                [](const Array3 &acc, const auto &tuple) {
                    using namespace util;
                    const double sigma_p = thrust::get<0>(tuple);
                    const double h_p = thrust::get<1>(tuple);
                    const auto &sigmaPQPerPlane = thrust::get<2>(tuple);
                    const auto &segmentDistancePerPlane = thrust::get<3>(tuple);
                    const auto &transcendentalExpressionsPerPlane = thrust::get<4>(tuple);
                    const std::pair<double, std::array<double, 3>> &singularitiesPerPlane = thrust::get<5>(
                            tuple);
                    const Array3 &Np = thrust::get<6>(tuple);

                    auto sum1Start = thrust::make_zip_iterator(thrust::make_tuple(
                            sigmaPQPerPlane.begin(),
                            segmentDistancePerPlane.begin(),
                            transcendentalExpressionsPerPlane.begin()));

                    auto sum1End = thrust::make_zip_iterator(thrust::make_tuple(
                            sigmaPQPerPlane.end(),
                            segmentDistancePerPlane.end(),
                            transcendentalExpressionsPerPlane.end()));

                    const double sum1 = std::accumulate(
                            sum1Start, sum1End, 0.0,
                            [](double acc, const auto &tuple) {
                                const double &sigma_pq = thrust::get<0>(tuple);
                                const double &h_pq = thrust::get<1>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<2>(tuple);
                                return acc + sigma_pq * h_pq * transcendentalExpressions.ln;
                            });

                    auto sum2Start = thrust::make_zip_iterator(thrust::make_tuple(
                            sigmaPQPerPlane.begin(),
                            transcendentalExpressionsPerPlane.begin()));

                    auto sum2End = thrust::make_zip_iterator(thrust::make_tuple(
                            sigmaPQPerPlane.end(),
                            transcendentalExpressionsPerPlane.end()));

                    const double sum2 = std::accumulate(
                            sum2Start, sum2End, 0.0,
                            [](double acc, const auto &tuple) {
                                const double &sigma_pq = thrust::get<0>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(tuple);
                                return acc + sigma_pq * transcendentalExpressions.an;
                            });

                    return acc + (Np * (sum1 + h_p * sum2 + singularitiesPerPlane.first));
                });

        V2 = abs(V2 * (util::GRAVITATIONAL_CONSTANT * _density));

        SPDLOG_INFO("Vx= {}", V2[0]);
        SPDLOG_INFO("Vy= {}", V2[1]);
        SPDLOG_INFO("Vz= {}", V2[2]);

        /*
         * Calculation of Vxx, Vyy, Vzz, Vxy, Vxz, Vyz
         */

        auto firstV3 = thrust::make_zip_iterator(thrust::make_tuple(planeNormalOrientation.begin(),
                                                                    planeDistances.begin(),
                                                                    segmentNormalOrientation.begin(),
                                                                    segmentDistances.begin(),
                                                                    transcendentalExpressions.begin(),
                                                                    singularities.begin(),
                                                                    planeUnitNormals.begin(),
                                                                    segmentUnitNormals.begin()));

        auto lastV3 = thrust::make_zip_iterator(thrust::make_tuple(planeNormalOrientation.end(),
                                                                   planeDistances.end(),
                                                                   segmentNormalOrientation.end(),
                                                                   segmentDistances.end(),
                                                                   transcendentalExpressions.end(),
                                                                   singularities.end(),
                                                                   planeUnitNormals.end(),
                                                                   segmentUnitNormals.end()));

        std::array<double, 6> V3 = std::accumulate(
                firstV3, lastV3, std::array<double, 6>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                [](const std::array<double, 6> &acc, const auto &tuple) {
                    using namespace util;
                    const double sigma_p = thrust::get<0>(tuple);
                    const double h_p = thrust::get<1>(tuple);
                    const auto &sigmaPQPerPlane = thrust::get<2>(tuple);
                    const auto &segmentDistancePerPlane = thrust::get<3>(tuple);
                    const auto &transcendentalExpressionsPerPlane = thrust::get<4>(tuple);
                    const std::pair<double, std::array<double, 3>> &singularitiesPerPlane = thrust::get<5>(
                            tuple);
                    const Array3 &Np = thrust::get<6>(tuple);
                    const auto &npqPerPlane = thrust::get<7>(tuple);


                    auto sum1Start = thrust::make_zip_iterator(thrust::make_tuple(
                            npqPerPlane.begin(),
                            transcendentalExpressionsPerPlane.begin()));

                    auto sum1End = thrust::make_zip_iterator(thrust::make_tuple(
                            npqPerPlane.end(),
                            transcendentalExpressionsPerPlane.end()));

                    const Array3 sum1 = std::accumulate(
                            sum1Start, sum1End, Array3{0.0, 0.0, 0.0},
                            [](const Array3 &acc, const auto &tuple) {
                                const Array3 &npq = thrust::get<0>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(tuple);
                                return acc + npq * transcendentalExpressions.ln;
                            });

                    auto sum2Start = thrust::make_zip_iterator(thrust::make_tuple(
                            sigmaPQPerPlane.begin(),
                            transcendentalExpressionsPerPlane.begin()));

                    auto sum2End = thrust::make_zip_iterator(thrust::make_tuple(
                            sigmaPQPerPlane.end(),
                            transcendentalExpressionsPerPlane.end()));

                    const double sum2 = std::accumulate(
                            sum2Start, sum2End, 0.0,
                            [](double acc, const auto &tuple) {
                                const double &sigma_pq = thrust::get<0>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(tuple);
                                return acc + sigma_pq * transcendentalExpressions.an;
                            });

                    const Array3 subSum = (sum1 + (Np * (sigma_p * sum2))) + singularitiesPerPlane.second;

                    const Array3 first = Np * subSum;
                    const Array3 reorderedNp = {Np[0], Np[0], Np[1]};
                    const Array3 reorderedSubSum = {subSum[1], subSum[2], subSum[2]};
                    const Array3 second = reorderedNp * reorderedSubSum;

                    return acc + concat(first, second);
                });

        V3 = V3 * (util::GRAVITATIONAL_CONSTANT * _density);

        SPDLOG_INFO("Vxx= {}", V3[0]);
        SPDLOG_INFO("Vyy= {}", V3[1]);
        SPDLOG_INFO("Vzz= {}", V3[2]);
        SPDLOG_INFO("Vxy= {}", V3[3]);
        SPDLOG_INFO("Vxz= {}", V3[4]);
        SPDLOG_INFO("Vyz= {}", V3[5]);


    }

    std::vector<Array3Triplet> GravityModel::calculateGij() {
        std::vector<Array3Triplet> g{_polyhedron.countFaces()};
        std::transform(_polyhedron.getFaces().cbegin(), _polyhedron.getFaces().cend(), g.begin(),
                       [&](const auto &face) -> Array3Triplet {
                           using util::operator-;
                           const auto &node0 = _polyhedron.getNode(face[0]);
                           const auto &node1 = _polyhedron.getNode(face[1]);
                           const auto &node2 = _polyhedron.getNode(face[2]);
                           return {node1 - node0, node2 - node1, node0 - node2};
                       });
        return g;
    }

    std::vector<Array3> GravityModel::calculatePlaneUnitNormals(const std::vector<Array3Triplet> &segmentVectors) {
        std::vector<Array3> planeUnitNormals{segmentVectors.size()};
        //Calculate N_i as (G_i1 * G_i2) / |G_i1 * G_i2| with * being the cross product
        std::transform(segmentVectors.cbegin(), segmentVectors.cend(), planeUnitNormals.begin(),
                       [](const auto &gi) -> Array3 {
                           using namespace util;
                           const Array3 crossProduct = cross(gi[0], gi[1]);
                           const double norm = euclideanNorm(crossProduct);
                           return crossProduct / norm;
                       });
        return planeUnitNormals;
    }

    std::vector<Array3Triplet>
    GravityModel::calculateSegmentUnitNormals(const std::vector<Array3Triplet> &segmentVectors,
                                              const std::vector<Array3> &planeUnitNormals) {
        std::vector<Array3Triplet> segmentUnitNormals{segmentVectors.size()};
        //Calculate n_ij as (G_ij * N_i) / |G_ig * N_i| with * being the cross product
        //Outer "loop" over G_i (running i) and N_i calculating n_i
        std::transform(segmentVectors.cbegin(), segmentVectors.cend(), planeUnitNormals.cbegin(),
                       segmentUnitNormals.begin(),
                       [](const Array3Triplet &gi, const Array3 &Ni) {
                           Array3Triplet ni{};
                           //Inner "loop" over G_ij (fixed i, running j) with parameter N_i calculating n_ij
                           std::transform(gi.cbegin(), gi.end(), ni.begin(),
                                          [&Ni](const auto &gij) -> Array3 {
                                              using namespace util;
                                              const Array3 crossProduct = cross(gij, Ni);
                                              const double norm = euclideanNorm(crossProduct);
                                              return crossProduct / norm;
                                          });
                           return ni;
                       });
        return segmentUnitNormals;
    }

    std::vector<double>
    GravityModel::calculatePlaneNormalOrientations(const std::vector<Array3> &planeUnitNormals) {
        std::vector<double> planeNormalOrientations(planeUnitNormals.size(), 0.0);
        //Calculate N_i * -G_i1 where * is the dot product and then use the inverted sgn
        std::transform(planeUnitNormals.cbegin(), planeUnitNormals.cend(), _polyhedron.getFaces().begin(),
                       planeNormalOrientations.begin(),
                       [&](const Array3 &ni, const std::array<size_t, 3> &gi) {
                           using namespace util;
                           //The first vertices' coordinates of the given face consisting of G_i's
                           const auto &Gi1 = _polyhedron.getNode(gi[0]);
                           //We abstain on the double multiplication with -1 in the line above and beyond since two
                           //times multiplying with -1 equals no change
                           return sgn(dot(ni, Gi1), util::EPSILON);
                       });
        return planeNormalOrientations;
    }


    std::vector<HessianPlane> GravityModel::calculateFacesToHessianPlanes(const Array3 &p) {
        std::vector<HessianPlane> hessianPlanes{_polyhedron.countFaces()};
        //Calculate for each face/ plane/ triangle (here) the Hessian Plane
        std::transform(_polyhedron.getFaces().cbegin(), _polyhedron.getFaces().cend(), hessianPlanes.begin(),
                       [&](const auto &face) -> HessianPlane {
                           using namespace util;
                           const auto &node0 = _polyhedron.getNode(face[0]);
                           const auto &node1 = _polyhedron.getNode(face[1]);
                           const auto &node2 = _polyhedron.getNode(face[2]);
                           //The three vertices put up the plane, p is the origin of the reference system default 0,0,0
                           return computeHessianPlane(node0, node1, node2, p);
                       });
        return hessianPlanes;
    }

    HessianPlane GravityModel::computeHessianPlane(const Array3 &p, const Array3 &q,
                                                   const Array3 &r, const Array3 &origin) {
        using namespace util;
        const auto crossProduct = cross(p - q, p - r);
        const auto res = (origin - p) * crossProduct;
        const double d = res[0] + res[1] + res[2];

        return {crossProduct[0], crossProduct[1], crossProduct[2], d};
    }

    std::vector<double> GravityModel::calculatePlaneDistances(const std::vector<HessianPlane> &plane) {
        std::vector<double> planeDistances(plane.size(), 0.0);
        //For each plane calculate h_p as D/sqrt(A^2 + B^2 + C^2)
        std::transform(plane.cbegin(), plane.cend(), planeDistances.begin(),
                       [](const HessianPlane &plane) -> double {
                           return std::abs(
                                   plane.d / std::sqrt(plane.a * plane.a + plane.b * plane.b + plane.c * plane.c));
                       });
        return planeDistances;
    }

    std::vector<Array3> GravityModel::calculateOrthogonalProjectionPointsOnPlane(
            const std::vector<HessianPlane> &hessianPlanes,
            const std::vector<Array3> &planeUnitNormals,
            const std::vector<double> &planeDistances) {
        std::vector<Array3> orthogonalProjectionPointsOfP{planeUnitNormals.size()};

        //Zip the three required arguments together: Plane normal N_i, Plane Distance h_i and the Hessian Form
        auto first = thrust::make_zip_iterator(thrust::make_tuple(planeUnitNormals.begin(),
                                                                  planeDistances.begin(),
                                                                  hessianPlanes.begin()));

        auto last = thrust::make_zip_iterator(thrust::make_tuple(planeUnitNormals.end(),
                                                                 planeDistances.end(),
                                                                 hessianPlanes.end()));

        thrust::transform(first, last, orthogonalProjectionPointsOfP.begin(), [](const auto &tuple) {
            using namespace util;
            const Array3 &Ni = thrust::get<0>(tuple);
            const double hi = thrust::get<1>(tuple);
            const HessianPlane &plane = thrust::get<2>(tuple);

            //Calculate the projection point by (22) P'_ = N_i / norm(N_i) * h_i
            // norm(N_i) is always 1 since N_i is a "normed" vector --> we do not need this division
            Array3 orthogonalProjectionPoint = Ni * hi;

            //Calculate the sign of the projections points x, y, z coordinates and apply it
            //if -D/A > 0 --> D/A < 0 --> everything is fine, no change
            //if -D/A < 0 --> D/A > 0 --> change sign if Ni is positive, else no change
            Array3 intersections = {plane.a == 0.0 ? 0.0 : plane.d / plane.a,
                                    plane.b == 0.0 ? 0.0 : plane.d / plane.b,
                                    plane.c == 0.0 ? 0.0 : plane.d / plane.c};

            for (unsigned int index = 0; index < 3; ++index) {
                if (intersections[index] < 0) {
                    orthogonalProjectionPoint[index] = std::abs(orthogonalProjectionPoint[index]);
                } else {
                    if (Ni[index] > 0) {
                        orthogonalProjectionPoint[index] = -1.0 * orthogonalProjectionPoint[index];
                    }
                    orthogonalProjectionPoint[index] = orthogonalProjectionPoint[index];
                }
            }
            return orthogonalProjectionPoint;
        });
        return orthogonalProjectionPointsOfP;
    }

    std::vector<Array3> GravityModel::calculateSegmentNormalOrientations(
            const std::vector<Array3Triplet> &segmentUnitNormals,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane) {
        std::vector<Array3> segmentNormalOrientations{segmentUnitNormals.size()};

        std::vector<Array3Triplet> x{segmentUnitNormals.size()};

        //First part of equation (23):
        //Calculate x_P' - x_ij^1 (x_P' is the projectionPoint and x_ij^1 is the first vertices of one segment,
        //i.e. the coordinates of the training-planes' nodes
        //The result is saved in x
        std::transform(orthogonalProjectionPointsOnPlane.cbegin(), orthogonalProjectionPointsOnPlane.cend(),
                       _polyhedron.getFaces().cbegin(), x.begin(),
                       [&](const Array3 &projectionPoint, const std::array<size_t, 3> &face)
                               -> Array3Triplet {
                           using util::operator-;
                           const auto &node0 = _polyhedron.getNode(face[0]);
                           const auto &node1 = _polyhedron.getNode(face[1]);
                           const auto &node2 = _polyhedron.getNode(face[2]);
                           return {projectionPoint - node0, projectionPoint - node1, projectionPoint - node2};
                       });
        //The second part of equation (23)
        //Calculate n_ij * x_ij with * being the dot product and use the inverted sgn to determine the value of sigma_pq
        //running over n_i and x_i (running i)
        std::transform(segmentUnitNormals.cbegin(), segmentUnitNormals.cend(), x.cbegin(),
                       segmentNormalOrientations.begin(),
                       [](const Array3Triplet &ni, const Array3Triplet &xi) {
                           //running over n_ij and x_ij (fixed i, running j)
                           std::array<double, 3> sigmaPQ{};
                           std::transform(ni.cbegin(), ni.cend(), xi.cbegin(), sigmaPQ.begin(),
                                          [](const Array3 &nij, const Array3 &xij) {
                                              using namespace util;
                                              return sgn((dot(nij, xij)), util::EPSILON) * -1.0;
                                          });
                           return sigmaPQ;
                       });
        return segmentNormalOrientations;
    }

    std::vector<Array3Triplet> GravityModel::calculateOrthogonalProjectionPointsOnSegments(
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<Array3> &segmentNormalOrientation) {
        std::vector<Array3Triplet> orthogonalProjectionPointsOfPPrime{orthogonalProjectionPointsOnPlane.size()};

        //Zip the three required arguments together: P' for every plane, sigma_pq for every segment, the faces
        auto first = thrust::make_zip_iterator(thrust::make_tuple(orthogonalProjectionPointsOnPlane.begin(),
                                                                  segmentNormalOrientation.begin(),
                                                                  _polyhedron.getFaces().begin()));
        auto last = thrust::make_zip_iterator(thrust::make_tuple(orthogonalProjectionPointsOnPlane.end(),
                                                                 segmentNormalOrientation.end(),
                                                                 _polyhedron.getFaces().end()));

        //The outer loop with the running i --> the planes
        thrust::transform(first, last, orthogonalProjectionPointsOfPPrime.begin(), [&](const auto &tuple) {
            //P' for plane i, sigma_pq[i] with fixed i, the nodes making up plane i
            const auto &pPrime = thrust::get<0>(tuple);
            const auto &sigmaP = thrust::get<1>(tuple);
            const auto &face = thrust::get<2>(tuple);

            auto counterJ = thrust::counting_iterator<unsigned int>(0);
            Array3Triplet pDoublePrime{};

            //The inner loop with fixed i, running over the j --> the segments of a plane
            thrust::transform(sigmaP.begin(), sigmaP.end(), counterJ, pDoublePrime.begin(),
                              [&](const auto &sigmaPQ, unsigned int j) {
                                  //We actually only accept +0.0 or -0.0, so the equal comparison is ok
                                  if (sigmaPQ == 0.0) {
                                      //Geometrically trivial case, in neither of the half space --> already on segment
                                      return pPrime;
                                  } else {
                                      //In one of the half space, calculate the projection point P'' for the segment
                                      //with the endpoints v1 and v2
                                      const auto &v1 = _polyhedron.getNode(face[j]);
                                      const auto &v2 = _polyhedron.getNode(face[(j + 1) % 3]);
                                      return calculateOrthogonalProjectionOnSegment(v1, v2, pPrime);
                                  }
                              });
            return pDoublePrime;
        });

        return orthogonalProjectionPointsOfPPrime;
    }

    Array3 GravityModel::calculateOrthogonalProjectionOnSegment(const Array3 &v1, const Array3 &v2,
                                                                const Array3 &pPrime) {
        using namespace util;
        //Preparing our the planes/ equations in matrix form
        const Array3 matrixRow1 = v2 - v1;
        const Array3 matrixRow2 = cross(v1 - pPrime, matrixRow1);
        const Array3 matrixRow3 = cross(matrixRow2, matrixRow1);
        const Array3 d = {dot(matrixRow1, pPrime), dot(matrixRow2, pPrime), dot(matrixRow3, v1)};
        Matrix<double, 3, 3> columnMatrix = transpose(Matrix<double, 3, 3>{matrixRow1, matrixRow2, matrixRow3});
        //Calculation and solving the equations of above
        const double determinant = det(columnMatrix);
        return Array3{
                det(Matrix<double, 3, 3>{d, columnMatrix[1], columnMatrix[2]}),
                det(Matrix<double, 3, 3>{columnMatrix[0], d, columnMatrix[2]}),
                det(Matrix<double, 3, 3>{columnMatrix[0], columnMatrix[1], d})
        } / determinant;
    }

    std::vector<Array3> GravityModel::calculateSegmentDistances(
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment) {
        std::vector<Array3> segmentDistances{orthogonalProjectionPointsOnPlane.size()};
        //The outer loop with the running i --> iterating over planes (P'_i and P''_i are the parameters of the lambda)
        std::transform(orthogonalProjectionPointsOnPlane.cbegin(), orthogonalProjectionPointsOnPlane.cend(),
                       orthogonalProjectionPointsOnSegment.cbegin(), segmentDistances.begin(),
                       [](const Array3 pPrime, const Array3Triplet &pDoublePrimes) {
                           std::array<double, 3> hp{};
                           //The inner loop with the running j --> iterating over the segments
                           //Using the values P'_i and P''_ij for the calculation of the distance
                           std::transform(pDoublePrimes.cbegin(), pDoublePrimes.cend(), hp.begin(),
                                          [&pPrime](const Array3 &pDoublePrime) {
                                              using namespace util;
                                              return euclideanNorm(pDoublePrime - pPrime);
                                          });
                           return hp;
                       });
        return segmentDistances;
    }

    std::vector<std::array<Distance, 3>> GravityModel::calculateDistances(
            const std::vector<Array3Triplet> &gij,
            const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment) {
        std::vector<std::array<Distance, 3>> distances{gij.size()};

        //Zip the three required arguments together: G_ij for every segment, P'' for every segment
        auto first = thrust::make_zip_iterator(thrust::make_tuple(gij.begin(),
                                                                  orthogonalProjectionPointsOnSegment.begin(),
                                                                  _polyhedron.getFaces().begin()));
        auto last = thrust::make_zip_iterator(thrust::make_tuple(gij.end(),
                                                                 orthogonalProjectionPointsOnSegment.end(),
                                                                 _polyhedron.getFaces().end()));

        thrust::transform(first, last, distances.begin(), [&](const auto &tuple) {
            const auto &gi = thrust::get<0>(tuple);
            const auto &pDoublePrimePerPlane = thrust::get<1>(tuple);
            const auto &face = thrust::get<2>(tuple);

            std::array<Distance, 3> distancesArray{};
            auto counterJ = thrust::counting_iterator<unsigned int>(0);

            auto first = thrust::make_zip_iterator(thrust::make_tuple(gi.begin(),
                                                                      pDoublePrimePerPlane.begin(),
                                                                      counterJ));
            auto last = thrust::make_zip_iterator(thrust::make_tuple(gi.end(),
                                                                     pDoublePrimePerPlane.end(),
                                                                     counterJ + 3));

            thrust::transform(first, last, distancesArray.begin(), [&](const auto &tuple) {
                using namespace util;
                Distance distance{};
                const Array3 &gijVector = thrust::get<0>(tuple);
                const Array3 &pDoublePrime = thrust::get<1>(tuple);
                const unsigned int j = thrust::get<2>(tuple);

                //Get the vertices (endpoints) of this segment
                const auto &v1 = _polyhedron.getNode(face[j]);
                const auto &v2 = _polyhedron.getNode(face[(j + 1) % 3]);

                //Calculate the 3D distances between P (0, 0, 0) and the segment endpoints v1 and v2
                distance.l1 = euclideanNorm(v1);
                distance.l2 = euclideanNorm(v2);
                //Calculate the 1D distances between P'' (every segment has its own) and the segment endpoints v1 and v2
                distance.s1 = euclideanNorm(pDoublePrime - v1);
                distance.s2 = euclideanNorm(pDoublePrime - v2);

                //Change the sign depending on certain conditions
                //1D and 3D distance are small
                if (std::abs(distance.s1 - distance.l1) < 1e-10 && std::abs(distance.s2 - distance.l2) < 1e-10) {
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
                    const double norm = euclideanNorm(gijVector);
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
            return distancesArray;
        });
        return distances;
    }

    std::vector<std::array<TranscendentalExpression, 3>>
    GravityModel::calculateTranscendentalExpressions(const std::vector<std::array<Distance, 3>> &distances,
                                                     const std::vector<double> &planeDistances,
                                                     const std::vector<Array3> &segmentDistances,
                                                     const std::vector<Array3> &segmentNormalOrientation,
                                                     const std::vector<Array3> &orthogonalProjectionPointsOnPlane) {
        //The result of this functions
        std::vector<std::array<TranscendentalExpression, 3>> transcendentalExpressions{distances.size()};

        //Zip iterator consisting of  3D and 1D distances l1/l2 and s1/2 | h_p | h_pq | sigma_pq | P'_p | faces
        //TODO Add h_pq
        auto first = thrust::make_zip_iterator(thrust::make_tuple(distances.begin(),
                                                                  planeDistances.begin(),
                                                                  segmentDistances.begin(),
                                                                  segmentNormalOrientation.begin(),
                                                                  orthogonalProjectionPointsOnPlane.begin(),
                                                                  _polyhedron.getFaces().begin()));
        auto last = thrust::make_zip_iterator(thrust::make_tuple(distances.end(),
                                                                 planeDistances.end(),
                                                                 segmentDistances.end(),
                                                                 segmentNormalOrientation.end(),
                                                                 orthogonalProjectionPointsOnPlane.end(),
                                                                 _polyhedron.getFaces().end()));

        thrust::transform(first, last, transcendentalExpressions.begin(), [&](const auto &tuple) {
            const auto &distancesPerPlane = thrust::get<0>(tuple);
            const double hp = thrust::get<1>(tuple);
            const auto &segmentDistancesPerPlane = thrust::get<2>(tuple);
            const auto &segmentNormalOrientationPerPlane = thrust::get<3>(tuple);
            const Array3 &pPrime = thrust::get<4>(tuple);
            const auto &face = thrust::get<5>(tuple);

            auto counterJ = thrust::counting_iterator<unsigned int>(0);


            //Zip iterator consisting of 3D and 1D distances l1/l2 and s1/2 for this plane | sigma_pq for this plane
            auto first = thrust::make_zip_iterator(thrust::make_tuple(distancesPerPlane.begin(),
                                                                      segmentDistancesPerPlane.begin(),
                                                                      segmentNormalOrientationPerPlane.begin()));
            auto last = thrust::make_zip_iterator(thrust::make_tuple(distancesPerPlane.end(),
                                                                     segmentDistancesPerPlane.end(),
                                                                     segmentNormalOrientationPerPlane.end()));

            //Result for this plane
            std::array<TranscendentalExpression, 3> transcendentalExpressionsPerPlane{};

            thrust::transform(first, last, counterJ, transcendentalExpressionsPerPlane.begin(),
                              [&](const auto &tuple, const unsigned int j) {
                                  using namespace util;
                                  const Distance &distance = thrust::get<0>(tuple);
                                  const double hpq = thrust::get<1>(tuple);
                                  const double sigmaPQ = thrust::get<2>(tuple);

                                  //Result for this segment
                                  TranscendentalExpression transcendentalExpressionPerSegment{};

                                  //Vertices (endpoints) of this segment
                                  const Array3 &v1 = _polyhedron.getNode(face[(j + 1) % 3]);
                                  const Array3 &v2 = _polyhedron.getNode(face[j]);

                                  //Compute LN_pq according to (14)
                                  //If either sigmaPQ has no sign AND either of the distances of P' to the two
                                  //segment endpoints is zero OR the 1D and 3D distances are below some threshold
                                  //then LN_pq is zero, too
                                  if ((sigmaPQ == 0.0 &&
                                       (euclideanNorm(pPrime - v1) == 0.0 || euclideanNorm(pPrime - v2) == 0.0)) ||
                                      (distance.s1 + distance.s2 < 1e-10 && distance.l1 + distance.l2 < 1e-10)) {
                                      transcendentalExpressionPerSegment.ln = 0.0;
                                  } else {
                                      transcendentalExpressionPerSegment.ln =
                                              std::log((distance.s2 + distance.l2) / (distance.s1 + distance.l1));
                                  }

                                  //Compute AN_pq according to (15)
                                  //If h_p or h_pq is zero then AN_pq is zero, too
                                  if (hp == 0 || hpq == 0) {
                                      transcendentalExpressionPerSegment.an = 0.0;
                                  } else {
                                      transcendentalExpressionPerSegment.an =
                                              std::atan((hp * distance.s2) / (hpq * distance.l2)) -
                                              std::atan((hp * distance.s1) / (hpq * distance.l1));
                                  }

                                  return transcendentalExpressionPerSegment;
                              });

            return transcendentalExpressionsPerPlane;

        });

        return transcendentalExpressions;
    }

    std::vector<std::pair<double, Array3>> GravityModel::calculateSingularityTerms(
            const std::vector<Array3Triplet> &gijVectors,
            const std::vector<Array3> &segmentNormalOrientation,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<double> &planeDistances,
            const std::vector<double> &planeNormalOrientation,
            const std::vector<Array3> &planeUnitNormals) {
        //The result
        std::vector<std::pair<double, Array3>> singularities{planeDistances.size()};

        //Zip iterator consisting of G_ij vectors | sigma_pq | faces | P' | h_p
        auto first = thrust::make_zip_iterator(thrust::make_tuple(gijVectors.begin(),
                                                                  segmentNormalOrientation.begin(),
                                                                  _polyhedron.getFaces().begin(),
                                                                  orthogonalProjectionPointsOnPlane.begin(),
                                                                  planeDistances.begin(),
                                                                  planeNormalOrientation.begin(),
                                                                  planeUnitNormals.begin()));
        auto last = thrust::make_zip_iterator(thrust::make_tuple(gijVectors.end(),
                                                                 segmentNormalOrientation.end(),
                                                                 _polyhedron.getFaces().end(),
                                                                 orthogonalProjectionPointsOnPlane.end(),
                                                                 planeDistances.end(),
                                                                 planeNormalOrientation.end(),
                                                                 planeUnitNormals.end()));

        thrust::transform(first, last, singularities.begin(), [&](const auto &tuple) {
            const auto &gijVectorsPerPlane = thrust::get<0>(tuple);
            const auto segmentNormalOrientationPerPlane = thrust::get<1>(tuple);
            const auto &face = thrust::get<2>(tuple);
            const Array3 &pPrime = thrust::get<3>(tuple);
            const double hp = thrust::get<4>(tuple);
            const double sigmaP = thrust::get<5>(tuple);
            const Array3 &Np = thrust::get<6>(tuple);

            //1. case: If all sigma_pq for a given plane p are 1.0 then P' lies inside the plane S_p
            if (std::all_of(segmentNormalOrientationPerPlane.cbegin(), segmentNormalOrientationPerPlane.cend(),
                            [](const double sigma) { return sigma == 1.0; })) {
                using namespace util;
                return std::make_pair(-1.0 * util::PI2 * hp, Np / euclideanNorm(Np) * -1.0 * util::PI2 * sigmaP);
            }
            //2. case: If sigma_pq == 0 AND norm(P' - v1) < norm(G_ij) && norm(P' - v2) < norm(G_ij) with G_ij
            // as the vector of v1 and v2
            // then P' is located on one line segment G_p of plane p, but not on any of its vertices
            auto counterJ = thrust::counting_iterator<unsigned int>(0);
            auto secondCaseBegin = thrust::make_zip_iterator(thrust::make_tuple(
                    gijVectorsPerPlane.begin(),
                    segmentNormalOrientationPerPlane.begin(),
                    counterJ));
            auto secondCaseEnd = thrust::make_zip_iterator(thrust::make_tuple(
                    gijVectorsPerPlane.end(),
                    segmentNormalOrientationPerPlane.end(),
                    counterJ + 3));
            if (std::any_of(secondCaseBegin, secondCaseEnd, [&](const auto &tuple) {
                using namespace util;
                const Array3 &gij = thrust::get<0>(tuple);
                const double sigmaPQ = thrust::get<1>(tuple);
                const unsigned int j = thrust::get<2>(tuple);

                if (sigmaPQ != 0.0) {
                    return false;
                }

                const Array3 &v1 = _polyhedron.getNode(face[(j + 1) % 3]);
                const Array3 &v2 = _polyhedron.getNode(face[j]);
                const double gijNorm = euclideanNorm(gij);
                return euclideanNorm(pPrime - v1) < gijNorm && euclideanNorm(pPrime - v2) < gijNorm;
            })) {
                using namespace util;
                return std::make_pair(-1.0 * util::PI * hp, Np / euclideanNorm(Np) * -1.0 * util::PI * sigmaP);
            }
            //3. case If sigma_pq == 0 AND norm(P' - v1) < 0 || norm(P' - v2) < 0
            // then P' is located at one of G_p's vertices
            auto counterJ3 = thrust::counting_iterator<int>(0);
            auto thirdCaseBegin = thrust::make_zip_iterator(thrust::make_tuple(
                    segmentNormalOrientationPerPlane.begin(),
                    counterJ3));
            auto thirdCaseEnd = thrust::make_zip_iterator(thrust::make_tuple(
                    segmentNormalOrientationPerPlane.end(),
                    counterJ3 + 3));
            double e1;
            double e2;
            unsigned int j;
            if (std::any_of(thirdCaseBegin, thirdCaseEnd, [&](const auto &tuple) {
                using namespace util;
                const double sigmaPQ = thrust::get<0>(tuple);
                j = thrust::get<1>(tuple);

                if (sigmaPQ != 0.0) {
                    return false;
                }

                const Array3 &v1 = _polyhedron.getNode(face[(j + 1) % 3]);
                const Array3 &v2 = _polyhedron.getNode(face[j]);
                e1 = euclideanNorm(pPrime - v1);
                e2 = euclideanNorm(pPrime - v2);
                return e1 == 0.0 || e2 == 0.0;
            })) {
                using namespace util;
                const Array3 &g1 = e1 == 0.0 ? gijVectorsPerPlane[j] : gijVectorsPerPlane[(j - 1 + 3) % 3];
                const Array3 &g2 = e1 == 0.0 ? gijVectorsPerPlane[(j + 1) % 3] : gijVectorsPerPlane[j];
                const double gdot = dot(g1 * -1.0, g2);
                const double theta =
                        gdot == 0.0 ? util::PI_2 : std::acos(gdot / (euclideanNorm(g1) * euclideanNorm(g2)));
                return std::make_pair(-1.0 * theta * hp, Np / euclideanNorm(Np) * -1.0 * theta * sigmaP);
            }

            //4. case Otherwise P' is located outside the plane S_p and then the singularity equals zero
            return std::make_pair(0.0, Array3{0, 0, 0});
        });

        return singularities;
    }


}
