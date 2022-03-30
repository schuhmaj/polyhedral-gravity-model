#include "GravityModel.h"

namespace polyhedralGravity {

    void GravityModel::calculate() {
        SPDLOG_INFO("Calculate...");
        auto x = calculateGij();
    }

    CartesianSegmentPropertyVector GravityModel::calculateGij() {
        CartesianSegmentPropertyVector g{_polyhedron.countFaces()};
        std::transform(_polyhedron.getFaces().cbegin(), _polyhedron.getFaces().cend(), g.begin(),
                       [&](const auto &face) -> CartesianSegmentProperty {
                           using util::operator-;
                           const auto &node0 = _polyhedron.getNode(face[0]);
                           const auto &node1 = _polyhedron.getNode(face[1]);
                           const auto &node2 = _polyhedron.getNode(face[2]);
                           return {node1 - node0, node2 - node1, node0 - node2};
                       });
        return g;
    }

    CartesianPlanePropertyVector GravityModel::calculatePlaneUnitNormals(const CartesianSegmentPropertyVector &g) {
        CartesianPlanePropertyVector planeUnitNormals{g.size()};
        //Calculate N_i as (G_i1 * G_i2) / |G_i1 * G_i2| with * being the cross product
        std::transform(g.cbegin(), g.cend(), planeUnitNormals.begin(), [](const auto &gi) -> Cartesian {
            using namespace util;
            const Cartesian crossProduct = cross(gi[0], gi[1]);
            const double norm = euclideanNorm(crossProduct);
            return crossProduct / norm;
        });
        return planeUnitNormals;
    }

    CartesianSegmentPropertyVector GravityModel::calculateSegmentUnitNormals(const CartesianSegmentPropertyVector &g,
                                                                             const CartesianPlanePropertyVector &planeUnitNormals) {
        CartesianSegmentPropertyVector segmentUnitNormals{g.size()};
        //Calculate n_ij as (G_ij * N_i) / |G_ig * N_i| with * being the cross product
        //Outer "loop" over G_i (running i) and N_i calculating n_i
        std::transform(g.cbegin(), g.cend(), planeUnitNormals.cbegin(), segmentUnitNormals.begin(),
                       [](const CartesianSegmentProperty &gi, const Cartesian &Ni) {
                           CartesianSegmentProperty ni{};
                           //Inner "loop" over G_ij (fixed i, running j) with parameter N_i calculating n_ij
                           std::transform(gi.cbegin(), gi.end(), ni.begin(),
                                          [&Ni](const auto &gij) -> Cartesian {
                                              using namespace util;
                                              const Cartesian crossProduct = cross(gij, Ni);
                                              const double norm = euclideanNorm(crossProduct);
                                              return crossProduct / norm;
                                          });
                           return ni;
                       });
        return segmentUnitNormals;
    }

    PlanePropertyVector
    GravityModel::calculatePlaneNormalOrientations(const CartesianPlanePropertyVector &planeUnitNormals) {
        PlanePropertyVector planeNormalOrientations(planeUnitNormals.size(), 0.0);
        //Calculate N_i * -G_i1 where * is the dot product and then use the inverted sgn
        std::transform(planeUnitNormals.cbegin(), planeUnitNormals.cend(), _polyhedron.getFaces().begin(),
                       planeNormalOrientations.begin(),
                       [&](const Cartesian &ni, const std::array<size_t, 3> &gi) {
                           using namespace util;
                           //The first vertices' coordinates of the given face consisting of G_i's
                           const auto &Gi1 = _polyhedron.getNode(gi[0]);
                           //We abstain on the double multiplication with -1 in the line above and beyond since two
                           //times multiplying with -1 equals no change
                           return sgn(dot(ni, Gi1), util::epsilon);
                       });
        return planeNormalOrientations;
    }


    std::vector<HessianPlane> GravityModel::calculateFacesToHessianPlanes(const Cartesian &p) {
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

    HessianPlane GravityModel::computeHessianPlane(const Cartesian &p, const Cartesian &q,
                                                   const Cartesian &r, const Cartesian &origin) {
        using namespace util;
        const auto crossProduct = cross(p - q, p - r);
        const auto res = (origin - p) * crossProduct;
        const double d = res[0] + res[1] + res[2];

        return {crossProduct[0], crossProduct[1], crossProduct[2], d};
    }

    PlanePropertyVector GravityModel::calculatePlaneDistances(const std::vector<HessianPlane> &plane) {
        PlanePropertyVector planeDistances(plane.size(), 0.0);
        //For each plane calculate h_p as D/sqrt(A^2 + B^2 + C^2)
        std::transform(plane.cbegin(), plane.cend(), planeDistances.begin(),
                       [](const HessianPlane &plane) -> double {
                           return std::abs(
                                   plane.d / std::sqrt(plane.a * plane.a + plane.b * plane.b + plane.c * plane.c));
                       });
        return planeDistances;
    }

    CartesianPlanePropertyVector GravityModel::calculateOrthogonalProjectionPointsOnPlane(
            const std::vector<HessianPlane> &hessianPlanes,
            const CartesianPlanePropertyVector &planeUnitNormals,
            const PlanePropertyVector &planeDistances) {
        CartesianPlanePropertyVector orthogonalProjectionPointsOfP{planeUnitNormals.size()};

        //Zip the three required arguments together: Plane normal N_i, Plane Distance h_i and the Hessian Form
        auto first = thrust::make_zip_iterator(thrust::make_tuple(planeUnitNormals.begin(),
                                                                  planeDistances.begin(),
                                                                  hessianPlanes.begin()));

        auto last = thrust::make_zip_iterator(thrust::make_tuple(planeUnitNormals.end(),
                                                                 planeDistances.end(),
                                                                 hessianPlanes.end()));

        thrust::transform(first, last, orthogonalProjectionPointsOfP.begin(), [](const auto &tuple) {
            using namespace util;
            const Cartesian &Ni = thrust::get<0>(tuple);
            const double hi = thrust::get<1>(tuple);
            const HessianPlane &plane = thrust::get<2>(tuple);

            //Calculate the projection point by (22) P'_ = N_i / norm(N_i) * h_i
            const Cartesian directionCosine = Ni / euclideanNorm(Ni);
            Cartesian orthogonalProjectionPoint = abs(directionCosine * hi);

            //Calculate the sign of the projections points x, y, z coordinates and apply it
            //if -D/A > 0 --> D/A < 0 --> everything is fine, no change
            //if -D/A < 0 --> D/A > 0 --> change sign if Ni is positive, else no change
            orthogonalProjectionPoint[0] *= plane.a == 0.0 ? 1.0 : plane.d / plane.a > 0.0 && Ni[0] < 0.0 ? -1.0 : 1.0;
            orthogonalProjectionPoint[1] *= plane.b == 0.0 ? 1.0 : plane.d / plane.b > 0.0 && Ni[1] < 0.0 ? -1.0 : 1.0;
            orthogonalProjectionPoint[2] *= plane.c == 0.0 ? 1.0 : plane.d / plane.c > 0.0 && Ni[2] < 0.0 ? -1.0 : 1.0;

            return orthogonalProjectionPoint;
        });
        return orthogonalProjectionPointsOfP;
    }

    SegmentPropertyVector GravityModel::calculateSegmentNormalOrientations(
            const CartesianSegmentPropertyVector &segmentUnitNormals,
            const CartesianPlanePropertyVector &orthogonalProjectionPointsOnPlane) {
        SegmentPropertyVector segmentNormalOrientations{segmentUnitNormals.size()};

        CartesianSegmentPropertyVector x{segmentUnitNormals.size()};

        //First part of equation (23):
        //Calculate x_P' - x_ij^1 (x_P' is the projectionPoint and x_ij^1 is the first vertices of one segment,
        //i.e. the coordinates of the training-planes' nodes
        //The result is saved in x
        std::transform(orthogonalProjectionPointsOnPlane.cbegin(), orthogonalProjectionPointsOnPlane.cend(),
                       _polyhedron.getFaces().cbegin(), x.begin(),
                       [&](const Cartesian &projectionPoint, const std::array<size_t, 3> &face)
                               -> CartesianSegmentProperty {
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
                       [](const CartesianSegmentProperty &ni, const CartesianSegmentProperty &xi) {
                           //running over n_ij and x_ij (fixed i, running j)
                           std::array<double, 3> sigmaPQ{};
                           std::transform(ni.cbegin(), ni.cend(), xi.cbegin(), sigmaPQ.begin(),
                                          [](const Cartesian &nij, const Cartesian &xij) {
                                              using namespace util;
                                              return sgn((dot(nij, xij)), util::epsilon) * -1.0;
                                          });
                           return sigmaPQ;
                       });
        return segmentNormalOrientations;
    }

    CartesianSegmentPropertyVector GravityModel::calculateOrthogonalProjectionPointsOnSegments(
            const CartesianPlanePropertyVector &orthogonalProjectionPointsOnPlane,
            const SegmentPropertyVector &segmentNormalOrientation) {
        CartesianSegmentPropertyVector orthogonalProjectionPointsOfPPrime{orthogonalProjectionPointsOnPlane.size()};

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
            CartesianSegmentProperty pDoublePrime{};

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

    Cartesian GravityModel::calculateOrthogonalProjectionOnSegment(const Cartesian &v1, const Cartesian &v2,
                                                                   const Cartesian &pPrime) {
        using namespace util;
        Cartesian pDoublePrime{};
        //Preparing our the planes/ equations in matrix form
        const Cartesian matrixRow1 = v2 - v1;
        const Cartesian matrixRow2 = cross(v1 - pPrime, matrixRow1);
        const Cartesian matrixRow3 = cross(matrixRow2, matrixRow1);
        const Cartesian d = {dot(matrixRow1, pPrime), dot(matrixRow2, pPrime), dot(matrixRow3, v1)};
        Matrix<double, 3, 3> columnMatrix = transpose(Matrix<double, 3, 3>{matrixRow1, matrixRow2, matrixRow3});
        //Calculation and solving the equations of above
        const double determinant = det(columnMatrix);
        pDoublePrime[0] = det(Matrix<double, 3, 3>{d, columnMatrix[1], columnMatrix[2]});
        pDoublePrime[1] = det(Matrix<double, 3, 3>{columnMatrix[0], d, columnMatrix[2]});
        pDoublePrime[2] = det(Matrix<double, 3, 3>{columnMatrix[0], columnMatrix[1], d});
        return pDoublePrime / determinant;
    }

    SegmentPropertyVector GravityModel::calculateSegmentDistances(
            const CartesianPlanePropertyVector &orthogonalProjectionPointsOnPlane,
            const CartesianSegmentPropertyVector &orthogonalProjectionPointsOnSegment) {
        SegmentPropertyVector segmentDistances{orthogonalProjectionPointsOnPlane.size()};
        //The outer loop with the running i --> iterating over planes (P'_i and P''_i are the parameters of the lambda)
        std::transform(orthogonalProjectionPointsOnPlane.cbegin(), orthogonalProjectionPointsOnPlane.cend(),
                       orthogonalProjectionPointsOnSegment.cbegin(), segmentDistances.begin(),
                       [](const Cartesian pPrime, const CartesianSegmentProperty &pDoublePrimes) {
                           std::array<double, 3> hp{};
                           //The inner loop with the running j --> iterating over the segments
                           //Using the values P'_i and P''_ij for the calculation of the distance
                           std::transform(pDoublePrimes.cbegin(), pDoublePrimes.cend(), hp.begin(),
                                          [&pPrime](const Cartesian &pDoublePrime) {
                                              using namespace util;
                                              return euclideanNorm(pDoublePrime - pPrime);
                                          });
                           return hp;
                       });
        return segmentDistances;
    }

    std::vector<std::array<Distance, 3>> GravityModel::calculateDistances(
            const CartesianSegmentPropertyVector &gij,
            const CartesianSegmentPropertyVector &orthogonalProjectionPointsOnSegment) {
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
                const Cartesian &gijVector = thrust::get<0>(tuple);
                const Cartesian &pDoublePrime = thrust::get<1>(tuple);
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
                                                     const PlanePropertyVector &planeDistances,
                                                     const SegmentPropertyVector &segmentDistances,
                                                     const SegmentPropertyVector &segmentNormalOrientation,
                                                     const CartesianPlanePropertyVector &orthogonalProjectionPointsOnPlane) {
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
            const Cartesian &pPrime = thrust::get<4>(tuple);
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
                                  const Cartesian &v1 = _polyhedron.getNode(face[j]);
                                  const Cartesian &v2 = _polyhedron.getNode(face[(j + 1) % 3]);

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

    PlanePropertyVector GravityModel::calculateAlphaSingularityTerms(
            const CartesianSegmentPropertyVector &gij,
            const SegmentPropertyVector &segmentNormalOrientation,
            const CartesianPlanePropertyVector &orthogonalProjectionPointsOnPlane,
            const PlanePropertyVector &planeDistances) {
        //The result
        PlanePropertyVector alphaSingularity(planeDistances.size(), 0.0);

        //Zip iterator consisting of G_ij vectors | sigma_pq | faces | P' | h_p
        auto first = thrust::make_zip_iterator(thrust::make_tuple(gij.begin(),
                                                                  segmentNormalOrientation.begin(),
                                                                  _polyhedron.getFaces().begin(),
                                                                  orthogonalProjectionPointsOnPlane.begin(),
                                                                  planeDistances.begin()));
        auto last = thrust::make_zip_iterator(thrust::make_tuple(gij.end(),
                                                                 segmentNormalOrientation.end(),
                                                                 _polyhedron.getFaces().end(),
                                                                 orthogonalProjectionPointsOnPlane.end(),
                                                                 planeDistances.end()));

        thrust::transform(first, last, alphaSingularity.begin(), [&](const auto &tuple) {
            double alphaSingularityForPlane{0.0};
            const auto &gijPerPlane = thrust::get<0>(tuple);
            const auto segmentNormalOrientationPerPlane = thrust::get<1>(tuple);
            const auto &face = thrust::get<2>(tuple);
            const Cartesian &pPrime = thrust::get<3>(tuple);
            const double hp = thrust::get<4>(tuple);


            return alphaSingularityForPlane;
        });

        return alphaSingularity;
    }

    CartesianPlanePropertyVector GravityModel::calculateBetaSingularityTerms(
            const CartesianSegmentPropertyVector &gij,
            const SegmentPropertyVector &segmentNormalOrientation,
            const CartesianPlanePropertyVector &orthogonalProjectionPointsOnPlane,
            const PlanePropertyVector &planeDistances,
            const PlanePropertyVector &planeNormalOrientation,
            const CartesianPlanePropertyVector &planeUnitNormals) {

        return polyhedralGravity::CartesianPlanePropertyVector();
    }


}