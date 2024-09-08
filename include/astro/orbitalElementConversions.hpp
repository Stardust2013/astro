/*
 * Copyright (c) 2014-2022 Kartik Kumar (me@kartikkumar.com)
 * Copyright (c) 2014-2016 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2016 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2016 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "astro/stateVectorIndices.hpp"

namespace astro
{

//! Convert Cartesian elements to Keplerian elements.
/*!
 * Converts a given set of Cartesian elements (position, velocity) to classical (osculating)
 * Keplerian elements. See Vallado (2007) for a derivation of the conversion.
 *
 * The tolerance is set to a default value. It should not be changed unless required for specific
 * scenarios. Below this tolerance value for eccentricity and inclination, the orbit is considered
 * to be a limit case.
 *
 * Essentially, special solutions are then used for parabolic, circular inclined,
 * non-circular equatorial, and circular equatorial orbits. These special solutions are
 * required because of singularities in the classical Keplerian elements. If high precision is
 * required near these singularities, users are encouraged to consider using other elements, such
 * as Modified Equinoctial Elements (MEE). It should be noted that MEE also suffer from
 * singularities, but not for zero eccentricity and inclination.
 *
 * WARNING: If eccentricity is 1.0 within tolerance, keplerianElements(0) = semi-latus rectum,
 *          since the orbit is parabolic.
 * WARNING: If eccentricity is 0.0 within tolerance, argument of periapsis is set to NaN, since the
 *          orbit is circular.
 * WARNING: If inclination is 0.0 within tolerance, longitude of ascending node is set to NaN, since
 *          since the orbit is equatorial.
 * WARNING: If inclination is 0.0 within tolerance, longitude of ascending node is set to NaN, since
 *          since the orbit is equatorial.
 * WARNING: If eccentricity = 0.0 within tolerance, inclination = non-zero
 *          keplerianElements(5) = argument of latitude
 * WARNING: If eccentricity = non-zero, inclination = 0.0 within tolerance
 *          keplerianElements(5) = true longitude of periapsis
 * WARNING: If eccentricity = 0.0 within tolerance, inclination = 0.0 within tolerance
 *          keplerianElements(5) = true latitude
 *
 * @tparam  Real                    Real type
 * @tparam  Vector                  Vector type
 * @param   cartesianElements       Vector containing Cartesian elements                        <br>
 *                                  N.B.: Order of elements is important!                       <br>
 *                                  cartesianElements(0) = x-position coordinate  [m]           <br>
 *                                  cartesianElements(1) = y-position coordinate  [m]           <br>
 *                                  cartesianElements(2) = z-position coordinate  [m]           <br>
 *                                  cartesianElements(3) = x-velocity coordinate  [m/s]         <br>
 *                                  cartesianElements(4) = y-velocity coordinate  [m/s]         <br>
 *                                  cartesianElements(5) = z-velocity coordinate  [m/s]
 * @param   gravitationalParameter  Gravitational parameter of central body       [m^3 s^-2]
 * @param   tolerance               Tolerance used to check for limit cases
 *                                  (eccentricity, inclination)
 * @return                          Converted vector of Keplerian elements                      <br>
 *                                  N.B.: Order of elements is important!                       <br>
 *                                  keplerianElements(0) = semiMajorAxis          [m]           <br>
 *                                  keplerianElements(1) = eccentricity           [-]           <br>
 *                                  keplerianElements(2) = inclination            [rad]         <br>
 *                                  keplerianElements(3) = argument of periapsis  [rad]         <br>
 *                                  keplerianElements(4) = longitude of
 *                                                         ascending node         [rad]         <br>
 *                                  keplerianElements(5) = true anomaly           [rad]
 */
template <typename Real, typename Vector6>
Vector6 convertCartesianToKeplerianElements(
    const Vector6& cartesianElements,
    const Real gravitationalParameter,
    const Real tolerance = 10.0 * std::numeric_limits<Real>::epsilon())
{
    assert(cartesianElements.size() == 6);
    assert(gravitationalParameter > 0.0);

    const Real pi = 3.14159265358979323846;
    const Real gravitationalParameterInverse = 1.0 / gravitationalParameter;

    typedef std::vector<Real> Vector;
    Vector6 keplerianElements = cartesianElements;

    // Cartesian position
    const Vector position = Vector({cartesianElements[xPositionIndex],
                                    cartesianElements[yPositionIndex],
                                    cartesianElements[zPositionIndex]});
    const Real positionNormSquared = position[0] * position[0]
                                     + position[1] * position[1]
                                     + position[2] * position[2];
    const Real positionNorm = std::sqrt(positionNormSquared);
    const Real positionNormInverse = 1.0 / positionNorm;

    // Cartesian velocity
    const Vector velocity = Vector({cartesianElements[xVelocityIndex],
                                    cartesianElements[yVelocityIndex],
                                    cartesianElements[zVelocityIndex]});
    const Real velocityNormSquared = velocity[0] * velocity[0]
                                     + velocity[1] * velocity[1]
                                     + velocity[2] * velocity[2];
    const Real velocityNorm = std::sqrt(velocityNormSquared);

    // Angular momentum
    const Vector angularMomentum
        = Vector({position[1] * velocity[2] - position[2] * velocity[1],
                    position[2] * velocity[0] - position[0] * velocity[2],
                    position[0] * velocity[1] - position[1] * velocity[0]});
    const Real angularMomentumNormSquared = angularMomentum[0] * angularMomentum[0]
                                            + angularMomentum[1] * angularMomentum[1]
                                            + angularMomentum[2] * angularMomentum[2];
    const Real angularMomentumNorm = std::sqrt(angularMomentumNormSquared);
    const Vector angularMomentumUnitVector
        = Vector({angularMomentum[0] / angularMomentumNorm,
                  angularMomentum[1] / angularMomentumNorm,
                  angularMomentum[2] / angularMomentumNorm});

    // Eccentricity
    const Real eccentricityVectorFirsTermMultiplier
        = velocityNormSquared * gravitationalParameterInverse - positionNormInverse;

    const Real positionDotVelocity = position[0] * velocity[0]
                                     + position[1] * velocity[1]
                                     + position[2] * velocity[2];

    const Real positionDotVelocityScaled = positionDotVelocity / gravitationalParameter;

    const Vector eccentricityVector
    = Vector({eccentricityVectorFirsTermMultiplier * position[0]
                - positionDotVelocity * gravitationalParameterInverse * velocity[0],
              eccentricityVectorFirsTermMultiplier * position[1]
                - positionDotVelocity * gravitationalParameterInverse * velocity[1],
              eccentricityVectorFirsTermMultiplier * position[2]
                - positionDotVelocity * gravitationalParameterInverse * velocity[2]});

    const Real eccentricityNormSquared = eccentricityVector[0] * eccentricityVector[0]
                                            + eccentricityVector[1] * eccentricityVector[1]
                                            + eccentricityVector[2] * eccentricityVector[2];

    const Real eccentricity = std::sqrt(eccentricityNormSquared);
    keplerianElements[1] = eccentricity;

    // Specifc total energy using the vis-viva equation
    const Real specificTotalEnergy
        = 0.5 * velocityNormSquared - gravitationalParameter / positionNorm;

    // Semi-major axis
    Real semiMajorAxis = std::numeric_limits<Real>::quiet_NaN();
    Real semiLatusRectum = std::numeric_limits<Real>::quiet_NaN();
    if (std::fabs(eccentricity - 1.0) > tolerance)
    {
        semiMajorAxis = -0.5 * gravitationalParameter / specificTotalEnergy;
        semiLatusRectum = semiMajorAxis * (1.0 - eccentricityNormSquared);
        keplerianElements[0] = semiMajorAxis;
    }
    else
    {
        semiMajorAxis = std::numeric_limits<Real>::infinity();
        semiLatusRectum = angularMomentumNormSquared * gravitationalParameterInverse;
        keplerianElements[0] = semiLatusRectum;
    }
    assert(semiMajorAxis > 0.0);

    // Inclination
    const Real inclination = std::acos(angularMomentumUnitVector[2]);
    keplerianElements[2] = inclination;

    // Ascending node vector
    const Vector ascendingNodeVector
        = Vector({-angularMomentum[1], angularMomentum[0], 0.0});

    const Real ascendingNodeVectorNormSquared
        = ascendingNodeVector[0] * ascendingNodeVector[0]
          + ascendingNodeVector[1] * ascendingNodeVector[1]
          + ascendingNodeVector[2] * ascendingNodeVector[2];

    const Real ascendingNodeVectorNorm = std::sqrt(ascendingNodeVectorNormSquared);

    const Vector ascendingNodeUnitVector
        = Vector({ascendingNodeVector[0] / ascendingNodeVectorNorm,
                    ascendingNodeVector[1] / ascendingNodeVectorNorm,
                    ascendingNodeVector[2] / ascendingNodeVectorNorm });

    // Longitude of ascending node
    Real longitudeOfAscendingNode = std::acos(ascendingNodeUnitVector[0]);

    if (ascendingNodeVector[1] < 0.0)
    {
        longitudeOfAscendingNode = 2.0 * pi - longitudeOfAscendingNode;
    }

    keplerianElements[4] = longitudeOfAscendingNode;

    // Argument of periapsis
    const Real ascendingNodeVectorDotEccentricityVector
        = ascendingNodeVector[0] * eccentricityVector[0]
          + ascendingNodeVector[1] * eccentricityVector[1]
          + ascendingNodeVector[2] * eccentricityVector[2];

    Real argumentOfPeriapsis = std::acos(ascendingNodeVectorDotEccentricityVector
                                         / (ascendingNodeVectorNorm * eccentricity));

    if (eccentricityVector[2] < 0.0)
    {
        argumentOfPeriapsis = 2.0 * pi - argumentOfPeriapsis;
    }

    keplerianElements[3] = argumentOfPeriapsis;

    // True anomaly
    const Real eccentricityVectorDotPosition = eccentricityVector[0] * position[0]
                                               + eccentricityVector[1] * position[1]
                                               + eccentricityVector[2] * position[2];

    Real trueAnomaly = std::acos(eccentricityVectorDotPosition / (eccentricity * positionNorm));

    if (positionDotVelocity < 0.0)
    {
        trueAnomaly = 2.0 * pi - trueAnomaly;
    }

    keplerianElements[5] = trueAnomaly;

    // True longitude of periapsis
    Real trueLongitudeOfPeriapsis = std::acos(eccentricityVector[0] / eccentricity);

    if (eccentricityVector[1] < tolerance)
    {
            trueLongitudeOfPeriapsis = 2.0 * pi - trueLongitudeOfPeriapsis;
    }

    // Special case: elliptical, equatorial
    if(std::fabs(eccentricity) > tolerance && std::fabs(inclination) < tolerance)
    {
            longitudeOfAscendingNode = std::numeric_limits<Real>::quiet_NaN();
            keplerianElements[4] = trueLongitudeOfPeriapsis;
    }

    // Argument of latitude
    const Real ascendingNodeVectorDotPosition = ascendingNodeVector[0] * position[0]
                                                + ascendingNodeVector[1] * position[1]
                                                + ascendingNodeVector[2] * position[2];

    Real argumentOfLatitude = std::acos(ascendingNodeVectorDotPosition
                                        / (ascendingNodeVectorNorm * positionNorm));

    if (position[2] < 0.0)
    {
        argumentOfLatitude = 2.0 * pi - argumentOfLatitude;
    }

    // Special case: circular, inclined
    if(std::fabs(eccentricity) < tolerance && std::fabs(inclination) > tolerance)
    {
        argumentOfPeriapsis = std::numeric_limits<Real>::quiet_NaN();
        keplerianElements[3] = argumentOfLatitude;
    }

    // True longitude
    Real trueLongitude = std::acos(position[0] / positionNorm);
    if (position[1] < 0.0)
    {
        trueLongitude = 2.0 * pi - trueLongitude;
    }

    // // Special case: circular, equatorial
    if (std::fabs(eccentricity) < tolerance && std::fabs(inclination) < tolerance)
    {
        argumentOfPeriapsis = std::numeric_limits<Real>::quiet_NaN();
        longitudeOfAscendingNode = std::numeric_limits<Real>::quiet_NaN();
        keplerianElements[5] = trueLongitude;
    }

    return keplerianElements;
}

//! Convert Keplerian elements to Cartesian elements.
/*!
 * Converts a given set of Keplerian (osculating) elements to Cartesian elements (position,
 * velocity). See Chobotov (2006) for a full derivation of the conversion.
 *
 * WARNING: If eccentricity is 1.0 within tolerance, the user should provide
 *          keplerianElements(0) = semi-latus rectum, since the orbit is parabolic.
 *
 * @tparam  Real                    Real type
 * @tparam  Vector                  Vector type
 * @param   keplerianElements       Vector containing Keplerian elemenets                       <br>
 *                                  N.B.: Order of elements is important!                       <br>
 *                                  keplerianElements(0) = semiMajorAxis          [m]           <br>
 *                                  keplerianElements(1) = eccentricity           [-]           <br>
 *                                  keplerianElements(2) = inclination            [rad]         <br>
 *                                  keplerianElements(3) = argument of periapsis  [rad]         <br>
 *                                  keplerianElements(4) = longitude of           [rad]
 *                                                          ascending node        [rad]         <br>
 *                                  keplerianElements(5) = true anomaly           [rad]
 * @param   gravitationalParameter  Gravitational parameter of central body       [m^3 s^-2]
 * @param   tolerance               Tolerance used to check for limit case of eccentricity
 * @return                          Converted vector of Cartesian elements                      <br>
 *                                  N.B.: Order of elements is important!                       <br>
 *                                  cartesianElements(0) = x-position coordinate  [m]           <br>
 *                                  cartesianElements(1) = y-position coordinate  [m]           <br>
 *                                  cartesianElements(2) = z-position coordinate  [m]           <br>
 *                                  cartesianElements(3) = x-velocity coordinate  [m/s]         <br>
 *                                  cartesianElements(4) = y-velocity coordinate  [m/s]         <br>
 *                                  cartesianElements(5) = z-velocity coordinate  [m/s]
 */
template <typename Real, typename Vector6>
Vector6 convertKeplerianToCartesianElements(
    const Vector6& keplerianElements, const Real gravitationalParameter,
    const Real tolerance = 10.0 * std::numeric_limits<Real>::epsilon())
{
    Vector6 cartesianElements = keplerianElements;

    const Real semiMajorAxis                    = keplerianElements[semiMajorAxisIndex];
    const Real eccentricity                     = keplerianElements[eccentricityIndex];
    const Real inclination                      = keplerianElements[inclinationIndex];
    const Real argumentOfPeriapsis              = keplerianElements[argumentOfPeriapsisIndex];
    const Real longitudeOfAscendingNode         = keplerianElements[longitudeOfAscendingNodeIndex];
    const Real trueAnomaly                      = keplerianElements[trueAnomalyIndex];

    // Pre-compute sines and cosines of angles for efficient computation.
    const Real cosineOfInclination              = std::cos(inclination);
    const Real sineOfInclination                = std::sin(inclination);
    const Real cosineOfArgumentOfPeriapsis      = std::cos(argumentOfPeriapsis);
    const Real sineOfArgumentOfPeriapsis        = std::sin(argumentOfPeriapsis);
    const Real cosineOfLongitudeOfAscendingNode = std::cos(longitudeOfAscendingNode);
    const Real sineOfLongitudeOfAscendingNode   = std::sin(longitudeOfAscendingNode);
    const Real cosineOfTrueAnomaly              = std::cos(trueAnomaly);
    const Real sineOfTrueAnomaly                = std::sin(trueAnomaly);

    // Compute semi-latus rectum in the case the orbit is not a parabola.
    Real semiLatusRectum = 0.0;
    if (std::fabs(eccentricity - 1.0) > tolerance)
    {
        semiLatusRectum = semiMajorAxis * (1.0 - eccentricity * eccentricity);
    }

    // Else set the semi-latus rectum as the first element in the vector of Keplerian elements.
    else
    {
        semiLatusRectum = keplerianElements[0];
    }

    // Compute the magnitude of the orbital radius, measured from the focal point.
    const Real radiusMagnitude = semiLatusRectum / (1.0 + eccentricity * cosineOfTrueAnomaly);

    // Define position and velocity in the perifocal coordinate system.
    const Real xPositionPerifocal = radiusMagnitude * cosineOfTrueAnomaly;
    const Real yPositionPerifocal = radiusMagnitude * sineOfTrueAnomaly;
    const Real xVelocityPerifocal
        = -std::sqrt(gravitationalParameter / semiLatusRectum) * sineOfTrueAnomaly;
    const Real yVelocityPerifocal
        = std::sqrt(gravitationalParameter / semiLatusRectum)
          * (eccentricity + cosineOfTrueAnomaly);

    // Compute scalar components of rotation matrix to rotate from periforcal to Earth-Centered
    // Inertial (ECI) frame.
    const Real rotationMatrixComponent11
        = (cosineOfLongitudeOfAscendingNode * cosineOfArgumentOfPeriapsis
            -sineOfLongitudeOfAscendingNode * sineOfArgumentOfPeriapsis * cosineOfInclination);
    const Real rotationMatrixComponent12
        = (-cosineOfLongitudeOfAscendingNode * sineOfArgumentOfPeriapsis
            -sineOfLongitudeOfAscendingNode * cosineOfArgumentOfPeriapsis * cosineOfInclination);

    const Real rotationMatrixComponent21
        = (sineOfLongitudeOfAscendingNode * cosineOfArgumentOfPeriapsis
            + cosineOfLongitudeOfAscendingNode * sineOfArgumentOfPeriapsis * cosineOfInclination);
    const Real rotationMatrixComponent22
        = (-sineOfLongitudeOfAscendingNode * sineOfArgumentOfPeriapsis
            + cosineOfLongitudeOfAscendingNode * cosineOfArgumentOfPeriapsis * cosineOfInclination);

    const Real rotationMatrixComponent31 = (sineOfArgumentOfPeriapsis * sineOfInclination);
    const Real rotationMatrixComponent32 = (cosineOfArgumentOfPeriapsis * sineOfInclination);

    // Compute Cartesian position and velocities.
    cartesianElements[xPositionIndex] = rotationMatrixComponent11 * xPositionPerifocal
                                          + rotationMatrixComponent12 * yPositionPerifocal;

    cartesianElements[yPositionIndex] = rotationMatrixComponent21 * xPositionPerifocal
                                          + rotationMatrixComponent22 * yPositionPerifocal;

    cartesianElements[zPositionIndex] = rotationMatrixComponent31 * xPositionPerifocal
                                          + rotationMatrixComponent32 * yPositionPerifocal;

    cartesianElements[xVelocityIndex] = rotationMatrixComponent11 * xVelocityPerifocal
                                          + rotationMatrixComponent12 * yVelocityPerifocal;

    cartesianElements[yVelocityIndex] = rotationMatrixComponent21 * xVelocityPerifocal
                                          + rotationMatrixComponent22 * yVelocityPerifocal;

    cartesianElements[zVelocityIndex] = rotationMatrixComponent31 * xVelocityPerifocal
                                          + rotationMatrixComponent32 * yVelocityPerifocal;

    return cartesianElements;
}

//! Convert true anomaly to elliptical eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical orbits (0 <= eccentricity < 1.0).
 *
 * This implementation performs a check on the eccentricity and throws an error if it is
 * non-elliptical, i.e., eccentricity < 0.0 and eccentricity >= 1.0.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam Real          Real type
 * @param  trueAnomaly   True anomaly                  [rad]
 * @param  eccentricity  Eccentricity                  [-]
 * @return               Elliptical eccentric anomaly  [rad]
 */
template <typename Real>
Real convertTrueAnomalyToEllipticalEccentricAnomaly(const Real trueAnomaly,
                                                    const Real eccentricity)
{
    assert(eccentricity >= 0.0 && eccentricity < 1.0);

    // Compute sine and cosine of eccentric anomaly.
    const Real sineOfEccentricAnomaly
        = std::sqrt(1.0 - std::pow(eccentricity, 2.0))
          * std::sin(trueAnomaly) / (1.0 + eccentricity * std::cos(trueAnomaly));
    const Real cosineOfEccentricAnomaly
        = (eccentricity + std::cos(trueAnomaly))
          / (1.0 + eccentricity * std::cos(trueAnomaly));

    // Return elliptical eccentric anomaly.
    return std::atan2(sineOfEccentricAnomaly, cosineOfEccentricAnomaly);
}

//! Convert true anomaly to hyperbolic eccentric anomaly.
/*!
 * Converts true anomaly to hyperbolic eccentric anomaly for hyperbolic orbits
 * (eccentricity > 1.0).
 *
 * This implementation performs a check on the eccentricity and throws an error if it is
 * non-hyperbolic, i.e., eccentricity >= 1.0.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam Real          Real type
 * @param  trueAnomaly   True anomaly                  [rad]
 * @param  eccentricity  Eccentricity                  [-]
 * @return               Hyperbolic eccentric anomaly  [rad]
 */
template <typename Real>
Real convertTrueAnomalyToHyperbolicEccentricAnomaly(const Real trueAnomaly,
                                                    const Real eccentricity)
{
    assert(eccentricity > 1.0);

    // Compute hyperbolic sine and hyperbolic cosine of hyperbolic eccentric anomaly.
    const Real hyperbolicSineOfHyperbolicEccentricAnomaly
        = std::sqrt(std::pow(eccentricity, 2.0) - 1.0)
          * std::sin(trueAnomaly) / (1.0 + std::cos(trueAnomaly));

    const Real hyperbolicCosineOfHyperbolicEccentricAnomaly
        = (std::cos(trueAnomaly) + eccentricity) / (1.0 + std::cos(trueAnomaly));

    // Return hyperbolic eccentric anomaly.
    // The inverse hyperbolic tangent is computed here manually, since the atanh() function is not
    // available in older C++ compilers.
    const Real angle
        = hyperbolicSineOfHyperbolicEccentricAnomaly
          / hyperbolicCosineOfHyperbolicEccentricAnomaly;
    return 0.5 * (std::log(1.0 + angle) - std::log(1.0 - angle));
}

//! Convert true anomaly to eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical and hyperbolic orbits
 * (eccentricity < 1.0 && eccentricity > 1.0). This function is essentially a wrapper for
 * functions that treat each case. It should be used in cases where the eccentricity of the orbit
 * is not known a priori.
 *
 * This implementation performs a check on the eccentricity and throws an error for
 * eccentricity < 0.0 and parabolic orbits, which have not been implemented.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @sa convertTrueAnomalyToEllipticalEccentricAnomaly,
 *     convertTrueAnomalyToHyperbolicEccentricAnomaly
 * @tparam Real          Real type
 * @param  trueAnomaly   True anomaly       [rad]
 * @param  eccentricity  Eccentricity       [-]
 * @return               Eccentric anomaly  [rad]
 */
template <typename Real>
Real convertTrueAnomalyToEccentricAnomaly(const Real trueAnomaly, const Real eccentricity)
{
    assert(eccentricity >= 0.0
        && (std::fabs(eccentricity - 1.0) > std::numeric_limits<Real>::epsilon()));

    Real eccentricAnomaly = 0.0;

    // Check if orbit is elliptical and compute eccentric anomaly.
    if (eccentricity >= 0.0 && eccentricity < 1.0)
    {
        eccentricAnomaly
            = convertTrueAnomalyToEllipticalEccentricAnomaly(trueAnomaly, eccentricity);
   }

    // Check if orbit is hyperbolic and compute eccentric anomaly.
    else if (eccentricity > 1.0)
    {
        eccentricAnomaly
            = convertTrueAnomalyToHyperbolicEccentricAnomaly(trueAnomaly, eccentricity);
   }

    return eccentricAnomaly;
}

//! Convert elliptical eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical orbits (0 <= eccentricity < 1.0).
 *
 * This implementation performs a check on the eccentricity and throws an error if it is
 * non-elliptical, i.e., eccentricity < 0.0 and eccentricity >= 1.0.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam Real                        Real type
 * @param  ellipticalEccentricAnomaly  Elliptical eccentric anomaly  [rad]
 * @param  eccentricity                Eccentricity                  [-]
 * @return                             Mean anomaly                  [rad]
 */
template <typename Real>
Real convertEllipticalEccentricAnomalyToMeanAnomaly(const Real ellipticalEccentricAnomaly,
                                                    const Real eccentricity)
{
    assert(eccentricity >= 0.0 && eccentricity < 1.0);

    return ellipticalEccentricAnomaly - eccentricity * std::sin(ellipticalEccentricAnomaly);
}

//! Convert hyperbolic eccentric anomaly to mean anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to mean anomaly for hyperbolic orbits
 * (eccentricity > 1.0).
 *
 * This implementation performs a check on the eccentricity and throws an error if it is
 * non-hyperbolic, i.e., eccentricity >= 1.0.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam Real                        Real type
 * @param  hyperbolicEccentricAnomaly  Hyperbolic eccentric anomaly  [rad]
 * @param  eccentricity                Eccentricity                  [-]
 * @return                             Mean anomaly                  [rad]
 */
template <typename Real>
Real convertHyperbolicEccentricAnomalyToMeanAnomaly(
    const Real hyperbolicEccentricAnomaly, const Real eccentricity)
{
    assert(eccentricity > 1.0);

    return eccentricity * std::sinh(hyperbolicEccentricAnomaly) - hyperbolicEccentricAnomaly;
}

//! Convert eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical and hyperbolic orbits
 * (eccentricity < 1.0 && eccentricity > 1.0). This function is essentially a wrapper for
 * functions that treat each case. It should be used in cases where the eccentricity of the orbit
 * is not known a priori.
 *
 * Currently, this implementation performs a check on the eccentricity and throws an error for
 * eccentricity < 0.0 and parabolic orbits, which have not been implemented.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @sa convertEllipticalEccentricAnomalyToMeanAnomaly,
 *     convertHyperbolicEccentricAnomalyToMeanAnomaly
 * @tparam Real              Real type
 * @param  eccentricity      Eccentricity       [-]
 * @param  eccentricAnomaly  Eccentric anomaly  [rad]
 * @return                   Mean anomaly       [rad]
 */
template <typename Real>
Real convertEccentricAnomalyToMeanAnomaly(
    const Real eccentricAnomaly, const Real eccentricity)
{
    assert(eccentricity >= 0.0
        && std::fabs(eccentricity - 1.0) > std::numeric_limits<Real>::epsilon());

    Real meanAnomaly = 0.0;

    // Check if orbit is elliptical and compute mean anomaly.
    if (eccentricity >= 0.0 && eccentricity < 1.0)
    {
        meanAnomaly
            = convertEllipticalEccentricAnomalyToMeanAnomaly(eccentricAnomaly, eccentricity);
   }

    // Check if orbit is hyperbolic and compute mean anomaly.
    else if (eccentricity > 1.0)
    {
        meanAnomaly
            = convertHyperbolicEccentricAnomalyToMeanAnomaly(eccentricAnomaly, eccentricity);
   }

    return meanAnomaly;
}

//! Convert elliptical eccentric anomaly to true anomaly.
/*!
 * Converts eccentric anomaly to true anomaly for elliptical orbits (0 <= eccentricity < 1.0).
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam  Real                      Real type
 * @param   ellipticEccentricAnomaly  Elliptical eccentric anomaly  [rad]
 * @param   eccentricity              Eccentricity                  [-]
 * @return                            True anomaly                  [rad]
 */
template <typename Real>
Real convertEllipticalEccentricAnomalyToTrueAnomaly(const Real ellipticEccentricAnomaly,
                                                    const Real eccentricity)
{
    const Real sineOfTrueAnomaly = std::sqrt(1.0 - eccentricity * eccentricity)
                                   * std::sin(ellipticEccentricAnomaly)
                                   / (1.0 - eccentricity * std::cos(ellipticEccentricAnomaly));

    const Real cosineOfTrueAnomaly
        = (std::cos(ellipticEccentricAnomaly) - eccentricity)
          / (1.0 - eccentricity * std::cos(ellipticEccentricAnomaly));

    return std::atan2(sineOfTrueAnomaly, cosineOfTrueAnomaly);
}

//! Convert hyperbolic eccentric anomaly to true anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to true anomaly for hyperbolic orbits
 * (eccentricity > 1.0). The equations used can be found in (Chobotov, 2002).
 *
 * @tparam  Real                        Real type
 * @param   hyperbolicEccentricAnomaly  Hyperbolic eccentric anomaly  [rad]
 * @param   eccentricity                Eccentricity                  [-]
 * @return                              True anomaly                  [rad]
 */
template <typename Real>
Real convertHyperbolicEccentricAnomalyToTrueAnomaly(const Real hyperbolicEccentricAnomaly,
                                                    const Real eccentricity)
{
    const Real sineOfTrueAnomaly
        = std::sqrt(eccentricity * eccentricity - 1.0)
          * std::sinh(hyperbolicEccentricAnomaly)
          / (eccentricity * std::cosh(hyperbolicEccentricAnomaly) - 1.0);

    const Real cosineOfTrueAnomaly
        = (eccentricity - std::cosh(hyperbolicEccentricAnomaly))
          / (eccentricity * std::cosh(hyperbolicEccentricAnomaly) - 1.0);

    return std::atan2(sineOfTrueAnomaly, cosineOfTrueAnomaly);
}

//! Convert eccentric anomaly to true anomaly.
/*!
 * Converts eccentric anomaly to true anomaly for elliptical and hyperbolic orbits
 * (eccentricity < 1.0 && eccentricity > 1.0). This function is essentially a wrapper for
 * functions that treat each case. It should be used in cases where the eccentricity of the orbit
 * is not known a priori.
 *
 * Currently, this implementation performs a check on the eccentricity and throws an error for
 * eccentricity < 0.0 and parabolic orbits, which have not been implemented.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @sa convertEllipticalEccentricAnomalyToTrueAnomaly,
 *     convertHyperbolicEccentricAnomalyToTrueAnomaly
 * @tparam  Real              Real type
 * @param   eccentricAnomaly  Eccentric anomaly  [rad]
 * @param   eccentricity      Eccentricity       [-]
 * @return                    True anomaly       [rad]
 */
template <typename Real>
Real convertEccentricAnomalyToTrueAnomaly(const Real eccentricAnomaly, const Real eccentricity)
{
    assert(eccentricity >= 0.0
        && std::fabs(eccentricity - 1.0) > std::numeric_limits<Real>::epsilon());

    Real trueAnomaly = 0.0;

    // Check if orbit is elliptical and compute true anomaly.
    if (eccentricity < 1.0)
    {
        trueAnomaly = convertEllipticalEccentricAnomalyToTrueAnomaly(eccentricAnomaly,
                                                                     eccentricity);
   }

    else if (eccentricity > 1.0)
    {
        trueAnomaly = convertHyperbolicEccentricAnomalyToTrueAnomaly(eccentricAnomaly,
                                                                     eccentricity);
   }

    return trueAnomaly;
}

//! Compute Kepler function for elliptical orbits.
/*!
 * Computes Kepler function, given as:
 *
 * \f[
 *      f(E) = E - e * sin(E) - M
 * \f]
 *
 * for elliptical orbits, where \f$E\f$ is the eccentric anomaly, \f$e\f$ is the
 * eccentricity, \f$M\f$ is the mean anomaly.
 *
 * This function can be used for root-finding for mean-to-eccentric anomaly conversion.
 *
 * All elliptical eccentricities >= 0.0 and < 1.0 are valid.
 *
 * @sa convertEllipticalMeanAnomalyToEccentricAnomaly
 * @tparam    Real              Real type
 * @param     eccentricAnomaly  Eccentric anomaly                                 [rad]
 * @param     eccentricity      Eccentricity                                      [-]
 * @param     meanAnomaly       Mean anomaly                                      [rad]
 * @return                      Kepler equation value for given elliptical orbit  [rad]
 */
template <typename Real>
Real computeEllipticalKeplerFunction(const Real eccentricAnomaly,
                                     const Real eccentricity,
                                     const Real meanAnomaly)
{
    return eccentricAnomaly - eccentricity * std::sin(eccentricAnomaly) - meanAnomaly;
}

//! Compute 1st-derivative of Kepler function for elliptical orbits.
/*!
 * Computes the 1st-derivative of Kepler function, given as:
 *
 * \f[
 *      \frac{df(E)} {dE} = 1 - e * cos(E)
 * \f]
 *
 * for elliptical orbits, where \f$E\f$ is the eccentric anomaly, and \f$e\f$ is the
 * eccentricity.
 *
 * All elliptical eccentricities >= 0.0 and < 1.0 are valid.
 *
 * @tparam    Real              Real type
 * @param     eccentricAnomaly  Eccentric anomaly                                            [rad]
 * @param     eccentricity      Eccentricity                                                 [-]
 * @return                      First-derivative of Kepler's function for elliptical orbits  [rad]
 */
template <typename Real>
Real computeFirstDerivativeEllipticalKeplerFunction(const Real eccentricAnomaly,
                                                    const Real eccentricity)
{
    return 1.0 - eccentricity * std::cos(eccentricAnomaly);
}

//! Convert elliptical mean anomaly to eccentric anomaly.
/*!
 * Converts mean anomaly to eccentric anomaly for elliptical orbits,
 * for all eccentricities >= 0.0 and < (1.0 - tolerance).
 *
 * If the conversion fails or the eccentricity falls outside the valid range, then a runtime
 * exception is thrown.
 *
 * The tolerance (1.0e-11) for the upper bound of the eccentricty range is based on the fact that
 * the root-finding algorithm struggles with near-parabolic orbits. This is an experimental finding
 * based on extensive research into the behaviour of the root-finding process (Musegaas, 2012).
 *
 * Also, note that the mean anomaly is automatically transformed to fit within the 0 to 2.0pi range.
 *
 * @sa computeEllipticalKeplerFunction, computeFirstDerivativeEllipticalKeplerFunction
 * @tparam    Real                  Real number type
 * @tparam    Integer               Integer type
 * @param     eccentricity          Eccentricity                                               [-]
 * @param     meanAnomaly           Mean anomaly                                               [rad]
 * @param     rootFindingTolerance  Stopping condition tolernace for Newton-Raphson algorithm  [rad]
 * @param     maximumIterations     Maximum iteration for Newton-Raphson algorithm             [-]
 * @return                          Eccentric anomaly                                          [rad]
 */
template <typename Real, typename Integer>
Real convertEllipticalMeanAnomalyToEccentricAnomaly(
    const Real      eccentricity,
    const Real      meanAnomaly,
    const Real      rootFindingTolerance = 1.0e-3 * std::numeric_limits<Real>::epsilon(),
    const Integer   maximumIterations = 100)
{
    assert(eccentricity >= 0.0 && eccentricity < (1.0 - 1.0e-11));

    const Real pi = 3.14159265358979323846;

    // Set mean anomaly to domain between 0 and 2pi.
    Real meanAnomalyShifted = std::fmod(meanAnomaly, 2.0 * pi);
    if (meanAnomalyShifted <0.0)
    {
        meanAnomalyShifted += 2.0 * pi;
    }

    Real eccentricAnomaly = std::numeric_limits<Real>::quiet_NaN();

    // Set the initial guess for the eccentric anomaly.
    // !!!!!!!!!!!!!     IMPORTANT     !!!!!!!!!!!!!
    // If this scheme is changed, please run a very extensive test suite. The Newton-Raphson
    // root-finding scheme tends to be chaotic for very specific combinations of mean anomaly
    // and eccentricity. Various random tests of 100,000,000 samples were done to verify the
    // default scheme for the initial guess. (Musegaas, 2012).
    Real initialGuess = 0.0;
    if (meanAnomalyShifted > pi)
    {
        initialGuess =  meanAnomalyShifted - eccentricity;
    }
    else
    {
        initialGuess = meanAnomalyShifted + eccentricity;
    }

    // Execute Newton-Raphson root-finding algorithm.
    eccentricAnomaly = initialGuess;
    for (int i = 0; i < maximumIterations + 1; ++i)
    {
        if (i == maximumIterations)
        {
            throw std::runtime_error(
                "ERROR: Maximum iterations for Newton-Raphson root-finding exceeded!");
        }

        const Real nextEccentricAnomaly
            = eccentricAnomaly
              - computeEllipticalKeplerFunction(eccentricAnomaly,
                                                eccentricity,
                                                meanAnomalyShifted)
              / computeFirstDerivativeEllipticalKeplerFunction(eccentricAnomaly,
                                                               eccentricity);

        const Real eccentricAnomalyDifference = eccentricAnomaly - nextEccentricAnomaly;

        eccentricAnomaly = nextEccentricAnomaly;

        if (eccentricAnomalyDifference < rootFindingTolerance)
        {
            break;
        }
    }

    // Return eccentric anomaly.
    return eccentricAnomaly;
}

//! Convert elliptical mean anomaly to eccentric anomaly using a binary search.
/*!
 * Converts mean anomaly to eccentric anomaly for elliptical orbits,
 * for all eccentricities >= 0.0 and < 1.0.
 *
 * If the eccentricity falls outside the valid range, then a runtime
 * exception is thrown.
 *
 * A binary search with 100 iterations should yield about 32 bits of precision.
 *
 * Also, note that the mean anomaly is automatically transformed to fit within the 0 to 2.0pi range.
 *
 * @tparam    Real                  Real number type
 * @tparam    Integer               Integer type
 * @param     eccentricity          Eccentricity                                               [-]
 * @param     meanAnomaly           Mean anomaly                                               [rad]
 * @param     iterations            Binary search iterations                                   [-]
 * @return                          Eccentric anomaly                                          [rad]
 */
template <typename Real, typename Integer>
Real convertEllipticalMeanAnomalyToEccentricAnomalyBS(
    const Real      eccentricity,
    const Real      meanAnomaly,
    const Integer   iterations = 100)
{
    assert(eccentricity >= 0.0 && eccentricity < 1.0);

    const Real pi = 3.14159265358979323846;

    // Set mean anomaly to domain between 0 and 2pi.
    Real meanAnomalyShifted = std::fmod(meanAnomaly, 2.0 * pi);
    if (meanAnomalyShifted <0.0)
    {
        meanAnomalyShifted += 2.0 * pi;
    }

    Real M = meanAnomalyShifted;
    Real E = eccentricity;
    Real E0;
    Real D;

    // TODO: sign(M). This is used twice in the algorithm, and could be extracted to a separate function.
    Real F = (Real(0) < M) - (M < Real(0));
    M = fabs(M) / (2.0 * pi);
    M = (M - (int)M) * 2.0 * pi * F;
    if (M < 0)
    {
        M = M + 2.0 * pi;
    }

    F = 1.0;
    if (M > pi)
    {
        F = -1.0;
        M = 2.0 * pi - M;
    }

    E0 = pi * 0.5;
    D = pi * 0.25;
    for (int J = 0; J < iterations; J++)
    {
        Real M1 = E0 - E * sin(E0);
        Real val = M - M1;
        E0 = E0 + D * ((Real(0) < val) - (val < Real(0)));
        D *= 0.5;
    }

    E0 *= F;

    return E0;
}

} // namespace astro

/*!
 * References
 *  Vallado, D.A.. Fundamentals of Astrodynamics and Applications. Third Edition, Microcosm Press,
 *      2007.
 *  Musegaas, P. Optimization of Space Trajectories Including Multiple Gravity Assists and Deep
 *      Space Maneuvers. MSc thesis, Delft University of Technology, 2013.
 *  Meeus, J. Astronomical Algorithms. Second Edition, Willmann-Bell, Inc. 1991.
 */
