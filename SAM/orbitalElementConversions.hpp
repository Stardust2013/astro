/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#ifndef SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP
#define SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP

#include <cmath>
 
#include <SML/sml.hpp>

namespace sam
{

//! Cartesian element array indices.
enum CartesianElementIndices
{
    xPositionIndex = 0,
    yPositionIndex = 1,
    zPositionIndex = 2,
    xVelocityIndex = 3,
    yVelocityIndex = 4,
    zVelocityIndex = 5
};

//! Cartesian element array indices.
enum KeplerianElementIndices
{
    semiMajorAxisIndex = 0,
    semiLatusRectumIndex = 0,    
    eccentricityIndex = 1,
    inclinationIndex = 2,
    argumentOfPeriapsisIndex = 3,
    longitudeOfAscendingNodeIndex = 4,
    trueAnomalyIndex = 5,
    meanAnomalyIndex = 5
};

//! Convert Cartesian elements to Keplerian elements.
/*!
 * Converts a given set of Cartesian elements (position, velocity) to classical (osculating) 
 * Keplerian elements. See Chobotov (2006) for a full derivation of the conversion.
 *
 * The tolerance is given a default value. It should not be changed unless required for specific 
 * scenarios. Below this tolerance value for eccentricity and inclination, the orbit is considered 
 * to be a limit case. Essentially, special solutions are then used for parabolic, circular 
 * inclined, non-circular equatorial, and circular equatorial orbits. These special solutions are 
 * required because of singularities in the classical Keplerian elements. If  high precision is 
 * required near these singularities, users are encouraged to consider using other elements, such 
 * as Modified Equinoctial Elements (MEE). It should be noted that MEE also suffer from 
 * singularities, but not for zero eccentricity and inclination.
 * 
 * WARNING: If eccentricity is 1.0 within 1.0e-15, keplerianElements( 0 ) = semi-latus rectum,
 *          since the orbit is parabolic.
 * WARNING: If eccentricity is 0.0 within 1.0e-15, argument of periapsis is set to 0.0, since the
 *          orbit is circular.
 * WARNING: If inclination is 0.0 within 1.0e-15, longitude of ascending node is set to 0.0, since 
 *          the orbit is equatorial.
 *
 * @tparam  Real              Real type
 * @tparam  Vector6           6-Vector type
 * @tparam  Vector3           3-Vector type
 * @param   cartesianElements Vector containing Cartesian elements. 
 *                            N.B.: Order of elements is important!
 *                            cartesianElements( 0 ) = x-position coordinate [m]
 *                            cartesianElements( 1 ) = y-position coordinate [m]
 *                            cartesianElements( 2 ) = z-position coordinate [m]
 *                            cartesianElements( 3 ) = x-velocity coordinate [m/s]
 *                            cartesianElements( 4 ) = y-velocity coordinate [m/s]
 *                            cartesianElements( 5 ) = z-velocity coordinate [m/s]
 * @param   gravitationalParameter Gravitational parameter of central body [m^3 s^-2]
 * @param   tolerance         Tolerance used to check for limit cases 
 *                            (zero eccentricity, inclination)
 * @return  Converted vector of Keplerian elements.
 *          N.B.: Order of elements is important! 
 *          keplerianElements( 0 ) = semiMajorAxis                [m]
 *          keplerianElements( 1 ) = eccentricity                 [-]
 *          keplerianElements( 2 ) = inclination                  [rad]
 *          keplerianElements( 3 ) = argument of periapsis        [rad]
 *          keplerianElements( 4 ) = longitude of ascending node  [rad]
 *          keplerianElements( 5 ) = true anomaly                 [rad]
 */
template< typename Real, typename Vector6, typename Vector3 >
Vector6 convertCartesianToKeplerianElements(
    const Vector6& cartesianElements, const Real gravitationalParameter,
    const Real tolerance = 10.0 * std::numeric_limits< Real >::epsilon( ) )
{
    Vector6 keplerianElements( 6 );

    // Set position and velocity vectors.
    Vector3 position( 3 );
    position[ 0 ] = cartesianElements[ 0 ];
    position[ 1 ] = cartesianElements[ 1 ];
    position[ 2 ] = cartesianElements[ 2 ];

    Vector3 velocity( 3 );
    velocity[ 0 ] = cartesianElements[ 3 ];
    velocity[ 1 ] = cartesianElements[ 4 ];
    velocity[ 2 ] = cartesianElements[ 5 ];     

    // Compute orbital angular momentum vector.
    const Vector3 angularMomentum( sml::cross( position, velocity ) );

    // Compute semi-latus rectum.
    const Real semiLatusRectum 
        = sml::squaredNorm< Real >( angularMomentum ) / gravitationalParameter;

    // Compute unit vector to ascending node.
    Vector3 ascendingNodeUnitVector 
        = sml::normalize< Real >( 
            sml::cross( sml::getZUnitVector< Vector3 >( ), 
                        sml::normalize< Real >( angularMomentum ) ) );        

    // Compute eccentricity vector.
    Vector3 eccentricityVector 
        = sml::add( sml::multiply( sml::cross( velocity, angularMomentum ), 
                                   1.0 / gravitationalParameter ),
                    sml::multiply( sml::normalize< Real >( position ), -1.0 ) );

    // Store eccentricity.
    keplerianElements[ eccentricityIndex ] = sml::norm< Real >( eccentricityVector );        

    // Compute and store semi-major axis.
    // Check if orbit is parabolic. If it is, store the semi-latus rectum instead of the
    // semi-major axis.
    if ( std::fabs( keplerianElements[ eccentricityIndex ] - 1.0 ) < tolerance )
    {
        keplerianElements[ semiLatusRectumIndex ] = semiLatusRectum;
    }

    // Else the orbit is either elliptical or hyperbolic, so store the semi-major axis.
    else
    {
        keplerianElements[ semiMajorAxisIndex ] 
            = semiLatusRectum / ( 1.0 - keplerianElements[ eccentricityIndex ]
                                  * keplerianElements[ eccentricityIndex ] );
    }

    // Compute and store inclination.
    keplerianElements[ inclinationIndex ]
        = std::acos( angularMomentum[ zPositionIndex ] 
                     / sml::norm< Real >( angularMomentum ) );    

    // Compute and store longitude of ascending node.
    // Define the quadrant condition for the argument of perigee.
    Real argumentOfPeriapsisQuandrantCondition = eccentricityVector[ zPositionIndex ];

    // Check if the orbit is equatorial. If it is, set the vector to the line of nodes to the
    // x-axis.
    if ( std::fabs( keplerianElements[ inclinationIndex ] ) < tolerance )
    {
        ascendingNodeUnitVector = sml::getXUnitVector< Vector3 >( );

        // If the orbit is equatorial, eccentricityVector_z is zero, therefore the quadrant
        // condition is taken to be the y-component, eccentricityVector_y.
        argumentOfPeriapsisQuandrantCondition = eccentricityVector[ yPositionIndex ];
    }

    // Compute and store the resulting longitude of ascending node.
    keplerianElements[ longitudeOfAscendingNodeIndex ] 
        = std::acos( ascendingNodeUnitVector[ xPositionIndex] );

    // Check if the quandrant is correct.
    if ( ascendingNodeUnitVector[ yPositionIndex ] < 0.0 )
    {
        keplerianElements[ longitudeOfAscendingNodeIndex ] 
            = 2.0 * sml::SML_PI - keplerianElements[ longitudeOfAscendingNodeIndex ];
    }

    // Compute and store argument of periapsis.
    // Define the quadrant condition for the true anomaly.
    Real trueAnomalyQuandrantCondition = sml::dot< Real >( position, velocity );

    // Check if the orbit is circular. If it is, set the eccentricity vector to unit vector
    // pointing to the ascending node, i.e. set the argument of periapsis to zero.
    if ( std::fabs( keplerianElements[ eccentricityIndex ] ) < tolerance )
    {
        eccentricityVector = ascendingNodeUnitVector;

        keplerianElements[ argumentOfPeriapsisIndex ] = 0.0;

        // Check if orbit is also equatorial and set true anomaly quandrant check condition
        // accordingly.
        if ( ascendingNodeUnitVector == sml::getXUnitVector< Vector3 >( ) )
        {
            // If the orbit is circular, dot( position, velocity ) = 0, therefore this value
            // cannot be used as a quadrant condition. Moreover, if the orbit is equatorial,
            // position_z is also zero and therefore the quadrant condition is taken to be the
            // y-component, position_y.
            trueAnomalyQuandrantCondition = position[ yPositionIndex ];
        }

        else
        {
            // If the orbit is circular, dot( position, velocity ) = 0, therefore the quadrant
            // condition is taken to be the z-component of the position, position_z.
            trueAnomalyQuandrantCondition = position[ zPositionIndex ];
        }
    }

    // Else, compute the argument of periapsis as the angle between the eccentricity vector and
    // the unit vector to the ascending node.
    else
    {
        keplerianElements[ argumentOfPeriapsisIndex ]
            = std::acos( sml::dot< Real >( sml::normalize< Real >( eccentricityVector ), ascendingNodeUnitVector ) );

        // Check if the quadrant is correct.
        if ( argumentOfPeriapsisQuandrantCondition < 0.0 )
        {
            keplerianElements[ argumentOfPeriapsisIndex ] 
                = 2.0 * sml::SML_PI - keplerianElements[ argumentOfPeriapsisIndex ];
        }
    }    

    // Compute dot-product of position and eccentricity vectors.
    Real dotProductPositionAndEccentricityVectors
        = sml::dot< Real >( sml::normalize< Real >( position ), 
                            sml::normalize< Real >(eccentricityVector ) );

    // Check if the dot-product is one of the limiting cases: 0.0 or 1.0
    // (within prescribed tolerance).
    if ( std::fabs( 1.0 - dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors = 1.0;
    }

    if ( std::fabs( dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors = 0.0;
    }

    // Compute and store true anomaly.
    keplerianElements[ trueAnomalyIndex ] 
        = std::acos( dotProductPositionAndEccentricityVectors );    

    // Check if the quandrant is correct.
    if ( trueAnomalyQuandrantCondition < 0.0 )
    {
        keplerianElements[ trueAnomalyIndex ]
            = 2.0 * sml::SML_PI - keplerianElements[ trueAnomalyIndex ];
    }

    return keplerianElements;
}

// convert true anomaly to eccentric anomaly

// convert eccentric anomaly to mean anomaly

} // namespace sam

#endif // SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP
