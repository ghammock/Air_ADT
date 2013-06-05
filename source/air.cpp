/******************************************************************************
||  air.cpp      (implementation file)                                       ||
||===========================================================================||
||                                                                           ||
||    Author: Gary Hammock                                                   ||
||    Creation Date:  2010-02-08                                             ||
||    Last Edit Date: 2013-06-04                                             ||
||                                                                           ||
||===========================================================================||
||  DESCRIPTION                                                              ||
||===========================================================================||
||    This abstract data type is used to calculate the properties of         ||
||    equilibrium air for pressures between 1E-4 atm and 1E2 atm and for     ||
||    temperatures between 0 K and 30000 K.                                  ||
||                                                                           ||
||    The reference stated evaluates the properties in SI units.             ||
||                                                                           ||
||===========================================================================||
||  CODE REQUIREMENTS                                                        ||
||===========================================================================||
||    airCoefficients.h                                                      ||
||    air.h                                                                  ||
||                                                                           ||
||===========================================================================||
||  REFERENCES                                                               ||
||===========================================================================||
||    Gupta, R., K. Lee, R. Thompson, J. Yos.  "Calculations and Curve Fits  ||
||        of Thermodynamic and Transport Properties for Equilibrium Air to   ||
||        30000 K".  NASA Reference Publication 1260.  October 1991.         ||
||                                                                           ||
||    Moran, Michael J., Howard N. Shapiro.  "Fundamentals of Engineering    ||
||        Thermodynamics".  5th Edition.  John Wiley and Sons.  Hoboken, NJ, ||
||        2004.  ISBN 0-471-27471-2.                                         ||
||                                                                           ||
||    Liepmann, H. W., A. Roshko.  "Elements of Gasdynamics".  Dover         ||
||        Publications.  2001.  (original copyright: New York.  John Wiley & ||
||        Sons.  1957.)  ISBN 978-0-486-41963-3.                             ||
||                                                                           ||
||===========================================================================||
||  LICENSE    (MIT/X11 License)                                             ||
||===========================================================================||
||    Copyright (C) 2013 Gary Hammock                                        ||
||                                                                           ||
||    Permission is hereby granted, free of charge, to any person obtaining  ||
||    a copy of this software and associated documentation files (the        ||
||    "Software"), to deal in the Software without restriction, including    ||
||    without limitation the rights to use, copy, modify, merge, publish,    ||
||    distribute, sublicense, and/or sell copies of the Software, and to     ||
||    permit persons to whom the Software is furnished to do so, subject to  ||
||    the following conditions:                                              ||
||                                                                           ||
||    The above copyright notice and this permission notice shall be         ||
||    included in all copies or substantial portions of the Software.        ||
||                                                                           ||
||    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,        ||
||    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF     ||
||    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. ||
||    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY   ||
||    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,   ||
||    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      ||
||    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                 ||
||                                                                           ||
******************************************************************************/

/**
 *  @file air.cpp
 *  @author Gary Hammock, PE
 *  @date 2013-06-04
*/

#include "air.h"
#include "airCoefficients.h"

// Define the universal gas constant:
//    8.314462175 kJ/kgmol-K
const double Air::_R_univ = 8.314462175;

/******************************************************
**           Constructors / Destructors              **
******************************************************/

/** Default constructor.  */
Air::Air()
  : _temperature(0.0), _pressure(0.0), _enthalpy(0.0), _intEnergy(0.0),
    _density(0.0), _cp(0.0), _gamma(0.0), _k(0.0), _pr(0.0), _mu(0.0),
    _nu(0.0), _comp(0.0), _gasConstant(0.0), _molarMass(0.0),
    _entropy(0.0), _soundSpeed(0.0), _refraction(0.0),
    _gibbsEnergy(0.0), _helmholtzEn(0.0), _chemPoten(0.0)
{}

/** Copy constructor.
 *
 *  @pre none.
 *  @post A new object is created from the copied values.
 *  @param copyFrom An Air object whose values are to be copied.
 *  @return none.
*/
Air::Air (const Air &copyFrom)
  : _temperature(copyFrom._temperature), _pressure(copyFrom._pressure),
    _enthalpy(copyFrom._enthalpy), _intEnergy(copyFrom._intEnergy),
    _density(copyFrom._density), _cp(copyFrom._cp),
    _gamma(copyFrom._gamma), _k(copyFrom._k), _pr(copyFrom._pr),
    _mu(copyFrom._mu), _nu(copyFrom._nu),
    _comp(copyFrom._comp), _gasConstant(copyFrom._gasConstant),
    _molarMass(copyFrom._molarMass), _entropy(copyFrom._entropy),
    _soundSpeed(copyFrom._soundSpeed), _refraction(copyFrom._refraction),
    _gibbsEnergy(copyFrom._gibbsEnergy), _helmholtzEn(copyFrom._helmholtzEn),
    _chemPoten(copyFrom._chemPoten)
{}

/** Initialization constructor.
 *
 *  @pre none.
 *  @post A new object is created using the input pressure and
 *        temperature to calculate the remaining thermodynamic
 *        properties of the state.
 *  @param pressure The pressure of the state in MPa.
 *  @param temperature The temperature of the state in K.
*/
Air::Air (double pressure, double temperature)
{  calculateProperties(pressure, temperature);  }

/** Default destructor.  */
Air::~Air() {}

/******************************************************
**               Accessors / Mutators                **
******************************************************/

////////////////////
//    Getters
////////////////////

/** Retrieve the temperature of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _temperature [units: K].
*/
double Air::getTemperature (void) const
{  return _temperature;  }

/** Retrieve the pressure of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _pressure [units: MPa].
*/
double Air::getPressure (void) const
{  return _pressure;  }

/** Retrieve the enthalpy of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _enthalpy [units: kJ/kg].
*/
double Air::getEnthalpy (void) const
{  return _enthalpy;  }

/** Retrieve the internal energy of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _intEnergy [units: kJ/kg].
*/
double Air::getInternalEnergy (void) const
{  return _intEnergy;  }

/** Retrieve the density of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _density [units: kg/m^3].
*/
double Air::getDensity (void) const
{  return _density;  }

/** Retrieve the isobaric specific heat of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _cp [units: kJ/kg-K].
*/
double Air::getSpecificHeat (void) const
{  return _cp;  }

/** Retrieve the ratio of specific heats of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _gamma [-dimensionless-].
*/
double Air::getGamma (void) const
{  return _gamma;  }

/** Retrieve the thermal conductivity of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _k [units: W/m-K].
*/
double Air::getThermalConductivity (void) const
{  return _k;  }

/** Retrieve the Prandtl number of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _pr [-dimensionless-].
*/
double Air::getPrandtlNumber (void) const
{  return _pr;  }

/** Retrieve the dynamic viscosity of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _mu [units: g/cm-s].
*/
double Air::getDynamicViscosity (void) const
{  return _mu;  }

/** Retrieve the kinematic viscosity of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _nu [units: m^2/s].
*/
double Air::getKinematicViscosity (void) const
{  return _nu;  }

/** Retrieve the compressibility factor of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _comp [-dimensionless-].
*/
double Air::getCompressibilityFactor (void) const
{  return _comp;  }

/** Retrieve the gas constant of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _gasConstant [-dimensionless-].
*/
double Air::getGasConstant (void) const
{  return _gasConstant;  }

/** Retrieve the molar mass of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _molarMass [units: kg/kgmol].
*/
double Air::getMolarMass (void) const
{  return _molarMass;  }

/** Retrieve the specific entropy of the state (note: in this
 *  implementation, entropy is a derived quantity).
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _entropy [units: kJ/kg-K]
*/
double Air::getEntropy (void) const
{  return _entropy;  }

/** Retrieve the calculated speed of sound of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _soundSpeed [units: m/s]
*/
double Air::getSoundSpeed (void) const
{  return _soundSpeed;  }

/** Retrieve the calculated index of refraction of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _refraction [-dimensionless-]
*/
double Air::getRefractionIndex (void) const
{  return _refraction;  }

/** Retrieve the specific Gibbs free energy (enthalpy) of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _gibbsEnergy [units: kJ/kg]
*/
double Air::getGibbsFreeEnergy (void) const
{  return _gibbsEnergy;  }

/** Retrieve the specific Helmholtz free energy of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _helmholtzEn [units: kJ/kg]
*/
double Air::getHelmholtzFreeEnergy (void) const
{  return _helmholtzEn;  }

/** Retrieve the chemical potential of the state.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _chemPoten [units: kJ/kgmol].
*/
double Air::getChemicalPotential (void) const
{  return _chemPoten;  }

////////////////////
//    Setters
////////////////////

/** Reset all the values to zero.
 *
 *  @pre The object is instantiated.
 *  @post All of the values are reinitialized to zero.
 *  @return none.
*/
void Air::reset (void)
{
    _temperature = 0.0;  // Air temperature [units: K]
    _pressure    = 0.0;  // Air pressure [units: MPa]
    _enthalpy    = 0.0;  // Air enthalpy [units: kJ/kg]
    _intEnergy   = 0.0;  // Specific internal energy [units: kJ/kg]
    _density     = 0.0;  // Air density [units: kg/m^3]
    _cp          = 0.0;  // Specific heat [units: kJ/kg-K]
    _gamma       = 0.0;  // Ratio of specific heats [-dimensionless-]
    _k           = 0.0;  // Thermal conductivity [units: W/m-K]
    _pr          = 0.0;  // Prandtl number [-dimensionless-]
    _mu          = 0.0;  // Dynamic viscosity [g/cm-s]
    _nu          = 0.0;  // Kinematic viscosity [m^2/s]
    _comp        = 0.0;  // Compressibility factor [-dimensionless-]
    _gasConstant = 0.0;  // Specific gas constant [units: kJ/kg-K]
    _molarMass   = 0.0;  // The substance molar mass [units: kg/kgmol]
    _entropy     = 0.0;  // Air specific entropy [units: kJ/kg-K]
    _soundSpeed  = 0.0;  // The speed of sound of air [units: m/s]
    _refraction  = 0.0;  // The index of refraction of air [-dimensionless-]
    _gibbsEnergy = 0.0;  // The specific Gibbs free energy [units: kJ/kg]
    _helmholtzEn = 0.0;  // The specific Helmholtz free energy [units: kJ/kg]
    _chemPoten   = 0.0;  // The chemical potential of air [units: kJ/kgmol]

    return;
}

/** Calculate the properties of air at the given pressure and temperature.
 *
 *  @pre The object is instantiated.
 *  @post The properties are calculated with values
 *        stored in the appropriate variables.
 *  @param pressure The air pressure of the state (in MPa).
 *  @param temperature The air temperature of the state (in K).
 *  @return true The calculation was performed successfully.
 *  @return false The calculation could not be performed.
*/
bool Air::calculateProperties (double pressure, double temperature)
{
    // Store the input temperature value [units: K]
    _temperature = temperature;

    // The calculation assumes pressure in units of atm,
    // so we need to convert MPa to atm to do a check for validity.
    //   0.101325 = conversion factor MPa -> atm
    _pressure = pressure / 0.101325;

    // Check that the pressure and temperature values
    // are in the appropriate ranges.
    //
    // 1E-4 <= _pressure <= 100 atm
    // 0 <= _temperature <= 30,000 K
    if (   (_pressure < 1E-4)   || (_pressure > 100.0)
        || (_temperature < 0.0) || (_temperature > 30000.0)
        )
    {
        _pressure = _pressure * 0.101325;
        return false;
    }

    // Convert the pressure back to MPa from atm.
    //   0.101325 = conversion factor MPa -> atm
    _pressure = _pressure * 0.101325;

    _enthalpy = _calculateEnthalpy(pressure, temperature);  // Units: kJ/kg
    _cp = _calculateSpecificHeat(pressure, temperature);    // Units: kJ/kg-K
    _k = _calculateThermalCond(pressure, temperature);      // Units: W/m-K
    _mu = _calculateViscosity(pressure, temperature);       // Units: kg/m-s
    _comp = _calculateCompFactor(pressure, temperature);    // -dimensionless-

    ////////////////////////////////////
    // Calculate the remaining thermodynamic properties.
    ////////////////////////////////////

    // Store the molar mass of air [units: kg/kgmol]
    _molarMass = 28.96755;

    // Store the air gas constant in SI units [Units: kJ/(kg-K)]
    //    8.314 = Universal gas constant [units: kJ/kgmol-K]
    _gasConstant = _R_univ / _molarMass;

    // Calculate gamma based on the cp value [-dimensionless-]
    _gamma = _cp / (_cp - _gasConstant);

    // Calculate the density [Units: kg/m^3]
    //    1000.0 = convert MPa -> kPa
    _density = (_pressure * 1000.0) / (_comp * _gasConstant * _temperature);

    // Calculate the internal energy based on the thermodynamic relation:
    //    h = u + p/rho
    //    [units: kJ/kg]
    //    1000.0 = convert MPa -> kPa
    _intEnergy = _enthalpy - (_pressure * 1000.0 / _density);

    // Calculate the thermal conductivity [-dimensionless-]
    //    1000.0 = convert kJ -> J. (So that J/s = W)
    _pr = _mu * _cp * 1000.0 / _k;

    // Calculate the kinematic viscosity [units: m^2/s]
    _nu = _mu / _density;

    // Calculate the entropy of the state [units: kJ/kg-K]
    _entropy = _calculateEntropy();

    // Calculate the speed of sound [units: m/s]
    //    1000.0 = convert kJ -> J
    _soundSpeed = sqrt(_gamma * _gasConstant * _temperature * 1000.0);

    // Calculate the refractive index [-dimensionless-]
    _refraction = _calculateRefractionIndex();

    // Calculate the specific Gibbs free energy (enthalpy) using the relation:
    //    G = H - TS
    //    [units: kJ/kg]
    _gibbsEnergy = _enthalpy - (_temperature * _entropy);

    // Calculate the specific Helmholtz free energy using the relation:
    //    F = U - TS
    //    [units: kJ/kg]
    _helmholtzEn = _intEnergy - (_temperature * _entropy);

    // Calculate the chemical potential of the state which is defined as
    // the Gibbs function (total) divided by the molar amount of substance.
    // (Which is also the specific Gibbs function multiplied by the molar
    // mass of the substance).  [units: kJ/kgmol]
    //
    //    ch = G / n = gm / n = gM
    _chemPoten = _gibbsEnergy * _molarMass;

    return true;
}

/** Calculate the properties of air at the given pressure and enthalpy.
 *
 *  @pre The object is instantiated.
 *  @post The properties are calculated with values
 *        stored in the appropriate variables.
 *  @param pressure The air pressure of the state (in MPa).
 *  @param enthalpy The air enthalpy of the state (in kJ/kg).
 *  @return true The calculation was performed successfully.
 *  @return false The calculation could not be performed.
*/
bool Air::calculateProps_PH (double pressure, double enthalpy)
{
    // Store the input pressure.  [units: MPa]
    _pressure = pressure;

    // Guess an initial temperature based on a nominal STP
    // specific heat (1.005 kJ/kg-K) and the input enthalpy.
    double temperature = enthalpy / 1.005;  // units: K

    // Enthalpy convergence tolerance [units: kJ/kg]
    const double tolerance = 1E-4;

    double calcH = 0.0,      // The calculated enthalpy [units: kJ/kg]
           deltaT = 1000.0,  // Iterative temperature interval [units: K]
           T_under = 0.0,    // Stored temp. from prev. low iteration [K]
           h_under = 0.0;    // Stored enthalpy from prev. low iter. [kJ/kg]

    bool converged = false;

    while (!converged)
    {
        // Calculate the enthalpy based on the guessed temperature
        //    [units: kJ/kg]
        calcH = _calculateEnthalpy(pressure, temperature);

        // Assume convergence if the enthalpy difference is
        //    sufficiently small or if the temperature
        //    differential is also small (this alleviates
        //    jump discontinuities in the curve fit equations).
        if ( (fabs(calcH - enthalpy) < tolerance) || (deltaT < 1E-2) )
            converged = true;

        else
        {
            if ((calcH - enthalpy) > 0)
            {
                deltaT *= 0.5;

                // Use a linear interpolation routine to speed convergence
                temperature = (((temperature - T_under) / (calcH - h_under))
                                  * (enthalpy - h_under)) + T_under;
            }
            else
            {
                T_under = temperature;
                h_under = calcH;

                temperature += deltaT;
            }
        }
    }

    return calculateProperties(pressure, temperature);
}

/******************************************************
**                 Helper Methods                    **
******************************************************/

/** Logarithmicaly interpolate the thermodynamic-property
 *  surface (phi, T) between pressure contours.
 *
 *  @pre The object is instantiated and p1 < _pressure < p2.
 *  @post none.
 *  @param p1 The smaller pressure order-of-magnitude.
 *  @param p2 The larger pressure order-of-magnitude.
 *  @param phi_1 The thermodynamic property at p1.
 *  @param phi_2 The thermodynamic property at p2.
 *  @return The interpolated thermodynamic property.
*/
double Air::_interpolate (double p1, double p2,
                          double phi_1, double phi_2) const
{
    double propValue,
           p = _pressure / 0.101325;

    propValue = (((log10(phi_2) - log10(phi_1)) / (log10(p2) - log10(p1)))
                           * (log10(p) - log10(p1))) + log10(phi_1);

    propValue = pow(10.0, propValue);

    return propValue;
}

/** Determine the index of the enthalpy coefficient array
 *  based on the pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The order-of-magnitude of the pressure of interest.
 *  @param temperature The temperature of the state in K.
 *  @return An index into the _h_coeffs array.
*/
uint32 Air::_get_h_row (double pressure, double temperature) const
{
    uint32 temperatureOffset,
           pressureOffset,
           hOffset;

    // The initial index begins with 0
    if ((pressure / 1E-4) < 10.0)
        pressureOffset = 0;

    // There are 6 entries in the 10^-4 regime,
    //    so we begin the 10^-3 index shifted by 6
    else if ((pressure / 1E-3) < 10.0)
        pressureOffset = 6;

    // There are 6 entries in the 10^-3 regime,
    //    so we begin the 10^-2 index shifted by 12
    else if ((pressure / 1E-2) < 10.0)
        pressureOffset = 12;

    // There are 5 entries in the 10^-2 regime,
    //    so we begin the 10^-1 index shifted by 17
    else if ((pressure / 1E-1) < 10.0)
        pressureOffset = 17;

    // There are 4 entries in the 10^-1 regime,
    //    so we begin the 10^0 index shifted by 21
    else if ((pressure / 1E0) < 10.0)
        pressureOffset = 21;

    // There are 4 entries in the 10^0 regime,
    //    so we begin the 10^1 index shifted by 25
    else if ((pressure / 1E1) < 10.0)
        pressureOffset = 25;

    // There are 4 entries in the 10^1 regime,
    //    so we begin the 10^2 index shifted by 29
    else
        pressureOffset = 29;

    ///////////////////////////////////
    // Pressure ~ 10^-4 atm
    ////////////////////
    if (pressureOffset == 0)
    {
             if (temperature <  2250.0) temperatureOffset = 0;
        else if (temperature <  4250.0) temperatureOffset = 1;
        else if (temperature <  6750.0) temperatureOffset = 2;
        else if (temperature < 10750.0) temperatureOffset = 3;
        else if (temperature < 17750.0) temperatureOffset = 4;
        else                            temperatureOffset = 5;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-3 atm
    ////////////////////
    else if (pressureOffset == 6)
    {
             if (temperature <  2250.0) temperatureOffset = 0;
        else if (temperature <  4250.0) temperatureOffset = 1;
        else if (temperature <  6750.0) temperatureOffset = 2;
        else if (temperature < 11750.0) temperatureOffset = 3;
        else if (temperature < 18750.0) temperatureOffset = 4;
        else                            temperatureOffset = 5;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-2 atm
    ////////////////////
    else if (pressureOffset == 12)
    {
             if (temperature <  2750.0) temperatureOffset = 0;
        else if (temperature <  5250.0) temperatureOffset = 1;
        else if (temperature <  9750.0) temperatureOffset = 2;
        else if (temperature < 17750.0) temperatureOffset = 3;
        else                            temperatureOffset = 4;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-1 atm
    ////////////////////
    else if (pressureOffset == 17)
    {
             if (temperature <  3250.0) temperatureOffset = 0;
        else if (temperature <  6250.0) temperatureOffset = 1;
        else if (temperature < 15250.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    ///////////////////////////////////
    // Pressure ~ 10^0 atm
    ////////////////////
    else if (pressureOffset == 21)
    {
             if (temperature <  3750.0) temperatureOffset = 0;
        else if (temperature <  8250.0) temperatureOffset = 1;
        else if (temperature < 17750.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    ///////////////////////////////////
    // Pressure ~ 10^1 atm
    ////////////////////
    else if (pressureOffset == 25)
    {
             if (temperature <  4250.0) temperatureOffset = 0;
        else if (temperature <  9250.0) temperatureOffset = 1;
        else if (temperature < 18750.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    ///////////////////////////////////
    // Pressure ~ 10^2 atm
    ////////////////////
    else
    {
             if (temperature <  6250.0) temperatureOffset = 0;
        else if (temperature < 12750.0) temperatureOffset = 1;
        else                            temperatureOffset = 2;
    }

    hOffset = pressureOffset + temperatureOffset;

    return hOffset;
}

/** Determine the index of the specific heat coefficient array
 *  based on the pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The order-of-magnitude of the pressure of interest.
 *  @param temperature The temperature of the state in K.
 *  @return An index into the _cp_coeffs array.
*/
uint32 Air::_get_cp_row (double pressure, double temperature) const
{
    uint32 temperatureOffset,
           pressureOffset,
           cpOffset;

    // The initial index begins with 0
    if ((pressure / 1E-4) < 10.0)
        pressureOffset = 0;

    // There are 9 entries in the 10^-4 regime,
    //    so we begin the 10^-3 index shifted by 9
    else if ((pressure / 1E-3) < 10.0)
        pressureOffset = 9;

    // There are 8 entries in the 10^-3 regime,
    //    so we begin the 10^-2 index shifted by 17
    else if ((pressure / 1E-2) < 10.0)
        pressureOffset = 17;

    // There are 7 entries in the 10^-2 regime,
    //    so we begin the 10^-1 index shifted by 24
    else if ((pressure / 1E-1) < 10.0)
        pressureOffset = 24;

    // There are 8 entries in the 10^-1 regime,
    //    so we begin the 10^0 index shifted by 32
    else if ((pressure / 1E0) < 10.0)
        pressureOffset = 32;

    // There are 7 entries in the 10^0 regime,
    //    so we begin the 10^1 index shifted by 39
    else if ((pressure / 1E1) < 10.0)
        pressureOffset = 39;

    // There are 7 entries in the 10^1 regime,
    //    so we begin the 10^2 index shifted by 46
    else
        pressureOffset = 46;

    ///////////////////////////////////
    // Pressure ~ 10^-4 atm
    ////////////////////
    if (pressureOffset == 0)
    {
             if (temperature <  1250.0) temperatureOffset = 0;
        else if (temperature <  1750.0) temperatureOffset = 1;
        else if (temperature <  2750.0) temperatureOffset = 2;
        else if (temperature <  4750.0) temperatureOffset = 3;
        else if (temperature <  6250.0) temperatureOffset = 4;
        else if (temperature <  9750.0) temperatureOffset = 5;
        else if (temperature < 14250.0) temperatureOffset = 6;
        else if (temperature < 19750.0) temperatureOffset = 7;
        else                            temperatureOffset = 8;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-3 atm
    ////////////////////
    else if (pressureOffset == 9)
    {
             if (temperature <  1250.0) temperatureOffset = 0;
        else if (temperature <  2250.0) temperatureOffset = 1;
        else if (temperature <  3750.0) temperatureOffset = 2;
        else if (temperature <  5250.0) temperatureOffset = 3;
        else if (temperature <  7250.0) temperatureOffset = 4;
        else if (temperature < 10750.0) temperatureOffset = 5;
        else if (temperature < 17250.0) temperatureOffset = 6;
        else                            temperatureOffset = 7;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-2 atm
    ////////////////////
    else if (pressureOffset == 17)
    {
             if (temperature <  1750.0) temperatureOffset = 0;
        else if (temperature <  2750.0) temperatureOffset = 1;
        else if (temperature <  4750.0) temperatureOffset = 2;
        else if (temperature <  6750.0) temperatureOffset = 3;
        else if (temperature < 12750.0) temperatureOffset = 4;
        else if (temperature < 19750.0) temperatureOffset = 5;
        else                            temperatureOffset = 6;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-1 atm
    ////////////////////
    else if (pressureOffset == 24)
    {
             if (temperature <  1750.0) temperatureOffset = 0;
        else if (temperature <  2750.0) temperatureOffset = 1;
        else if (temperature <  4250.0) temperatureOffset = 2;
        else if (temperature <  6750.0) temperatureOffset = 3;
        else if (temperature <  9750.0) temperatureOffset = 4;
        else if (temperature < 15750.0) temperatureOffset = 5;
        else if (temperature < 21500.0) temperatureOffset = 6;
        else                            temperatureOffset = 7;
    }

    ///////////////////////////////////
    // Pressure ~ 10^0 atm
    ////////////////////
    else if (pressureOffset == 32)
    {
             if (temperature <  1750.0) temperatureOffset = 0;
        else if (temperature <  3250.0) temperatureOffset = 1;
        else if (temperature <  4750.0) temperatureOffset = 2;
        else if (temperature <  7750.0) temperatureOffset = 3;
        else if (temperature < 11750.0) temperatureOffset = 4;
        else if (temperature < 20500.0) temperatureOffset = 5;
        else                            temperatureOffset = 6;
    }

    ///////////////////////////////////
    // Pressure ~ 10^1 atm
    ////////////////////
    else if (pressureOffset == 39)
    {
             if (temperature <  1750.0) temperatureOffset = 0;
        else if (temperature <  3250.0) temperatureOffset = 1;
        else if (temperature <  5750.0) temperatureOffset = 2;
        else if (temperature <  9250.0) temperatureOffset = 3;
        else if (temperature < 13750.0) temperatureOffset = 4;
        else if (temperature < 22500.0) temperatureOffset = 5;
        else                            temperatureOffset = 6;
    }

    ///////////////////////////////////
    // Pressure ~ 10^2 atm
    ////////////////////
    else
    {
             if (temperature <  1750.0) temperatureOffset = 0;
        else if (temperature <  3750.0) temperatureOffset = 1;
        else if (temperature <  6750.0) temperatureOffset = 2;
        else if (temperature < 10750.0) temperatureOffset = 3;
        else if (temperature < 17750.0) temperatureOffset = 4;
        else                            temperatureOffset = 5;
    }

    cpOffset = pressureOffset + temperatureOffset;

    return cpOffset;
}

/** Determine the index of the thermal conductivity coefficient array
 *  based on the pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The order-of-magnitude of the pressure of interest.
 *  @param temperature The temperature of the state in K.
 *  @return An index into the _k_coeffs array.
*/
uint32 Air::_get_k_row (double pressure, double temperature) const
{
    uint32 temperatureOffset,
           pressureOffset,
           kOffset;

    // The initial index begins with 0
    if ((pressure / 1E-4) < 10.0)
        pressureOffset = 0;

    // There are 7 entries in the 10^-4 regime,
    //    so we begin the 10^-3 index shifted by 7
    else if ((pressure / 1E-3) < 10.0)
        pressureOffset = 7;

    // There are 7 entries in the 10^-3 regime,
    //    so we begin the 10^-2 index shifted by 14
    else if ((pressure / 1E-2) < 10.0)
        pressureOffset = 14;

    // There are 7 entries in the 10^-2 regime,
    //    so we begin the 10^-1 index shifted by 21
    else if ((pressure / 1E-1) < 10.0)
        pressureOffset = 21;

    // There are 6 entries in the 10^-1 regime,
    //    so we begin the 10^0 index shifted by 27
    else if ((pressure / 1E0) < 10.0)
        pressureOffset = 27;

    // There are 6 entries in the 10^0 regime,
    //    so we begin the 10^1 index shifted by 33
    else if ((pressure / 1E1) < 10.0)
        pressureOffset = 33;

    // There are 5 entries in the 10^1 regime,
    //    so we begin the 10^2 index shifted by 38
    else
        pressureOffset = 38;

    ///////////////////////////////////
    // Pressure ~ 10^-4 atm
    ////////////////////
    if (pressureOffset == 0)
    {
             if (temperature <  1750.0) temperatureOffset = 0;
        else if (temperature <  2750.0) temperatureOffset = 1;
        else if (temperature <  4750.0) temperatureOffset = 2;
        else if (temperature <  6250.0) temperatureOffset = 3;
        else if (temperature < 10250.0) temperatureOffset = 4;
        else if (temperature < 17750.0) temperatureOffset = 5;
        else                            temperatureOffset = 6;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-3 atm
    ////////////////////
    else if (pressureOffset == 7)
    {
             if (temperature <  1750.0) temperatureOffset = 0;
        else if (temperature <  2750.0) temperatureOffset = 1;
        else if (temperature <  4750.0) temperatureOffset = 2;
        else if (temperature <  6250.0) temperatureOffset = 3;
        else if (temperature < 11250.0) temperatureOffset = 4;
        else if (temperature < 18250.0) temperatureOffset = 5;
        else                            temperatureOffset = 6;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-2 atm
    ////////////////////
    else if (pressureOffset == 14)
    {
             if (temperature <  2250.0) temperatureOffset = 0;
        else if (temperature <  3250.0) temperatureOffset = 1;
        else if (temperature <  5750.0) temperatureOffset = 2;
        else if (temperature <  7750.0) temperatureOffset = 3;
        else if (temperature < 12750.0) temperatureOffset = 4;
        else if (temperature < 18750.0) temperatureOffset = 5;
        else                            temperatureOffset = 6;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-1 atm
    ////////////////////
    else if (pressureOffset == 21)
    {
             if (temperature <  2250.0) temperatureOffset = 0;
        else if (temperature <  4250.0) temperatureOffset = 1;
        else if (temperature <  6750.0) temperatureOffset = 2;
        else if (temperature <  9250.0) temperatureOffset = 3;
        else if (temperature < 16750.0) temperatureOffset = 4;
        else                            temperatureOffset = 5;
    }

    ///////////////////////////////////
    // Pressure ~ 10^0 atm
    ////////////////////
    else if (pressureOffset == 27)
    {
             if (temperature <  2250.0) temperatureOffset = 0;
        else if (temperature <  4250.0) temperatureOffset = 1;
        else if (temperature <  7750.0) temperatureOffset = 2;
        else if (temperature < 10750.0) temperatureOffset = 3;
        else if (temperature < 19250.0) temperatureOffset = 4;
        else                            temperatureOffset = 5;
    }

    ///////////////////////////////////
    // Pressure ~ 10^1 atm
    ////////////////////
    else if (pressureOffset == 33)
    {
             if (temperature <  3250.0) temperatureOffset = 0;
        else if (temperature <  5250.0) temperatureOffset = 1;
        else if (temperature <  8750.0) temperatureOffset = 2;
        else if (temperature < 13750.0) temperatureOffset = 3;
        else                            temperatureOffset = 4;
    }

    ///////////////////////////////////
    // Pressure ~ 10^2 atm
    ////////////////////
    else
    {
             if (temperature <  3750.0) temperatureOffset = 0;
        else if (temperature <  6250.0) temperatureOffset = 1;
        else if (temperature < 10750.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    kOffset = pressureOffset + temperatureOffset;

    return kOffset;
}

/** Determine the index of the viscosity coefficient array
 *  based on the pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The order-of-magnitude of the pressure of interest.
 *  @param temperature The temperature of the state in K.
 *  @return An index into the _mu_coeffs array.
*/
uint32 Air::_get_mu_row (double pressure, double temperature) const
{
    uint32 temperatureOffset,
           pressureOffset,
           muOffset;

    // The initial index begins with 0
    if ((pressure / 1E-4) < 10.0)
        pressureOffset = 0;

    // There are 4 entries in the 10^-4 regime,
    //    so we begin the 10^-3 index shifted by 4
    else if ((pressure / 1E-3) < 10.0)
        pressureOffset = 4;

    // There are 4 entries in the 10^-3 regime,
    //    so we begin the 10^-2 index shifted by 8
    else if ((pressure / 1E-2) < 10.0)
        pressureOffset = 8;

    // There are 4 entries in the 10^-2 regime,
    //    so we begin the 10^-1 index shifted by 12
    else if ((pressure / 1E-1) < 10.0)
        pressureOffset = 12;

    // There are 4 entries in the 10^-1 regime,
    //    so we begin the 10^0 index shifted by 16
    else if ((pressure / 1E0) < 10.0)
        pressureOffset = 16;

    // There are 3 entries in the 10^0 regime,
    //    so we begin the 10^1 index shifted by 19
    else if ((pressure / 1E1) < 10.0)
        pressureOffset = 19;

    // There are 3 entries in the 10^1 regime,
    //    so we begin the 10^2 index shifted by 22
    else
        pressureOffset = 22;

    ///////////////////////////////////
    // Pressure ~ 10^-4 atm
    ////////////////////
    if (pressureOffset == 0)
    {
             if (temperature <  7750.0) temperatureOffset = 0;
        else if (temperature < 10750.0) temperatureOffset = 1;
        else if (temperature < 16750.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-3 atm
    ////////////////////
    else if (pressureOffset == 4)
    {
             if (temperature <  8250.0) temperatureOffset = 0;
        else if (temperature < 12250.0) temperatureOffset = 1;
        else if (temperature < 18750.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-2 atm
    ////////////////////
    else if (pressureOffset == 8)
    {
             if (temperature <  8750.0) temperatureOffset = 0;
        else if (temperature < 14250.0) temperatureOffset = 1;
        else if (temperature < 19750.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-1 atm
    ////////////////////
    else if (pressureOffset == 12)
    {
             if (temperature <  9750.0) temperatureOffset = 0;
        else if (temperature < 16750.0) temperatureOffset = 1;
        else if (temperature < 24500.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    ///////////////////////////////////
    // Pressure ~ 10^0 atm
    ////////////////////
    else if (pressureOffset == 16)
    {
             if (temperature < 11250.0) temperatureOffset = 0;
        else if (temperature < 19750.0) temperatureOffset = 1;
        else                            temperatureOffset = 2;
    }

    ///////////////////////////////////
    // Pressure ~ 10^1 atm
    ////////////////////
    else if (pressureOffset == 19)
    {
             if (temperature < 12750.0) temperatureOffset = 0;
        else if (temperature < 21500.0) temperatureOffset = 1;
        else                            temperatureOffset = 2;
    }


    ///////////////////////////////////
    // Pressure ~ 10^2 atm
    ////////////////////
    else
    {
             if (temperature < 15250.0) temperatureOffset = 0;
        else                            temperatureOffset = 1;
    }


    muOffset = pressureOffset + temperatureOffset;

    return muOffset;
}

/** Determine the index of the compressibility coefficient array
 *  based on the pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The order-of-magnitude of the pressure of interest.
 *  @param temperature The temperature of the state in K.
 *  @return An index into the _z_coeffs array.
*/
uint32 Air::_get_z_row (double pressure, double temperature) const
{
    uint32 temperatureOffset,
           pressureOffset,
           zOffset;

    // The initial index begins with 0
    if ((pressure / 1E-4) < 10.0)
        pressureOffset = 0;

    // There are 5 entries in the 10^-4 regime,
    //    so we begin the 10^-3 index shifted by 5
    else if ((pressure / 1E-3) < 10.0)
        pressureOffset = 5;

    // There are 5 entries in the 10^-3 regime,
    //    so we begin the 10^-2 index shifted by 10
    else if ((pressure / 1E-2) < 10.0)
        pressureOffset = 10;

    // There are 5 entries in the 10^-2 regime,
    //    so we begin the 10^-1 index shifted by 15
    else if ((pressure / 1E-1) < 10.0)
        pressureOffset = 15;

    // There are 5 entries in the 10^-1 regime,
    //    so we begin the 10^0 index shifted by 20
    else if ((pressure / 1E0) < 10.0)
        pressureOffset = 20;

    // There are 5 entries in the 10^0 regime,
    //    so we begin the 10^1 index shifted by 25
    else if ((pressure / 1E1) < 10.0)
        pressureOffset = 25;

    // There are 4 entries in the 10^1 regime,
    //    so we begin the 10^2 index shifted by 29
    else
        pressureOffset = 29;

    ///////////////////////////////////
    // Pressure ~ 10^-4 atm
    ////////////////////
    if (pressureOffset == 0)
    {
             if (temperature <  2750.0) temperatureOffset = 0;
        else if (temperature <  5750.0) temperatureOffset = 1;
        else if (temperature <  8750.0) temperatureOffset = 2;
        else if (temperature < 17750.0) temperatureOffset = 3;
        else                            temperatureOffset = 4;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-3 atm
    ////////////////////
    else if (pressureOffset == 5)
    {
             if (temperature <  3250.0) temperatureOffset = 0;
        else if (temperature <  6750.0) temperatureOffset = 1;
        else if (temperature <  9750.0) temperatureOffset = 2;
        else if (temperature < 19750.0) temperatureOffset = 3;
        else                            temperatureOffset = 4;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-2 atm
    ////////////////////
    else if (pressureOffset == 10)
    {
             if (temperature <  3250.0) temperatureOffset = 0;
        else if (temperature <  7250.0) temperatureOffset = 1;
        else if (temperature < 11750.0) temperatureOffset = 2;
        else if (temperature < 21500.0) temperatureOffset = 3;
        else                            temperatureOffset = 4;
    }

    ///////////////////////////////////
    // Pressure ~ 10^-1 atm
    ////////////////////
    else if (pressureOffset == 15)
    {
             if (temperature <  3750.0) temperatureOffset = 0;
        else if (temperature <  8250.0) temperatureOffset = 1;
        else if (temperature < 13750.0) temperatureOffset = 2;
        else if (temperature < 23500.0) temperatureOffset = 3;
        else                            temperatureOffset = 4;
    }

    ///////////////////////////////////
    // Pressure ~ 10^0 atm
    ////////////////////
    else if (pressureOffset == 20)
    {
             if (temperature <  5750.0) temperatureOffset = 0;
        else if (temperature <  9250.0) temperatureOffset = 1;
        else if (temperature < 15750.0) temperatureOffset = 2;
        else if (temperature < 23500.0) temperatureOffset = 3;
        else                            temperatureOffset = 4;
    }

    ///////////////////////////////////
    // Pressure ~ 10^1 atm
    ////////////////////
    else if (pressureOffset == 25)
    {
             if (temperature <  5750.0) temperatureOffset = 0;
        else if (temperature <  9750.0) temperatureOffset = 1;
        else if (temperature < 17250.0) temperatureOffset = 2;
        else                            temperatureOffset = 3;
    }

    ///////////////////////////////////
    // Pressure ~ 10^2 atm
    ////////////////////
    else
    {
             if (temperature <  8750.0) temperatureOffset = 0;
        else if (temperature < 17750.0) temperatureOffset = 1;
        else                            temperatureOffset = 2;
    }

    zOffset = pressureOffset + temperatureOffset;

    return zOffset;
}

/** Determine the upper and lower order of magnitudes
 *  based on the input pressure.
 *
 *  @pre none.
 *  @post none.
 *  @param pressure The pressure of interest in atm.
 *  @param pLower The lower order of magnitude pressure in atm.
 *  @param pUpper The upper order of magnitude pressure in atm.
 *  @return none.
*/
void Air::_getPressureOM (double pressure,
                          double &pLower, double &pUpper) const
{
    if ((pressure / 1E-4) < 10.0)
    {
        pLower = 1E-4;
        pUpper = 1E-3;
    }
    else if ((pressure / 1E-3) < 10.0)
    {
        pLower = 1E-3;
        pUpper = 1E-2;
    }
    else if ((pressure / 1E-2) < 10.0)
    {
        pLower = 1E-2;
        pUpper = 1E-1;
    }
    else if ((pressure / 1E-1) < 10.0)
    {
        pLower = 1E-1;
        pUpper = 1E0;
    }
    else if ((pressure / 1E0) < 10.0)
    {
        pLower = 1E0;
        pUpper = 1E1;
    }
    else
    {
        pLower = 1E1;
        pUpper = 1E2;
    }

    return;
}

/** Calculate the enthalpy using the input pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The pressure of the state in units of MPa.
 *  @param temperature The temperature of the state in K.
 *  @return The calculated enthalpy in units of kJ/kg.
*/
double Air::_calculateEnthalpy (double pressure, double temperature) const
{
    double enth1;

    double p1,  // Smaller order-of-magnitude [units: atm]
           p2;  // Larger order-of-magnitude  [units: atm]

    double x;     // Independent variable for curve fit relations.
    double h[2];  // order-of-magnitude values for enthalpy.

    uint32 index[2];  // Stores the indices into the coefficient tables.

    // The reference states that for temperatures below 500 K,
    // simpler relations may be used to generate properties.
    if (temperature <= 500.0)
        enth1 = 0.24E-3 * temperature;  // Units: kcal/g

    else
    {
        // The calculation assumes pressure in units of atm,
        // so we need to convert MPa to atm.
        //   0.101325 = conversion factor MPa -> atm
        pressure /= 0.101325;

        // Determine the correct order-of-magnitude values
        // for evaluating the properties
        _getPressureOM(pressure, p1, p2);

        // The enthalpy, specific heat, and thermal conductivity
        // curve fits use the natural log of temperature as the
        // independent variable.
        x = log(temperature / 10000.0);

        // Based on the magnitude of the input pressure
        // and the temperature range, return an enthalpy
        // table location index
        index[0] = _get_h_row(p1, temperature);
        index[1] = _get_h_row(p2, temperature);

        for (uint32 i = 0; i < 2; ++i)
        {
            h[i] =  (_h_coeffs[index[i]][0] * pow(x, 4.0))
                  + (_h_coeffs[index[i]][1] * pow(x, 3.0))
                  + (_h_coeffs[index[i]][2] * pow(x, 2.0))
                  + (_h_coeffs[index[i]][3] * x)
                  +  _h_coeffs[index[i]][4];

            h[i] = exp(h[i]);  // [Units: kcal/g]
        }

        // Evaluate the properties by using log-linear interpolation
        // between the pressure intervals specified
        enth1 = _interpolate(p1, p2, h[0], h[1]);   // units: kcal/g
    }

    // Convert enthalpy from kcal/g -> kJ/kg
    //    1000.0   = convert g -> kg
    //    1000.0   = convert kcal -> cal
    //    238.8459 = convert cal -> kJ
    double enthalpy = enth1 * 1000.0 * 1000.0 / 238.8459;

    return enthalpy;
}

/** Calculate the specific heat using the input pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The pressure of the state in units of MPa.
 *  @param temperature The temperature of the state in K.
 *  @return The calculated specific heat in units of kJ/kg-K.
*/
double Air::_calculateSpecificHeat (double pressure, double temperature) const
{
    double cp1;

    double p1,  // Smaller order-of-magnitude [units: atm]
           p2;  // Larger order-of-magnitude  [units: atm]

    double x;      // Independent variable for curve fit relations.
    double cp[2];  // order-of-magnitude values for cp.

    uint32 index[2];  // Stores the indices into the coefficient tables.

    // The reference states that for temperatures below 500 K,
    // simpler relations may be used to generate properties.
    if (temperature <= 500.0)
    {
        // The reference states that for temperatures below
        // 500 K, use the constant value cp = 0.24 cal/(g-K)
        cp1 = 0.24;  // Units: cal/(g-K)
    }

    else
    {
        // The calculation assumes pressure in units of atm,
        // so we need to convert MPa to atm.
        //   0.101325 = conversion factor MPa -> atm
        pressure /= 0.101325;

        // Determine the correct order-of-magnitude values
        // for evaluating the properties
        _getPressureOM(pressure, p1, p2);

        // The enthalpy, specific heat, and thermal conductivity
        // curve fits use the natural log of temperature as the
        // independent variable.
        x = log(temperature / 10000.0);

        // Based on the magnitude of the input pressure
        // and the temperature range, return a specific heat
        // table location index
        index[0] = _get_cp_row(p1, temperature);
        index[1] = _get_cp_row(p2, temperature);

        for (uint32 i = 0; i < 2; ++i)
        {
            cp[i] =  (_cp_coeffs[index[i]][0] * pow(x, 4.0))
                   + (_cp_coeffs[index[i]][1] * pow(x, 3.0))
                   + (_cp_coeffs[index[i]][2] * pow(x, 2.0))
                   + (_cp_coeffs[index[i]][3] * x)
                   +  _cp_coeffs[index[i]][4];

            cp[i] = exp(cp[i]);  // [Units: cal/(g-K)]
        }

        // Evaluate the properties by using log-linear interpolation
        // between the pressure intervals specified
        cp1 = _interpolate(p1, p2, cp[0], cp[1]); // units: cal/g-K
    }

    // Convert specific heat from cal/g-K to kJ/kg-K
    //    1000.0   = convert g -> kg
    //    238.8459 = convert cal -> kJ
    double specificHeat = cp1 * 1000.0 / 238.8459;

    return specificHeat;
}

/** Calculate the thermal cond. using the input pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The pressure of the state in units of MPa.
 *  @param temperature The temperature of the state in K.
 *  @return The calculated thermal cond. in units of W/m-K.
*/
double Air::_calculateThermalCond (double pressure, double temperature) const
{
    double k1;

    double p1,  // Smaller order-of-magnitude [units: atm]
           p2;  // Larger order-of-magnitude  [units: atm]

    double x;     // Independent variable for curve fit relations.
    double k[2];  // order-of-magnitude values for thermal cond.

    uint32 index[2];  // Stores the indices into the coefficient tables.

    // The reference states that for temperatures below 500 K,
    // simpler relations may be used to generate properties.
    if (_temperature <= 500.0)
    {
        // If the temperature is less than 500 K, use Sutherland's
        // thermal conductivity law...  [Units: cal/(cm-s-K)]
        k1 = 5.9776E-6 * (pow(_temperature, 1.5) / (_temperature + 194.4));
    }

    else
    {
        // The calculation assumes pressure in units of atm,
        // so we need to convert MPa to atm.
        //   0.101325 = conversion factor MPa -> atm
        pressure /= 0.101325;

        // Determine the correct order-of-magnitude values
        // for evaluating the properties
        _getPressureOM(pressure, p1, p2);

        // The enthalpy, specific heat, and thermal conductivity
        // curve fits use the natural log of temperature as the
        // independent variable.
        x = log(temperature / 10000.0);

        // Based on the magnitude of the input pressure
        // and the temperature range, return a therm. cond.
        // table location index
        index[0] = _get_k_row(p1, temperature);
        index[1] = _get_k_row(p2, temperature);

        for (uint32 i = 0; i < 2; ++i)
        {
            k[i] =  (_k_coeffs[index[i]][0] * pow(x, 4.0))
                  + (_k_coeffs[index[i]][1] * pow(x, 3.0))
                  + (_k_coeffs[index[i]][2] * pow(x, 2.0))
                  + (_k_coeffs[index[i]][3] * x)
                  +  _k_coeffs[index[i]][4];

            k[i] = exp(k[i]);  // [Units: cal/(cm-s-K)]
        }

        // Evaluate the properties by using log-linear interpolation
        // between the pressure intervals specified
        k1 = _interpolate(p1, p2, k[0], k[1]);   // units: cal/cm-s-K
    }

    // Convert thermal conductivity from cal/cm-s-K to W/m-K
    //    100.0     = convert cm -> m
    //    0.2388459 = convert cal -> J
    //    1 J/s     = 1 W
    double thermalCond = k1 * 100.0 / 0.2388459;

    return thermalCond;
}

/** Calculate the compressibility using the input pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The pressure of the state in units of MPa.
 *  @param temperature The temperature of the state in K.
 *  @return The calculated compressibility factor.
*/
double Air::_calculateCompFactor (double pressure, double temperature) const
{
    double comp;

    double p1,  // Smaller order-of-magnitude [units: atm]
           p2;  // Larger order-of-magnitude  [units: atm]

    double x;     // Independent variable for curve fit relations.
    double z[2];  // order-of-magnitude values for compressibility.

    uint32 index[2];  // Stores the indices into the coefficient tables.

    // The reference states that for temperatures below 500 K,
    // simpler relations may be used to generate properties.
    if (_temperature <= 500.0)
    {
        // The reference states that for temperatures below
        // 500 K, use the compressibility factor of unity
        comp = 1.0;  // -dimensionless-
    }

    else
    {
        // The calculation assumes pressure in units of atm,
        // so we need to convert MPa to atm.
        //   0.101325 = conversion factor MPa -> atm
        pressure /= 0.101325;

        // Determine the correct order-of-magnitude values
        // for evaluating the properties
        _getPressureOM(pressure, p1, p2);

        // The compressibility factor and viscosity curve fits
        // use a scaled temperature as the independent value.
        x = _temperature / 1000.0;

        // Based on the magnitude of the input pressure
        // and the temperature range, return a comp. factor
        // table location index
        index[0] = _get_z_row(p1, temperature);
        index[1] = _get_z_row(p2, temperature);

        for (uint32 i = 0; i < 2; ++i)
        {
            z[i] =   _z_coeffs[index[i]][0]
                  + (_z_coeffs[index[i]][1] * x)
                  + (_z_coeffs[index[i]][2] * pow(x, 2.0))
                  + (_z_coeffs[index[i]][3] * pow(x, 3.0))
                  + (_z_coeffs[index[i]][4] * pow(x, 4.0));
        }

        // Evaluate the properties by using log-linear interpolation
        // between the pressure intervals specified
        comp = _interpolate(p1, p2, z[0], z[1]);   // -dimensionless-
    }

    return comp;
}

/** Calculate the viscosity using the input pressure and temperature.
 *
 *  @pre The object is instantiated and
 *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
 *  @post none.
 *  @param pressure The pressure of the state in units of MPa.
 *  @param temperature The temperature of the state in K.
 *  @return The calculated viscosity in units of kg/m-s.
*/
double Air::_calculateViscosity (double pressure, double temperature) const
{
    double mu1;

    double p1,  // Smaller order-of-magnitude [units: atm]
           p2;  // Larger order-of-magnitude  [units: atm]

    double x;      // Independent variable for curve fit relations.
    double mu[2];  // order-of-magnitude values for viscosity.

    uint32 index[2];  // Stores the indices into the coefficient tables.

    // The reference states that for temperatures below 500 K,
    // simpler relations may be used to generate properties.
    if (_temperature <= 500.0)
    {
        // If the computed film temperature is less than 500 K,
        // we compute the viscosity using Sutherland's Viscosity Law...
        // [Units: poise (g/cm-s)]
        mu1 = 1.4584E-5 * (pow(_temperature, 1.5) / (_temperature + 110.33));
    }

    else
    {
        // The calculation assumes pressure in units of atm,
        // so we need to convert MPa to atm.
        //   0.101325 = conversion factor MPa -> atm
        pressure /= 0.101325;

        // Determine the correct order-of-magnitude values
        // for evaluating the properties
        _getPressureOM(pressure, p1, p2);

        // The compressibility factor and viscosity curve fits
        // use a scaled temperature as the independent value.
        x = _temperature / 1000.0;

        // Based on the magnitude of the input pressure
        // and the temperature range, return a viscosity
        // table location index
        index[0] = _get_mu_row(p1, temperature);
        index[1] = _get_mu_row(p2, temperature);

        for (uint32 i = 0; i < 2; ++i)
        {
            mu[i] =   _mu_coeffs[index[i]][0]
                  + (_mu_coeffs[index[i]][1] * x)
                  + (_mu_coeffs[index[i]][2] * pow(x, 2.0))
                  + (_mu_coeffs[index[i]][3] * pow(x, 3.0))
                  + (_mu_coeffs[index[i]][4] * pow(x, 4.0))
                  + (_mu_coeffs[index[i]][5] * pow(x, 5.0));
        }

        // Evaluate the properties by using log-linear interpolation
        // between the pressure intervals specified
        mu1 = _interpolate(p1, p2, mu[0], mu[1]); // units: poise
    }

    // Convert viscosity from poise to kg/m-s
    //    1000.0 = convert g -> kg
    //    100.0  = convert cm -> m
    double viscosity = mu1 * 100.0 / 1000.0;

    return viscosity;
}

/** Calculate the entropy of the state using the stored pressure,
 *  temperature, specific heat, and gas constant.
 *
 *  @pre The object is instantited and the values for _pressure,
 *       _temperature, _gasConstant, and _cp have been calculated.
 *  @post none.
 *  @return A double precision value containing the calculated
 *          entropy value [units: kJ/kg-K].
*/
double Air::_calculateEntropy (void) const
{
    /******************************************************
    **  REFERENCES                                       **
    ** ------------------------------------------------- **
    **  1.) Moran, Michael J., Howard N. Shapiro.        **
    **      "Fundamentals of Engineering                 **
    **      Thermodynamics".  5th Edition.  John Wiley   **
    **      and Sons.  Hoboken, NJ,  2004.               **
    **      ISBN 0-471-27471-2.                          **
    **                                                   **
    ******************************************************/

    double entropy;

    // Reference temperature, pressure, and entropy values.
    // (From Table A-22, Ref 1.)
    static const double T0 = 300.0,     // Ref. temperature, 300 K.
                        p0 = 0.101325,  // Ref. pressure 0.101325 MPa (1 atm).
                        s0 = 1.70203;   // Ref. entropy [units: kJ/kg-K].

    // Calculate the entropy of the gas state.  This equation assumes a
    // thermally and calorically perfect gas, (it is a limitation) but
    // it's a pretty good assumption for most air cases.
    //
    //    (From Equation 6.23, Ref 1.)

    entropy =   (_cp * log(_temperature / T0))        // Caloric component.
              - (_gasConstant * log(_pressure / p0))  // enthalpic component.
              + s0;                                   // Reference offset.

    return entropy;
}

/** Calculate the refractive index of air using the calculated density.
 *
 *  @pre The object is instantiated and the value for _density
 *       has been computed.
 *  @post none.
 *  @return A double precision value for the calculated index of
 *          refraction of air [-dimensionless-].
*/
double Air::_calculateRefractionIndex (void) const
{
    /******************************************************
    **  REFERENCES                                       **
    ** ------------------------------------------------- **
    **  1.) Liepmann, H. W., A. Roshko.  "Elements of    **
    **      Gasdynamics".  Dover Publications.  2001.    **
    **      (original copyright: New York.  John Wiley & **
    **      Sons.  1957.)  ISBN 978-0-486-41963-3.       **
    **                                                   **
    ******************************************************/

    double refIndex;

    static const double beta = 0.000292,            // Ref 1.
                        rho0 = 1.2925694365458342;  // @ STP [units: kg/m^3].

    // Calculate the refractive index. [-dimensionless-].
    refIndex = 1.0 + (beta * _density / rho0);

    return refIndex;
}