/******************************************************************************
||  air.h      (definition file)                                             ||
||===========================================================================||
||                                                                           ||
||    Author: Gary Hammock                                                   ||
||    Creation Date:  2010-02-08                                             ||
||    Last Edit Date: 2014-01-30                                             ||
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
||    air.cpp                                                                ||
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
||    Incropera, Frank P.  David P. DeWitt.  "Fundamentals of Heat and Mass  ||
||        Transfer."  5th Edition.  John Wiley and Sons.  Hoboken, NJ, 2002. ||
||        ISBN 0-471-38650-2.                                                ||
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
 *  @file air.h
 *  @author Gary Hammock, PE
 *  @date 2014-01-30
*/

#ifndef _GH_DEF_AIR_H
#define _GH_DEF_AIR_H

// We need the math library for log() and exp()
#ifndef _CMATH_
#include <cmath>
#endif

///////////////////////////////////////
//    Simplifying type definitions
///////////////////////////////////////
typedef unsigned int uint32;

/**
 *  @class Air An ADT to calculate and store the thermodynamic and
 *         transport properties of equilibrium air.
*/
class Air
{
  public:
    /******************************************************
    **           Constructors / Destructors              **
    ******************************************************/

    /** Default constructor.  */
    Air();

    /** Copy constructor.
     *
     *  @pre none.
     *  @post A new object is created from the copied values.
     *  @param copyFrom An Air object whose values are to be copied.
     *  @return none.
    */
    Air (const Air &copyFrom);

    /** Initialization constructor.
     *
     *  @pre none.
     *  @post A new object is created using the input pressure and
     *        temperature to calculate the remaining thermodynamic
     *        properties of the state.
     *  @param pressure The pressure of the state in MPa.
     *  @param temperature The temperature of the state in K.
    */
    Air (double pressure, double temperature);

    /** Default destructor.  */
    ~Air();

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
    double getTemperature (void) const;

    /** Retrieve the pressure of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _pressure [units: MPa].
    */
    double getPressure (void) const;

    /** Retrieve the enthalpy of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _enthalpy [units: kJ/kg].
    */
    double getEnthalpy (void) const;

    /** Retrieve the internal energy of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _intEnergy [units: kJ/kg].
    */
    double getInternalEnergy (void) const;

    /** Retrieve the density of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _density [units: kg/m^3].
    */
    double getDensity (void) const;

    /** Retrieve the isobaric specific heat of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _cp [units: kJ/kg-K].
    */
    double getSpecificHeat (void) const;

    /** Retrieve the ratio of specific heats of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _gamma [-dimensionless-].
    */
    double getGamma (void) const;

    /** Retrieve the thermal conductivity of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _k [units: W/m-K].
    */
    double getThermalConductivity (void) const;

    /** Retrieve the Prandtl number of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _pr [-dimensionless-].
    */
    double getPrandtlNumber (void) const;

    /** Retrieve the dynamic viscosity of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _mu [units: kg/m-s].
    */
    double getDynamicViscosity (void) const;

    /** Retrieve the kinematic viscosity of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _nu [units: m^2/s].
    */
    double getKinematicViscosity (void) const;

    /** Retrieve the compressibility factor of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _comp [-dimensionless-].
    */
    double getCompressibilityFactor (void) const;

    /** Retrieve the gas constant of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _gasConstant [units: kJ/kg-K].
    */
    double getGasConstant (void) const;

    /** Retrieve the molar mass of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _molarMass [units: kg/kgmol].
    */
    double getMolarMass (void) const;

    /** Retrieve the specific entropy of the state (note: in this
     *  implementation, entropy is a derived quantity).
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _entropy [units: kJ/kg-K]
    */
    double getEntropy (void) const;

    /** Retrieve the calculated speed of sound of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _soundSpeed [units: m/s]
    */
    double getSoundSpeed (void) const;

    /** Retrieve the calculated index of refraction of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _refraction [-dimensionless-]
    */
    double getRefractionIndex (void) const;

    /** Retrieve the specific Gibbs free energy (enthalpy) of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _gibbsEnergy [units: kJ/kg]
    */
    double getGibbsFreeEnergy (void) const;

    /** Retrieve the specific Helmholtz free energy of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _helmholtzEn [units: kJ/kg]
    */
    double getHelmholtzFreeEnergy (void) const;

    /** Retrieve the chemical potential of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _chemPoten [units: kJ/kgmol].
    */
    double getChemicalPotential (void) const;

    /** Retrieve the Schmidt number of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _schmidt [-dimensionless-].
    */
    double getSchmidtNumber (void) const;

    /** Retrieve the Lewis number of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _lewis [-dimensionless-].
    */
    double getLewisNumber (void) const;

    /** Retrieve the thermal diffusivity of the state.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _thermalDiff [units: m^s/s].
    */
    double getThermalDiffusivity (void) const;

    ////////////////////
    //    Setters
    ////////////////////

    /** Reset all the values to zero.
     *
     *  @pre The object is instantiated.
     *  @post All of the values are reinitialized to zero.
     *  @return none.
    */
    void reset (void);

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
    bool calculateProperties (double pressure, double temperature);

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
    bool calculateProps_PH (double pressure, double enthalpy);

  private:
    /******************************************************
    **                     Members                       **
    ******************************************************/
    double _temperature,  // Air temperature [units: K]
           _pressure,     // Air pressure [units: MPa]
           _enthalpy,     // Air enthalpy [units: kJ/kg]
           _intEnergy,    // Specific internal energy [units: kJ/kg]
           _density,      // Air density [units: kg/m^3]
           _cp,           // Specific heat [units: kJ/kg-K]
           _gamma,        // Ratio of specific heats [-dimensionless-]
           _k,            // Thermal conductivity [units: W/m-K]
           _pr,           // Prandtl number [-dimensionless-]
           _mu,           // Dynamic viscosity [kg/m-s]
           _nu,           // Kinematic viscosity [m^2/s]
           _comp,         // Compressibility factor [-dimensionless-]
           _gasConstant,  // Specific gas constant [units: kJ/kg-K]
           _molarMass,    // The substance molar mass [units: kg/kgmol]
           _entropy,      // Air specific entropy [units: kJ/kg-K]
           _soundSpeed,   // The speed of sound of air [units: m/s]
           _refraction,   // The index of refraction of air [-dimensionless-]
           _gibbsEnergy,  // The specific Gibbs free energy [units: kJ/kg].
           _helmholtzEn,  // The specific Helmholtz free energy [units: kJ/kg]
           _chemPoten,    // The chemical potential of air [units: kJ/kgmol]
           _schmidt,      // The Schmidt number (Sc) [-dimensionless-]
           _lewis,        // The Lewis number (Le) [-dimensionless-]
           _thermalDiff;  // Thermal diffusivity (alpha) [units: m^2/s]

    static const double _R_univ;  // Universal gas constant [units: kJ/kgmol-K]

    // These coefficients are initialized in "airCoefficients.h"
    static const double _h_coeffs[32][5];
    static const double _cp_coeffs[52][5];
    static const double _k_coeffs[42][5];
    static const double _mu_coeffs[24][6];
    static const double _z_coeffs[32][5];

    /******************************************************
    **                 Helper Methods                    **
    ******************************************************/

    /** Logarithmicaly interpolate the thermodynamic-property
     *  surface (phi, T) between pressure contours.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @param p1 The smaller pressure order-of-magnitude.
     *  @param p2 The larger pressure order-of-magnitude.
     *  @param phi_1 The thermodynamic property at p1.
     *  @param phi_2 The thermodynamic property at p2.
     *  @return The interpolated thermodynamic property.
    */
    double _interpolate (double p1, double p2,
                         double phi_1, double phi_2) const;

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
    uint32 _get_h_row (double pressure, double temperature) const;

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
    uint32 _get_cp_row (double pressure, double temperature) const;

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
    uint32 _get_k_row (double pressure, double temperature) const;

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
    uint32 _get_mu_row (double pressure, double temperature) const;

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
    uint32 _get_z_row (double pressure, double temperature) const;

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
    void _getPressureOM (double pressure,
                         double &pLower, double &pUpper) const;

    /** Calculate the enthalpy using the input pressure and temperature.
     *
     *  @pre The object is instantiated and
     *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
     *  @post none.
     *  @param pressure The pressure of the state in units of MPa.
     *  @param temperature The temperature of the state in K.
     *  @return The calculated enthalpy in units of kJ/kg.
    */
    double _calculateEnthalpy (double pressure, double temperature) const;

    /** Calculate the specific heat using the input pressure and temperature.
     *
     *  @pre The object is instantiated and
     *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
     *  @post none.
     *  @param pressure The pressure of the state in units of MPa.
     *  @param temperature The temperature of the state in K.
     *  @return The calculated specific heat in units of kJ/kg-K.
    */
    double _calculateSpecificHeat (double pressure, double temperature) const;

    /** Calculate the thermal cond. using the input pressure and temperature.
     *
     *  @pre The object is instantiated and
     *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
     *  @post none.
     *  @param pressure The pressure of the state in units of MPa.
     *  @param temperature The temperature of the state in K.
     *  @return The calculated thermal cond. in units of W/m-K.
    */
    double _calculateThermalCond (double pressure, double temperature) const;

    /** Calculate the compressibility using the input pressure and temperature.
     *
     *  @pre The object is instantiated and
     *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
     *  @post none.
     *  @param pressure The pressure of the state in units of MPa.
     *  @param temperature The temperature of the state in K.
     *  @return The calculated compressibility factor.
    */
    double _calculateCompFactor (double pressure, double temperature) const;

    /** Calculate the viscosity using the input pressure and temperature.
     *
     *  @pre The object is instantiated and
     *       10E-4 <= P <= 100 atm, 0 <= T <= 30000 K.
     *  @post none.
     *  @param pressure The pressure of the state in units of MPa.
     *  @param temperature The temperature of the state in K.
     *  @return The calculated viscosity in units of kg/m-s.
    */
    double _calculateViscosity (double pressure, double temperature) const;

    /** Calculate the entropy of the state using the stored pressure,
     *  temperature, specific heat, and gas constant.
     *
     *  @pre The object is instantited and the values for _pressure,
     *       _temperature, _gasConstant, and _cp have been calculated.
     *  @post none.
     *  @return A double precision value containing the calculated
     *          entropy value [units: kJ/kg-K].
    */
    double _calculateEntropy (void) const;

    /** Calculate the refractive index of air using the calculated density.
     *
     *  @pre The object is instantiated and the value for _density
     *       has been computed.
     *  @post none.
     *  @return A double precision value for the calculated index of
     *          refraction of air [-dimensionless-].
    */
    double _calculateRefractionIndex (void) const;

};  // end class Air

#endif