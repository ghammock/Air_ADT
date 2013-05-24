Air ADT
Author: Gary Hammock, PE

================================================================================
                                DESCRIPTION
================================================================================

This C++ abstract data type is used to calculate the thermodynamic and
transport properties of equilibrium air in the pressure range of
1E-4 <= P <= 100 atm and with a temperature range of 0 <= T < 30000 K.

The data are to be read and stored in typical SI units as follows:
    Temperature, K
    Pressure, MPa
    Enthalpy, kJ/kg
    Density, kg/m^3
    Specific Heat, kJ/kg-K
    Ratio of Specific Heats (Gamma), -dimensionless-
    Thermal Conductivity, W/m-K
    Prandtl Number, -dimensionless-
    Dynamic Viscosity, kg/m-s
    Kinematic Viscosity, m^2/s
    Compressibility, -dimensionless-
    Gas Constant, kJ/kg-K
    Speed of Sound, m/s

================================================================================
                               IMPLEMENTATION
================================================================================

To create a new Air object:

                   Air state1;

To create a new initialized Air object using pressure and temperature as input:

                   Air state2(pressure, temperature);

To calculate the properties (state lookup) of air using pressure and
temperature:

                   Air state1.calculateProperties(pressure, temperature);

To calculate the properties (state lookup) of air using pressure and enthalpy:

                   Air state1.calculateProps_PH(pressure, enthalpy);

List of accessor methods used by the ADT:
    double getTemperature (void)
    double getPressure (void)
    double getEnthalpy (void)
    double getDensity (void)
    double getSpecificHeat (void)
    double getGamma (void)
    double getThermalConductivity (void)
    double getPrandtlNumber (void)
    double getDynamicViscosity (void)
    double getKinematicViscosity (void)
    double getCompressibilityFactor (void)
    double getGasConstant (void)
    double getSoundSpeed (void)
    
================================================================================
                              DESIRED UPDATES
================================================================================
I have a legitimate need to have extended thermophysical properties through
1000 atm (10 MPa).  If anyone has generated the curve fit coefficients, please
consider including them (or contacting me).  I have Fought's paper containing
the enthalpy coefficients, but specific heat and thermal conductivity are 
the high priority coefficient sets.

================================================================================
                                 REFERENCES
================================================================================
1.) Gupta, R., K. Lee, R. Thompson, J. Yos.  "Calculations and Curve Fits of
      Thermodynamic and Transport Properties for Equilibrium Air to 30000 K".
      NASA Reference Publication 1260.  October 1991.