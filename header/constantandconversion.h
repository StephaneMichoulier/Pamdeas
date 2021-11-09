#ifndef CONSTANTANDCONVERSION_H_INCLUDED
#define CONSTANTANDCONVERSION_H_INCLUDED

#include <iostream>

using namespace std;


#define         au 149597870700.   // Astronomical unit [m]
#define       msol 1.98847e30      // Solar mass [kg]
#define          G 6.67430e-11     // Gravitational constant [m³/kg/s²]
#define     rossby 3.              // Rossby number 
#define kboltzmann 1.3806488e-23   // Boltzmann constant kb [J/K]
#define   mgasmean 3.852450297e-27 // Mean molecular mass of gas [kg]
#define   sigmamol 2.0e-19         // Molecular cross-section of gas in PPDs [m²]
#define  M_SQRT_PI 1.7724538509055 // Square root of PI, undefined in <cmath>
#define     cratio -0.5801454844   // Common ratio for a power related to h&s regime
#define      b_oku 0.15            // Parameter b (Okuzumi et al. 2012)
#define maxpacking 0.74048         // Max packing for hexagonal close packing 
#define     deltaR 1e-5            // DeltaR to compute derivative d/dr [m]
#define       year 31556952.       // Year [s]


/* ------------------------- CONVERSIONS ------------------------- */

// Convert AU in meter
inline double AUtoMeter(const double& R)
{   return R*au;    }

// Convert meter in AU
inline double MeterToAU(const double& R)
{   return R/au;    }

// Convert Msol in kg
inline double MsolToKg(const double& mass)
{   return mass*msol;   }

// Convert kg in Msol
inline double KgToMsol(const double& mass)
{   return mass/msol;   }

// Convert time from second to year
inline double SecToYear(const double& time)
{   return time/year;   }

// Convert time from year to second
inline double YearToSec(const double& time)
{   return time*year;   }

#endif // CONSTANTANDCONVERSION_H_INCLUDED