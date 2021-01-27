#ifndef POROSITY_H_INCLUDED
#define POROSITY_H_INCLUDED

#include <iostream>

using namespace std;


/* ------------------------ SHARED FUNCTIONS ------------------------ */

// A parameter shared by all transition masses and sizes
double CParam(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& eroll,// ->
              const double& a0, const double& rhos, const double& alpha);

// Return the porosity = 1 - filling factor
double Porosity(const double& phifinal);


/* ------------------------ BOUNCE ------------------------*/

// Compute the probability to bounce
double ProbaBounce(const double& phi, const double& philimbounce, const double& vrel, const double& vstick, const double& vend);

// Compute the variation in volume after bouncing [m³]
double VarVolumeBounce(const double& phi, const double& philim, const double& coeffrest, const double& ekin,// ->
                       const double& volume, const double& Yd0, const double& Ydpower);

// Compute the filling factor after bouncing
double PhiBounce(const double& sizei, const double& phii, const double& rhos, const double& philim,// ->
                 const double& vrel, const double& vstick, const double& vyield, const double& ncoll,// ->
                 const double& Yd0, const double& Ydpower);


/* ------------------------ MASS MODEL------------------------ */

/* ------------------------ OKUZUMI GROWTH MODEL ------------------------ */
              
// Transition mass M1/m0 (Phi h&s = Phi Ep-St<1)
double M1(const double& cparam);

// Transition mass M2/m0 (Phi h&s = Phi St-St<1)
double M2(const double& cparam, const double& rhog, const double& cg, const double& a0);

// Transition mass M3/m0 (Phi Ep-St<1 = Phi St-St<1)
double M3(const double& m1_on_m0, const double& m2_on_m0);
double M3(const double& cparam, const double& rhog, const double& cg, const double& a0);

// Transition mass M4/m0 (Phi Ep-St<1 = Phi St=1)
double M4(const double& R, const double cparam, const double& Mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos);

// Transition mass M5/m0 (Phi St-St<1 = Phi St=1)
double M5(const double& R, const double cparam, const double& Mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos);

// Filling factor du to growth
double PhiMGr(const double& massf, const double& massi, const double& phii,// ->
              const double& rhos, const double& eroll, const double& vrel);

// Filling factor due to collision [Okuzumi]
double PhiMColl(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& st,// ->
                const double& massf, const double eroll, const double& a0, const double& rhos, const double& alpha);


/* ------------------------ KATAOKA ------------------------*/

// Filling factor due to gas compression
double PhiMGas(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& deltav,// ->
               const double& massf, const double& rhos, const double& eroll, const double& a0, const int& iregime);

// Filling factor due to self-gravity
double PhiMGrav(const double& massf, const double& rhos, const double& a0, const double& eroll);


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

// Return the minimum filling factor between PhiMColl, PhiMGas and PhiMgrav
double PhiMinMColGasGrav(const double& R, const double& Mstar, const double& rhog, const double cg, const double& deltav,// ->
                         const double st, const double& massf, const double& rhos, const double& eroll, const double& a0,// ->
                         const double& alpha, const int& iregime);

// Return the minimum filling factor between phiMGr, PhiMinMColGasGrav and PhiMBounce
double PhiMFinal(const double& R, const double& Mstar, const double& rhog, const double& cg, const double deltav, const double st,// ->
                 const double& massf, const double& massi, const double& phii, const double& a0, const double& rhos, const double& eroll,// ->
                 const double& alpha, const double& ncoll, const int& ifrag, const int& ibounce, const double& vfrag, const double& vrel,// ->
                 const double& vstick, const double& probabounce, const double& philim, const double& Yd0, const double& Ydpower,// ->
                 const int& iregime);


/* ------------------------ SIZE MODEL ------------------------ */

/* ------------------------ OKUZUMI GROWTH MODEL ------------------------ */

// Transition size S1/a0 (Phi h&s = Phi Ep-St<1)
double S1(const double& cparam);

// Transition size S2/a0 (Phi h&s = Phi St-St<1)
double S2(const double& cparam, const double& rhog, const double& cg, const double& a0);

// Transition size S3/a0 (Phi Ep-St<1 = Phi St-St<1)
double S3(const double& rhog, const double& cg, const double& a0);

// Transition size S4/a0 (Phi Ep-St<1 = Phi St=1)
double S4(const double& R, const double cparam, const double& Mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos);

// Transition size S5/a0 (Phi St-St<1 = Phi St=1)
double S5(const double& R, const double cparam, const double& Mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos);

// Filling factor du to growth
double PhiSGr(const double& sizef, const double& sizei, const double& phii,// ->
              const double& rhos, const double& eroll, const double& vrel, double&phipow);

// Filling factor due to collision [Okuzumi]
double PhiSColl(const double& R, const double& Mstar, const double& rhog, const double& cg, const double st,// ->
                const double& sizef, const double& eroll, const double& a0, const double& rhos, const double& alpha,// ->
                const int& iregime, double& phipow);


/* ------------------------ KATAOKA ------------------------*/

// Filling factor due to gas compression
double PhiSGas(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& deltav,// ->
               const double& sizef, const double& eroll, const double& a0, const int& iregime, double& phipow);

// Filling factor due to self-gravity
double PhiSGrav(const double& sizef, const double& rhos, const double& a0, const double& eroll);


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

// Return the min filling factor between PhiSColl, PhiSGas and PhiSgrav
double PhiMinSColGasGrav(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& deltav,// ->
                         const double st, const double& sizef, const double& rhos, const double& eroll, const double& a0,// ->
                         const double& alpha, const int& iregime, double& phipow);

// Return the min filling factor between phiSGr and PhiMinSColGasGrav
double PhiSFinal(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& deltav, const double st,// ->
                 const double& sizef, const double& sizei, const double& phii, const double& a0, const double& rhos, const double& eroll,// -> 
                 const double& alpha, const int& ifrag, const double& vfrag, const double& vrel, const int& iregime, double& phipow);

#endif // POROSITY_H_INCLUDED