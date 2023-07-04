#ifndef POROSITY_H_INCLUDED
#define POROSITY_H_INCLUDED

#include <iostream>

using namespace std;


/* ------------------------ SHARED FUNCTIONS ------------------------ */

// A parameter shared by all transition masses and sizes
double CParam(const double& R, const double& mstar, const double& rhog, const double& cg, const double& eroll,//->
              const double& a0, const double& rhos, const double& alpha);
                 
/* ------------------------ BOUNCE ------------------------*/

// Compute the probability to bounce
double ProbaBounce(const double& filfac, const double& filfacbnc, const double& vrel, const double& vstick, const double& vend);

// Compute the variation in volume after bouncing [m³]
//double VarVolumeBounce(const double& filfac, const double& filfaclim, const double& coeffrest, const double& ekin,//->
//                       const double& volume, const double& Yd0, const double& Ydpower);

// Compute the filling factor after bouncing
double FilFacBounce(const double& sizei, const double& filfaci, const double& rhos,//->
                    const double& vrel, const double& vstick, const double& vyield, const double& ncoll,//->
                    const double& eroll,const double& a0);

/* ------------------------ MASS MODEL------------------------ */

/* ------------------------ OKUZUMI GROWTH MODEL ------------------------ */
              
// Transition mass M1/m0 (FilFac h&s = FilFac Ep-St<1)
double M1(const double& cparam);

// Transition mass M2/m0 (FilFac h&s = FilFac St-St<1)
double M2(const double& cparam, const double& rhog, const double& cg, const double& a0);

// Transition mass M3/m0 (FilFac Ep-St<1 = FilFac St-St<1)
double M3(const double& m1_on_m0, const double& m2_on_m0);
double M3(const double& cparam, const double& rhog, const double& cg, const double& a0);

// Transition mass M4/m0 (FilFac Ep-St<1 = FilFac St=1)
double M4(const double& R, const double cparam, const double& mstar, const double& rhog,//->
          const double& cg, const double& a0, const double& rhos);

// Transition mass M5/m0 (FilFac St-St<1 = FilFac St=1)
double M5(const double& R, const double cparam, const double& mstar, const double& rhog,//->
          const double& cg, const double& a0, const double& rhos);

// Filling factor du to growth
double FilFacMGr(const double& massf, const double& massi, const double& filfaci,//->
                 const double& rhos, const double& eroll, const double& vrel);

// Filling factor due to collision, mfrac is massf/m0 [Okuzumi]
double FilFacMColl(const double& R, const double& mstar, const double& rhog, const double& cg, const double& st,//->
                   const double& mfrac, const double eroll, const double& a0, const double& rhos, const double& alpha,//->
                   int& porreg);


/* ------------------------ KATAOKA ------------------------*/

// Filling factor due to gas compression
double FilFacMGas(const double& R, const double& mstar, const double& deltav, const double& st, const double& massf, const double& rhos,//-> 
                  const double& eroll, const double& a0, const double& m0);//, const int& dragreg, const double& rhog, const double& cg);

// Filling factor due to self-gravity
double FilFacMGrav(const double& massf, const double& rhos, const double& eroll, const double& a0, const double& m0);


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

// Return the minimum filling factor between FilFacMColl, FilFacMGas and FilFacMgrav
double FilFacMinMColGasGrav(const double& R, const double& mstar, const double& rhog, const double cg, const double& deltav,//->
                            const double st, const double& massf, const double& rhos, const double& eroll, const double& a0,//->
                            const double& m0, const double& alpha, int& porreg);

// Return the minimum filling factor between filfacMGr, FilFacMinMColGasGrav and FilFacMBounce
double FilFacMFinal(const double& R, const double& mstar, const double& rhog, const double& cg, const double deltav, const double st,//->
                    const double& massf, const double& massi, const double& filfaci, const double& a0, const double& rhos, const double& eroll,//->
                    const double& alpha, const double& ncoll, const int& ifrag, const int& ibounce, const double& vfrag, const int& icomp,//->
                    const double& vrel, const double& vstick, const double& probabounce, const double& filfaclim, const double& Yd0,//-> 
                    const double& Ydpower, int& porreg);


/* ------------------------ SIZE MODEL ------------------------ */

/* ------------------------ OKUZUMI GROWTH MODEL ------------------------ */

// Transition size S1/a0 (FilFac h&s = FilFac Ep-St<1)
double S1(const double& cparam);

// Transition size S2/a0 (FilFac h&s = FilFac St-St<1)
double S2(const double& cparam, const double& rhog, const double& cg, const double& a0);

// Transition size S3/a0 (FilFac Ep-St<1 = FilFac St-St<1)
double S3(const double& rhog, const double& cg, const double& a0);

// Transition size S4/a0 (FilFac Ep-St<1 = FilFac St=1)
double S4(const double& R, const double cparam, const double& mstar, const double& rhog,//->
          const double& cg, const double& a0, const double& rhos);

// Transition size S5/a0 (FilFac St-St<1 = FilFac St=1)
double S5(const double& R, const double cparam, const double& mstar, const double& rhog,//->
          const double& cg, const double& a0, const double& rhos);

// Filling factor du to growth
double FilFacSGr(const double& sizef, const double& sizei, const double& filfaci,//->
              const double& rhos, const double& eroll, const double& vrel, double&filfacpow);

// Filling factor due to collision [Okuzumi]
double FilFacSColl(const double& R, const double& mstar, const double& rhog, const double& cg, const double st, const double& sfrac,//->
                   const double& eroll, const double& a0, const double& rhos, const double& alpha, int& porreg, double& filfacpow);


/* ------------------------ KATAOKA ------------------------*/

// Filling factor due to gas compression
double FilFacSGas(const double& R, const double& mstar, const double& deltav, const double& st,const double& sizef, const double& eroll,//->
                  const double& a0, const double& rhos, const int& dragreg, double& filfacpow);//, const double& rhog, const double& cg);

// Filling factor due to self-gravity
double FilFacSGrav(const double& sizef, const double& rhos, const double& a0, const double& eroll);


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

// Return the min filling factor between FilFacSColl, FilFacSGas and FilFacSgrav
double FilFacMinSColGasGrav(const double& R, const double& mstar, const double& rhog, const double& cg, const double& deltav,//->
                            const double st, const double& sizef, const double& rhos, const double& eroll, const double& a0,//->
                            const double& alpha, const int& dragreg, int& porreg, double& filfacpow);

// Return the min filling factor between filfacSGr and FilFacMinSColGasGrav
double FilFacSFinal(const double& R, const double& mstar, const double& rhog, const double& cg, const double& deltav, const double st,//->
                    const double& sizef, const double& sizei, const double& filfaci, const double& a0, const double& rhos, const double& eroll,//-> 
                    const double& alpha, const int& ifrag, const double& vfrag, const double& vrel, const int& dragreg, int& porreg, double& filfacpow);

#endif // POROSITY_H_INCLUDED