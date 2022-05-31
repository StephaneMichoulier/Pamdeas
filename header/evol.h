#ifndef EVOL_H_INCLUDED
#define EVOL_H_INCLUDED

#include <iostream>

using namespace std;


/* ------------------------ TIME STEP ------------------------ */

// Compute a dt as a fraction of orbital time at radius R
double KeplerDt(double& omegakfraction, const double omegak);

// Compute an adaptative dt to ensure stability
double AdaptativeDt(const double& time, const double& timeend, const int& massorsize, const double& vargrowth,// ->
                    const double& R, const double& dvargrowthdt, const double& dRdt);


/* ------------------------ RADIAL DRIFT ------------------------*/

double DRDt(const double& R, const double& mstar, double p, double q, const double& rhog, const double& cg, const double& R0,// -> 
            const double& sigma0, const double& hg0, const double& dustfrac, const double& st, const double& alpha, const int& ibr,// ->
            const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight);


/* ------------------------  GROWTH-FRAG-BOUNCE dm/dt ------------------------ */

// Compute the final dm/dt between growth, frag, bounce [kg/s]
double DmDt(const double& size, const double& rhog, const double& dustfrac, const double& vrel, const int& ifrag,// ->
            const int& ieros, const double& veros, const int& ibounce, const double& vfrag, const double& vstick,// ->
            const double& probabounce);


/* ------------------------  GROWTH-FRAG ds/dt ------------------------ */

// Compute the final ds/dt between growth, frag [kg/s], Bounce not included !
double DsDt(const double& filfac, const double& rhog, const double& rhos, const double& dustfrac,// ->
            const double& vrel, const int& ifrag, const double& vfrag, const double& filfacpow);


/* ------------------------ CHARACTERISTIC TIMES ------------------------*/

// Mean drift time [s]
double Tdrift(const double& R, const double& drdt);

// Mean growth time [s]
double Tgrowth(const double& var, const double& dvardt);

// Compute the mean collision time [s]
double Tcoll(const double& size, const double& rhog, const double& filfac, const double& rhos,// ->
             const double& dustfrac, const double& vrel);

// Compute the number of collision during dt
double Ncoll(const double& dt, const double& tcoll);

#endif // EVOL_H_INCLUDED
