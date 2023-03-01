#ifndef DISRUPTION_H_INCLUDED
#define DISRUPTION_H_INCLUDED

#include <iostream>

using namespace std;


/* ------------------------- DISPLAY ------------------------- */

// Display disrupt param of a grain
void DisplayDisruptParam(const double& size, const double& filfac, const double& rhos,const double& freqspin, const double& tensilestress);


/* ------------------------- DISRUPTION ------------------------- */

// Compute the steady-state angular velocity [rad/s]
double FreqSpin(const double& size, const double& deltav, const double& gammaft);

// Compute the tensile stress due to the centrifugal force [Pa]
double TensileStess(const double& size, const double& filfac, const double& rhos, const double& freqspin);
double TensileStess(const double& size, const double& filfac, const double& rhos, const double& deltav, const double& gammaft);

// Compare tensile stress and tensile strength 
bool Disruptlim(const double& size, const double& filfac, const double& rhos, const double& deltav, const double& gammaft,//->
                const double& weirmod, const double& esurf, const double& a0, const int& disrupteq);

// Compute the mass & filling factor of the disrupted grain
void Disrupt(double& mass, double& filfac, const double& sizeini, const double& rhos,const double& R, const double& mstar, //->
             const double& rhog, const double cg, const double& deltav, const double st, const double& eroll, //->
             const double& a0, const double& alpha, int& porreg);
             
#endif // DISRUPTION_H_INCLUDED
