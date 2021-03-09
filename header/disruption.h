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
bool Disrupt(const double& size, const double& filfac, const double& rhos, const double& deltav, const double& gammaft,// -> 
             const double& esurf, const double& a0);
             
#endif // DISRUPTION_H_INCLUDED
