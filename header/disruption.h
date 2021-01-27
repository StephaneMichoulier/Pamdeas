#ifndef DISRUPTION_H_INCLUDED
#define DISRUPTION_H_INCLUDED

#include <iostream>

using namespace std;

double FreqSpin(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st, 
                const double& phi, const double& size, const double& gammaft);

double TensileStess(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st,
                    const double& phi, const double& size, const double& rhos, const double& gammaft);

bool Disrupt(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st,
             const double& phi, const double& size, const double& rhos, const double& gammaft, const double& esurf,
             const double& a0);

#endif // DISRUPTION_H_INCLUDED
