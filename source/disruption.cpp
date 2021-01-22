#include <iostream>
#include <cmath>

#include "../header/general_functions.h"

using namespace std;

double FreqSpin(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st, 
                const double& phi, const double& size, const double& gammaft)
{
    return gammaft*phi*DeltaV(R,Mstar,p,q,cg,st)/9./size;
}

double TensileStess(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st, 
                    const double& phi, const double& size, const double& rhos, const double& gammaft)
{
    double freqspin = FreqSpin(gammaft,phi,size);
    return rhos*size*size*freqspin*freqspin/4.
}

bool Disrupt(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st, 
             const double& phi, const double& size, const double& rhos, const double& gammaft, const double& esurf, 
             const double& a0)
{
    double tensilestress = TensileStess(gammaft,phi,size,rhos);
    double maxtensiletress = 0.6*pow(phi,1.8)*esurf/a0;

    if (maxtensiletress <= tensilestress)
        return true;
    else
        return false;
}