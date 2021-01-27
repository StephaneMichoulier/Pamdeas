#include <iostream>
#include <cmath>

#include "../header/general_functions.h"
#include "../header/disruption.h"

using namespace std;

double FreqSpin(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st,
                const double& phi, const double& size, const double& gammaft)
{
    double deltav = st*st/(1.+st*st)*DeltaV(R,Mstar,p,q,cg,st);
    return 5.*gammaft*deltav/3./size;
}

double TensileStess(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st,
                    const double& phi, const double& size, const double& rhos, const double& gammaft)
{
    double freqspin = FreqSpin(R,Mstar,p,q,cg,st,phi,size,gammaft);
    return rhos*phi*size*size*freqspin*freqspin/4.;
}

bool Disrupt(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st,
             const double& phi, const double& size, const double& rhos, const double& gammaft, const double& esurf,
             const double& a0)
{
    double freqspin = FreqSpin(R,Mstar,p,q,cg,st,phi,size,gammaft);
    double tensilestress = TensileStess(R,Mstar,p,q,cg,st,phi,size,rhos,gammaft);
    double maxtensiletress = 0.6*pow(phi,1.8)*esurf/a0;

    if (maxtensiletress <= tensilestress)
    {    cout << endl
             << "mass: " << 4.*M_PI/3.*rhos*phi*size*size*size*1000. << endl
             << "phi: " << phi*10000. << endl
             << "size: " << size << endl
             << "freq: " << 2.*M_PI/freqspin/60. << endl
             << "tensilestress: "<< tensilestress << endl
             << "maxtensiletress: " << maxtensiletress << endl;
             return true;}
    else
        return false;
}
