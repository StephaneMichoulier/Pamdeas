#include <iostream>
#include <cmath>

#include "../header/general_functions.h"
#include "../header/disruption.h"

using namespace std;

double FreqSpin(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st, 
                const double& phi, const double& size, const double& gammaft)
{
    double deltav = DeltaV(R,Mstar,p,q,cg,st);
    /*cout << endl << "size: " << size << endl 
         << "phi: " << phi << endl 
         << "Freq: " << 2.*M_PI/(3600*5.*gammaft*phi*deltav/3./size)/60. << endl;*/
    return 3600*5.*gammaft*phi*deltav/3./size;
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
    double tensilestress = TensileStess(R,Mstar,p,q,cg,st,phi,size,rhos,gammaft);
    double maxtensiletress = 6e5*pow(phi,1.8)*esurf/0.1*(0.1e-6/a0);

    //cout << "tensilestress: " << tensilestress << endl;
    //cout << "maxtensiletress: " << maxtensiletress << endl;

    if (maxtensiletress <= tensilestress)
        return true;
    else
        return false;
}