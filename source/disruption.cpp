#include <iostream>
#include <cmath>

#include "../header/general_functions.h"
#include "../header/disruption.h"

using namespace std;


/* ------------------------- DISRUPTION ------------------------- */

double FreqSpin(const double& size, const double& deltav, const double& gammaft)
{
    return 5.*gammaft*deltav/3./size;
}

double TensileStess(const double& size, const double& phi, const double& rhos, const double& deltav, const double& gammaft)
{
    double freqspin = FreqSpin(size,deltav,gammaft);
    return rhos*phi*size*size*freqspin*freqspin/4.;
}

bool Disrupt(const double& size, const double& phi, const double& rhos, const double& deltav, const double& gammaft,// -> 
             const double& esurf, const double& a0)
{
    double freqspin = FreqSpin(size,deltav,gammaft);
    double tensilestress = TensileStess(size,phi,rhos,deltav,gammaft);
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
