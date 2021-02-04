#include <iostream>
#include <cmath>

#include "../header/disruption.h"

using namespace std;


/* ------------------------- DISRUPTION ------------------------- */

double FreqSpin(const double& size, const double& deltav, const double& gammaft)
{   return 5.*gammaft*deltav/3./size;   }

double TensileStess(const double& size, const double& phi, const double& rhos, const double& deltav, const double& gammaft)
{
    double freqspin = FreqSpin(size,deltav,gammaft);
    return rhos*phi*size*size*freqspin*freqspin/4.;
}

bool Disrupt(const double& size, const double& phi, const double& rhos, const double& deltav, const double& gammaft,// -> 
             const double& esurf, const double& a0, const double& st)
{
    double freqspin = FreqSpin(size,deltav,gammaft);
    double tensilestress = TensileStess(size,phi,rhos,deltav,gammaft);
    double maxtensilestress = 0.6*pow(phi,1.8)*esurf/a0;

    if (maxtensilestress <= tensilestress)
    {    
        //display values
        cout << endl
             << "mass (g): " << 4.*M_PI/3.*rhos*phi*size*size*size*1000. << endl
             << "phi*1e4 (phi*1e4): " << phi*10000. << endl
             << "size (m): " << size << endl
             << "St: " << st << endl
             << "wc (r/s): " << freqspin << endl
             << "2pi/wc (min): " << 2.*M_PI/freqspin/60. << endl
             << "a*wc (m/s): " << size*freqspin << endl
             << "tensilestress (Pa): "<< tensilestress << endl
             << "maxtensiletress (Pa): " << maxtensilestress << endl;

        return true;}
    else
        return false;
}
