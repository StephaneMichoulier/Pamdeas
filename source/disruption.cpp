#include <iostream>
#include <cmath>

#include "../header/disruption.h"

using namespace std;


/* ------------------------- DISRUPTION ------------------------- */

void DisplayDisruptParam(const double& size, const double& filfac, const double& rhos,const double& freqspin, const double& tensilestress)
{
    cout << endl;
    cout << "mass (g): " << 4.*M_PI/3.*rhos*filfac*size*size*size*1000. << endl;
    cout << "filfac: " << filfac << endl;
    cout << "size (m): " << size << endl;
    cout << "2pi/wc (min): " << 2.*M_PI/freqspin/60. << endl;
    cout << "tensilestress (Pa): "<< tensilestress << endl;
}


/* ------------------------- DISRUPTION ------------------------- */

double FreqSpin(const double& size, const double& deltav, const double& gammaft)
{   return 5.*gammaft*deltav/3./size;   }

double TensileStess(const double& size, const double& filfac, const double& rhos, const double& freqspin)
{   return rhos*filfac*size*size*freqspin*freqspin/4.;  }

double TensileStess(const double& size, const double& filfac, const double& rhos, const double& deltav, const double& gammaft)
{
    double freqspin = FreqSpin(size,deltav,gammaft);
    return rhos*filfac*size*size*freqspin*freqspin/4.;
}

bool Disrupt(const double& size, const double& filfac, const double& rhos, const double& deltav, const double& gammaft,// -> 
             const double& esurf, const double& a0)
{
    double freqspin = FreqSpin(size,deltav,gammaft);
    double tensilestress = TensileStess(size,filfac,rhos,freqspin);
    double maxtensilestress = 0.6*pow(filfac,1.8)*esurf/a0;

    if (maxtensilestress <= tensilestress)
    {   //display values when grain is disrupted
        DisplayDisruptParam(size,filfac,rhos,freqspin,tensilestress);
        return true;
    }
    else return false;
}
