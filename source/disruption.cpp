#include <iostream>
#include <cmath>
#include <random>

#include "../header/disruption.h"
#include "../header/dust.h"
#include "../header/porosity.h"

using namespace std;


/* ------------------------- DISRUPTION ------------------------- */

void DisplayDisruptParam(const double& size, const double& filfac, const double& rhos,const double& freqspin, const double& tensilestress)
{
    cout << endl;
    cout << "mass (g): " << 4.*M_PI/3.*rhos*filfac*size*size*size*1000. << endl;
    cout << "filfac: " << filfac << endl;
    cout << "size (m): " << size << endl;
    cout << "2pi/wc (min): " << 2.*M_PI/freqspin/60. << endl;
    cout << "tensilestress (Pa): " << tensilestress << endl << endl;
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

bool Disruptlim(const double& size, const double& filfac, const double& rhos, const double& deltav, const double& gammaft,// -> 
                const double& weirmod, const double& esurf, const double& a0, const int& disrupteq)
{
    double freqspin = FreqSpin(size,deltav,gammaft);
    double tensilestress = TensileStess(size,filfac,rhos,freqspin);
    double maxtensilestress;

    if (disrupteq == 0)
    {   // Tatsuuma et al. 2019
        maxtensilestress = 0.6*pow(filfac,1.8)*esurf/a0;
    }
    else
    {   // Kimura et al. 2020
        double vol = GrainVolumeSize(size,filfac,rhos);
        maxtensilestress = 8000*(esurf/0.1)*pow(a0*1e6/0.1,3./weirmod-1)*pow(filfac/0.1,1.5-1./weirmod)*exp(0.24*(filfac/0.1 - 1.))*pow(vol*1e18/686,-1./weirmod);
    }

    if (maxtensilestress <= tensilestress)
    {   //display values when grain is disrupted
        //DisplayDisruptParam(size,filfac,rhos,freqspin,tensilestress);
        return true;
    }
    else return false;
}

void Disrupt(double& mass, double& filfac, const double& sizeini, const double& rhos,const double& R, const double& mstar, // ->
             const double& rhog, const double cg, const double& deltav, const double st, const double& eroll, // ->
             const double& a0, const double& alpha, int& porreg)
{
    double maxmass = log10(mass/2.);
    double minmass = log10(GrainMass(sizeini,filfac,rhos));
    double curmass = log10(mass);
    double randmass;
    if (maxmass > minmass)
    {   randmass =((maxmass- minmass) * ((float)rand() / RAND_MAX)) + minmass;}
    else
    {   if (curmass > minmass)  randmass = minmass;
        else randmass = curmass;
    }
    mass = pow(10.,randmass);
    if (mass < GrainMass(a0,1.,rhos)) mass=GrainMass(a0,1.,rhos);
    double filfacmin = FilFacMinMColGasGrav(R,mstar,rhog,cg,deltav,st,mass,rhos,eroll,a0,GrainMass(a0,1.,rhos),alpha,porreg);
    if (filfacmin > filfac)    filfac = filfacmin;
}
