#include <iostream>
#include <cmath>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"

using namespace std;


/* ------------------------- GRAIN MASS & SIZE ------------------------- */

double GrainMass(const double& size, const double& phi, const double& rhos)
{   return 4.*M_PI*size*size*size*rhos*phi/3.;  }

double GrainCrossSection(const double& size)
{   return M_PI*size*size;   }

double GrainMassToSize(const double& mass, const double& phi, const double& rhos)
{   return pow(3.*mass/(4.*M_PI*rhos*phi),1./3.);   }

double GrainVolumeMass(const double& mass, const double& phi, const double& rhos)
{
    double size = GrainMassToSize(mass,phi,rhos);
    return 4.*M_PI*size*size*size/3.;
}

double GrainVolumeSize(const double& size, const double& phi, const double& rhos)
{   return 4.*M_PI*size*size*size/3.;   }


/* ------------------------- AERODYNAMICAL PARAMETERS ------------------------- */

double St(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& size,// ->
          const double& phi, const double& rhos, int& iregime)
{
    double st;

    if (TransRegEpSt(rhog,cg,size) < 1.)
    {
        iregime = 1;
        st = rhos*phi*size*OmegaK(R,Mstar)/(cg*rhog);

        if (st >= 1.)
        {   iregime = 3;   }
    }
    else
    {
        iregime = 2;
        st = rhos*phi*size*size*OmegaK(R,Mstar)/(4.5*NuMolGas(rhog,cg)*rhog);

        if (st >= 1.)
        {   iregime = 4;    }
    }

    return st;
}


/* ------------------------ ENERGIES & VELOCITIES ------------------------ */

double DeltaV2(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st)
{   return 0.5*(p+0.5*q+1.5)*cg*cg/Vk(R,Mstar);   }

double Vrel(const double& cg, const double& st, const double& alpha)
{   return sqrt(2.*M_SQRT2*alpha*Rossby*st)*cg/(1.+st); }

double Ekin(const double& mass, const double& vrel)
{   return 0.25*mass*vrel*vrel;     }

double Ekin(const double& size, const double& phi, const double& rhos, const double& vrel)
{   return M_PI*rhos*phi*size*size*size*vrel*vrel/3.;   }

double YoungMod(const double& phi, const double& youngmod0)
{
    if (phi >= 0.4)
    {   return youngmod0*pow(phi,2.41);    }
    else
    {   return youngmod0*pow(0.4,2.41);    }
}

double Eroll(const double& a0, const double& esurf, const double& youngmod0)
{   return 302.455974078*pow(pow(esurf,5.)*pow(a0,4.)/(youngmod0*youngmod0),1./3.);   }

double Yd(const double& phi, const double& philim, const double& Yd0, const double& Ydpower)
{
    if (phi >= philim)
    {   return Yd0*pow(phi,Ydpower); }
    else
    {   return Yd0*pow(philim,Ydpower);  }
}

double Vstick(const double& size, const double& phi, const double& rhos, const double& esurf, const double& youngmod0)
{
    double mass = GrainMass(size,phi,rhos);
    return 4.23*pow(pow(esurf,5.)*pow(size,4.)/(mass*mass*mass*YoungMod(phi,youngmod0)*YoungMod(phi,youngmod0)),1./6.);
}

double Vyield(const double& size, const double& phi, const double& rhos, const double& esurf, const double& youngmod0)
{   return 10.*Vstick(size,phi,rhos,esurf,youngmod0);  }

double Vyield(const double& vstick)
{   return 10.*vstick;  }

double Vend(const double& size, const double& phi, const double& rhos, const double& esurf, const double& youngmod0)
{   return 24343220.*Vstick(size,phi,rhos,esurf,youngmod0);    }     // constant taken from the initial code PACED

double Vend(const double& vstick)
{   return 24343220.*vstick;   }     // constant taken from the initial code PACED

double Vfrag(const double& phi, const double& philim, const double& vfragi, const int& constvfrag)
{
    if (constvfrag == 1)
    {   return vfragi;  }
    else
    {                     // Need to look at the model for variable vfrag
        if (phi>=0.6)
        {
            return vfragi*sqrt(pow(100.0,1.5*(1-sqrt(phi)))*pow(phi,3.75*sqrt(phi)-1.0));
        }
        else
        {
            if (phi>=philim)
            {
                return vfragi*sqrt(pow(100.0,3*(0.5-(1/2.58)))*pow(phi,4.92/2.58));
            }
            else
            {
                return vfragi*sqrt(pow(100.0,3*(0.5-(1/2.58)))*pow(philim,4.92/2.58));
            }
        }
    }
}

double CoeffRest(const double& vrel, const double& vstick, const double& vyield)
{
    if (vrel <= vstick)
    {   return 0.;    }
    else
    {
        if (vrel <= vyield)
        {   return sqrt(1.-(vstick*vstick/(vrel*vrel)));  }
        else
        {
            double vyieldonvrel=vyield/vrel;
            double vstickonvrel=vstick/vrel;
            return sqrt(1.2*sqrt(3.)*(1.-(vyieldonvrel*vyieldonvrel/6.))*
                   sqrt(1./(1.+2.*sqrt((1.2/(vyieldonvrel*vyieldonvrel))-0.2)))-(vstickonvrel*vstickonvrel));
        }
    }
}