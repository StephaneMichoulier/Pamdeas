#include <iostream>
#include <cmath>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"

using namespace std;


/* ------------------------- GRAIN MASS & SIZE ------------------------- */

double GrainMass(const double& size, const double& filfac, const double& rhos)
{   return 4.*M_PI*size*size*size*rhos*filfac/3.;  }

double GrainCrossSection(const double& size)
{   return M_PI*size*size;   }

double GrainMassToSize(const double& mass, const double& filfac, const double& rhos)
{   return pow(3.*mass/(4.*M_PI*rhos*filfac),1./3.);   }

double GrainVolumeMass(const double& mass, const double& filfac, const double& rhos)
{
    double size = GrainMassToSize(mass,filfac,rhos);
    return 4.*M_PI*size*size*size/3.;
}

double GrainVolumeSize(const double& size, const double& filfac, const double& rhos)
{   return 4.*M_PI*size*size*size/3.;   }


/* ------------------------- AERODYNAMICAL PARAMETERS ------------------------- */

double St(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& size,// ->
          const double& filfac, const double& rhos, const double& deltav, int& iregime)
{
    double st;

    if (TransRegEpSt(rhog,cg,size) < 1.)
    {   
        iregime = 1;
        st = rhos*filfac*size/(cg*rhog);
    }
    else
    {   
        double Re = 2.*size*deltav/NuMolGas(rhog,cg);
        if (Re < 1.)
        {
            iregime = 2;
            st = rhos*filfac*size*size/(4.5*NuMolGas(rhog,cg)*rhog);   
        }
        else if (Re < 800.)
        {
            iregime = 3;
            st = pow(Re,0.6)*rhos*filfac*size/(9.*rhog*deltav);
        }
        else
        {   
            iregime = 4;
            st = rhos*filfac*size/(0.165*rhog*deltav);    //0.165 = 1.32/8
        }
    }
    return st*Omegak(R,Mstar);
}


/* ------------------------ ENERGIES & VELOCITIES ------------------------ */

//double DeltaV2(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st)
//{   return 0.5*(p+0.5*q+1.5)*cg*cg/Vk(R,Mstar);   }

double Vrel(const double& cg, const double& st, const double& alpha)
{   return sqrt(2.*M_SQRT2*alpha*rossby*st)*cg/(1.+st); }

double Ekin(const double& mass, const double& vrel)
{   return 0.25*mass*vrel*vrel;     }

double Ekin(const double& size, const double& filfac, const double& rhos, const double& vrel)
{   return M_PI*rhos*filfac*size*size*size*vrel*vrel/3.;   }

double YoungMod(const double& filfac, const double& youngmod0)
{
    if (filfac >= 0.4)
    {   return youngmod0*pow(filfac,2.41);    }
    else
    {   return youngmod0*pow(0.4,2.41);    }
}

double Eroll(const double& a0, const double& esurf, const double& youngmod0)
{   
    //return 6.*M_PI*M_PI*esurf*a0*8.*1e-10;
    return 302.455974078*pow(pow(esurf,5.)*pow(a0,4.)/(youngmod0*youngmod0),1./3.);   
}

double Yd(const double& filfac, const double& filfaclim, const double& Yd0, const double& Ydpower)
{
    if (filfac >= filfaclim)
    {   return Yd0*pow(filfac,Ydpower); }
    else
    {   return Yd0*pow(filfaclim,Ydpower);  }
}

double Vstick(const double& size, const double& filfac, const double& rhos, const double& esurf, const double& youngmod0)
{
    double mass = GrainMass(size,filfac,rhos);
    return 4.23*pow(pow(esurf,5.)*pow(size,4.)/(mass*mass*mass*YoungMod(filfac,youngmod0)*YoungMod(filfac,youngmod0)),1./6.);
}

double Vyield(const double& size, const double& filfac, const double& rhos, const double& esurf, const double& youngmod0)
{   return 10.*Vstick(size,filfac,rhos,esurf,youngmod0);  }

double Vyield(const double& vstick)
{   return 10.*vstick;  }

double Vend(const double& size, const double& filfac, const double& rhos, const double& esurf, const double& youngmod0)
{   return 24343220.*Vstick(size,filfac,rhos,esurf,youngmod0);    }     // constant taken from the initial code PACED

double Vend(const double& vstick)
{   return 24343220.*vstick;   }     // constant taken from the initial code PACED

double Vfrag(const double& filfac, const double& filfaclim, const double& vfragi, const int& constvfrag)
{
    if (constvfrag == 1)
    {   return vfragi;  }
    else
    {                     // Need to look at the model for variable vfrag
        if (filfac>=0.6)
        {
            return vfragi*sqrt(pow(100.,1.5*(1-sqrt(filfac)))*pow(filfac,3.75*sqrt(filfac)-1.0));
        }
        else
        {
            if (filfac>=filfaclim)
            {
                return vfragi*sqrt(pow(100.,3*(0.5-(1/2.58)))*pow(filfac,4.92/2.58));
            }
            else
            {
                return vfragi*sqrt(pow(100.,3*(0.5-(1/2.58)))*pow(filfaclim,4.92/2.58));
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
