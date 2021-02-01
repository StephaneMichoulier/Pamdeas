#include <iostream>
#include <cmath>

#include "../header/constant.h"
#include "../header/general_functions.h"

using namespace std;


/* ------------------------- CONVERSIONS ------------------------- */

double AUtoMeter(const double& R)
{   return R*AU;    }

double MeterToAU(const double& R)
{   return R/AU;    }

double MsolToKg(const double& mass)
{   return mass*Msol;   }

double KgToMsol(const double& mass)
{   return mass/Msol;   }

double SecToYear(const double& time)
{   return time/365.25/3600./24.;   }

double YearToSec(const double& time)
{   return time*365.25*3600.*24.;   }


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


/* ------------------------- DISK QUANTITIES ------------------------- */

double Sigma0(const double& Rin, const double& Rout, const double& R0, const double& Mdisk, const double& p,
              const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight)
{
    double sigma0;

    if (p != 2.)
    {   sigma0 = AUtoMeter(R0)*AUtoMeter(R0)/(2.-p)*(pow(Rout/R0,2.-p)-pow(Rin/R0,2.-p));  }
    else
    {   sigma0 = AUtoMeter(R0)*AUtoMeter(R0)*log(Rout/Rin);    }

    if (ibump == 1)
    {
        double bumpin = (Rin-Rbump)/(M_SQRT2*bumpwidth);
        double bumpout = (Rout-Rbump)/(M_SQRT2*bumpwidth);

        sigma0 += bumpheight*AUtoMeter(bumpwidth)*(M_SQRT_PI/M_SQRT2*AUtoMeter(Rbump)*(erf(bumpout)-erf(bumpin))
        +AUtoMeter(bumpwidth)*(exp(-bumpin*bumpin)-exp(-bumpout*bumpout)));
    }

    return (Mdisk/2./M_PI)/sigma0;
}

double Sigma(const double& R, const double& p, const double& R0, const double& sigma0,
             const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight)
{
    double sigma = sigma0*pow(R/R0,-p);

    if (ibump == 1)
    {   sigma += bumpheight*sigma0*exp(-(R-Rbump)*(R-Rbump)/(2.*bumpwidth*bumpwidth));    }

    return sigma;
}

double Hg(const double& R, const double& q, const double& R0, const double& Hg0)
{   return Hg0*pow(R/R0,(3.-q)*0.5); }

double OmegaK(const double& R, const double& Mstar)
{   return sqrt(G*Mstar/AUtoMeter(R))/AUtoMeter(R);    }

double Vk(const double& R, const double& Mstar)
{   return OmegaK(R,Mstar)*AUtoMeter(R);    }

double Gravity(const double& R, double& Mstar)
{   return G*Mstar/(AUtoMeter(R)*AUtoMeter(R));    }

double Cg(const double& R, const double& Mstar, const double& Hg)
{   return AUtoMeter(Hg)*OmegaK(R,Mstar);   }

double Cg(const double& R, const double& Mstar, const double& q, const double& R0, const double& Hg0)
{   return AUtoMeter(Hg(R,q,R0,Hg0))*OmegaK(R,Mstar);   }

double T(const double& R, const double& q, const double& R0, const double& cg)
{   return cg*cg*Mgasmean/Kboltzmann;  }

double T(const double& R, const double& Mstar, const double& q, const double& R0, const double& Hg0)
{   return Cg(R,Mstar,q,R0,Hg0)*Mgasmean/Kboltzmann;  }

double Rhog(const double& sigma, const double& Hg)
{   return sigma/(AUtoMeter(Hg)*M_SQRT2*M_SQRT_PI);    }

double Rhog(const double& R, const double& p, const double& q, const double& sigma0, const double& R0,//->
            const double& Hg0,const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight)
{
    return Sigma(R,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight)/(AUtoMeter(Hg(R,q,R0,Hg0))*M_SQRT2*M_SQRT_PI);
}

double Pg(const double& rhog, const double& cg)
{   return rhog*cg*cg;    }

double Pg(const double& R, const double& Mstar, const double& p, const double& q,
          const double& sigma0, const double& R0, const double& Hg0,const int& ibump, const double& Rbump,//->
          const double& bumpwidth, const double& bumpheight)
{
    double cg = Cg(R,Mstar,q,R0,Hg0);
    return Rhog(R,p,q,sigma0,R0,Hg0,ibump,Rbump,bumpwidth,bumpheight)*cg*cg;
}

double NuMolGas(const double& rhog, const double& cg)
{   return 5.*M_SQRT_PI*Mgasmean*cg/(64.*Sigmamol*rhog);   }

double NuTurbGas(const double& R, const double& Mstar, const double& alpha, const double& cg)
{   return alpha*cg*cg/OmegaK(R,Mstar);   }

double NuTurbGas(const double& R, const double& Mstar, const double& q, const double& R0,//->
                 const double& Hg0, const double& alpha)
{
    double cg = Cg(R,Mstar,q,R0,Hg0);
    return alpha*cg*cg/OmegaK(R,Mstar);
}

double Lambda(const double& rhog, const double& cg)
{   return 2.*NuMolGas(rhog,cg)/cg;     }

double TransRegEpSt(const double& rhog, const double& cg, const double& size)
{   return size/(2.25*Lambda(rhog,cg));  }


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

double DustFrac(const double& dustfrac0, const double& dustfracmax, const double& R, const double& Rbump,
                const double& bumpwidth, const int& ibump)
{
    if (ibump == 0)
    {   return dustfrac0;   }
    else
    {   return dustfrac0 + ibump*(dustfracmax-dustfrac0)*exp(-(R-Rbump)*(R-Rbump)/(2.*bumpwidth*bumpwidth));    }
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
{   //return 302.455974078*pow(pow(esurf,5.)*pow(a0,4.)/(youngmod0*youngmod0),1./3.);   
    return 6.*M_PI*M_PI*esurf*a0*8.e-10;
}

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
