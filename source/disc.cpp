#include <iostream>
#include <cmath>

#include "../header/constantandconversion.h"
#include "../header/disc.h"

using namespace std;


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

double DustFrac(const double& dustfrac0, const double& dustfracmax, const double& R, const double& Rbump,
                const double& bumpwidth, const int& ibump)
{
    if (ibump == 0)
    {   return dustfrac0;   }
    else
    {   return dustfrac0 + ibump*(dustfracmax-dustfrac0)*exp(-(R-Rbump)*(R-Rbump)/(2.*bumpwidth*bumpwidth));    }
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
