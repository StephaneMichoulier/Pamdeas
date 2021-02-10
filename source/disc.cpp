#include <iostream>
#include <cmath>

#include "../header/constantandconversion.h"
#include "../header/disc.h"

using namespace std;


/* ------------------------- DISK QUANTITIES ------------------------- */

double Sigma0(const double& Rin, const double& Rout, const double& R0, const double& mdisk, const double& p,
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

    return (mdisk/2./M_PI)/sigma0;
}

double Sigma(const double& R, const double& p, const double& R0, const double& sigma0,
             const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight)
{
    double sigma = sigma0*pow(R/R0,-p);

    if (ibump == 1)
    {   sigma += bumpheight*sigma0*exp(-(R-Rbump)*(R-Rbump)/(2.*bumpwidth*bumpwidth));    }

    return sigma;
}

double Hg(const double& R, const double& q, const double& R0, const double& hg0)
{   return hg0*pow(R/R0,(3.-q)*0.5); }

double Omegak(const double& R, const double& mstar)
{   return sqrt(G*mstar/AUtoMeter(R))/AUtoMeter(R);    }

double Vk(const double& R, const double& mstar)
{   return Omegak(R,mstar)*AUtoMeter(R);    }

double Gravity(const double& R, double& mstar)
{   return G*mstar/(AUtoMeter(R)*AUtoMeter(R));    }

double Cg(const double& R, const double& mstar, const double& hg)
{   return AUtoMeter(hg)*Omegak(R,mstar);   }

double Cg(const double& R, const double& mstar, const double& q, const double& R0, const double& hg0)
{   return AUtoMeter(Hg(R,q,R0,hg0))*Omegak(R,mstar);   }

double T(const double& R, const double& q, const double& R0, const double& cg)
{   return cg*cg*mgasmean/kboltzmann;  }

double T(const double& R, const double& mstar, const double& q, const double& R0, const double& hg0)
{   return Cg(R,mstar,q,R0,hg0)*mgasmean/kboltzmann;  }

double Rhog(const double& sigma, const double& hg)
{   return sigma/(AUtoMeter(hg)*M_SQRT2*M_SQRT_PI);    }

double Rhog(const double& R, const double& p, const double& q, const double& sigma0, const double& R0,//->
            const double& hg0,const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight)
{
    return Sigma(R,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight)/(AUtoMeter(Hg(R,q,R0,hg0))*M_SQRT2*M_SQRT_PI);
}

double Pg(const double& rhog, const double& cg)
{   return rhog*cg*cg;    }

double Pg(const double& R, const double& mstar, const double& p, const double& q,
          const double& sigma0, const double& R0, const double& hg0,const int& ibump, const double& Rbump,//->
          const double& bumpwidth, const double& bumpheight)
{
    double cg = Cg(R,mstar,q,R0,hg0);
    return Rhog(R,p,q,sigma0,R0,hg0,ibump,Rbump,bumpwidth,bumpheight)*cg*cg;
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
{   return 5.*M_SQRT_PI*mgasmean*cg/(64.*sigmamol*rhog);   }

double NuTurbGas(const double& R, const double& mstar, const double& alpha, const double& cg)
{   return alpha*cg*cg/Omegak(R,mstar);   }

double NuTurbGas(const double& R, const double& mstar, const double& q, const double& R0,//->
                 const double& hg0, const double& alpha)
{
    double cg = Cg(R,mstar,q,R0,hg0);
    return alpha*cg*cg/Omegak(R,mstar);
}

double Lambda(const double& rhog, const double& cg)
{   return 2.*NuMolGas(rhog,cg)/cg;     }

double TransRegEpSt(const double& rhog, const double& cg, const double& size)
{   return size/(2.25*Lambda(rhog,cg));  }