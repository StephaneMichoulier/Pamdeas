#include <iostream>
#include <cmath>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/evol.h"

using namespace std;

/* ------------------------  TIME STEP ------------------------ */

double KeplerDt(double& omegakfraction, const double omegak)
{
    return SecToYear(omegakfraction*2.*M_PI/omegak);
}

double AdaptativeDt(const double& time, const double& timeend, const int& massorsize, const double& vargrowth,// -> 
                    const double& R, const double& dvargrowthdt, const double& dRdt)
{
    double dt = 1.;
    double limup = 0.15;
    double limdown = 0.075;
    double CFLcdt;
    double Cgrowth = abs(vargrowth/dvargrowthdt);    // CFL from growth rate
    double Cdrift = abs(R/dRdt);                     // CFL from radial drift rate

    if (Cgrowth != 0 && Cdrift != 0)
    {
        if (Cgrowth < Cdrift)       CFLcdt = Cgrowth;
        else                        CFLcdt = Cdrift;
    }
    else if (Cgrowth != 0 && Cdrift == 0)       CFLcdt = Cgrowth;
    else if (Cgrowth == 0 && Cdrift != 0)       CFLcdt = Cdrift;
    else
    {
        cerr << "Error: time step iteration impossible" << endl;
        exit(1);
    }

    if (massorsize == 1)
    {   limup /= 3.;   limdown /= 3.;   }

    CFLcdt = SecToYear(CFLcdt);
    limup *= CFLcdt;
    limdown *= CFLcdt;
    do
    {
        if (dt > limup)     dt /= 2.;
        if (dt < limdown)   dt *= 2.;

    } while (limup < dt || limdown >  dt);

    if (time+dt > timeend)
    {   dt = timeend - time;   }

    return dt;
}


/* ------------------------ RADIAL DRIFT ------------------------*/

double DRDt(const double& R, const double& mstar, double p, double q, const double& rhog, const double& cg, const double& R0,// -> 
            const double& sigma0, const double& hg0, const double& dustfrac, const double& st, const double& alpha, const int& ibr,// ->
            const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight)
{
    double dfbr = 1.;   //dust fraction back reaction
    if (ibr == 1)   dfbr += dustfrac;

    double vdrift = VDrift(R,mstar,p,q,rhog,cg,R0,sigma0,hg0,ibump,Rbump,bumpwidth,bumpheight);
    double vvisc = VVisc(R,mstar,p,q,rhog,cg,R0,sigma0,hg0,alpha,ibump,Rbump,bumpwidth,bumpheight);

    vdrift *= st/(dfbr*dfbr+st*st);
    vvisc *= dfbr/(dfbr*dfbr+st*st);

    return MeterToAU(vdrift+vvisc);
}


/* ------------------------  GROWTH-FRAG-BOUNCE dm/dt ------------------------ */

double DmDt(const double& size, const double& rhog, const double& dustfrac, const double& vrel, const int& ifrag,// ->
            const int& ibounce, const double& vfrag, const double& vstick, const double& probabounce)
{
    double dmdt = 4.*M_PI*dustfrac*rhog*size*size*vrel;

    switch (ibounce)
    {
        case (0):
        {
            if (vrel >= vfrag && ifrag > 0) switch (ifrag)
            {
                case (1):
                {
                    dmdt *= -1.;
                    break;
                }
                case (2):
                {
                    dmdt *= -vrel*vrel/(vrel*vrel+vfrag*vfrag);
                    break;
                }
            }
            break;
        }
        case (1):
        {
            if (vrel >= vfrag && ifrag > 0)    switch (ifrag)
            {   
                case (1):
                {
                    dmdt *= -1.;
                    break;
                }
                case (2):
                {
                    dmdt *= -vrel*vrel/(vrel*vrel+vfrag*vfrag);
                    break;
                }
            }
            else
            {
                if (vrel > vstick)
                    {
                        dmdt *= probabounce;
                        break;
                    }
            }
            break;
        }
    }
    return dmdt;
}


/* ------------------------  GROWTH-FRAG-BOUNCE ds/dt ------------------------ */

double DsDt(const double& filfac, const double& rhog, const double& rhos, const double& dustfrac,// ->
            const double& vrel, const int& ifrag, const double& vfrag, const double& filfacpow)
{
    double dsdt = dustfrac*rhog*vrel/(filfac*rhos*(1.+filfacpow/3.));

    if (vrel >= vfrag && ifrag > 0) switch (ifrag)
    {
        case (1):
        {
            dsdt *= -1.;
            break;
        }
        case (2):
        {
            dsdt *= -vrel*vrel/(vrel*vrel+vfrag*vfrag);
            break;
        }
    }
    return dsdt;
}


/* ------------------------ CHARACTERISTIC TIMES ------------------------*/

double Tdrift(const double& R, const double& drdt)
{   return abs(AUtoMeter(R)/drdt);  }

double Tgrowth(const double& var, const double& dvardt)
{   return abs(var/dvardt);  }

double Tcoll(const double& size, const double& rhog, const double& filfac, const double& rhos,// ->
             const double& dustfrac, const double& vrel)
{
    return rhos*filfac*size/(3.*dustfrac*rhog*vrel);
}

double Ncoll(const double& dt, const double& tcoll)
{   return YearToSec(dt)/tcoll;   }
