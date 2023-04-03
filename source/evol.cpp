#include <iostream>
#include <cmath>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"
#include "../header/evol.h"

using namespace std;

/* ------------------------  TIME STEP ------------------------ */

double KeplerDt(double& omegakfraction, const double omegak)
{
    return SecToYear(omegakfraction*2.*M_PI/omegak);
}

double AdaptativeDt(const double& time, const double& timeend, const int& massorsize, const double& vargrowth,//-> 
                    const double& R, const double& dvargrowthdt, const double& dRdt)
{
    double dt;
    double dtgrowth = 0.05*abs(vargrowth/dvargrowthdt);    // dt from growth rate
    double dtdrift =  0.3*abs(R/dRdt);                     // dt from radial drift rate

    if (massorsize == 1)    dtgrowth /= 4.;

    if (dtgrowth != 0 && dtdrift != 0)
    {
        if (dtgrowth < dtdrift)     dt = dtgrowth;
        else                        dt = dtdrift;
    }
    else if (dtgrowth != 0 && dtdrift == 0)       dt = dtgrowth;
    else if (dtgrowth == 0 && dtdrift != 0)       dt = dtdrift;
    else
    {
        cerr << "Error: time step iteration impossible" << endl;
        exit(1);
    }

    dt = SecToYear(dt);

    //Prevent problem with large grains with large dt but small tend to get some data
    if (dt > timeend/100.)
    {   dt = timeend/100.;   }

    if (time+dt > timeend)
    {   dt = timeend - time;   }

    return dt;
}


/* ------------------------ RADIAL DRIFT ------------------------*/

double DRDt(const double& R, const double& Rin, const double& mstar, double p, double q, const double& rhog, const double& cg, const double& R0,//-> 
            const double& sigma0, const double& hg0, const double& dustfrac, const double& st, const double& alpha, const int& ibr,//->
            const int& ismooth, const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight)
{
    double dfbr = 1.;   //dust fraction back reaction
    if (ibr == 1)   dfbr += dustfrac;

    double vdrift = VDrift(R,Rin,mstar,p,q,rhog,cg,R0,sigma0,hg0,ismooth,ibump,Rbump,bumpwidth,bumpheight);
    double vvisc = VVisc(R,Rin,mstar,p,q,rhog,cg,R0,sigma0,hg0,alpha,ismooth,ibump,Rbump,bumpwidth,bumpheight);

    vdrift *= st/(dfbr*dfbr+st*st);
    vvisc *= dfbr/(dfbr*dfbr+st*st);

    return MeterToAU(vdrift+vvisc);
}


/* ------------------------  GROWTH-FRAG-BOUNCE dm/dt ------------------------ */

double DmDt(const double& size, const double& rhog, const double rhos, const double& dustfrac, const double& vrel,//->
            const int& ifrag, const int& ieros, const double& ejectasize, const double& cohacc, const int& ibounce,//->
            const double& vfrag, const double& vstick, const double deltav, const double& probabounce)
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
    if (ieros == 1)
    {   if (deltav > 0.00142460/sqrt(rhog*ejectasize))  dmdt -= rhog*pow(deltav,3)*size*GrainMass(ejectasize,1,rhos)/cohacc/ejectasize;
    }
    return dmdt;
}


/* ------------------------  GROWTH-FRAG ds/dt ------------------------ */

double DsDt(const double& size, const double& filfac, const double& rhog, const double& rhos,//-> 
            const double& dustfrac,const double& vrel, const int& ifrag, const int& ieros, const double& ejectasize,//->
            const double& cohacc, const double deltav, const double& vfrag, const double& filfacpow)
{
    double dsdt = dustfrac*rhog*vrel/(filfac*rhos);
    if (filfac != 1) dsdt /= 1.+filfacpow/3.;
    

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
    if (ieros == 1)
    {
        if (deltav > 0.00142461/sqrt(rhog*ejectasize))   dsdt -= rhog*pow(deltav,3)*ejectasize*ejectasize/3./size/cohacc;
    }

    return dsdt;
}


/* ------------------------ CHARACTERISTIC TIMES ------------------------*/

double Tdrift(const double& R, const double& drdt)
{   return abs(AUtoMeter(R)/drdt);  }

double Tgrowth(const double& var, const double& dvardt)
{   return abs(var/dvardt);  }

double Tcoll(const double& size, const double& rhog, const double& filfac, const double& rhos,//->
             const double& dustfrac, const double& vrel)
{
    return rhos*filfac*size/(3.*dustfrac*rhog*vrel);
}

double Ncoll(const double& dt, const double& tcoll)
{   return YearToSec(dt)/tcoll;   }
