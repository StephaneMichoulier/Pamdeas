#include <iostream>
#include <cmath>

#include "../header/constant.h"
#include "../header/general_functions.h"
#include "../header/evol.h"

using namespace std;

/* ------------------------  TIME STEP ------------------------ */

double KepletDt(double& omegakfraction, const double omegak)
{
    return SecToYear(omegakfraction*2.*M_PI/omegak);
}

double AdaptativeDt(const double& time, const double& timeend, const int& massorsize, const int& ibump,// ->
                    const double& vargrowth, const double& R, const double& dvargrowthdt, const double& dRdt)
{
    double dt = 1.;
    double limup = 0.15;
    double limdown = 0.075;
    double Ccdt;
    double C1 = abs(vargrowth/dvargrowthdt);
    double C2 = abs(R/dRdt);

    if (C1 != 0 && C2 != 0)
    {
        if (C1 < C2)        Ccdt = C1;
        else                Ccdt = C2;
    }
    else if (C1 != 0 && C2 == 0)       Ccdt = C1;
    else if (C1 == 0 && C2 != 0)       Ccdt = C2;
    else
    {
        cout << "Error, time step iteration impossible" << endl;
        exit(1);
    }

    if (ibump == 1)
    {   limup *= 2.;   limdown *= 2.;   }

    if (massorsize == 1)
    {   limup /= 3.;   limdown /= 3.;   }

    Ccdt = SecToYear(Ccdt);
    limup *= Ccdt;
    limdown *= Ccdt;
    do
    {
        if (dt > limup)     dt /= 2.;
        if (dt < limdown)   dt *= 2.;

    } while (limup < dt || limdown >  dt);

    if (time+dt > timeend)
    {   dt = timeend - time;   }

    return dt;
}


/* ------------------------ DRIFT ------------------------*/

double DRDt(const double& R, const double& Mstar, double p, double q, const double& hg, const double& rhog, const double& cg,// ->
            const double& sigma0, const double& R0, const double& Hg0, const double& dustfrac, const double& st, const int& ibr,// ->
            const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight)
{
    double deriver;
    double vdrift;
    double dfbr = 0;

    deriver = Pg(R+DeltaR,Mstar,p,q,sigma0,R0,Hg0,ibump,Rbump,bumpwidth,bumpheight);
    deriver -= Pg(R-DeltaR,Mstar,p,q,sigma0,R0,Hg0,ibump,Rbump,bumpwidth,bumpheight);
    deriver /= 2.*AUtoMeter(DeltaR);

    if (ibr == 1)   dfbr = dustfrac;

    vdrift = (AUtoMeter(hg)*AUtoMeter(hg))/(AUtoMeter(R)*Pg(rhog,cg))*deriver;
    vdrift *= Vk(R,Mstar)*st/((1.+dfbr)*(1.+dfbr)+st*st);

    return MeterToAU(vdrift);
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
            if ((vrel < vfrag && ifrag > 0) || ifrag == 0)
            {
                    break;
            }
            else switch (ifrag)
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
            if ((vrel < vfrag && ifrag > 0) || ifrag == 0)
            {
                if (vrel <= vstick)
                {
                    break;
                }
                else
                {
                    dmdt *= probabounce;
                    break;
                }
            }
            else switch (ifrag)
            {   case (1):
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
    }
    return dmdt;
}


/* ------------------------  GROWTH-FRAG-BOUNCE ds/dt ------------------------ */

double DsDt(const double& phi, const double& rhog, const double& rhos, const double& dustfrac,// ->
            const double& vrel, const int& ifrag, const double& vfrag, const double& phipow)
{
    double dsdt = dustfrac*rhog*vrel/(phi*rhos*(1.+phipow/3.));

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

double Tcoll(const double& size, const double& rhog, const double& phi, const double& rhos,// ->
             const double& dustfrac, const double& vrel)
{
    return rhos*phi*size/(3.*dustfrac*rhog*vrel);
}

double Ncoll(const double& dt, const double& tcoll)
{   return YearToSec(dt)/tcoll;   }
