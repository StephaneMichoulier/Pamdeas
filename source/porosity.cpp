#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"
#include "../header/porosity.h"

using namespace std;


/* ------------------------ SHARED FUNCTIONS ------------------------ */

double CParam(const double& R, const double& mstar, const double& rhog, const double& cg, const double& eroll,// ->
              const double& a0, const double& rhos, const double& alpha)
{
    return (243.*M_SQRT2*M_PI/15625.)*(rossby*alpha*pow(a0,4.)*rhos*rhos*cg*Omegak(R,mstar)/(rhog*b_oku*eroll));
}


/* ------------------------ BOUNCE ------------------------*/

double ProbaBounce(const double& filfac, const double& filfacbnc, const double& vrel, const double& vstick, const double& vend)
{
    if (filfac >= filfacbnc)
    {
        if (vrel <= vstick) 
        {   return 1.;  }
        else
        {
            if (vrel <= vend)
            {   return (log(vrel)-log(vend))/(log(vstick)-log(vend));   }
            else
            {   return 0.;  }
        }
    }
    else 
    {   return 1.;  }
}

double VarVolumeBounce(const double& filfac, const double& filfaclim, const double& coeffrest, const double& ekin,// ->
                       const double& volume, const double& Yd0, const double& Ydpower)
{
    double compactvol = (1.-coeffrest*coeffrest)*ekin/Yd(filfac,filfaclim,Yd0,Ydpower);

    if (compactvol >= volume)   
    {   compactvol = volume;    }
    
    return compactvol;
}

double FilFacBounce(const double& sizei, const double& filfaci, const double& rhos, const double& filfaclim,// ->
                    const double& vrel, const double& vstick, const double& vyield, const double& ncoll,// ->
                    const double& Yd0, const double& Ydpower)
{
    double coeffrest = CoeffRest(vrel,vstick,vyield);
    double volume = GrainVolumeSize(sizei,filfaci,rhos);
    double ekin = Ekin(sizei,filfaci,rhos,vrel);
    double varvolbounce = VarVolumeBounce(filfaci,filfaclim,coeffrest,ekin,volume,Yd0,Ydpower);
    double filfacbounce = filfaci*pow(1./(1.-(0.5*varvolbounce/volume)),ncoll);

    if (vrel <= vyield)
    {   return filfaci; }
    else
    {
        if (filfacbounce <= 1.)
        {   return filfacbounce;    }
        else
        {   return 1.;  }
    }
}


/* ------------------------ MASS MODEL------------------------ */

/* ------------------------ OKUZUMI GROWTH MODEL------------------------ */

double M1(const double& cparam)
{   return pow(cparam/(2.*(pow(2.,0.075)-1.)),0.375/(cratio+0.125)); }

double M2(const double& cparam, const double& rhog, const double& cg, const double& a0)
{   return pow(cparam*cg*a0/(9.*NuMolGas(rhog,cg)*(pow(2.,0.2)-1.)),1./(3.*cratio));   }

double M3(const double& m1_on_m0, const double& m2_on_m0)
{   return pow(m1_on_m0,8.*cratio+1.)*pow(m2_on_m0,-8.*cratio); }

double M3(const double& cparam, const double& rhog, const double& cg, const double& a0)
{
    return pow(cparam/(2.*(pow(2.,0.075)-1.)),3.)*pow(cparam*cg*a0/(9.*NuMolGas(rhog,cg)*(pow(2.,0.2)-1.)),-8./3.);
}

double M4(const double& R, const double cparam, const double& mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos)
{
    return pow(rhog*cg/(rhos*a0*Omegak(R,mstar)),4.)*2.*(pow(2.,0.075)-1.)/cparam;
}

double M5(const double& R, const double cparam, const double& mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos)
{
    double nu = NuMolGas(rhog,cg);
    return pow(9.*nu*rhog/(2.*rhos*a0*a0*Omegak(R,mstar)),1.5)*pow(9.*nu*(pow(2.,0.2)-1.)/(cg*a0*cparam),1./6.);
}

double FilFacMGr(const double& massf, const double& massi, const double& filfaci, const double& rhos, const double& eroll, const double& vrel)
{
    double power;
    if (Ekin(massi,vrel)/(3.*b_oku*eroll) <= 1.)
    {   power = cratio;   }
    else
    {   power = -0.2; }

    return filfaci*pow(massf/massi,power);
}

/*double FilFacMColl(const double& R, const double& mstar, const double& rhog, const double& cg, const double& massf,// ->
                const double eroll, const double& a0, const double& rhos, const double& alpha)
{
    // here, m1 to m5 and the associated functions M1 to M5 are normalized by m0 !
    double cparam = CParam(R,mstar,rhog,cg,eroll,a0,rhos,alpha);
    double m = massf/GrainMass(a0,1.,rhos);
    double m1 = M1(cparam);
    double m2 = M2(cparam,rhog,cg,a0);
    double m3 = M3(m1,m2);
    double m4 = M4(R,cparam,mstar,rhog,cg,a0,rhos);
    double m5 = M5(R,cparam,mstar,rhog,cg,a0,rhos);
    double filfaccoll;

    if (m1 <= m2)
    {
        if (m <= m1)
        {
            filfaccoll = pow(m,cratio);   // FilFac h&s
        }
        else
        {
            if (m3 <= m4)
            {
                if (m <= m3)
                {
                    filfaccoll = pow(m1,cratio+0.125)*pow(m,-0.125);   // FilFac Ep-St<1
                }
                else
                {
                    if (m <= m5)
                    {
                        filfaccoll = pow(m2,cratio);   // FilFac St-St<1
                    }
                    else
                    {
                        filfaccoll = pow(m2,cratio)*pow(m5/m,0.2);  // FilFac St-St>1
                    }
                }
            }
            else
            {
                if (m <= m4)
                {
                    filfaccoll = pow(m1,cratio+0.125)*pow(m,-0.125);   // FilFac Ep-St<1
                }
                else
                {
                    filfaccoll = pow(m1,cratio+0.125)*pow(m4,-0.125)*pow(m4/m,0.2);   // FilFac Ep-St>1
                }
            }
        }
    }
    else
    {   if (m <= m2)
        {
            filfaccoll = pow(m,cratio);   // FilFac h&s
        }
        else
        {
            if (m <= m5)
            {
                filfaccoll = pow(m2,cratio);   // FilFac St-St<1
            }
            else
            {
                filfaccoll = pow(m2,cratio)*pow(m5/m,0.2);  // FilFac St-St>1
            }
        }
    }
    return filfaccoll;
}*/

double FilFacMColl(const double& R, const double& mstar, const double& rhog, const double& cg, const double& st,// ->
                   const double& mfrac, const double eroll, const double& a0, const double& rhos, const double& alpha)
{
    // here, m1 to m5 and the associated functions M1 to M5 are normalized by m0 !
    double cparam = CParam(R,mstar,rhog,cg,eroll,a0,rhos,alpha);
    double m1 = M1(cparam);
    double m2 = M2(cparam,rhog,cg,a0);
    double m4;
    double m5;
    double filfaccoll;

    if (st < 1.)
    {
        if (m1 <= m2)
        {
            if (mfrac <= m1)
            {   filfaccoll = pow(mfrac,cratio);    }   // FilFac h&s
            else
            {   if (mfrac <= M3(m1,m2))
                {
                    filfaccoll = pow(m1,cratio+0.125)/pow(mfrac,0.125);   // FilFac Ep-St<1
                }
                else
                {
                    filfaccoll = pow(m2,cratio);   // FilFac St-St<1
                }
            }
        }
        else
        {
            if (mfrac <= m2)
            {   filfaccoll = pow(mfrac,cratio);    }   // FilFac h&s
            else
            {   filfaccoll = pow(m2,cratio);   }   // FilFac St-St<1
        }
    }
    else
    {
        m4 = M4(R,cparam,mstar,rhog,cg,a0,rhos);
        m5 = M5(R,cparam,mstar,rhog,cg,a0,rhos);

        if (m4 <= m5)    filfaccoll = pow(m1,cratio+0.125)*pow(m4,-0.125)*pow(m4/mfrac,0.2);   // FilFac Ep-St>1
        else             filfaccoll = pow(m2,cratio)*pow(m5/mfrac,0.2);  // FilFac St-St>1
    }
    return filfaccoll;
}


/* ------------------------ KATAOKA ------------------------*/

// Original algorithm
/*double FilFacMGas(const double& R, const double& mstar, const double& p, const double& q, const double& rhog,// ->
               const double& cg, const double& massf, const double& rhos, const double& eroll, const double& a0)
{
    double m0 = GrainMass(a0,1.,rhos);
    double filfaclimSt = pow(m0*rhog*DeltaV(R,mstar,p,q,cg)*cg/(M_PI*eroll*rhos),1./3.);

    if  (TransRegEpSt(rhog,cg,GrainMassToSize(massf,filfaclimSt,rhos)) < 1.)
    {   return filfaclimSt;    }
    else
    {   return pow(6.*a0*a0*DeltaV(R,mstar,p,q,cg)*NuMolGas(rhog,cg)*rhog/(eroll),0.375)*pow(massf/m0,-0.125); }
}*/

double FilFacMGas(const double& R, const double& mstar, const double& deltav, const double& st, const double& massf, const double& rhos,// -> 
                  const double& eroll, const double& a0, const double& m0)//, const int& iregime, const double& rhog, const double& cg)
{
    /*double filfacgas;
    switch (iregime)
    {
        case (1):
        {
            filfacgas = pow(m0*rhog*deltav*cg/(M_PI*eroll*rhos),1./3.);
            break;
        }
        case (2):
        {
            filfacgas = pow(6.*a0*a0*deltav*NuMolGas(rhog,cg)*rhog/(eroll),3./8.)*pow(massf/m0,-0.125);
            break;
        }
        case (3):
        {
            filfacgas = pow((9.*m0*deltav*deltav*rhog)/(M_PI*eroll*rhos),5./14.)*pow(NuMolGas(rhog,cg)/(2.*a0*deltav),3./14.)/pow(massf/m0,1./14.);
            break;
        }
        case (4):
        {
            filfacgas = pow(0.165*m0*rhog*deltav*deltav/(M_PI*eroll*rhos),1./3.);
            break;
        }
    }
    return filfacgas;*/
    return pow((m0*a0*deltav*Omegak(R,mstar))/(eroll*M_PI*st),3./7.)*pow(massf/m0,1./7.);
}

double FilFacMGrav(const double& massf, const double& rhos, const double& eroll, const double& a0, const double& m0)
{    return pow(G*m0*m0/(eroll*M_PI*a0),0.6)*pow(massf/m0,0.4) ;    }


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

double FilFacMinMColGasGrav(const double& R, const double& mstar, const double& rhog, const double cg, const double& deltav,// ->
                            const double st, const double& massf, const double& rhos, const double& eroll, const double& a0,// ->
                            const double& m0, const double& alpha)//, const int& iregime)
{ 
    double filfaccoll = FilFacMColl(R,mstar,rhog,cg,st,massf/m0,eroll,a0,rhos,alpha);
    double filfacgas  = FilFacMGas(R,mstar,deltav,st,massf,rhos,eroll,a0,m0);
    double filfacgrav = FilFacMGrav(massf,rhos,eroll,a0,m0);

    double array[] = {filfaccoll,filfacgas,filfacgrav};
    double filfacmin = *max_element(array,array+3);

    return filfacmin;
}

double FilFacMFinal(const double& R, const double& mstar, const double& rhog, const double& cg, const double deltav, const double st,// ->
                    const double& massf, const double& massi, const double& filfaci, const double& a0, const double& rhos, const double& eroll,// ->
                    const double& alpha, const double& ncoll, const int& ifrag, const int& ibounce, const double& vfrag, const double& vrel,// ->
                    const double& vstick, const double& probabounce, const double& filfaclim, const double& Yd0, const double& Ydpower)
                    //const int& iregime)
{
    double filfacf = 1.;
    double m0 = GrainMass(a0,1.,rhos);
    double sizei = GrainMassToSize(massi,filfaci,rhos);
    double filfacmincollgasgrav = FilFacMinMColGasGrav(R,mstar,rhog,cg,deltav,st,massf,rhos,eroll,a0,m0,alpha);
    double filfacgr = FilFacMGr(massf,massi,filfaci,rhos,eroll,vrel);

    if (massf > m0)
    {
        switch (ibounce)
        {   case (0):
            {
                if ((vrel < vfrag && ifrag > 0) || (ifrag == 0))
                {
                    filfacf = filfacgr;
                    break;
                }
                else
                {
                    filfacf = filfaci;
                    break;
                }
            }
            case (1):
            {
                if ((vrel < vfrag && ifrag > 0) || (ifrag == 0))
                {
                    if (vrel <= vstick)
                    {
                        filfacf = filfacgr;
                        break;
                    }
                    else
                    {
                        double filfacbounce = FilFacBounce(sizei,filfaci,rhos,filfaclim,vrel,vstick,Vyield(vstick),ncoll,Yd0,Ydpower);

                        if (vrel < Vend(vstick))
                        {
                            if (filfacgr < filfacmincollgasgrav)
                            {   filfacf = filfacmincollgasgrav;   }
                            else
                            {   filfacf = filfacgr;   }

                            filfacf = filfacf*probabounce+(1.-probabounce)*filfacbounce;
                            break;
                        }
                        else
                        {
                            filfacf = filfacbounce;
                            break;
                        }
                    }
                }
                else
                {
                    filfacf = filfaci;
                    break;
                }
            }
        }
        if (filfacf < filfacmincollgasgrav)
        {   filfacf = filfacmincollgasgrav;   }
        
        if (filfacf > 1.)
        {   filfacf = 1.;    }

        return filfacf;
    }
    else 
    {   return filfacf; }
}


/* ------------------------ SIZE MODEL ------------------------ */

/* ------------------------ OKUZUMI GROWTH MODEL------------------------ */

double S1(const double& cparam)
{   return pow(cparam/(2.*(pow(2.,0.075)-1.)),(1.-cratio)/(1.+8.*cratio));  }

double S2(const double& cparam, const double& rhog, const double& cg, const double& a0)
{   return pow(cparam*cg*a0/(9.*NuMolGas(rhog,cg)*(pow(2.,0.2)-1.)),(1.-cratio)/(9.*cratio));  }

double S3(const double& rhog, const double& cg, const double& a0)
{   return (9.*NuMolGas(rhog,cg)/(2.*cg*a0)*((pow(2.,0.2)-1.)/(pow(2.,0.075)-1.)));    }

double S4(const double& R, const double cparam, const double& mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos)
{
    return pow(rhog*cg/(rhos*a0*Omegak(R,mstar)),1.5)*sqrt(2.*(pow(2.,0.075)-1.)/cparam);
}

double S5(const double& R, const double cparam, const double& mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos)
{
    double nu = NuMolGas(rhog,cg);
    return sqrt(9.*nu*rhog/(2.*rhos*a0*a0*Omegak(R,mstar)))*pow(9.*nu*(pow(2.,0.2)-1.)/(cg*a0*cparam),1./6.);
}

double FilFacSGr(const double& sizef, const double& sizei, const double& filfaci,// ->
              const double& rhos, const double& eroll, const double& vrel)
{
    double power;
    if (Ekin(sizei,filfaci,rhos,vrel)/(3.*b_oku*eroll) <= 1.)
    {   power = 3.*cratio/(1.-cratio);  }
    else
    {   power = -0.5;   }

    return filfaci*pow(sizef/sizei,power);
}

double FilFacSColl(const double& R, const double& mstar, const double& rhog, const double& cg, const double st, const double& sfrac,// ->
                   const double& eroll, const double& a0, const double& rhos, const double& alpha, const int& iregime, double& filfacpow)
{
    // here, s1 to s5 and the associated functions S1 to S5 are normalized by a0 !
    double cparam = CParam(R,mstar,rhog,cg,eroll,a0,rhos,alpha);
    double s1 = S1(cparam);
    double s2 = S2(cparam,rhog,cg,a0);
    double s4 = S4(R,cparam,mstar,rhog,cg,a0,rhos);
    double s5 = S5(R,cparam,mstar,rhog,cg,a0,rhos);
    double filfaccoll;

    if (st < 1.)
    {   if (s1 <= s2)
        {
            if (sfrac <= s1)
            {
                filfacpow = 3.*cratio/(1.-cratio);
                filfaccoll = pow(sfrac,filfacpow);   // FilFac h&s
            }
            else
            {   if (sfrac <= S3(rhog,cg,a0))
                {
                    filfacpow = -1./3.;
                    filfaccoll = pow(s1,(1.+8.*cratio)/(3.*(1.-cratio)))*pow(sfrac,filfacpow);   // FilFac Ep-St<1
                }
                else
                {   filfacpow = 0.;
                    filfaccoll = pow(s2,3.*cratio/(1.-cratio));   // FilFac St-St<1

                }
            }
        }
        else
        {
            if (sfrac <= s2)
            {
                filfacpow = 3.*cratio/(1.-cratio);
                filfaccoll = pow(sfrac,filfacpow);   // FilFac h&s
            }
            else
            {
                filfacpow = 0.;
                filfaccoll = pow(s2,3.*cratio/(1.-cratio));   // FilFac St-St<1
            }
        }
    }
    else
    {
        filfaccoll = pow((cparam*rhog*cg)/(2.*rhos*a0*Omegak(R,mstar)),0.25)/sqrt(sfrac);
        if (s4 <= s5)    filfaccoll *= pow(pow(2.,0.075)-1.,-0.25);
        else             filfaccoll *= pow(pow(2.,0.2)-1.,-0.25);
        filfacpow = -0.5;
    }
    return filfaccoll;
}

/* ------------------------ KATAOKA ------------------------*/

double FilFacSGas(const double& R, const double& mstar, const double& deltav, const double& st,const double& sizef, const double& eroll,// ->
                  const double& a0, const double& rhos, const int& iregime, double& filfacpow)//, const double& rhog, const double& cg)
{
    //double filfacgas;
    switch (iregime)
    {
        case (1):
        {
            filfacpow = 0.;
            //filfacgas = pow(4.*a0*a0*a0*deltav*rhog*cg/(3.*eroll),1./3.);
            break;
        }
        case (2):
        {
            filfacpow = -1./3.;
            //filfacgas = pow(6.*a0*a0*a0*deltav*NuMolGas(rhog,cg)*rhog/(sizef*eroll),1./3.);
            break;
        }
        case (3):
        {
            filfacpow = -0.2;
            //filfacgas = pow((12.*a0*a0*a0*deltav*deltav*rhog)/eroll,1./3.)*pow(NuMolGas(rhog,cg)/(2.*a0*deltav),0.2)/pow(sizef/a0,0.2);
            break;
        }
        case (4):
        {
            filfacpow = 0.;
            //filfacgas = pow(0.22*a0*a0*a0*rhog*deltav*deltav/eroll,1./3.);
            break;
        }
    }
    //return filfacgas;
    return pow(4.*a0*a0*a0*rhos*deltav*sizef*Omegak(R,mstar)/(3.*eroll*st),0.5);
}

double FilFacSGrav(const double& sizef, const double& rhos, const double& a0, const double& eroll)
{   return 16.*M_PI*G*rhos*rhos*a0*a0*a0*sizef*sizef/(9.*eroll);  }


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

double FilFacMinSColGasGrav(const double& R, const double& mstar, const double& rhog, const double& cg, const double& deltav,// ->
                            const double st, const double& sizef, const double& rhos, const double& eroll, const double& a0,// ->
                            const double& alpha, const int& iregime, double& filfacpow)
{
    double filfacmin;
    double filfacpowgas;
    int maxindex;
    double filfaccoll = FilFacSColl(R,mstar,rhog,cg,st,sizef/a0,eroll,a0,rhos,alpha,iregime,filfacpow);
    double filfacgas  = FilFacSGas(R,mstar,deltav,st,sizef,eroll,a0,rhos,iregime,filfacpowgas);
    double filfacgrav = FilFacSGrav(sizef,rhos,a0,eroll);

    double array[] = {filfaccoll,filfacgas,filfacgrav};
    double arrayindex[] = {filfacpow,filfacpowgas,2.};

    maxindex = distance(array,max_element(array,array+3));
    filfacmin = array[maxindex];
    filfacpow = arrayindex[maxindex];

    return filfacmin;
}

double FilFacSFinal(const double& R, const double& mstar, const double& rhog, const double& cg, const double& deltav, const double st,// ->
                    const double& sizef, const double& sizei, const double& filfaci, const double& a0, const double& rhos, const double& eroll,// -> 
                    const double& alpha, const int& ifrag, const double& vfrag, const double& vrel, const int& iregime, double& filfacpow)
{
    double filfacf = 1.;
    double filfacmincollgasgrav = FilFacMinSColGasGrav(R,mstar,rhog,cg,deltav,st,sizef,rhos,eroll,a0,alpha,iregime,filfacpow);
    double filfacgr = FilFacSGr(sizef,sizei,filfaci,rhos,eroll,vrel);

    if (sizef > a0)
    {
        if ((vrel < vfrag && ifrag > 0) || (ifrag == 0))
        {   filfacf = filfacgr; }
        else
        {   filfacf = filfaci;  }

        if (filfacf <= filfacmincollgasgrav)
        {   filfacf = filfacmincollgasgrav;     }

        if (filfacf > 1.)
        {
            filfacf = 1;
            filfacpow = 3.*cratio/(1.-cratio);
        }

        return filfacf;
    }
    else 
    {   return filfacf; }

}
