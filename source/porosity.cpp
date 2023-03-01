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

double CParam(const double& R, const double& mstar, const double& rhog, const double& cg, const double& eroll,//->
              const double& a0, const double& rhos, const double& alpha)
{
    return (243.*M_SQRT2*M_PI/15625.)*(rossby*alpha*pow(a0,4.)*rhos*rhos*cg*Omegak(R,mstar)/(rhog*b_oku*eroll));
}

/* ------------------------ BOUNCE ------------------------*/

double ProbaBounce(const double& filfac, const double& filfacbnc, const double& vrel, const double& vstick, const double& vend)
{
    if (filfac >= filfacbnc)
    {
        if (vrel < vstick) 
        {   return 1.;  }
        else
        {
            if (vrel < vend)
            {   return (log(vrel)-log(vend))/(log(vstick)-log(vend));   }
            else
            {   return 0.;  }
        }
    }
    else 
    {   return 1.;  }
}

double VarVolumeBounce(const double& filfac, const double& filfaclim, const double& coeffrest, const double& ekin,//->
                       const double& volume, const double& Yd0, const double& Ydpower)
{
    double compactvol = (1.-coeffrest*coeffrest)*ekin/Yd(filfac,filfaclim,Yd0,Ydpower);
    
    if (compactvol >= volume)   
    {   compactvol = volume;    }
    
    return compactvol;
}

double FilFacBounce(const double& sizei, const double& filfaci, const double& rhos, const double& filfaclim,//->
                    const double& vrel, const double& vstick, const double& vyield, const double& ncoll,//->
                    const double& Yd0, const double& Ydpower)
{
    double coeffrest = CoeffRest(vrel,vstick,vyield);
    double volume = GrainVolumeSize(sizei,filfaci,rhos);
    double ekin = Ekin(sizei,filfaci,rhos,vrel);
    double varvolbounce = VarVolumeBounce(filfaci,filfaclim,coeffrest,ekin,volume,Yd0,Ydpower);
    double filfacbounce = filfaci*pow(1./(1.-0.5*(varvolbounce/volume)),ncoll);

    if (vrel < vyield)
    {   return filfaci; }
    else
    {
        if (filfacbounce <= 1.)
        {   return filfacbounce;    }
        else
        {   return 1.;  }
    }
}

double FilFacBouncebis(const double& sizei, const double& filfaci, const double& rhos,//->
                    const double& vrel, const double& vstick, const double& vyield, const double& ncoll,//->
                    const double& eroll,const double& a0)
{
    double coeffrest = CoeffRest(vrel,vstick,vyield);
    double volume = GrainVolumeSize(sizei,filfaci,rhos);
    double ekin = Ekin(sizei,filfaci,rhos,vrel);
    double varvolbounce  = (1.-coeffrest*coeffrest)*ekin/(eroll*pow(filfaci/(maxpacking-filfaci)/a0,3));
    if (varvolbounce >= volume)   
    {   varvolbounce = volume;    }
    double filfacbounce = filfaci*pow(1./(1.-0.5*(varvolbounce/volume)),ncoll);

    if (vrel < vyield)
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

double M4(const double& R, const double cparam, const double& mstar, const double& rhog,//->
          const double& cg, const double& a0, const double& rhos)
{
    return pow(rhog*cg/(rhos*a0*Omegak(R,mstar)),4.)*2.*(pow(2.,0.075)-1.)/cparam;
}

double M5(const double& R, const double cparam, const double& mstar, const double& rhog,//->
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

double FilFacMColl(const double& R, const double& mstar, const double& rhog, const double& cg, const double& st,//->
                   const double& mfrac, const double eroll, const double& a0, const double& rhos, const double& alpha,//->
                   int& porreg)
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
            {   
                filfaccoll = pow(mfrac,cratio);         // FilFac h&s
                porreg = 1;

            }   
            else
            {   if (mfrac <= M3(m1,m2))
                {
                    filfaccoll = pow(m1,cratio+0.125)/pow(mfrac,0.125);   // FilFac Ep-St<1
                    porreg = 2;
                }
                else
                {
                    filfaccoll = pow(m2,cratio);   // FilFac St-St<1
                    porreg = 3;
                }
            }
        }
        else
        {
            if (mfrac <= m2)
            {   
                filfaccoll = pow(mfrac,cratio);     // FilFac h&s
                porreg = 1;
            }   
            else
            {   
                filfaccoll = pow(m2,cratio);      // FilFac St-St<1
                porreg = 3;
            }
        }
    }
    else
    {
        m4 = M4(R,cparam,mstar,rhog,cg,a0,rhos);
        m5 = M5(R,cparam,mstar,rhog,cg,a0,rhos);

        if (m4 <= m5)    
        {
            filfaccoll = pow(m1,cratio+0.125)*pow(m4,-0.125)*pow(m4/mfrac,0.2);   // FilFac Ep-St>1
            porreg = 4;
        }
        else             
        {
            filfaccoll = pow(m2,cratio)*pow(m5/mfrac,0.2);  // FilFac St-St>1
            porreg = 5;
        }
    }
    return filfaccoll;
}


/* ------------------------ KATAOKA ------------------------*/

double FilFacMGas(const double& R, const double& mstar, const double& deltav, const double& st, const double& massf, const double& rhos,//-> 
                  const double& eroll, const double& a0, const double& m0)//, const int& dragreg, const double& rhog, const double& cg)
{
    return pow((m0*a0*deltav*Omegak(R,mstar))/(eroll*M_PI*st),3./7.)*pow(massf/m0,1./7.);
}

double FilFacMGrav(const double& massf, const double& rhos, const double& eroll, const double& a0, const double& m0)
{    return pow(G*m0*m0/(eroll*M_PI*a0),0.6)*pow(massf/m0,0.4) ;    }


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

double FilFacMinMColGasGrav(const double& R, const double& mstar, const double& rhog, const double cg, const double& deltav,//->
                            const double st, const double& massf, const double& rhos, const double& eroll, const double& a0,//->
                            const double& m0, const double& alpha, int& porreg)
{ 
    int maxindex;
    double filfacmin;
    double filfaccoll = FilFacMColl(R,mstar,rhog,cg,st,massf/m0,eroll,a0,rhos,alpha,porreg);
    double filfacgas  = FilFacMGas(R,mstar,deltav,st,massf,rhos,eroll,a0,m0);
    double filfacgrav = FilFacMGrav(massf,rhos,eroll,a0,m0);

    double filfacminarray[] = {filfaccoll,filfacgas,filfacgrav};
    int porregindex[] = {porreg,6,7};

    maxindex = distance(filfacminarray,max_element(filfacminarray,filfacminarray+3));
    filfacmin = filfacminarray[maxindex];
    porreg = porregindex[maxindex];

    return filfacmin;
}

double FilFacMFinal(const double& R, const double& mstar, const double& rhog, const double& cg, const double deltav, const double st,//->
                    const double& massf, const double& massi, const double& filfaci, const double& a0, const double& rhos, const double& eroll,//->
                    const double& alpha, const double& ncoll, const int& ifrag, const int& ibounce, const double& vfrag, const int& icomp,//->
                    const double& vrel, const double& vstick, const double& probabounce, const double& filfaclim, const double& Yd0,//-> 
                    const double& Ydpower, int& porreg)
{
    double filfacf = 1.;
    double m0 = GrainMass(a0,1.,rhos);
    double sizei = GrainMassToSize(massi,filfaci,rhos);
    double filfacmincollgasgrav = FilFacMinMColGasGrav(R,mstar,rhog,cg,deltav,st,massf,rhos,eroll,a0,m0,alpha,porreg);
    double filfacgr = FilFacMGr(massf,massi,filfaci,rhos,eroll,vrel);
    double filfacmax = maxpacking + (1.-maxpacking)*(m0/massf);
    int model = 5;

    if (massf > m0)
    {
        switch (ibounce)
        {   case (0):
            {
                if ((vrel < vfrag && ifrag > 0) || (ifrag == 0))
                {
                    filfacf = filfacgr;
                
                    //In developement
                    
                    /*if (icomp == 1)
                    {   double vrelvf = vrel/vfrag;
                        if (vrelvf < 0.1)
                        {   filfacf = filfacgr;  }
                        else if (vrelvf <= 0.3)
                        {   filfacf = filfaci * (2.9*vrelvf*pow(filfacf,-0.2) +1); }
                        else    
                        {   filfacf = filfaci * (2.9*0.3*pow(filfacgr,-0.2)*exp(1.5*(0.3-vrelvf)) + 1);  }
                    }*/
                    break;
                }
                else
                { 
                    if (icomp == 1)
                    {   
                        switch (model)
                        {
                            case (0): // Fit1ncoll
                            {
                                double vrelvf = vrel/vfrag;
                                double xi =(3*0.3*pow(filfaci,-0.2)*exp(1.5*(0.3-vrelvf))+1);
                                filfacf = filfaci*pow(xi,ncoll);
                                break;
                            }
                            case (1):// fit2ncoll
                            {
                                double vrelvf = vrel/vfrag;
                                double xi2 = 27.*pow(filfaci,-0.2)*pow(vrelvf,1.5)/(2.*exp(4.*vrelvf)-1) +1;
                                filfacf = filfaci*pow(xi2,ncoll);
                                break;
                            }
                            case (2): // Garcia 
                            {    
                                double Vi = GrainVolumeSize(sizei,filfaci,rhos);
                                double deltaVol = VarVolumeBounce(filfaci,filfaclim,0,Ekin(massi,vrel),Vi,Yd0,Ydpower);
                                filfacf = filfaci*pow(1/(1-0.5*(deltaVol/Vi)),ncoll);
                                break;
                            }
                            case (3):// Fit1 + ms
                            {
                                double vrelvf = vrel/vfrag;
                                double xi =(3*0.3*pow(filfaci,-0.2)*exp(1.5*(0.3-vrelvf))+1);
                                if (massf != massi)
                                {   
                                filfacf = filfaci * pow( massf/massi, -log(xi)/log(1.+vrelvf*vrelvf) ); 
                                }
                                else
                                {   
                                    filfacf = filfaci * pow(xi, ncoll/1.4);   
                                }
                                break;
                            }
                            case (4):// Fit1 + garcia
                            {
                                double vrelvf = vrel/vfrag;
                                double Vi = GrainVolumeSize(sizei,filfaci,rhos);
                                double xi =(3*0.3*pow(filfaci,-0.2)*exp(1.5*(0.3-vrelvf))+1);
                                double deltaVol = Vi - massf/massi*Vi/xi;
                                filfacf = filfaci*pow(1./(1.-(deltaVol/Vi)),ncoll);
                                break;
                            }
                            case (5): // Garcia + mod compressive strenght kataoka
                            {    
                                double Vi = GrainVolumeSize(sizei,filfaci,rhos);
                                double Ebreak = 1.5*48./302.46*eroll;
                                //double Ebreak = 6*272.21/302.46*eroll;
                                double Ecomp = Ekin(massi,vrel) - (2.*massi-massf)*Ebreak/m0;

                                if (Ecomp<0) Ecomp=0;

                                double deltaVol = (Ecomp/(eroll*pow(filfaci/(maxpacking-filfaci)/a0,3)));
                                if (deltaVol > Vi)  deltaVol = Vi;
                                   

                                filfacf = filfaci*pow( 1./ (1. - 0.5*exp(1-pow(vrel/vfrag,2))*(deltaVol/Vi) ),ncoll);
                               // cout << "ratiofil = " << filfacf/filfaci << endl;
                                break;
                            }
                        }
                    }
                    else    filfacf = filfaci;
                    break;
                }
            }
            case (1):
            {
                if ((vrel < vfrag && ifrag > 0) || (ifrag == 0))
                {
                    if (vrel < vstick)
                    {
                        filfacf = filfacgr;
                        break;
                    }
                    else
                    {   
                        int bouncemod = 1;
                        double filfacbounce = 0;

                        if (bouncemod == 0)
                        {
                            filfacbounce = FilFacBounce(sizei,filfaci,rhos,filfaclim,vrel,vstick,Vyield(vstick),ncoll,Yd0,Ydpower);
                        }
                        else
                        {
                            filfacbounce = FilFacBouncebis(sizei,filfaci,rhos,vrel,vstick,Vyield(vstick),ncoll,eroll,a0);
                        }

                        if (vrel < Vend(vstick))
                        {   
                            if (filfacgr < filfacmincollgasgrav)
                            {   
                                filfacf = filfacmincollgasgrav;   
                            }
                            else
                            {   
                                filfacf = filfacgr;  
                            }
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
                    if (icomp == 1)
                    {   
                        switch (model)
                        {
                            case (0): // Fit1ncoll
                            {
                                double vrelvf = vrel/vfrag;
                                double xi =(3*0.3*pow(filfaci,-0.2)*exp(1.5*(0.3-vrelvf))+1);
                                filfacf = filfaci*pow(xi,ncoll);
                                break;
                            }
                            case (1):// fit2ncoll
                            {
                                double vrelvf = vrel/vfrag;
                                double xi2 = 27.*pow(filfaci,-0.2)*pow(vrelvf,1.5)/(2.*exp(4.*vrelvf)-1) +1;
                                filfacf = filfaci*pow(xi2,ncoll);
                                break;
                            }
                            case (2): // Garcia 
                            {    
                                double Vi = GrainVolumeSize(sizei,filfaci,rhos);
                                double deltaVol = VarVolumeBounce(filfaci,filfaclim,0,Ekin(massi,vrel),Vi,Yd0,Ydpower);
                                filfacf = filfaci*pow(1/(1-0.5*(deltaVol/Vi)),ncoll);
                                break;
                            }
                            case (3):// Fit1 + ms
                            {
                                double vrelvf = vrel/vfrag;
                                double xi =(3*0.3*pow(filfaci,-0.2)*exp(1.5*(0.3-vrelvf))+1);
                                if (massf != massi)
                                {   
                                filfacf = filfaci * pow( massf/massi, -log(xi)/log(1.+vrelvf*vrelvf) ); 
                                }
                                else
                                {   
                                    filfacf = filfaci * pow(xi, ncoll/1.4);   
                                }
                                break;
                            }
                            case (4):// Fit1 + garcia
                            {
                                double vrelvf = vrel/vfrag;
                                double Vi = GrainVolumeSize(sizei,filfaci,rhos);
                                double xi =(3*0.3*pow(filfaci,-0.2)*exp(1.5*(0.3-vrelvf))+1);
                                double deltaVol = Vi - massf/massi*Vi/xi;
                                filfacf = filfaci*pow(1./(1.-(deltaVol/Vi)),ncoll);
                                break;
                            }
                            case (5): // Garcia + mod compressive strenght kataoka
                            {    
                                double Vi = GrainVolumeSize(sizei,filfaci,rhos);
                                double Ebreak = 1.5*48./302.46*eroll;
                                double Ecomp = Ekin(massi,vrel) - (2.*massi-massf)*Ebreak/m0;

                                if (Ecomp<0) Ecomp=0;

                                double deltaVol = (Ecomp/(eroll*pow(filfaci/(maxpacking-filfaci)/a0,3)));
                                if (deltaVol > Vi)  deltaVol = Vi;
                                   

                                filfacf = filfaci*pow( 1./ (1. - 0.5*exp(1-pow(vrel/vfrag,2))*(deltaVol/Vi) ),ncoll);
                                break;
                            }
                        }
                    }
                    else    filfacf = filfaci;
                    break;
                }
            }
        }
        if (filfacf < filfacmincollgasgrav)
        {   filfacf = filfacmincollgasgrav;   }
        
        if (filfacf > filfacmax)
        {   filfacf = filfacmax;    }

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

double S4(const double& R, const double cparam, const double& mstar, const double& rhog,//->
          const double& cg, const double& a0, const double& rhos)
{
    return pow(rhog*cg/(rhos*a0*Omegak(R,mstar)),1.5)*sqrt(2.*(pow(2.,0.075)-1.)/cparam);
}

double S5(const double& R, const double cparam, const double& mstar, const double& rhog,//->
          const double& cg, const double& a0, const double& rhos)
{
    double nu = NuMolGas(rhog,cg);
    return sqrt(9.*nu*rhog/(2.*rhos*a0*a0*Omegak(R,mstar)))*pow(9.*nu*(pow(2.,0.2)-1.)/(cg*a0*cparam),1./6.);
}

double FilFacSGr(const double& sizef, const double& sizei, const double& filfaci,//->
              const double& rhos, const double& eroll, const double& vrel)
{
    double power;
    if (Ekin(sizei,filfaci,rhos,vrel)/(3.*b_oku*eroll) <= 1.)
    {   power = 3.*cratio/(1.-cratio);  }
    else
    {   power = -0.5;   }

    return filfaci*pow(sizef/sizei,power);
}

double FilFacSColl(const double& R, const double& mstar, const double& rhog, const double& cg, const double st, const double& sfrac,//->
                   const double& eroll, const double& a0, const double& rhos, const double& alpha, int& porreg, double& filfacpow)
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
                porreg = 1;
            }
            else
            {   if (sfrac <= S3(rhog,cg,a0))
                {
                    filfacpow = -1./3.;
                    filfaccoll = pow(s1,(1.+8.*cratio)/(3.*(1.-cratio)))*pow(sfrac,filfacpow);   // FilFac Ep-St<1
                    porreg = 2;
                }
                else
                {   
                    filfacpow = 0.;
                    filfaccoll = pow(s2,3.*cratio/(1.-cratio));   // FilFac St-St<1
                    porreg = 3;
                }
            }
        }
        else
        {
            if (sfrac <= s2)
            {
                filfacpow = 3.*cratio/(1.-cratio);
                filfaccoll = pow(sfrac,filfacpow);   // FilFac h&s
                porreg = 1;
            }
            else
            {
                filfacpow = 0.;
                filfaccoll = pow(s2,3.*cratio/(1.-cratio));   // FilFac St-St<1
                porreg = 3;
            }
        }
    }
    else
    {
        filfaccoll = pow((cparam*rhog*cg)/(2.*rhos*a0*Omegak(R,mstar)),0.25)/sqrt(sfrac);
        if (s4 <= s5)    
        {   
            filfaccoll *= pow(pow(2.,0.075)-1.,-0.25);      // FilFac Ep-St>1
            porreg = 4;
        }
        else             
        {   
            filfaccoll *= pow(pow(2.,0.2)-1.,-0.25);       // FilFac St-St>1
            porreg = 5;
        }
        filfacpow = -0.5;
    }
    return filfaccoll;
}

/* ------------------------ KATAOKA ------------------------*/

double FilFacSGas(const double& R, const double& mstar, const double& deltav, const double& st,const double& sizef, const double& eroll,//->
                  const double& a0, const double& rhos, const int& dragreg, double& filfacpow)//, const double& rhog, const double& cg)
{
    //double filfacgas;
    switch (dragreg)
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

double FilFacMinSColGasGrav(const double& R, const double& mstar, const double& rhog, const double& cg, const double& deltav,//->
                            const double st, const double& sizef, const double& rhos, const double& eroll, const double& a0,//->
                            const double& alpha, const int& dragreg, int& porreg, double& filfacpow)
{
    double filfacmin;
    double filfacpowgas;
    int maxindex;
    double filfaccoll = FilFacSColl(R,mstar,rhog,cg,st,sizef/a0,eroll,a0,rhos,alpha,porreg,filfacpow);
    double filfacgas  = FilFacSGas(R,mstar,deltav,st,sizef,eroll,a0,rhos,dragreg,filfacpowgas);
    double filfacgrav = FilFacSGrav(sizef,rhos,a0,eroll);

    double filfacminarray[] = {filfaccoll,filfacgas,filfacgrav};
    double filfacpowindex[] = {filfacpow,filfacpowgas,2.};
    int porregindex[] = {porreg,6,7};

    maxindex = distance(filfacminarray,max_element(filfacminarray,filfacminarray+3));
    filfacmin = filfacminarray[maxindex];
    filfacpow = filfacpowindex[maxindex];
    porreg = porregindex[maxindex];

    return filfacmin;
}

double FilFacSFinal(const double& R, const double& mstar, const double& rhog, const double& cg, const double& deltav, const double st,//->
                    const double& sizef, const double& sizei, const double& filfaci, const double& a0, const double& rhos, const double& eroll,//-> 
                    const double& alpha, const int& ifrag, const double& vfrag, const double& vrel, const int& dragreg, int& porreg, double& filfacpow)
{
    double filfacf = 1.;
    double filfacmincollgasgrav = FilFacMinSColGasGrav(R,mstar,rhog,cg,deltav,st,sizef,rhos,eroll,a0,alpha,dragreg,porreg,filfacpow);
    double filfacgr = FilFacSGr(sizef,sizei,filfaci,rhos,eroll,vrel);
    double filfacmax;

    if (sizef > a0)
    {
        if ((vrel < vfrag && ifrag > 0) || (ifrag == 0))
        {   filfacf = filfacgr; }
        else
        {   filfacf = filfaci;  
        
        // modify for grain compaction
        }

        if (filfacf <= filfacmincollgasgrav)
        {   filfacf = filfacmincollgasgrav;     }


        filfacmax = 0.5*maxpacking * (1 + sqrt(1 + 4*(1-maxpacking)/maxpacking/maxpacking * pow(sizef/a0,-3)));
        if (filfacf > filfacmax)
        {
            filfacf = filfacmax;
            filfacpow = 3.*cratio/(1.-cratio);
        }

        return filfacf;
    }
    else 
    {   return filfacf; }

}
