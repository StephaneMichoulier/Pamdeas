#include <iostream>
#include <cmath>

#include "../header/constant.h"
#include "../header/general_functions.h"
#include "../header/porosity.h"

using namespace std;


/* ------------------------ SHARED FUNCTIONS ------------------------ */

double CParam(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& eroll,// ->
              const double& a0, const double& rhos, const double& alpha)
{
    return (243.*M_SQRT2*M_PI/15625.)*(Rossby*alpha*pow(a0,4.)*rhos*rhos*cg*OmegaK(R,Mstar)/(rhog*b_oku*eroll));
}

double Porosity(const double& phifinal)
{   return 1.-phifinal; }


/* ------------------------ BOUNCE ------------------------*/

double ProbaBounce(const double& phi, const double& philimbounce, const double& vrel, const double& vstick, const double& vend)
{
    if (phi >= philimbounce)
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

double VarVolumeBounce(const double& phi, const double& philim, const double& coeffrest, const double& ekin,// ->
                       const double& volume, const double& Yd0, const double& Ydpower)
{
    double compactvol = (1.-coeffrest*coeffrest)*ekin/Yd(phi,philim,Yd0,Ydpower);

    if (compactvol <= volume)
    {   return compactvol;    }
    else
    {   return volume;  }
}

double PhiBounce(const double& sizei, const double& phii, const double& rhos, const double& philim,// ->
                 const double& vrel, const double& vstick, const double& vyield, const double& ncoll,// ->
                 const double& Yd0, const double& Ydpower)
{
    double coeffrest = CoeffRest(vrel,vstick,vyield);
    double volume = GrainVolumeSize(sizei,phii,rhos);
    double ekin = Ekin(sizei,phii,rhos,vrel);
    double varvolbounce = VarVolumeBounce(phii,philim,coeffrest,ekin,volume,Yd0,Ydpower);
    double phibounce = phii*pow(1./(1.-(0.5*varvolbounce/volume)),ncoll);

    if (vrel <= vyield)
    {   return phii;    }
    else
    {
        if (phibounce <= 1.)
        {   return phibounce;   }
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

double M4(const double& R, const double cparam, const double& Mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos)
{
    return pow(rhog*cg/(rhos*a0*OmegaK(R,Mstar)),4.)*2.*(pow(2.,0.075)-1.)/cparam;
}

double M5(const double& R, const double cparam, const double& Mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos)
{
    double nu = NuMolGas(rhog,cg);
    return pow(9.*nu*rhog/(2.*rhos*a0*a0*OmegaK(R,Mstar)),1.5)*pow(9.*nu*(pow(2.,0.2)-1.)/(cg*a0*cparam),1./6.);
}

double PhiMGr(const double& massf, const double& massi, const double& phii, const double& rhos, const double& eroll, const double& vrel)
{
    double power;
    if (Ekin(massi,vrel)/(3.*b_oku*eroll) <= 1.)
    {   power = cratio;   }
    else
    {   power = -0.2; }

    return phii*pow(massf/massi,power);
}

/*double PhiMColl(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& massf,// ->
                const double eroll, const double& a0, const double& rhos, const double& alpha)
{
    // here, m1 to m5 and the associated functions M1 to M5 are normalized by m0 !
    double cparam = CParam(R,Mstar,rhog,cg,eroll,a0,rhos,alpha);
    double m = massf/GrainMass(a0,1.,rhos);
    double m1 = M1(cparam);
    double m2 = M2(cparam,rhog,cg,a0);
    double m3 = M3(m1,m2);
    double m4 = M4(R,cparam,Mstar,rhog,cg,a0,rhos);
    double m5 = M5(R,cparam,Mstar,rhog,cg,a0,rhos);
    double phicoll;

    if (m1 <= m2)
    {
        if (m <= m1)
        {
            phicoll = pow(m,cratio);   // Phi h&s
        }
        else
        {
            if (m3 <= m4)
            {
                if (m <= m3)
                {
                    phicoll = pow(m1,cratio+0.125)*pow(m,-0.125);   // Phi Ep-St<1
                }
                else
                {
                    if (m <= m5)
                    {
                        phicoll = pow(m2,cratio);   // Phi St-St<1
                    }
                    else
                    {
                        phicoll = pow(m2,cratio)*pow(m5/m,0.2);  // Phi St-St>1
                    }
                }
            }
            else
            {
                if (m <= m4)
                {
                    phicoll = pow(m1,cratio+0.125)*pow(m,-0.125);   // Phi Ep-St<1
                }
                else
                {
                    phicoll = pow(m1,cratio+0.125)*pow(m4,-0.125)*pow(m4/m,0.2);   // Phi Ep-St>1
                }
            }
        }
    }
    else
    {   if (m <= m2)
        {
            phicoll = pow(m,cratio);   // Phi h&s
        }
        else
        {
            if (m <= m5)
            {
                phicoll = pow(m2,cratio);   // Phi St-St<1
            }
            else
            {
                phicoll = pow(m2,cratio)*pow(m5/m,0.2);  // Phi St-St>1
            }
        }
    }
    return phicoll;
}*/

double PhiMColl(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& st,// ->
                const double& massf, const double eroll, const double& a0, const double& rhos, const double& alpha)
{
    // here, m1 to m5 and the associated functions M1 to M5 are normalized by m0 !
    double cparam = CParam(R,Mstar,rhog,cg,eroll,a0,rhos,alpha);
    double m = massf/GrainMass(a0,1.,rhos);
    double m1 = M1(cparam);
    double m2 = M2(cparam,rhog,cg,a0);
    double m4;
    double m5;
    double phicoll;

    if (st < 1.)
    {
        if (m1 <= m2)
        {
            if (m <= m1)
            {   phicoll = pow(m,cratio);    }   // Phi h&s
            else
            {   if (m <= M3(m1,m2))
                {
                    phicoll = pow(m1,cratio+0.125)*pow(m,-0.125);   // Phi Ep-St<1
                }
                else
                {
                    phicoll = pow(m2,cratio);   // Phi St-St<1
                }
            }
        }
        else
        {
            if (m <= m2)
            {   phicoll = pow(m,cratio);    }   // Phi h&s
            else
            {   phicoll = pow(m2,cratio);   }   // Phi St-St<1
        }
    }
    else
    {
        m4 = M4(R,cparam,Mstar,rhog,cg,a0,rhos);
        m5 = M5(R,cparam,Mstar,rhog,cg,a0,rhos);

        if (m4 <= m5)    phicoll = pow(m1,cratio+0.125)*pow(m4,-0.125)*pow(m4/m,0.2);   // Phi Ep-St>1
        else             phicoll = pow(m2,cratio)*pow(m5/m,0.2);  // Phi St-St>1
    }
    return phicoll;
}


/* ------------------------ KATAOKA ------------------------*/

/*double PhiMGas(const double& R, const double& Mstar, const double& p, const double& q, const double& rhog,// ->
               const double& cg, const double& massf, const double& rhos, const double& eroll, const double& a0)
{
    double m0 = GrainMass(a0,1.,rhos);
    double philimSt = pow(m0*rhog*DeltaV(R,Mstar,p,q,cg)*cg/(M_PI*eroll*rhos),1./3.);

    if  (TransRegEpSt(rhog,cg,GrainMassToSize(massf,philimSt,rhos)) < 1.)
    {   return philimSt;    }
    else
    {   return pow(6.*a0*a0*DeltaV(R,Mstar,p,q,cg)*NuMolGas(rhog,cg)*rhog/(eroll),0.375)*pow(massf/m0,-0.125); }
}*/

double PhiMGas(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& deltav, const double& st,// ->
               const double& massf, const double& rhos, const double& eroll, const double& a0, const int& iregime)
{

    double m0 = GrainMass(a0,1.,rhos);
    double phigas;
    switch (iregime)
    {
        case (1): case (3):
        {
            phigas = pow(m0*rhog*deltav*cg/(M_PI*eroll*rhos),1./3.);
            break;
        }
        case (2): case (4):
        {
            phigas = pow(6.*a0*a0*deltav*NuMolGas(rhog,cg)*rhog/(eroll),3./8.)*pow(massf/m0,-0.125);
            break;
        }
    }
    return phigas;
    //return pow((m0*a0*deltav*OmegaK(R,Mstar))/(eroll*M_PI*st),3./7.)*pow(massf/m0,1./7.);

}

double PhiMGrav(const double& massf, const double& rhos, const double& a0, const double& eroll)
{
    double m0 = GrainMass(a0,1.,rhos);
    return pow(G*m0*m0/(eroll*M_PI*a0),0.6)*pow(massf/m0,0.4) ;
}


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

double PhiMinMColGasGrav(const double& R, const double& Mstar, const double& rhog, const double cg, const double& deltav,// ->
                         const double st, const double& massf, const double& rhos, const double& eroll, const double& a0,// ->
                         const double& alpha, const int& iregime)
{
    double phimin;
    double phicoll = PhiMColl(R,Mstar,rhog,cg,st,massf,eroll,a0,rhos,alpha);
    double phigas  = PhiMGas(R,Mstar,rhog,cg,deltav,st,massf,rhos,eroll,a0,iregime);
    double phigrav = PhiMGrav(massf,rhos,a0,eroll);

    if (phigas < phicoll)
    {
        if (phigrav < phicoll)
        {   phimin = phicoll;   }
        else
        {   phimin = phigrav;   }
    }
    else
    {
        if (phigrav < phigas)
        {   phimin = phigas;    }
        else
        {   phimin = phigrav;   }
    }

    if (phimin > 1.)
    {   phimin = 1.;    }

    return phimin;
}

double PhiMFinal(const double& R, const double& Mstar, const double& rhog, const double& cg, const double deltav, const double st,// ->
                 const double& massf, const double& massi, const double& phii, const double& a0, const double& rhos, const double& eroll,// ->
                 const double& alpha, const double& ncoll, const int& ifrag, const int& ibounce, const double& vfrag, const double& vrel,// ->
                 const double& vstick, const double& probabounce, const double& philim, const double& Yd0, const double& Ydpower,// ->
                 const int& iregime)
{
    double phif = 1.;
    double sizei = GrainMassToSize(massi,phii,rhos);
    double phimincollgasgrav = PhiMinMColGasGrav(R,Mstar,rhog,cg,deltav,st,massf,rhos,eroll,a0,alpha,iregime);
    double phigr = PhiMGr(massf,massi,phii,rhos,eroll,vrel);

    if (sizei > a0)
    {
        switch (ibounce)
        {   case (0):
            {
                if ((vrel < vfrag && ifrag > 0) || (ifrag == 0 ))
                {
                    phif = phigr;
                    break;
                }
                else
                {
                    phif = phii;
                    break;
                }
            }
            case (1):
            {
                if ((vrel < vfrag && ifrag > 0) || (ifrag == 0 ))
                {
                    if (vrel <= vstick)
                    {
                        phif = phigr;
                        break;
                    }
                    else
                    {
                        double phibounce = PhiBounce(sizei,phii,rhos,philim,vrel,vstick,Vyield(vstick),ncoll,Yd0,Ydpower);

                        if (vrel < Vend(vstick))
                        {
                            if (phigr < phimincollgasgrav)
                            {   phif = phimincollgasgrav;   }
                            else
                            {   phif = phigr;   }

                            phif = phif*probabounce+(1.-probabounce)*phibounce;
                            break;
                        }
                        else
                        {
                            phif = phibounce;
                            break;
                        }
                    }
                }
                else
                {
                    phif = phii;
                    break;
                }
            }
        }
        if (phif < phimincollgasgrav)
        {   return phimincollgasgrav;   }
        else
        {   return phif;    }
    }
    else
    {  return phif; }
}


/* ------------------------ SIZE MODEL ------------------------ */

/* ------------------------ OKUZUMI GROWTH MODEL------------------------ */

double S1(const double& cparam)
{   return pow(cparam/(2.*(pow(2.,0.075)-1.)),(1.-cratio)/(1.+8.*cratio));  }

double S2(const double& cparam, const double& rhog, const double& cg, const double& a0)
{   return pow(cparam*cg*a0/(9.*NuMolGas(rhog,cg)*(pow(2.,0.2)-1.)),(1.-cratio)/(9.*cratio));  }

double S3(const double& rhog, const double& cg, const double& a0)
{   return (9.*NuMolGas(rhog,cg)/(2.*cg*a0)*((pow(2.,0.2)-1.)/(pow(2.,0.075)-1.)));    }

double S4(const double& R, const double cparam, const double& Mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos)
{
    return pow(rhog*cg/(rhos*a0*OmegaK(R,Mstar)),1.5)*sqrt(2.*(pow(2.,0.075)-1.)/cparam);
}

double S5(const double& R, const double cparam, const double& Mstar, const double& rhog,// ->
          const double& cg, const double& a0, const double& rhos)
{
    double nu = NuMolGas(rhog,cg);
    return sqrt(9.*nu*rhog/(2.*rhos*a0*a0*OmegaK(R,Mstar)))*pow(9.*nu*(pow(2.,0.2)-1.)/(cg*a0*cparam),1./6.);
}

double PhiSGr(const double& sizef, const double& sizei, const double& phii,// ->
              const double& rhos, const double& eroll, const double& vrel)
{
    double power;
    if (Ekin(sizei,phii,rhos,vrel)/(3.*b_oku*eroll) <= 1.)
    {   power = 3.*cratio/(1.-cratio);  }
    else
    {   power = -0.5;   }

    return phii*pow(sizef/sizei,power);
}

/*double PhiSColl(const double& R, const double& Mstar, const double& rhog, const double& cg, const double st, const double& sizef,// ->
                const double& eroll, const double& a0, const double& rhos, const double& alpha, const int& iregime, double& phipow)
{
    // here, s1 to s5 and the associated functions S1 to S5 are normalized by a0 !
    double cparam = CParam(R,Mstar,rhog,cg,eroll,a0,rhos,alpha);
    double s = sizef/a0;
    double s1 = S1(cparam);
    double s2 = S2(cparam,rhog,cg,a0);
    double s3 = S3(rhog,cg,a0);
    double s4 = S4(R,cparam,Mstar,rhog,cg,a0,rhos);
    double s5 = S5(R,cparam,Mstar,rhog,cg,a0,rhos);
    double phicoll;

    if (s1 <= s2)
    {
        if (s <= s1)
        {
            phipow = 3.*cratio/(1.-cratio);
            phicoll = pow(s,phipow);   // Phi h&s
        }
        else
        {
            if (s3 <= s4)
            {
                if (s <= s3)
                {
                    phipow = -1./3.;
                    phicoll = pow(s1,(1.+8.*cratio)/(3.*(1.-cratio)))*pow(s,phipow);   // Phi Ep-St<1
                }
                else
                {
                    if (s <= s5)
                    {
                        phipow = 0.;
                        phicoll = pow(s2,3.*cratio/(1.-cratio));   // Phi St-St<1
                    }
                    else
                    {
                        phicoll = pow(s2,3.*cratio/(1.-cratio))*sqrt(s5/s);  // Phi St-St>1
                        phipow = -0.5;
                    }
                }
            }
            else
            {
                if (s <= s4)
                {
                    phipow = -1./3.;
                    phicoll = pow(s1,(1.+8.*cratio)/(3.*(1.-cratio)))*pow(s,phipow);   // Phi Ep-St<1
                }
                else
                {
                    phicoll = pow(s1,(1.+8.*cratio)/(3.*(1.-cratio)))*pow(s4,-1./3.)*sqrt(s4/s);   // Phi Ep-St>1
                    phipow = -0.5;
                }
            }
        }
    }
    else
    {   if (s <= s2)
        {
            phipow = 3.*cratio/(1.-cratio);
            phicoll = pow(s,phipow);   // Phi h&s
        }
        else
        {
            if (s <= s5)
            {
                phipow = 0.;
                phicoll = pow(s2,3.*cratio/(1.-cratio));   // Phi St-St<1
            }
            else
            {
                phicoll = pow(s2,3.*cratio/(1.-cratio))*sqrt(s5/s);  // Phi St-St>1
                phipow = -0.5;
            }
        }
    }
    return phicoll;
}*/

double PhiSColl(const double& R, const double& Mstar, const double& rhog, const double& cg, const double st, const double& sizef,// ->
                const double& eroll, const double& a0, const double& rhos, const double& alpha, const int& iregime, double& phipow)
{
    // here, s1 to s5 and the associated functions S1 to S5 are normalized by a0 !
    double cparam = CParam(R,Mstar,rhog,cg,eroll,a0,rhos,alpha);
    double s = sizef/a0;
    double s1 = S1(cparam);
    double s2 = S2(cparam,rhog,cg,a0);
    double s4 = S4(R,cparam,Mstar,rhog,cg,a0,rhos);
    double s5 = S5(R,cparam,Mstar,rhog,cg,a0,rhos);
    double phicoll;

    if (st < 1.)
    {   if (s1 <= s2)
        {
            if (s <= s1)
            {
                phipow = 3.*cratio/(1.-cratio);
                phicoll = pow(s,phipow);   // Phi h&s
            }
            else
            {   if (s <= S3(rhog,cg,a0))
                {
                    phipow = -1./3.;
                    phicoll = pow(s1,(1.+8.*cratio)/(3.*(1.-cratio)))*pow(s,phipow);   // Phi Ep-St<1
                }
                else
                {   phipow = 0.;
                    phicoll = pow(s2,3.*cratio/(1.-cratio));   // Phi St-St<1

                }
            }
        }
        else
        {
            if (s <= s2)
            {
                phipow = 3.*cratio/(1.-cratio);
                phicoll = pow(s,phipow);   // Phi h&s
            }
            else
            {
                phipow = 0.;
                phicoll = pow(s2,3.*cratio/(1.-cratio));   // Phi St-St<1
            }
        }
    }
    else
    {
        phicoll = pow((cparam*rhog*cg)/(2*rhos*a0*OmegaK(R,Mstar)),0.25)/sqrt(s);
        if (s4 <= s5)    phicoll *= pow(pow(2.,0.075)-1.,-0.25);
        else             phicoll *= pow(pow(2.,0.2)-1.,-0.25);
        phipow = -0.5;
    }
    return phicoll;
}

/* ------------------------ KATAOKA ------------------------*/

double PhiSGas(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& deltav,// ->
               const double& sizef, const double& eroll, const double& a0, const int& iregime, double& phipow)
{
    double phigas;

    switch (iregime)
    {
        case (1): case (3):
        {
            phipow = 0.;
            phigas = pow(4.*a0*a0*a0*deltav*rhog*cg/(3.*eroll),1./3.);
            break;
        }
        case (2): case (4):
        {
            phipow = -1./3.;
            phigas = pow(6.*a0*a0*a0*deltav*NuMolGas(rhog,cg)*rhog/(sizef*eroll),1./3.);
            break;
        }
    }
    return phigas;
}

double PhiSGrav(const double& sizef, const double& rhos, const double& a0, const double& eroll)
{   return 16.*M_PI*G*rhos*rhos*a0*a0*a0*sizef*sizef/(9.*eroll);  }


/* ------------------------ FINAL FILLING FACTOR ------------------------*/

double PhiMinSColGasGrav(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& deltav,// ->
                         const double st, const double& sizef, const double& rhos, const double& eroll, const double& a0,// ->
                         const double& alpha, const int& iregime, double& phipow)
{
    double phimin;
    double powgas;
    double phicoll = PhiSColl(R,Mstar,rhog,cg,st,sizef,eroll,a0,rhos,alpha,iregime,phipow);
    double phigas  = PhiSGas(R,Mstar,rhog,cg,deltav,sizef,eroll,a0,iregime,powgas);
    double phigrav = PhiSGrav(sizef,rhos,a0,eroll);

    if (phigas < phicoll)
    {
        if (phigrav < phicoll)
        {   phimin = phicoll;   }
        else
        {
            phimin = phigrav;
            phipow = powgas;
        }
    }
    else
    {
        if (phigrav < phigas)
        {
            phimin = phigas;
            phipow = powgas;
        }
        else
        {
            phimin = phigrav;
            phipow = 2.;
        }
    }

    if (phimin > 1.)
    {
        phimin = 1.;
    }
    return phimin;
}

double PhiSFinal(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& deltav, const double st,// ->
                 const double& sizef, const double& sizei, const double& phii, const double& a0, const double& rhos, const double& eroll,// -> 
                 const double& alpha, const int& ifrag, const double& vfrag, const double& vrel, const int& iregime, double& phipow)
{
    double phif = 1.;
    double phimincollgasgrav = PhiMinSColGasGrav(R,Mstar,rhog,cg,deltav,st,sizef,rhos,eroll,a0,alpha,iregime,phipow);
    double phigr = PhiSGr(sizef,sizei,phii,rhos,eroll,vrel);

    if ((vrel < vfrag && ifrag > 0) || (ifrag == 0 ))
    {
        phif = phigr;
    }
    else
    {
        phif = phii;
    }

    if (phif <= phimincollgasgrav)
    {   return phimincollgasgrav;     }
    else
    {   return phif;    }
}
