#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"
#include "../header/porosity.h"
#include "../header/evol.h"
#include "../header/disruption.h"
#include "../header/readwrite.h"

using namespace std;

// Presentation
void Presentation()
{
    cout << "\e[1m" << endl;
    cout << "!+-----                 -----+!" << endl << endl;
    cout << "      Welcome in EDEN !" << endl << endl;
    cout << "!+-----                 -----+!" << endl;
    cout << "\e[0m" << endl << endl;
}

// Terminal waiting animation
void Animation(const double& time, const double timeend, const double& RfRin, const double& sizelimsize, const bool& disrupted, const string& outputfile)
{
    cout << "\rProgression: " << setprecision(4) << time*100./(timeend) <<" %     "<< flush;

    if (time>=timeend || RfRin < 1. || sizelimsize > 1. || disrupted == true)
    {
        cout << "\rProgression: " << time*100./(timeend) <<" %             ";

        if (RfRin < 1.)
        {   cout << endl << "R < Rin: particle accreted               "<<flush; }
        if (sizelimsize > 1.)
        {   cout << endl << "Particle size bigger than max size       "<<flush;}
        if (disrupted == true)
        {   cout << endl << "Particle disrupted by spinning motion    "<<flush;}

        cout <<endl <<outputfile << " written" <<endl << endl;
    }
}

// Running time
void Runningtime(const double& runningtime)
{   cout << endl << "Running time: " << runningtime << " s" << endl;  }

// End
void End()
{   cout << "\e[1m" << endl;
    cout << "Job done ! " << endl;
    cout << "\e[0m" << endl;
}

/* ------------------------ MAIN ------------------------*/

int main()
{
    /*------------------------ INITIALS PARAMETERS FROM THE USER ------------------------*/

    int    massorsize;      // choose between dm/dt or ds/dt
    double tend;            // End time of the simulation
    int    stepmethod;      // Choose a time-stepping method
    double step;            // time step related to stepmethod
    int    profile;         // Compute disk profiles
    double Mstar;           // Star mass
    double Mdisk;           // Disk mass
    double Rin;             // Inner radius
    double Rout;            // Outer radius
    double R0;              // Reference radius
    double Rbump;           // Radius of the pressure bump (when ibump enabled)
    double dustfrac0;       // Initial dust to gas ratio at reference radius
    double dustfracmax;     // Max dust to gas ratio possible (when ibackreaction enabled)
    double H0R0;            // H/R at reference radius
    double p;               // p index
    double q;               // q index
    double alpha;           // Alpha turbulence Shakura & Sunyaev
    double bumpwidth;       // half width at half maximum in AU (when ibump enabled)
    double bumpheight;      // fraction of surface density
    double sizeini;         // Initial size
    double phiini;          // Initial filling factor
    double a0;              // Monomer size
    double rhos;            // Dust monomer density
    double youngmod0;       // Young Modulus of grains
    double esurf;           // Surface energy of grains
    double Yd0;             // Dynamic compression resistance constant
    double Ydpower;         // Dynamic compression resistance power law
    int    iporosity;       // Porous or compact grain option
    int    idrift;          // Drift option
    int    ibounce;         // Bounce option
    int    ifrag;           // Fragmention option
    int    ibump;           // Pressure bump option
    int    ibr;             // Back-reaction option
    int    idisrupt;        // Disruption by spinning motion
    double gammaft;         // Force-to-torque efficiency
    double vfragi;          // Initial fragmentation threshold (when fragmentation enabled)
    int    constvfrag;      // Constant vfrag option (when fragmentation enabled)
    double philim;          // Filling factor dynamic compression resistance limit (when fragmentation enabled)
    double philimbounce;    // Filling factor bounce limit (when bounce enabled)
    double limsize;         // Limit to the max size
    int    ngrains;         // Number of grains
    vector <double> Rini;   // Initials radii

    /*------------------------ INITIALS PARAMETERS TO COMPUTE ------------------------*/

    double h0;              // Disc height at R0
    double hg;              // Disc height at R
    double rhog0;           // Gas density at R0
    double rhog;            // Gas density at R
    double cg0;             // Gas sound speed at R0
    double cg;              // Gas sound speed at R
    double sigma0;          // Gas surface density at R0
    double sigma;           // Gas surface density at R
    double st;              // Stoke number
    double dustfrac;        // dust to gas ratio at R
    double vrel;            // Relative velocity between grains
    double vfrag = 0;       // Fragmentation threshold
    double vstick = 0;      // Sticking velocity
    double probabounce = 0; // Probability for a grain to bounce
    double ncoll = 0;       // number of collision per dt
    double eroll;           // Rolling energy


    /*------------------------ LOOP PARAMETERS ------------------------*/

    double massi;           // Mass before one loop
    double sizei;           // Size before one loop
    double phii;            // Filling factor before one loop
    double Ri;              // Radius before one loop
    double massf;           // Mass after one loop
    double sizef;           // Size after one loop
    double phif;            // Filling factor after one loop
    double Rf;              // Radius after one loop
    double dt;              // time step for loop
    double t = 0;           // time for loop
    double dmdt = 0;        // Variation of mass per unit of time
    double dsdt = 0;        // Variation of size per unit of time
    double drdt = 0;        // Variation of radius per unit of time (drift velocity)
    double deltav;          // Velocity difference between dust and gas
    int    ireg;            // Regime of the dust grain (growth, h&s, Ep-St<1, frag etc)
    double phipow;          // Power of the dominante filling factor to remove dsdt degeneracy
    bool disrupted = false; // Is grain disrupted by spinning motion

    double Rprofile;        // Radius to compute disk profiles

    string outputfile;      // Variable for the name of the output file

    clock_t t1, t2 = 0;     // sets of two variable to determine running time

    t1 = clock();


    /*------------------------ READ INPUT FILE ------------------------*/

    ReadFile(massorsize,tend,stepmethod,step,profile,Mstar,Mdisk,Rin,Rout,R0,dustfrac0,H0R0,p,q,alpha,iporosity,sizeini,
             phiini,a0,rhos,youngmod0,esurf,Yd0,Ydpower,idrift,ibounce,idisrupt,ifrag,ibr,ibump,gammaft,vfragi,constvfrag,
             philim,philimbounce,limsize,Rbump,dustfracmax,bumpwidth,bumpheight,ngrains,Rini);

    Presentation();

    /*------------------------ INITIALIZATION OF PARAMETERS AT R0 ------------------------*/


    h0 = H0R0*R0;
    sigma0 = Sigma0(Rin,Rout,R0,Mdisk,p,ibump,Rbump,bumpwidth,bumpheight);
    cg0 = Cg(R0,Mstar,h0);
    rhog0 = Rhog(sigma0,h0);
    eroll = Eroll(a0,esurf,youngmod0);

    /*------------------------ COMPUTE DISK QUANTITIES PROFILES ------------------------*/

    if (profile == 1)
    {
        ofstream writebump;
        writebump.open("diskprofiles.out");
        Rprofile = Rin;
        for (Rprofile = Rin; Rprofile <= Rout; Rprofile += 0.01)
        {   hg = Hg(Rprofile,q,R0,h0);
            cg = Cg(Rprofile,Mstar,hg);
            sigma = Sigma(Rprofile,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
            dustfrac = DustFrac(dustfrac0,dustfracmax,Rprofile,Rbump,bumpwidth,ibump);
            rhog = Rhog(sigma,hg);
            WriteProfileFile(writebump,Rprofile,hg,cg ,sigma,rhog,dustfrac,Pg(rhog,cg),T(Rprofile,q,R0,cg));
            Animation(Rprofile,Rout-0.01,2,0,false,"diskprofiles.out");
        }
        WriteProfileColumns();
    }

    /*------------------------ TIME LOOP FOR SIMULATION------------------------*/

    for (int j = 0; j < ngrains; j++)
    {
        // Set output file name
        ofstream writer;
        outputfile = FileName(massorsize,Rini[j],iporosity);
        writer.open(outputfile.c_str());

        // Initialize loop parameters
        dt = step;
        sizei = sizeini;
        sizef = sizeini;
        phii = phiini;
        phif = phiini;
        massi = GrainMass(sizeini,phiini,rhos);
        massf = GrainMass(sizeini,phiini,rhos);
        Ri = Rini[j];
        Rf = Rini[j];
        phipow = 3*cratio/(1-cratio);

        // Compute gas and dust quantities at t=0
        hg = Hg(Rf,q,R0,h0);
        sigma = Sigma(Rf,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
        dustfrac = DustFrac(dustfrac0,dustfracmax,Rf,Rbump,bumpwidth,ibump);
        cg = Cg(Rf,Mstar,hg);
        rhog = Rhog(sigma,hg);
        st = St(Ri,Mstar,rhog,cg,sizef,phiini,rhos,ireg);
        vrel = Vrel(cg,st,alpha);

        // Write in output file quantities at t=0
        WriteOutputFile(writer,t,Ri,massi,phii,sizei,st,cg,sigma,rhog,dustfrac,vrel,OmegaK(Ri,Mstar),0,0,ireg);

        // This is where the loop begin
        switch (massorsize)
        {
            case(0):
            {
                while(Rf > Rin && t < tend && sizef < limsize && disrupted != true)
                {
                    // Compute additionnal quantities for ifrag = 1 or 2, ibounce = 1 and idrift = 1
                    if (ifrag > 0)  vfrag = Vfrag(phif,philim,vfragi,constvfrag);

                    // Compute additionnal quantities for ibounce = 1
                    if (ibounce == 1)
                    {
                        vstick = Vstick(sizef,phif,rhos,esurf,youngmod0);
                        probabounce = ProbaBounce(phif,philimbounce,vrel,vstick,Vend(vstick));
                    }

                    // Compute additionnal quantities for idrift = 1: vdrift=drdt
                    if (idrift == 1)
                    {   drdt = DRDt(Rf,Mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,Rbump,bumpwidth,bumpheight); }

                    // Compute deltav
                    deltav = DeltaV(Rf,Mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,idrift,Rbump,bumpwidth,bumpheight);

                    // Compute dm/dt
                    dmdt = DmDt(sizef,rhog,dustfrac,vrel,ifrag,ibounce,vfrag,vstick,probabounce);

                    // Compute new dt and new time t
                    switch (stepmethod)
                    {   case (0):   break;
                        case (1):
                        {
                            dt = KepletDt(step,OmegaK(Rf,Mstar));
                            break;
                        }
                        case (2):
                        {   
                            dt = AdaptativeDt(t,tend,massorsize,ibump,massf,Rf,dmdt,drdt);
                            break;
                        }
                    }
                    t += dt;

                    // Compute new radius after dt if idrift=1
                    if (idrift == 1)    Rf += drdt*YearToSec(dt);

                    // Compute new mass after dt
                    massf += dmdt*YearToSec(dt);

                    // Check if new mass < monomer mass
                    if (massf < GrainMass(a0,1.,rhos))
                    {   massf = GrainMass(a0,1.,rhos);  }

                    // Compute filling factor for porous grains
                    if (iporosity == 1)
                    {
                        // Compute additionnal quantities for ibounce = 1
                        if (ibounce == 1)   ncoll = Ncoll(dt,Tcoll(sizef,rhog,phif,rhos,dustfrac,vrel));

                        // Compute the new filling factor after dt
                        phif = PhiMFinal(Ri,Mstar,rhog,cg,DeltaV2(Ri,Mstar,p,q,cg,st),st,massf,massi,phii,a0,rhos,eroll,alpha,ncoll,ifrag,
                                         ibounce,vfrag,vrel,vstick,probabounce,philim,Yd0,Ydpower,ireg);
                    }

                    // Compute gas and dust quantities at time t and radius R for a new loop

                    if (idrift == 1)   Ri = Rf;
                    massi = massf;
                    phii = phif;
                    sizef = GrainMassToSize(massf,phif,rhos);

                    hg = Hg(Rf,q,R0,h0);
                    sigma = Sigma(Rf,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
                    dustfrac = DustFrac(dustfrac0,dustfracmax,Rf,Rbump,bumpwidth,ibump);
                    cg = Cg(Rf,Mstar,hg);
                    rhog = Rhog(sigma,hg);
                    st = St(Rf,Mstar,rhog,cg,sizef,phif,rhos,ireg);
                    vrel = Vrel(cg,st,alpha);

                    // Write in output file quantities at time t
                    WriteOutputFile(writer,t,Rf,massf,phif,sizef,st,cg,sigma,rhog,dustfrac,vrel,OmegaK(Rf,Mstar),drdt,dmdt,ireg);

                    // Disruption by spinning motion
                    if (idisrupt == true)
                    {   disrupted = Disrupt(sizef,phif,rhos,deltav,gammaft,esurf,a0,st);   }
                    
                    // Waiting animation
                    Animation(t,tend,Rf/Rin,sizef/limsize,disrupted,outputfile);
                }
                break;
            }
            case (1):
            {
                while(Rf > Rin && t < tend && sizef < limsize)
                {
                    // Compute additionnal quantities for ifrag = 1 or 2, ibounce = 1 and idrift = 1
                    if (ifrag > 0)  vfrag = Vfrag(phif,philim,vfragi,constvfrag);

                    // Compute additionnal quantities for idrift = 1: vdrift=drdt
                    if (idrift == 1)
                    {   drdt = DRDt(Rf,Mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,Rbump,bumpwidth,bumpheight);  }

                    // Compute deltav
                    deltav = DeltaV(Rf,Mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,idrift,Rbump,bumpwidth,bumpheight);

                    // Compute ds/dt
                    dsdt = DsDt(phif,rhog,rhos,dustfrac,vrel,ifrag,vfrag,phipow);

                    // Compute new dt and new time t
                    switch (stepmethod)
                    {   case (0):   break;
                        case (1):
                        {
                            dt = KepletDt(step,OmegaK(Rf,Mstar));
                            break;
                        }
                        case (2):
                        {
                            dt = AdaptativeDt(t,tend,massorsize,ibump,sizef,Rf,dsdt,drdt);
                            break;
                        }
                    }
                    t += dt;

                    // Compute new radius after dt if idrift=1
                    if (idrift == 1)    Rf += drdt*YearToSec(dt);

                    // Compute new size after dt
                    sizef += dsdt*YearToSec(dt);

                    // Check if new size < monomer size
                    if (sizef < a0)  sizef = a0;

                    // Compute filling factor for porous grains
                    if (iporosity == 1)
                    {
                        // Compute the new filling factor after dt
                        phif = PhiSFinal(Ri,Mstar,rhog,cg,deltav,st,sizef,sizei,phii,a0,rhos,eroll,alpha,
                                         ifrag,vfrag,vrel,ireg,phipow);
                    }

                    // Compute gas and dust quantities at time t and radius R for a new loop
                    if (idrift == 1)   Ri = Rf;
                    sizei = sizef;
                    phii = phif;
                    massf = GrainMass(sizef,phif,rhos);

                    hg = Hg(Rf,q,R0,h0);
                    sigma = Sigma(Rf,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
                    dustfrac = DustFrac(dustfrac0,dustfracmax,Rf,Rbump,bumpwidth,ibump);
                    cg = Cg(Rf,Mstar,hg);
                    rhog = Rhog(sigma,hg);
                    st = St(Rf,Mstar,rhog,cg,sizef,phif,rhos,ireg);
                    vrel = Vrel(cg,st,alpha);

                    // Write in output file quantities at time t
                    WriteOutputFile(writer,t,Rf,massf,phif,sizef,st,cg,sigma,rhog,dustfrac,vrel,OmegaK(Rf,Mstar),drdt,dsdt,ireg);

                    // Disruption by spinning motion
                    if (idisrupt == true)
                    {   disrupted = Disrupt(sizef,phif,rhos,deltav,gammaft,esurf,a0,st);   }

                    // Waiting animation
                    Animation(t,tend,Rf/Rin,sizef/limsize,disrupted,outputfile);
                }
            break;
            }
        }
        // Reinitialize time for a new particle
        t = 0;
        dt = 1;
        if (idisrupt == true) disrupted = false;
        writer.close();
    }

    // Compute running time
    t2 = clock();
    Runningtime((t2-t1)/(1.*CLOCKS_PER_SEC));

    // Write initials conditions
    WriteInitFile(massorsize,tend,stepmethod,step,Mstar,Mdisk,Rin,Rout,R0,Rbump,dustfrac0,H0R0,p,q,alpha,iporosity,sizeini,
                  phiini,a0,rhos,idrift,ibounce,idisrupt,ifrag,ibr,ibump,vfragi,ngrains,sigma0,rhog0,cg0,
                  (t2-t1)/(1.*CLOCKS_PER_SEC));
    WriteOutputColumns(massorsize);
    End();

    return 0;
}
