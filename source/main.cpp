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


/* ------------------------ TERMINAL DISPLAY ------------------------*/

// Presentation animation
void PresentationAnim()
{
    cout << "\e[1m\n" ;
    cout << "!+-----                 -----+!\n";
    cout << "      Welcome in EDEN !\n";
    cout << "!+-----                 -----+!\n";
    cout << "\e[0m\n";
}

// Terminal progression / waiting animation
void ProgressionAnim(const double& iteratedvalue, const double endvalue, const int& istate, const string& outputfile)
{
    cout << "\rProgression: " << setprecision(4) << iteratedvalue*100./(endvalue) << " %     " << flush;

    if (iteratedvalue >= endvalue || istate != 0)
    {
        cout << "\rProgression: " << iteratedvalue*100./(endvalue) << " %             ";

        switch (istate)
        {   
            case(1):    
            {   
                cout << "\nGrain size bigger than max size       " << flush;
                break;
            }     
            case(2):    
            {
                cout << "\nR < Rin: grain accreted               " << flush;
                break;
            }
            case(3):    
            {   
                cout << "\nGrain disrupted by spinning motion    " << flush;
                break;
            }
        }

        cout << "\n" << outputfile << " written\n\n";
    }
}

// Running time
void Runningtime(const double& runningtime)
{   cout << "\nRunning time: " << runningtime << " s\n";  }

// End animation
void EndAnim()
{   cout << "\e[1m\n";
    cout << "Job done !\n";
    cout << "\e[0m\n";
}


/* ------------------------------------------------------*/
/* ------------------------ MAIN ------------------------*/
/* ------------------------------------------------------*/

int main(int argc, char* argv[])
{   
    if (argc == 1){
    /*------------------------ INITIALS PARAMETERS FROM THE USER ------------------------*/

    int    massorsize;      // choose between dm/dt or ds/dt
    double tend;            // End time of the simulation   [yr]
    int    stepmethod;      // Choose a time-stepping method
    double step;            // time step related to stepmethod
    int    profile;         // Compute disk profiles
    double mstar;           // Star mass    [Msol]
    double mdisk;           // Disk mass    [Msol]
    double Rin;             // Inner radius [AU]
    double Rout;            // Outer radius [AU]
    double R0;              // Reference radius [AU]
    double Rbump;           // Radius of the pressure bump (when ibump enabled) [AU]
    double dustfrac0;       // Initial dust to gas ratio at reference radius
    double dustfracmax;     // Max dust to gas ratio possible (when ibackreaction enabled)
    double h0R0;            // H/R at reference radius
    double p;               // p index
    double q;               // q index
    double alpha;           // Alpha turbulence Shakura & Sunyaev
    double bumpwidth;       // half width at half maximum (when ibump enabled)  [AU]
    double bumpheight;      // fraction of surface density
    double sizeini;         // Initial size [m]
    double filfacini;       // Initial filling factor
    double a0;              // Monomer size [m]
    double rhos;            // Dust monomer density [kg/m³]
    double youngmod0;       // Young Modulus of grains [PA]
    double esurf;           // Surface energy of grains [J/m²]
    double Yd0;             // Dynamic compression resistance constant  [Pa]
    double Ydpower;         // Dynamic compression resistance power law
    int    iporosity;       // Porous or compact grain option
    int    idrift;          // Drift option
    int    ibounce;         // Bounce option
    int    ifrag;           // Fragmention option
    int    ibump;           // Pressure bump option
    int    ibr;             // Back-reaction option
    int    idisrupt;        // Disruption by spinning motion
    double gammaft;         // Force-to-torque efficiency
    double vfragi;          // Initial fragmentation threshold (when fragmentation enabled) [m/s]
    int    constvfrag;      // Constant vfrag option (when fragmentation enabled)
    double filfaclim;       // Filling factor dynamic compression resistance limit (when fragmentation enabled)
    double filfacbnc;       // Filling factor bounce limit (when bounce enabled)
    double limsize;         // Limit to the max size    [m]
    int    ngrains;         // Number of grains
    vector <double> Rini;   // Initials radii

    /*------------------------ INITIALS PARAMETERS TO COMPUTE ------------------------*/

    double h0;              // Disc height at R0    [AU]
    double hg;              // Disc height at R     [AU]
    double rhog0;           // Gas density at R0    [kg/m³]
    double rhog;            // Gas density at R     [kg/m³]
    double cg0;             // Gas sound speed at R0    [m/s]
    double cg;              // Gas sound speed at R     [m/s]
    double sigma0;          // Gas surface density at R0    [kg/m²]
    double sigma;           // Gas surface density at R     [kg/m²]
    double st;              // Stoke number
    double dustfrac;        // dust to gas ratio at R
    double vrel;            // Relative velocity between grains [m/s]
    double vfrag;           // Fragmentation threshold  [m/s]
    double vstick;          // Sticking velocity    [m/s]
    double probabounce;     // Probability for a grain to bounce
    double ncoll;           // number of collision per dt   [s⁻¹]
    double eroll;           // Rolling energy   [J]


    /*------------------------ LOOP PARAMETERS ------------------------*/

    double massi;           // Mass before one loop [kg]
    double sizei;           // Size before one loop [m]
    double filfaci;         // Filling factor before one loop
    double Ri;              // Radius before one loop   [AU]
    double massf;           // Mass after one loop  [kg]
    double sizef;           // Size after one loop  [m]
    double filfacf;         // Filling factor after one loop
    double Rf;              // Radius after one loop    [AU]
    double dt;              // time step for loop   [yr]
    double t;               // time for loop    [ye]
    double dmdt;            // Variation of mass per unit of time   [kg/s]
    double dsdt;            // Variation of size per unit of time   [m/s]
    double drdt;            // Variation of radius per unit of time (drift velocity)    [AU/s]
    double deltav;          // Velocity difference between dust and gas     [m/s]
    int    ireg;            // Regime of the dust grain (growth, h&s, Ep-St<1, frag etc)
    double filfacpow;       // Power of the dominante filling factor to remove dsdt degeneracy
    bool   disrupted;       // Is grain disrupted by spinning motion
    vector <int> istate;    // State of the grain: 0=alive, 1=maxsize, 2=accreted,3=disrupted 

    double Rprofile;        // Radius to compute disk profiles

    string outputfile;      // Variable for the name of the output file

    clock_t t1, t2 = 0.;    // sets of two variable to determine running time

    t1 = clock();

    PresentationAnim();


    /*------------------------ READ INPUT FILE ------------------------*/

    cout << "Reading input file\n";

    ReadFile(massorsize,tend,stepmethod,step,profile,mstar,mdisk,Rin,Rout,R0,dustfrac0,h0R0,p,q,alpha,iporosity,sizeini,
             filfacini,a0,rhos,youngmod0,esurf,Yd0,Ydpower,idrift,ibounce,idisrupt,ifrag,ibr,ibump,gammaft,vfragi,constvfrag,
             filfaclim,filfacbnc,limsize,Rbump,dustfracmax,bumpwidth,bumpheight,ngrains,Rini,istate);

    cout << "\r\rInput file read" << flush << "\n\n";

    /*------------------------ INITIALIZATION OF PARAMETERS AT R0 ------------------------*/

    h0 = h0R0*R0;
    sigma0 = Sigma0(Rin,Rout,R0,mdisk,p,ibump,Rbump,bumpwidth,bumpheight);
    cg0 = Cg(R0,mstar,h0);
    rhog0 = Rhog(sigma0,h0);
    eroll = Eroll(a0,esurf,youngmod0);


    /*------------------------ COMPUTE DISK QUANTITIES PROFILES ------------------------*/

    if (profile == 1)
    {
        cout << "Compute disk profiles\n";
        ofstream writebump;
        writebump.open("disk_profiles.out");
        Rprofile = Rin;
        for (Rprofile = Rin; Rprofile <= Rout; Rprofile += 0.01)
        {   hg = Hg(Rprofile,q,R0,h0);
            cg = Cg(Rprofile,mstar,hg);
            sigma = Sigma(Rprofile,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
            dustfrac = DustFrac(dustfrac0,dustfracmax,Rprofile,Rbump,bumpwidth,ibump);
            rhog = Rhog(sigma,hg);
            WriteProfileFile(writebump,Rprofile,hg,cg ,sigma,rhog,dustfrac,Pg(rhog,cg),T(Rprofile,q,R0,cg));
            ProgressionAnim(Rprofile,Rout-0.01,0,"diskprofiles.out");
        }
        WriteProfileHeader();
        writebump.close();
        cout << "\n";
    }


    /*------------------------ OPEN FILE FOR DISRUPTION PARAM ------------------------*/

    ofstream writerdisruption;
    if (idisrupt == 1)  
    {   
        writerdisruption.open(DisruptFileName(massorsize).c_str());
        WriteDisruptHeader();
    }


    /*------------------------ TIME LOOP FOR SIMULATION ------------------------*/
   
    cout << "\e[1m" << "Simulation starting:" << "\e[0m\n\n";

    for (int j = 0; j < ngrains; j++)
    {
        // Set output file name
        ofstream writer;
        outputfile = OutputFileName(massorsize,Rini[j],iporosity);
        writer.open(outputfile.c_str());

        // Initialize loop parameters
        if (stepmethod == 0) dt = step;
        t = 0;
        sizei = sizeini;
        sizef = sizeini;
        filfaci = filfacini;
        filfacf = filfacini;
        massi = GrainMass(sizeini,filfacini,rhos);
        massf = GrainMass(sizeini,filfacini,rhos);
        Ri = Rini[j];
        Rf = Rini[j];
        filfacpow = 3.*cratio/(1.-cratio);
        disrupted = false;

        // Compute gas and dust quantities at t=0
        hg = Hg(Rf,q,R0,h0);
        sigma = Sigma(Rf,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
        dustfrac = DustFrac(dustfrac0,dustfracmax,Rf,Rbump,bumpwidth,ibump);
        cg = Cg(Rf,mstar,hg);
        rhog = Rhog(sigma,hg);
        st = St(Ri,mstar,rhog,cg,sizei,filfacini,rhos,0.,ireg); //we assume deltav very small for small grain strongly coupled with gas
        deltav = DeltaV(Rf,mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,idrift,Rbump,bumpwidth,bumpheight);
        vrel = Vrel(cg,st,alpha);

        // Write in output file quantities at t=0
        WriteOutputFile(writer,t,Ri,massi,filfaci,sizei,st,cg,sigma,rhog,dustfrac,vrel,Omegak(Ri,mstar),0.,0.,ireg);

        // This is where the loop begin
        switch (massorsize)
        {
            case(0):
            {
                while(t < tend && istate[j] == 0)
                {
                    // Compute additionnal quantities for ifrag = 1 or 2, ibounce = 1 and idrift = 1
                    if (ifrag > 0)  vfrag = Vfrag(filfaci,filfaclim,vfragi,constvfrag);

                    // Compute additionnal quantities for ibounce = 1
                    if (ibounce == 1)
                    {
                        vstick = Vstick(sizef,filfaci,rhos,esurf,youngmod0);
                        probabounce = ProbaBounce(filfaci,filfacbnc,vrel,vstick,Vend(vstick));
                    }

                    // Compute additionnal quantities for idrift = 1: vdrift=drdt
                    if (idrift == 1)
                    {   drdt = DRDt(Ri,mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,Rbump,bumpwidth,bumpheight); }

                    // Compute dm/dt
                    dmdt = DmDt(sizef,rhog,dustfrac,vrel,ifrag,ibounce,vfrag,vstick,probabounce);

                    // Compute new dt and new time t
                    switch (stepmethod)
                    {   
                        case (1):
                        {
                            dt = KeplerDt(step,Omegak(Ri,mstar));
                            break;
                        }
                        case (2):
                        {   
                            dt = AdaptativeDt(t,tend,massorsize,massi,Ri,dmdt,drdt);
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
                        if (ibounce == 1)   ncoll = Ncoll(dt,Tcoll(sizef,rhog,filfaci,rhos,dustfrac,vrel));

                        // Compute the new filling factor after dt
                        filfacf = FilFacMFinal(Ri,mstar,rhog,cg,deltav,st,massf,massi,filfaci,a0,rhos,eroll,alpha,ncoll,ifrag,
                                         ibounce,vfrag,vrel,vstick,probabounce,filfaclim,Yd0,Ydpower);
                    }

                    // Compute gas and dust quantities at time t and radius R for a new loop

                    if (idrift == 1)   Ri = Rf;
                    massi = massf;
                    filfaci = filfacf;
                    sizef = GrainMassToSize(massf,filfacf,rhos);

                    hg = Hg(Rf,q,R0,h0);
                    sigma = Sigma(Rf,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
                    dustfrac = DustFrac(dustfrac0,dustfracmax,Rf,Rbump,bumpwidth,ibump);
                    cg = Cg(Rf,mstar,hg);
                    rhog = Rhog(sigma,hg);
                    st = St(Rf,mstar,rhog,cg,sizef,filfacf,rhos,deltav,ireg);
                    deltav = DeltaV(Rf,mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,idrift,Rbump,bumpwidth,bumpheight);
                    vrel = Vrel(cg,st,alpha);

                    // Disruption by spinning motion
                    if (idisrupt == 1)
                    {   disrupted = Disrupt(sizef,filfacf,rhos,deltav,gammaft,esurf,a0);   }

                    //State of the grain
                    State(istate[j],Rf/Rin,sizef/limsize,disrupted);

                    if (istate[j] == 3)
                    {   
                        WriteDisruptFile(writerdisruption,Rf,massf,filfacf,sizef,st,vrel,FreqSpin(sizef,deltav,gammaft),
                        TensileStess(sizef,filfacf,rhos,deltav,gammaft),gammaft,alpha,a0);
                    }

                    // Write in output file quantities at time t
                    WriteOutputFile(writer,t,Rf,massf,filfacf,sizef,st,cg,sigma,rhog,dustfrac,vrel,Omegak(Rf,mstar),drdt,dmdt,ireg);

                    // Waiting animation
                    ProgressionAnim(t,tend,istate[j],outputfile);
                }
                break;
            }
            case (1):
            {
                while(t < tend && istate[j] == 0)
                {
                    // Compute additionnal quantities for ifrag = 1 or 2 
                    if (ifrag > 0)  vfrag = Vfrag(filfaci,filfaclim,vfragi,constvfrag);

                    // Compute additionnal quantities for idrift = 1: vdrift=drdt
                    if (idrift == 1)
                    {   drdt = DRDt(Ri,mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,Rbump,bumpwidth,bumpheight);  }

                    // Compute ds/dt
                    dsdt = DsDt(filfaci,rhog,rhos,dustfrac,vrel,ifrag,vfrag,filfacpow);

                    // Compute new dt and new time t
                    switch (stepmethod)
                    {
                        case (1):
                        {
                            dt = KeplerDt(step,Omegak(Ri,mstar));
                            break;
                        }
                        case (2):
                        {
                            dt = AdaptativeDt(t,tend,massorsize,sizef,Ri,dsdt,drdt);
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
                        filfacf = FilFacSFinal(Ri,mstar,rhog,cg,deltav,st,sizef,sizei,filfaci,a0,rhos,eroll,alpha,
                                         ifrag,vfrag,vrel,ireg,filfacpow);
                    }

                    // Compute gas and dust quantities at time t and radius R for a new loop
                    if (idrift == 1)   Ri = Rf;
                    sizei = sizef;
                    filfaci = filfacf;
                    massf = GrainMass(sizef,filfacf,rhos);

                    hg = Hg(Rf,q,R0,h0);
                    sigma = Sigma(Rf,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
                    dustfrac = DustFrac(dustfrac0,dustfracmax,Rf,Rbump,bumpwidth,ibump);
                    cg = Cg(Rf,mstar,hg);
                    rhog = Rhog(sigma,hg);
                    st = St(Rf,mstar,rhog,cg,sizef,filfacf,rhos,deltav,ireg);
                    deltav = DeltaV(Rf,mstar,p,q,rhog,cg,R0,sigma0,h0,dustfrac,st,alpha,ibr,ibump,idrift,Rbump,bumpwidth,bumpheight);
                    vrel = Vrel(cg,st,alpha);

                    // Disruption by spinning motion
                    if (idisrupt == 1)
                    {   disrupted = Disrupt(sizef,filfacf,rhos,deltav,gammaft,esurf,a0);   }

                   //State of the grain
                    State(istate[j],Rf/Rin,sizef/limsize,disrupted);

                    if (istate[j] == 3)
                    {   
                        WriteDisruptFile(writerdisruption,Rf,massf,filfacf,sizef,st,vrel,FreqSpin(sizef,deltav,gammaft),
                        TensileStess(sizef,filfacf,rhos,deltav,gammaft),gammaft,alpha,a0);
                    }

                    // Write in output file quantities at time t
                    WriteOutputFile(writer,t,Rf,massf,filfacf,sizef,st,cg,sigma,rhog,dustfrac,vrel,Omegak(Rf,mstar),drdt,dsdt,ireg);

                    // Waiting animation
                    ProgressionAnim(t,tend,istate[j],outputfile);
                }
            break;
            }
        }
        writer.close();
    }
    writerdisruption.close();

    // Compute running time
    t2 = clock();
    Runningtime((t2-t1)/(1.*CLOCKS_PER_SEC));

    // Write initials conditions
    WriteInitFile(massorsize,tend,stepmethod,step,mstar,mdisk,Rin,Rout,R0,Rbump,dustfrac0,h0R0,p,q,alpha,iporosity,sizeini,
                  filfacini,a0,rhos,idrift,ibounce,ifrag,ibr,ibump,idisrupt,vfragi,ngrains,sigma0,rhog0,cg0,istate,
                  (t2-t1)/(1.*CLOCKS_PER_SEC));
    WriteOutputHeader(massorsize);
    EndAnim();

    }


    /*------------------------ WRITE INPUT FILE ------------------------*/

    else if (argc == 2)
    {
        if (strcmp(argv[1], "setup") == 0)
        {   
            WriteInputFile();
            cout << "Input file has been written: input.in" << endl;
        }
        else
        {   cerr << "Error: unknown argument passed" << endl;   }
    }
    else
    {   cerr << "Error: to many arguments passed" << endl;  }

    return 0;
}