#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <chrono>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"
#include "../header/porosity.h"
#include "../header/evol.h"
#include "../header/disruption.h"
#include "../header/readwrite.h"

using namespace std;


/* ------------------------ TERMINAL DISPLAY ------------------------*/

void checksize(double& mass, double& filfac,const double& rhos, double& sdustmin)
{
    // Used to test functions for Phantom when grainsizemin != a_0
    double sdust = GrainMassToSize(mass,filfac,rhos);
    if (sdust < sdustmin)
    mass = mass * pow(sdustmin/sdust,3);
}

// Presentation animation
void PresentationAnim()
{
    cout << "\e[1m" << endl;
    cout << "!+-----                 -----+!" << endl;
    cout << "      Welcome in PAMPEAS !" << endl;
    cout << "!+-----                 -----+!" << endl;
    cout << "\e[0m" << endl;
}

// Terminal progression / waiting animation
void ProgressionAnim(const double& iteratedvalue, const double endvalue, const int& istate, const string& outputfile)
{
    cout << "\rProgression: " << setprecision(4) << iteratedvalue*100./(endvalue) << " %     " << flush;

    if (iteratedvalue >= endvalue || istate != 0)
    {
        cout << "\rProgression: " << iteratedvalue*100./(endvalue) << " %             " << endl;

        switch (istate)
        {   
            case(1):    
            {   
                cout << "Grain size bigger than max size       " << endl;
                break;
            }     
            case(2):    
            {
                cout << "R < Rin: grain accreted               " << endl;
                break;
            }
            case(3):    
            {   
                cout << "Grain disrupted by spinning motion    " << endl;
                break;
            }
        }

        cout << outputfile << " written\n" << endl;
    }
}

// Wall time
void Walltime(const double& walltime)
{   cout << "\nWall time: " << walltime << " s" << endl;    }

// End animation
void EndAnim()
{   cout << "\e[1m" << endl;
    cout << "Job done !" << endl;
    cout << "\e[0m" << endl;
}


/* ------------------------------------------------------*/
/* -------------------- MAIN ROUTINE --------------------*/
/* ------------------------------------------------------*/
    
void Pamdeas(const string& input)
{
    /*------------------------ INITIALS PARAMETERS FROM THE USER ------------------------*/

    int    massorsize;      // choose between dm/dt or ds/dt
    double tend;            // End time of the simulation [yr]
    int    stepmethod;      // Choose a time-stepping method
    double step;            // time step related to stepmethod
    int    profile;         // Compute disc profiles
    double mstar;           // Star mass [Msol]
    double mdisc;           // Disc mass [Msol]
    double sigma0;          // Gas surface density at R0 [kg/m²]
    double Rin;             // Inner radius [AU]
    double Rout;            // Outer radius [AU]
    double R0;              // Reference radius [AU]
    double Rbump;           // Radius of the pressure bump (when ibump is enabled) [AU]
    double Rsnow;           // Radius of the snow line (when isnow is enabled)
    double dustfrac0;       // Initial dust to gas ratio at reference radius
    double dustfracmax;     // Max dust to gas ratio possible (when ibackreaction is enabled)
    double hg0R0;            // H/R at reference radius
    double T0;              // Temperature at R0 [K]
    double p;               // p index
    double q;               // q index
    double alpha;           // Alpha turbulence Shakura & Sunyaev
    double bumpwidth;       // half width at half maximum (when ibump is enabled) [AU]
    double bumpheight;      // fraction of surface density
    double sizeini;         // Initial size [m]
    double a0;              // Monomer size [m]
    double rhos;            // Dust monomer density [kg/m³]
    double youngmod0;       // Young Modulus of grains [PA]
    double esurf;           // Surface energy of grains [J/m²]
    double Yd0;             // Dynamic compression resistance constant [Pa]
    double Ydpower;         // Dynamic compression resistance power law
    int    isetdens;        // Set density profile option
    int    isettemp;        // Set temperature profile option
    int    iporosity;       // Porous or compact grain option
    int    idrift;          // Drift option
    int    ibounce;         // Bounce option
    int    ifrag;           // Fragmention option
    int    ieros;           // Erosion option
    int    icomp;           // Compaction option
    int    ibump;           // Pressure bump option
    int    ibr;             // Back-reaction option
    int    idisrupt;        // Disruption by spinning motion
    int    isnow;           // Snow line option
    double vfragi;          // Initial fragmentation threshold (when ifrag is enabled) [m/s]
    double vfragin;         // Inward fragmentation threshold (when ifrag & isnow is enabled) [m/s]
    double vfragout;        // Outward fragmentation threshold (when ifrag & isnow is enabled) [m/s]
    double verosi;          // Initial erosion threshold (when ieros is enabled) [vfrag]
    int    constvfrag;      // Constant vfrag option (when fragmentation is enabled)
    double filfaclim;       // Filling factor dynamic compression resistance limit (when ifrag is enabled)
    double filfacbnc;       // Filling factor bounce limit (when ibounce is enabled)
    double maxsize;         // Limit to the max size [m]
    double gammaft;         // Force-to-torque efficiency
    int    disrupteq;       // Disruption equation (0=Tatsuuma2019, 1=Kimura2020)
    double weirmod;         // Weirbull modulus if disrupteq=1 (Kimura et al. 2020)
    int    ngrains;         // Number of grains
    vector <double> Rini;   // Initials radii

    /*------------------------ INITIALS PARAMETERS TO COMPUTE ------------------------*/

    double hg0;             // Disc height at R0 [AU]
    double hg;              // Disc height at R [AU]
    double rhog0;           // Gas density at R0 [kg/m³]
    double rhog;            // Gas density at R [kg/m³]
    double cg0;             // Gas sound speed at R0 [m/s]
    double cg;              // Gas sound speed at R [m/s]
    double sigma;           // Gas surface density at R [kg/m²]
    double st;              // Stoke number
    double dustfrac;        // dust to gas ratio at R
    double vrel;            // Relative velocity between grains [m/s]
    double vfrag;           // Fragmentation threshold [m/s]
    double veros;           // Initial erosion threshold [m/s]
    double vstick;          // Sticking velocity [m/s]
    double probabounce;     // Probability for a grain to bounce
    double ncoll;           // number of collision per dt [s⁻¹]
    double eroll;           // Rolling energy [J]


    /*------------------------ LOOP PARAMETERS ------------------------*/

    double massi;           // Mass before one loop [kg]
    double sizei;           // Size before one loop [m]
    double filfaci;         // Filling factor before one loop
    double Ri;              // Radius before one loop [AU]
    double massf;           // Mass after one loop [kg]
    double sizef;           // Size after one loop [m]
    double filfacf;         // Filling factor after one loop
    double Rf;              // Radius after one loop [AU]
    double dt;              // time step for loop [yr]
    double t;               // time for loop [yr]
    double tlastwrite;      // time since last write in output to not write ridiculous amount of unecessary data [yr]
    double dmdt;            // Variation of mass per unit of time [kg/s]
    double dsdt;            // Variation of size per unit of time [m/s]
    double drdt;            // Variation of radius per unit of time (drift velocity) [AU/s]
    double deltav;          // Velocity difference between dust and gas [m/s]
    int    dragreg;         // Drag regime of the dust grain (1=Ep, 2=St lin, 3=St nlin, 4=St quad)
    int    porreg;          // Porous expansion/compression regime (1=h&s, 2=EpSt<1, 3=StSt<1, 4=EpSt>1, 5=StSt>1, 6=Gas, 7=Grav)
    double filfacpow;       // Power of the dominante filling factor to remove dsdt degeneracy
    bool   disrupted;       // Is grain disrupted by spinning motion
    vector <int> istate;    // State of the grain: 0=alive, 1=maxsize, 2=accreted, 3=disrupted 

    double Rprofile;        // Radius to compute disc profiles

    string outputfile;      // Variable for the name of the output file

    auto wtbegin = chrono::high_resolution_clock::now(); // start counter for walltime

    PresentationAnim();


    /*------------------------ READ INPUT FILE ------------------------*/

    cout << "Reading input file" << endl;

    ReadFile(massorsize,tend,stepmethod,step,profile,isetdens,isettemp,Rin,Rout,R0,mstar,mdisc,sigma0,hg0R0,T0,dustfrac0,p,q,alpha,ibr,ibump,Rbump,
             dustfracmax,bumpwidth,bumpheight,iporosity,sizeini,a0,rhos,idrift,ibounce,idisrupt,ifrag,vfragi,ieros,verosi,icomp,maxsize,isnow,
             Rsnow,vfragin,vfragout,youngmod0,esurf,Yd0,Ydpower,constvfrag,filfaclim,filfacbnc,gammaft,disrupteq,weirmod,ngrains,Rini,istate,input);

    cout << "Input file read\n" << endl;

    if (ngrains == 0 && profile == 0)
    {   cout << "Too bad, nothing to compute" << endl;  }


    /*------------------------ INITIALIZATION OF PARAMETERS AT R0 ------------------------*/

    if (isetdens == 0)  sigma0 = Sigma0(Rin,Rout,R0,mdisc,p,ibump,Rbump,bumpwidth,bumpheight);
    else                mdisc = Mdisc(Rin,Rout,R0,sigma0,p,ibump,Rbump,bumpwidth,bumpheight);

    if (isettemp == 0)  
    {   
        hg0 = hg0R0*R0;
        T0 = T(R0,mstar,q,R0,hg0);
    }
    else                
    {   hg0 = Hg(R0,mstar,Cg(T0));   }

    cg0 = Cg(R0,mstar,hg0);
    rhog0 = Rhog(sigma0,hg0);
    eroll = Eroll(a0,esurf,youngmod0);


    /*------------------------ COMPUTE DISK QUANTITIES PROFILES ------------------------*/

    if (profile == 1)
    {
        cout << "Compute disc profiles" << endl;
        ofstream writebump;
        writebump.open("disc_profiles.out");
        WriteProfileHeader(writebump);
        Rprofile = Rin;
        for (Rprofile = Rin; Rprofile <= Rout; Rprofile += 0.01)
        {   hg = Hg(Rprofile,q,R0,hg0);
            cg = Cg(Rprofile,mstar,hg);
            sigma = Sigma(Rprofile,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
            dustfrac = DustFrac(dustfrac0,dustfracmax,Rprofile,Rbump,bumpwidth,ibump);
            rhog = Rhog(sigma,hg);
            WriteProfileFile(writebump,Rprofile,hg,cg ,sigma,rhog,dustfrac,Pg(rhog,cg),T(cg));
            ProgressionAnim(Rprofile,Rout-0.01,0,"disc_profiles.out");
        }
        writebump.close();
    }


    /*------------------------ OPEN FILE FOR DISRUPTION PARAM ------------------------*/

    ofstream writerdisruption;
    if (idisrupt == 1 && ngrains > 0)  
    {   
        writerdisruption.open(DisruptFileName(massorsize).c_str());
        WriteDisruptHeader(writerdisruption);
    }


    /*------------------------ TIME LOOP FOR SIMULATION ------------------------*/
   
    if (ngrains > 0)
    {   cout << "\n\e[1m" << "Simulation starting:" << "\e[0m\n" << endl;   }

    for (int j = 0; j < ngrains; j++)
    {
        // Set output file name
        ofstream writer;
        outputfile = OutputFileName(massorsize,Rini[j],iporosity);
        writer.open(outputfile.c_str());

        // Initialize loop parameters
        if (stepmethod == 0) dt = step;
        t = 0;
        tlastwrite = -1;
        sizei = sizeini;
        sizef = sizeini;
        Ri = Rini[j];
        Rf = Rini[j];
        filfacpow = 3.*cratio/(1.-cratio);
        disrupted = false;

        // Compute gas and dust quantities at t=0
        hg = Hg(Ri,q,R0,hg0);
        sigma = Sigma(Ri,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
        dustfrac = DustFrac(dustfrac0,dustfracmax,Ri,Rbump,bumpwidth,ibump);
        cg = Cg(Ri,mstar,hg);
        rhog = Rhog(sigma,hg);

        if(sizei > a0 && iporosity == 1)
        {   filfaci = FilFacSColl(Ri,mstar,rhog,cg,0.,sizeini/a0,eroll,a0,rhos,alpha,porreg,filfacpow); }
        else
        {   filfaci = 1.;
            porreg = 1;
        }
        
        filfacf = filfaci;
        massi = GrainMass(sizei,filfaci,rhos);
        massf = GrainMass(sizei,filfaci,rhos);
        st = St(Ri,mstar,rhog,cg,sizei,filfaci,rhos,0.,dragreg); //we assume deltav very small for small grain strongly coupled with gas in the epstein regime
        deltav = DeltaV(Ri,mstar,p,q,rhog,cg,R0,sigma0,hg0,dustfrac,st,alpha,ibr,ibump,idrift,Rbump,bumpwidth,bumpheight);
        vrel = Vrel(cg,st,alpha,deltav);
        //vfrag = vfragi; // need to be initialized for erosion even if ifrag=0

        // Write in output file quantities at t=0
        WriteOutputHeader(writer,massorsize);
        WriteOutputFile(writer,t,Ri,massi,filfaci,sizei,st,cg,sigma,rhog,dustfrac,vrel,Omegak(Ri,mstar),0.,0.,dragreg,porreg);

        // This is where the loop begin
        switch (massorsize)
        {
            case(0):
            {
                while(t < tend && istate[j] == 0)
                {
                    // Compute additionnal quantities for ifrag = 1 or 2, ibounce = 1 and idrift = 1
                    if (ifrag > 0)  
                    {   vfrag = Vfrag(Ri,isnow,Rsnow,filfaci,filfaclim,vfragi,vfragin,vfragout,constvfrag);    }

                    // Compute additionnal quantities for ibounce = 1
                    if (ibounce == 1)
                    {
                        vstick = Vstick(sizef,filfaci,rhos,esurf,youngmod0);
                        probabounce = ProbaBounce(filfaci,filfacbnc,vrel,vstick,Vend(vstick));
                    }

                    // Compute additionnal quantities for idrift = 1: vdrift=drdt
                    if (idrift == 1)
                    {   drdt = DRDt(Ri,mstar,p,q,rhog,cg,R0,sigma0,hg0,dustfrac,st,alpha,ibr,ibump,Rbump,bumpwidth,bumpheight); }

                    // Compute dm/dt
                    dmdt = DmDt(sizef,rhog,dustfrac,vrel,ifrag,ieros,veros,ibounce,vfrag,vstick,probabounce);

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
                    //checksize(massf,filfaci,rhos,sizeini);   // Used to check compatibility with Phantom

                    // Compute filling factor for porous grains
                    if (iporosity == 1)
                    {
                        ncoll = Ncoll(dt,Tcoll(sizef,rhog,filfaci,rhos,dustfrac,vrel));
                        // Compute the new filling factor after dt
                        filfacf = FilFacMFinal(Ri,mstar,rhog,cg,deltav,st,massf,massi,filfaci,a0,rhos,eroll,alpha,ncoll,ifrag,
                                               ibounce,vfrag,icomp,vrel,vstick,probabounce,filfaclim,Yd0,Ydpower,porreg);
                    }

                    // Compute gas and dust quantities at time t and radius R for a new loop
                    if (idrift == 1)   Ri = Rf;
                    massi = massf;
                    filfaci = filfacf;
                    sizef = GrainMassToSize(massf,filfacf,rhos);

                    // Disruption by spinning motion
                    if (idisrupt == 1)
                    {   disrupted = Disruptlim(sizef,filfacf,rhos,deltav,gammaft,weirmod,esurf,a0,disrupteq);   }

                    //State of the grain
                    State(istate[j],Rf/Rin,sizef/maxsize,disrupted);

                    // Write parameters when grain is disrupted
                    if (istate[j] == 3)
                    {   
                        WriteDisruptFile(writerdisruption,Rf,massf,filfacf,sizef,st,vrel,FreqSpin(sizef,deltav,gammaft),
                        TensileStess(sizef,filfacf,rhos,deltav,gammaft),gammaft,alpha,a0);

                        //If disruption do not stop simulation
                        if (maxsize != -1)
                        {
                            Disrupt(massf,filfacf,sizeini,rhos,Rf,mstar,rhog,cg,deltav,st,eroll,a0,alpha,porreg);
                            massi = massf;
                            filfaci=filfacf;
                            sizef = GrainMassToSize(massf,filfacf,rhos);
                            istate[j] = 0;
                        }
                    }

                    // Update hydro quantities
                    hg = Hg(Rf,q,R0,hg0);
                    sigma = Sigma(Rf,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
                    dustfrac = DustFrac(dustfrac0,dustfracmax,Rf,Rbump,bumpwidth,ibump);
                    cg = Cg(Rf,mstar,hg);
                    rhog = Rhog(sigma,hg);
                    st = St(Rf,mstar,rhog,cg,sizef,filfacf,rhos,deltav,dragreg);
                    deltav = DeltaV(Rf,mstar,p,q,rhog,cg,R0,sigma0,hg0,dustfrac,st,alpha,ibr,ibump,idrift,Rbump,bumpwidth,bumpheight);
                    vrel = Vrel(cg,st,alpha,deltav);

                    // Write in output file quantities at time t
                    if ( t-tlastwrite > 0.01 || (t == tend && istate[j] == 0 ))
                    {   
                        WriteOutputFile(writer,t,Rf,massf,filfacf,sizef,st,cg,sigma,rhog,dustfrac,vrel,Omegak(Rf,mstar),drdt,dmdt,dragreg,porreg);
                        tlastwrite = t;
                    }

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
                    if (ifrag > 0)
                    {   vfrag = Vfrag(Ri,isnow,Rsnow,filfaci,filfaclim,vfragi,vfragin,vfragout,constvfrag);    }

                    // Compute additionnal quantities for idrift = 1: vdrift=drdt
                    if (idrift == 1)
                    {   drdt = DRDt(Ri,mstar,p,q,rhog,cg,R0,sigma0,hg0,dustfrac,st,alpha,ibr,ibump,Rbump,bumpwidth,bumpheight);  }

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
                        // ncoll for grain compaction
                        ncoll = Ncoll(dt,Tcoll(sizef,rhog,filfaci,rhos,dustfrac,vrel));

                        // Compute the new filling factor after dt
                        filfacf = FilFacSFinal(Ri,mstar,rhog,cg,deltav,st,sizef,sizei,filfaci,a0,rhos,eroll,alpha,
                                         ifrag,vfrag,vrel,dragreg,porreg,filfacpow);
                    }

                    // Compute gas and dust quantities at time t and radius R for a new loop
                    if (idrift == 1)   Ri = Rf;
                    sizei = sizef;
                    filfaci = filfacf;
                    massf = GrainMass(sizef,filfacf,rhos);

                    // Disruption by spinning motion
                    if (idisrupt == 1)
                    {   disrupted = Disruptlim(sizef,filfacf,rhos,deltav,gammaft,weirmod,esurf,a0,disrupteq);   }

                    //State of the grain
                    State(istate[j],Rf/Rin,sizef/maxsize,disrupted);

                    // Write parameters when grain is disrupted
                    if (istate[j] == 3)
                    {
                        WriteDisruptFile(writerdisruption,Rf,massf,filfacf,sizef,st,vrel,FreqSpin(sizef,deltav,gammaft),
                        TensileStess(sizef,filfacf,rhos,deltav,gammaft),gammaft,alpha,a0);
                    }

                    // Update hydro quantities
                    hg = Hg(Rf,q,R0,hg0);
                    sigma = Sigma(Rf,p,R0,sigma0,ibump,Rbump,bumpwidth,bumpheight);
                    dustfrac = DustFrac(dustfrac0,dustfracmax,Rf,Rbump,bumpwidth,ibump);
                    cg = Cg(Rf,mstar,hg);
                    rhog = Rhog(sigma,hg);
                    st = St(Rf,mstar,rhog,cg,sizef,filfacf,rhos,deltav,dragreg);
                    deltav = DeltaV(Rf,mstar,p,q,rhog,cg,R0,sigma0,hg0,dustfrac,st,alpha,ibr,ibump,idrift,Rbump,bumpwidth,bumpheight);
                    vrel = Vrel(cg,st,alpha,deltav);

                    // Write in output file quantities at time t, 
                    if ( t-tlastwrite > 0.01 || (t == tend && istate[j] == 0 ))
                    {   
                        WriteOutputFile(writer,t,Rf,massf,filfacf,sizef,st,cg,sigma,rhog,dustfrac,vrel,Omegak(Rf,mstar),drdt,dsdt,dragreg,porreg);
                        tlastwrite = t;
                    }

                    // Waiting animation
                    ProgressionAnim(t,tend,istate[j],outputfile);
                }
            break;
            }
        }
        writer.close();
    }
    writerdisruption.close();

    // Compute wall time
    auto wtend = std::chrono::high_resolution_clock::now();
    double wt = (chrono::duration_cast<chrono::nanoseconds>(wtend-wtbegin)).count()* 1e-9;
    Walltime(wt);

    // Write initials conditions
    WriteInitFile(massorsize,tend,stepmethod,step,isetdens,isettemp,Rin,Rout,R0,mstar,mdisc,sigma0,hg0,T0,dustfrac0,rhog0,cg0,p,q,alpha,ibr,
                  ibump,Rbump,iporosity,sizeini,a0,rhos,idrift,ibounce,idisrupt,ifrag,vfragi,ieros,verosi,icomp,gammaft,disrupteq,isnow,
                  vfragin,vfragout,Rsnow,ngrains,Rini,istate,wt);

    EndAnim();

    }


/*------------------------ MAIN ------------------------*/

int main(int argc, char* argv[])
{   
    if (argc == 1)
    {   
        Pamdeas("input.in");    
    }
    else if (argc == 2)
    {
        if (strcmp(argv[1], "setup") == 0)
        {   
            WriteInputFile();
            cout << "Input file has been written: input.in" << endl;
        }
        else
        {   Pamdeas(string(argv[1]));   }
    }
    else
    {   cerr << "Error: to many arguments passed" << endl;  }

    return 0;
}