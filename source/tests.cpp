#include <iostream>
#include <cmath>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"
#include "../header/evol.h"
/*

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
    double ejectasize;      // Size of ejecta [m]
    double cohacc;          // Strength of the cohesive acceleration in kg/s^2
    double a0;              // Monomer size [m]
    double rhos;            // Dust monomer density [kg/m³]
    double youngmod0;       // Young Modulus of grains [PA]
    double esurf;           // Surface energy of grains [J/m²]
    double Yd0;             // Dynamic compression resistance constant [Pa]
    double Ydpower;         // Dynamic compression resistance power law
    int    isetdens;        // Set density profile option
    int    isettemp;        // Set temperature profile option
    int    ismooth;         // Smooth inner disc
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
    int    constvfrag;      // Constant vfrag option (when fragmentation is enabled)
    double filfaclim;       // Filling factor dynamic compression resistance limit (when ifrag is enabled)
    double filfacbnc;       // Filling factor bounce limit (when ibounce is enabled)
    double maxsize;         // Limit to the max size [m]
    double gammaft;         // Force-to-torque efficiency
    int    disrupteq;       // Disruption equation (0=Tatsuuma2019, 1=Kimura2020)
    double weirmod;         // Weirbull modulus if disrupteq=1 (Kimura et al. 2020)
    int    ngrains;         // Number of grains
    vector <double> Rini;   // Initials radii
*/
    /*------------------------ INITIALS PARAMETERS TO COMPUTE ------------------------*/
/*
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
    double vstick;          // Sticking velocity [m/s]
    double probabounce;     // Probability for a grain to bounce
    double ncoll;           // number of collision per dt [s⁻¹]
    double eroll;           // Rolling energy [J]

*/
    /*------------------------ LOOP PARAMETERS ------------------------*/
/*
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


void Test_dustgrowth()*/