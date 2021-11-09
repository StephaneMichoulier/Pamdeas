#ifndef DUST_H_INCLUDED
#define DUST_H_INCLUDED

#include <iostream>

using namespace std;


/* ------------------------ GRAIN MASS & SIZE ----------------------------- */

// Compute grain mass [kg]
double GrainMass(const double& size, const double& filfac, const double& rhos);

// Compute the cross section of a grain [m²]
double GrainCrossSection(const double& size);

// Compute grain size from mass, density and filfac [m³]
double GrainMassToSize(const double& mass, const double& filfac, const double& rhos);

// Compute grain volume from mass, density and filfac [m³]
double GrainVolumeMass(const double& mass, const double& filfac, const double& rhos);

// Compute grain volume from size, density and filfac [m³]
double GrainVolumeSize(const double& size, const double& filfac, const double& rhos);


/* ------------------------- AERODYNAMICAL PARAMETERS ------------------------- */

// Compute the stokes number, return also in which regime you are: Ep-St<1 = 1, St-St<1 = 2, Ep-St>1 = 3, St-St>1 = 4
double St(const double& R, const double& mstar, const double& rhog, const double& cg, const double& size,// ->
          const double& filfac, const double& rhos, const double& deltav, int& iregime);

// Accurate velocity difference between gas and dust [m/s]
double DeltaV(const double& R, const double& mstar, double p, double q, const double& rhog, const double& cg, const double& R0,// -> 
              const double& sigma0, const double& hg0, const double& dustfrac, const double& st, const double& alpha, const int& ibr,// ->
              const int& ibump, const int& idrift, const double& Rbump, const double& bumpwidth, const double& bumpheight);

/* ------------------------ ENERGIES & VELOCITIES ------------------------ */

// Compute the relative velocity between two grains [m/s]
double Vrel(const double& cg, const double& st, const double& alpha, const double& deltav);

// Compute the kinetic energy of 2 grains at the impact [J]
double Ekin(const double& mass, const double& vrel);
double Ekin(const double& size, const double& filfac, const double& rhos, const double& vrel);

// Young Modulus (Shimaki 2010) [Pa]
double YoungMod(const double& filfac, const double& youngmod0);

// Rolling Energy [J]
double Eroll(const double& a0, const double& esurf, const double& youngmod0);

// Dynamic compression resistance (Shimaki and Arakawa, 2012b) [Pa]
double Yd(const double& filfac, const double& filfaclim, const double& Yd0, const double& Ydpower);

// Velocity limit between full sticking regime and partial sticking + bouncing regime [m/s]
double Vstick(const double& size, const double& filfac, const double& rhos, const double& esurf, const double& youngmod0);

// Velocity limit between elastic and inelastic bouncing regime [m/s]
double Vyield(const double& size, const double& filfac, const double& rhos, const double& esurf, const double& youngmod0);
double Vyield(const double& vstick);

// Velocity limit between partial sticking + bouncing regime and full bouncing regime [m/s]
double Vend(const double& size, const double& filfac, const double& rhos, const double& esurf, const double& youngmod0);
double Vend(const double& vstick);

// Fragmentation threshold velocity [m/s]
double Vfrag(const double& R, const int& isnow, const double& Rsnow, const double& filfac, const double& filfaclim,// ->
             const double& vfragi, const double& vfragin, const double& vfragout, const int& constvfrag);

// Restitution coefficient e
double CoeffRest(const double& vrel, const double& vstick, const double& vyield);


/* ------------------------ GRAIN STATE ------------------------ */

// Determine if a grain is alive, reaches maxsize, accreted or disrupted
void State(int& istate, const double& RfRin, const double& sizelimsize, const bool& disrupted);

#endif // DUST_H_INCLUDED