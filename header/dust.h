#ifndef DUST_H_INCLUDED
#define DUST_H_INCLUDED

#include <iostream>

using namespace std;


/* ------------------------ GRAIN MASS & SIZE ----------------------------- */

// Compute grain mass [kg]
double GrainMass(const double& size, const double& phi, const double& rhos);

// Compute the cross section of a grain [m²]
double GrainCrossSection(const double& size);

// Compute grain size from mass, density and filfac [m³]
double GrainMassToSize(const double& mass, const double& phi, const double& rhos);

// Compute grain volume from mass, density and filfac [m³]
double GrainVolumeMass(const double& mass, const double& phi, const double& rhos);

// Compute grain volume from size, density and filfac [m³]
double GrainVolumeSize(const double& size, const double& phi, const double& rhos);


/* ------------------------- AERODYNAMICAL PARAMETERS ------------------------- */

// Compute the stokes number, return also in which regime you are: Ep-St<1 = 1, St-St<1 = 2, Ep-St>1 = 3, St-St>1 = 4
double St(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& size,// ->
          const double& phi, const double& rhos, const double& deltav, int& iregime);


/* ------------------------ ENERGIES & VELOCITIES ------------------------ */

// Simple presciption for velocity difference between gas and dust [m/s]
double DeltaV2(const double& R, const double& Mstar, const double& p, const double& q, const double& cg, const double& st);

// Compute the relative velocity between two grains [m/s]
double Vrel(const double& cg, const double& st, const double& alpha);

// Compute the kinetic energy of 2 grains at the impact [J]
double Ekin(const double& mass, const double& vrel);
double Ekin(const double& size, const double& phi, const double& rhos, const double& vrel);

// Young Modulus (Shimaki 2010) [Pa]
double YoungMod(const double& phi, const double& youngmod0);

// Rolling Energy [J]
double Eroll(const double& a0, const double& esurf, const double& youngmod0);

// Dynamic compression resistance (Shimaki and Arakawa, 2012b) [Pa]
double Yd(const double& phi, const double& philim, const double& Yd0, const double& Ydpower);

// Velocity limit between full sticking regime and partial sticking + bouncing regime [m/s]
double Vstick(const double& size, const double& phi, const double& rhos, const double& esurf, const double& youngmod0);

// Velocity limit between elastic and inelastic bouncing regime [m/s]
double Vyield(const double& size, const double& phi, const double& rhos, const double& esurf, const double& youngmod0);
double Vyield(const double& vstick);

// Velocity limit between partial sticking + bouncing regime and full bouncing regime [m/s]
double Vend(const double& size, const double& phi, const double& rhos, const double& esurf, const double& youngmod0);
double Vend(const double& vstick);

// Fragmentation threshold velocity [m/s]
double Vfrag(const double& phi, const double& philim, const double& vfragi, const int& constvfrag);

// Restitution coefficient e
double CoeffRest(const double& vrel, const double& vstick, const double& vyield);

#endif // DUST_H_INCLUDED