#ifndef GENERAL_FUNCTIONS_H_INCLUDED
#define GENERAL_FUNCTIONS_H_INCLUDED

#include <iostream>

using namespace std;


/* ------------------------- CONVERSIONS ------------------------- */

// Convert AU in meter
double AUtoMeter(const double& R); 

// Convert meter in AU
double MeterToAU(const double& R);

// Convert Msol in kg
double MsolToKg(const double& mass);

// Convert kg in Msol
double KgToMsol(const double& mass);

// Convert time from second to year
double SecToYear(const double& time);

// Convert time from year to second
double YearToSec(const double& time);


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


/* ------------------------- DISK QUANTITIES ------------------------- */

// Gas surface density at reference radius R0, Mdisk in Msol, Sigma [kg/m²]
double Sigma0(const double& Rin, const double& Rout, const double& R0, const double& Mdisk, const double& p,
              const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight);

// Gas surface density [kg/m³]
double Sigma(const double& R, const double& p, const double& R0, const double& sigma0,
             const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight);

// Disk height [AU] 
double Hg(const double& R, const double& q, const double& R0, const double& Hg0);

// Keplerian orbital frequency [s⁻¹]
double OmegaK(const double& R, const double& Mstar); 

// Keplerian orbital velocity [m/s]
double Vk(const double& R, const double& Mstar);

// Compute the gravity field intensity at a certain radius R [m/s²]
double Gravity(const double& R, const double& Mstar); 

// Gas sound speed at R [m/s]
double Cg(const double& R, const double& Mstar, const double& Hg);
double Cg(const double& R, const double& Mstar, const double& q, const double& R0, const double& Hg0);

// Gas temperature at R [K]
double T(const double& R, const double& q, const double& R0, const double& cg);
double T(const double& R, const double& Mstar, const double& q, const double& R0, const double& Hg0);

// Gas density at R [kg/m³]
double Rhog(const double& sigma, const double& Hg);
double Rhog(const double& R, const double& p, const double& q, const double& sigma0, const double& R0,
            const double& Hg0,const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight);

// Gas pressure at R [Pa]
double Pg(const double& rhog, const double& cg);
double Pg(const double& R, const double& Mstar, const double& p, const double& q,//->
          const double& sigma0, const double& R0, const double& Hg0,const int& ibump, const double& Rbump,//-> 
          const double& bumpwidth, const double& bumpheight);

// ν_mol : Gas molecular kinematic viscosity [m²/s]
double NuMolGas(const double& rhog, const double& cg);

// ν_t : Gas turbulent kinematic viscosity [m²/s]
double NuTurbGas(const double& R, const double& Mstar, const double& alpha, const double& cg);
double NuTurbGas(const double& R, const double& Mstar, const double& q, const double& R0,//->
                 const double& Hg0, const double& alpha);

// Mean free path λ of the gas [m]
double Lambda(const double& rhog, const double& cg); 

// Compute the limit between Epstein and Stokes regime: if TransRegEpSt<1 -> Epstein, if TransRegEpSt>1 -> Stokes
double TransRegEpSt(const double& rhog, const double& cg, const double& size);


/* ------------------------- AERODYNAMICAL PARAMETERS ------------------------- */

// Compute the stokes number, return also in which regime you are: Ep-St<1 = 1, St-St<1 = 2, Ep-St>1 = 3, St-St>1 = 4
double St(const double& R, const double& Mstar, const double& rhog, const double& cg, const double& size,// ->
          const double& phi, const double& rhos, int& iregime);

// Compute dust to gas ratio (dustfrac(R)) if ibump=1
double DustFrac(const double& dustfrac0, const double& dustfracmax, const double& R, const double& Rbump,//-> 
                const double& bumpwidth, const int& ibump);

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

#endif // GENERAL_FUNCTIONS_H_INCLUDED