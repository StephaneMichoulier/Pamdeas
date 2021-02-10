#ifndef DISC_H_INCLUDED
#define DISC_H_INCLUDED

#include <iostream>

using namespace std;


/* ------------------------- DISK QUANTITIES ------------------------- */

// Gas surface density at reference radius R0, Mdisk in Msol, Sigma [kg/m²]
double Sigma0(const double& Rin, const double& Rout, const double& R0, const double& mdisk, const double& p,
              const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight);

// Gas surface density [kg/m³]
double Sigma(const double& R, const double& p, const double& R0, const double& sigma0,
             const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight);

// Disk height [AU] 
double Hg(const double& R, const double& q, const double& R0, const double& hg0);

// Keplerian orbital frequency [s⁻¹]
double Omegak(const double& R, const double& Mstar); 

// Keplerian orbital velocity [m/s]
double Vk(const double& R, const double& Mstar);

// Compute the gravity field intensity at a certain radius R [m/s²]
double Gravity(const double& R, const double& Mstar); 

// Gas sound speed at R [m/s]
double Cg(const double& R, const double& Mstar, const double& hg);
double Cg(const double& R, const double& Mstar, const double& q, const double& R0, const double& hg0);

// Gas temperature at R [K]
double T(const double& R, const double& q, const double& R0, const double& cg);
double T(const double& R, const double& Mstar, const double& q, const double& R0, const double& hg0);

// Gas density at R [kg/m³]
double Rhog(const double& sigma, const double& hg);
double Rhog(const double& R, const double& p, const double& q, const double& sigma0, const double& R0,
            const double& hg0,const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight);

// Gas pressure at R [Pa]
double Pg(const double& rhog, const double& cg);
double Pg(const double& R, const double& Mstar, const double& p, const double& q,//->
          const double& sigma0, const double& R0, const double& hg0,const int& ibump, const double& Rbump,//-> 
          const double& bumpwidth, const double& bumpheight);

// Compute dust to gas ratio (dustfrac(R)) if ibump=1
double DustFrac(const double& dustfrac0, const double& dustfracmax, const double& R, const double& Rbump,//-> 
                const double& bumpwidth, const int& ibump);

// ν_mol : Gas molecular kinematic viscosity [m²/s]
double NuMolGas(const double& rhog, const double& cg);

// ν_t : Gas turbulent kinematic viscosity [m²/s]
double NuTurbGas(const double& R, const double& Mstar, const double& alpha, const double& cg);
double NuTurbGas(const double& R, const double& Mstar, const double& q, const double& R0,//->
                 const double& hg0, const double& alpha);

// Mean free path λ of the gas [m]
double Lambda(const double& rhog, const double& cg); 

// Compute the limit between Epstein and Stokes regime: if TransRegEpSt<1 -> Epstein, if TransRegEpSt>1 -> Stokes
double TransRegEpSt(const double& rhog, const double& cg, const double& size);


#endif // DISC_H_INCLUDED