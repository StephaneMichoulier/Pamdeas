#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>


#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"
#include "../header/evol.h"
#include "../header/readwrite.h"

using namespace std;


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* ---------- YOU SHOULD NOT MODIFY THIS FILE DESIGN FOR TESTS ---------*/
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */


bool CompTol(const double& x, const double& y, const double tolerance = 0.01)
{
    // compare 2 values within tolerance
    if(abs(x - y)/y < tolerance)  return true; //they are same
    else                          return false;//they are not same
}

void TestDiscInitParam(double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& ngrains)
{
    tend = 0;       stepmethod = 0;       step = 0.01;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 0.5;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;
    ngrains = 0;
    mstar = MsolToKg(mstar);
}

void TestGrowthMassParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//-<
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 0;     
    tend = 1e4;     stepmethod = 2;     step = 0.01;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 0.5;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 0;
    sizeini = 1e-5; a0 = 1e-6;          rhos = 1000;
    idrift = 0;     ibounce = 0;        idisrupt = 0;       ifrag = 0;          vfragi = 0;
    ieros = 0;      icomp = 0;
    maxsize = -1;
    isnow = 0;      constvfrag = 1;
    ngrains = 1;
    Rini.resize(ngrains);
    istate.resize(ngrains);
    Rini[0] = 10.; 
    istate[0] = 0;
    mstar = MsolToKg(mstar);
}

void TestGrowthSizeParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//-<
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 1;     
    tend = 1e4;     stepmethod = 2;     step = 0.01;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 0.5;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 0;
    sizeini = 1e-5; a0 = 1e-6;          rhos = 1000;
    idrift = 0;     ibounce = 0;        idisrupt = 0;       ifrag = 0;          vfragi = 0;
    ieros = 0;      icomp = 0;
    maxsize = -1;
    isnow = 0;      constvfrag = 1;
    ngrains = 1;
    Rini.resize(ngrains);
    istate.resize(ngrains);
    Rini[0] = 10.; 
    istate[0] = 0;
    mstar = MsolToKg(mstar);
}

void TestDriftParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//-<
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 0;
    tend = 1e6;     stepmethod = 2;     step = 0.01;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 0.5;
    alpha = 1e-5;
    ibr = 0;        ibump = 0;          iporosity = 0;
    sizeini = 1;   a0 = 1e-6;          rhos = 1000;
    idrift = 1;     ibounce = 0;        idisrupt = 0;       ifrag = 0;          vfragi = 0;
    ieros = 0;      icomp = 0;
    maxsize = -1;
    isnow = 0;      constvfrag = 1;
    ngrains = 1;
    Rini.resize(ngrains);
    istate.resize(ngrains);
    Rini[0] = 100.;
    istate[0] = 0;
    mstar = MsolToKg(mstar);
}

void TestPorosityGrowthMassParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp, double& youngmod0,//->
                double& esurf, double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 0;
    tend = 1.4e5;     stepmethod = 2;     step = 0.01;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 0.5;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 1;
    sizeini = 2e-7; a0 = 2e-7;          rhos = 2700;
    idrift = 1;     ibounce = 0;        idisrupt = 0;       ifrag = 2;          vfragi = 10;
    esurf = 20e-2;  youngmod0 = 72e9;
    ieros = 0;      icomp = 0;
    maxsize = -1;
    isnow = 0;      constvfrag = 1;
    ngrains = 1;
    Rini.resize(ngrains);
    istate.resize(ngrains);
    Rini[0] = 100.;
    istate[0] = 0;
    mstar = MsolToKg(mstar);
}

void TestPorosityGrowthSizeParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp, double& youngmod0,//->
                double& esurf, double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 1;
    tend = 1.4e5;     stepmethod = 2;     step = 0.01;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 0.5;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 1;
    sizeini = 2e-7; a0 = 2e-7;          rhos = 2700;
    idrift = 1;     ibounce = 0;        idisrupt = 0;       ifrag = 2;          vfragi = 10;
    esurf = 20e-2;  youngmod0 = 72e9;
    ieros = 0;      icomp = 0;
    maxsize = -1;
    isnow = 0;      constvfrag = 1;
    ngrains = 1;
    Rini.resize(ngrains);
    istate.resize(ngrains);
    Rini[0] = 100.;
    istate[0] = 0;
    mstar = MsolToKg(mstar);
}

void TestPorosityAllMassParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//->
                double& youngmod0, double& esurf, double& filfaclim, double& maxsize, int& isnow, int& constvfrag,//->
                int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 0;
    tend = 1.4e5;     stepmethod = 2;     step = 0.01;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 0.5;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 1;
    sizeini = 2e-7; a0 = 2e-7;          rhos = 2700;
    idrift = 1;     ibounce = 1;        idisrupt = 0;       ifrag = 2;          vfragi = 10;
    esurf = 20e-2;  youngmod0 = 72e9;   filfaclim = 0.3;
    ieros = 0;      icomp = 1;
    maxsize = -1;
    isnow = 0;      constvfrag = 1;
    ngrains = 1;
    Rini.resize(ngrains);
    istate.resize(ngrains);
    Rini[0] = 100.;
    istate[0] = 0;
    mstar = MsolToKg(mstar);
}

void TestDisruptParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//->
                double& youngmod0, double& esurf, double& gammaft, int& disrupteq, double& maxsize,//->
                int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 0;
    tend = 1e5;     stepmethod = 2;     step = 0.01;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 0.5;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 1;
    sizeini = 2e-7; a0 = 2e-7;          rhos = 2700;
    idrift = 1;     ibounce = 0;        idisrupt = 1;       ifrag = 0;          vfragi = 10;
    esurf = 20e-2;  youngmod0 = 72e9;
    gammaft = 0.1;  disrupteq = 0;
    ieros = 0;      icomp = 0;
    maxsize = -1;
    isnow = 0;      constvfrag = 1;
    ngrains = 1;
    Rini.resize(ngrains);
    istate.resize(ngrains);
    Rini[0] = 100.;
    istate[0] = 0;
    mstar = MsolToKg(mstar);
}

void TestInitDisc(const double& R, const double& sigma, const double& hg, const double& cg, const double& rhog,//->
                  const double& T, const double& P, const double& mdisc, const double& vk, const double& numol,//->
                  const double& nuturb, const double& gaspath, bool& ifailinit)
{
    //At 10 au, value expected
    double hgexp = 0.88914;
    double cgexp =  837.4708;
    double sigmaexp = 142.841;
    double rhogexp = 4.28418e-10;
    double Texp = 195.701;
    double Pexp = 0.00030047397;
    double mdiscexp = 0.01;
    double vkexp = 9418.889;
    double numolexp = 5214.03681;
    double nuturbexp = 1.1139483e+11;
    double gaspathexp = 12.451866;

    // Compute with computed values
    if (CompTol(sigma,sigmaexp,1e-5) == false)
    {
        ifailinit = true;
        cout  << "Gas surface density: " << sigma << ", expected: " << sigmaexp << " kg/m²" << endl;
    }
    if (CompTol(hg,hgexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Scale height: " << hg/R << ", expected: " << hgexp/R << endl;
    }
    if (CompTol(cg,cgexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Gas sound speed: " << cg << ", expected: " << cgexp << " m/s" << endl;
    }
    if (CompTol(rhog,rhogexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Gas density: " << rhog << ", expected: " << rhogexp <<" kg/m³" << endl;
    } 
    if (CompTol(T,Texp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Gas temperature: " << T << ", expected: " << Texp << " K" << endl;
    }
    if (CompTol(P,Pexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Gas Pressure: " << P << ", expected: " << Pexp << " Pa" << endl;
    }
    if (CompTol(vk,vkexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Orbital velocity: " << vk << ", expected: " << vkexp << " m/s" << endl;
    }
    if (CompTol(numol,numolexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Gas molecular viscosity: " << numol << ", expected: " << numolexp <<" m²/s" << endl;
    }
    if (CompTol(nuturb,nuturbexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Gas turbulent viscosity: " << nuturb << ", expected: " << nuturbexp << " m²/s" << endl;
    }
    if (CompTol(gaspath,gaspathexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Mean free path λ of gas: " << gaspath << ", expected: " << gaspathexp << " m" << endl;
    }
    if (CompTol(KgToMsol(mdisc),mdiscexp,1e-5) == false)
    {
        ifailinit = true;
        cout << "Disc mass: " << KgToMsol(mdisc) << ", expected: " << mdiscexp << " Msol" << endl;
    }
}

void TestAllDiscConfig(const double& Rin, const double& mstar, const double& mdisc, const double& p, const double& q,//->
                       const double& sigma0, const double& R0, const double& hg0, const double& alpha, const int& ismooth,//-> 
                       const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight,//->
                       double& hgtest, double& sigmatest, double& rhogtest, double& cgtest, double& tgtest, double& pgtest,//->
                       double& vktest, double& numoltest, double& nuturbtest, double& gaspathtest, const int& whichtest, bool& ifailinit)
{   
    //compute sigma once (1 equ)
    sigmatest   = Sigma(10,Rin,p,R0,sigma0,ismooth,ibump,Rbump,bumpwidth,bumpheight);

    // test first set of equations to initialise a disc
    hgtest      = Hg(10,q,R0,hg0);
    cgtest      = Cg(10,mstar,hgtest);
    tgtest      = T(cgtest);
    rhogtest    = Rhog(sigmatest,hgtest);
    pgtest      = Pg(rhogtest,cgtest);
    vktest      = Vk(10,mstar);
    numoltest   = NuMolGas(rhogtest,cgtest);
    if (whichtest != 4) nuturbtest  = NuTurbGas(10,mstar,alpha,cgtest);
    else nuturbtest  = NuTurbGas(10,mstar,alpha*100,cgtest);            //Change alpha for Drift test to pass, alpha different
    gaspathtest = Lambda(rhogtest,cgtest);
    TestInitDisc(10,sigmatest,hgtest,cgtest,rhogtest,tgtest,pgtest,mdisc,vktest,numoltest,nuturbtest,gaspathtest,ifailinit);
    
    if (whichtest == 1)
    {
        if (ifailinit == true)  cout << "!!!  Disc config 1 initialisation failed, check test output !!!" << endl;
        else
        {
            // test second set of equations to initialise a disc
            hgtest      = Hg(10,mstar,cgtest);
            tgtest      = T(10,mstar,q,R0,hg0);
            cgtest      = Cg(tgtest);
            rhogtest    = Rhog(sigmatest,hgtest);
            pgtest      = Pg(rhogtest,cgtest);
            vktest      = Vk(10,mstar);
            numoltest   = NuMolGas(rhogtest,cgtest);
            nuturbtest  = NuTurbGas(10,mstar,alpha,cgtest);
            gaspathtest = Lambda(rhogtest,cgtest);
            TestInitDisc(10,sigmatest,hgtest,cgtest,rhogtest,tgtest,pgtest,mdisc,vktest,numoltest,nuturbtest,gaspathtest,ifailinit);
        }
        if (ifailinit == true) cout << "!!!  Disc config 2 initialisation failed, check test output !!!" << endl;
        else{
            // test third set of equations to initialise a disc
            hgtest      = Hg(10,q,R0,hg0);
            cgtest      = Cg(10,mstar,q,R0,hg0);
            tgtest      = T(cgtest);
            rhogtest    = Rhog(10,Rin,p,q,sigma0,R0,hg0,ismooth,ibump,Rbump,bumpwidth,bumpheight);
            pgtest      = Pg(10,Rin,mstar,p,q,sigma0,R0,hg0,ismooth,ibump,Rbump,bumpwidth,bumpheight);
            vktest      = Vk(10,mstar);
            numoltest   = NuMolGas(rhogtest,cgtest);
            nuturbtest  = NuTurbGas(10,mstar,q,R0,hg0,alpha);
            gaspathtest = Lambda(rhogtest,cgtest);
            TestInitDisc(10,sigmatest,hgtest,cgtest,rhogtest,tgtest,pgtest,mdisc,vktest,numoltest,nuturbtest,gaspathtest,ifailinit);
        }
        if (ifailinit == true) cout << "!!!  Disc config 3 initialisation failed, check test output !!!" << endl;
    }
}


void TestGrowthCompare(const double& time, const double& size, const double& st, const double& st0, const double& rhog,//->
                       const double& cg, const double& dustfrac, const double& alpha, const double& rhos, const double& omegak,//->
                       bool& ifailtest, const bool& verbosetest, ofstream& outputfile)
{
    // Compute parameters for analytic solution for a stationnary disc
    double tau = (1.+dustfrac)/(sqrt(2*M_SQRT2*alpha*rossby)*omegak*dustfrac);
    double Time = YearToSec(time)/tau + 2*sqrt(st0)*(1+st0/3.);
    double sigma = pow(8 + 9*Time*Time + 3*Time*sqrt(16+9*Time*Time),1./3.);
    double stcomp = (sigma/2. + 2./sigma -2);
    double sizecomp = rhog*cg/rhos/omegak*(sigma*sigma - 4*sigma + 4)/(2.*sigma);

    // Compute with computed values
    if (CompTol(st,stcomp,5e-2) == false)
    {
        ifailtest = true;
        cout  << "Stoke number: " << st << ", expected: " << stcomp << endl;
    }
    if (CompTol(size,sizecomp,5e-2) == false)
    {
        ifailtest = true;
        cout  << "Size: " << size << ", expected: " << sizecomp << " m" << endl;
    }
    // Write in file if asked
    if (verbosetest)    TestGrowthOutputfile(outputfile,time,st,stcomp,size,sizecomp);
}

void TestPorosityCompare(const double& size, const double& filfac, bool& ifailtest)
{
    double filfacexp = 0.000658745;
    double sizeexp = 0.516713;
           
    // Compute with computed values
    if (CompTol(filfac,filfac,1e-3) == false)
    {
        ifailtest = true;
        cout  << "Filling factor: " << filfac << ", expected: " << filfacexp << endl;
    }
    if (CompTol(size,sizeexp,1e-3) == false)
    {
        ifailtest = true;
        cout  << "Size: " << size << ", expected: " << sizeexp << " m" << endl;
    }
}

void TestPorosityAllCompare(const double& size, const double& filfac, bool& ifailtest)
{
    double filfacexp = 0.424881;
    double sizeexp = 0.00257567;
    // Compute with computed values
    if (CompTol(filfac,filfac,1e-3) == false)
    {
        ifailtest = true;
        cout  << "Filling factor: " << filfac << ", expected: " << filfacexp << endl;
    }
    if (CompTol(size,sizeexp,1e-3) == false)
    {
        ifailtest = true;
        cout  << "Size: " << size << ", expected: " << sizeexp << " m" << endl;
    }
}

void TestDisruptCompare(const double& freqspin, const double& tensilestess, bool& ifailtest)
{
    double freqspinexp = 4.60644;
    double tensilestessexp = 2.52714;
    // Compute with computed values
    if (CompTol(freqspin,freqspinexp,1e-5) == false)
    {
        ifailtest = true;
        cout  << "Spin frequency: " << freqspin << ", expected: " << freqspinexp << endl;
    }
    if (CompTol(tensilestess,tensilestessexp,1e-5) == false)
    {
        ifailtest = true;
        cout  << "Tensile Stress: " << tensilestess << ", expected: " << tensilestessexp << " m" << endl;
    }
}

void TestDriftCompare(const double& time, const double& R, const double& p, const double& q, const double& cg, const double& vk, const double& st, const double drdt,//->
                      const double& deltav, bool& ifailtest, const bool& verbosetest, ofstream& outputfile)
{
    // Compute parameters for analytic solution for a stationnary disc
    double etavk = 0.5*(p + 0.5*q + 1.5)*cg*cg/vk; //eta*vk
    double drdtcomp = -2*st*etavk/(1+st*st);
    double dvorb = st*st*etavk/(1+st*st);
    double deltavcomp = sqrt(drdtcomp*drdtcomp + dvorb*dvorb);

    // Compute with computed values
    if (CompTol(drdt,drdtcomp,1e-5) == false)
    {
        ifailtest = true;
        cout  << "vdrift: " << drdt << ", expected: " << drdtcomp << " m/s"<< endl;
    }
    if (CompTol(deltav,deltavcomp,1e-4) == false)
    {
        ifailtest = true;
        cout  << "Deltav: " << deltav << ", expected: " << deltav << " m/s" << endl;
    }
    // Write in file if asked
    if (verbosetest)    TestDriftOutputfile(outputfile,time,R,drdt,drdtcomp,deltav,deltavcomp);

}

void WriteTestsResultsFiles(const double& tend, const double& Rtest, const double& mdisc, const double& sigmatest,//->
                    const double& hgtest, const double& tgtest, const double& pgtest, const double& rhogtest,const double& cgtest,//->
                    double& vktest, double& numoltest, double& nuturbtest, double& gaspathtest,const double& p, const double& q,//->
                    const double& alpha, const double& sizeini, const double& a0,const double& rhos, const int& idrift, const double& icomp,//->
                    const int& ifrag, const double& vfragi, const int& iporosity, const int& ibounce, const int& idisrupt, const double& gammaft,//->
                    const int& disrupteq, const double& walltime, const string& test, const bool& ifailinit, const bool& ifailtest)
{
    ofstream writerdoc;
    string dash;
    string mod;

    //At 10 au
    double hgexp = 0.88914;
    double cgexp =  837.4708;
    double sigmaexp = 142.841;
    double rhogexp = 4.28418e-10;
    double Texp = 195.701;
    double Pexp = 0.00030047397;
    double mdiscexp = 0.01;
    double vkexp = 9418.889;
    double numolexp = 5214.03681;
    double nuturbexp = 1.1139483e+11;
    if (test == "TestDrift") nuturbexp *= 0.01;
    double gaspathexp = 12.451866;

    writerdoc.open(test+".txt");

	writerdoc << endl<< "     Wall Time : " << walltime << " s" << endl;
    writerdoc << endl;

    if (Rtest < 10)        dash = "";
    else if (Rtest >= 10)  dash = "-";
    else if (Rtest >= 100) dash = "--";
    else                dash = "---";
    writerdoc << setprecision(4);
    writerdoc << "------------------------------------------------------" << dash << endl;
	writerdoc << "! ---------- Initials Conditions at " << Rtest << " AU ---------- !" << endl;
    writerdoc << "------------------------------------------------------" << dash << endl;
    writerdoc << endl;
	writerdoc << "    Gas surface density: " << sigmatest << ", expected: " << sigmaexp << " kg/m²" << endl;
	writerdoc << "           Scale height: " << hgtest/Rtest << ", expected: " << hgexp/Rtest << endl;
	writerdoc << "        Gas sound speed: " << cgtest << ", expected: " << cgexp << " m/s" << endl;
	writerdoc << "            Gas density: " << rhogtest << ", expected: " << rhogexp <<" kg/m³" << endl;
	writerdoc << "        Gas temperature: " << tgtest << ", expected: " << Texp << " K" << endl;
    writerdoc << "           Gas Pressure: " << pgtest << ", expected: " << Pexp << " Pa" << endl;
    writerdoc << "       Orbital velocity: " << vktest << ", expected: " << vkexp << " m/s" << endl;
	writerdoc << "Gas molecular viscosity: " << numoltest << ", expected: " << numolexp <<" m²/s" << endl;
	writerdoc << "Gas turbulent viscosity: " << nuturbtest << ", expected: " << nuturbexp << " m²/s" << endl;
    writerdoc << "Mean free path λ of gas: " << gaspathtest << ", expected: " << gaspathexp << " m" << endl;
    writerdoc << "              Disc mass: " << KgToMsol(mdisc) << ", expected: " << mdiscexp << " Msol" << endl << endl;
    if (ifailinit == false)
    {
        writerdoc << "  --> Disc correctly initialised " << endl;
    }
    else
    {
        writerdoc << " !!! Disc initialisation failed !!! " << endl;
    }
    writerdoc << endl;
    writerdoc << endl;
    writerdoc << "------------------------------------------" << endl;
    writerdoc << "! ---------- Input Parameters ---------- !" << endl;
    writerdoc << "------------------------------------------" << endl;
    writerdoc << endl;

    writerdoc << "                p index: " << p << endl;
    writerdoc << "                q index: " << q << endl;
    writerdoc << "       Alpha turbulence: " << alpha << endl;
    writerdoc << endl;

    if (test != "TestDiscinit")
    {
        writerdoc << "         Simulated time: " << tend << " yrs" << endl << endl;

        writerdoc << "   -------- Dust properties --------" << endl;
        writerdoc << endl;
        writerdoc << "                  Grain: ";

        if (iporosity == 0)  writerdoc << "compact" << endl;
        else                 writerdoc << "porous" << endl;

        writerdoc << "     Initial grain size: " << sizeini << " m" << endl;
        if (iporosity == 1) writerdoc << "           Monomer size: " << a0 << " m" << endl;
        writerdoc << " Dust intrinsic density: " << rhos << " kg/m³" << endl;

        writerdoc << endl;
        writerdoc << "   -------- Grain options --------" << endl<< endl;
        writerdoc << "           Radial drift: ";

        if (idrift == 0)    writerdoc << "no" << endl;
        else                writerdoc << "yes" << endl;

        writerdoc << "          Fragmentation: ";
        switch (ifrag)
        {
            case (0):
            {
                writerdoc << "no" << endl;
                break;
            }
            case (1):
            {
                writerdoc << "Hard model" << endl;
                break;
            }
            case (2):
            {
                writerdoc << "Smooth model" << endl;
                break;
            }
        }

        if (ifrag !=0)
        {
            writerdoc << "Fragmentation threshold: " << vfragi << " m/s" << endl;
        }

        if (iporosity == 1)
        {   writerdoc << " Compaction during frag: ";

            if (icomp  == 0)    writerdoc << "no" << endl;
            else                writerdoc << "yes" << endl;
        }

        if (iporosity == 1)
        {
            writerdoc << "                 Bounce: ";

            if (ibounce == 0)   writerdoc << "no" << endl;
            else                writerdoc << "yes" << endl;
        }

        writerdoc << "             Disruption: ";

        if (idisrupt == 0)  writerdoc << "no" << endl;
        else                writerdoc << "yes" << endl;

        if (idisrupt == 1 && iporosity == 1)
        {
            writerdoc << "    Force-to-torque eff: " << gammaft << endl;
            if (disrupteq  == 0)    writerdoc << "Disruption equation: Tatsuuma & al. 2019" << endl;
            else                    writerdoc << "Disruption equation: Kimura & al. 2020" << endl;
        }
        writerdoc.close();
        if (ifailinit == true) exit(1);
    }
}

void VerboseTest(bool& verbosetest)
{
	char answer='n'; //To verify what is entered
    //display message
    do{
    cout<<"Verbose test (y,n)? : ";
    cin >> answer;
    if (cin.fail()){    answer = 'n';   }
    cin.clear();
    cin.ignore(100,'\n');
    }while (answer!='Y' && answer!='y' && answer!='N' && answer!='n');

    if (answer=='Y' || answer=='y') verbosetest = true;
    else    verbosetest = false;
}