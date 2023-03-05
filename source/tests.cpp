#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>


#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"
#include "../header/evol.h"

using namespace std;


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* ---------- YOU SHOULD NOT MODIFY THIS FILE DESIGN FOR TESTS ---------*/
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */


bool CompTol(const double& x, const double& y, float tolerance = 0.01)
{
   if(abs(x - y) < tolerance)
      return true; //they are same
      return false; //they are not same
}

void TestGrowthFragMassParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//-<
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 0;     
    tend = 1e3;     stepmethod = 0;     step = 1;
    profile = 0;    isetdens = 1;       isettemp = 0;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01; hg0R0 = 0.05;
    p = 1;          q = 1;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 0;
    sizeini = 1e2;  a0 = 1e-6;          rhos = 300;
    idrift = 0;     ibounce = 0;        idisrupt = 0;       ifrag = 1;          vfragi = 1000;
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

void TestGrowthFragSizeParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//-<
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    massorsize = 1;     
    tend = 1e3;     stepmethod = 0;     step = 1;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01;
    p = 1;          q = 1;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 0;
    sizeini = 1e2;  a0 = 1e-6;          rhos = 300;
    idrift = 0;     ibounce = 0;        idisrupt = 0;       ifrag = 1;          vfragi = 1000;
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
    massorsize = 1;     
    tend = 1e3;     stepmethod = 0;     step = 1;
    profile = 0;    isetdens = 1;       isettemp = 1;       ismooth = 0;
    Rin = 1;        Rout = 100;         R0 = 1;
    mstar = 1;      sigma0 = 1428.41;   T0 = 618.861;       dustfrac0 = 0.01;
    p = 1;          q = 1;
    alpha = 1e-3;
    ibr = 0;        ibump = 0;          iporosity = 0;
    sizeini = 1e-1; a0 = 1e-6;          rhos = 300;
    idrift = 1;     ibounce = 0;        idisrupt = 0;       ifrag = 1;          vfragi = 1000;
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


void TestInitDisc(const double& R, const double& sigma, const double& hg, const double& cg, const double& rhog,//->
                  const double& T, const double& P, const double& mdisc, const double& vk, const double& numol,//->
                  const double& nuturb, const double& gaspath, bool& ifailinit)
{
    //At 10 au
    double hgRexp = 0.05;
    double cgexp = 470.944;
    double sigmaexp = 142.841;
    double rhogexp = 7.61847e-10;
    double Texp = 61.8861; 
    double Pexp = 0.00016897;
    double mdiscexp = 0.01;
    double vkexp = 9418.889;
    double numolexp = 1648.8235;
    double nuturbexp = 3.5226e+10;
    double gaspathexp = 7.0022;

    if (CompTol(sigma,sigmaexp,1e-3) == false)
    {
        ifailinit = true;
        cout  << "Gas surface density: " << sigma << ", expected: " << sigmaexp << " kg/m²" << endl;
    }
    if (CompTol(hg/R,hgRexp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Scale height: " << hg/R << ", expected: " << hgRexp << endl;
    }
    if (CompTol(cg,cgexp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Gas sound speed: " << cg << ", expected: " << cgexp << " m/s" << endl;
    }
    if (CompTol(rhog,rhogexp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Gas density: " << rhog << ", expected: " << rhogexp <<" kg/m³" << endl;
    } 
    if (CompTol(T,Texp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Gas temperature: " << T << ", expected: " << Texp << " K" << endl;
    }
    if (CompTol(P,Pexp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Gas Pressure: " << P << ", expected: " << Pexp << " Pa" << endl;
    }
    if (CompTol(vk,vkexp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Orbital velocity: " << vk << ", expected: " << vkexp << " m/s" << endl;
    }
    if (CompTol(numol,numolexp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Gas molecular viscosity: " << numol << ", expected: " << numolexp <<" m²/s" << endl;
    }
    if (CompTol(nuturb/1e10,nuturbexp/1e10,1e-3) == false)
    {
        ifailinit = true;
        cout << "Gas turbulent viscosity: " << nuturb << ", expected: " << nuturbexp << " m²/s" << endl;
    }
    if (CompTol(gaspath,gaspathexp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Mean free path λ of gas: " << gaspath << ", expected: " << gaspathexp << " m" << endl;
    }
    if (CompTol(KgToMsol(mdisc),mdiscexp,1e-3) == false)
    {
        ifailinit = true;
        cout << "Disc mass: " << KgToMsol(mdisc) << ", expected: " << mdiscexp << " Msol" << endl;
    }

}

void TestAllDiscConfig(const double& Rin, const double& mstar, const double& mdisc, const double& p, const double& q,//->
                       const double& sigma0, const double& R0, const double& hg0, const double& alpha, const int& ismooth,//-> 
                       const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight,//->
                       double& hgtest, double& sigmatest, double& rhogtest, double& cgtest, double& tgtest, double& pgtest,//->
                       double& vktest, double& numoltest, double& nuturbtest, double& gaspathtest, bool& ifailinit)
{   
    sigmatest   = Sigma(10,Rin,p,R0,sigma0,ismooth,ibump,Rbump,bumpwidth,bumpheight);

    hgtest      = Hg(10,q,R0,hg0);
    cgtest      = Cg(10,mstar,hgtest);
    tgtest      = T(cgtest);
    rhogtest    = Rhog(sigmatest,hgtest);
    pgtest      = Pg(rhogtest,cgtest);
    vktest      = Vk(10,mstar);
    numoltest   = NuMolGas(rhogtest,cgtest);
    nuturbtest  = NuTurbGas(10,mstar,alpha,cgtest);
    gaspathtest = Lambda(rhogtest,cgtest);
    TestInitDisc(10,sigmatest,hgtest,cgtest,rhogtest,tgtest,pgtest,mdisc,vktest,numoltest,nuturbtest,gaspathtest,ifailinit);

    if (ifailinit == true)  cout << "!!!  Disc config 1 initialisation failed, check test output !!!" << endl;
    else
    {
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

void TestGrowthCompare()
{


}





void WriteTestsResultsFiles(const int& massorsize, const double& tend, const double& Rtest, const double& mdisc, const double& sigmatest,//->
                    const double& hgtest, const double& tgtest, const double& pgtest, const double& rhogtest,const double& cgtest,//->
                    double& vktest, double& numoltest, double& nuturbtest, double& gaspathtest,//->
                    const double& p, const double& q, const double& alpha, const double& sizeini, const double& a0,const double& rhos,//->
                    const int& idrift, const int& ifrag, const double& vfragi, const double& walltime, const string& test,//-> 
                    const bool& ifailinit, const bool& ifailtest)
{
    ofstream writerdoc;
    string dash;
    string mod;

    //At 10 au
    double hgRexp = 0.05;
    double cgexp = 470.944;
    double sigmaexp = 142.841;
    double rhogexp = 7.61847e-10;
    double Texp = 61.8861; 
    double Pexp = 0.00016897;
    double mdiscexp = 0.01;
    double vkexp = 9418.889;
    double numolexp = 1648.8235;
    double nuturbexp = 3.5226e+10;
    double gaspathexp = 7.0022;

    if (massorsize == 0)    mod = "M";
    else                    mod = "S";
    
    writerdoc.open(test+"_"+mod+".txt");

	writerdoc << endl<< "     Wall Time : " << walltime << " s" << endl;
    writerdoc << endl;

    if (Rtest < 10)        dash = "";
    else if (Rtest >= 10)  dash = "-";
    else if (Rtest >= 100) dash = "--";
    else                dash = "---";
    writerdoc << setprecision(10);
    writerdoc << "------------------------------------------------------" << dash << endl;
	writerdoc << "! ---------- Initials Conditions at " << Rtest << " AU ---------- !" << endl;
    writerdoc << "------------------------------------------------------" << dash << endl;
    writerdoc << endl;
	writerdoc << "    Gas surface density: " << sigmatest << ", expected: " << sigmaexp << " kg/m²" << endl;
	writerdoc << "           Scale height: " << hgtest/Rtest << ", expected: " << hgRexp << endl;
	writerdoc << "        Gas sound speed: " << cgtest << ", expected: " << cgexp << " m/s" << endl;
	writerdoc << "            Gas density: " << rhogtest << ", expected: " << rhogexp <<" kg/m³" << endl;
	writerdoc << "        Gas temperature: " << tgtest << ", expected: " << Texp << " K" << endl;
    writerdoc << "           Gas Pressure: " << pgtest << ", expected: " << Pexp << " K" << endl;
    writerdoc << "       Orbital velocity: " << vktest << ", expected: " << vkexp << " m/s" << endl;
	writerdoc << "Gas molecular viscosity: " << numoltest << ", expected: " << numolexp <<" m²/s" << endl;
	writerdoc << "Gas turbulent viscosity: " << nuturbtest << ", expected: " << nuturbexp << " m²/s" << endl;
    writerdoc << "Mean free path λ of gas: " << gaspathtest << ", expected: " << gaspathexp << " m" << endl;
    writerdoc << "              Disc mass: " << KgToMsol(mdisc) << ", expected: " << mdiscexp << " Msol" << endl;
    if (ifailinit == false)
    {
        writerdoc << "  Disc correctly initialised " << endl;
    }
    else
    {
        writerdoc << "  Disc initialisation failed" << endl;
    }
    writerdoc << endl;
    writerdoc << endl;
    writerdoc << "------------------------------------------" << endl;
    writerdoc << "! ---------- Input Parameters ---------- !" << endl;
    writerdoc << "------------------------------------------" << endl;
    writerdoc << endl;
    writerdoc << "             Model used: ";

    if (massorsize == 0)     writerdoc << "dm/dt" << endl;
    else                     writerdoc << "ds/dt" << endl;

    writerdoc << "         Simulated time: " << tend << " yrs" << endl;

    writerdoc << "                p index: " << p << endl;
    writerdoc << "                q index: " << q << endl;
    writerdoc << "       Alpha turbulence: " << alpha << endl;
    writerdoc << endl;

    writerdoc << "   -------- Dust properties --------" << endl;
    writerdoc << endl;
    writerdoc << "                  Grain: ";

    writerdoc << "     Initial grain size: " << sizeini << " m" << endl;
    writerdoc << " Dust intrinsic density: " << rhos << " kg/m³" << endl;

    writerdoc << endl;
    writerdoc << "   -------- Grain options --------" << endl<< endl;
    writerdoc << "           Radial drift: ";

    if (idrift == 0)    writerdoc << "no" << endl;
    else                writerdoc << "yes" << endl;

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

    writerdoc.close();
    if (ifailinit == true) exit(1);
}