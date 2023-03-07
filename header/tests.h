#ifndef TESTS_H_INCLUDED
#define TESTS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../header/constantandconversion.h"
#include "../header/disc.h"
#include "../header/dust.h"
#include "../header/evol.h"

using namespace std;


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* ---------- YOU SHOULD NOT MODIFY THIS FILE DESIGN FOR TESTS ---------*/
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

// compare 2 values within tolerance
bool CompTol(const double& x, const double& y, const double tolerance = 0.01);

// Init parameters for Disc init test
void TestDiscInitParam(double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& ngrains);

// Init parameters for Growth Mass test
void TestGrowthMassParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);

// Init parameters for Growth Size test
void TestGrowthSizeParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);

// Init parameters for Drift test
void TestDriftParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);

// Init parameters for Porosity growth Mass test
void TestPorosityGrowthMassParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp, double& youngmod0,//->
                double& esurf, double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);

// Init parameters for Porosity growth Size test
void TestPorosityGrowthSizeParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp, double& youngmod0,//->
                double& esurf, double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);

// Init parameters for Porosity growth+frag+comp+bounce Mass test
void TestPorosityAllMassParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//->
                double& youngmod0, double& esurf, double& filfaclim, double& maxsize, int& isnow, int& constvfrag,//->
                int& ngrains, vector <double>& Rini, vector <int>& istate);

// Init parameters for Disruption test
void TestDisruptParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//->
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//->
                double& youngmod0, double& esurf, double& gammaft, int& disrupteq, double& maxsize,//->
                int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);

// Test the initialisation of a given disc by comparing computed values with reference values at a given radius != Ref
void TestInitDisc(const double& R, const double& sigma, const double& hg, const double& cg, const double& rhog,//->
                  const double& T, const double& P, const double& mdisc, const double& vk, const double& numol,//->
                  const double& nuturb, const double& gaspath, const int& whichtest, bool& ifailinit);

// Test the initialisation of a given disc by using all possible equations available
void TestAllDiscConfig(const double& Rin, const double& mstar, const double& mdisc, const double& p, const double& q,//->
                       const double& sigma0, const double& R0, const double& hg0, const double& alpha, const int& ismooth,//-> 
                       const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight,//->
                       double& hgtest, double& sigmatest, double& rhogtest, double& cgtest, double& tgtest, double& pgtest,//->
                       double& vktest, double& numoltest, double& nuturbtest, double& gaspathtest, const int& whichtest, bool& ifailinit);

// Test the evalutation of the growth rate with courant condition with respect to the analytic solution
void TestGrowthCompare(const double& time, const double& size, const double& st, const double& st0, const double& rhog,//->
                       const double& cg, const double& dustfrac, const double& alpha, const double& rhos, const double& omegak,//->
                       bool& ifailtest, const bool& verbosetest, ofstream& outputfile);

// Test the evalutation of the drift rate with courant condition with respect to the analytic solution
void TestDriftCompare(const double& time, const double& R, const double& p, const double& q, const double& cg, const double& vk, const double& st, const double drdt,//->
                      const double& deltav, bool& ifailtest, const bool& verbosetest, ofstream& outputfile);


// Test the porosity and size evaluation with growth and frag by comparing computed values with reference values
void TestPorosityCompare(const double& size, const double& filfac, bool& ifailtest);

// Test the porosity and size evaluation with growth, frag, bounce and comp by comparing computed values with reference values
void TestPorosityAllCompare(const double& size, const double& filfac, bool& ifailtest);

// Test the disruption evalutation by comparing computed values with reference values
void TestDisruptCompare(const double& freqspin, const double& tensilestess, bool& ifailtest);

// Write summary of the test
void WriteTestsResultsFiles(const double& tend, const double& Rtest, const double& mdisc, const double& sigmatest,//->
                    const double& hgtest, const double& tgtest, const double& pgtest, const double& rhogtest,const double& cgtest,//->
                    double& vktest, double& numoltest, double& nuturbtest, double& gaspathtest,const double& p, const double& q,//->
                    const double& alpha, const double& sizeini, const double& a0,const double& rhos, const int& idrift, const double& icomp,//->
                    const int& ifrag, const double& vfragi, const int& iporosity, const int& ibounce, const int& idisrupt, const double& gammaft,//->
                    const int& disrupteq, const double& walltime, const string& test, const bool& ifailinit, const bool& ifailtest);

// ask to write output from test or not
void VerboseTest(bool& verbosetest);

#endif // TESTS_H_INCLUDED
