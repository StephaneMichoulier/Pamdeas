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


bool CompTol(const double& x, const double& y, float tolerance = 0.01);

void TestGrowthFragMassParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//-<
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);

void TestGrowthFragSizeParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//-<
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);

void TestDriftParam(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
                double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
                double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, int& iporosity, double& sizeini, double& a0, double& rhos,//-<
                int& idrift, int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, int& icomp,//-> 
                double& maxsize, int& isnow, int& constvfrag, int& ngrains, vector <double>& Rini, vector <int>& istate);


void TestInitDisc(const double& R, const double& sigma, const double& hg, const double& cg, const double& rhog,//->
                  const double& T, const double& P, const double& mdisc, const double& vk, const double& numol,//->
                  const double& nuturb, const double& gaspath, bool& ifailinit);

void TestAllDiscConfig(const double& Rin, const double& mstar, const double& mdisc, const double& p, const double& q,//->
                       const double& sigma0, const double& R0, const double& hg0, const double& alpha, const int& ismooth,//-> 
                       const int& ibump, const double& Rbump, const double& bumpwidth, const double& bumpheight,//->
                       double& hgtest, double& sigmatest, double& rhogtest, double& cgtest, double& tgtest, double& pgtest,//->
                       double& vktest, double& numoltest, double& nuturbtest, double& gaspathtest, bool& ifailinit);

void TestGrowthCompare();

// Write summary of the test
void WriteTestsResultsFiles(const int& massorsize, const double& tend, const double& Rtest, const double& mdisc, const double& sigmatest,//->
                    const double& hgtest, const double& tgtest, const double& pgtest, const double& rhogtest,const double& cgtest,//->
                    double& vktest, double& numoltest, double& nuturbtest, double& gaspathtest,//->
                    const double& p, const double& q, const double& alpha, const double& sizeini, const double& a0,const double& rhos,//->
                    const int& idrift, const int& ifrag, const double& vfragi, const double& walltime, const string& test,//-> 
                    const bool& ifailinit, const bool& ifailtest);

#endif // TESTS_H_INCLUDED
