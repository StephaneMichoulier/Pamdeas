#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "../header/constantandconversion.h"
#include "../header/readwrite.h"
#include "../header/disc.h"

using namespace std;


/* ------------------------ READING ------------------------*/

void ReadFile(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp,// ->
              double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& h0R0, double& T0,// ->
              double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, double& Rbump, double& dustfracmax,// ->
              double& bumpwidth, double& bumpheight, int& iporosity, double& sizeini, double& a0, double& rhos, int& idrift,// ->
              int& ibounce, int& idisrupt, int& ifrag, double& vfragi, double& gammaft, double& limsize, double& youngmod0,// ->
              double& esurf, double& Yd0, double& Ydpower, int& constvfrag, double& filfaclim, double& filfacbnc,// ->
              int& ngrains, vector <double>& Rini, vector <int>& istate)
{
    ifstream Reader("input.in");
    if (!Reader)
   	{   
        cerr << "Error: input.in is missing" << endl;
        exit(1);
    }

    ReadVoid(Reader, 6);     CheckType(Reader, massorsize, "massorsize");
    ReadVoid(Reader, 10);    CheckType(Reader, tend, "tend");
    ReadVoid(Reader, 5);     CheckType(Reader, stepmethod, "stepmethod");
    ReadVoid(Reader, 9);     CheckType(Reader, step, "step");
    ReadVoid(Reader, 18);    CheckType(Reader, profile, "profile");

    ReadVoid(Reader, 7);     CheckType(Reader, isetdens, "isetdens");
    ReadVoid(Reader, 8);     CheckType(Reader, isettemp, "isettemp");

    ReadVoid(Reader,11);     CheckType(Reader, Rin, "Rin");
    ReadVoid(Reader, 5);     CheckType(Reader, Rout, "Rout");
    ReadVoid(Reader, 5);     CheckType(Reader, R0, "R0");
    ReadVoid(Reader, 5);     CheckType(Reader, mstar, "mstar");
    ReadVoid(Reader, 5);     CheckType(Reader, mdisc, "mdisc");
    ReadVoid(Reader, 5);     CheckType(Reader, sigma0, "sigma0");
    ReadVoid(Reader, 7);     CheckType(Reader, h0R0, "h0R0");
    ReadVoid(Reader, 5);     CheckType(Reader, T0, "T0");
    ReadVoid(Reader, 5);     CheckType(Reader, dustfrac0, "dustfrac0");
    ReadVoid(Reader, 9);     CheckType(Reader, p, "p");
    ReadVoid(Reader, 3);     CheckType(Reader, q, "q");
    ReadVoid(Reader, 2);     CheckType(Reader, alpha, "alpha");
    ReadVoid(Reader, 6);     CheckType(Reader, ibr, "ibr");
    ReadVoid(Reader, 5);     CheckType(Reader, ibump, "ibump");

    ReadVoid(Reader, 14);    CheckType(Reader, Rbump, "Rbump");
    ReadVoid(Reader, 6);     CheckType(Reader, dustfracmax, "dustfracmax");
    ReadVoid(Reader, 5);     CheckType(Reader, bumpwidth, "bumpwidth");
    ReadVoid(Reader, 9);     CheckType(Reader, bumpheight, "bumpheight");

    ReadVoid(Reader, 10);    CheckType(Reader, iporosity, "iporosity");
    ReadVoid(Reader, 4);     CheckType(Reader, sizeini, "sizeini");
    ReadVoid(Reader, 5);     CheckType(Reader, a0, "a0");
    ReadVoid(Reader, 5);     CheckType(Reader, rhos, "rhos");

    ReadVoid(Reader, 8);     CheckType(Reader, idrift, "idrift");
    ReadVoid(Reader, 5);     CheckType(Reader, ibounce, "ibounce");
    ReadVoid(Reader, 5);     CheckType(Reader, idisrupt, "idisrupt");
    ReadVoid(Reader, 6);     CheckType(Reader, ifrag, "ifrag");
    ReadVoid(Reader, 6);     CheckType(Reader, vfragi, "vfragi");
    ReadVoid(Reader, 5);     CheckType(Reader, gammaft, "gammaft");
    ReadVoid(Reader, 9);     CheckType(Reader, limsize, "limsize");

    ReadVoid(Reader, 15);     CheckType(Reader, youngmod0, "youngmod0");
    ReadVoid(Reader, 11);    CheckType(Reader, esurf, "esurf");
    ReadVoid(Reader, 13);    CheckType(Reader, Yd0, "Yd0");
    ReadVoid(Reader, 9);     CheckType(Reader, Ydpower, "Ydpower");
    ReadVoid(Reader, 10);    CheckType(Reader, constvfrag, "constvfrag");
    ReadVoid(Reader, 7);     CheckType(Reader, filfaclim, "filfaclim");
    ReadVoid(Reader, 8);     CheckType(Reader, filfacbnc, "filfacbnc");

    ReadVoid(Reader, 13);    CheckType(Reader, ngrains, "ngrain");

    ReadVoid(Reader, 7);

    if (ngrains > 0)
    {
        Rini.resize(ngrains);

        for (int i = 0; i < ngrains; i++)
        {
            ReadVoid(Reader, 2); CheckType(Reader, Rini[i],"Rini for grain "+to_string(i+1));
        }
    }

    CheckData(massorsize,tend,stepmethod,step,profile,isetdens,isettemp,Rin,Rout,R0,mstar,mdisc,sigma0,h0R0,T0,dustfrac0,p,q,alpha,ibr,ibump,
              Rbump,dustfracmax,bumpwidth,bumpheight,iporosity,sizeini,a0,rhos,idrift,ibounce,idisrupt,ifrag,vfragi,gammaft,limsize,
              youngmod0,esurf,Yd0,Ydpower,constvfrag,filfaclim,filfacbnc,ngrains,Rini,istate);
    
    // Initialize some parameters depending on option
    istate.resize(ngrains);
    for (int i = 0; i < ngrains; i++)   istate[i] = 0;

    if (stepmethod == 2)    step = 1.;

    mstar = MsolToKg(mstar);
    mdisc = MsolToKg(mdisc);

	Reader.close();
}

void CheckData(const int& massorsize, const double& tend, const int& stepmethod, const double& step, const int& profile, const int& isetdens,// ->
               const int& isettemp, const double& Rin, const double& Rout, const double& R0, const double& mstar, const double& mdisc,// ->
               const double& sigma0, const double& h0R0, const double& T0, const double& dustfrac0, const double& p, const double& q,// ->
               const double& alpha, const int& ibr, const int& ibump, const double& Rbump, const double& dustfracmax, const double& bumpwidth,// ->
               const double& bumpheight, const int& iporosity, const double& sizeini, const double& a0, const double& rhos, const int& idrift,// ->
               const int& ibounce, const int& idisrupt, const int& ifrag, const double& vfragi, const double& gammaft, const double& limsize,// ->
               const double& youngmod0, const double& esurf, const double& Yd0, const double& Ydpower, const int& constvfrag,// ->
               const double& filfaclim, const double& filfacbnc, const int& ngrains, const vector <double>& Rini, const vector <int>& istate)
{
    bool error = false;

    // Wrong values -> error
    if (massorsize != 0 && massorsize != 1)     ErrorValue(error, "Error", "massorsize != 0 & != 1");
    if (tend <= 0)                              ErrorValue(error, "Error", "tend <= 0");
    if (stepmethod != 0 && stepmethod != 1 && stepmethod != 2)
                                                ErrorValue(error, "Error", "stepmethod != 0 & != 1 & != 2");
    if (stepmethod != 2)
    {
        if (step <= 0)                          ErrorValue(error, "Error", "step <= 0");
    }

    if (profile != 0 && profile != 1)           ErrorValue(error, "Error", "profile != 0 & != 1");
    if (isetdens!= 0 && isetdens != 1)          ErrorValue(error, "Error", "isetdens!= 0 & != 1");
    if (isettemp!= 0 && isettemp != 1)          ErrorValue(error, "Error", "isettemp!= 0 & != 1");

    if (Rin <= 0)                               ErrorValue(error, "Error", "Rin <= 0");
    if (Rout <= 0)                              ErrorValue(error, "Error", "Rout <= 0");
    if (Rin >= Rout)                            ErrorValue(error, "Error", "Rin >= Rout");
    if (R0 <= 0)                                ErrorValue(error, "Error", "R0 <= 0");
    if (mstar <= 0)                             ErrorValue(error, "Error", "mstar <= 0");

    if (isetdens == 0)
    {   if (mdisc <= 0)                         ErrorValue(error, "Error", "mdisc <= 0");    }
    else
    {   if (sigma0 <= 0)                        ErrorValue(error, "Error", "sigma0 <= 0");    }

    if (isettemp == 0)
    {   if (h0R0 <= 0)                          ErrorValue(error, "Error", "H0R0 <= 0");    }
    else
    {   if (T0 <= 0)                            ErrorValue(error, "Error", "T0 <= 0");    }

    if (dustfrac0 <= 0 || dustfrac0 >= 1)
    {
        if (dustfrac0 >= 1 )                    ErrorValue(error, "Error", "dustfrac0 >= 1");
        else                                    ErrorValue(error, "Error", "dustfrac0 <= 0");
    }

    if (alpha <= 0)                             ErrorValue(error, "Error", "alpha <= 0");
    if (ibr != 0 && ibr != 1)                   ErrorValue(error, "Error", "ibr != 0 & != 1");
    if (ibump != 0 && ibump != 1)               ErrorValue(error, "Error", "ibump != 0 & != 1");

    if (Rbump < Rin || Rbump > Rout)
    {
        if (Rbump < Rin )                       ErrorValue(error, "Error", "Rbump < Rin");
        else                                    ErrorValue(error, "Error", "Rbump > Rout");
    }

    if (dustfracmax < dustfrac0 || dustfracmax >= 1)
    {
        if (dustfracmax >= 1 )                  ErrorValue(error, "Error", "dustfracmax >= 1");
        else                                    ErrorValue(error, "Error", "dustfracmax < dustfrac0");
    }

    if (bumpwidth <= 0)                         ErrorValue(error, "Error", "bumpwidth <= 0");
    if (bumpheight <= 0)                        ErrorValue(error, "Error", "bumpheight <= 0");
    if (iporosity != 0 && iporosity != 1)       ErrorValue(error, "Error", "iporosity != 0 & != 1");
    if (sizeini <= 0)                           ErrorValue(error, "Error", "sizeini <= 0");

    if (a0 <= 0 || a0 > sizeini)
    {
        if(a0 <= 0)                             ErrorValue(error, "Error", "a0 <= 0");
        else                                    ErrorValue(error, "Error", "a0 > sizeini");
    }

    if (rhos <= 0)                              ErrorValue(error, "Error", "rhos <= 0");

    if (idrift != 0 && idrift != 1)             ErrorValue(error, "Error", "idrift != 0 & != 1");
    if (ibounce != 0 && ibounce != 1)           ErrorValue(error, "Error", "ibounce != 0 & != 1");
    if (idisrupt != 0 && idisrupt != 1)         ErrorValue(error, "Error", "idisrupt != 0 & != 1");
    if (ifrag != 0 && ifrag != 1 && ifrag != 2) ErrorValue(error, "Error", "ifrag != 0 & != 1 & != 2");
    if (vfragi < 0)                             ErrorValue(error, "Error", "vfragi < 0");

    if (gammaft <= 0 || gammaft > 1)
    {
        if(gammaft > 1)                         ErrorValue(error, "Error", "gammaft > 1");
        else                                    ErrorValue(error, "Error", "gammaft <= 0");
    }

    if (limsize <= sizeini)                     ErrorValue(error, "Error", "limsize <= sizeini");

    if (youngmod0 <= 0)                         ErrorValue(error, "Error", "Youngmod0 <= 0");
    if (esurf <= 0)                             ErrorValue(error, "Error", "Esurf <= 0");
    if (Yd0 <= 0)                               ErrorValue(error, "Error", "Yd0 <= 0");
    if (constvfrag != 0 && constvfrag != 1)     ErrorValue(error, "Error", "constvfrag != 0 & != 1");

    if (filfaclim < 0 || filfaclim > 1)
    {
        if(filfaclim > 1)                       ErrorValue(error, "Error", "filfaclim > 1");
        else                                    ErrorValue(error, "Error", "filfaclim < 0");
    }

    if (filfacbnc < 0 || filfacbnc > 1)
    {
        if(filfacbnc > 1)                       ErrorValue(error, "Error", "filfacbnc > 1");
        else                                    ErrorValue(error, "Error", "filfacbnc < 0");
    }

    if (ngrains < 0)                            ErrorValue(error, "Error", "ngrain < 0");
    else
    {
        if (error == false && ngrains > 0)
        {   for (int i=0; i<ngrains; i++)
            {   
                if (Rini[i] < Rin)              ErrorValue(error, "Error", "Rini < Rin for grain " + to_string(i+1));
                if (Rini[i] > Rout)             ErrorValue(error, "Error", "Rini > Rout for grain " + to_string(i+1));
            }
        }
    }

    // Incompatible option

    if (ibounce == 1 && massorsize == 1)        ErrorValue(error, "Incompatible", "bounce is not included with ds/dt model");


    // Stop programm if error

    if (error == true)                          exit(1);

}


/* ------------------------ WRITING ------------------------*/


void WriteInputFile()
{
    ofstream writerinput;

	writerinput.open("input.in");

    writerinput << endl;
    writerinput << "#-Use dm/dt or ds/dt" << endl;
    writerinput << "massorsize = 0          >(dm/dt = 0, ds/dt = 1)" << endl;
    writerinput << endl;
    writerinput << "#-Time options" << endl;
    writerinput << "      tend = 1.0e6      >End time (yrs)" << endl;
    writerinput << "stepmethod = 2          >(Fixe=0, Fraction of orbital period=1, Adaptative dt=2)" << endl;
    writerinput << "      step = 1          >(stepmethod=0 -> step in yrs, stepmethod=1 -> step in fraction of orbital time)" << endl;
    writerinput << endl;
    writerinput << "#-Compute Disc profiles" << endl;
    writerinput << "  profile  = 0          >(0=no, 1=yes)" << endl;
    writerinput << endl;
    writerinput << "#-Set disc profiles" << endl;
    writerinput << "  isetdens = 0          >Set Density profile with (0=mdisc, 1=sigma0)" << endl;
    writerinput << "  isettemp = 0          >Set temperature profile with (0=H0/R0, 1=T0)" << endl;
    writerinput << endl;
    writerinput << "#-Gas disc properties" << endl;
    writerinput << endl;
    writerinput << "      Rin  = 3.         >Inner radius (AU)" << endl;
    writerinput << "      Rout = 300.       >Outer radius (AU)" << endl;
    writerinput << "      Rref = 100.       >Reference radius (AU)" << endl;
    writerinput << "     mstar = 1.         >Star mass (Msol)" << endl;
    writerinput << "     mdisc = 0.02       >Disc mass (Msol)" << endl;
    writerinput << "    sigma0 = 9.1        >Surface density at Rref (kg/m²)" << endl;
    writerinput << "     H0/R0 = 0.05       >H/R at Rref" << endl;
    writerinput << "        T0 = 6.2        >Temperature at Rref" << endl;
    writerinput << " dustfrac0 = 0.01       >Dust to gas ratio at Rref" << endl;
    writerinput << "   p index = 1.5        " << endl;
    writerinput << "   q index = 0.75       " << endl;
    writerinput << "     alpha = 1.e-2      >Turbulence Shakura & Sunyaev" << endl;
    writerinput << "       ibr = 0          >Back-reaction (0=no, 1=yes)" << endl;
    writerinput << "     ibump = 0          >Pressure bump (0=no, 1=yes)" << endl;
    writerinput << endl;
    writerinput << "#-Pressure bump option, Available if ibump = 1" << endl;
    writerinput << endl;
    writerinput << "      Rbump = 6.5       >Pressure bump radius (AU)" << endl;
    writerinput << "dustfracmax = 0.1       >Max dustfrac possible" << endl;
    writerinput << "  bumpwidth = 1.1       >Bump half width at half maximum (AU)" << endl;
    writerinput << " bumpheight = 300       >Bump height (Surface density at Rref)" << endl;
    writerinput << endl;
    writerinput << "#-Dust properties" << endl;
    writerinput << endl;
    writerinput << " iporosity = 1          >(compact=0, porous=1)" << endl;
    writerinput << "   sizeini = 1.0e-7     >Initial size (m)" << endl;
    writerinput << "        a0 = 1.0e-7     >Monomer size (m)" << endl;
    writerinput << "      rhos = 917        >Dust intrinsic density (kg/m³)" << endl;
    writerinput << endl;
    writerinput << "#-Grain options" << endl;
    writerinput << endl;
    writerinput << "    idrift = 0          >Drift (0=no, 1=yes)" << endl;
    writerinput << "   ibounce = 0          >Bounce (0=no, 1=yes)" << endl;
    writerinput << "  idisrupt = 0          >Rotationnal disruption (0=no, 1=yes)" << endl;
    writerinput << "     ifrag = 0          >Fragmentation (0=no, 1=hardfrag, 2=smoothfrag)" << endl;
    writerinput << "     vfrag = 15         >Fragmentation threshold (m/s)" << endl;
    writerinput << "   gammaft = 0.1        >Force-to-torque efficiency (Disruption) [Tatsuuma et al. 2021]" << endl;
    writerinput << "   maxsize = 1.0e3      >Maximum size to stop the simulation" << endl;
    writerinput << endl;
    writerinput << "#-Porosity properties, Available if iporosity = 1" << endl;
    writerinput << endl;
    writerinput << " Youngmod0 = 9.4e9      >Young Modulus for ice [Yamamoto et al. 2014] (Pa)" << endl;
    writerinput << "     Esurf = 7.3e-2     >Surface energy for ice grains J/m² [Yamamoto et al. 2014] (J/m²)" << endl;
    writerinput << "       Yd0 = 9.8e6      >Dynamic compression resistance constant for ice (Pa)" << endl;
    writerinput << "   Ydpower = 4          >Dynamic compression resistance power for ice [Mellor, 1975]" << endl;
    writerinput << "constvfrag = 1          >Constant fragmentation threshold (0=no, 1=yes)" << endl;
    writerinput << " filfaclim = 0.01       >Filling factor dynamic compression resistance limit" << endl;
    writerinput << " filfacbnc = 0.3        >Filling factor bounce limit" << endl;
    writerinput << endl;
    writerinput << "#-Number of dust grains and initial radii" << endl;
    writerinput << endl;
    writerinput << "   ngrains = 10         >Number of dust grains" << endl;
    writerinput << endl;
    writerinput << "   Initial radii (AU)" << endl;
    writerinput << "        R01 = 5" << endl;
    writerinput << "        R02 = 10" << endl;
    writerinput << "        R03 = 20" << endl;
    writerinput << "        R04 = 50" << endl;
    writerinput << "        R05 = 75" << endl;
    writerinput << "        R06 = 100" << endl;
    writerinput << "        R07 = 150" << endl;
    writerinput << "        R08 = 200" << endl;
    writerinput << "        R09 = 250" << endl;
    writerinput << "        R10 = 300" << endl;
    writerinput << endl;

    writerinput.close();
}

void WriteProfileFile(ofstream& outputprofile, const double& Rprofile, const double& hg, const double& cg, const double& sigma,
                      const double& rhog, const double& dustfrac, const double& pg, const double& T)
{
    WriteValue(outputprofile,  6, 6, Rprofile);
    WriteValue(outputprofile, 10, 6, hg);
    WriteValue(outputprofile, 10, 6, cg);
    WriteValue(outputprofile, 10, 6, sigma);
    WriteValue(outputprofile, 12, 6, rhog);
    WriteValue(outputprofile, 10, 6, dustfrac);
    WriteValue(outputprofile, 12, 6, pg);
    WriteValue(outputprofile, 10, 6, T);
    outputprofile << "\n";
}

void WriteProfileHeader()
{
    ofstream writecol("disc_profiles_header.txt");
    writecol << "R(AU)\n" << "Hg(AU)\n" << "cg(m/s)\n" << "sigma(kg/m^2)\n" << "rhog(kg/m^3)\n" << "dustfrac\n" << "Pg(Pa)\n" << "T(K)\n";
    writecol.close();
}

void WriteOutputFile(ofstream& outputfile, const double& t, const double& Rf, const double& massf, const double& filfacf,// ->
                     const double& sizef, const double& st, const double& cg, const double& sigma,// ->
                     const double& rhog, const double& dustfrac, const double& vrel, const double& omegak,// ->
                     const double& drdt, const double& dvardt, const int& iregime)
{
    WriteValue(outputfile, 11, 6, t);
    WriteValue(outputfile, 10, 6, Rf);
    WriteValue(outputfile, 12, 6, massf);
    WriteValue(outputfile, 12, 6, filfacf);
    WriteValue(outputfile, 12, 6, sizef);
    WriteValue(outputfile, 12, 6, st);
    WriteValue(outputfile, 10, 6, cg);
    WriteValue(outputfile, 10, 6, sigma);
    WriteValue(outputfile, 12, 6, rhog);
    WriteValue(outputfile, 10, 6, dustfrac);
    WriteValue(outputfile, 12, 6, vrel);
    WriteValue(outputfile, 12, 6, omegak);
    WriteValue(outputfile, 13, 6, drdt);
    WriteValue(outputfile, 12, 6, dvardt);
    WriteValue(outputfile,  2, 1, iregime);
    outputfile << "\n";
}

void WriteOutputHeader(const double& massorsize)
{
    ofstream writecol("output_header.txt");
    writecol << "t(yr)\n" << "R(AU)\n" << "mass(kg)\n" << "filfac\n" << "size(m)\n" << "St\n" << "cg(m/s)\n" << "sigma(kg/m^2)\n"
             << "rhog(kg/m^3)\n" << "dustfrac\n" << "vrel(m/s)\n" << "omegak(1/s)\n" << "drdt(au/s)\n";

    if (massorsize == 0)
    {   writecol << "dmdt(kg/s)\n";  }
    else
    {   writecol << "dsdt(m/s)\n";  }

    writecol << "drag_regime\n";
    writecol.close();
}

void WriteDisruptFile(ofstream& outputfile, const double& R, const double& massf, const double& filfacf, const double& sizef, 
                      const double& st, const double& vrel, const double& freqspin, const double& tensilestress,
                      const double& gammaft, const double& alpha, const double& a0)
{
    WriteValue(outputfile,  6, 4, R);
    WriteValue(outputfile,  6, 4, gammaft);
    WriteValue(outputfile,  6, 4, a0);
    WriteValue(outputfile,  6, 4, alpha);
    WriteValue(outputfile, 10, 6, vrel);
    WriteValue(outputfile, 12, 6, massf);
    WriteValue(outputfile, 12, 6, filfacf);
    WriteValue(outputfile, 12, 6, sizef);
    WriteValue(outputfile, 12, 6, st);
    WriteValue(outputfile, 10, 6, freqspin);
    WriteValue(outputfile, 10, 6, tensilestress);
    outputfile << "\n";
}

void WriteDisruptHeader()
{
    ofstream writecol("disrupt_param_header.txt");
    writecol << "R(AU)\n" << "gammaft\n" << "a0(m)\n" << "alpha\n" << "vrel(m/s)\n" << "mass(kg)\n" << "filfac\n" << "size(m)\n"
             << "St\n" << "freqspin(rad/s)\n" << "tensilestress(Pa)\n";

    writecol.close();
}

void WriteInitFile(const int& massorsize, const double& tend, const int& stepmethod, const double& dt, const int& isetdens,// ->
                   const int& isettemp, const double& Rin, const double& Rout, const double& R0, const double& mstar, const double& mdisc,// ->
                   const double& sigma0, const double& h0, const double& T0, const double& dustfrac0, const double& rhog0, const double& cg0,// ->
                   const double& p, const double& q, const double& alpha, const int& ibr, const int& ibump, const double& Rbump,// ->
                   const int& iporosity, const double& sizeini, const double& a0, const double& rhos, const int& idrift, const int& ibounce,// ->
                   const int& idisrupt, const int& ifrag, const double& vfragi, const double& gammaft, const int& ngrains,// ->
                   const vector <double>& Rini, const  vector <int>& istate, const double& runningtime)
{
    ofstream writerdoc;
    string dash;
    string mod;

    if (massorsize == 0)    mod = "M";
    else                    mod = "S";
    
    writerdoc.open("results_setup_"+mod+".txt");

	writerdoc << endl<< "     Running Time : " << runningtime << " s" << endl;
    writerdoc << endl;

    if (R0 < 10)        dash = "";
    else if (R0 >= 10)  dash = "-";
    else if (R0 >= 100) dash = "--";
    else                dash = "---";

    writerdoc << "------------------------------------------------------" << dash << endl;
	writerdoc << "! ---------- Initials Conditions at " << R0 << " AU ---------- !" << endl;
    writerdoc << "------------------------------------------------------" << dash << endl;
    writerdoc << endl;
	writerdoc << "    Gas surface density: " << sigma0 << " kg/m²" << endl;
	writerdoc << "           Scale height: " << h0 << " AU" << endl;
	writerdoc << "        Gas sound speed: " << cg0 << " m/s" << endl;
	writerdoc << "            Gas density: " << rhog0 << " kg/m³" << endl;
	writerdoc << "        Gas temperature: " << T0 << " K" << endl;
	writerdoc << "           Gas pressure: " << Pg(rhog0,cg0) << " Pa" << endl;
	writerdoc << "     Gas mean free path: " << Lambda(rhog0,cg0) << " m" << endl;
    writerdoc << endl;
    writerdoc << endl;
    writerdoc << "------------------------------------------" << endl;
    writerdoc << "! ---------- Input Parameters ---------- !" << endl;
    writerdoc << "------------------------------------------" << endl;
    writerdoc << endl;
    writerdoc << "             Model used: ";

    if (massorsize == 0)     writerdoc << "dm/dt" << endl;
    else                     writerdoc << "ds/dt" << endl;

    writerdoc << "         Dust particles: " << ngrains << endl;
    writerdoc << endl;
    writerdoc << "         Simulated time: " << tend << " yrs" << endl;

    writerdoc << "   Time-stepping method: ";
    switch (stepmethod)
    {
        case (0):
        {
            writerdoc << "Fixe" << endl;
            writerdoc << "                     dt: " << dt << " yrs" << endl;
            break;
        }
        case (1):
        {
            writerdoc << "Fraction of orbital period" << endl;
            writerdoc << "                     dt: " << dt << endl;
            break;
        }
        case (2):
        {
            writerdoc << "Adaptative" << endl;
            break;
        }
    }
    writerdoc << endl;
    writerdoc << "   -------- Gas disc properties --------" << endl;
    writerdoc << endl;
    writerdoc << "           Inner radius: " << Rin << " AU" << endl;
    writerdoc << "           Outer radius: " << Rout << " AU" << endl;
    writerdoc << "       Reference radius: " << R0 << " AU" << endl;
    writerdoc << "              Star mass: " << KgToMsol(mstar) << " Msol" << endl;
    writerdoc << "              Disc mass: " << KgToMsol(mdisc) << " Msol" << endl;
    if (isetdens == 1)
    {   writerdoc << "  Surface density at R0: " << sigma0 << " kg/m²" << endl; }
    if (isettemp == 0)
    {   writerdoc << "H/R at reference radius: " << h0/R0 << endl;  }
    else
    {   writerdoc << "      Temperature at R0: " << T0 << " K" << endl; }
    writerdoc << "  Dust/gas ratio (Rref): " << dustfrac0 << endl;
    writerdoc << "                p index: " << p << endl;
    writerdoc << "                q index: " << q << endl;
    writerdoc << "       Alpha turbulence: " << alpha << endl;
    writerdoc << endl;

    writerdoc << "          Back-reaction: ";

    if (ibr == 0)       writerdoc << "no" << endl;
    else                writerdoc << "yes" << endl;

    writerdoc << "          Pressure bump: ";

    if (ibump == 0)     writerdoc << "no" << endl;
    else
    {
        writerdoc << "yes" << endl;
        writerdoc << "   Pressure bump radius: " << Rbump << " AU" << endl;
    }

    writerdoc << endl;
    writerdoc << "   -------- Dust properties --------" << endl;
    writerdoc << endl;
    writerdoc << "                  Grain: ";

    if (iporosity == 0)  writerdoc << "compact" << endl;
    else                 writerdoc << "porous" << endl;

    writerdoc << "     Initial grain size: " << sizeini << " m" << endl;
    writerdoc << "           Monomer size: " << a0 << " m" << endl;
    writerdoc << " Dust intrinsic density: " << rhos << " kg/m³" << endl;
    if (ifrag !=0)      writerdoc << "Fragmentation threshold: " << vfragi << " m/s" << endl;
    if (idisrupt == 1)  writerdoc << "    Force-to-torque eff: " << gammaft << endl;

    writerdoc << endl;
    writerdoc << "   -------- Grain options --------" << endl<< endl;
    writerdoc << "           Radial drift: ";

    if (idrift == 0)    writerdoc << "no" << endl;
    else                writerdoc << "yes" << endl;

    if (massorsize == 0)
    {   
        writerdoc << "                 Bounce: ";

        if (ibounce == 0)   writerdoc << "no" << endl;
        else                writerdoc << "yes" << endl;
    }

    writerdoc << "             Disruption: ";

    if (idisrupt == 0)  writerdoc << "no" << endl;
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
    writerdoc << endl;

    writerdoc << "   -------- Grain final state --------" << endl;
    writerdoc << endl;

    for (int i=0; i < ngrains; i++)
    {
        if (i+1 < 10) writerdoc << "            R0" << i+1; 
        else          writerdoc << "            R" << i+1;

        switch (istate[i])
        {   
            case(0):    
            {   
                writerdoc << " still alive" << endl;
                break;
            }
            case(1):    
            {
                writerdoc << " reached maximum size" << endl;
                break;
            }
            case(2):    
            {   
                writerdoc << " was accreted" << endl;
                break;
            }
            case(3):    
            {   
                writerdoc << " was disrupted" << endl;
                break;
            }
        }
    }

    writerdoc.close();
}

/* ------------------------ TOOLS ------------------------*/

void ErrorValue(bool& error, const string& typerror, const string& variable)
{
    cerr << typerror <<": " << variable << endl;
    error = true;
}

string ToStringWithPrecision(const double& value, const int& precision)
{
    ostringstream out;
    out << setprecision(precision) << value;
    return out.str();
}

string OutputFileName(const int& massorsize, const double& Rini, const int& iporosity)
{
    string ms;
    string porosity;

    if (iporosity == 1)  porosity = "por_";
    else                 porosity = "comp_";

    if (massorsize == 0) ms = "M_";
    else                 ms = "S_";

    return "output" + porosity + ms + ToStringWithPrecision(Rini,10).c_str() + ".out";
}

string DisruptFileName(const int& massorsize)
{
    string ms;

    if (massorsize == 0) ms = "M_";
    else                 ms = "S_";

    return "disrupt_param_" + ms + ".out";
}

void ReadVoid(ifstream& reader, int nbvoid)
{
    string blank;
    for (int i = 0; i < nbvoid; i++)
    {   reader >> blank;  }
}

