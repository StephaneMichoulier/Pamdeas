#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "../header/constant.h"
#include "../header/readwrite.h"
#include "../header/general_functions.h"

using namespace std;


/* ------------------------ READING ------------------------*/

void CheckData(const int& massorsize, const double& tend, const int& stepmethod, const double& step, const int& profile,// ->
               const double& Mstar, const double& Mdisk, const double& Rin, const double& Rout, const double& R0,// -> 
               const double& dustfrac0, const double& H0R0, const double& p, const double& q, const double& alpha, // ->
               const int& iporosity, const double& sizeini, const double& phiini, const double& a0, const double& rhos,// ->
               const double& youngmod0, const double& esurf, const double& Yd0, const double& Ydpower, const int& idrift,// -> 
               const int& ibounce, const int& ifrag, const int& ibr, const int& ibump, const double& vfragi,// -> 
               const int& constvfrag, const double& philim, const double& philimbounce, const double& limsize,// -> 
               const double& Rbump, const double& dustfracmax, const double& bumpwidth, const double& bumpheight,// ->
               const int& ngrains, const vector <double>& Rini)
{
    bool error = false;

    if (massorsize != 0 && massorsize != 1)     Error(error, "massorsize != 0 & != 1");
    
    if (tend <= 0)                              Error(error, "tend <= 0");

    if (stepmethod != 0 && stepmethod != 1 && stepmethod != 2)                       
                                                Error(error, "stepmethod != 0 & != 1 & != 2");

    if (stepmethod != 2)
    {
        if (step <= 0)                          Error(error, "step <= 0");
    }

    if (profile != 0 && profile != 1)           Error(error, "profile != 0 & != 1");

    if (Mstar <= 0)                             Error(error, "Mstar <= 0");

    if (Mdisk <= 0)                             Error(error, "Mdisk <= 0");

    if (Rin <= 0)                               Error(error, "Rin <= 0");

    if (Rout <= 0)                              Error(error, "Rout <= 0");

    if (Rin >= Rout)                            Error(error, "Rin >= Rout");

    if (R0 <= 0)                                Error(error, "R0 <= 0");

    if (dustfrac0 <= 0 || dustfrac0 >= 1)
    {   
        if (dustfrac0 >= 1 )                    Error(error, "dustfrac0 >= 1");
        else                                    Error(error, "dustfrac0 <= 0");
    }

    if (H0R0 <= 0)                              Error(error, "H0R0 <= 0");

    if (alpha <= 0)                             Error(error, "alpha <= 0");

    if (ibr != 0 && ibr != 1)                   Error(error, "ibr != 0 & != 1");

    if (ibump != 0 && ibump != 1)               Error(error, "ibump != 0 & != 1");

    if (iporosity != 0 && iporosity != 1)       Error(error, "iporosity != 0 & != 1");

    if (sizeini <= 0)                           Error(error, "sizeini <= 0");

    if (phiini <= 0 || phiini > 1)              
    {   
        if(phiini > 1)                          Error(error, "phiini > 1");
        else                                    Error(error, "phiini <= 0");
    }

    if (a0 <= 0 || a0 > sizeini)                                
    {   
        if(a0 <= 0)                             Error(error, "a0 <= 0");
        else                                    Error(error, "a0 > sizeini");
    }

    if (rhos <= 0)                              Error(error, "rhos <= 0");

    if (youngmod0 <= 0)                         Error(error, "Youngmod0 <= 0");

    if (esurf <= 0)                             Error(error, "Esurf <= 0");

    if (Yd0 <= 0)                               Error(error, "Yd0 <= 0");

    if (idrift != 0 && idrift != 1)             Error(error, "idrift != 0 & != 1");

    if (ibounce != 0 && ibounce != 1)           Error(error, "ibounce != 0 & != 1");

    if (ibounce == 1 && massorsize == 1)        Error(error, "bounce is not included with ds/dt model");

    if (ifrag != 0 && ifrag != 1 && ifrag != 2) Error(error, "ifrag != 0 & != 1 & != 2");

    if (vfragi < 0)                             Error(error, "vfragi < 0");

    if (constvfrag != 0 && constvfrag != 1)     Error(error, "constvfrag != 0 & != 1");

    if (philim < 0 || philim > 1)           
    {   
        if(philim > 1)                          Error(error, "philim > 1");
        else                                    Error(error, "philim < 0");
    }

    if (philimbounce < 0 || philimbounce > 1)   
    {   
        if(philimbounce > 1)                    Error(error, "philimbounce > 1");
        else                                    Error(error, "philimbounce < 0");
    }

    if (limsize <= sizeini)                     Error(error, "limsize <= sizeini");

    if (Rbump < Rin || Rbump > Rout)
    {   
        if (Rbump < Rin )                       Error(error, "Rbump < Rin");
        else                                    Error(error, "Rbump > Rout");
    }

    if (dustfracmax < dustfrac0 || dustfracmax >= 1)
    {   
        if (dustfracmax >= 1 )                  Error(error, "dustfracmax >= 1");
        else                                    Error(error, "dustfracmax < dustfrac0");
    }

    if (bumpwidth <= 0)                         Error(error, "bumpwidth <= 0");

    if (bumpheight <= 0)                        Error(error, "bumpheight <= 0");

    if (ngrains < 0)                            Error(error, "ngrain < 0");
    else
    {
        if (error == false)
        {   for (int i=0; i<ngrains; i++)
            {
                if (Rini[i] < Rin)              Error(error, "Rini < Rin for particle " + to_string(i+1));
                if (Rini[i] > Rout)             Error(error, "Rini > Rout for particle " + to_string(i+1));
            }
        }
    }

    if (error == true)                          exit(1);

} 


void ReadFile(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, double& Mstar, double& Mdisk,// ->
              double& Rin, double& Rout, double& R0, double& dustfrac0, double& H0R0, double& p, double& q,// ->
              double& alpha, int& iporosity, double& sizeini, double& phiini, double& a0, double& rhos, double& youngmod0,// ->
              double& esurf, double& Yd0, double& Ydpower, int& idrift, int& ibounce, int& ifrag, int& ibr, int& ibump,// ->
              double& vfragi, int& constvfrag, double& philim, double& philimbounce, double& limsize, double& Rbump,// ->
              double& dustfracmax, double& bumpwidth, double& bumpheight, int& ngrains, vector <double>& Rini)
{   
    ifstream Reader("input.in");
    if (!Reader)
   	{   cout << "Fatal error, input.in is missing " << endl;
        WriteInputFile();
        cout << "Input file has been written: input.in" << endl;
        exit(1);
    }

    ReadVoid(Reader,6);     Reader >> massorsize;
    ReadVoid(Reader,10);    Reader >> tend;
    ReadVoid(Reader,5);     Reader >> stepmethod;
    ReadVoid(Reader,9);     Reader >> step;
    ReadVoid(Reader,17);    Reader >> profile;
    ReadVoid(Reader,7);     Reader >> Mstar; 
    ReadVoid(Reader,5);     Reader >> Mdisk;
    ReadVoid(Reader,5);     Reader >> Rin; 
    ReadVoid(Reader,5);     Reader >> Rout;
    ReadVoid(Reader,5);     Reader >> R0; 
    ReadVoid(Reader,5);     Reader >> dustfrac0; 
    ReadVoid(Reader,8);     Reader >> H0R0; 
    ReadVoid(Reader,6);     Reader >> p;
    ReadVoid(Reader,3);     Reader >> q;
    ReadVoid(Reader,2);     Reader >> alpha; 
    ReadVoid(Reader,6);     Reader >> ibr;
    ReadVoid(Reader,5);     Reader >> ibump; 

    
    ReadVoid(Reader,8);     Reader >> iporosity; 
    ReadVoid(Reader,4);     Reader >> sizeini; 
    ReadVoid(Reader,5);     Reader >> a0; 
    ReadVoid(Reader,5);     Reader >> rhos; 
    
    ReadVoid(Reader,7);     Reader >> idrift;
    ReadVoid(Reader,5);     Reader >> ibounce; 
    ReadVoid(Reader,5);     Reader >> ifrag; 
    ReadVoid(Reader,6);     Reader >> vfragi; 
    ReadVoid(Reader,5);     Reader >> limsize;

    ReadVoid(Reader,15);    Reader >> phiini; 
    ReadVoid(Reader,5);     Reader >> youngmod0;
    ReadVoid(Reader,11);    Reader >> esurf;
    ReadVoid(Reader,13);    Reader >> Yd0; 
    ReadVoid(Reader,9);     Reader >> Ydpower;
    ReadVoid(Reader,9);     Reader >> constvfrag;
    ReadVoid(Reader,7);     Reader >> philim;
    ReadVoid(Reader,8);     Reader >> philimbounce; 

    ReadVoid(Reader,14);    Reader >> Rbump;
    ReadVoid(Reader,6);     Reader >> dustfracmax;
    ReadVoid(Reader,5);     Reader >> bumpwidth;
    ReadVoid(Reader,9);     Reader >> bumpheight;

    ReadVoid(Reader,15);    Reader >>ngrains;

    ReadVoid(Reader,7);

    Rini.resize(ngrains);

    for (int i = 0; i < ngrains; i++)
    {
        ReadVoid(Reader,2); Reader >> Rini[i];
    }

    CheckData(massorsize,tend,stepmethod,step,profile,Mstar,Mdisk,Rin,Rout,R0,dustfrac0,H0R0,p,q,alpha,iporosity,sizeini,
              phiini,a0,rhos,youngmod0,esurf,Yd0,Ydpower,idrift,ibounce,ifrag,ibr,ibump,vfragi,constvfrag,philim,
              philimbounce,limsize,Rbump,dustfracmax,bumpwidth,bumpheight,ngrains,Rini);

    // Convert values to SI units
    Mstar = MsolToKg(Mstar);
    Mdisk *= Mstar;

    if (stepmethod == 2)    step = 1.;

	Reader.close();
}


/* ------------------------ WRITING ------------------------*/


void WriteInputFile()
{
    ofstream writerinput;

	writerinput.open("input.in");

    writerinput << endl;
    writerinput << "#-Use dm/dt or ds/dt" << endl;
    writerinput << "massorsize = 0          >(dm/dt = 0, ds/dt = 1)" << endl << endl;;

    writerinput << "#-Time options" << endl;
    writerinput << "      tend = 1.0e6      >End time (yrs)" << endl;
    writerinput << "stepmethod = 1          >(Fixe=0, Fraction of orbital period=1, Adaptative dt=2)" << endl;
    writerinput << "      step = 1          >(stepmethod=0 -> step in yrs, stepmethod=1 -> step in fraction of orbital time)" << endl << endl;  

    writerinput << "#-Disk profiles" << endl;
    writerinput << "     profile  = 0       >(0=no, 1=yes)" << endl << endl;

    writerinput << "#-Gas disk properties" << endl << endl;

    writerinput << "     Mstar = 1.         >Star mass (Msol)" << endl;
    writerinput << "     Mdisk = 0.01       >Disk mass (Mstar)" << endl;
    writerinput << "      Rin  = 3.         >Inner radius (AU)" << endl;
    writerinput << "      Rout = 300.       >Outer radius (AU)" << endl;
    writerinput << "      Rref = 100.       >Reference radius (AU)" << endl;	
    writerinput << " dustfrac0 = 0.01       >Dust to gas ratio at Rref" << endl;
    writerinput << "     H0/R0 = 0.05       >H/R at Rref" << endl;
    writerinput << "   p index = 1.5        " << endl;
    writerinput << "   q index = 0.75       " << endl;
    writerinput << "     alpha = 0.01       >Turbulence Shakura & Sunyaev" << endl;
    writerinput << "       ibr = 0          >Back-reaction (0=no, 1=yes)" << endl; 
    writerinput << "     ibump = 0          >Pressure bump (0=no, 1=yes)" << endl << endl;

    writerinput << "#-Dust properties" << endl << endl;

    writerinput << " iporosity = 1          >(compact=0, porous=1)" << endl;
    writerinput << "   sizeini = 1.0e-7     >Initial size (m)" << endl;
    writerinput << "        a0 = 1.0e-7     >Monomer size (m)" << endl;
    writerinput << "      rhod = 917        >Dust density (kg/m³)" << endl << endl;

    writerinput << "#-Grain options" << endl << endl;
    writerinput << "    idrift = 0          >Drift (0=no, 1=yes)" << endl;
    writerinput << "   ibounce = 0          >Bounce (0=no, 1=yes)" << endl;
    writerinput << "     ifrag = 0          >Fragmentation (0=no, 1=hardfrag, 2=smoothfrag)" << endl;
    writerinput << "     vfrag = 15         >Fragmentation threshold (m/s)" << endl;
    writerinput << "   maxsize = 1.0e3      >Maximum size to stop the simulation" << endl << endl;

    writerinput << "#-Porosity properties, Available if iporosity = 1" << endl << endl;
    writerinput << "    phiini = 1          >Initial filling factor" << endl;
    writerinput << " Youngmod0 = 9.4e9      >Young Modulus for ice [Yamamoto et al. 2014] (Pa)" << endl;
    writerinput << "     Esurf = 7.3e-2     >Suface energy for ice grains J/m² [Yamamoto et al. 2014] (Pa)" << endl;
    writerinput << "       Yd0 = 9.8e6      >Dynamic compression resistance constant for ice (Pa)" << endl;
    writerinput << "   Ydpower = 4          >Dynamic compression resistance for ice [Mellor, 1975]" << endl;
    writerinput << "constvfrag = 1          >Constant fragmentation threshold (0=no, 1=yes)" << endl;
    writerinput << "    philim = 0.01       >Filling factor dynamic compression resistance limit" << endl;
    writerinput << "   philimb = 0.3        >Filling factor bounce limit" << endl << endl;

    writerinput << "#-Pressure bump option, Available if ibump = 1" << endl << endl;
    writerinput << "      Rbump = 6.5       >Pressure bump radius (AU)" << endl;
    writerinput << "dustfracmax = 0.1       >Max dustfrac possible" << endl;
    writerinput << "  bumpwidth = 1.1       >Bump half width at half maximum (AU)" << endl;
    writerinput << " bumpheight = 300       >Bump height (Surface density at Rref)" << endl << endl;

    writerinput << "#-Number of dust grains and initial radii" << endl << endl;

    writerinput << "   ngrains = 10         >Number of dust grains" << endl << endl;
	
    writerinput << "   Initial radii (AU)" << endl;
    writerinput << "       R01 = 5" << endl;
    writerinput << "       R02 = 10" << endl;
    writerinput << "       R03 = 20" << endl;
    writerinput << "       R04 = 50" << endl;
    writerinput << "       R05 = 75" << endl;
    writerinput << "       R06 = 100" << endl;
    writerinput << "       R07 = 150" << endl;
    writerinput << "       R08 = 200" << endl;
    writerinput << "       R09 = 250" << endl;
    writerinput << "       R10 = 300" << endl << endl;

writerinput.close();
}

void WriteProfileFile(ofstream& outputprofile, const double& Rprofile, const double& hg, const double& cg, const double& sigma,
                      const double& rhog, const double& dustfrac, const double& Pg, const double& T)
{
    WriteValue(outputprofile,10,6,Rprofile);
    WriteValue(outputprofile,10,6,hg);
    WriteValue(outputprofile,10,6,cg);
    WriteValue(outputprofile,10,6,sigma);
    WriteValue(outputprofile,10,6,rhog);
    WriteValue(outputprofile,10,6,dustfrac);
    WriteValue(outputprofile,10,6,Pg);
    WriteValue(outputprofile,10,6,T);
    outputprofile << endl;
}

void WriteProfileColumns()
{
    ofstream writecol("diskprofilescolumns.txt");
    writecol << "R" << endl
             << "Hg" << endl
             << "cg" << endl
             << "sigma" << endl
             << "rhog" << endl
             << "dustfrac" << endl
             << "Pg" << endl
             << "T" << endl;
}

void WriteOutputFile(ofstream& outputfile, const double& t, const double& Rf, const double& massf, const double& phif,// ->
                     const double& sizef, const double& St, const double& cg, const double& sigma,// ->
                     const double& rhog, const double& dustfrac, const double& vrel, const double& omegak,// ->
                     const double& drdt, const double& dvardt, const int& iregime)
{   
    WriteValue(outputfile,10,6,t);
    WriteValue(outputfile,10,6,Rf);
    WriteValue(outputfile,10,6,massf);
    WriteValue(outputfile,10,6,phif);
    WriteValue(outputfile,10,6,sizef);
    WriteValue(outputfile,10,6,St);
    WriteValue(outputfile,10,6,cg);
    WriteValue(outputfile,10,6,sigma);
    WriteValue(outputfile,10,6,rhog);
    WriteValue(outputfile,10,6,dustfrac);
    WriteValue(outputfile,10,6,vrel);
    WriteValue(outputfile,10,6,omegak);
    WriteValue(outputfile,10,6,drdt);
    WriteValue(outputfile,10,6,dvardt);
    WriteValue(outputfile,1,1,iregime);
    outputfile << endl;
}

void WriteOutputColumns(const double& massorsize)
{
    ofstream writecol("outputcolumns.txt");
    writecol << "t" << endl
             << "R" << endl
             << "mass" << endl
             << "phi" << endl
             << "size" << endl
             << "St" << endl
             << "cg" << endl
             << "sigma" << endl
             << "rhog" << endl
             << "dustfrac" << endl
             << "vrel" << endl
             << "omegaK"<< endl
             << "drdt" << endl;

    if (massorsize == 0)
    {   writecol << "dmdt" << endl;  }
    else
    {   writecol << "dsdt" << endl;  }
    
    writecol << "regime" << endl;
}

void WriteInitFile(const int& massorsize, const double& tend, const int& stepmethod, const double& dt, const double& Mstar,// ->
                   const double& Mdisk, const double& Rin, const double& Rout, const double& R0, const double& Rbump,// -> 
                   const double& dustfrac0, const double& H0R0, const double& p, const double& q, const double& alpha,// ->
                   const int& iporosity, const double& sizeini, const double& phiini, const double& a0, const double& rhos,// -> 
                   const int& idrift, const int& ibounce, const int& ifrag, const int& ibr, const int& ibump,// ->
                   const double& vfragi, const int& ngrains, const double& sigma0, const double& rhog0,// ->
                   const double& cg0, const double& runningtime)
{
    ofstream writerdoc;
    string dash;

	writerdoc.open("Initials_Conditions_" + to_string_with_precision(R0,10) + "AU.txt");

	writerdoc << endl<< "     Running Time : " << runningtime << " s" << endl << endl;

    if (R0 < 10)        dash = "";
    else if (R0 >= 10)  dash = "-";
    else if (R0 >= 100) dash = "--";
    else                dash = "---";

    writerdoc << "------------------------------------------------------" << dash << endl;
	writerdoc << "! ---------- Initials Conditions at " << R0 << " AU ---------- !" << endl;
    writerdoc << "------------------------------------------------------" << dash << endl << endl;
	writerdoc << "    Gas surface density: " << sigma0 << " kg/m²" << endl;
	writerdoc << "           Scale height: " << H0R0*R0 << " AU" << endl;
	writerdoc << "        Gas sound speed: " << cg0 << " m/s" << endl;
	writerdoc << "            Gas density: " << rhog0 << " kg/m³" << endl;
	writerdoc << "        Gas temperature: " << T(R0,q,R0,cg0) << " K" << endl;
	writerdoc << "           Gas pressure: " << Pg(rhog0,cg0) << " Pa" << endl;
	writerdoc << "     Gas mean free path: " << Lambda(rhog0,cg0) << " m" << endl << endl << endl;

    writerdoc << "------------------------------------------" << endl;
    writerdoc << "! ---------- Input Parameters ---------- !" << endl;
    writerdoc << "------------------------------------------" << endl << endl;
    writerdoc << "             Model used: ";

    if (massorsize == 0)     writerdoc << "dm/dt" << endl;
    else                     writerdoc << "ds/dt" << endl;

    writerdoc << "         Dust particles: " << ngrains << endl<< endl;

    writerdoc << "         Simulated time: " << tend << " yrs" << endl;

    writerdoc << "   Time-stepping method: ";
    switch (stepmethod)
    {
        case (0):
        {   
            writerdoc << "Fixe" << endl;  
            writerdoc << "                     dt: " << dt << " yrs" << endl << endl;
            break; 
        }
        case (1):    
        {   
            writerdoc << "Fraction of orbital period" << endl;  
            writerdoc << "                     dt: " << dt << endl << endl;
            break;  
        }
        case (2):    
        {
            writerdoc << "Adaptative" << endl << endl;    
            break; 
        }
    } 

    writerdoc << "   -------- Gas disk properties --------" << endl<< endl;
    writerdoc << "              Star mass: " << KgToMsol(Mstar) << " Msol" << endl;
    writerdoc << "              Disk mass: " << KgToMsol(Mdisk) << " Msol" << endl;
    writerdoc << "           Inner radius: " << Rin << " AU" << endl;
    writerdoc << "           Outer radius: " << Rout << " AU" << endl;
    writerdoc << "       Reference radius: " << R0 << " AU" << endl;
    writerdoc << "H/R at reference radius: " << H0R0 << endl;
    writerdoc << "  Dust/gas ratio (Rref): " << dustfrac0 << endl;
    writerdoc << "                p index: " << p << endl;
    writerdoc << "                q index: " << q << endl;
    writerdoc << "       Alpha turbulence: " << alpha << endl << endl;

    writerdoc << "          Back-reaction: ";

    if (ibr == 0)       writerdoc << "no" << endl;
    else                writerdoc << "yes" << endl;

    writerdoc << "          Pressure bump: ";

    if (ibump == 0)     writerdoc << "no" << endl;
    else
    {            
        writerdoc << "yes" << endl;
        writerdoc << "   Pressure bump radius: " << Rbump << " AU" << endl << endl;
    }

    writerdoc << endl;
    writerdoc << "   -------- Dust properties --------" << endl << endl;

    writerdoc << "                  Grain: ";

    if (iporosity == 0)  writerdoc << "compact" << endl;
    else                 writerdoc << "porous" << endl;

    writerdoc << "     Initial grain size: " << sizeini << " m" << endl;
    if (iporosity == 0) writerdoc << " Initial filling factor: " << phiini << endl;
    writerdoc << "           Monomer size: " << a0 << " m" << endl;
    writerdoc << " Dust intrinsic density: " << rhos << " kg/m³" << endl;
    if (ifrag !=0)      writerdoc << "Fragmentation threshold: " << vfragi << " m/s" << endl;

    writerdoc << endl;
    writerdoc << "   -------- Grain options --------" << endl<< endl;
    writerdoc << "           Radial drift: ";

    if (idrift == 0)    writerdoc << "no" << endl;
    else                writerdoc << "yes" << endl;

    writerdoc << "                 Bounce: ";

    if (ibounce == 0)   writerdoc << "no" << endl;
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
    
    writerdoc.close();
}

/* ------------------------ TOOLS ------------------------*/

void Error(bool& error, const string& variable)
{
    cout << "Error, " << variable << endl;
    error = true;
}

string to_string_with_precision(const double& value, const int& precision)
{
    ostringstream out;
    out << setprecision(precision) << value;
    return out.str();
}

string FileName(const int& massorsize, const double& Rini, const int& iporosity)
{
    string ms;
    string porosity;

    if (iporosity == 1)  porosity = "por_";
    else                 porosity = "comp_";

    if (massorsize == 0) ms = "M_";
    else                 ms = "S_";

    return "output" + porosity + ms + to_string_with_precision(Rini,10).c_str() + ".out";
}

void WriteValue(ostream& writer, const int& width, const double& precision, const double& value)
{
    writer << setw(width)<< setprecision(precision) << value << " ";
}

void ReadVoid (ifstream& reader,int nbvoid)
{
    string blank;
    for (int i = 0; i < nbvoid; i++)
    {   reader>>blank;  }
}

