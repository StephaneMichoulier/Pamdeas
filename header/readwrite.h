#ifndef READWRITE_H_INCLUDED
#define READWRITE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


/* ------------------------ READING ------------------------*/

// Read the input file
void ReadFile(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, int& isetdens, int& isettemp, int& ismooth,//->
              double& Rin, double& Rout, double& R0, double& mstar, double& mdisc, double& sigma0, double& hg0R0, double& T0,//->
              double& dustfrac0, double& p, double& q, double& alpha, int& ibr, int& ibump, double& Rbump, double& dustfracmax,//->
              double& bumpwidth, double& bumpheight, int& iporosity, double& sizeini, double& a0, double& rhos, int& idrift,//->
              int& ibounce, int& idisrupt, int& ifrag, double& vfragi, int& ieros, double& ejectasize, double & cohacc, int& icomp,//->
              double& maxsize, int& isnow, double& Rsnow, double& vfragin, double& vfragout, double& youngmod0, double& esurf, double& Yd0,//->
              double& Ydpower, int& constvfrag, double& filfaclim, double& filfacbnc, double& gammaft, int& disrupteq, double& weirmod,//->
              int& ngrains, vector <double>& Rini, vector <int>& istate, const string& input);

// Check data from the input file
void CheckData(const int& massorsize, const double& tend, const int& stepmethod, const double& step, const int& profile,//->
               const int& isetdens, const int& isettemp, const int& ismooth, const double& Rin, const double& Rout, const double& R0,//->
               const double& mstar, const double& mdisc, const double& sigma0, const double& hg0R0, const double& T0, const double& dustfrac0,//->
               const double& p, const double& q, const double& alpha, const int& ibr, const int& ibump, const double& Rbump,//->
               const double& dustfracmax, const double& bumpwidth, const double& bumpheight, const int& iporosity, const double& sizeini,//->
               const double& a0, const double& rhos, const int& idrift, const int& ibounce, const int& idisrupt, const int& ifrag,//->
               const double& vfragi, const int& ieros, const double& ejectasize, const double& cohacc, const int& icomp, const double& maxsize,//-> 
               const int& isnow, const double& Rsnow, const double& vfragin, const double& vfragout, const double& youngmod0, const double& esurf,//->
               const double& Yd0, const double& Ydpower, const int& constvfrag, const double& filfaclim, const double& filfacbnc,//->
               const double& gammaft, const int& disrupteq, const double& weirmod, const int& ngrains, const vector <double>& Rini,//->
               const vector <int>& istate);
               

/* ------------------------ WRITING ------------------------*/

// Write the input file if it does not exist
void WriteInputFile();

// Write disc profiles columns
void WriteProfileHeader(ofstream& outputprofile);

// Write disc profiles data
void WriteProfileFile(ofstream& outputprofile, const double& Rprofile, const double& hg, const double& cg, const double& sigma,//->
                      const double& rhog, const double& dustfrac, const double& pg, const double& T);


// Write output columns
void WriteOutputHeader(ofstream& outputfile, const double& massorsize);

// Write computed data in output files
void WriteOutputFile(ofstream& outputfile, const double& t, const double& Rf, const double& massf, const double& filfacf,//->
                     const double& sizef, const double& st, const double& cg, const double& sigma, const double& rhog,//->
                     const double& dustfrac, const double& vrel, const double& deltav, const double& omegak,//->
                     const double& drdt, const double& dvardt, const int& dragreg, const int& porreg);

// Write disrupt files columns
void WriteDisruptHeader(ofstream& outputfile);

// Write data in the disrupt file
void WriteDisruptFile(ofstream& outputfile, const double& R, const double& massf, const double& filfacf, const double& sizef, 
                      const double& st, const double& vrel, const double& freqspin, const double& tensilestress,
                      const double& gammaft, const double& alpha, const double& a0);


// Write the initials conditions
void WriteResultsFile(const int& massorsize, const double& tend, const int& stepmethod, const double& dt, const int& isetdens, const int& isettemp,//->
                   const int& ismooth, const double& Rin, const double& Rout, const double& R0, const double& mstar, const double& mdisc,//->
                   const double& sigma0, const double& hg0, const double& T0, const double& dustfrac0, const double& rhog0, const double& cg0,//->
                   const double& p, const double& q, const double& alpha, const int& ibr, const int& ibump, const double& Rbump,//->
                   const int& iporosity, const double& sizeini, const double& a0, const double& rhos, const int& idrift, const int& ibounce,//->
                   const int& idisrupt, const int& ifrag, const double& vfragi, const int& ieros, const double& ejectasize, const int& icomp,//->
                   const double& gammaft, const int& disrupteq, const int& isnow, const double& vfragin, const double& vfragout,//->
                   const double& Rsnow, const int& ngrains, const vector <double>& Rini, const  vector <int>& istate, const double& walltime);              


/* ------------------------ TOOLS ------------------------*/

// Print error
void ErrorValue(bool& error, const string& typerror, const string& variable);

// Check type of input values, can only be done in the header file
template <typename T>
void CheckType(ifstream& inputreader, T& outputvalue, const string& name)
{
    string inputvalue;
    inputreader >> inputvalue;

    if (typeid(T) == typeid(double))
    {   
        try
        {
             outputvalue = stod(inputvalue);
        }
        catch(exception& e)
        {
            cerr << "Error: " << name << " is not of type <double>" << endl;
            exit(1);
        }
    }
    if (typeid(T) == typeid(int))
    {   
        try
        {
            if ((!stod(inputvalue) && inputvalue != "0" && inputvalue != "0.") || (stod(inputvalue)-stoi(inputvalue) != 0))
            {   
                cerr << "Error: " << name << " is not of type <int>" << endl;
                exit(1);
            }
            else
            {   
                outputvalue = stoi(inputvalue); 
            }
        }
        catch(exception& e)
        {
            cerr << "Error: " << name << " is not of type <int>" << endl;
            exit(1);
        }
    }
}

// Print type = double without decimals
string ToStringWithPrecision(const double& value, const int& precision);

// Set output files name
string OutputFileName(const int& massorsize, const double& Rini, const int& iporosity);

// Set disrupt file name
string DisruptFileName(const int& massorsize);

// Write a value with a given column width and precision
template <typename T>
void WriteValue(ostream& writer, const int& width, const double& precision, const T& value)
{
    writer << setw(width) << setprecision(precision) << value << " ";
}

// Read blank space in the input file
void ReadVoid(ifstream& reader, int nbvoid);

#endif // READWRITE_H_INCLUDED
