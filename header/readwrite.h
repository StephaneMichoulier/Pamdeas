#ifndef READWRITE_H_INCLUDED
#define READWRITE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


/* ------------------------ READING ------------------------*/

// Read the input file
void ReadFile(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, double& mstar, double& mdisk,// ->
              double& Rin, double& Rout, double& R0, double& dustfrac0, double& H0R0, double& p, double& q, double& alpha,// ->
              int& iporosity, double& sizeini, double& filfacini, double& a0, double& rhos, double& youngmod0, double& esurf,// ->
              double& Yd0, double& Ydpower, int& idrift, int& ibounce, int& idisrupt, int& ifrag, int& ibr, int& ibump,// ->
              double& gammaft, double& vfragi, int& constvfrag, double& filfaclim, double& filfacbnc, double& limsize,// ->
              double& Rbump, double& dustfracmax, double& bumpwidth, double& bumpheight, int& ngrains, vector <double>& Rini,// ->
              vector <int>& istate);

// Check data from the input file
void CheckData(const int& massorsize, const double& tend, const int& stepmethod, const double& step, const int& profile,// ->
               const double& mstar, const double& mdisk, const double& Rin, const double& Rout, const double& R0,// ->
               const double& dustfrac0, const double& H0R0, const double& p, const double& q, const double& alpha, // ->
               const int& iporosity, const double& sizeini, const double& filfacini, const double& a0, const double& rhos,// ->
               const double& youngmod0, const double& esurf, const double& Yd0, const double& Ydpower, const int& idrift,// ->
               const int& ibounce, const int& idisrupt, const int& ifrag, const int& ibr, const int& ibump,// ->
               const double& gammaft, const double& vfragi, const int& constvfrag, const double& filfaclim,// ->
               const double& filfacbnc, const double& limsize, const double& Rbump, const double& dustfracmax,// ->
               const double& bumpwidth, const double& bumpheight, const int& ngrains, const vector <double>& Rini);
               

/* ------------------------ WRITING ------------------------*/

// Write the input file if it does not exist
void WriteInputFile();

// Write disk profiles data
void WriteProfileFile(ofstream& outputfile, const double& Rprofile, const double& hg, const double& cg, const double& sigma,// ->
                      const double& rhog, const double& dustfrac, const double& pg, const double& T);

// Write disk profiles columns
void WriteProfileHeader();

// Write computed data in output files
void WriteOutputFile(ofstream& outputfile, const double& t, const double& Rf, const double& massf, const double& filfacf,// ->
                     const double& sizef, const double& St, const double& cg, const double& sigma,// ->
                     const double& rhog, const double& dustfrac, const double& vrel, const double& omegak,// ->
                     const double& drdt, const double& dvardt, const int& iregime);

void WriteOutputHeader(const double& massorsize);

// Write data in the disrupt file
void WriteDisruptFile(ofstream& outputfile, const double& R, const double& massf, const double& filfacf, const double& sizef, 
                      const double& st, const double& vrel, const double& freqspin, const double& tensilestress,
                      const double& gammaft, const double& alpha, const double& a0);

// Write disrupt files columns
void WriteDisruptHeader();

// Write the initials conditions
void WriteInitFile(const int& massorsize, const double& tend, const int& stepmethod, const double& dt, const double& mstar,// ->
                   const double& mdisk, const double& Rin, const double& Rout, const double& R0, const double& Rbump,// ->
                   const double& dustfrac0, const double& H0R0, const double& p, const double& q, const double& alpha,// ->
                   const int& iporosity, const double& sizeini, const double& filfacini, const double& a0, const double& rhos,// ->
                   const int& idrift, const int& ibounce, const int& ifrag, const int& ibr, const int& ibump,// ->
                   const int& idisrupt, const double& vfragi, const int& ngrains, const double& sigma0, const double& rhog0,// ->
                   const double& cg0, const vector <int>& istate, const double& runningtime);               


/* ------------------------ TOOLS ------------------------*/

// Print error
void ErrorValue(bool& error, const string& variable);

// Check type of input values
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
