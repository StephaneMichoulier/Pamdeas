#ifndef READWRITE_H_INCLUDED
#define READWRITE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


/* ------------------------ READING ------------------------*/

// Check data from the input file
void CheckData(const int& massorsize, const double& tend, const int& stepmethod, const double& step, const int& profile,// ->
               const double& Mstar, const double& Mdisk, const double& Rin, const double& Rout, const double& R0,// -> 
               const double& dustfrac0, const double& H0R0, const double& p, const double& q, const double& alpha, // ->
               const int& iporosity, const double& sizeini, const double& phiini, const double& a0, const double& rhos,// ->
               const double& youngmod0, const double& esurf, const double& Yd0, const double& Ydpower, const int& idrift,// -> 
               const int& ibounce, const int& disrupt, const int& ifrag, const int& ibr, const int& ibump,// ->
               const double& vfragi, const int& constvfrag, const double& philim, const double& philimbounce,// ->
               const double& limsize, const double& Rbump, const double& dustfracmax, const double& bumpwidth,// ->
               const double& bumpheight, const int& ngrains, const vector <double>& Rini);

// Read the input file
void ReadFile(int& massorsize, double& tend, int& stepmethod, double& step, int& profile, double& Mstar, double& Mdisk,// ->
              double& Rin, double& Rout, double& R0, double& dustfrac0, double& H0R0, double& p, double& q,// ->
              double& alpha, int& iporosity, double& sizeini, double& phiini, double& a0, double& rhos, double& youngmod0,// ->
              double& esurf, double& Yd0, double& Ydpower, int& idrift, int& ibounce, int& idisrupt, int& ifrag, int& ibr,// ->
              int& ibump, double& vfragi, int& constvfrag, double& philim, double& philimbounce, double& limsize,// -> 
              double& Rbump, double& dustfracmax, double& bumpwidth, double& bumpheight, int& ngrains, vector <double>& Rini);


/* ------------------------ WRITING ------------------------*/

// Write the input file if it does not exist
void WriteInputFile();

// Write disk profiles data
void WriteProfileFile(ofstream& outputfile, const double& Rprofile, const double& hg, const double& cg, const double& sigma,// ->
                      const double& rhog, const double& dustfrac, const double& Pg, const double& T);

// Write disk profiles columns
void WriteProfileColumns();

// Write computed data in the output files
void WriteOutputFile(ofstream& outputfile, const double& t, const double& Rf, const double& massf, const double& phif,// ->
                     const double& sizef, const double& St, const double& cg, const double& sigma,// ->
                     const double& rhog, const double& dustfrac, const double& vrel, const double& omegak,// ->
                     const double& drdt, const double& dvardt, const int& iregime);

// Write output files columns
void WriteOutputColumns(const double& massorsize);

// Write the initials conditions
void WriteInitFile(const int& massorsize, const double& tend, const int& stepmethod, const double& dt, const double& Mstar,// ->
                   const double& Mdisk, const double& Rin, const double& Rout, const double& R0, const double& Rbump,// -> 
                   const double& dustfrac0, const double& H0R0, const double& p, const double& q, const double& alpha,// ->
                   const int& iporosity, const double& sizeini, const double& phiini, const double& a0, const double& rhos,// -> 
                   const int& idrift, const int& ibounce, const int& ifrag, const int& ibr, const int& ibump,// ->
                   const int& idisrupt, const double& vfragi, const int& ngrains, const double& sigma0, const double& rhog0,// ->
                   const double& cg0, const double& runningtime);               


/* ------------------------ TOOLS ------------------------*/

// Print error
void Error(bool& error, const string& variable);

// Print type = double without decimals
string to_string_with_precision(const double& value, const int& precision);

// Set the output file name
string FileName(const int& massorsize, const double& Rini, const int& iporosity);

// Write a value with a given column width and precision
void WriteValue(ostream& writer, const int& width, const double& precision, const double& value);

// Read blank space in the input file
void ReadVoid (ifstream& reader,int nbvoid);

#endif // READWRITE_H_INCLUDED