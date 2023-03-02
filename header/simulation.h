#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

/* ------------------------ TERMINAL DISPLAY ------------------------*/

void checksize(double& mass, double& filfac,const double& rhos, double& sdustmin);

// Presentation animation
void PresentationAnim();

// Terminal progression / waiting animation
void ProgressionAnim(const double& iteratedvalue, const double endvalue, const int& istate, const string& outputfile);

// Wall time
void Walltime(const double& walltime);

// End animation
void EndAnim();

/* ------------------------------------------------------*/
/* -------------------- MAIN ROUTINE --------------------*/
/* ------------------------------------------------------*/
    
void Pamdeas(const string& input);

