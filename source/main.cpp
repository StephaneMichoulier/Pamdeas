#include <iostream>
#include <iomanip>

#include "../header/simulation.h"
#include "../header/readwrite.h"

using namespace std;


/*------------------------ MAIN ------------------------*/

int main(int argc, char* argv[])
{   
    if (argc == 1)
    {   
        Pamdeas("input.in","../input/","No");
    }
    else if (argc == 2)
    {
        if (strcmp(argv[1],"setup") == 0)
        {   
            WriteInputFile();
            cout << "Input file has been written: ./input/input.in" << endl;
        }
        else if (strcmp(argv[1], "test") == 0)
        {   
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;

            cout << "->Running Discinit test" << endl;
            Pamdeas("","","TestDiscinit");
            cout << "Discinit test complete" << endl << endl;

            cout << "->Running GrowthSize test" << endl;
            Pamdeas("","","TestGrowthSize");
            cout << "GrowthSize test complete" << endl << endl;

            cout << "->Running GrowthMass test" << endl;
            Pamdeas("","","TestGrowthMass");
            cout << "GrowthMass test complete" << endl << endl;

            cout << "->Running Drift test" << endl;
            Pamdeas("","","TestDrift");
            cout << "Drift test complete" << endl << endl;

            cout << "->Running PorosityGrowthSize test" << endl;
            Pamdeas("","","TestPorosityGrowthSize");
            cout << "PorosityGrowthSize test complete" << endl << endl;

            cout << "->Running PorosityGrowthMass test" << endl;
            Pamdeas("","","TestPorosityGrowthMass");
            cout << "PorosityGrowthMass test complete" << endl << endl;

            cout << "->Running Disrupt test" << endl;
            Pamdeas("","","TestDisrupt");
            cout << "Disrupt test complete" << endl;
        }
        else if (strcmp(argv[1], "testdisc") == 0)
        {
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running Discinit test" << endl;
            Pamdeas("","","TestDiscinit");
            cout << "Discinit test complete" << endl << endl;
        }
        else if (strcmp(argv[1], "testgrowth") == 0)
        {
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running GrowthSize test" << endl;
            Pamdeas("","","TestGrowthSize");
            cout << "GrowthSize test complete" << endl << endl;

            cout << "->Running GrowthMass test" << endl;
            Pamdeas("","","TestGrowthMass");
            cout << "GrowthMass test complete" << endl << endl;
        }
        else if (strcmp(argv[1], "testdrift") == 0)
        {
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running Drift test" << endl;
            Pamdeas("","","TestDrift");
            cout << "Drift test complete" << endl;
        }
        else if (strcmp(argv[1], "testporosity") == 0)
        {
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running PorosityGrowthSize test" << endl;
            Pamdeas("","","TestPorosityGrowthSize");
            cout << "PorosityGrowthSize test complete" << endl << endl;

            cout << "->Running PorosityGrowthMass test" << endl;
            Pamdeas("","","TestPorosityGrowthMass");
            cout << "PorosityGrowthMass test complete" << endl << endl;

            cout << "->Running TestPorosityAllMass test" << endl;
            Pamdeas("","","TestPorosityAllMass");
            cout << "TestPorosityAllMass test complete" << endl << endl;
        }
        else if (strcmp(argv[1], "testdisrupt") == 0)
        {
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running Disrupt test" << endl;
            Pamdeas("","","TestDisrupt");
            cout << "Disrupt test complete" << endl;
        }
        else
        {   Pamdeas(string(argv[1]),"../input/","No");
        }
    }
    else if (argc == 3)
    {
        Pamdeas(string(argv[1]),"../"+string(argv[2])+"/","No");
    }
    else
    {   cerr << "Error: to many arguments passed" << endl;  }
    return 0;
}