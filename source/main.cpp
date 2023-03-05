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
        Pamdeas("input.in","No");
    }
    else if (argc == 2)
    {
        if (strcmp(argv[1], "setup") == 0)
        {   
            WriteInputFile();
            cout << "Input file has been written: input.in" << endl;
        }
        else if (strcmp(argv[1], "test") == 0)
        {   
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running Discinit test" << endl;
            Pamdeas("","TestDiscinit");
            cout << "Discinit test complete" << endl << endl;
            cout << "->Running GrowthFragMass test" << endl;
            Pamdeas("","TestGrowthFragMass");
            cout << "GrowthFragMass test complete" << endl << endl;
            cout << "->Running GrowthFragSize test" << endl;
            Pamdeas("","TestGrowthFragSize");
            cout << "GrowthFragSize test complete" << endl << endl;
            cout << "->Running Drift test" << endl;
            Pamdeas("","TestDrift");
            cout << "Drift test complete" << endl;
        }
        else if (strcmp(argv[1], "testdisc") == 0)
        {
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running Discinit test" << endl;
            Pamdeas("","TestDiscinit");
            cout << "Discinit test complete" << endl << endl;
        }
        else if (strcmp(argv[1], "testgrowth") == 0)
        {
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running GrowthFragMass test" << endl;
            Pamdeas("","TestGrowthFragMass");
            cout << "GrowthFragMass test complete" << endl << endl;
            cout << "->Running GrowthFragSize test" << endl;
            Pamdeas("","TestGrowthFragSize");
            cout << "GrowthFragSize test complete" << endl << endl;
        }
        else if (strcmp(argv[1], "testdrift") == 0)
        {
            cout << "\e[1mStarting Pamdeas test suite:\e[0m" << endl << endl;
            cout << "->Running Drift test" << endl;
            Pamdeas("","TestDrift");
            cout << "Drift test complete" << endl;
        }
        else
        {   Pamdeas(string(argv[1]),"No");   }
    }
    else
    {   cerr << "Error: to many arguments passed" << endl;  }

    return 0;
}