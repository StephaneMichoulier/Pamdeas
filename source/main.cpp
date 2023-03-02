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
        Pamdeas("input.in");    
    }
    else if (argc == 2)
    {
        if (strcmp(argv[1], "setup") == 0)
        {   
            WriteInputFile();
            cout << "Input file has been written: input.in" << endl;
        }
        else
        {   Pamdeas(string(argv[1]));   }
    }
    else
    {   cerr << "Error: to many arguments passed" << endl;  }

    return 0;
}