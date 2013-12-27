#include "main.h"

using namespace std;


int main(int argc, char *argv[])
{
    //string input_file = "/home/runneburger/Models/param.txt"; //pour test via CodeBlocks
    //string output_file = "";

    string input_file = argv[1];      // pour passage en ligne de commande :
    //string output_file = argv[2];

    ParameterSet param(input_file);
    OutputFormat::SetSummary(cout);

    Random::initialize();
    Architecture::initialize(param);
    Fitness::initialize(param);
    Environment::initialize(param);

    Population pop(param);
    int maxgen = param.getpar(SIMUL_GENER)->GetInt();
    for (int generation = 1; generation <= maxgen; generation++)
    {
        Fitness::update_generation(generation);
        if ((generation == 1) || (generation == maxgen) || (generation % param.getpar(SIMUL_OUTPUT)->GetInt() == 0))
        {
            pop.write();
        }
        Population offsp = pop.reproduce();
        if (generation < maxgen)
            pop = offsp;
    }
    int extragen = param.getpar(SIMUL_EXTRA)->GetInt();
    if (extragen > 0)
    {
        Population orig = pop;
        Fitness::update_extra(+1.0);
        pop.update();

        for (int generation = 1; generation <= extragen; generation++)
        {
            // No need to call Fitness::update_generation
            Population offsp = pop.reproduce();
            pop.write();
            if (generation < extragen)
                pop = offsp;
        }
        pop = orig;
        Fitness::update_extra(-1.0);
        pop.update();
        for (int generation = 1; generation <= extragen; generation++)
        {
            // No need to call Fitness::update_generation
            Population offsp = pop.reproduce();
            pop.write();
            if (generation < extragen)
                pop = offsp;
        }
    }

}

