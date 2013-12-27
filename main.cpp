
#include <iostream>
#include <fstream>
#include <string>

#include <boost/program_options.hpp>

#include "main.h"
#include "Parconst.h"

using namespace std;

int main(int argc, char *argv[])
{	
	string input_file;
	string output_file;	
	
	// Option parser
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages") 
      ("parameter,p", po::value<string>(&input_file), "Parameter file") 
      ("output,o", po::value<string>(&output_file), "Output file");
      
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm); //read the command line options
	notify(vm);
	    
    
    if (vm.count("help")) {
		cout << "Command line help" << endl;
		cout << desc << endl;
		return(0); // The program ends here
	}
	
	if (!vm.count("parameter")) {
		cerr << "A parameter file must be provided" << endl;
		cerr << desc << endl;
		return(1);
	}
	
	ostream* pt_output = &cout; // Default: the output goes to std::cout
	ofstream file_out;
	
	if (vm.count("output")) {
		file_out.open(output_file.c_str());
		pt_output = &file_out;
	} 
	
	OutputFormat::SetSummary(*pt_output);
	
	cerr << input_file << endl;
	ParameterSet param(input_file);
    
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

	file_out.close(); // This probably does not harm if the file is not open at all
	return(0);
}

