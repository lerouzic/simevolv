// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Parameters.h"
#include "Parconst.h"
#include "Architecture.h"
#include "Fitness.h"
#include "Environment.h"
#include "Population.h"
#include "Random.h"

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;



int main(int argc, char *argv[])
{
	string input_file;
	string archi_file;
	string output_file;
	long int seed;

	// Option parser
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Print help messages")
      ("parameter,p", po::value<string>(&input_file), "Parameter file")
      ("architecture,a", po::value<string>(&archi_file), "Architecture file")
      ("output,o", po::value<string>(&output_file), "Output file")
      ("seed,s", po::value<long int>(&seed), "Seed for the random number generator")
      ("template,t", "Print a template for the parameter file")
      ("parcheck,c", "Warns about inconsistencies in the parameter file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm); //read the command line options
	notify(vm);

	// First thing to do: set the output
	ostream* pt_output = &cout; 		/* Default: the output goes to std::cout */
	ofstream file_out;
	
	if (vm.count("output")) 
	{
		file_out.open(output_file.c_str());
		pt_output = &file_out;
	}
	
	if (vm.count("help")) 
    {
		*pt_output << "Command line help" << endl;
		*pt_output << desc << endl;
		return(EXIT_SUCCESS); // The program ends here
	}

	if (vm.count("template")) 
	{
		ParameterSet pp;
		pp.write(*pt_output);
		return(EXIT_SUCCESS);
	}

	if (!vm.count("parameter")) 
	{
		cerr << "A parameter file must be provided" << endl;
		cerr << desc << endl;
		return(EXIT_FAILURE);
	}

    if (vm.count("seed")) 
    {
		Random::initialize(seed);
	} else {
		Random::initialize();
	}

	ParameterSet param(input_file);

	if (vm.count("architecture")) 
	{
		Architecture::initialize(archi_file);
	} else {
		Architecture::initialize(param);
	}
	
	unsigned int global_generation = 1;
	
	// This is needed for the initial population
	// (explaining the unfortunate code duplication later)
	Fitness::initialize(param);
	Environment::initialize(param);
	
    Population pop(param);
    
	string next_parfile = "";    
    
	do { // while next_parfile == ""
		if (param.exists(FILE_NEXTPAR)) {
			if (next_parfile == param.getpar(FILE_NEXTPAR)->GetString()) {
				// Otherwise it is very easy to loop forever!
				next_parfile = "";
			} else {
				next_parfile = param.getpar(FILE_NEXTPAR)->GetString();
			}
		}
	
		// Partially duplicated code -- probably harmless, but not elegant
		if (global_generation > 1) {
			Fitness::initialize(param);
			Environment::initialize(param);
			Architecture::update_param(param);		
			pop.update_param(param);
		}

		unsigned int maxgen = param.getpar(SIMUL_GENER)->GetInt();
		unsigned int intervgen = param.getpar(SIMUL_OUTPUT)->GetInt();
		    
		for (unsigned int local_generation = 1; local_generation <= maxgen; local_generation++, global_generation++)
		{
			Fitness::fluctuate(local_generation);
			if ((global_generation == 1) || (global_generation == maxgen) || (global_generation % intervgen == 0))
			{
				pop.write(*pt_output, global_generation);
			}
        
			if ((local_generation < maxgen) || (next_parfile != "")) {
				// no need to compute a new population if the simulation is over
				Population offsp = pop.reproduce(param.getpar(INIT_PSIZE)->GetInt());
				pop = offsp;
			}
	        // Quick and dirty: send some debugging/detailed information into a specific file
	        /*if (generation==maxgen) {
				ofstream tmp_out("debug.txt");
				ostream & out = tmp_out;
				pop.write_debug(out);
				tmp_out.close();
			}*/
		}
		if (next_parfile != "") {
			param.read(next_parfile);
		}
    } while (next_parfile != "");
    
    if (vm.count("parcheck")) 
    {
		param.warning_unused();
		param.warning_multicalls();
	}
	
	Architecture::terminate();    // This allows to write the architecture content in a file if necessary
	file_out.close();	// This probably does not harm if the file is not open
	return(EXIT_SUCCESS);
}
