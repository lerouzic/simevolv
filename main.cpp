// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "main.h"

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;



int main(int argc, char *argv[])
{
	string input_file;
	string output_file;
	long int seed;

	// Option parser
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Print help messages")
      ("parameter,p", po::value<string>(&input_file), "Parameter file")
      ("output,o", po::value<string>(&output_file), "Output file")
      ("seed,s", po::value<long int>(&seed), "Seed for the random number generator")
      ("template,t", "Print a template for the parameter file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm); //read the command line options
	notify(vm);

	// First thing to do: set the output
	ostream* pt_output = &cout; 	/* Default: the output goes to std::cout */
	ofstream file_out;
	if (vm.count("output")) {
		file_out.open(output_file.c_str());
		pt_output = &file_out;
	}
	OutputFormat::SetSummary(*pt_output);

    if (vm.count("help")) {
		*pt_output << "Command line help" << endl;
		*pt_output << desc << endl;
		return(0); // The program ends here
	}

	if (vm.count("template")) {
		ParameterSet pp;
		pp.write(*pt_output);
		return(0);
	}

	if (!vm.count("parameter")) {
		cerr << "A parameter file must be provided" << endl;
		cerr << desc << endl;
		return(1);
	}

	ParameterSet param(input_file);

    if (vm.count("seed")) {
		Random::initialize(seed);
	} else {
		Random::initialize();
	}

    Architecture::initialize(param);
    Fitness::initialize(param);
    Environment::initialize(param);

    Population pop(param);
    int maxgen = param.getpar(SIMUL_GENER)->GetInt();
    for (int generation = 1; generation <= maxgen; generation++)
    {
        Fitness::fluctuate(generation);
        if ((generation == 1) || (generation == maxgen) || (generation % param.getpar(SIMUL_OUTPUT)->GetInt() == 0))
        {
            pop.write(generation);
        }
        Population offsp = pop.reproduce();
        if (generation < maxgen)
            pop = offsp;
    }
    
	file_out.close(); // This probably does not harm if the file is not open
	return(0);
}

