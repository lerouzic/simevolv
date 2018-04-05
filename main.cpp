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
	string iarchi_file;
    string oarchi_file;
    string ipop_file;
    string opop_file;
	string output_file;
	long int seed;

	// Option parser
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Print help messages")
      ("parameter,p", po::value<string>(&input_file), "Parameter file")
      ("iarchitecture,a", po::value<string>(&iarchi_file), "Input architecture file")
      ("oarchitecture,b", po::value<string>(&oarchi_file), "Output architecture file")
      ("ipopulation,P", po::value<string>(&ipop_file), "Input population file")
      ("opopulation,Q", po::value<string>(&opop_file), "Output population file")
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

	if (vm.count("iarchitecture")) 
	{
		Architecture::load(iarchi_file);
	} else {
		Architecture::initialize(param);
	}
	
	unsigned int global_generation = 0;
	
    Population pop(param); // Generates the genotypes based on the rules defined in the parameter set
    
    if (vm.count("ipopulation")) {
        #ifdef SERIALIZATION_TEXT 
            ifstream popin(ipop_file.c_str(), ios::in);
            popin >> pop; // Copies everything, but phenotypes and fitnesses will be updated. 
            pop.update_param(param); // warning: some rates are stored in the population class, they should be updated. 
            popin.close();
        #else
            assert("Compile the program with a SERIALIZATION flag before using the \"Input population\" option");
        #endif
    }
    
	string next_parfile = string("");
	bool continue_simulation = false;
    
	do { // while continue_simulation

        Fitness::initialize(param);
        Environment::initialize(param);
        Architecture::update_param(param);		

        // The population needs to be up-to-date against the parameter file. 
        pop.update_param(param);

		unsigned int current_pargen = param.getpar(SIMUL_GENER)->GetInt();
		unsigned int intervgen = param.getpar(SIMUL_OUTPUT)->GetInt();
		unsigned int maxgen = param.exists(SIMUL_MAXGEN) ? 
            param.getpar(SIMUL_MAXGEN)->GetInt() : global_generation + current_pargen;
		
        // current_pargen: number of generations to run with the current parameter file
        // intervgen     : number of generations between output events
        // maxgen        : total number of generations before stopping the simulation.
        //                 if not specified in the parameter file, runs to the end of the current parameter
        
        
		// Two counts: 
        //   - inner_generation corresponds to the inner loop (current parameter file)
        //   - global generation corresponds to the outer loop 
        //     (number of generations since the beginning of the simulation)
        // The simulation stops if global_generation == maxgen (max number of generations reached)
        // or if inner_generation == current_pargen (and no NEXTPAR file is provided). 
		for (unsigned int inner_generation = 0; 
            (inner_generation < current_pargen) && (global_generation < maxgen); 
						inner_generation++, global_generation++)
		{
            // At the beginning of the loop, the population is only defined by its genotype.
            pop.update_phenotype();
            pop.update_fitness();
            
            // Step 1: Output if necessary. If inner_generation == 0, this is the initial population
            // in the context of the new parameter file. 
			if ((global_generation == 0) || (global_generation % intervgen == 0))
			{
				pop.write(*pt_output, global_generation);
			}
        
            // Step 2: Run a generation (except at the last generation of the inner loop)
            // note: pop = pop.reproduce(...) does not work
            Population offspring = pop.reproduce(param.getpar(INIT_PSIZE)->GetInt());
            pop = offspring;

		} // end of inner loop
		
        continue_simulation = (param.exists(FILE_NEXTPAR) && (global_generation <= maxgen)); 
        
		if (continue_simulation) {
			next_parfile = param.getpar(FILE_NEXTPAR)->GetString();
            param.erase(FILE_NEXTPAR); // we do not want to keep this if FILE_NEXTPAR is missing in the file
			param.read(next_parfile);
            if (param.exists(FILE_NEXTPAR) 
                && (param.getpar(FILE_NEXTPAR)->GetString() == next_parfile) 
                && (!param.exists(SIMUL_MAXGEN))) {
                // Avoid infinite loop
                continue_simulation = false;
            }
        }
        
    } while (continue_simulation);
    
    // A last output is necessary at the end of the simulation
    pop.update_phenotype();
    pop.update_fitness();    
    pop.write(*pt_output, global_generation);
    
    if (vm.count("parcheck")) 
    {
		param.warning_unused();
		param.warning_multicalls();
	}
	
	if (vm.count("oarchitecture")) { // serializes the current architecture
        Architecture::save(oarchi_file);     
    }
    
    if (vm.count("opopulation")) { // serializes the current population
    #ifdef SERIALIZATION_TEXT                    
        ofstream popout(opop_file.c_str());
        popout << pop;
        popout.close();
    #else
        assert("Compile the program with a SERIALIZATION flag before using the \"Output population\" option");
    #endif
    }
    
	file_out.close();	// This probably does not harm if the file is not open
	return(EXIT_SUCCESS);
}
