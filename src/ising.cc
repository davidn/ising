/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
/*
 * main.cc
 * Copyright (C) David Newgas 2010 <dn271@cam.ac.uk>
 * 
 * ising is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ising is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <ctime>
#include "../config.h"
#include "lattice.h"
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

/* Print some documentation. */
void usage()
{
	cout << "Usage: " PACKAGE " [options]\n"\
		"Simulate a 2D ferromagnet using a monte-carlo ising model.\n\n"\
		"  -T temp          Temperature to simulate.\n"\
		"  -J J             Set interation parameter to J.\n"\
		"  -H H             Set external Magnetic field to H.\n"\
		"  -t time          Run each temperature run for time iterations.\n"\
		"  -s size          Use a size by size lattice.\n"\
		"  -d file.gif      Draw progress to progress.gif.\n"\
		"  -S               Run the simulation slowly, so progression can be seen.\n"\
		"  -h               Show this help text.\n";
	exit(1);
}

int main(int argc, char ** argv)
{
	bool automatic=true,slow=false,do_output=false;
	int opt, size = DEFAULT_SIZE, iterations, counter, stabilised=0, condition, change;
	const char *state_filename=NULL;
	double J=DEFAULT_J, muH=DEFAULT_muH, kT=DEFAULT_kT;
#ifdef OUTPUT_DOTS
	fstream state_out;
#elif OUTPUT_GNUPLOT
	Gnuplot state_out;
	//ostream &state_out = cout;
#endif
	/* Read command line options. */
	while ((opt = getopt(argc,argv,"St:s:J:H:dT:h?")) != -1)
	{
		switch(opt)
		{
			case 's':
				size = atoi(optarg);
				if (size < 1)
			{
				cout << "Lattice too small!" << endl;
				exit(2);
			}
				break;
			case 'T':
				kT= atof(optarg);
				if (kT < 0)
			{
				cout << "Temperature below 0 is meaningless!" << endl;
				exit(2);
			}
				break;
			case 'S':
				slow = true;
				break;
			case 't':
				automatic = false;  // If a time is specfied don't do auto end detection
				iterations = atoi(optarg);
				if (iterations < 1)
			{
				cout << "Must have at least 1 step!" << endl;
				exit(2);
			}
				break;
			case 'J':
				J = atof(optarg);
				break;
			case 'H':
				muH = atof(optarg);
				break;
			case 'd':
				do_output = true;
				state_filename = optarg;
				break;
			case 'h':
			case '?':
			default:
				usage();
		}
	}

	if (do_output)
	{
#ifdef OUTPUT_DOTS
		if (strcmp(state_filename,"-")==0)
			state_filename = "/dev/stdout";
		if (state_filename != NULL)
		{
			state_out.open(state_filename, fstream::trunc | fstream::out);
			if(!state_out)
			{
				cerr << "Error opening output" << endl;
			}
		}
#elif OUTPUT_GNUPLOT
		state_out << "set term gif animate" <<(slow? " delay 15" : "") << " \n"\
			"unset key\n"\
			"set view map\n"\
			"set output 'progression.gif'\n"\
			"unset tics\n"\
			"set cbtics (-1, 1)\n"\
			"set cblabel 'Spin Direction'\n";
#endif
	}

	/* Lets have a different seed each time (Ignoring errors)*/
	srand(time(NULL));
	
	/* Initialise a lattice object */
	Lattice lattice = Lattice(size,J,muH,kT);
	lattice.randomise();
	if (do_output)
	{
		state_out << lattice;
	}

	/* This is the expected variance of the system in equilibrium.  It assumes
	 * the energy requirement for a single point to move from equilibrium is 
	 * 4J 
	 */
	condition = exp(-4*J/kT) * size * size;
	/* We need this to be nonzero to ensure we terminate*/
	if (condition==0)
		condition = 1;
	/* counter is just keeping track of number of steps, ending is more complex */
	for (counter = 0;;counter ++)
	{
		/* Step the lattice */
		change = lattice.step();
		if (automatic)
		{
			/* Our automatic mode terminates when the square of the change in 
			 * net spin (ie variance) is less than our precomputed condition for
			 * a given (5) number of iterations. */
			if(change*change < condition)
				stabilised ++;
			else
				stabilised = 0;
			if (stabilised > 5)
				break;
		}
		else
		/* If we're not in automatic mode, just end after time iterations.*/
			if (counter >= iterations)
				break;
#ifdef OUTPUT_DOTS
		/* slow allows realtime viewing of the equilbriation. */
		if (slow && do_output)
			sleep(1);
#endif
		/* Print the state after each iteration. */
		if (do_output)
		{
#ifdef OUTPUT_DOTS
			state_out << '\f' << lattice;
#elif OUTPUT_GNUPLOT
			state_out << lattice;
#endif
		}
	}
	cout << kT/J << ' ' << lattice.M() << ' ' << lattice.E() << ' ' << counter << endl;
	return 0;
}
