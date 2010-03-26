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
		"  -d               Draw progress.\n"\
		"  -S               Run the simulation slowly, so progression can be seen.\n"\
		"  -h               Show this help text.\n";
	exit(1);
}

int main(int argc, char ** argv)
{
	bool automatic=true,slow=false,do_output=false;
	int opt, size = DEFAULT_SIZE, time, counter;
	const char *state_filename="/dev/null";
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
				break;
			case 'S':
				slow = true;
				break;
			case 't':
				automatic = false;  // If a time is specfied don't do auto end detection
				time = atoi(optarg);
				if (time < 1)
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
				state_filename = "/dev/stderr";
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
		/* Open files for writing.  TODO: check for errors */
		if (strcmp(state_filename,"-")==0)
			state_filename = "/dev/stdout";
		if (state_filename != NULL)
			state_out.open(state_filename, fstream::trunc | fstream::out);
#elif OUTPUT_GNUPLOT
		state_out << "set term gif animate\n"\
			"unset key\n"\
			"set view map\n"\
			"set output 'progression.gif'\n"\
			"unset tics\n"\
			"set cbtics (-1, 1)\n"\
			"set cblabel 'Spin Direction'\n";
#endif
	}
	/* Initialise a lattice object */
	Lattice lattice = Lattice(size,J,muH,kT);
	lattice.randomise();
	if (do_output)
	{
		state_out << lattice;
	}
	/* counter is just keeping track of number of steps, ending is more complex */
	for (counter = 0;;counter ++)
	{
		/* This always calls lattice.step().  We will then end the loop if 
		 * we are in automatic mode and the net spin change < 0.1*size. This
		 * condition is a naive attempt at detecting equilibrium.*/
		if ((lattice.step() < 0.1 * size ) && automatic)
			break;
		/* If we're not in automatic mode, just end after time iterations.*/
		if ((!automatic) && (counter >= time))
			break;
		/* slow allows realtime viewing of the equilbriation */
		if (slow)
			sleep(1);
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
	cout << lattice.kT/lattice.J << ' ' << lattice.M() << ' ' << lattice.E() << ' ' << counter << endl;
	return 0;
}
