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

#define DEFAULT_SIZE 100

#define OUTPUT_DOTS

#include <iostream>
#include <fstream>
#include <unistd.h>
#include "lattice.h"
#include "../config.h"
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

/* Print some documentation. */
void usage()
{
	cout << "Usage: " PACKAGE " [options]\n"\
		"Simulate a 2D ferromagnet using a monte-carlo ising model.\n\n"\
		"  -S               Run the simulation slowly, so progression can be seen.\n"\
		"  -t time          Run each temperature run for time iterations.\n"\
		"  -s size          Use a size by size lattice.\n"\
		"  -J J             Set interation parameter to J.\n"\
		"  -H H             Set external Magnetic field to H.\n"\
		"  -o output.dat    Save the data to output.dat.\n"\
		"  -d               Draw progress.\n"\
		"  -1 temp          Lowest temperature to simulate.\n"\
		"  -2 temp          Highest temperature to simulate.\n"\
		"  -n number        Number of temperature steps to take.\n"\
		"  -V               Be verbose.\n"\
		"  -h               Show this help text.\n";
	exit(1);
}

int main(int argc, char ** argv)
{
	bool automatic=true,slow=false,verbose=false;
	int opt, change, size = DEFAULT_SIZE, time, num_temps=1, counter;
	const char *state_filename="/dev/null", *output_filename="ising.dat";
	double J=DEFAULT_J, muH=DEFAULT_muH, kTfrom=DEFAULT_kT, kTto=DEFAULT_kT;
	Gnuplot gp;
	fstream state_out, output;
	/* Read command line options. */
	while ((opt = getopt(argc,argv,"St:s:J:H:o:d1:2:n:Vh?")) != -1)
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
			case 'n':
				num_temps = atoi(optarg);
				if (num_temps < 1)
				{
					cout << "Must have at least 1 temperature!" << endl;
					exit(2);
				}
				break;
			case '1':
				kTfrom = atof(optarg);
				if (kTfrom < 0)
				{
					cout << "Temperature below 0 is meaningless!" << endl;
					exit(2);
				}
				break;
			case '2':
				kTto = atof(optarg);
				if (kTto < 0)
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
			case 'o':
				output_filename = optarg;
				break;
			case 'g':
				state_filename = "-";
				break;
			case 'V':
				verbose = true;
				break;
			case 'h':
			case '?':
			default:
				usage();
		}
	}
	/* We need from and to in the right order or later a loop will go on forever */
	if (kTfrom>kTto)
	{
		double kTtmp = kTto;
		kTto = kTfrom;
		kTfrom = kTtmp;
		cerr << "Swapping temperatures 1 and 2." << endl;
	}

	/* Open files for writing.  TODO: check for errors */
	if (strcmp(output_filename,"-")==0)
		output_filename = "/dev/stdout";
	if (output_filename != NULL)
		output.open(output_filename, fstream::trunc | fstream::out);

	if (strcmp(state_filename,"-")==0)
		state_filename = "/dev/stdout";
	if (state_filename != NULL)
		state_out.open(state_filename, fstream::trunc | fstream::out);

	/* Initialise a lattice object */
	Lattice lattice = Lattice(size,J,muH);

	/* Do a run at each temperature. We detect kTfrom==kTto later*/
	for (lattice.kT=kTfrom; lattice.kT <= kTto; lattice.kT += (kTto-kTfrom)/num_temps)
	{
		if (verbose)
			cout << "Temperature: " << lattice.kT << endl;
		/* We want to start with a random lattice, and have it printed.*/
		lattice.randomise();
		state_out << lattice;
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
			state_out << '\f' << lattice;
		}
		/* After each temperature run, print the temperature (units of Kb/J),
		 * magnetization, energy and number of iterations taken. */
		output << lattice.kT/lattice.J << " " << lattice.M() << " " << lattice.E() << " " << counter << endl;
		/* if kTto==kTfrom the kT increment is zero, and we therefore must only
		 * do one run, so break in this case. */
		if (kTto == kTfrom)
			break;
		/* Print some fairly usless info.*/
		if (verbose)
			cout << "  Took " << counter << " steps" << endl;
	}
	/* Close the output to ensure it is flushed for gnuplot to read. */
	output.close();
	/* Tell gnuplot to draw us a graph.*/
	gp << \
		"set term png size 1024,768\n"\
		"set output 'graph.png'\n"\
		"set xlabel 'kT/J'\n"\
		"set ylabel 'Magnetization'\n"\
		"set xtics rotate by -45 add ('Tc(Onsager)' 2.269)\n"\
		"set y2label 'Energy'\n"\
		"set ytics nomirror\n"\
		"set y2tics\n"\
		"plot '"<<output_filename << "' u 1:(abs($2)) t 'Magnetization',"\
		"'"<<output_filename << "' u 1:3 t 'Energy' axes x1y2\n";
	
	return 0;
}
