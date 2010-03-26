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

#include <cstdio>
#include <unistd.h>
#include <sys/wait.h>
#include <libgen.h>
#include <cerrno>
#include <vector>
#include "../config.h"
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

/* Print some documentation. */
void usage()
{
	cout << "Usage: " PACKAGE " [options]\n"\
		"Simulate a 2D ferromagnet using a monte-carlo ising model.\n\n"\
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
	int opt, num_temps=1,i=0,parallel=1;
	vector<pid_t> children;
	pid_t tempid;
	const char *output_filename="ising.dat";
	double kTfrom=DEFAULT_kT, kTto=DEFAULT_kT,kT;
	vector<char *> parameters;
	Gnuplot gp;
	FILE * output;
	parameters.push_back("");
	/* Read command line options. */
	while ((opt = getopt(argc,argv,"t:s:J:H:o:1:2:j:n:h?")) != -1)
	{
		switch(opt)
		{
			case 's':
				if (atoi(optarg) < 1)
				{
					cout << "Lattice too small!" << endl;
					exit(2);
				}
				parameters.push_back("-s");
				parameters.push_back(optarg);
				break;
			case 'n':
				num_temps = atoi(optarg);
				if (num_temps < 1)
				{
					cout << "Must have at least 1 temperature!" << endl;
					exit(2);
				}
				break;
			case 'j':
				parallel = atoi(optarg);
				if (num_temps < 1)
				{
					cout << "Must have at least 1 subprocess!" << endl;
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
			case 't':
				if (atoi(optarg) < 1)
				{
					cout << "Must have at least 1 step!" << endl;
					exit(2);
				}
				parameters.push_back("-t");
				parameters.push_back(optarg);
				break;
			case 'J':
				parameters.push_back("-J");
				parameters.push_back(optarg);
				break;
			case 'H':
				parameters.push_back("-H");
				parameters.push_back(optarg);
				break;
			case 'o':
				output_filename = optarg;
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
		output = fopen(output_filename, "w");

	char command[255];
	command[0] = '\0';
	strncat(command,dirname(argv[0]),255);
	strncat(command,"/ising",255);
	
	/* Do a run at each temperature. We detect kTfrom==kTto later*/
	for (kT=kTfrom; kT <= kTto; kT += (kTto-kTfrom)/num_temps)
	{
		if (children.size() >= parallel)
			children.erase(find(children.begin(),children.end(),wait(NULL)));
		if((tempid=fork())==0)
		{
			freopen(output_filename,"a",stdout);
			char tempstr[31];
			sprintf(tempstr,"%.15E",kT);
			parameters.push_back("-T");
			parameters.push_back(tempstr);
			parameters.push_back(NULL);
			execv(command,&parameters[0]);
			perror("Running subprocess");
			exit(1);
		}
		children.push_back(tempid);
		/* if kTto==kTfrom the kT increment is zero, and we therefore must only
		 * do one run, so break in this case. */
		if (kTto == kTfrom)
			break;
	}

	for(vector<pid_t>::iterator it = children.begin(); it != children.end(); ++it)
		waitpid(*it,NULL,0);
	
	/* Close the output to ensure it is flushed for gnuplot to read. */
	fclose(output);
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