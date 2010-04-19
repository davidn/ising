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

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <cstdio>
#include <unistd.h>
#include <sys/wait.h>
#include <libgen.h>
#include <cerrno>
#include <vector>
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

static const char* progname;

/* Print some documentation. */
void usage()
{
	cout << "Usage: " <<progname<< " [options]\n"\
		"Simulate a 2D ferromagnet using a monte-carlo ising model.\n\n"\
		"  -t time          Run each temperature run for time iterations.\n"\
		"  -s size          Use a size by size lattice.\n"\
		"  -j n             Use n processes (best set to number of CPUs)\n"\
		"  -J J             Set interation parameter to J.\n"\
		"  -H H             Set external Magnetic field to H.\n"\
		"  -o output.dat    Save the data to output.dat.\n"\
		"  -g graph.png     Draw a graph of magnetisation and energy to graph.png.\n"\
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
	const char *output_filename="ising.dat", *graph_filename=NULL;
	double kTfrom=DEFAULT_kT, kTto=DEFAULT_kT,kT;
	vector<char *> parameters;
	FILE * output;
	char command[255];
	command[0] = '\0';
	strncat(command,dirname(argv[0]),255);
	strncat(command,"/ising",255);
	parameters.push_back(command);
	/* Read command line options. */
	while ((opt = getopt(argc,argv,"g:t:s:J:H:o:1:2:j:n:h?")) != -1)
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
			case 'g':
				graph_filename = optarg;
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

	if ((output = fopen(output_filename, "w")) == NULL)
	{
		perror("Could not open output file");
		exit(1);
	}
	
	/* Do a run at each temperature. We detect kTfrom==kTto later*/
	for (kT=kTfrom; kT <= kTto; kT += (kTto-kTfrom)/num_temps)
	{
		if (children.size() >= parallel)
		{
			tempid = wait(NULL);
			if (tempid == -1)
			{
				perror("Could not wait for child");
				exit(1);
			}
			children.erase(find(children.begin(),children.end(),tempid));
		}
		if((tempid=fork())==0)
		{
			freopen(output_filename,"a",stdout);
			char tempstr[31];
			sprintf(tempstr,"%.15E",kT);
			parameters.push_back("-T");
			parameters.push_back(tempstr);
			parameters.push_back(NULL);
			execv(command,&parameters[0]);
			perror("Could not run subprocess");
			exit(1);
		}
		if (tempid == -1)
		{
			perror("Could not fork()");
			exit(1);
		}
		children.push_back(tempid);
		/* if kTto==kTfrom the kT increment is zero, and we therefore must only
		 * do one run, so break in this case. */
		if (kTto == kTfrom)
			break;
	}

	for(vector<pid_t>::iterator it = children.begin(); it != children.end(); ++it)
	{
		if(waitpid(*it,NULL,0)==-1)
		{
			perror("Could not collect a child (ignoring)");
		}
	}
	
	/* Close the output to ensure it is flushed for gnuplot to read. */
	if(fclose(output)==EOF)
	{
		perror("Could not close ouput file (ignoring)");
	}
	/* Tell gnuplot to draw us a graph.*/
	if (graph_filename != NULL)
	{
		Gnuplot gp;
		gp << \
			"set term png size 1024,768\n"\
			"set output '" << graph_filename << "'\n"\
			"set xlabel 'kT/J'\n"\
			"set ylabel 'Magnetization'\n"\
			"set xtics rotate by -45 add ('Tc(Onsager)' 2.269185314)\n"\
			"set y2label 'Energy'\n"\
			"set ytics nomirror\n"\
			"set y2tics\n"\
			"plot '"<<output_filename << "' u 1:(abs($2)) t 'Magnetization',"\
			"'"<<output_filename << "' u 1:3 t 'Energy' axes x1y2\n";
	}
	
	return 0;
}