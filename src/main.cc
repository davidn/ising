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
#define DEFAULT_kT 1
#define DEFAULT_muH 0
#define DEFAULT_J 1

#define OUTPUT_DOTS

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include "../config.h"
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

#define RAND() (rand()/RAND_MAX)
#define boltzmann(energy) exp(-energy/kT)

double J = DEFAULT_J;
double muH = DEFAULT_muH;
double kT;

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

class Lattice : public vector<vector<char> >
{
	public:
		Lattice(size_type num);
		void randomise();
		int step();
		double M();
		double E();
	private:
		size_type sz;
};

Lattice::Lattice(size_type num): sz(num)
{
	this->resize(num, vector<char>::vector(num));
}

void Lattice::randomise()
{
	for(vector<vector<char> >::iterator it = this->begin(); it<this->end(); ++it)
	{
		for(vector<char>::iterator at = it->begin(); at<it->end(); ++at)
		{
			*at = (rand() > RAND_MAX/2) ? -1 : 1;
		}
	}
}

/* Step the lattice, return the net change in total spins */
int Lattice::step()
{
	int change = 0,i,j;
	double d_E;
	for (int num=0; num < sz*sz; num++)
	{
		i= rand() * sz / RAND_MAX;
		j= rand() * sz / RAND_MAX;
		/* Using 
		 d_E = J * Sum over neighbours ( initial*neighbour - final*neighbour) + muH * ( initial - final )
		 Note that final = - initial, giving
		 d_E = 2 J * Sum over neighbours(initial*neighbour) + 2 muH * initial*/
		d_E = 2 * J * this->at(i)[j] * ( this->at((i+1)%sz)[j]
		                                +this->at((i-1)%sz)[j]
		                                +this->at(i)[(j+1)%sz]
		                                +this->at(i)[(j-1)%sz] );
		/*d_E increases if spin is initially positive and goes to negative*/
		d_E += 2 * muH * this->at(i)[j];
		if (d_E < 0 || boltzmann(d_E) > RAND())
		{
			change += 2 * (this->at(i)[j] = - this->at(i)[j]);
		}
		else
		{
			this->at(i)[j] = this->at(i)[j];
		}
	}
	return change;
}

double Lattice::E()
{
	double Eret = 0;
	for (int i=0; i < sz; i++)
	{
		for (int j = 0; j < sz; j++)
		{
			Eret -= J * this->at(i)[j] * ( this->at((i+1)%sz)[j]
			                                +this->at((i-1)%sz)[j]
			                                +this->at(i)[(j+1)%sz]
			                                +this->at(i)[(j-1)%sz] );
			Eret -= muH * this->at(i)[j];
		}
	}
	return Eret / (sz * sz);
}

double Lattice::M()
{
	double Mret = 0;
	for(vector<vector<char> >::iterator it = this->begin(); it<this->end(); ++it)
	{
		for(vector<char>::iterator at = it->begin(); at<it->end(); ++at)
		{
			Mret += *at;
		}
	}
	return Mret / (sz * sz);
}

ostream & operator<<(std::ostream &os, Lattice &lattice)
{
#ifdef OUTPUT_DOTS
	for(vector<vector<char> >::iterator it = lattice.begin(); it<lattice.end(); ++it)
	{
		for(vector<char>::iterator at = it->begin(); at < it->end(); ++at)
		{
			os << ((*at > 0) ? '.' : 'O');
		}
		os << endl;
	}
#else
#endif
}

int main(int argc, char ** argv)
{
	bool automatic=true,slow=false,verbose=false;
	int opt, change, size = DEFAULT_SIZE, time, num_temps=1, counter;
	const char *state_filename="/dev/null", *output_filename="ising.dat";
	double kTfrom=DEFAULT_kT, kTto=DEFAULT_kT;
	Gnuplot gp;
	fstream state_out, output;
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
				automatic = false;
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

	if (kTfrom>kTto)
	{
		double kTtmp = kTto;
		kTto = kTfrom;
		kTfrom = kTtmp;
		cerr << "Swapping temperatures 1 and 2." << endl;
	}
	
	if (strcmp(output_filename,"-")==0)
		output_filename = "/dev/stdout";
	if (output_filename != NULL)
		output.open(output_filename, fstream::trunc | fstream::out);

	if (strcmp(state_filename,"-")==0)
		state_filename = "/dev/stdout";
	if (state_filename != NULL)
		state_out.open(state_filename, fstream::trunc | fstream::out);
	
	Lattice lattice = Lattice(size);

	for (kT=kTfrom; kT <= kTto; kT += (kTto-kTfrom)/num_temps)
	{
		lattice.randomise();
		if (verbose)
			cout << "Temperature: " << kT << endl;
		state_out << lattice;
		for (counter = 0;;counter ++)
		{
			if ((lattice.step() < 0.1 * size ) && automatic)
				break;
			if ((!automatic) && (counter >= time))
				break;
			if (slow)
				sleep(1);
			state_out << endl << lattice;
		}
		output << kT << " " << lattice.M() << " " << lattice.E() << " " << counter << endl;
		if (kTto == kTfrom)
		{
			break;
		}
		if (verbose)
			cout << "  Took " << counter << " steps" << endl;
	}
	output.close();
	gp << \
		"set term png size 1024,768\n"\
		"set output 'graph.png'\n"\
		"set xlabel 'kT/J'\n"\
		"set ylabel 'Magnetization'\n"\
		"set y2label 'Energy'\n"\
		"set ytics nomirror\n"\
		"set y2tics\n"\
		"plot '"<<output_filename << "' u 1:(abs($2)) t 'Magnetization',"\
		"'"<<output_filename << "' u 1:3 t 'Energy' axes x1y2\n";
	
	return 0;
}
