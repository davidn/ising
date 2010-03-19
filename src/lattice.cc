/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
/*
 * ising
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

#include <vector>
#include <cstdlib>
#include <cmath>
#include "lattice.h"

using namespace std;

Lattice::Lattice(size_type num, double J, double muH, double kT)
:sz(num), kT(kT), J(J), muH(muH), vector<vector<char> >(num, vector<char>::vector(num))
{
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
		if (d_E < 0 || exp(-d_E/kT) > rand()/RAND_MAX)
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

