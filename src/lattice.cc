/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <vector>
#include <cstdlib>
#include <cmath>
#include "lattice.h"

using namespace std;

/* exp lookup stores the value of exp(-d_E/kT) for each possible arrangement
 * of the current lattice point and it's neighbours.  The least significant
 * bit is the current direction of the current point (1=+ve, 0=-ve). The other
 * three bits are the number of positive going neighbours.
 *
 * (i & 0x1)?1:-1 is thus the spin of the current point
 * ((i & 0xe)-4) is thus the net spin of the neigbours
 */
#define MY_SPIN(i) ((i & 0x1)?1:-1)
#define NEIGH_SPIN(i) ((i & 0xe)-4)

/* Initialise a num x num lattice */
Lattice::Lattice(size_type num, double J, double muH, double kT)
:sz(num), kT(kT), J(J), muH(muH), vector<vector<char> >(num, vector<char>(num))
{
	this->initalise_exp_lookup ();
}

/* Randomise the lattice to +/-1 only, uniformly spread */
void Lattice::randomise()
{
	for(auto it: this)
	{
		for(auto& at: it)
		{
			at = (rand() > RAND_MAX/2) ? -1 : 1;
		}
	}
}

/* Step the lattice, return the net change in total spins */
int Lattice::step()
{
	/* change is a running total of net spin change, useful for equilibrium
	 * detection.*/
	int change = 0,i,j;
	unsigned char lookup;
	/* Do sz*sz random points rather than all... this prevents the entire
     * lattice being flipped at high temperature.  Perhaps we should just do 1
     * point and have Lattice.step() called more times? */
	for (int num=0; num < sz*sz; num++)
	{
		/* We only need 1 call to rand() if sz < 2^15 which is probably safe to
		   assume as a 32768x32768 grid would take a LONG time to run through */
		int rand_bits = rand(); //31 bits of randomness
		i= (rand_bits & 0x7FFF) * sz / 0x8000; // lowest 15 bits
		j= (rand_bits>>15 & 0x7FFF) * sz / 0x8000; //next 15 bits
		lookup = (*this)[(i+1)%sz][j]
			+(*this)[(i-1)%sz][j]
			+(*this)[i][(j+1)%sz]
			+(*this)[i][(j-1)%sz] + 4
			+ ((*this)[i][j] > 0 ? 1:0 );
		if (exp_lookup[lookup] * RAND_MAX > rand())
			change += 2 * ((*this)[i][j] = - (*this)[i][j]);
		/* NB: we flip the spin           ^^^    here  */
	}
	return change;
}

/* return energy density of the lattice. */
double Lattice::E()
{
	double Eret = 0;
	for (int i=0; i < sz; i++)
	{
		for (int j = 0; j < sz; j++)
		{
			/*Only sum in 1 direction per dimension to prevent double counting*/
			Eret -= J * (*this)[i][j] * ( (*this)[(i+1)%sz][j]
			                             +(*this)[i][(j+1)%sz] );
			Eret -= muH * (*this)[i][j];
		}
	}
	return Eret / (sz * sz);
}

/* return magnetization of the lattice */
double Lattice::M()
{
	double Mret = 0;
	for(auto it: this)
	{
		for(auto at: it)
		{
			Mret += at;
		}
	}
	return Mret / (sz * sz);
}

void Lattice::initalise_exp_lookup()
{
	unsigned char i;
	double d_E;
	for (i=0;i<10;i++)
	{
		/* Using 
		 d_E = J * Sum over neighbours ( initial*neighbour - final*neighbour) + muH * ( initial - final )
		 Note that final = - initial, giving
		 d_E = 2 J * initail * Sum of neighbours + 2 muH * initial*/
		d_E = 2 * muH * MY_SPIN(i)
			+ 2 * J * MY_SPIN(i) * NEIGH_SPIN(i);
		exp_lookup[i] = d_E < 0 ? 1 : exp(-d_E/kT);
	}
}

vector<vector<char> >::iterator & begin(Lattice* const &lattice)
{
	return begin(lattice);
}

vector<vector<char> >::iterator & end(Lattice* const &lattice)
{
	return end(lattice);
}

/* Allow easy lattice printing */
ostream & operator<<(std::ostream &os, Lattice &lattice)
{
	/* OUTPUT_DOTS mode is for console viewing... it prints . and 0 depending
     * on spin direction.*/
#ifdef OUTPUT_DOTS
	for(it: lattice)
	{
		for(at: it)
		{
			os << ((*at > 0) ? '.' : 'O');
		}
		os << endl;
	}
#elif OUTPUT_GNUPLOT
	os << "splot '-' matrix with image\n";
	for(auto it: lattice)
	{
		for(auto at: it)
		{
			os << (int)at << ' ';
		}
		os << endl;
	}
	os << "e\ne\n";
#endif
		
	return os;
}

