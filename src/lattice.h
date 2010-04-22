/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */

#ifndef _LATTICE_H_
#define _LATTICE_H_

#include <vector>
#include <ostream>

class Lattice : public std::vector<std::vector<char> > 
{
public:
	Lattice(size_type, double J=DEFAULT_J, double muH=DEFAULT_muH, double kT=DEFAULT_kT);
	void randomise();
	int step();
	double M();
	double E();
private:
	double kT;
	double J;
	double muH;
	size_type sz;
	double exp_lookup[10];
	void initalise_exp_lookup ();
};

std::ostream & operator<<(std::ostream &os, Lattice &lattice);

#endif // _LATTICE_H_
