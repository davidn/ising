/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */

#ifndef _LATTICE_H_
#define _LATTICE_H_

#include <vector>
#include <ostream>

class Lattice
{
public:
	Lattice(size_t, double J=DEFAULT_J, double muH=DEFAULT_muH, double kT=DEFAULT_kT);
	void randomise();
	int step();
	double M();
	double E();
private:
	std::vector<std::vector<char> > lattice;
	double kT;
	double J;
	double muH;
	size_t sz;
	double exp_lookup[10];
	void initalise_exp_lookup ();

	friend std::ostream & operator<<(std::ostream &os, Lattice &lattice);
};

std::ostream & operator<<(std::ostream &os, Lattice &lattice);

#endif // _LATTICE_H_
