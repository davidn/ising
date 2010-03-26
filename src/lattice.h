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
	double kT;
	double J;
	double muH;
private:
	size_type sz;
};

std::ostream & operator<<(std::ostream &os, Lattice &lattice);

#endif // _LATTICE_H_
