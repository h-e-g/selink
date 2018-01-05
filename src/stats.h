/* Copyright 2012 GEH Institut Pasteur */
/* Copyright 2014 Pierre Boutillier */
/* This file is part of Selink. */

/* Selink is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* Selink is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with Selink.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef __STATS_H_
#define __STATS_H_

#include "data_acquisition.h"

double SSD_WP (unsigned long c, unsigned int n);
double SSD_T(unsigned long c1, unsigned int n1,
	     unsigned long c2, unsigned int n2);

double sum_w_from_type(int type, unsigned int popsize);
double omega_from_type(int type, unsigned int popsize, unsigned int i);

double Haplotype_diversitycompute(unsigned int *occurences, unsigned int pop_size);

/* compute polymorphism on the given window from nb_a2 */
unsigned int Polymorphism_counter(unsigned long *nb_a2, unsigned int popsize,
				  unsigned long start, unsigned int len);

/* Update allele frequency on the given window
   BEWARE: it requires the invarient old_start <= new_start <= old_stop <= new_stop.
*/
unsigned long update_AFS(unsigned int *xi, unsigned long *nb_a2_holder, snp_t **snp,
			 unsigned int pop_len, unsigned long max_xi,
			 unsigned long old_start, unsigned long new_start,
			 unsigned long old_stop, unsigned long new_stop);

double *Omega_reader(char *filename, unsigned int snp_nb);

double Theta_compute(int type, unsigned int *ksi, unsigned int popsize);

double Theta2_fromS(unsigned int popsize, unsigned int snpnb);

/* cache harmonic sums from 0 -> pop size */
/* avoids to recompute HS for every value */
double *Cache_HarmonicSums(unsigned long n) ;

/* compute equation (9) from Achaz paper */
/* double Neutrality_test(unsigned int *ksi, unsigned int popsize, double theta1, double theta2);*/
void Neutrality_test_precomputation (unsigned int n, int w1,  int w2,
				     double **W, double *alpha, double *beta);
double Neutrality_test(unsigned int *xi, unsigned long n, double *W, double alpha,
		       double beta, double theta, double theta2);

double *HARMONICSUMS_CACHED;

#endif
