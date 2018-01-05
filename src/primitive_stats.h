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

#ifndef __PRIMITIVE_STATS_H_
#define __PRIMITIVE_STATS_H_

#include "data_acquisition.h"

double **create_double_array_for_results(samples *smpls, unsigned long nlocis);
unsigned long **create_long_array_for_results(samples *smpls, unsigned long nlocis);
void free_double_results_array(samples *smpls, double** t);
void free_long_results_array(samples *smpls, unsigned long** t);

void fill_half_ihh_pi_divs_and_nbd_arrays_one_pop
(int verbose, double min_value, int rev, unsigned long pop_id,
 unsigned int pop_len, samples *smpls, snp_t **snp,
 unsigned long nlocis, double *ihh_a1, double *ihh_a2,
 unsigned long *nb_a2_holder, double *pi_a1, double *pi_a2,
 unsigned long *diversities, unsigned long *haplo_sums);

#endif
