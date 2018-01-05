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

#define NA "NA"

void fprintf_if_number(FILE *f, const char *format,double x);

void print_beginning_of_header(FILE *OUT) ;

void print_output_header(FILE *OUT, char *pop_name, unsigned int pop_len,
			 int ocomp, char *ocustom, int snp_stats) ;

void print_excluded_header(FILE *OUT) ;

void print_beginning_of_line(FILE *OUT, unsigned long core, unsigned int winlen_left,
			     unsigned int winlen_right, snp_t **snp, int a2_is_ancestral,
			     int a2_is_MAJ) ;
void print_one_snp_stats(FILE *OUT, FILE *EXCLUDED,  snp_t **snp,
			 unsigned long core, int a2_is_ancestral, int a2_is_MAJ,
			 unsigned int pop_len, unsigned long *nb_a2_holder,
			 double *ihh_a1, double *ihh_a2, int out_of_window,
			 double *pi_a1, double *pi_a2) ;
void print_not_computed_line(FILE *OUT, int ocomp, char *ocustom, int snp_stats) ;

void print_files_one_pop(int verbose,int ocomp, char *ocustom, int snp_stats,
			 FILE *OUT, FILE *EXCLUDED, char *pop_name,
			 unsigned int *xi, snp_t **snp, unsigned long nlocis,
			 unsigned int pop_len, unsigned long *nb_a2_holder,
			 unsigned long *diversities, unsigned long *haplo_sums,
			 double *ihh_a1, double *ihh_a2, double *pi_a1, double *pi_a2) ;
// --- interpop ---
void print_FST_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
		       unsigned int *pop_sizes, unsigned long **nb_a2_holder,
		       unsigned int nb_pop);
/** guillaume: global Fst */
void print_GLOB_FST_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
		       unsigned int *pop_sizes, unsigned long **nb_a2_holder,
		       unsigned int nb_pop);
void print_Xpehh_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
			 double **ihh_a1, double **ihh_a2, unsigned int nb_pop) ;
/** guillaume: Delta IHH  */
void print_Deltaihh_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
			 double **ihh_a1, double **ihh_a2, unsigned int nb_pop) ;
/** guillaume: Delta DAF  */
void print_DeltaDAF_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
		       unsigned int *pop_sizes, unsigned long **nb_a2_holder,
		       unsigned int nb_pop) ;


void print_interpop_header(FILE *OUT, samples *smpls, int snp_stats) ;
void print_not_computed_interpop(FILE *OUT, unsigned int nb_pop, int snp_stats,
				 unsigned long core, snp_t **snp,
				 int is_a2_ancestral) ;
void print_interpop_file(int verbose, char *outfile, samples *smpls, snp_t **snp,
			 unsigned long nlocis,  unsigned int* pop_sizes,
			 double **ihh_a1, double **ihh_a2,
			 unsigned long **nb_a2_holder, int snp_stats);
