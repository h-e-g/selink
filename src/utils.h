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

#ifndef __UTILS_H_
#define __UTILS_H_

#include "data_acquisition.h"
#define MAYBE(tab,ind) ((tab)?(tab)[ind]:NULL)

/* toggle sequence to binary representation ::  inplace modification */
int Seq2bin(char *ref, char  **seq);

/* Reads the SNP positions selected from the positions file */
snp_pos *Positions_reader (char *pos_file, int num_chr);

/* Best SNP between two defining the ends of the optimal window size */
unsigned long Best_snp (snp_t **snp_def, unsigned long snp_ref_pos,
			unsigned long num_bef, unsigned long num_aft,
			unsigned int winsize);

/* Searches the SNP stop of the window */
unsigned long Snp_search_right (snp_t **snp_def, unsigned long rank_snp_ref,
				unsigned int winsize, unsigned long nlocis,
				int use_best);

/* Searches the SNP start of the window */
unsigned long Snp_search_left (snp_t **snp_def, unsigned long rank_snp_ref,
			       unsigned int winsize, int use_best);

/* Windows definition by SNP */
unsigned long Windowmaker_by_snp (snp_t **snp_def, unsigned int win_len,
				  unsigned long nlocis);

/* Windows definition by size */
unsigned long Windowmaker_by_size (snp_t **snp_def, unsigned int winsize,
				   unsigned long nlocis, int parameter);

/* Windows filter by pos */
unsigned long Windowmaker_for_pos (snp_t **snp_def, snp_pos *pos_file,
				   unsigned int winsize, unsigned long nlocis,
				   int parameter);

// Apply the heuristic to determine the ancestral if it is unknown
int is_a2_ancestral (snp_t *snp, unsigned long nb_a2_holder, unsigned int nb_hap) ;
int is_a2_MAJ (unsigned long nb_a2_holder, unsigned int nb_hap) ;
int is_ancestral_defined (snp_t *snp ) ;

#endif
