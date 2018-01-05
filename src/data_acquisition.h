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

#ifndef __DATA_ACQUISITION_H_
#define __DATA_ACQUISITION_H_

#include "minilib.h"

// A sample cell
typedef struct sample {
  char *id;
  unsigned int pop_id;
  int sex;
  unsigned long rank; //in the phased file
  char *haplo1;
  char *haplo2;
}sample_t;

// All the samples and information about populations
typedef struct {
  unsigned int pop_number;
  unsigned int size;
  char **pops;
  unsigned long max_rank;
  unsigned long next_snp_of_sample_offset;
  sample_t** samples;
} samples;

// A snp cell
typedef struct snp {
  char *id;
  unsigned int pos;
  char a1;
  char a2;
  char ancestor;
  char _dummy_padding;

  unsigned int winlen_left; //Number of SNPs to the left of the core
  unsigned int winlen_right; //Number of SNPs to the right of the core
  list *win_starting_here; //a list of pair of long sorted by there second projection.
  // First long is the pos of the core of the window
  // Second long is the length of the window

} snp_t;

// Struct used when a file gives the snp to take into account
// Absent pos are still taken into account for stats of
// present pos the neighborough
typedef struct snp_data {
  char *id; //identifier of the SNP
  unsigned int pos; //SNP position
  int chr; //Number of the chromosome associated with the corresponding SNP position
} snp_pos;

typedef char DtString[BUFSIZ]; //Take maximum characters as possible, DtString is a table of char

#define A1 '0'
#define A2 '1'

//Free a sample cell
void free_sample (sample_t *s) ;
//Create a dummy cell
sample_t *dummy_sample(void) ;

//Free the samples
void free_samples (samples *s) ;

//Free a snp
void free_snpcol (unsigned long nlocis, snp_t **col);

// Read a legend file
conj *Legend_reader(char *prefix) ;

// Read and check a .hap file if transp or a .phased file else
void Init_Phased_reader(char *prefix, int transp, unsigned long snp_nb,
			samples *samples) ;

// Alloc and feed snp_res sample_col with
// prefix{.sample/.legend/(.hap if transp .phased else)} according to popfile
unsigned long Tripli_read(char *prefix, char *popfile_name,
			 snp_t ***snp_res, samples **sample_col, int transp);
#endif
