/* Copyright 2012 GEH Institut Pasteur */
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

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/param.h>
#include <stdbool.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <err.h>

#if HAVE_PTHREAD_H
#include <pthread.h>
#endif

#include "usage.h"
#include "utils.h"
#include "stats.h"
#include "primitive_stats.h"
#include "arp-read.h"
#include "printers.h"

#define EHHTHRESHOLD 0.05

int Go_compute (char *prefix, char *outfile, int arp_read, unsigned int winlen, unsigned int winsize, char *popfile_name, int ocomp, char *ocustom, int snp_stats, int interpop, int best, snp_pos *posfile, int transp, int verbose, unsigned long threads);

int main (int argc, char **argv) {
  
  /* variables and initialisation */
  
  //char *infile_name,  *outfile, *popfile_name, *omega_custom; 
  char *infile_name,  *outfile, *popfile_name, *omega_custom; //Updated by Choumouss
  char *prog; 
  int i, arp_read, verbose;
  int interpop, snp_stats, omega_comp;
  unsigned int winlen, winsize;
  //char **poplist;
  int best, hap_transp=0; //Choumouss
  char *posfile_name; //Choumouss
  int chr; //Choumouss
  unsigned long threads;
  snp_pos *posfile_data; //Choumouss
  

  /* default values */
  outfile = NULL;
  popfile_name = NULL;
  omega_custom = NULL;
  winlen = 0;
  arp_read = 0;
  omega_comp = 0;
  verbose = true;
  snp_stats = 0; //2^0 means ihs 2^1 means piA/piD 2^2 means K
  interpop = 0;
  /** best = 0; //Choumouss */
  /** guillaume; pour virer l'option -b et mettre best par défaut */
  best = 1; 
  posfile_name = NULL; //Choumouss
  chr = 0; //Choumouss
  posfile_data = NULL; //Choumouss
  //window_param = 0; //Choumouss
  winsize = 0; //Choumouss
  threads = 1;

  /* get progname */
  prog = basename(argv[0]);
  
  /* check syntax option on command line */
  i = 0;
  while((i = getopt_long (argc, argv, "acf:hqil:n:bP:C:Tt:o:spKvw", long_options, NULL)) != -1) {
	 switch(i) {

	 case 'a':
		arp_read = 1;
		break;

	 case 'c':
		omega_custom = optarg;
		break;

	 case 'h':
		Usage(prog, EXIT_SUCCESS);

	 case 'q':
	   verbose=false;
	   break;

	 case 'i':
		interpop = 1;
		break;
		
	 /** guillaume: for Fst only computed on every snp
	 case 'I':
		interpop_fst = 1;
		break;
         */ 
	 
	 case 'l':
	   winsize = ATOUI(optarg);
	   break;

	 case 'n':
	   winlen = ATOUI(optarg);
	   break;

	 case 'o':
		outfile = optarg;
		break;

	 case 'f':
		popfile_name = optarg;
		break;

	 case 's':
		snp_stats |= PRINT_IHS;
		break;
	 case 'p':
		snp_stats |= PRINT_DIND;
		break;
	 case 'K':
		snp_stats |= PRINT_K;
		break;

	 case 'v':
      (void)fprintf(stdout, "%s version %s %s\n", prog, PACKAGE, VERSION);
      return EXIT_SUCCESS;

	 case 'w':	  		
		omega_comp = 1;
		break;

	 /* Choumouss */
	 case 'b': //Best SNP ends the windows
	 	/** best = 1; */
  		/** guillaume; pour virer l'option -b et mettre best par défaut */ 
	 	best = 0;
	 	break;

	 case 'P': //SNP position filter
		posfile_name = optarg;
		break;

	 case 'C':
	 	chr = atoi(optarg); //Chromosome number
	 	break; 

	 case 't':
	   threads = strtoul(optarg,NULL,10);
	   break;

	 case 'T':
	   hap_transp = 1;
	   break;

	 default :
      Usage(prog, EXIT_FAILURE);
    }
  }

 /* check for mandatory parameters */
  if (argc - optind != 1) { Usage(prog, EXIT_FAILURE); }
  infile_name = argv[optind];

  if (posfile_name){
  	posfile_data = Positions_reader(posfile_name, chr);
  }

  Go_compute(infile_name, outfile, arp_read, winlen, winsize, popfile_name,
	     omega_comp, omega_custom, snp_stats, interpop, best,
	     posfile_data, hap_transp, verbose, threads);

  return EXIT_SUCCESS;
}

#if HAVE_PTHREAD_H
static pthread_mutex_t main_is_working = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t a_thread_has_finished = PTHREAD_COND_INITIALIZER;
#endif

typedef struct {
  int verbose;
  int ocomp;
  char *ocustom;
  int snp_stats;
  int _dummy_pad;
  char *outfile;
  unsigned long nlocis;
  unsigned int *pop_sizes;
  samples *smpls;
  snp_t **snp;
  double **ihh_a1;
  double **ihh_a2;
  double **pi_a1;
  double **pi_a2;
  unsigned long **nb_a2_holder;
  unsigned long **diversities;
  unsigned long **haplo_sums;
} all_datas;

static void* one_pop_computations(void *stuff) {
  conj *p = (conj*) stuff;
  all_datas *data = (all_datas*) proj1(p);
  unsigned long pop_id = (unsigned long) proj2(p);

  unsigned int pop_len;
  char *pop_name;
  unsigned int *xi; // allele frequency

  pop_name = data->smpls->pops[pop_id];
  pop_len = 2*data->pop_sizes[pop_id];

  fill_half_ihh_pi_divs_and_nbd_arrays_one_pop
    (data->verbose,EHHTHRESHOLD,false,pop_id,pop_len,data->smpls,data->snp,
     data->nlocis,MAYBE(data->ihh_a1,pop_id),MAYBE(data->ihh_a2,pop_id),NULL,
     MAYBE(data->pi_a1,pop_id),MAYBE(data->pi_a2,pop_id),
     MAYBE(data->diversities,pop_id),MAYBE(data->haplo_sums,pop_id));

  fill_half_ihh_pi_divs_and_nbd_arrays_one_pop
    (data->verbose,EHHTHRESHOLD,true,pop_id,pop_len,data->smpls,data->snp,
     data->nlocis,MAYBE(data->ihh_a1,pop_id),MAYBE(data->ihh_a2,pop_id),
     MAYBE(data->nb_a2_holder,pop_id), MAYBE(data->pi_a1,pop_id),
     MAYBE(data->pi_a2,pop_id),NULL,NULL);

  if (!(xi = calloc(pop_len+1,sizeof(unsigned int))))
    err(EXIT_FAILURE, "Unable to allocate array of allele_frequencies");

  // An output file for each population
  FILE *OUT= stdout;
  FILE *EXCLUDED = NULL;
  char name_file[FILENAME_MAX];
  char name_excludedFile[FILENAME_MAX];
  if (data->outfile != NULL){
    //adds the population name to the output file prefix
    sprintf(name_file, "%s_%s.out", data->outfile, pop_name);
    //adds the population name to the output file prefix
    sprintf(name_excludedFile, "%s_%s_excluded.out", data->outfile, pop_name);
    OUT = safe_fopen(name_file, "w");
    EXCLUDED = safe_fopen(name_excludedFile, "w");
  }

  if (data->verbose)
    fprintf(stderr, "computes and outputs results for pop %lu/%d: 00%%",
	    pop_id+1,data->smpls->pop_number);
  print_files_one_pop(data->verbose, data->ocomp, data->ocustom, data->snp_stats,
		      OUT, EXCLUDED, pop_name, xi, data->snp, data->nlocis, pop_len,
		      data->nb_a2_holder[pop_id], MAYBE(data->diversities,pop_id),
		      MAYBE(data->haplo_sums,pop_id), MAYBE(data->ihh_a1,pop_id),
		      MAYBE(data->ihh_a2,pop_id), MAYBE(data->pi_a1,pop_id),
		      MAYBE(data->pi_a2,pop_id));
  free(xi);
  if(OUT != stdout) safe_fclose(name_file,OUT);
  if (EXCLUDED) safe_fclose(name_excludedFile,EXCLUDED);
#if HAVE_PTHREAD_H
  pthread_cond_signal(&a_thread_has_finished);
#endif
  return NULL;
}

int Go_compute (char *prefix, char *outfile, int arp_read, unsigned int winlen, unsigned int winsize, char *popfile_name, int ocomp, char *ocustom, int snp_stats, int interpop, int best, snp_pos *posfile, int transp, int verbose, unsigned long threads){

  unsigned int i;
  unsigned long pop_id;
  all_datas data = { verbose, ocomp, ocustom, snp_stats, 0, outfile, 0, NULL, NULL,
		     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
  void *dummy;


	
  /**  debug */
  /** fprintf(stderr, "yo 1\n"); */

  if (verbose) fprintf(stderr, "read and check input files\n");
  if (arp_read) err(EXIT_FAILURE, "arp reader is unplugged, Sorry!\n");
  //sample_nb = Arp_reader(prefix, &snp, &samples, &ancestor);
  else {
    data.nlocis =
      Tripli_read(prefix, popfile_name, &(data.snp), &(data.smpls), transp);
  }

  if (verbose) fprintf(stderr, "got %u samples description\n", data.smpls->size);
  if (verbose) fprintf(stderr, "read %lu snp descriptions\n", data.nlocis);

  if (winsize){
    if (posfile) //if positions filter is requested
      Windowmaker_for_pos (data.snp, posfile, winsize, data.nlocis, best);
    else Windowmaker_by_size(data.snp, winsize, data.nlocis, best);
  }
  else {
    if (winlen) {
      if (winlen > data.nlocis)
	err(EXIT_FAILURE,
	    "Incompatible values win len (%u) vs nlocis (%lu)",winlen,data.nlocis);
      Windowmaker_by_snp(data.snp, winlen, data.nlocis);
    } else Windowmaker_by_snp(data.snp, MIN(data.nlocis,DEFAULT_WINLEN), data.nlocis);
  }

  if ((data.pop_sizes = calloc(data.smpls->pop_number,sizeof(int))) == NULL)
    err(EXIT_FAILURE, "malloc failed");
  for (i=0; i<data.smpls->size; i++) data.pop_sizes[data.smpls->samples[i]->pop_id]++;

  data.nb_a2_holder = create_long_array_for_results(data.smpls,data.nlocis);
  if (snp_stats & PRINT_K) {
    data.diversities = create_long_array_for_results(data.smpls,data.nlocis);
    data.haplo_sums = create_long_array_for_results(data.smpls,data.nlocis);
  }
  if (snp_stats & PRINT_DIND) {
    data.pi_a1 = create_double_array_for_results(data.smpls,data.nlocis);
    data.pi_a2 = create_double_array_for_results(data.smpls,data.nlocis);
  }
  if (snp_stats & PRINT_IHS) {
    data.ihh_a1 = create_double_array_for_results(data.smpls,data.nlocis);
    data.ihh_a2 = create_double_array_for_results(data.smpls,data.nlocis);
  }

  /* cache harmonic sum used in various stats */
  HARMONICSUMS_CACHED = Cache_HarmonicSums(data.nlocis);

#if HAVE_PTHREAD_H
  if (threads > 1) {
    pthread_t *worker;
    if ((worker = calloc(data.smpls->pop_number,sizeof(pthread_t))) == NULL)
      err(EXIT_FAILURE, "Unable to allocate array to store thrads.\n");;
    conj *p;
    if ((p = calloc(data.smpls->pop_number,sizeof(conj))) == NULL)
      err(EXIT_FAILURE, "Unable to allocate array to store thrads.\n");;

    for(pop_id = 0;pop_id<data.smpls->pop_number;pop_id++) {
      p[pop_id].fst = &data;
      p[pop_id].snd = (void *) pop_id;
    }

    for(pop_id = 0;pop_id<MIN(data.smpls->pop_number,threads);pop_id++)
      pthread_create(worker + pop_id,NULL,one_pop_computations,p+pop_id);
    for(;pop_id<data.smpls->pop_number;pop_id++) {
      pthread_cond_wait(&a_thread_has_finished,&main_is_working);
      pthread_create(worker + pop_id,NULL,one_pop_computations,p+pop_id);
    }
    for(pop_id = 0;pop_id<data.smpls->pop_number;pop_id++)
      pthread_join(worker[pop_id],&dummy);
    free(p);
    free(worker);
  }
  else {
#endif

    /** debug */
    /** fprintf(stderr, "yo 2\n"); */
    
    
    for(pop_id = 0;pop_id<data.smpls->pop_number;pop_id++) {
      conj p = { &data, (void *) pop_id};
      dummy = one_pop_computations(&p);
    }
#if HAVE_PTHREAD_H
  }
#endif

  /** debug */
  /** fprintf(stderr, "yo 3\n"); */
  if (interpop == 1)
    print_interpop_file(verbose, outfile, data.smpls, data.snp, data.nlocis,
			data.pop_sizes, data.ihh_a1, data.ihh_a2,
			data.nb_a2_holder, snp_stats);

#ifndef EFBUG
  free(data.pop_sizes);
  free_double_results_array(data.smpls,data.ihh_a1);
  free_double_results_array(data.smpls,data.ihh_a2);
  free_long_results_array(data.smpls,data.nb_a2_holder);
  free_long_results_array(data.smpls,data.diversities);
  free_long_results_array(data.smpls,data.haplo_sums);
  free_double_results_array(data.smpls,data.pi_a1);
  free_double_results_array(data.smpls,data.pi_a2);

  free_samples(data.smpls);

  free_snpcol(data.nlocis,data.snp);

  free(HARMONICSUMS_CACHED);
  free(posfile);
 #endif

  return 0;
}
