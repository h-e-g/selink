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

#include <config.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include "utils.h"
#include "stats.h"
#include "usage.h"
#include "printers.h"


void fprintf_if_number(FILE *f, const char *format, double x) {
  if (isnan(x) || isinf(x)) fprintf(f, "\t%s",NA);
  else fprintf(f, format, x);
}

void print_beginning_of_header(FILE *OUT) {
  fprintf(OUT, "#SNP\tstart\tstop\tcore_id\tcore_pos\tminor_allele\tancestral");
}

void print_output_header(FILE *OUT, char *pop_name, unsigned int pop_len,
			 int ocomp, char *ocustom, int snp_stats) {
  print_beginning_of_header(OUT);
  if (snp_stats & PRINT_K) fprintf(OUT, "\tK\tH");
  fprintf(OUT, "\tTheta_Pi");

  /* if (ocomp == 1) fprintf(OUT,"\tTheta_S\tTajima_D\tTajima_D(seqerr)\tFay&WuH\tFu&Li_D\tFu&Li_F"); **/
  if (ocomp == 1) fprintf(OUT,"\tTheta_S\tTajima_D\tFay&WuH");

  if (ocustom) fprintf(OUT, "\tTheta_custom\tNeutrality_custom1\tNeutrality_custom2");

  fprintf(OUT, "\tmaf\tdaf");
  if (snp_stats & PRINT_IHS) fprintf(OUT, "\tiHHa\tiHHd\tihs\tDeltat_iHH");
  if (snp_stats & PRINT_DIND) fprintf(OUT, "\tpiA\tpiD\tpiA/piD");
  fprintf(OUT, "\nPOPULATION: %s  (%d haplotypes)\n", pop_name, pop_len);
}

void print_excluded_header(FILE *OUT) {
  fprintf(OUT, "SNP_id\tPos_SNP\n");
}

void print_beginning_of_line(FILE *OUT, unsigned long core,unsigned int winlen_left,
			     unsigned int winlen_right, snp_t **snp,
			     int a2_is_ancestral, int a2_is_MAJ) {
  if ((winlen_left == UINT_MAX) || (winlen_right == UINT_MAX))
    fprintf(OUT, "%s\t--\t--",NA);
  else
    fprintf(OUT, "%d\t%u\t%u",
	    winlen_left+winlen_right + 1,snp[core-winlen_left]->pos,
	    snp[core+winlen_right]->pos);
  fprintf(OUT,"\t%s\t%u\t%c\t%c",
	  snp[core]->id, snp[core]->pos, a2_is_MAJ?snp[core]->a1:snp[core]->a2,
	  snp[core]->ancestor);
}

void print_one_snp_stats(FILE *OUT, FILE *EXCLUDED, snp_t **snp,
			 unsigned long core, int a2_is_ancestral, int a2_is_MAJ,
			 unsigned int pop_len, unsigned long *nb_a2_holder,
			 double *ihh_a1, double *ihh_a2, int out_of_window,
			 double *pi_a1, double *pi_a2) {
  double ihhA, ihhD, piA, piD;

  double maf = (a2_is_MAJ
		?pop_len - nb_a2_holder[core]
		:nb_a2_holder[core])/(double)pop_len;

  fprintf (OUT, "\t%1.3lf", maf);




  double daf = (a2_is_ancestral
		?pop_len - nb_a2_holder[core]
		:nb_a2_holder[core])/(double)pop_len;
  if( is_ancestral_defined(snp[core]) ){
  	fprintf (OUT, "\t%1.3lf", daf);
  }else{  
  	fprintf (OUT, "\t-");
  }
  
  
  

  if ((daf == 0.0)||(daf == 1.0)){ // Etienne, Choumouss
    if (ihh_a1) fprintf(OUT, "\t%s\t%s\t%s\t%s",NA,NA,NA,NA);
    if (pi_a1) fprintf(OUT, "\t%s\t%s\t%s",NA,NA,NA);
    if (EXCLUDED)
      fprintf(EXCLUDED, "%s\t%u\n", snp[core]->id, snp[core]->pos);
  }
  else {
    if (ihh_a1&&ihh_a2) {
      ihhA = a2_is_ancestral?ihh_a2[core]:ihh_a1[core];
      ihhD = a2_is_ancestral?ihh_a1[core]:ihh_a2[core];
      fprintf_if_number(OUT,"\t%1.3lf",ihhA);
      fprintf_if_number(OUT,"\t%1.3lf",ihhD);
      fprintf_if_number(OUT,"\t%1.5lf",log(ihhA/ihhD));
      fprintf_if_number(OUT,"\t%1.5lf",ihhA-ihhD);
    }
    if (pi_a1&&pi_a2) {
      if (out_of_window) fprintf(OUT,"\t%s\t%s\t%s",NA,NA,NA);
      else {
	piA = a2_is_ancestral?pi_a2[core]:pi_a1[core];
	piD = a2_is_ancestral?pi_a1[core]:pi_a2[core];
	fprintf_if_number(OUT,"\t%1.5lf",piA);
	fprintf_if_number(OUT,"\t%1.5lf",piD);
	fprintf_if_number(OUT,"\t%1.5lf",piA/piD);
      }
    }
  }
}

void print_not_computed_line(FILE *OUT, int ocomp, char *ocustom, int snp_stats) {
  if (snp_stats & PRINT_K) fprintf(OUT, "\t%s\t%s",NA,NA);
  fprintf(OUT, "\t%s",NA);
  /** if (ocomp == 1) fprintf(OUT, "\t%s\t%s\t%s\t%s\t%s\t%s",NA,NA,NA,NA,NA,NA);*/
  if (ocomp == 1) fprintf(OUT, "\t%s\t%s\t%s",NA,NA,NA);
  if (ocustom) fprintf(OUT, "\t%s\t%s\t%s",NA,NA,NA);
}

void print_files_one_pop(int verbose,int ocomp, char *ocustom, int snp_stats,
			 FILE *OUT, FILE *EXCLUDED, char *pop_name,
			 unsigned int *xi, snp_t **snp, unsigned long nlocis,
			 unsigned int pop_len, unsigned long *nb_a2_holder,
			 unsigned long *diversities, unsigned long *haplo_sums,
			 double *ihh_a1, double *ihh_a2,
			 double *pi_a1, double *pi_a2) {
  unsigned long core, start, stop, max_xi=0;
  double TD, FWH, FuLiD, FuLiF, theta_s, theta2, H;
  double TD_seqerr;
  unsigned int winlen_left, winlen_right, window_len;
  int a2_is_ancestral, a2_is_MAJ;

  double *W_TD=NULL, *W_FWH=NULL;
  double alpha_TD=0.0, beta_TD=0.0, alpha_FWH=0.0, beta_FWH=0.0;

  if (ocomp == 1) {
    Neutrality_test_precomputation(pop_len,2,1,&W_TD,&alpha_TD,&beta_TD);
    Neutrality_test_precomputation(pop_len,2,7,&W_FWH,&alpha_FWH,&beta_FWH);
  }

  print_output_header(OUT,pop_name,pop_len,ocomp,ocustom,snp_stats);
  if (EXCLUDED) print_excluded_header(EXCLUDED);

  start = stop = UINT_MAX;
  for (core = 0; core<nlocis; core++) { /* loop on window */
    winlen_left = snp[core]->winlen_left;
    winlen_right = snp[core]->winlen_right;

    a2_is_ancestral = is_a2_ancestral(snp[core],nb_a2_holder[core],pop_len);
    a2_is_MAJ = is_a2_MAJ(nb_a2_holder[core],pop_len);

    print_beginning_of_line(OUT,core,winlen_left,winlen_right,
			    snp,a2_is_ancestral,a2_is_MAJ);

    // Show progress on console
    if (verbose && ((nlocis > 100) && (core % (nlocis/100) == 0)))
      fprintf (stderr,"\b\b\b%02lu%%",(100*core)/nlocis);

    if ((winlen_left == UINT_MAX) || (winlen_right == UINT_MAX)
	|| ((window_len = winlen_right + winlen_left + 1) <= 1)) {
      print_not_computed_line(OUT,ocomp,ocustom,snp_stats);
    } else {
      if (start == UINT_MAX && stop == UINT_MAX)
	max_xi = update_AFS(xi,nb_a2_holder,snp, pop_len, max_xi, 0, 0,
			    core - winlen_left - 1, core + winlen_right);
      else
	max_xi = update_AFS(xi,nb_a2_holder,snp, pop_len, max_xi, start,
			    core - winlen_left, stop, core + winlen_right);

      start = core - winlen_left;
      stop = core + winlen_right;

      if (snp_stats & PRINT_K) {
	/*compute haplotype diversity */
	H = (pop_len - (haplo_sums[core] / ((double) pop_len)))
	  / (pop_len - 1);
	fprintf(OUT, "\t%ld\t%1.3lf", diversities[core], H);
      }
      /*compute theta Pi */
      fprintf(OUT, "\t%1.3lf", Theta_compute(2, xi, pop_len));

      if (ocomp == 1) {
	theta_s = Theta_compute(1, xi, pop_len);
	theta2 = Theta2_fromS(pop_len,
			      Polymorphism_counter(nb_a2_holder, pop_len,
						   start, window_len)
			      );
	TD = Neutrality_test(xi, max_xi, W_TD, alpha_TD, beta_TD, theta_s, theta2);
	//TD_seqerr = Neutrality_test(xi, pop_len, 9, 8, theta_s, theta2);
	TD_seqerr = 0 ;
	FWH = Neutrality_test(xi, max_xi, W_FWH, alpha_FWH, beta_FWH, theta_s, theta2);
	//FuLiD = Neutrality_test(xi, pop_len, 1, 3, theta_s, theta2);
	FuLiD = 0 ;
	//FuLiF = Neutrality_test(xi, pop_len, 2, 3, theta_s, theta2);
	FuLiF = 0 ;
	
	/** fprintf(OUT, "\t%1.4lf\t%1.5lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf", theta_s, TD, TD_seqerr, FWH, FuLiD, FuLiF);*/
	fprintf(OUT, "\t%1.4lf\t%1.4lf\t%1.4lf", theta_s, TD, FWH );
      }
    }
    
    
    print_one_snp_stats(OUT, EXCLUDED, snp, core, a2_is_ancestral, a2_is_MAJ,
			pop_len, nb_a2_holder, ihh_a1, ihh_a2,
			((winlen_left == UINT_MAX) || (winlen_right ==UINT_MAX)),
			pi_a1, pi_a2);
    fprintf(OUT, "\n");
  }

  if (ocomp == 1) {
    free(W_TD);
    free(W_FWH);
  }

  // Show progress on console
  if (verbose) fprintf (stderr,"\b\b\ball.\n");
}

// --- interpop ---

static void print_one_interpop_header(FILE *OUT, unsigned int nb_pop,
				      char **poplist, char *prefix) {
  unsigned int i, j;
  for (i=0; i < nb_pop - 1; i++)
    for (j=i+1;j < nb_pop; j++)
      fprintf(OUT, "\t%s_%s_%s", prefix, poplist[i], poplist[j]);
}

void print_interpop_header(FILE *OUT, samples *smpls, int snp_stats) {
  print_beginning_of_header(OUT);
  /** guillaume: correct buggs for -i function alone */
  fprintf(OUT, "\tGlob_Fst");
  print_one_interpop_header(OUT, smpls->pop_number, smpls->pops, "FST");
  print_one_interpop_header(OUT, smpls->pop_number, smpls->pops, "DeltaDAF");
  if (snp_stats & PRINT_IHS) {
    print_one_interpop_header(OUT, smpls->pop_number, smpls->pops, "Xpehh");
    print_one_interpop_header(OUT, smpls->pop_number, smpls->pops, "DeltaiHH");
  }
  fprintf(OUT, "\n");
}

static void print_one_not_computed_interpop(FILE *OUT, unsigned int nb_pop) {
  unsigned int i, j;
  for (i=0; i < nb_pop - 1; i++)
    for (j=i+1;j < nb_pop; j++)
      fprintf(OUT, "\t%s",NA);
}

void print_FST_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
		       unsigned int *pop_sizes, unsigned long **nb_a2_holder,
		       unsigned int nb_pop) {
  unsigned int i, j, n1, n2;
  unsigned long c1, c2;
  double ssd_wp1, ssd_wp2, ssd_t, nn, N, D, fst;

  for (i = 0; i < nb_pop - 1; i++) {
    c1 = a2_is_ancestral ? 2*pop_sizes[i] - nb_a2_holder[i][core] : nb_a2_holder[i][core];
    n1 = pop_sizes[i]*2;
    for (j = i+1; j < nb_pop; j++){
      c2 = a2_is_ancestral ? 2*pop_sizes[j] - nb_a2_holder[j][core] : nb_a2_holder[j][core];
      n2 = pop_sizes[j]*2;
      ssd_wp1 = SSD_WP(c1, n1);
      ssd_wp2 = SSD_WP(c2, n2);
      ssd_t = SSD_T(c1, n1, c2, n2);
      nn = (n1 + n2) - (((n1 * n1) + (n2 * n2)) / ((double) n1 + n2 ));
      N = ((n1 + n2 - 2) * ssd_t) - ((n1 + n2 - 1) * (ssd_wp1 + ssd_wp2));
      D = ((n1 + n2 -2) * ssd_t) - ((n1 + n2 - 1 -nn) * (ssd_wp1 + ssd_wp2));
      fst = N / D;

#ifdef VERBOSEDEBUG
      fprintf(stderr, ">> %i\t%d\t%1.3f\t%i\t%d\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n", c1[n], n1, ssd_wp1, c2[n], n2, ssd_wp2, ssd_t, nn, fst);
#endif

      fprintf_if_number(OUT,"\t%1.3lf",fst);
    }
  }
}
/** guillaume: global Fst */
void print_GLOB_FST_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
		       unsigned int *pop_sizes, unsigned long **nb_a2_holder,
		       unsigned int nb_pop) {
  unsigned int i=0, c=0, n=0, n_2=0, n_i=0;
  unsigned long c_i;
  double ssd_wp=0.0, ssd_t=0.0, n_c=0.0, b=0.0, N=0.0, D=0.0, fst=0.0;

  for (i = 0; i < nb_pop; i++) {
    c_i = a2_is_ancestral ? 2*pop_sizes[i] - nb_a2_holder[i][core] : nb_a2_holder[i][core];
    n_i = pop_sizes[i]*2;

    c += c_i; n += n_i; n_2 += n_i*n_i;

    ssd_wp+= (c_i * (n_i - c_i)) / ((double) n_i) ;
  }
  ssd_t=(c * (n - c)) / ((double) n);

  n_c=( n  -  ( n_2 / ((double) n) ) ) / ((double) nb_pop - 1) ;
  b=n_c*(nb_pop-1);

  N = ( ( n - nb_pop )*ssd_t  - ( n - 1 )*ssd_wp  );
  D = ( ( n - nb_pop )*ssd_t  - ( n - 1 - b )*ssd_wp );
  fst = N / D;

  fprintf_if_number(OUT,"\t%1.3lf",fst);
}

void print_Xpehh_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
			 double **ihh_a1, double **ihh_a2, unsigned int nb_pop) {
  unsigned int i, j;
  for (i=0; i < nb_pop - 1; i++)
    for (j = i + 1; j < nb_pop; j++)
      fprintf_if_number(OUT, "\t%1.4lf", log(a2_is_ancestral?ihh_a1[i][core] / ihh_a1[j][core]:ihh_a2[i][core] / ihh_a2[j][core]));
}

void print_Deltaihh_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
			    double **ihh_a1, double **ihh_a2, unsigned int nb_pop) {
  unsigned int i, j;
  for (i=0; i < nb_pop - 1; i++)
    for (j = i + 1; j < nb_pop; j++)
      fprintf(OUT, "\t%1.4lf", a2_is_ancestral?ihh_a1[i][core] - ihh_a1[j][core]:ihh_a2[i][core] - ihh_a2[j][core]);
}

void print_DeltaDAF_one_pos(FILE *OUT, unsigned long core, int a2_is_ancestral,
			    unsigned int *pop_sizes, unsigned long **nb_a2_holder,
			    unsigned int nb_pop) {
  unsigned int i, j;
  double p_i, p_j;
  for (i=0; i < nb_pop - 1; i++)
    for (j = i + 1; j < nb_pop; j++){
      p_i=nb_a2_holder[i][core] / ((double) 2*pop_sizes[i]);
      p_j=nb_a2_holder[j][core] / ((double) 2*pop_sizes[j]);
      fprintf(OUT, "\t%1.4lf", a2_is_ancestral?p_i - p_j:p_j - p_i);
    }
}

void print_interpop_file(int verbose, char *outfile, samples *smpls,
			 snp_t **snp, unsigned long nlocis,
			 unsigned int* pop_sizes, double **ihh_a1, double **ihh_a2,
			 unsigned long **nb_a2_holder, int snp_stats) {
  FILE *XPOUT = stdout;
  char *xpoutfile=NULL;
  int a2_is_ancestral, a2_is_MAJ;
  unsigned long core;
  unsigned int i, winlen_left, winlen_right, window_len;
  unsigned long total_nb_a2;

  if (outfile) {
    // generate xpehh outfile
    xpoutfile = append_extension(outfile, "interpop");
    XPOUT = safe_fopen(xpoutfile, "w");
  }

  if (verbose)
    fprintf(stderr, "computes and outputs results for interpop: 00%%");

  print_interpop_header(XPOUT, smpls,snp_stats);
  for (core = 0; core<nlocis; core++) { /* loop on window */

    // Show progress on console
    if (verbose && ((nlocis > 100) && (core % (nlocis/100) == 0)))
      fprintf (stderr,"\b\b\b%02lu%%",(100*core)/nlocis);

    winlen_left = snp[core]->winlen_left;
    winlen_right = snp[core]->winlen_right;

    total_nb_a2 = 0;
    for (i=0;i<smpls->pop_number;i++) total_nb_a2 += nb_a2_holder[i][core];
    a2_is_ancestral = is_a2_ancestral(snp[core],total_nb_a2,smpls->size*2);
    a2_is_MAJ = is_a2_MAJ(total_nb_a2,smpls->size*2);

    print_beginning_of_line(XPOUT,core,winlen_left,winlen_right,
			    snp,a2_is_ancestral,a2_is_MAJ);

    print_GLOB_FST_one_pos(XPOUT,core,a2_is_ancestral,
			   pop_sizes,nb_a2_holder,smpls->pop_number);
    print_FST_one_pos(XPOUT,core,a2_is_ancestral,
		      pop_sizes,nb_a2_holder,smpls->pop_number);
    print_DeltaDAF_one_pos(XPOUT,core,a2_is_ancestral,
			   pop_sizes,nb_a2_holder,smpls->pop_number);
    if (snp_stats & PRINT_IHS) {
      if ((winlen_left == UINT_MAX) || (winlen_right == UINT_MAX) || ((window_len = winlen_right + winlen_left + 1) < 1)) {
	print_one_not_computed_interpop(XPOUT, smpls->pop_number);
	print_one_not_computed_interpop(XPOUT, smpls->pop_number);
      } else {
	print_Xpehh_one_pos(XPOUT,core,a2_is_ancestral,
			    ihh_a1,ihh_a2,smpls->pop_number);
	/** guillaume Delta ihh */
	print_Deltaihh_one_pos(XPOUT,core,a2_is_ancestral,
			       ihh_a1,ihh_a2,smpls->pop_number);
      }
    }
    fprintf(XPOUT,"\n");
  }
  // Show progress on console
  if (verbose) fprintf (stderr,"\b\b\ball.\n");

  if (outfile) {
    safe_fclose(xpoutfile,XPOUT);
    free(xpoutfile);
  }
}
