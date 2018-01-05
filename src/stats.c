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

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include "data_acquisition.h"
#include "utils.h"
#include "stats.h"



#define BUFFLEN 100

/* compute harmonic sum see Fu equation (4) */
static inline double a_FU(unsigned int n);
/* compute beta n : see FU paper equation (5) */
static double Beta_FU(unsigned int i, unsigned int n);
/* compute sigma_ii see FU paper equation (2) */
static double Sigma_ii(unsigned int i, unsigned int n);
/* compute sigma_ij see FU paper equation (3) */
static double Sigma_ij(unsigned int i, unsigned int j, unsigned int n);


double SSD_WP (unsigned long c, unsigned int n) {
  return (c * (n - c)) / ((double) n);
}

double SSD_T(unsigned long c1, unsigned int n1,
	     unsigned long c2, unsigned int n2) {
  return ((c1 + c2) * ((n1 - c1) +  (n2 - c2))) / ((double) n1 + n2);
}



/* cache harmonic sums from 0 -> pop size */
/* avoids to recompute HS for every value */
double *Cache_HarmonicSums(unsigned long n){
  double *hs;
  unsigned long i;
  // prepare enough space to have HS values for n and n+1
  hs = malloc(sizeof(double) * (n+1+1));
  if (!hs) err(EXIT_FAILURE, "malloc failed");
  
  hs[0] = 0;
  for(i = 1; i <= n; i++) {
	 hs[i] = hs[i-1] + 1.0/i;
  }
  return hs;
}

/* returns harmonic sum for given value. */
static inline double a_FU( unsigned int n ){
  return HARMONICSUMS_CACHED[n-1];
}

static inline double Beta_FU(unsigned int i, unsigned int n) {
  return ( ((2.0 * n) / ((n-i+1)*(n-i))) * (a_FU(n+1) - a_FU(i)) )
    - (2.0/(n -i));
}
static double Sigma_ii(unsigned int i, unsigned int n) {
  double s_ii;
  /* n may be odd so avoid round problems raised by entire division  */
  if (i*2 < n) // i < n/2
    s_ii = Beta_FU(i+1, n);
  else if (2*i == n) //i == n/2
    s_ii =(2.0 * (a_FU(n) - a_FU(i)) / (n-i)) - ( 1.0/(i*i) );
  else // i > n/2
    s_ii = Beta_FU(i, n) - ( 1.0/(i*i) );
  return s_ii;
}

static double Sigma_ij(unsigned int j, unsigned int i, unsigned int n) {
  double s_ij;

  if (i+j < n)
    s_ij =(Beta_FU(i+1, n) - Beta_FU(i, n)) / 2.0;
  else if (i+j == n) {
    s_ij  = (a_FU(n) - a_FU(i)) / (n - i);
    s_ij += (a_FU(n) - a_FU(j)) / (n - j);
    s_ij -= ( Beta_FU(i, n) + Beta_FU(j+1, n) ) / 2.0;
    s_ij -= 1.0 / (i*j);
  }
  else { //i + j > n
    s_ij  =  (Beta_FU(j, n) - Beta_FU(j+1, n)) / 2.0;
    s_ij -=  1.0 / (i*j);
  }
  return s_ij;
}

double omega_from_type(int type, unsigned int popsize, unsigned int i) {
  switch (type) {
  case 1: // theta S
  case 4: // theta S-etha1
    return 1.0/i;
  case 2: // theta PI
  case 5: // theta PI-etha1
    return popsize - i;
  case 3: // theta etha1 since we use DAFS raher than AFS, singletons are mutation at i=1 and i=n-1
    if(i == 1){
    	return 1.0;
    }else if ( i == popsize - 1 ){
    	return 1.0;
    }else{
    	return 0.0;
    }
  case 6: // theta xi
    return 0.0;
  case 7: //  theta etha  == theta L
    return i;
  case 8: // theta S-xi
    if( i == 1 ){
    	return 0.0;
    }else if ( i == popsize - 1 ){
    	return 0.0;
    }else{
    	return 1.0/i;
    }
  case 9: // theta PI-xi
    if(i == 1){
    	return 0.0;
    }else if ( i == popsize - 1 ){
    	return 0.0;
    }else{
    	return popsize - i;
    }
  default :
    err(EXIT_FAILURE, "%d: unknow omega estimator", type);
  }
}

double sum_w_from_type(int type, unsigned int popsize) {
  switch (type) {
  case 1: // theta S
  case 4: // theta S-etha1
    return HARMONICSUMS_CACHED[popsize-1];
  case 2: // theta PI
  case 5: // theta PI-etha1
  case 7: //  theta etha  == theta L
    return (popsize * (popsize - 1)) / 2.0;
  case 3: // theta etha1 since we use DAFS raher than AFS, singletons are mutation at i=1 and i=n-1
    /** return 1.0; */
    return 2.0;
  case 6: // theta xi
    return 0.0;
  case 8: // theta S-xi
    return ( HARMONICSUMS_CACHED[popsize-2] - 1.0 );
  case 9: // theta PI-xi
    return  ( (popsize * (popsize - 1)) / 2.0  - (popsize - 1) - 1.0 );
  default :
    err(EXIT_FAILURE, "%d: unknow sum omega estimator", type);
  }
}

/* compute equation (9) from Achaz paper */
/** guillaume: modif to generalise this function*/
/* double Neutrality_test( unsigned int *xi, unsigned int n, double theta, double theta2) { */
void Neutrality_test_precomputation (unsigned int n, int w1,  int w2,
				     double **W, double *alpha, double *beta) {
  double tmp, Sum_w1, Sum_w2;
  unsigned int i, j;

  if ((*W = calloc(n,sizeof(double))) == NULL)
    err(EXIT_FAILURE, "malloc failed");

  Sum_w1 = sum_w_from_type(w1,n);
  Sum_w2 = sum_w_from_type(w2,n);

  for(i=1; i<n; i++)
    (*W)[i] = (omega_from_type(w1,n,i)/Sum_w1) - (omega_from_type(w2,n,i)/Sum_w2);

  *alpha=*beta=0.0;
  for(i=1; i < n ; i++){
    tmp = (*W)[i] * (*W)[i] * i;
    *alpha += tmp;
    *beta  +=  i*tmp * Sigma_ii(i, n);
    for( j=i+1; j< n ; j++){
      *beta  +=  2*i*j*(*W)[i]*(*W)[j] * Sigma_ij(i, j, n);
    }
  }
}

double Neutrality_test(unsigned int *xi, unsigned long n, double *W, double alpha,
		       double beta, double theta, double theta2) {
  double T_OMEGA = 0.0;
  unsigned long i;

  for(i=1; i <= n ; i++)
    T_OMEGA += (double) i * xi[i] * W[i];
  T_OMEGA /= sqrt( alpha*theta  + beta*theta2);

  return T_OMEGA;
}

/*
 * compute polymorphism on the given window (number of column with polymorphism)
 */
unsigned int Polymorphism_counter(unsigned long *nb_a2, unsigned int pop_len,
				  unsigned long start, unsigned int len) {
  unsigned int i, n = 0;

  for (i=0; i<len;i++)
    if ((nb_a2[start+i] != 0) && (nb_a2[start+i] != pop_len))
      n++;
  return n;
}

unsigned long update_AFS(unsigned int *xi, unsigned long *nb_a2_holder, snp_t **snp,
			 unsigned int pop_len, unsigned long max_xi,
			 unsigned long old_start, unsigned long new_start,
			 unsigned long old_stop, unsigned long new_stop) {
  unsigned long i, k;
  for (i=old_start; i < new_start; i++) {
    k = is_a2_ancestral(snp[i],nb_a2_holder[i],pop_len)
      ?pop_len - nb_a2_holder[i]
      :nb_a2_holder[i];
    if (k != 0) xi[k]--;
  }
  for (i=old_stop+1; i <= new_stop; i++) {
    k = is_a2_ancestral(snp[i],nb_a2_holder[i],pop_len)
      ?pop_len - nb_a2_holder[i]
      :nb_a2_holder[i];
    if (k != 0) {
      xi[k]++;
      max_xi = max_xi<k?k:max_xi;
    }
  }
  return max_xi;
}

double Theta_compute(int type, unsigned int *xi, unsigned int popsize){
  double sum_wixi = 0.0;
  unsigned int i;

  for (i = 1; i < popsize; i++) {
    sum_wixi += omega_from_type(type,popsize,i) * i * xi[i];
  }
  return sum_wixi / sum_w_from_type(type,popsize);
}

double Theta2_fromS(unsigned int popsize, unsigned int snpnb) {
  double an, bn;
  double S;
  double theta2;
  unsigned int i;

  S = (double) snpnb;
  an = HARMONICSUMS_CACHED[popsize-1];

  bn = 0.0;
  for (i = 1; i < popsize; i++) {
	 bn += 1.0/(i*i);
  }

  theta2 =( S * (S - 1.0 )) / (an*an + bn);

#ifdef VERBOSEDEBUG
  fprintf (stderr, "DEBUG: theta2_fromS: (posize = %d, snp nb = %d)\n", popsize, snpnb);
  fprintf(stderr, ">> S\tS-1\tan\tbn\t(S * S-1)/ (an*an + bn)\n");
  fprintf(stderr, ">> %1.1f\t%1.1f\t%1.3f\t%1.3f\t%1.3f", S, S-1, an, bn, theta2);
  fprintf(stderr, "\n\n");
#endif

  return theta2;
}

/* TODO : check that sum of values == 1 */
double *Omega_reader(char *filename, unsigned int pop_size) {

  unsigned int l, i;
  double *omega;
  FILE *IN;
  char *BUFF, *p, *q;
  size_t L;

  if((IN = fopen(filename, "r")) == NULL) { err(EXIT_FAILURE, "%s: open failed", filename);
  }

  if ((BUFF  = (char *)malloc(sizeof(char)*(BUFFLEN))) == NULL) err(EXIT_FAILURE, "malloc failed");

  if ((omega = malloc(pop_size*sizeof(double))) == NULL) err(EXIT_FAILURE, "malloc failed");

  L = BUFFLEN;
  p = BUFF;
  while (fgets(p, BUFFLEN, IN) != NULL) {
    /*check if all datas read */
    if ((strchr(p, '\n')) == NULL) { //realloc
      L += BUFFLEN;
      if ((BUFF = (char *)realloc(BUFF, L * sizeof(char))) == NULL) err(EXIT_FAILURE, "malloc failed");
      p = BUFF + strlen(BUFF);
      continue;
    }
  }
  
  /* split buffer */
  q = BUFF;
  L = strlen(BUFF);

  /* split buffer */
  while ((p = strsep(&q, " \t\n")) != NULL) { continue; }

  i = 0, l = 0;
  p = BUFF;
  while (l<L) {
    if (isdigit(*p)) {
      omega[i] = strtod(p,NULL);
      printf(">> %s %lf\n", p, omega[i]);
      l += strlen(p);
      p += strlen(p);
      i++;
      continue;
    }
    l++;
    p++;
  }
  if (i != pop_size) err(EXIT_FAILURE, "inconsistent values");

  if (fclose(IN) == EOF) err(EXIT_FAILURE, "%s: close failed", filename);
#ifndef EFBUG
	free(BUFF);
#endif

  return omega;
}
