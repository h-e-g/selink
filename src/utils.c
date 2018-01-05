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
#include <err.h>
#include <sys/param.h>
#include <ctype.h>
#include <string.h>

#include "data_acquisition.h"
#include "utils.h"

#define TABLELEN 100

/* compare sequence with reference (ancestor) and toggles to */
/* binary representation 0 = similarity 1 = difference */
int Seq2bin(char *ref, char **seq){
  int l;
  char *s, *r;
  l = 0;
  s = *seq;
  r = ref;

  while(*r) {
	 if (*s == *r) { *s++ = '0'; }
	 else  { *s++ = '1'; }
	 r++;
	 l++;
  }
  *s = '\0';
   return l;
}

unsigned long Windowmaker_by_snp (snp_t **snp,
				  unsigned int win_len, unsigned long nlocis) {
  unsigned long i, winmax, dummy;

  if (win_len % 2 == 0) win_len++;
  dummy = win_len;
  winmax = nlocis - win_len/2  ;

  if (win_len <= nlocis)
    for (i=win_len/2; i < winmax; i++) {
      snp[i]->winlen_left = win_len/2; //Number of SNPs to the left of the core
      snp[i]->winlen_right = win_len/2; //Number of SNPs to the right of the core
      snp[i-win_len/2]->win_starting_here =
	cons(pair((void *)i, (void *) dummy),
	     snp[i-win_len/2]->win_starting_here);
    }
  return win_len;
}

snp_pos *Positions_reader (char *pos_file, int num_chr){
  FILE *fichier_pos;
  char *p, *q;
  char *tmp[4], **t;
  snp_pos *pos_data, *dt;
  unsigned int n, m;

  if ((pos_data=(snp_pos *)malloc(TABLELEN*sizeof(snp_pos))) == NULL)
    err(EXIT_FAILURE, "malloc failed");
  fichier_pos = fopen(pos_file, "r");
  DtString ligne;
  dt = pos_data;
  n = 0;
  m = TABLELEN;
  while (fgets(ligne, TABLELEN, fichier_pos) != NULL){
  	q = ligne;
  	t = tmp;
	while ((p = strsep (&q, "\t")) != NULL){
		*t++= p;
	}
	if (atoi(tmp[2]) == num_chr){ //verify the chromosome name
 		dt[n].id = strdup(tmp[0]); //identifier
 		dt[n].pos = ATOUI(tmp[1]); //SNP position
 		dt[n].chr = atoi(tmp[2]); //Chromosome
 		n++;
 	}
	if (n >= m){
		m += TABLELEN;
		if (!(pos_data =(snp_pos *) realloc(pos_data, m*sizeof(snp_pos))))
		  err(EXIT_FAILURE, "realloc failed");
		dt = pos_data+n;
	}
  }
  dt[n].id = NULL;
  dt[n].pos = UINT_MAX;
  dt[n].chr = -1;

  return dt;
}

unsigned long Best_snp (snp_t **snp_def, unsigned long snp_ref_pos,
			unsigned long num_bef, unsigned long num_aft,
			unsigned int winsize) {
  unsigned long pos_target, num;
  snp_t *snp1_prox, *snp2_prox;
  snp1_prox = snp_def[num_bef];
  snp2_prox = snp_def[num_aft];
  pos_target = snp_ref_pos + (winsize/2);
  if ((snp2_prox->pos - pos_target) > (pos_target - snp1_prox->pos))
    num = num_bef;
  else num = num_aft;
  return num; //absolute position
}

// Determine last SNP (for stop): get the absolute position of the SNP
unsigned long Snp_search_right (snp_t **snp_def, unsigned long rank_snp_ref,
				unsigned int winsize, unsigned long nlocis,
				int use_best) {
  unsigned long num = ULONG_MAX;
  snp_t *snp, *snp_ref = snp_def[rank_snp_ref];
  unsigned long cpt = rank_snp_ref;
  snp = snp_def[cpt];
  while (((snp->pos - snp_ref->pos + 1) < winsize/2) && (cpt < nlocis - 1))
    snp = snp_def[++cpt];

  if ((snp->pos - snp_ref->pos + 1) == winsize/2)
    num = cpt;
  else
    if ((cpt != nlocis-1) || ((snp->pos - snp_ref->pos + 1) > winsize/2))
      { //return -1 if the border is reached
	if (use_best)
	  num = Best_snp (snp_def, snp_ref->pos, cpt-1, cpt, winsize);
	else num = cpt-1;
      }

  return num; //absolute position
}

// Determine last SNP (for start)
unsigned long Snp_search_left (snp_t **snp_def, unsigned long rank_snp_ref,
			       unsigned int winsize, int use_best) {
  unsigned long cpt, num=ULONG_MAX;
  snp_t *snp, *snp_ref = snp_def[rank_snp_ref];
  cpt = rank_snp_ref;
  snp = snp_def[cpt];

  while (((snp_ref->pos - snp->pos + 1) < winsize/2) && (cpt > 0))
    snp = snp_def[--cpt];

  if ((snp_ref->pos - snp->pos + 1) == winsize/2) num = cpt;
  else
    if ((cpt != 0) || ((snp_ref->pos - snp->pos + 1) > winsize/2))
      { //return -1 if the border is reached
	if (use_best)
	  num = Best_snp (snp_def, snp_ref->pos, cpt+1, cpt, winsize);
	else num = cpt+1;
      }

  return num;
}

unsigned long Windowmaker_for_pos (snp_t **snp_def, snp_pos *pos_file,
				   unsigned int winsize, unsigned long nlocis,
				   int use_best) {
  unsigned long i, l = 0;
  unsigned long start_num, stop_num;
  unsigned long pos, win_len, win_max_len=0;

  while(pos_file[l].pos != UINT_MAX){
    pos = pos_file[l].pos;
    i = 0;
    while (i < nlocis && snp_def[i]->pos != pos)
      i++;
    if (snp_def[i]->pos != pos) {
      printf("The position %lud in the position file does not exist as a SNP in this chromosome.\n", pos);
    } else {
      if ((i != 0) && (i == nlocis-1)
	  && ((start_num = Snp_search_left(snp_def, i, winsize, use_best))
	      != ULONG_MAX)
	  && ((stop_num = Snp_search_right(snp_def,i,winsize,nlocis,use_best))
	      != ULONG_MAX)) {
	win_len = stop_num - start_num + 1;
	win_max_len = MAX(win_len,win_max_len);

	//number of SNPs to the left of the core
	snp_def[i]->winlen_left = (unsigned int) (i - start_num);
	//number of SNPs to the right of the core
	snp_def[i]->winlen_right = (unsigned int) (stop_num - i);
	snp_def[start_num]->win_starting_here =
	  cons(pair((void *)i, (void *)win_len),
	       snp_def[start_num]->win_starting_here);
      }
    }
    l++;
  }
  // win_starting_here must be sorted in increasing length
  for (i=0;i<nlocis;i++)
    snp_def[i]->win_starting_here =
      destructive_rev_append(snp_def[i]->win_starting_here,nil());

  return win_max_len;
}


unsigned long Windowmaker_by_size (snp_t **snp_def, unsigned int winsize,
				   unsigned long nlocis, int use_best) {
  unsigned long start_num, stop_num;
  unsigned long i,win_len,win_max_len=0;

  // We start from the end because
  // win_starting_here must be sorted incrementally by length
  for ( i = nlocis - 2; i > 0; i-- ) {
    if ((i != 0) && (i != nlocis-1)
	&& ((start_num = Snp_search_left (snp_def, i, winsize, use_best))
	    != ULONG_MAX)
	&& ((stop_num = Snp_search_right (snp_def,i,winsize,nlocis,use_best))
	    != ULONG_MAX)) {
      win_len = stop_num - start_num + 1;
      win_max_len = MAX(win_len,win_max_len);

      //number of SNPs to the left of the central SNP
      snp_def[i]->winlen_left = (unsigned int) (i - start_num);
      //number of SNPs to the right of the central SNP
      snp_def[i]->winlen_right = (unsigned int) (stop_num - i);
      snp_def[start_num]->win_starting_here =
	cons(pair((void *)i, (void *)win_len),
	     snp_def[start_num]->win_starting_here);
    }
  }
  return win_max_len;
}

int is_a2_ancestral (snp_t *snp, unsigned long nb_a2_holder, unsigned int nb_hap) {
  return ((snp->a2 == snp->ancestor) ||
	  ((snp->a1 != snp->ancestor) && (2*nb_a2_holder > nb_hap)));
}

int is_a2_MAJ (unsigned long nb_a2_holder, unsigned int nb_hap) {
  return ( 2*nb_a2_holder > nb_hap );
}


int is_ancestral_defined (snp_t *snp ) {
  
        if( snp->ancestor == 'A' ){
	return 1;
  }else if( snp->ancestor == 'T' ){
  	return 1;  	
  }else if( snp->ancestor == 'G' ){
  	return 1;  	
  }else if( snp->ancestor == 'C' ){
  	return 1;  	
  }else{
  	return 0;
  } 
  
  
}
