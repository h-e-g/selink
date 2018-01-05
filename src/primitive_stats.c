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
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <err.h>

#include "tree.h"
#include "stats.h"
#include "primitive_stats.h"

static void one_next_nb_a2(void *val, void *el) {
  tree_t *n, *d;
  conj *p = (conj*) val;
  unsigned int *out = (unsigned int*) proj1(p);
  borned_stack *next_sk = (borned_stack*) proj2(p);

  n = (tree_t*) el;
  if (n) {
    tree_compute_node(n);
    d = n->val.node.a2;
    if (d) {
      tree_compute_node(d);
      push(next_sk,d);
      *out += d->val.node.nb;
    }
    if (n->val.node.a1)
      push(next_sk,n->val.node.a1);
  }
}

static unsigned int next_nb_a2(computation_pack *pack) {
  unsigned int out_nb_a2 = 0;
  conj p = { &out_nb_a2, pack->next_sk };
  iter_non_empty_stack(pack->current_sk,&p,one_next_nb_a2);
  reinit_stack(pack->current_sk);
  switch_stacks(pack);
  return out_nb_a2;
}

static void one_next_ehh_and_nb_a2(void *val, void *el) {
  tree_t *n, *d;
  conj *p = (conj*) val;
  conj *pp = (conj*) proj1(p);
  unsigned long *out_ehh_num = proj1(pp);
  unsigned int *out_nb_a2 = proj2(pp);
  borned_stack *next_sk = (borned_stack*) proj2(p);

  n = (tree_t*) el;
  if (n) {
    tree_compute_node(n);
    *out_ehh_num += n->val.node.nb * (n->val.node.nb - 1);
    //printf("%i ",*out_ehh_num);
    d = n->val.node.a2;
    if (d) {
      tree_compute_node(d);
      push(next_sk,d);
      *out_nb_a2 += d->val.node.nb;
    }
    if (n->val.node.a1)
      push(next_sk,n->val.node.a1);
  }
}

static void next_ehh_num_and_nb_a2
(computation_pack *pack, unsigned long *out_ehh_num, unsigned long *out_nb_a2) {
  *out_ehh_num = 0;
  *out_nb_a2 = 0;
  conj out = {out_ehh_num, out_nb_a2};
  conj val = {&out, pack->next_sk};

  iter_non_empty_stack(pack->current_sk,&val,one_next_ehh_and_nb_a2);
  reinit_stack(pack->current_sk);
  switch_stacks(pack);
  //printf("\n");
}

static void one_cross_one_hap(void *val, void *el) {
  conj *p = (conj*) val;
  unsigned long *nb_singleton = (unsigned long *) proj1(p);
  borned_stack *next_sk = (borned_stack*) proj2(p);

  tree_t *n = (tree_t*) el;

  if(n) {
    tree_compute_node(n);
    if (n->val.node.nb==1) (*nb_singleton)++;
    else {
      if (n->val.node.a2)
	push(next_sk,n->val.node.a2);
      if (n->val.node.a1)
	push(next_sk,n->val.node.a1);
    }
  }
}

static void cross_one_hap(computation_pack *pack, unsigned long *nb_singleton) {
  conj out = {nb_singleton, pack->next_sk};
  iter_non_empty_stack(pack->current_sk,&out,one_cross_one_hap);
  reinit_stack(pack->current_sk);
  switch_stacks(pack);
}

static void one_diversity_sum(void *val, void *el) {
  unsigned long *haplo_sum = (unsigned long*) val;

  tree_t *n = (tree_t*) el;

  if (n) {
    tree_compute_node(n);
    *haplo_sum += n->val.node.nb * n->val.node.nb;
  }
}

static unsigned long diversity_sum(computation_pack *pack,
				   unsigned long haplo_sum) {
  iter_non_empty_stack(pack->current_sk,&haplo_sum,one_diversity_sum);
  return haplo_sum;
}

static void prepare_primitive_computations(computation_pack *pack) {
  tree_compute_node(pack->current_tree);
  reinit_stack(pack->current_sk);
  push(pack->current_sk,pack->current_tree);
}

static void compute_diversities(computation_pack *pack, list *wins,
				unsigned long *pop_div,
				unsigned long *pop_haplo_sum, unsigned int i) {
  unsigned long tmp_core, tmp_len, nb_singleton=0;
  conj *p;

  if (pop_div&&pop_haplo_sum) {
    while(wins) {
      p = (conj*) wins->head;
      tmp_core = (unsigned long) proj1(p);
      tmp_len = (unsigned long) proj2(p);
      while((i<tmp_len)&&!(stack_is_empty_or_full(pack->current_sk))) {
	cross_one_hap (pack,&nb_singleton);
	i++;
      }
      if (stack_is_empty_or_full(pack->current_sk)) {
	pop_div[tmp_core] += nb_singleton;
	pop_haplo_sum[tmp_core] += nb_singleton;
      } else {
	pop_div[tmp_core] +=
	  nb_singleton + non_empty_stack_length(pack->current_sk);
	pop_haplo_sum[tmp_core] += diversity_sum(pack,nb_singleton);
      }
      wins = wins->tail;
    }
  }
}

static unsigned int half_ihh_and_pi(computation_pack *pack, double min_value,
				    unsigned int max_size, int rev,
				    unsigned long start, unsigned int pop_len,
				    tree_t *tree, snp_t **snp, double *out_ihh,
				    double *out_pi) {
  unsigned int i = 0, nb_carrier;
  unsigned long nb_a2, out_ihh_num=0;
  unsigned long prev_ehh_num, ehh_num, tmp_pi = 0.0;

  tree_compute_node(tree);
  reinit_stack(pack->current_sk);
  push(pack->current_sk,tree);
  if (tree) {
    if ((nb_carrier = tree->val.node.nb) > 1) {
      // compute ehh and nb_a2
      next_ehh_num_and_nb_a2(pack,&prev_ehh_num,&nb_a2);
      if (out_pi)
	tmp_pi += (nb_carrier-nb_a2/*=nb_a1*/)*nb_a2;
      if (prev_ehh_num != nb_carrier*(nb_carrier -1))
	err(1,"Fists ehh value is %f not 1.0\n",
	    ((double) prev_ehh_num)/nb_carrier*(nb_carrier -1));
      if (out_pi || out_ihh)
	while (!stack_is_empty_or_full(pack->current_sk) && i < max_size
	       && prev_ehh_num > (min_value-1.0/pop_len)*(nb_carrier*(nb_carrier-1))) {
	  next_ehh_num_and_nb_a2(pack,&ehh_num,&nb_a2);
	  if (out_pi)
	    tmp_pi += (nb_carrier-nb_a2/*=nb_a1*/)*nb_a2;
	  if (out_ihh)
	    out_ihh_num += (prev_ehh_num + ehh_num)
	      * (snp[rev?start-i:start+i+1]->pos - snp[rev?start-(i+1):start+i]->pos);
	  //if (rev)
	  //printf("in %i, ehh = %f. From now, ihh = %f\n",snp[start+(rev?-i:i+1)]->pos,ehh,*out_ihh);
	  prev_ehh_num = ehh_num;
	  i++;
	}
      if (out_ihh) {
	if (max_size==UINT_MAX&&stack_is_empty_or_full(pack->current_sk))
	  *out_ihh = NAN;
	else
	  *out_ihh += 0.5 * (((double) out_ihh_num)/(nb_carrier*(nb_carrier-1))) +
	    (1.0/pop_len -0.05)*(snp[rev?start:start+i]->pos - snp[rev?start-i:start]->pos);
      }
    }

    // compute end of nb_a2 when ehh is < min_value
    if (out_pi)
      for(;i<max_size;i++) {
	nb_a2 = next_nb_a2(pack);
	tmp_pi += (nb_carrier-nb_a2/*=nb_a1*/)*nb_a2;
      }
    if (out_pi)
      *out_pi += tmp_pi/sum_w_from_type(2,nb_carrier);
  }
  return i;
}

static void **create_array_for_results(samples *smpls, unsigned long nlocis,
				       size_t size) {
  void **res, *out;
  unsigned int pop_id;

  if ((res = malloc(sizeof(void*) * smpls->pop_number)) == NULL) err(EXIT_FAILURE, "malloc failed");

  for (pop_id=0; pop_id<smpls->pop_number; pop_id++) {
    if ((out = calloc(nlocis,size)) == NULL) err(EXIT_FAILURE, "malloc failed");
    res[pop_id] = out;
  }
  return res;
}

double **create_double_array_for_results(samples *smpls,unsigned long nlocis) {
  return (double **) create_array_for_results(smpls,nlocis,sizeof(double));
}

unsigned long **create_long_array_for_results(samples *smpls,unsigned long nlocis) {
  return (unsigned long**)
    create_array_for_results(smpls,nlocis,sizeof(unsigned long));
}

void free_double_results_array(samples *smpls,double **t) {
  unsigned int i;
  if (t) {
    for (i=0; i<smpls->pop_number; i++)
      free(t[i]);
    free(t);
  }
}

void free_long_results_array(samples *smpls, unsigned long **t) {
  unsigned int i;
  if (t) {
    for (i=0; i<smpls->pop_number; i++)
      free(t[i]);
    free(t);
  }
}

void fill_half_ihh_pi_divs_and_nbd_arrays_one_pop
(int verbose, double min_value, int rev, unsigned long pop_id,
 unsigned int pop_len,samples *smpls, snp_t **snp,
 unsigned long nlocis, double *ihhA, double *ihhD,
 unsigned long *nb_a2_holder, double *piA, double *piD,
 unsigned long *diversities, unsigned long *haplo_sums) {
  unsigned long i, pos;
  unsigned int depth=0;
  unsigned int max_size;
  tree_t *tree_der;
  computation_pack *pack;

  if (rev && (diversities||haplo_sums))
    err(EXIT_FAILURE, "computing diversities does not work in reverse order");

  //printf ("compute ihh for pop %s in direction %i\n",smpls->pops[pop_id],rev);
  pack = build_tree(rev,pop_id,smpls,nlocis);
  if (verbose)
    fprintf(stderr, "done window in %s order for pop %lu/%d: 00%%",
	    rev?"reverse":"regular",pop_id+1,smpls->pop_number);
  for(i=0; i<nlocis; i++) {
    if (verbose && ((nlocis > 100) && (i % (nlocis/100) == 0)))
      fprintf (stderr,"\b\b\b%02lu%%",(100*i)/nlocis);
    pos = rev?nlocis-i-1:i;
    if (!pack->current_tree) err(EXIT_FAILURE, "fill_half_ihh_and_nbd_arrays builds a tree too short\n");
    tree_compute_node(pack->current_tree);
    if (nb_a2_holder) {
      tree_der = pack->current_tree->val.node.a2;
      tree_compute_node(tree_der);
      nb_a2_holder[pos] = tree_der?tree_der->val.node.nb:0;
    }
    max_size = rev?snp[pos]->winlen_left:snp[pos]->winlen_right;

    if (pack->current_tree->val.node.a1
	&& pack->current_tree->val.node.a2) {//window is not skipped
      //if (rev)
      //printf(" A1 of core %i\n",pos);
      depth = half_ihh_and_pi(pack,min_value,max_size,rev,pos,pop_len,
			      pack->current_tree->val.node.a1,snp,
			      ihhA?ihhA+pos:NULL,
			      (max_size!=UINT_MAX&&piA)?piA+pos:NULL);
      if (max_size != UINT_MAX&&!rev&&diversities&&haplo_sums)
	compute_diversities(pack,snp[pos]->win_starting_here,
			    diversities,haplo_sums,depth+2);

      //if (rev)
      //printf(" A2 of core %i\n",pos);
      depth = half_ihh_and_pi(pack,min_value,max_size,rev,pos,pop_len,
			      pack->current_tree->val.node.a2,snp,
			      ihhD?ihhD+pos:NULL,
			      (max_size!=UINT_MAX&&piD)?piD+pos:NULL);
      if (max_size != UINT_MAX&&!rev&&diversities&&haplo_sums)
	compute_diversities(pack,snp[pos]->win_starting_here,
			    diversities,haplo_sums,depth+2);
    }
    if ((!pack->current_tree->val.node.a1
	 || !pack->current_tree->val.node.a2
	 || max_size == UINT_MAX)
	&& !rev&&diversities&&haplo_sums) {
	prepare_primitive_computations(pack);
	compute_diversities(pack,snp[pos]->win_starting_here,
			    diversities,haplo_sums,0);
      }
    tree_pop_root(pack);
  }
  tree_pop_root(pack);
  free_computation_pack(pack);
  if (verbose) fprintf (stderr,"\b\b\ball.\n");
}
