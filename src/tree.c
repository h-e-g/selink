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

#include <err.h>
#include <stdlib.h>

#include "data_acquisition.h"
#include "tree.h"

#define NEXT_POS(rev,offset,pos) ((rev)?(pos)-offset:(pos)+offset)

static tree_t *mk_empty() {
  return NULL;
}

static void tree_update_node(tree_t *node, tree_t *a1, tree_t *a2, unsigned int nb) {
  node->typ = NODE;
  node->val.node.a1 = a1;
  node->val.node.a2 = a2;
  node->val.node.nb = nb;
}

static tree_t *mk_stuck(int rev, unsigned long remains,
			unsigned long offset, list *el) {
  tree_t *n;
  if (!(n = malloc(sizeof(tree_t))))
    err(EXIT_FAILURE, "malloc of a tree failed");
  n->typ = STUCK;
  n->val.stuck.remains = remains;
  n->val.stuck.rev = rev;
  n->val.stuck.offset = offset;
  n->val.stuck.elements = el;
  return n;
}

void tree_compute_node(tree_t *node) {
  int rev;
  unsigned int nb = 0;
  unsigned long i, off;
  char *pos;
  list *l, *ll = nil(), *lr = nil();
  tree_t *a1 = mk_empty(), *a2 = mk_empty();
  if (node && (node->typ == STUCK)) {
    l = node->val.stuck.elements;
    rev = node->val.stuck.rev;
    off = node->val.stuck.offset;
    i = node->val.stuck.remains;
    if (i == 0) {
      while (l) {
	destructive_tail(&l);
	nb++;
      }
    } else {
      while (l) {
	pos = (char *) l->head;
	if (*pos == A1)
	  ll = cons(NEXT_POS(rev,off,pos),ll);
	else if (*pos == A2)
	  lr = cons(NEXT_POS(rev,off,pos),lr);
	else
	  err(EXIT_FAILURE,
	      "Violated invariant in tree:compute_node (rev=%d, remains=%lu)",
	      rev,i);
	nb++;
	destructive_tail(&l);
      }
      if (ll) a1 = mk_stuck(rev,i-1,off,ll);
      if (lr) a2 = mk_stuck(rev,i-1,off,lr);
    }
  tree_update_node(node,a1,a2,nb);
  }
}

static void tree_insert(int rev, char *label, unsigned long length,
	    unsigned long offset, tree_t **tree) {
  if (*tree) {
    if ((*tree)->typ == STUCK) {
      if ((length == (*tree)->val.stuck.remains) && (rev == (*tree)->val.stuck.rev))
	(*tree)->val.stuck.elements = cons(label,(*tree)->val.stuck.elements);
      else err(EXIT_FAILURE,"Violated invariant in tree:insert.");
    } else {
      (*tree)->val.node.nb++;
      if (*label == A1)
	tree_insert(rev,NEXT_POS(rev,offset,label),length-1,
	       offset,&((*tree)->val.node.a1));
      else if (*label == A2)
	tree_insert(rev,NEXT_POS(rev,offset,label),length-1,
	       offset,&((*tree)->val.node.a2));
      else err(EXIT_FAILURE,"Violated invarient in tree:insert");

    }
  } else *tree = mk_stuck(rev,length,offset,cons(label,nil()));
}

static tree_t *destructive_merge(tree_t *t1, tree_t *t2) {
  if (t1) {
    if (t2) {
      if ((t1->typ == STUCK) && (t2->typ == STUCK)) {
	if ((t1->val.stuck.rev == t2->val.stuck.rev) &&
	    (t1->val.stuck.remains == t2->val.stuck.remains))
	  t1->val.stuck.elements =
	    destructive_rev_append(t1->val.stuck.elements,t2->val.stuck.elements);
	else err(EXIT_FAILURE,"Violated invariant in tree:destructive_merge.");
      } else {
	tree_compute_node(t1);
	tree_compute_node(t2);
	t1->val.node.nb = t1->val.node.nb + t2->val.node.nb;
	t1->val.node.a1 = destructive_merge(t1->val.node.a1,t2->val.node.a1);
	t1->val.node.a2 = destructive_merge(t1->val.node.a2,t2->val.node.a2);
      }
      free(t2);
    }
    return t1;
  } else return t2;
}

void free_computation_pack(computation_pack *pack) {
  free_stack(pack->current_sk);
  free_stack(pack->next_sk);
  if (pack->current_tree)
    err(EXIT_FAILURE, "To free computation_pack contains an non empty tree\n");
  free(pack);
}

void tree_pop_root(computation_pack *pack) {
  tree_t *tmp=pack->current_tree;
  if (tmp) {
    tree_compute_node(tmp);
    pack->current_tree = destructive_merge(tmp->val.node.a1,tmp->val.node.a2);
    free(tmp);
  }
}

computation_pack *build_tree(int rev, unsigned long pop_id,
			     samples *smpls, unsigned long nlocis) {
  tree_t *out = mk_empty();
  computation_pack *pack;
  unsigned int i, nb=0;

  for (i = 0; i < smpls->size; i++)
    if (smpls->samples[i]->pop_id == pop_id) {
      tree_insert(rev,
	     (rev?(nlocis-1)*smpls->next_snp_of_sample_offset:0)
	     + smpls->samples[i]->haplo1,nlocis,
	     smpls->next_snp_of_sample_offset, &out);
      tree_insert(rev,
	     (rev?(nlocis-1)*smpls->next_snp_of_sample_offset:0)
	     + smpls->samples[i]->haplo2,nlocis,
	     smpls->next_snp_of_sample_offset, &out);
      nb++;
    }
  if (!(pack = malloc(sizeof(computation_pack))))
    err(EXIT_FAILURE, "malloc of a computation pack failed");

  pack->current_tree = out;
  pack->current_sk = create_stack(2*nb+1);
  pack->next_sk = create_stack(2*nb+1);

  return pack;
}

void switch_stacks(computation_pack *pack) {
  borned_stack *tmp = pack->current_sk;
  pack->current_sk=pack->next_sk;
  pack->next_sk=tmp;
}
