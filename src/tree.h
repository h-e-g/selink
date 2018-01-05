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

#include "data_acquisition.h"

typedef struct tree {
  enum { NODE = 1, STUCK = 2 } typ;
  int _dummy_pad;
  union {
    struct {
      unsigned long remains;
      unsigned long offset;
      int rev;
      int _dummy_pad;
      list *elements;
    } stuck;
    //Left is a1, right is a2
    struct {
      struct tree *a1;
      struct tree *a2;
      unsigned int nb;
      int _dummy_pad;
    } node;
  } val;
} tree_t;

void tree_compute_node(tree_t *node) ;

typedef struct {
  borned_stack *current_sk;
  borned_stack *next_sk;
  tree_t *current_tree;
} computation_pack;

void free_computation_pack(computation_pack *pack) ;
void switch_stacks(computation_pack *pack) ;

void tree_pop_root(computation_pack *pack) ;

computation_pack *build_tree(int rev, unsigned long pop_id,
			     samples *smpls, unsigned long nlocis) ;
