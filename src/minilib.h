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

#ifndef __MINILIB_H_
#define __MINILIB_H_

#include <stdio.h>

// --- Pair ---
typedef struct cconj {
  void *fst;
  void *snd;
} conj;

// Build a conj cell
conj *pair(void *fst, void *snd);

// Get projections
void *proj1(const conj *pair);
void *proj2(const conj *pair);

// --- List ---
typedef struct clist list;
struct clist {
  void* head;
  list* tail;
};

// Build the empty list
list* nil(void);
// Create a list cell
list* cons(void* head, list* tail);
// Replace a list by its tail, free the cell
void destructive_tail(list** l);
// Get the number of elements
unsigned long list_length(const list *l);

// Create an array from a list, reverse the element and free the list
void** destructive_rev_to_array(list* l, unsigned long list_length);
// Create an array from a list, free the list
void** destructive_to_array(list* l, unsigned long list_length);
// Append elements of l1 to l2 in reverse order
list* destructive_rev_append(list *l1, list * l2);

// Return the first element where f is true or NULL
void* conservative_filter(int (f)(const void*), const list *l);
// Associative list functions
// returns NULL if not present
void* conservative_assoc(int (cmp)(const void*,const void*),
			 const void *el, const list *l);
//Returns the rank of [el] in [l]. Add [el] at the end of [l] in place if [el]
//was not present in [l].
unsigned long imperative_find
(int (cmp)(const void *,const void*), void* (cp)(void *), void *el, list **l);

// --- FIFO stack with a maximum length ---
typedef struct cstack {
  unsigned long start;
  unsigned long stop;
  unsigned long size;
  void *data[];
} borned_stack;

// Create an empty stack
borned_stack* create_stack(unsigned long size);
// Push an element in sk
void push(borned_stack* sk, void* el);
// Remove and Return the first element of the stack
void* pop(borned_stack* sk);
// Free a stack. DO NOT free its elements.
void free_stack(borned_stack* sk);
// Reinit a stack to the empty one. DO NOT do anything with elements
void reinit_stack(borned_stack* sk);
int stack_is_empty_or_full(const borned_stack* sk);
// Get the number of element
unsigned long non_empty_stack_length(const borned_stack* sk);

// Apply [f] to [r] and all elements of sk (in the order of pop).
void iter_non_empty_stack(const borned_stack* sk, void *r, void (f)(void *,void *));

// --- File(name)s ---
//fopen but ensures that the file descriptor is not NULL
FILE *safe_fopen(char *filename, char *mode);
//fclose but exit if it fails
void safe_fclose(char *filename,FILE *file);
/* generate prefix.ext names */
char *append_extension(char *prefix, const char *ext);

// --- Dummy ---
/* right strip */
int Right_strip(char *target);

#define ATOUI(x) ((unsigned int) strtoul(x,NULL,10))
#endif
