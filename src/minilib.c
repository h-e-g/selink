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
#include <ctype.h>
#include <err.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "minilib.h"

conj *pair(void *fst, void *snd) {
  conj *cell = malloc(sizeof(conj));
  if (cell) {
    cell->fst = fst;
    cell->snd = snd;
  } else {
    err(EXIT_FAILURE, "unable to allocate a conj cell");
  }
  return cell;
}

void *proj1(const conj *pair) {
  if(pair) return pair->fst;
  else err(EXIT_FAILURE, "Tried to get proj1 of uninitialized conj");
}

void *proj2(const conj *pair) {
  if(pair) return pair->snd;
  else err(EXIT_FAILURE, "Tried to get proj2 of uninitialized conj");
}

list* nil() {
  return NULL;
}

list* cons(void* head, list* tail) {
  list* cell = malloc(sizeof(list));
  if (cell) {
    cell->head=head;
    cell->tail=tail;
  } else {
    err(EXIT_FAILURE, "unable to allocate a list cell");
  }
  return cell;
}

void destructive_tail(list** l) {
  list* l2;
  l2=*l;
  *l=(*l)->tail;
  free(l2);
}

void *conservative_filter(int (f)(const void*), const list *l) {
  const list *l2=l;
  while(l2&&f(l2->head)) l2 = l2->tail;
  if (l2) return (l2->head);
  else return NULL;
}

void *conservative_assoc(int (cmp)(const void*,const void*), const void *el, const list *l) {
  const list *l2=l;
  while(l2 && cmp(el,proj1(l2->head))) l2 = l2->tail;
  if (l2) return (proj2(l2->head));
  else return NULL;
}

void** destructive_rev_to_array(list* l, unsigned long length) {
  void** out=malloc(length*sizeof(void*));
  while(length>0&&l) {
    out[--length]=l->head;
    destructive_tail(&l);
  }
  if (length>0||l) err(EXIT_FAILURE, "length was incorrect in destructive_rev_to_array");
  else return out;
}

void** destructive_to_array(list* l, unsigned long length) {
  unsigned long i=0;
  void** out=malloc(length*sizeof(void*));
  while(i<length&&l) {
    out[i]=l->head;
    i++;
    destructive_tail(&l);
  }
  if (i<length||l) err(EXIT_FAILURE, "length was incorrect in destructive_to_array");
  else return out;
}

unsigned long imperative_find
(int (cmp)(const void *,const void*), void* (cp) (void *), void *el, list **l) {
  list *l2=*l;
  unsigned long i = 0;
  if (l2) {
    while (l2->tail && cmp(el,l2->head)) {
      l2 = l2->tail;
      i++;
    }
    if (cmp(el,l2->head)) {
      if (l2->tail) err(EXIT_FAILURE,"Big problem in imperative_find\n");
      l2->tail = cons(cp(el),nil());
      i++;
    }
  }
  else
    *l=cons(cp(el),nil());
  return i;
}

static unsigned long list_length_aux(const list *l, unsigned long i) {
  if (l) return list_length_aux(l->tail,i+1);
  else return i;
}

unsigned long list_length(const list *l) {
  return list_length_aux(l,0);
}

list *destructive_rev_append(list *l1, list *l2) {
  list *out = l1;

  while (l2) {
    out = cons(l2->head,out);
    destructive_tail(&l2);
  }
  return out;
}

void reinit_stack(borned_stack* sk) {
  sk->start = 0;
  sk->stop = 0;
}

borned_stack* create_stack(unsigned long size){
  borned_stack *sk =
    malloc(sizeof(borned_stack)+sizeof(void *)*size);
  if((size>0)&&sk) {
    reinit_stack(sk);
    sk->size = size;
  } else
    err(EXIT_FAILURE, "unable to allocate a borned stack");
  return sk;
}

void free_stack(borned_stack* sk) {
  if (sk)
    free(sk);
}

int stack_is_empty_or_full(const borned_stack* sk) {
  if (sk) return (sk->stop == sk->start);
  else return 0;
}

unsigned long non_empty_stack_length(const borned_stack* sk) {
  if (sk->stop > sk->start)
    return sk->stop - sk->start;
  else
    return sk->size - sk->start + sk->stop;
}

void push(borned_stack* sk, void* el){
  if (sk) {
    sk->data[sk->stop] = el;
    if (sk->stop == sk->size-1) sk->stop = 0;
    else sk->stop++;
  } else
    err(EXIT_FAILURE, "cannot push in an uninitialized stack");
}

void* pop(borned_stack* sk) {
  void *tmp;

  if (sk) {
    tmp = sk->data[sk->start];
    if (sk->start == sk->size - 1) sk->start = 0;
    else sk->start++;
  } else
    err(EXIT_FAILURE, "cannot pop from an uninitialized stack");
  return tmp;
}

void iter_non_empty_stack(const borned_stack* sk, void *r, void (f)(void *,void *)) {
  unsigned long i;
  if (sk) {
    if (sk->stop > sk->start) {
      for (i=sk->start;i<sk->stop;i++) f(r,sk->data[i]);
    } else {
      for (i=sk->start;i<sk->size;i++) f(r,sk->data[i]);
      for (i=0;i<sk->stop;i++) f(r,sk->data[i]);
    }
  }
}


//fopen but ensures that the file descriptor is not NULL
FILE *safe_fopen(char *filename, char *mode) {
  FILE *file = fopen(filename,mode);
  if (file) return file;
  else err(EXIT_FAILURE, "%s: open failed", filename);
}

//fclose but exits if it fails
void safe_fclose(char *filename,FILE *file) {
  if (fclose(file) == EOF)
    err(EXIT_FAILURE, "%s: close failed", filename);
}

/* generate file names based on prefix + extension */
char *append_extension(char *prefix, const char *ext) {
  char *res;
  if (asprintf(&res, "%s.%s", prefix, ext) < 0)
    err(EXIT_FAILURE, "extension generator");
  return res;
}

/* right strip target */
int Right_strip(char *target) {
  int success = false;
  char *p;
  if ((p = strchr(target, '\n')) != NULL){
    success = true;
    while ((p >= target) && (isspace(*p) || (*p == '\0'))) {
      *p ='\0';
      p--;
    }
  }
  return success;
}
