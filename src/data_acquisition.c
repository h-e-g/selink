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
#include <err.h>

#include <stdlib.h>
#include <string.h>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <limits.h>
#include <stdio.h>

#include "minilib.h"
#include "data_acquisition.h"

#define BUFFLEN 100

// Free a sample cell
void free_sample (sample_t *s) {
  free(s->id);
  free(s);
}

// Free the samples struct
void free_samples (samples *s) {
  unsigned int i;
  for (i=0;i<s->pop_number;i++)
    free(s->pops[i]);
  free(s->pops);
  for (i=0;i<s->size;i++)
    free_sample(s->samples[i]);
  free(s->samples);
  free(s);
}

/* Parse [filename] with respect to grammar */
/* {pop_grp_name {: comma_separated_list(pop_name)}\n} list*/
static conj *parse_pop_filter(char * filename) {
  size_t buff_size = BUFFLEN*sizeof(char);
  FILE *file=safe_fopen(filename,"r");
  char* buff;
  if (!(buff=malloc(buff_size))) err(EXIT_FAILURE, "malloc failed");
  list *out_grp=nil(), *out_name=nil();
  samples *out;
  if (!(out=malloc(sizeof(samples)))) err(EXIT_FAILURE, "malloc failed");

  char *line, *pop_grp, *pop, *tmp = NULL;
  unsigned int nb=0, *id;
  while (getline(&buff,&buff_size,file) != -1) {
    line=buff;
    id = malloc(sizeof(int));
    *id = nb;
    //parse the pop_grp
    pop_grp=strdup(strsep(&line,":\r\n"));
    Right_strip(pop_grp);
    out_name = cons(pop_grp, out_name);
    //if there is a colon, parse the pops to merge
    if (line[0]==':') {
      while ((pop = strsep(&line,",\r\n")) != NULL) {
	tmp=strdup(pop);
	Right_strip(tmp);//TODO at left too
	out_grp = cons(pair(tmp,id), out_grp);
      }
    }
    out_grp = cons(pair(strdup(pop_grp),id), out_grp);
    nb++;
  }
  safe_fclose(filename,file);
  free(buff);
  out->pop_number=nb;
  out->pops= (char **) destructive_rev_to_array(out_name,nb);
  return pair(out,out_grp);
}

sample_t *dummy_sample() {
  sample_t *sample=malloc(sizeof(sample_t));
  if (sample) {
    sample->id = NULL;
    sample->pop_id = UINT_MAX;
    sample->sex = 0;
    sample->rank= UINT_MAX;
    sample->haplo1 = sample->haplo2 = NULL;
  } else err(EXIT_FAILURE, "unable to allocate a sample struct");
  return sample;
}

static sample_t *select_pop(list *pop_filter, list **out_pops, sample_t *active_sample, char *element) {
  unsigned int *pop_id;

  if (pop_filter) {
    pop_id = (unsigned int*) conservative_assoc((int (*)(const void*, const void*)) strcmp,element,pop_filter);
    if (pop_id) {//Select the sample
      active_sample->pop_id = *pop_id;
      //ensure NULLity of haplos
      active_sample->haplo1 = NULL;
      active_sample->haplo2 = NULL;
      return active_sample;
    } else {//pop not in pop_filter, drop the sample
      free_sample(active_sample);
      return NULL;
    } } else {//no pop_filter, create on the fly pops
    active_sample->pop_id = (unsigned int)
      imperative_find((int (*)(const void*, const void*)) strcmp,
		      (void* (*)(void *)) strdup,element,out_pops);
    //ensure NULLity of haplos
    active_sample->haplo1 = NULL;
    active_sample->haplo2 = NULL;
    //incr list
    return active_sample;
  }
}

// parse [prefix].sample with
/* format : id  sex population */
/*        : seaprator tab */
/*        : sex either 1 or 2 : 1 male, 2 female */
// return them in a struct.
static void parse_samples (char *prefix, list *pop_filter, samples* out) {
  size_t buff_size = BUFFLEN*sizeof(char);
  char* buff;
  if (!(buff=malloc(buff_size))) err(EXIT_FAILURE, "malloc failed");
  char* filename = append_extension(prefix,"sample");
  FILE *file = safe_fopen(filename, "r");

  unsigned int sample_nb = 0, rank = 0;
  sample_t *active_sample;
  char *element, *line;
  list* out_list = nil(), *out_pops = nil();

  while (getline(&buff,&buff_size,file) != -1) {
    line = buff;
    //init cell
    active_sample=malloc(sizeof(sample_t));
    if (!active_sample) err(EXIT_FAILURE, "malloc failed in parse_samples");
    //set rank
    active_sample->rank = rank;
    //parse id
    element=strsep(&line," \t");
    if (element) active_sample->id = strdup(element);
    else err(EXIT_FAILURE, "Badly formed sample file");
    //parse sex
    element=strsep(&line," \t");
    if (strcmp(element,"1")) active_sample->sex = 1;
    else if (strcmp(element,"2")) active_sample->sex = 2;
    else err(EXIT_FAILURE, "Unknown sex \"%s\" in sample file",element);
    //parse pop
    element=strsep(&line," \t\r\n");
    if (element) {
      if (select_pop(pop_filter,&out_pops,active_sample,element)) {
	out_list = cons(active_sample, out_list);
	sample_nb++;
      }
    } else err(EXIT_FAILURE, "Badly formed sample file");
    rank++;
  }

  safe_fclose(filename,file);
  free(buff);
  free(filename);

  out->next_snp_of_sample_offset = 0;
  out->size = sample_nb;
  out->max_rank = rank;
  out->samples = (sample_t**) destructive_rev_to_array(out_list,out->size);

  if (pop_filter)
    while (pop_filter) {
      free(proj1(pop_filter->head));
      free(proj2(pop_filter->head));
      free(pop_filter->head);
      destructive_tail(&pop_filter);
    }
  else {
    out->pop_number = (unsigned int)
      list_length(out_pops);
    out->pops = (char **) destructive_to_array(out_pops,out->pop_number);
  }
}

/* reads one line of .legend files. returns a *snp_t or NULL if it fails */
/* descrition read */
/* format : id      pos     allele1      allele2    ancestral */
/*        : separator tab or spaces */
/* in case ancestor is not known, it supplys - as value */
static snp_t *parse_one_legend(FILE *file, char *buff, size_t buff_size) {
  long l, i;
  char *p, *q, *tmp[5];
  snp_t *c = NULL;
  while ((l = getline(&buff, &buff_size, file)) == 0);//Read until an none empty line
  if (l != -1) {//if end of file is not reached
    q = buff;
    i = 0;
    while ((p = strsep (&q, " \t")) != NULL) {
      tmp[i]= p;
      i++;
    }
    if (i > 5 || i < 4)
      err(EXIT_FAILURE, "incorrect legend definition: %s", buff);
    if ((c = malloc(sizeof(snp_t)))) {
      c->id = strdup(tmp[0]);
      c->pos = ATOUI(tmp[1]);
      c->a1 = *tmp[2];
      c->a2 = *tmp[3];
      /* is ancestor defined ? if yes grab the data else set it to default */
      if (i > 4) c->ancestor = *tmp[4];
      else c->ancestor = '-';
      //every information about windows are set to -1 for now
      c->winlen_left = UINT_MAX;
      c->winlen_right = UINT_MAX;
      c->win_starting_here = nil();
    } else err(EXIT_FAILURE, "unable to allocate a snp struct");
  }
  return c;
}

/* This fuction fails if lines have more then BUFFLEN characters. */
conj *Legend_reader(char *prefix) {
  size_t buff_size = BUFFLEN*sizeof(char);
  char *filename = append_extension(prefix, "legend");
  FILE *file = safe_fopen(filename,"r");
  char * buff;
  if (!(buff=malloc(buff_size))) err(EXIT_FAILURE, "malloc failed");

  snp_t *current;
  unsigned long n=0;
  list *out = nil();

  //Drop the header line
  if (getline(&buff, &buff_size, file) == -1) err(EXIT_FAILURE, "%s: empty file", filename );

  while ((current = parse_one_legend(file,buff,buff_size)) != NULL) {
    out = cons(current,out);
    n++;
  }

  safe_fclose(filename,file);
  free(buff);
  free(filename);

  return pair(destructive_rev_to_array(out,n),(void *) n);
}

// Check .hap or .phased file is well formed and init haplo positions and offset.
void Init_Phased_reader(char *prefix, int transp, unsigned long nb_snp, samples *samples) {
  char *filename, *addr;
  int fd, has_space, has_training_space = 0;
  unsigned long i, height;
  unsigned long j, length;
  long long off;
  struct stat sd;

  // Virtually load the file to memory
  if (transp) filename = append_extension(prefix, "hap");
  else  filename = append_extension(prefix, "phased");
  if ((fd = open(filename,O_RDONLY)) == -1)
    err(EXIT_FAILURE, "%s: open failed", filename);
  if ((fstat(fd,&sd)) == -1)
    err(EXIT_FAILURE, "%s: failed to get stats", filename);
  addr = mmap(NULL, sd.st_size, PROT_READ, MAP_SHARED, fd, 0);
  if (addr == MAP_FAILED)
    err(EXIT_FAILURE, "%s: failed to load in memory", filename);

  // Get its structure
  has_space = ((sd.st_size >= 2)&&(addr[1] == ' '));
  if (transp) {
    height = nb_snp;
    length = has_space?4*samples->max_rank:2*samples->max_rank+1;
  } else {
    height = 2*samples->max_rank;
    length = has_space?2*nb_snp:nb_snp+1;
  }
  if (has_space) {
    has_training_space = (sd.st_size >= length)&&(addr[length-1] == ' ');
    if (has_training_space) length++;
  }

  // Check its structure
  off = 0;
  for (i=0;i<height;i++) {
    if ((off >= sd.st_size)
	|| (*(addr+off) != A1 && *(addr+off) != A2))
      err(EXIT_FAILURE,"badly formed %s ligne %lu column 1: %c",
	  filename,i,*(addr+off));
    j=1; off++;
    while (j<length-1) {
      if (has_space) {
	if ((off >= sd.st_size) || (*(addr+off) != ' '))
	  err(EXIT_FAILURE,"badly formed %s ligne %lu column %lu: %c",
	      filename,i,j+1,*(addr+off+1));
	else { off++; j++; }
      }
      if ((off >= sd.st_size)
	  || (*(addr+off) != A1 && *(addr+off) != A2))
	err(EXIT_FAILURE,"badly formed %s ligne %lu column %lu: %c",
	    filename,i,j,*(addr+off));
      else { off++; j++; }
    }
    if (has_training_space) {
      if ((off >= sd.st_size) || (*(addr+off) != ' '))
	err(EXIT_FAILURE,"badly formed %s ligne %lu column %lu: %c",
	    filename,i,j+1,*(addr+off+1));
      else { off++; j++; }
    }
    if ((off >= sd.st_size) || (*(addr+off) != '\n'))
      err(EXIT_FAILURE,"badly formed %s ligne %lu column %lu",
	  filename,i,j);
    else { off++; j++; }
  }

  free(filename);

  // Set the beginning of snp and next_snp_of_sample_offset
  if (transp) {
    samples->next_snp_of_sample_offset = length;
    for (i=0;i<samples->size;i++) {
      samples->samples[i]->haplo1 =
	addr + (has_space?2:1)*2*samples->samples[i]->rank;
      samples->samples[i]->haplo2 =
	addr + (has_space?2:1)*(2*samples->samples[i]->rank + 1);
    }
  } else {
    samples->next_snp_of_sample_offset = has_space?2:1;
    for (i=0;i<samples->size;i++) {
      samples->samples[i]->haplo1 =
	addr + length*2*samples->samples[i]->rank;
      samples->samples[i]->haplo2 =
	addr + length*(2*samples->samples[i]->rank + 1);
    }
  }
}

/* wrapper to get datas from .legend, .phased, .sample (filtered by popfile)*/
unsigned long Tripli_read(char *prefix, char *popfile_name,
			 snp_t ***snp_res, samples **s_res, int transp){
  conj *legend, *pop_filter;
  unsigned long nlocis;
  samples *s;

  /* get .legend datas */
  legend = Legend_reader(prefix);
  nlocis = (unsigned long) proj2(legend);
  *snp_res = (snp_t **) proj1(legend);
  free(legend);
  /* get .sample datas */
  if (popfile_name) {
    pop_filter = parse_pop_filter(popfile_name);
    s = (samples *) proj1(pop_filter);
    parse_samples(prefix,proj2(pop_filter),s);
  } else {
    if (!(s=malloc(sizeof(samples)))) err(EXIT_FAILURE, "malloc failed");
    parse_samples(prefix,nil(),s);
  }
  /* get .phased datas */
  /* fill the sample holder with haplotypes */
  Init_Phased_reader(prefix, transp, nlocis, s);

  *s_res = s;

  return nlocis;
}

void free_snpcol (unsigned long nlocis, snp_t **col){
  unsigned long i;

  for (i=0;i<nlocis;i++) {
    while (col[i]->win_starting_here)
      {
	free(col[i]->win_starting_here->head);
	destructive_tail(&(col[i]->win_starting_here));
      }
    free(col[i]->id);
    free(col[i]);
  }
  free(col);
}
