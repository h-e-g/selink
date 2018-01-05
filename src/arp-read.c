#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#endif

#include <err.h>
#include <ctype.h>
#include <stdio.h>

#include "arp-read.h"
#include "utils.h"


#define BUFFLEN 100
#define NLOCISTAG  "NumLinkedLoci="
#define ANCESTORTAG "ANCESTOR_MUT="
#define SAMPLENBTAG "NbSamples="
#define SAMPLENAMETAG "SampleName="
#define SAMPLESIZETAG "SampleSize="
#define SAMPLEDATATAG "SampleData="
#define BLOCKCLOSE '}'

static unsigned long Arp_getintval(const char *s1, const char *s2) {
  char *p;
  unsigned long ret;
  ret = 0;
  if ((p = strcasestr(s1, s2)) != NULL) {
    p += strlen(s2);
    ret = ATOUI(p);
  }
  return ret;
}

static char *Arp_getstrval(const char *s1, const char *s2) {
  char *ret, *p;
  ret = NULL;
  
  if ((p = strcasestr(s1, s2)) != NULL){
    p += strlen(s2);
    ret = strdup(p);
  }
  return ret;
}

static unsigned int *Arp_getlocispos(char *buff, unsigned long nbval){

  unsigned int *pos; 
  unsigned long i;
  char *p, *q;
  size_t n, l, L;
  L = strlen(buff);
  /* split buffer */
  q = buff; 
  while ((p = strsep(&q, " \t\n")) != NULL) { continue; }
  
  if ((pos = malloc(sizeof(unsigned int) * nbval)) == NULL)
    err(EXIT_FAILURE, "malloc failed");

  /*extract values */
  i = l = 0;
  p = buff;
  while(l < L){
    if (*p) { //something to deal with
      pos[i] = ATOUI(p);
      i++;
      n = strlen(p);
      l += n;
      p += n; //jump
    }
    // anyway try to span 
    p++; l++;
  }
  if (i - nbval) err(EXIT_FAILURE, "arp positions inconsistent values"); 
  
  return pos;
}


static sample_t * Process_sample
(char *dat, unsigned int n, unsigned int hl, char *ref) {
  
  int i, haplo;
  char *p, *q, *s;
  size_t l;
  sample_t *sample;
  if ((sample = malloc(sizeof(sample_t))) == NULL) err(EXIT_FAILURE, "malloc failed");
  /* span empty char */
  p  = dat;
  while (*p && isblank(*p)) { p++; }
  s = p;
  l = strlen(p);
  
  /* which haplo we are dealing with */
  haplo = (n % 2) + 1;
  
  /* split datas */
  while ((q = strsep (&p, " \t")) != NULL) { }
  p =s;
  q = p + l;
  i = 0;
  if (haplo == 1) { // dealing with id val haplotype_1
	 while (p <q) {
		if (*p) {
		  switch (i) {
		  case 0: //get id
			 sample->id = strdup(p);
			 break;
		  case 1: //skip value
			 break;
		  case 2: // get haplotype
			 if (strlen(p) != hl) err(EXIT_FAILURE, "arp file inconsistent format");
			 Seq2bin(ref, &p); // convert to 0/1 according to ancestor definition
			 sample->haplo1 = strdup(p);
			 break;
		  default: // we expect only 3 values
			 err(EXIT_FAILURE, "arp file, inconsistent format");
		  }// end switch
		  i++;
		  p += strlen(p);
		}
		p++;
	 }
  }// end haplotype_1
  else { // haplotype_2
	 while (p <q) {
		if (*p) { 
		  switch (i) {
		  case 0:
			 if (strlen(p) != hl) err(EXIT_FAILURE, "arp file, inconsistent format");
			 Seq2bin(ref, &p);
			 sample->haplo2 = strdup(p);
			 break;
		  default: // we expect only 1 values
			 err(EXIT_FAILURE, "arp file, inconsistent format");
		  }
		  i++;
		  p += strlen(p);
		}
		p++;
	 }
  }// end haplotype_2
	 
  /* fill non available datas in .arp format, don't care doing it twice */
  sample->sex = 0;
  sample->pop_id = UINT_MAX;//TODO completly wrong
  return sample;

}


unsigned int Arp_reader(char *prefix, snp_t ***snp_res, sample_t ***sample_col) {
  
  FILE *IN;
  char *realfile;
  char *BUFF, *p;
  int success;
  unsigned int i, n, status, cntrl, sampleread,pos_in_sa=UINT_MAX;
  unsigned int nlocis=UINT_MAX, samplesize=UINT_MAX, totalsamples=0;
  unsigned int *pos = NULL;
  char *ancestor=NULL;
  size_t L;
  sample_t **sa = NULL;
  snp_t **snp_col, *snp; 
  

  /* generate .arp extension */
  realfile = append_extension(prefix, "arp");
  
  if ((IN = fopen(realfile, "r")) == NULL) err(EXIT_FAILURE, "%s: open failed", realfile);

  if ((BUFF=malloc(BUFFLEN)) == NULL) err(EXIT_FAILURE, "malloc failed"); 

  status = cntrl = n = sampleread = totalsamples =0;

  L = BUFFLEN;
  p = BUFF; 
  
  while (fgets(p, BUFFLEN, IN) != NULL) {
	 /* skip empty lines */
    if (*p == '\0') continue;
	 /* check if line completly read */
	 if (strchr(p, '\n') == NULL) {
		if ((BUFF = (char *)realloc(BUFF, (L+BUFFLEN) * sizeof(char))) == NULL) err(EXIT_FAILURE, "malloc failed"); 
		p = BUFF + strlen(BUFF);
		L += BUFFLEN;
		continue;
	 }
	 success = Right_strip(p);
	 p = BUFF;
    switch (status) {
    case 0: /* check for NlocisNumber definition  */
      if (strcasestr(p, NLOCISTAG) != NULL) { 
		  nlocis = Arp_getintval(p, NLOCISTAG);		 
		  status = 1; //jump ancestor acquisition
      }
      break;
		
    case 1: /*check for ancestor data */
      if (strcasestr(p, ANCESTORTAG) != NULL) { 
		  ancestor = Arp_getstrval(p, ANCESTORTAG);
		  /* check data coherency */
		  if (strlen(ancestor) != nlocis) err(EXIT_FAILURE, "%s: invalid format", realfile); 
		  status = 2; // jump to locis position acquisition
      }
		//because ancestor desc follow imediatly n locis description
      else err(EXIT_FAILURE, "%s: invalid format", realfile); 
      break;
      
    case 2: /* get locis pos following ancestor description */
      /* check if line fuly read, may be long */
      if (success == 0) { // we miss datas
		  continue;
      }
      else {
		  pos =Arp_getlocispos(BUFF, nlocis);
		  if ((snp_col = malloc(nlocis*sizeof(snp_t*))) == NULL) err(EXIT_FAILURE, "malloc failed");
		  i = 0;
		  while (i < nlocis) {
		    if ((snp = malloc(sizeof(snp_t))) == NULL) err(EXIT_FAILURE, "malloc failed");
		    asprintf(&(snp->id), "%u", pos[i]);
		    snp->pos = pos[i];
		    snp->ancestor = ancestor[i];
		    snp->a1 = snp->a2 = '\0';
		    snp_col[i]=snp;
		    i++;
		  }
		  *snp_res = snp_col;
		  status = 3; // jump to population number acquisition
      }
      break;
		
    case 3: /* get number of populations */
		if (strcasestr(p, SAMPLENBTAG) != NULL) { 
		  status = 4; // jump to sample name acquisition
		}
		break;
		  
    case 4: /* get sample name == population name description */
      if (strcasestr(p, SAMPLENAMETAG) != NULL) { 
		  status = 5; // jump to sample size acquiisition
      }
      break;

    case 5: /* get sample size */
		if (strcasestr(p, SAMPLESIZETAG) != NULL) { 
		  samplesize = Arp_getintval(p,SAMPLESIZETAG );
		  totalsamples += samplesize;
		  n = 0; 
		  // allocate sample holder according sample size
		  if (totalsamples == samplesize) {
			 if ((sa = malloc((samplesize+1)*sizeof(sample_t*))) == NULL)
				err(EXIT_FAILURE, "malloc failed");
			 pos_in_sa = 0;
		  }
		  else {
			 if ((sa = realloc(sa,(totalsamples+1)*sizeof(sample_t*))) == NULL)
				err(EXIT_FAILURE, "malloc failed"); 
			 pos_in_sa = totalsamples -  samplesize;
		  }
		  status = 6; // jump to sample block checking
		}
		break;

	 case 6: /* just check that sample data block exists  */
		if (strcasestr(p, SAMPLEDATATAG) != NULL) { 
		  cntrl = 1;
		  status = 7;
		} // jump to sample acquisition
		break;
		
	 case 7: /* get samples: should get 2n samplesize values */
		if (cntrl == 1 && (strchr(p, BLOCKCLOSE) != NULL)) {
		  if ((2*samplesize) != n) err(EXIT_FAILURE, "%s: inconsistent datas", realfile);
		  sampleread += samplesize;
		  status = 4; //  closing block found, start again for next sample batch
		}
		else {
		  sa[pos_in_sa] = Process_sample(p, n,  nlocis, ancestor);
		  if ((n%2) == 1) pos_in_sa++;  //store it.
		  n += 1;
		}
		break;
    } // switch
  }// end read

  /* close collection */
  sa[pos_in_sa] = dummy_sample();
  *sample_col = sa;

  if (fclose(IN) == EOF) err(EXIT_FAILURE, "%s: close failed", realfile);
  free(BUFF);
#ifndef EFBUG
  free(realfile);
  if (pos) free(pos);
#endif

  return totalsamples;
}
