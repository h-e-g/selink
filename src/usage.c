/* Copyright 2012 GEH Institut Pasteur */
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
#include <stdio.h>

#include "usage.h"

const struct option long_options[] = {
    { "arp",    no_argument,       NULL, 'a' },
    { "comp",    required_argument,       NULL, 'c' },
    { "help",    no_argument,       NULL, 'h' },
    { "quiet",    no_argument,       NULL, 'q' },
    { "inter",    no_argument,       NULL, 'i' },
    { "lenw",    required_argument,       NULL, 'l' },
    { "numb",    required_argument,       NULL, 'n' },
    { "outfile",    required_argument,       NULL, 'o' },
    { "filter",    required_argument,       NULL, 'f' },
    { "pi",    no_argument,       NULL, 'p' },
    { "ihs",    no_argument,       NULL, 's' },
    { "freq",    no_argument,       NULL, 'K' },
    { "version", no_argument,       NULL, 'v' },
    { "omega", no_argument,       NULL, 'w' },
    { "best", no_argument,       NULL, 'b' },
    { "pos",    required_argument, NULL, 'P' },
    { "chr",    required_argument, NULL, 'C' },
    { "thread",    required_argument, NULL, 't' },
    { "hapfile",    no_argument, NULL, 'T' },
    { NULL, 0, NULL, 0 }
};

void Usage(char *prog, int exitval) {
  FILE *f = stderr;
  (void)fprintf(f, "usage: %s [options] file_prefix\n", prog);
  (void)fprintf(f, "  -v, --version         ... Print version number and exit.\n");
  (void)fprintf(f, "  -h, --help            ... Print this message and exit.\n");
  (void)fprintf(f, "  -q, --quiet           ... Turn OFF messages on stderr.\n");
  (void)fprintf(f, "  -c, --comp    <file>  ... Load omega vectors from <file>.\n");
  /*  (void)fprintf(f, "  -I, --inter           ... Display Fst statistics only.\n"); **/
  (void)fprintf(f, "  -i, --inter           ... Turns inter population statistics ON.\n");
  (void)fprintf(f, "  -p, --pi              ... Display piA/piD statistics.\n");
  (void)fprintf(f, "  -s, --ihs             ... Display iHS statistics.\n");
  (void)fprintf(f, "  -w, --omega           ... Display AFS (omega) based statistics.\n");
  (void)fprintf(f, "  -K, --freq            ... Display haplotipic frequency.\n");
  (void)fprintf(f, "  -P, --pos     <file>  ... Use <file> as position of SNPs filter. (default all positions)\n");
  (void)fprintf(f, "  -C, --chr     <val>   ... Chromosome number (for position file).\n"); //Choumouss
  (void)fprintf(f, "  -l, --lenw    <val>   ... Size in base of analysis window.\n");
  (void)fprintf(f, "  -n, --numb    <val>   ... Size in SNP of analysis window (default %i).\n",DEFAULT_WINLEN);
  /*** (void)fprintf(f, "  -b, --best            ... Window size with best SNP method.\n"); */
  (void)fprintf(f, "  -o, --outfile <file>  ... Use <file> as output results. (default stdout)\n");
  (void)fprintf(f, "  -f, --filter  <file>  ... Use <file> as population filter with also merged populations in one. (default all population)\n");
  (void)fprintf(f, "  -T, --hapfile         ... Use .hap haplotypes file type.\n");
  (void)fprintf(f, "  -t, --thread  <val>   ... computes <val> pops in parallel.\n");
  (void)fprintf(f, "  -a, --arp             ... Toggles .arp input file reading ON.\n");
  exit(exitval);
}
