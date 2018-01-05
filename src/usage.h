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

#ifndef __USAGE_H_
#define __USAGE_H_

#include <getopt.h>

#define DEFAULT_WINLEN 1000
#define PRINT_IHS 1
#define PRINT_DIND 2
#define PRINT_K 4

//Management of the long options
const extern struct option long_options[24];

/* usage display */
void Usage(char *prog, int exitval) __attribute__ ((noreturn));
#endif
