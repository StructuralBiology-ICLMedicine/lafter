
/*                                                                         
 * Copyright 10/12/2017 - Dr. Christopher H. S. Aylett                     
 *                                                                         
 * This program is free software; you can redistribute it and/or modify    
 * it under the terms of version 3 of the GNU General Public License as    
 * published by the Free Software Foundation.                              
 *                                                                         
 * This program is distributed in the hope that it will be useful,         
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           
 * GNU General Public License for more details - YOU HAVE BEEN WARNED!     
 *                                                                         
 * Program: LAFTER V1.1                                                    
 *                                                                         
 * Authors: Chris Aylett                                                   
 *          Colin Palmer                                                   
 *                                                                         
 */

// Inclusions
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdint.h>
#include <pthread.h>
#include <unistd.h>
#include <complex.h>
#include <fftw3.h>

// Noise and signal power thread arguments structure
typedef struct{
  double          *in;
  double         *out;
  int32_t *neighbours;
  int32_t        size;
  int32_t        step;
  int32_t      thread;
} soft_arg;

// Filter map thread arguments structure
typedef struct{
  arguments            *args;
  fftw_complex           *in;
  fftw_complex          *out;
  double               hires;
  double                 dim;
  int32_t               size;
  int32_t          half_size;
  int32_t     full_half_size;
  int32_t           ups_size;
  int32_t      ups_half_size;
  int32_t ups_full_half_size;
  int32_t               step;
  int32_t             thread;
} ups_arg;

void soften_map_thread(soft_arg *arg);
// Soften map edges
// pthread function

void upsample_map_thread(ups_arg *arg);
// Soften map edges
// pthread function
