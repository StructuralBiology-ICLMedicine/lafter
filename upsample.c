
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

// Library header inclusion for linking                                  
#include "lafter.h"
#include "upsample.h"

// Soften map
void soften_map(double *in, double *out, int32_t size, int32_t nthreads){
  int32_t i, j;
  int32_t full = size * size * size;
  int32_t neighbours[6] = {1, -1, size, -size, size * size, -size * size};
  double *tmp;
  pthread_t threads[nthreads];
  soft_arg arg[nthreads];
  for (i = 0; i < nthreads; i++){
    arg[i].neighbours = neighbours;
    arg[i].size = full;
    arg[i].step = nthreads;
    arg[i].thread = i;
  }
  for (j = 0; j < 8; j++){
    // Start threads
    for (i = 0; i < nthreads; i++){
      arg[i].in = in;
      arg[i].out = out;
      if (pthread_create(&threads[i], NULL, (void*) soften_map_thread, &arg[i])){
        printf("\nThread initialisation failed!\n");
        fflush(stdout);
        exit(1);
      }
    }
    // Join threads
    for (i = 0; i < nthreads; i++){
      if (pthread_join(threads[i], NULL)){
        printf("\nThread failed during run!\n");
        fflush(stdout);
        exit(1);
      }
    }
    tmp = in;
    in = out;
    out = tmp;
  }
  return;
}

void soften_map_thread(soft_arg *arg){
  int32_t i, j;
  double mean;
  for (i = arg->thread; i < arg->size; i += arg->step){
    if (arg->in[i] != 0.0){
      arg->out[i] = arg->in[i];
      continue;
    }
    mean = 0;
    for (j = 0; j < 6; j++){
      mean += arg->in[(arg->size + i + arg->neighbours[j]) % (arg->size)];
    }
    arg->out[i] = mean / 6.0;
  }
  return;
}

// Write out upsampled mrc file combining in1/2 - sharpening specified by command line args
void write_upsampled(fftw_complex *in, double res, r_mrc *mask, char *name, arguments *args, int32_t size, int32_t nthreads){

  int32_t i;

  // Allocate memory for maps
  fftw_complex *ks = fftw_alloc_complex(size * args->ups * size * args->ups * ((size * args->ups) / 2 + 1));
  double *rs = fftw_alloc_real(size * args->ups * size * args->ups * size * args->ups);

  // Make FFTW plans
  fftw_plan fft_ks_rs = fftw_plan_dft_c2r_3d(size * args->ups, size * args->ups, size * args->ups, ks, rs, FFTW_ESTIMATE);

  // Zero fill maps
  memset(ks, 0, size * args->ups * size * args->ups * ((size * args->ups) / 2 + 1) * sizeof(fftw_complex));
  memset(rs, 0, size * args->ups * size * args->ups * size * args->ups * sizeof(double));

  // Start threads
  double hires = res * res;
  double dim = (double) size;
  pthread_t threads[nthreads];
  ups_arg arg[nthreads];
  for (i = 0; i < nthreads; i++){
    arg[i].args = args;
    arg[i].in = in;
    arg[i].out = ks;
    arg[i].hires = hires;
    arg[i].dim = dim;
    arg[i].size = size;
    arg[i].half_size = size / 2 + 1;
    arg[i].full_half_size = size * arg[i].half_size;
    arg[i].ups_size = size * args->ups;
    arg[i].ups_half_size = arg[i].ups_size / 2 + 1;
    arg[i].ups_full_half_size = arg[i].ups_size * arg[i].ups_half_size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) upsample_map_thread, &arg[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
  }
  fftw_execute(fft_ks_rs);
  write_mrc(mask, rs, name, size * args->ups);
  fftw_free(ks);
  fftw_free(rs);
  return;
}

void upsample_map_thread(ups_arg *arg){
  // Butterworth lowpass with sharpening from in to out
  double norms, kd, jd, id;
  int32_t in_i, out_i, i, j, k, jj, kk;
  for(int32_t k = 0; k < arg->size; k++){
    kk = k > arg->half_size ? k - arg->size : k;
    kd = ((double) kk) / arg->dim;
    kd = kd * kd;
    for(int32_t j = 0; j < arg->size; j++){
      jj = j > arg->half_size ? j - arg->size : j;
      jd = ((double) jj) / arg->dim;
      jd = jd * jd;
      for(int32_t i = arg->thread; i < arg->half_size; i += arg->step){
        id = ((double) i) / arg->dim;
        id = id * id;
        norms = kd + jd + id;
        out_i = ((arg->ups_size + kk) % arg->ups_size) * arg->ups_full_half_size + ((arg->ups_size + jj) % arg->ups_size) * arg->ups_half_size + i;
        in_i = k * arg->full_half_size + j * arg->half_size + i;
        arg->out[out_i] = arg->in[in_i] * (pow(2.0, (8.0 * norms * arg->args->sharp)) / sqrt(1.0 + pow((norms / arg->hires), 8.0)));
      }
    }
  }
  return;
}
