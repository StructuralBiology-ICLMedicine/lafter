
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
#include "fouriers.h"

// Add FFT in to FFT out
void add_fft(fftw_complex *in, fftw_complex *out, int32_t full, int32_t nthreads){
  int32_t size = full * full * (full / 2 + 1), i;
  pthread_t threads[nthreads];
  add_fft_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in = in;
    arg[i].out = out;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) add_fft_thread, &arg[i])){
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
  return;
}

void add_fft_thread(add_fft_arg *arg){
  for(int32_t index = arg->thread; index < arg->size; index += arg->step){
    arg->out[index] += arg->in[index];
  }
  return;
}

// Calculate FSC over map
double calc_fsc(fftw_complex *half1, fftw_complex *half2, int32_t full, int32_t nthreads){
  int32_t size = full * full * (full / 2 + 1), i;
  pthread_t threads[nthreads];
  calc_fsc_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].half1 = half1;
    arg[i].half2 = half2;
    arg[i].numerator = 0.0;
    arg[i].denomin_2 = 0.0;
    arg[i].denomin_1 = 0.0;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) calc_fsc_thread, &arg[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  long double numerator = 0.0;
  long double denomin_1 = 0.0;
  long double denomin_2 = 0.0;
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
    numerator += arg[i].numerator;
    denomin_1 += arg[i].denomin_1;
    denomin_2 += arg[i].denomin_2;
  }
  return (double) (numerator / sqrtl(fabsl(denomin_1 * denomin_2)));
}

void calc_fsc_thread(calc_fsc_arg *arg){
  for(int32_t index = arg->thread; index < arg->size; index += arg->step){
    arg->numerator += creal(arg->half1[index] * conj(arg->half2[index]));
    arg->denomin_1 += creal(arg->half1[index] * conj(arg->half1[index]));
    arg->denomin_2 += creal(arg->half2[index] * conj(arg->half2[index]));
  }
  return;
}

// Apply bandpass to in and writes to out
void bandpass_filter(fftw_complex *in, fftw_complex *out, list *node, int32_t full, int32_t nthreads){
  double hires = node->res + node->stp;
  double lores = node->res;
  hires = hires * hires;
  lores = lores * lores;
  double dim = (double) full;
  int32_t size = (full / 2) + 1, i;
  int32_t full_size = full * size;
  pthread_t threads[nthreads];
  filter_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in = in;
    arg[i].out = out;
    arg[i].hires = hires;
    arg[i].lores = lores;
    arg[i].dim = dim;
    arg[i].full_size = full_size;
    arg[i].full = full;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) bandpass_filter_thread, &arg[i])){
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
  return;
}

void bandpass_filter_thread(filter_arg* arg){
  double norms, kd, jd, id;
  int32_t index;
  for(int _k = 0, k = 0; _k < arg->full; _k++, k = (_k < arg->size) ? _k : _k - arg->full){
    kd = ((double) k) / arg->dim;
    for(int _j = 0, j = 0; _j < arg->full; _j++, j = (_j < arg->size) ? _j : _j - arg->full){
      jd = ((double) j) / arg->dim;
      for(int _i = arg->thread, i = arg->thread; _i < arg->size; _i += arg->step, i = _i){
        id = ((double) i) / arg->dim;
        norms = kd * kd + jd * jd + id * id;
        index = _k * arg->full_size + _j * arg->size + _i;
        arg->out[index] = arg->in[index] * (sqrt(1.0 / (1.0 + pow((norms / arg->hires), 8.0))) - sqrt(1.0 / (1.0 + pow((norms / arg->lores), 8.0))));
      }
    }
  }
  return;
}

// Butterworth lowpass from in to out
void lowpass_filter(fftw_complex *in, fftw_complex *out, list *node, int32_t full, int32_t nthreads){
  double hires = node->res + node->stp;
  hires = hires * hires;
  double dim = (double) full;
  int32_t size = (full / 2) + 1, i;
  int32_t full_size = full * size;
  pthread_t threads[nthreads];
  filter_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in = in;
    arg[i].out = out;
    arg[i].hires = hires;
    arg[i].dim = dim;
    arg[i].full_size = full_size;
    arg[i].full = full;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) lowpass_filter_thread, &arg[i])){
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
  return;
}

void lowpass_filter_thread(filter_arg *arg){
  double norms, kd, jd, id;
  int32_t index;
  for(int _k = 0, k = 0; _k < arg->full; _k++, k = (_k < arg->size) ? _k : _k - arg->full){
    kd = ((double) k) / arg->dim;
    for(int _j = 0, j = 0; _j < arg->full; _j++, j = (_j < arg->size) ? _j : _j - arg->full){
      jd = ((double) j) / arg->dim;
      for(int _i = arg->thread, i = arg->thread; _i < arg->size; _i += arg->step, i = _i){
        id = ((double) i) / arg->dim;
        norms = kd * kd + jd * jd + id * id;
        index = _k * arg->full_size + _j * arg->size + _i;
        arg->out[index] = arg->in[index] * sqrt(1.0 / (1.0 + pow((norms / arg->hires), 8.0)));
      }
    }
  }
  return;
}
