
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
#include "truncate.h"

// Updates out if in1/2 over noise - returns fractional recovery
double truncate_map(double *in1, double *in2, double *out, r_mrc *mask, list *node, arguments *args, int32_t size, int32_t nthreads){
  int32_t i, m, n, full = size * size * size;
  double cor, cur;
  // Calculate max noise
  pthread_t threads[nthreads];
  max_arg arg1[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg1[i].mask = mask;
    arg1[i].in1 = in1;
    arg1[i].in2 = in2;
    arg1[i].noise = 0.0;
    arg1[i].count = 0.0;
    arg1[i].size = full;
    arg1[i].step = nthreads;
    arg1[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) calc_max_noise_thread, &arg1[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  double count = 0.0;
  double noise = 0.0;
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
    count += arg1[i].count;
    if (noise < arg1[i].noise){
      noise = arg1[i].noise;
    }
  }
  if (node->nxt == NULL && args->ovfit != 1.0){
    if (args->ovfit == 0.0){
      args->ovfit = 1.0;
    }
    m = 0;
    do {
      n = 0;
      for (i = 0; i < full; i += 1){
        if (mask->data[i] < 0.99){
          continue;
        }
        cur = in1[i] + in2[i];
        cor = cur * cur;
        if (cor > args->ovfit * noise){
          n++;
        }
      }
      if (((double) n) / count > 0.025){
        m++;
        args->ovfit *= 2.0;
      } else {
        printf("\t ESTIMATED OVERFITTING %i-fold \n\n", m);
        break;
      }
    } while (1);
  }
  noise *= args->ovfit;
  // Pass through signal greater than noise
  ass_vox_arg arg2[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg2[i].in1 = in1;
    arg2[i].in2 = in2;
    arg2[i].out = out;
    arg2[i].noise = noise;
    arg2[i].rcv = 0.0;
    arg2[i].size = full;
    arg2[i].step = nthreads;
    arg2[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) assign_voxels_thread, &arg2[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  double rcv = 0.0;
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
    rcv += arg2[i].rcv;
  }
  return rcv / count;
}

void calc_max_noise_thread(max_arg *arg){
  int32_t i;
  double cor, cur;
  for (i = arg->thread; i < arg->size; i += arg->step){
    // Normalise input transforms first
    arg->in1[i] = arg->in1[i] / arg->size;
    arg->in2[i] = arg->in2[i] / arg->size;
    // Do not calculate statistics from voxels outside the mask
    if (arg->mask->data[i] < 0.99){
      continue;
    }
    arg->count += 1.0;
    cur = arg->in1[i] - arg->in2[i];
    cor = cur * cur;
    if (cor > arg->noise){
      arg->noise = cor;
    }
  }
  return;
}

void assign_voxels_thread(ass_vox_arg *arg){
  int32_t i;
  double cor, cur;
  for (i = arg->thread; i < arg->size; i += arg->step){
    if (fabs(arg->out[i]) > 0.0){
      arg->rcv += 1.0;
      continue;
    }
    cur = arg->in1[i] + arg->in2[i];
    cor = cur * cur;
    if (cor > arg->noise){
      arg->out[i] = cur;
      arg->rcv += 1.0;
    }
  }
  return;
}
