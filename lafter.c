
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

// Main algorithm function
int main(int argc, char **argv){

  // Get arguments
  arguments *args = parse_args(argc, argv);
  int32_t nthread = get_num_jobs();
  
  // Read MRC inputs
  r_mrc *vol1 = read_mrc(args->vol1);
  r_mrc *vol2 = read_mrc(args->vol2);
  r_mrc *mask;
  if (args->rad == 0.0){
    mask = read_mrc(args->mask);
  } else {
    mask = make_msk(vol1, args, nthread);
  }
  double apix = vol1->length_xyz[0] / (float) vol1->n_xyz[0];

  // Check map sizes and CPUs
  int32_t xyz = mask->n_crs[0];
  if (xyz != vol1->n_crs[0] || xyz != vol2->n_crs[0]){
    printf("\n\t MAPS MUST BE THE SAME SIZE! \n");
    return 1;
  }
  size_t r_st = xyz * xyz * xyz * sizeof(double);
  size_t k_st = xyz * xyz * ((xyz / 2) + 1) * sizeof(fftw_complex);

  // FFTW set-up
  printf("\n\t Setting up threads and maps\n");
  printf("\n\t Using %i threads. If you want to override this, set the OMP_NUM_THREADS environment variable.\n", nthread);
  fftw_init_threads();
  fftw_plan_with_nthreads(nthread);

  // Allocate memory for maps
  double *ri1 = fftw_malloc(r_st);
  double *ri2 = fftw_malloc(r_st);
  double *ro1 = fftw_malloc(r_st);
  double *ro2 = fftw_malloc(r_st);
  fftw_complex *ki1 = fftw_malloc(k_st);
  fftw_complex *ki2 = fftw_malloc(k_st);
  fftw_complex *ko1 = fftw_malloc(k_st);
  fftw_complex *ko2 = fftw_malloc(k_st);
  
  // Make FFTW plans
  printf("\n\t FFTW doing its thing - ");
  fflush(stdout);
  fftw_plan fft_ro1_ki1 = fftw_plan_dft_r2c_3d(xyz, xyz, xyz, ro1, ki1, FFTW_MEASURE);
  printf("#");
  fflush(stdout);
  fftw_plan fft_ro2_ki2 = fftw_plan_dft_r2c_3d(xyz, xyz, xyz, ro2, ki2, FFTW_ESTIMATE);
  printf("#");
  fflush(stdout);
  fftw_plan fft_ko1_ri1 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ko1, ri1, FFTW_MEASURE);
  printf("#");
  fflush(stdout);
  fftw_plan fft_ko2_ri2 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ko2, ri2, FFTW_ESTIMATE);
  printf("#\n");
  fflush(stdout);

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);
  memset(ri1, 0, r_st);
  memset(ri2, 0, r_st);
  memset(ko1, 0, k_st);
  memset(ko2, 0, k_st);
  memset(ki1, 0, k_st);
  memset(ki2, 0, k_st);

  // Copy data into place
  add_map(vol1, ro1, nthread);
  add_map(vol2, ro2, nthread);

  // Apply masks in situ
  apply_mask(mask, ro1, nthread);
  apply_mask(mask, ro2, nthread);

  // Execute forward transform
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // Zero centre
  ki1[0] = 0.0 + 0.0I;
  ki2[0] = 0.0 + 0.0I;

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);

  // Initialise list
  list head;
  head.res = 1e-10;
  head.stp = 0.025;
  head.prv = NULL;
  head.nxt = NULL;
  head.crf = 0.00;
  head.fsc = 1.00;

  list *tail = &head;

  // Loop parameterisation
  double max_res = 0.45;
  double mean_p  = 1.00;

  // Noise suppression loop
  printf("\n\t Suppressing noise -- Pass 1 \n");
  printf("\n\t # Resolution is reported in Ångströms [Å] everywhere it is quoted ");
  printf("\n\t # MeanProb records the estimated probability voxels are not noise ");
  printf("\n\t # FSC indicates the Fourier Shell Correlation between half sets -\n\n");
  fflush(stdout);
  do {

    bandpass_filter(ki1, ko1, tail, xyz, nthread);
    bandpass_filter(ki2, ko2, tail, xyz, nthread);

    tail->fsc = calc_fsc(ko1, ko2, xyz, nthread);
    tail->crf = sqrt(fabs((2.0 * tail->fsc) / (1.0 + tail->fsc)));

    fftw_execute(fft_ko1_ri1);
    fftw_execute(fft_ko2_ri2);

    mean_p = suppress_noise(ri1, ri2, ro1, ro2, mask, tail, xyz, nthread);

    printf("\t Resolution = %12.6f | MeanProb = %12.6f | FSC = %12.6f\n", apix / (tail->res + tail->stp), mean_p, tail->fsc);
    fflush(stdout);

#ifdef DEBUG
    if (args->cut == -1.0){
      mean_p = 0.333;
    }
#endif

    if ((tail->res + tail->stp) > max_res || tail->fsc < args->cut){
      if (args->ovfit != 1.0 && tail->stp < 0.0625){
        max_res = tail->res;
        tail = end_list(tail);
        continue;
      } else {
        max_res = tail->res;
        break;
      }
    }

    tail = extend_list(tail, mean_p);

  } while (1);

  // ri1, ri2, ko2 are not needed for a while
  fftw_free(ri1);
  fftw_free(ri2);
  fftw_free(ko2);

  // Apply masks in situ
  apply_mask(mask, ro1, nthread);
  apply_mask(mask, ro2, nthread);

  // Back-transform noise-suppressed maps
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // ro1, ro2 are not needed for a while
  fftw_free(ro1);
  fftw_free(ro2);

  // Zero centre
  ki1[0] = 0.0 + 0.0I;
  ki2[0] = 0.0 + 0.0I;

  // Output result
  printf("\n\t Writing noise suppressed MRC file\n");
  fflush(stdout);
  char *name1 = "noise_suppressed.mrc";
  memset(ko1, 0, k_st);
  add_fft(ki1, ko1, xyz, nthread);
  add_fft(ki2, ko1, xyz, nthread);
  write_upsampled(ko1, max_res, mask, name1, args, xyz, nthread);

  // ro1 is needed again now, along with its FFT plan
  ro1 = fftw_malloc(r_st);
  fft_ro1_ki1 = fftw_plan_dft_r2c_3d(xyz, xyz, xyz, ro1, ki1, FFTW_ESTIMATE);
  memset(ro1, 0, r_st);

  // ri1, ri2, ko2 are needed again now, along with their FFT plans
  ri1 = fftw_malloc(r_st);
  ri2 = fftw_malloc(r_st);
  ko2 = fftw_malloc(k_st);
  fft_ko1_ri1 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ko1, ri1, FFTW_ESTIMATE);
  fft_ko2_ri2 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ko2, ri2, FFTW_ESTIMATE);

  // Truncate by SNR
  printf("\n\t De-noising volume -- Pass 2 \n");
  printf("\n\t # Recovery indicates the fraction of the mask recovered by the current resolution");
  printf("\n\t # This should reach at least 1.0 but will preferably end up considerably higher -\n\n");
  fflush(stdout);
  do {

    lowpass_filter(ki1, ko1, tail, xyz, nthread);
    lowpass_filter(ki2, ko2, tail, xyz, nthread);

    fftw_execute(fft_ko1_ri1);
    fftw_execute(fft_ko2_ri2);

    mean_p = truncate_map(ri1, ri2, ro1, mask, tail, args, xyz, nthread);

    printf("\t Resolution = %12.6f | Recovery = %12.6f\n", apix / (tail->res + tail->stp), mean_p);
    fflush(stdout);

    if (tail->prv == NULL){
      break;
    } else{
      tail = tail->prv;
    }

  } while (1);

  // ri1, ri2 are not needed again
  fftw_free(ri1);
  fftw_free(ri2);

  // ko1, ko2 are not needed for a while
  fftw_free(ko1);
  fftw_free(ko2);

  // ro2 is needed again now, along with its FFT plan
  ro2 = fftw_malloc(r_st);
  fft_ro2_ki2 = fftw_plan_dft_r2c_3d(xyz, xyz, xyz, ro2, ki2, FFTW_ESTIMATE);
  memset(ro2, 0, r_st);

  // Soften map ro1
  soften_map(ro1, ro2, xyz, nthread);

  // Sum half-maps and compare to lafter
  memset(ro2, 0, r_st);
  add_map(vol1, ro2, nthread);
  add_map(vol2, ro2, nthread);

  // Apply masks in situ
  apply_mask(mask, ro1, nthread);
  apply_mask(mask, ro2, nthread);

  // Forward transform
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // ro2 is not needed again
  fftw_free(ro2);

  // Output final volume
  printf("\n\t Outputing noise truncated MRC file\n");
  fflush(stdout);
  char *name2 = "LAFTER_filtered.mrc";
  write_upsampled(ki1, max_res, mask, name2, args, xyz, nthread);

  // ro1 is not needed again
  fftw_free(ro1);

  // ko1, ko2 are needed again now
  ko1 = fftw_malloc(k_st);
  ko2 = fftw_malloc(k_st);

  // Output quality curves
  int n;
  double diff, dist = 0.0, rmsd = 0.0;
  printf("\n\t Comparing Half set Cref to the LAFTER-Sum cross-FSC \n");
  printf("\n\t # These curves should be very close until cut-off: if they differ greatly treat the results with care - \n");
  do {

    bandpass_filter(ki1, ko1, tail, xyz, nthread);
    bandpass_filter(ki2, ko2, tail, xyz, nthread);

    tail->fsc = calc_fsc(ko1, ko2, xyz, nthread);

    printf("\n\t Resolution = %12.6f | Cref = %12.6f - ", apix / (tail->res + tail->stp), tail->crf);
    for (n = 0; n < (int)(50 * tail->crf); n++){
      printf(">");
      fflush(stdout);
    }

    printf("\n\t                           | xFSC = %12.6f - ", tail->fsc);
    for (n = 0; n < (int)(50 * tail->fsc); n++){
      printf("#");
      fflush(stdout);
    }

    diff = (tail->crf - tail->fsc);
    rmsd += diff * diff * tail->stp;
    dist += tail->stp;

    if ((tail->res + tail->stp) > max_res){
      break;
    }

    printf("\n\t                           |");
    tail = tail->nxt;

  } while (1);

  // Report RMSD between Cref and xFSC
  printf("\n\n\n\t RMSD between the half-set Cref and LAFTER-Sum cross-FSC  -  %3.6f", sqrt(rmsd / dist));
  fflush(stdout);

  // Over and out...
  printf("\n\n\n\t ++++ ++++ That's All Folks! ++++ ++++ \n\n\n");

  return 0;
}
