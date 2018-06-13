
# LAFTER [ Local Agreement Filter for Transmission EM Reconstructions ]
__KAILASH RAMLAUL & CHRISTOPHER H S AYLETT - ICL SSB - 12/06/2018__

## TLDR:
- LAFTER is a local filter for single particle TEM reconstructions
- LAFTER requires FFTW3 compiled with POSIX threads for compilation
- LAFTER is a C99 program (the compilation script is included)
- LAFTER prints all its available options to stdout when called
- LAFTER requires independent half-volumes and mask from refinement
- If no mask used give a radius (--rad #) to which noise remains
- LAFTER reports the Cref and LAFTER-sum xFSC curves as validation
- LAFTER outputs a twofold upsampled MRC volume for interpretation
- The output of LAFTER is INCOMPATIBLE with coordinate refinement - overfitting will occur if this is tried as the maps are denoised

## LAFTER Algorithm
- LAFTER requires the independent unfiltered halfmaps for agreement
	 estimation, and any and all masks used for the refinement process
	 to identify regions where the noise has been down-weighted.
- LAFTER performs two passes over the volume, the first from low to
	 high resolution for noise minimisation, and the second from high
	 to low, to filter each voxel at the point of agreement.
- Each successive resolution shell is incorporated using an adaptive
	 step size proportional to both the current resolution, and the SNR
	 between halves, minimising the step size at low resolution, where
	 the signal varies strongly, and near the resolution-limit, as the
	 noise eclipses the signal.
- The incorporation of new resolution shells terminates when the FSC
	 between the masked halfmaps reaches either 0.143 (default), or the
	 user provided FSC cut-off ( --fsc fsc_cut-off e.g. 0.5 ).
- The final density maps are then explicitly lowpass filtered at the
	 resolution cut-off to truncate the resolution at this point.
- The working density map is updated at each resolution shell at the
	 Gaussian estimate of the noise distribution, the variance of which
	 is calculated from the differences between each pair of halfmaps.
- The first pass over increasing resolution suppresses noise in each
	 halfmap by weighting each voxel in each shell by the probability
	 that it is part of the noise distribution. This can be thought of
	 as “adaptive masking”.
- During a second pass the degree of overfitting is estimated at the
	 resolution limit and the noise cut-off raised proportionally. Then
	 the two noise-weighted half-maps are combined at decreasing total
	 resolution and each voxel lowpass filtered at the maximum value of
	 the observed noise distribution to attempt to avoid the carry-over
	 of any noise into the final volume.

## Calculation of CRef and comparison to LAFTER-density FSC
- The theoretical value of CRef is calculated between the unfiltered
halfsets according to the equation of Henderson and Rosenthal:
       	   double Cref = sqrt((2 * FSC) / (1 + FSC));
- Here the FSC refers to the half-set FSC. Cref is compared to the
cross FSC (xFSC) between the summed dataset and the output of the LAFTER algorithm.
- All FSCs are calculated with masked densities in every instance.

## Installation and License
- LAFTER is a performance optimised C program using FFTW3 in Fourier transformation (Frigo and Johnson, 2005) for speed and portability
- LAFTER is open source and is made available under the GNU public license, which should be included in any package.

## Interpreting LAFTER output
- In each cycle LAFTER outputs the resolution in reciprocal voxels, Nyquist being equal to 0.5.
- On the first cycle LAFTER reports MeanProb - the mean probability that voxels within the mask represent noise for that shell and FSC - the FSC between the two masked, unfiltered halfmaps.
- On the second cycle LAFTER reports Recovery - the fraction of the voxels recovered with respect to the number within the mask.
- LAFTER performs a last pass to plot Cref and the LAFTER-sum cross- FSC (xFSC) to stdout.
- These curves should ideally mirror one another however some slight variation particularly at the resolution limit is expected.

       ++++ ++++ Good luck & happy cryo-EM-ing - Kailash and ChrisA ++++ ++++
