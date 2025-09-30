#!/bin/bash

# Author: Radu Cimpeanu
# Date: 30/09/2025

# Additional resolution levels, drop velocities and pool velocities can be added below

for LEVEL in 10; do
	for UDROP in 0.603; do
		for UPOOL in 0.15; do
		
			# create new output folder 
			mkdir BounceOutput-Udrop$UDROP-Upool$UPOOL-Level$LEVEL
			
			# everything from -L$BASILISK/gl onward is just for visualisation (but needs a correct installation)
			qcc -autolink -O2 -Wall -fopenmp -g -Wall -pipe -D_FORTIFY_SOURCE=2 -o BounceOutput-Udrop$UDROP-Upool$UPOOL-Level$LEVEL/3D_Impact 3D_Impact_Bounce.c -L$BASILISK/gl -lglutils -lfb_tiny -lGLU -lGLEW -lGL -lX11 -lm
			
			# enter folder
			cd BounceOutput-Udrop$UDROP-Upool$UPOOL-Level$LEVEL
			mkdir Slices
			mkdir Interfaces
			mkdir Animations

			# execute code given resources
			export OMP_NUM_THREADS=8
			# currently the parameters are:
			# 1: maxlevel 
			# 2: impingement angle
			# 3: initial drop (vertical) velocity - dimensional, m/s
			# 4: pool (horizontal) velocity - dimensional, m/s
			# 5: drop radius - dimensional, m
			# 6: pool depth - dimensional, m
			# 7: domain size non-dimensional, relative to radius
			./3D_Impact $LEVEL 90.0 $UDROP $UPOOL 0.23e-3 2.30e-3 20.0 
			
			# exit folder      
			cd ..    

		done
	done		
done

       
