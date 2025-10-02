# Bouncing droplets onto a moving pool
Direct numerical simulation code infrastructure for three-dimensional drop impact onto a moving pool, supporting collaborative work with the Harris Lab at Brown. The code complements the manuscript available on arXiv.

## Installation
* The code relies on [Basilisk](<http://basilisk.fr/>) to model the Navier-Stokes equations. See the [installation page](<http://basilisk.fr/src/INSTALL>) for instructions. 
* Full visualisation capabilities have been used in order to generate animations. These may be switched off depending on the local architecture.
* The two-phase non-coalescing fluid volume implementation by V. Sanjay available [here](https://github.com/VatsalSy/Lifting-a-sessile-drop/blob/master/CaseI/two-phaseDOD.h) has been successfully employed in this study to limit numerical artifacts during contact time.
  
## Running the code
Once the Basilisk structure is in place, the driver code here is built in order to navigate parameter sweeps in resolution level, drop velocity $U_{\textrm{drop}}$ and pool velocity $U_{\textrm{pool}}$ , with one of each values added to the run_master_example.sh for brevity. Other parameters can be varied through this shell script, with both physical and computational handles provided. 

The code can be executed by simply executing this shell script via *sh run_master_example.sh* inside a terminal. Output will then be produced within a foldering structure that consists of summary DNS execution information, mass conservation and VOF data, interface coordinates and animations, which can be used for further post-processing.

## Supplementary movies
Supplementary material referred to in the companion manuscript is made available in a [separate subfolder](https://github.com/rcsc-group/BouncingDropletsMovingPool3D/tree/main/SupplementaryMovies) (alongside descriptive captions), providing visualisations of both real-life and numerical experiments.

