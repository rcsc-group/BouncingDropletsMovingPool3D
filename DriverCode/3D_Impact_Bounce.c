// Author: Radu Cimpeanu
// Date: 30/09/2025

#include <stdlib.h>

#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "grid/octree.h"
#include "navier-stokes/centered.h"  // solve NS equations
#include "two-phaseDOD.h"            // prevent drop-bath coalescence
#include "tension.h"                 // include surf tension between phases
#include "vof.h"                     // solve using VoF method
#include "fractions.h"               // initially define fraction of fluid phases
#include "view.h"                    // need to make the animations
#include "tag.h"                     // helps track droplet properties
#include "draw.h"		      // visualisation helper

#include "tracer.h"

// Dimensional quantities (to be passed as arguments in main below)
double rhol = 873.0; // liquid phase density (kg/m^3)
double rhoa = 1.21;  // gas phase density (kg/m^3)

double mul = 2.0e-3; // liquid dynamic viscosity (kg/ms)
double mua = 1.8e-5; // gas dynamic viscosity(kg/ms)

double sig = 18.7e-3; // surface tension (N/m)

double R0 = 0.33e-3;  // drop radius (m)

double Pool_depth = 1.0e-3; // default pool depth (m)

double Udrop = 0.5; // default drop velocity (m/s)

double Upool = 0.0; // default drop velocity (m/s)

double Gacc = 9.81; // gravitational acceleration (m/s^2)

// Dimensionless specifications
double Re = 1.0; // Reynolds number
double Fr = 1.0; // Froude number
double We = 1.0; // Weber number

int maxlevel = 7;           // default maximum refinement level
double impact_angle = 90.0; // impact angle (oblique impacts)

double t_end = 20.0001;       // max time (ensuring the final time is exported)
double domainSize = 8.0;      // dimensionless domain size
double southPoleHeight = 1.0; // dimensionless distance of drop south pole from pool surface

scalar omega[];    // vorticity
scalar velfield[]; // velocity field norm

vector h[];
double theta0 = 90;

face vector av[];

scalar drop_tracer[], pool_tracer[];
scalar * tracers = {drop_tracer, pool_tracer};

// Output folder operations
FILE * fp_interface;
FILE * fp_stats;
FILE * fp_pressure;
FILE * fp_droplets;

/* ====== */
/* Set Up */
/* ====== */
int main(int argc, char * argv[]) {

  maxlevel = atoi(argv[1]);     // prescribed maximum resolution level
  impact_angle = atof(argv[2]); // prescribed drop impact angle
  Udrop = atof(argv[3]);        // prescribed drop impact velocity
  Upool = atof(argv[4]);        // prescribed pool horizontal velocity
  R0 = atof(argv[5]);           // prescribed drop radius
  Pool_depth = atof(argv[6]);   // prescribed pool depth
  domainSize = atof(argv[7]);   // prescribed computational domain size
  
  Re = (rhol * Udrop * R0 / mul);  
  Fr = (Udrop / sqrt(Gacc * R0));
  We = (rhol * Udrop * Udrop * R0 / sig);

  size(domainSize); 
  origin(0.0, 0.0, -(domainSize/2.0));

  a = av;
  
  mu1 = 1./Re;           // fluid viscosity
  mu2 = mu1/(mul/mua);   // gas viscosity
  
  rho1 = 1.;             // fluid density
  rho2 = 1./(rhol/rhoa); // gas density
  
  f1.sigma = 1./We;      // surface tension
  f2.sigma = 1./We;      // surface tension

  init_grid(128);

    {
      char name[200];
      sprintf(name, "loginterface.dat");
      fp_interface = fopen(name, "w");
    }

    {
      char name[200];
      sprintf(name, "logstats.dat");
      fp_stats = fopen(name, "w");
    }

    {
      char name[200];
      sprintf(name, "logpressure.dat");
      fp_pressure = fopen(name, "w");
    }
    
    {
      char name[200];
      sprintf(name, "logdroplets.dat");
      fp_droplets = fopen(name, "w");
    }

    DT = 1.0e-2; 
    NITERMIN = 1; 
    NITERMAX = 200; 
    TOLERANCE = 1e-3; 
 
    run(); 

    fclose(fp_interface);
    fclose(fp_stats);
    fclose(fp_pressure);
    fclose(fp_droplets);
}


/* =================== */
/* Boundary Conditions */
/* =================== */

// outflow at top 
u.n[top] = neumann(0);
u.t[top] = neumann(0);
p[top] = neumann(0);
pf[top] = neumann(0);

// impermeability and free slip at bottom
u.n[bottom] = dirichlet(0.);
u.t[bottom] = neumann(0.);
p[bottom] = dirichlet(0.);
pf[bottom] = dirichlet(0.);

// inflow at front
// the air velocity is also set to the pool velocity value to model a converged background flow
u.n[front] = dirichlet((Upool/Udrop));
u.t[front] = dirichlet(0.);

// outflow at back
u.n[back] = neumann(0.);
p[back] = neumann(0.);
pf[back] = neumann(0.);

// impermeability and free slip at bottom left and right
u.n[left] = dirichlet(0.);
u.t[left] = neumann(0.);

u.n[right] = dirichlet(0.);
u.t[right] = neumann(0.);

event acceleration (i++) {
  foreach_face(x)  
    av.x[] += 0.0;
  foreach_face(y)  
    av.y[] -= 1./(Fr*Fr);
  foreach_face(z)  
    av.z[] += 0.0;
}

/* ================== */
/* Initial Conditions */
/* ================== */
event init(t=0.) {

  // Double-check dimensionless groupings
  fprintf(stdout, "Reynolds number Re = %0.6f \n", Re);
  fflush(stdout);
  fprintf(stdout, "Weber number We = %0.6f \n", We);
  fflush(stdout);
  fprintf(stdout, "Froude number Fr = %0.6f \n", Fr);
  fflush(stdout);
  fprintf(stdout, "Density ratio = %0.6f \n", (rhol/rhoa));
  fflush(stdout);
  fprintf(stdout, "Viscosity ratio = %0.6f \n", (mul/mua));
  fflush(stdout);

  double dropr = 1.;

  double dropx = 0.;
  double dropy = (Pool_depth/R0) + dropr + southPoleHeight;
  double dropz = 0.;
  if (impact_angle < 90.0) // check for normal impact
  	dropz = -0.016666*impact_angle + 1.5; // general formula to maintain the south pole at the given height for non-normal impact, linearly from z = 1.25 for angle = 15 deg. to z = 0 for angle = 90 deg.

  refine((((sq(x - dropx) + sq(y - dropy) + sq(z - dropz) < sq(dropr*1.025)) && \
          ((sq(x - dropx) + sq(y - dropy) + sq(z - dropz) > sq(dropr*0.975)))))
          && level < 9);
  refine(((y > ((Pool_depth/R0) - 0.025)) && (y < ((Pool_depth/R0) + 0.025))) && level < 9);
  
  
  fraction(f1, sq(dropr) - sq(x - dropx) - sq(y - dropy) - sq(z - dropz)); // The drop
  fraction(f2, (Pool_depth/R0) - y); // The pool

  fraction(drop_tracer, sq(dropr) - sq(x - dropx) - sq(y - dropy) - sq(z - dropz)); // Drop tracer
  fraction(pool_tracer, (Pool_depth/R0) - y); // Pool tracer

  DT = 1.0e-2;
  NITERMIN = 1;
  NITERMAX = 200; 
  TOLERANCE = 1e-3;

  /* Set initial velocity */
  foreach() {
    // if within the drop (and slightly outside for smoothness) set relevant velocities
    if (sq(x - dropx) + sq(y - dropy) + sq(z - dropz) < 1.05*sq(dropr)) {
	u.x[] =  0.0;      
	u.y[] = -sin(M_PI*(impact_angle/180.0));
	u.z[] = -cos(M_PI*(impact_angle/180.0)); 
    } else {
	u.x[] =  0.0;      
	u.y[] =  0.0;
	u.z[] =  (Upool/Udrop); 
    }
  }
}


/* ========== */
/* Adapt Grid */
/* ========== */
event adapt (i++) {
  // Compute vorticity
  vorticity(u, omega);

  foreach()
    velfield[] = pow(u.x[]*u.x[]+u.y[]*u.y[]+u.z[]*u.z[], 0.5);

  // Refine based on interface locations and velocity changes
  adapt_wavelet((scalar*){f1,f2,u}, (double[]){1e-6, 1e-6, 1e-2, 1e-2, 1e-2}, maxlevel, maxlevel-4);

  // Unrefine far away from the impact site
  unrefine( (sq(x) + sq(z) > sq(5.0)) && level > (maxlevel-2));

}

event loginterface (t = 0.0; t += 0.1) {

    scalar posX[],posY[],posZ[];
    position (f1, posX, {1,0,0}); // (1,0,0) indicates the unit vector in the x-direction
    position (f1, posY, {0,1,0}); // (0,1,0) indicates the unit vector in the y-direction
    position (f1, posZ, {0,0,1}); // (0,1,0) indicates the unit vector in the y-direction

    fprintf(fp_interface, "%i %g %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f \n", i, t, statsf(f1).sum, statsf(posX).min, statsf(posX).max, statsf(posY).min, statsf(posY).max, statsf(posZ).min, statsf(posZ).max);
    fflush(fp_interface);
}

event logstats (t = 0.0; t += 0.1; t <= t_end) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}

// Fluid volume metrics
event droplets (t += 0.1)
{
  scalar m[];
  foreach()
    m[] = f1[] > 1e-3;
  int n = tag (m);

  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f1[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f1[]*p.x;
    }

  #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  for (int j = 0; j < n; j++)
    fprintf (fp_droplets, "%d %g %d %g %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j]);
  fflush (fp_droplets);
}

event small_droplet_removal (i++) {
    // Removes any small droplets that have formed, that are smaller than a specific size
    remove_droplets(f1, 8);       // Removes droplets of diameter 8 cells or less
    remove_droplets(f1, 8, true); // Removes bubbles of diameter 8 cells or less
    remove_droplets(f2, 8);       // Removes droplets of diameter 8 cells or less
    remove_droplets(f2, 8, true); // Removes bubbles of diameter 8 cells or less
}


event saveInterfaces (t += 0.1) {

    char nameInterfacesDrop[200];
    char nameInterfacesPool[200];

    sprintf(nameInterfacesDrop,"Interfaces/interfaceDrop-%0.1f.dat",t);
    sprintf(nameInterfacesPool,"Interfaces/interfacePool-%0.1f.dat",t);

    FILE * fpDrop = fopen(nameInterfacesDrop, "w");
    FILE * fpPool = fopen(nameInterfacesPool, "w");
	
    output_facets (f1, fpDrop);
    output_facets (f2, fpPool);

    fclose(fpDrop);
    fclose(fpPool);
}



/* ========== */
/* Animations */
/* ========== */
// Declare a few helper storage variables first
scalar liquidregions[], liquidregionsequal[], velnorm[], gaspressure[];

event movies (t += 0.05)
{

  char timestring[100];
  scalar omega[];
  vorticity (u, omega);
  
    foreach(){	
        liquidregionsequal[] = f1[] + f2[];
  	liquidregions[] = 1.0 - f1[] - 0.5*f2[];
  	velnorm[] = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
	gaspressure[] = (1.0 - f1[] - f2[])*p[];	
  }

  fprintf(fp_pressure, "%g %g %g %g %g %g %g\n", t, normf(p).avg, normf(p).rms, normf(p).max, normf(gaspressure).avg, normf(gaspressure).rms, normf(gaspressure).max);
  fflush(fp_pressure);
  
  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("u.x", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Vel_Ux.mp4");

  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("u.y", map = cool_warm, n = {1,0,0}, min = -1.0, max = 0.5);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Vel_Uy.mp4");
  
  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("u.z", map = cool_warm, n = {1,0,0}, min = -0.1, max = 0.6);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Vel_Uz.mp4");  
  
  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("p", map = cool_warm, n = {1,0,0}, min = -0.25, max = 1.25);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Pressure.mp4");  

  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("gaspressure", map = cool_warm, n = {1,0,0}, min = -0.25, max = 0.75);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/GasPressure.mp4");  
  
  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("liquidregions", n = {1,0,0}, map = cool_warm, min = -0.5, max = 2.5);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Liquids.mp4");  
  
  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  squares ("liquidregions", n = {1,0,0}, map = cool_warm, min = -0.5, max = 2.5);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Liquids_NoGrid.mp4");  


  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Velocity.mp4");  
  
  view (fov = 22.0, camera = "left", tx = -0.01, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("omega", n = {1,0,0}, map = cool_warm, min = -0.5, max = 0.5);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Vorticity.mp4");  
  
  view (fov = 22.0, camera = "right", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  cells(n = {1,0,0});
  isosurface("f1", 0.5, fc = {0.8,0.8,1.0});
  isosurface("f2", 0.5, fc = {1.0,1.0,1.0});
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Velocity_Front_All.mp4");
  
  view (fov = 22.0, camera = "right", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  isosurface("f1", 0.5, fc = {0.8,0.8,1.0});
  isosurface("f2", 0.5, fc = {1.0,1.0,1.0});
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Velocity_Front_All_NoGrid.mp4");
  
  view (fov = 22.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  isosurface("f1", 0.5, fc = {0.8,0.8,1.0});
  isosurface("f2", 0.5, fc = {1.0,1.0,1.0});
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Velocity_Back_All.mp4");
  
  view (fov = 22.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  isosurface("f1", 0.5, fc = {0.8,0.8,1.0});
  isosurface("f2", 0.5, fc = {1.0,1.0,1.0});
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Animations/Velocity_Back_All_NoGrid.mp4");

  view (fov = 22.0, camera = "right", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  cells(n = {1,0,0});
  isosurface("f1", 0.5, fc = {0.8,0.8,1.0});
  isosurface("f2", 0.5, fc = {1.0,1.0,1.0});
  save ("Animations/Velocity_Front_All_Bottom.mp4");
  
  view (fov = 35.0, camera = "iso", tx = -0.35, ty = -0.05, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {0,0,1}, alpha = -(domainSize/2.0));
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {0,1,0});
  cells(n = {0,0,1}, alpha = -(domainSize/2.0));
  cells(n = {1,0,0});
  isosurface("f1", 0.5, fc = {1.0,1.0,1.0});
  isosurface("f2", 0.5, fc = {1.0,1.0,1.0});
  save ("Animations/Iso_All_Far.mp4");

  view (fov = 25.0, camera = "iso", tx = -0.35, ty = -0.05, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {0,0,1}, alpha = -(domainSize/2.0));
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {0,1,0});
  cells(n = {0,0,1}, alpha = -(domainSize/2.0));
  cells(n = {1,0,0});
  isosurface("f1", 0.5, fc = {1.0,1.0,1.0});
  isosurface("f2", 0.5, fc = {1.0,1.0,1.0});
  save ("Animations/Iso_All_Close.mp4");

}
