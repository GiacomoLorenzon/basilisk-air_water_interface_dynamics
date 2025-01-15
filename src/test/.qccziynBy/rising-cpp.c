@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#define BGHOSTS 2
#include "common.h"
#include "grid/quadtree.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "rising.c"
/**
# Rising bubble

A two-dimensional bubble is released in a rectangular box and raises
under the influence of buoyancy. This test case was proposed by
[Hysing et al, 2009](/src/references.bib#hysing2009) (see also [the
FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html)).

We solve the incompressible, variable-density, Navier--Stokes
equations with interfaces and surface tension. We can solve either the
axisymmetric or planar version. We can used standard or "reduced"
gravity. */
#if AXIS //true
#include "axi.h" // fixme: does not run with -catch
#endif

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#if false
#include "view.h"
#endif

#if REDUCED //false
/** Re-express gravity in two-phase flows as an interfacial force.
 with G the gravity vector and Z an optional reference level.*/
#include "reduced.h"
#endif

#ifndef LEVEL
#define LEVEL 9
#endif

#define FINAL_T 6.0001
#define LOGFILE "rising_output_9_2/log_file.txt"
#define INTERFACE "rising_output_9_2/interface.txt"
#define OUTPUTFIELD "rising_output_9_2/output_field_bolla-%2.1f.txt"
#define OUTPUTFIELDFUORI "rising_output_9_2/output_field-%2.1f.txt"
#define OUTPUTLOGCONV "rising_output_9_2/log_convergence.txt"
#define VELGRIDX "rising_output_9_2/velocity_grid_x-%2.1f.txt"
#define VELGRIDY "rising_output_9_2/velocity_grid_y-%2.1f.txt"
#define PGRID "rising_output_9_2/pressure_grid-%2.1f.txt"
#define VORTGRID "rising_output_9_2/vorticty_grid-%2.1f.txt"

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */
u.t[right] = dirichlet(0);
u.t[left] = dirichlet(0);

/**
We make sure there is no flow through the top and bottom boundary,
otherwise the compatibility condition for the Poisson equation can be
violated. */
uf.n[bottom] = 0.;
uf.n[top] = 0.;

int main()
{
  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  $256\times 64$ grid points. */
  size(2);
  init_grid(1 << LEVEL);

  /**
  Hysing et al. consider two cases (1 and 2), with the densities, dynamic
  viscosities and surface tension of fluid 1 and 2 given below. */
  rho1 = 1000., mu1 = 10.;
    
  // CASE 1
  // rho2 = 100., mu2 = 1., f.sigma = 24.5;
  // CASE2
  rho2 = 1., mu2 = 0.1, f.sigma = 1.96;
    
  // Consider lowering the Tolerance
  TOLERANCE = 1e-4;
#if REDUCED //false
  G.x = -0.98;
  Z.x = 1.;
#endif
  run();
}

// Initial conditions
event init(t = 0)
{
#if false // !MAC
    restore (file = "lid-restore.dump");
#endif
  /**
  The domain is a rectangle. We only simulate half the bubble. */
  mask(y > 0.5 ? top : none);
  /**
  The bubble is centered on (0.5,0) and has a radius of 0.25. */
  fraction(f, sq(x - 0.5) + sq(y) - sq(0.25));
}

/**
We add the acceleration of gravity. */
#if !REDUCED // true
event acceleration(i++)
{
  face vector av = a;
  foreach_face(x)
      av.x[] -= 0.98;
}
#endif

/**
A utility function to check the convergence of the multigrid
solvers. */
#if false
void mg_print(mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
  {
    printf("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
           mg.resb > 0 ? exp(log(mg.resb / mg.resa) / mg.i) : 0.,
           mg.nrelax);
  }
}

event log_conv(i++)
{
    char nameF[80] = OUTPUTLOGCONV;
    static FILE * fp10;
    
    if (i==0)
    {
        // write file
        fp10 = fopen (nameF, "w");
          mg_print(mgp);
          mg_print(mgpf);
          mg_print(mgu);
        fprintf (fp10, "\n");
        fclose(fp10);
    }
    else // than append data
    {
        // append file
        fp10 = fopen (nameF, "a");
          mg_print(mgp);
          mg_print(mgpf);
          mg_print(mgu);
        fprintf (fp10, "\n");
        fclose(fp10);
    }
}
#endif

/**
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */
event logfile(i++)
{
  // Declare and initialize the name of the file
  char nameF[80] = LOGFILE;
  static FILE * fp1;
    
  double xb = 0., vb = 0., sb = 0.;
  foreach (reduction(+
                     : xb) reduction(+
                                     : vb) reduction(+
                                                     : sb))
  {
    // dv() is the volume of the grid element
    // f = 0 bolla
    double dv = (1. - f[]) * dv();
    xb += x * dv; // posiz centro massa
    vb += u.x[] * dv; // velocita centro massa
    sb += dv; // area bolla
  }
    // At the first iteration write the file
    if (i==0)
    {
        // write file
        fp1 = fopen (nameF, "w");
        // N.B. perf: performance data. See below
        // http://basilisk.fr/src/navier-stokes/perfs.h
        fprintf (fp1, "%g %g %g %g %g %g %g",
                t, sb, xb / sb, vb / sb, dt, perf.t, perf.speed);
        fprintf(fp1, "\r\n");
        fclose(fp1);
    }
    else // than append data
    {
        // append file
        fp1 = fopen (nameF, "a");
        fprintf (fp1, "%g %g %g %g %g %g %g",
                t, sb, xb / sb, vb / sb, dt, perf.t, perf.speed);
        fprintf(fp1, "\r\n");
        fclose(fp1);
    }
}

/**
 We output the pressure field, the velocity field and the vorticty*/
event output_field_bolla(t = 0.0; t += 0.1; t <= FINAL_T)
{
    // Initialize the file-name with time-indexing
    char nameF[80];
    sprintf (nameF, OUTPUTFIELD, t);
    
    // Compute vorticty
    scalar vort[];
    vorticity (u, vort);
    
    static FILE * fp2;
    
     if (i == 0) {
       fp2 = fopen (nameF, "w");
       foreach()
       {
         fprintf (fp2, "%f %f %f %f %f %f\r\n",x,y,u.x[]*(1. - f[]) ,u.y[]*(1. - f[]), p[]*(1. - f[]), vort[]*(1. - f[]));
       }
       fclose(fp2);
     } else {
       fp2 = fopen (nameF, "a");
       foreach()
       {
         fprintf (fp2, "%f %f %f %f %f %f\r\n",x,y,u.x[]*(1. - f[]),u.y[]*(1. - f[]), p[]*(1. - f[]), vort[]*(1. - f[]));
       }
       fclose(fp2);
   }
    
    /*
     * or
           output_field({p}, fp2, linear = true);
     */
}

event output_field_fuori(t = 0.0; t += 0.1; t <= FINAL_T)
{
    // Initialize the file-name with time-indexing
    char nameF[80];
    sprintf (nameF, OUTPUTFIELDFUORI, t);
    
    // Compute vorticty
    scalar vort[];
    vorticity (u, vort);
    
    static FILE * fp2;
    
     if (i == 0) {
       fp2 = fopen (nameF, "w");
       foreach()
       {
         fprintf (fp2, "%f %f %f %f %f %f\r\n",x,y,u.x[]*(f[]) ,u.y[]*(f[]), p[]*(f[]), vort[]*(f[]));
       }
       fclose(fp2);
     } else {
       fp2 = fopen (nameF, "a");
       foreach()
       {
         fprintf (fp2, "%f %f %f %f %f %f\r\n",x,y,u.x[]*(f[]),u.y[]*(f[]), p[]*(f[]), vort[]*(f[]));
       }
       fclose(fp2);
   }
    
    /*
     * or
           output_field({p}, fp2, linear = true);
     */
}


#if false
/**
We generate an animation using Basilisk View. */
event movie (t += 0.1)
{
  scalar omega[];
  vorticity (u, omega);
  view (tx = -0.5);
  clear();
  draw_vof ("f");
  squares ("omega", bilinear = true,  = 10);
  box ();
  static FILE * fp3 = popen ("ppm2mp4 > movie.mp4", "w");
  save (fp3 = fp3);
}
#endif

/**
We output the shape of the bubble every $t=0.1s$ */
event interface(t = 0.0; t += 0.1; t <= FINAL_T)
{
    char nameF[80] = INTERFACE;
    static FILE * fp5;
    
    if(i==0)
    {
      fp5 = fopen (nameF, "w");
      output_facets(f, fp5);
    } else {
      fp5 = fopen (nameF, "a");
      output_facets(f, fp5);
    }
  fprintf(fp5, "\n=====================\n");
  fclose(fp5);
}

/**
If gfsview is installed on the system, we can also visualise the
simulation as it proceeds. */
#if false
event gfsview (i += 10) {
  static FILE * fp4 = popen("gfsview2D rising.gfv", "w");
  scalar vort[];
  vorticity (u, vort);
  output_gfs (fp4);
}
#endif

/**
Write datafiles. To write ".gfs" files, out_gfs (file = nameGfs, t=1); is used. For writing data files read and
 processed by view, "dump (file = nameBiew); can be used.*/
event snapshot (t = 0; t += 0.1)
{
  dump (file = "dump");
  char nameGfs[80];
  sprintf (nameGfs, "rising_intermediateGfs/snapshot-%2.1f.gfs", t);
  output_gfs (file = nameGfs, t = t);
  char nameBview[80];
  sprintf (nameBview, "rising_intermediateBview/snapshot-%2.1f", t);
  dump (file = nameBview);
}

#if true //ADAPT
/**
Wavelet method for AMR based on the maximum allowable error. In the syntax {f,u}
 defines the variables based on which refinmenet is to be done. In the second term {..., ..., ...} enlist the maximum error allowed.*/
event adapt(i++)
{
  adapt_wavelet({f, u}, (double[]){5e-4, 1e-3, 1e-3}, LEVEL);
}
#endif

/**
## Results

The final shape of the bubble is compared to that obtained with the
MooNMD Lagrangian solver (see [the FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble/bubble_verification.html))
at the highest resolution. We also display the shape of the
axisymmetric version of the test. The axisymmetric bubble moves much
faster.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 1.
set term push
set term @SVG size 640,320
set size ratio -1
set grid
plot [][0:0.4]'../c1g3l4s.txt' u 2:($1-0.5) w l t 'MooNMD', \
              'log' u 1:2 w l t 'Basilisk', \
              '../rising-axi/log' u 1:2 w l t 'Basilisk (axisymmetric)'
~~~

For test case 2, the mesh in Basilisk is too coarse to accurately
resolve the skirt.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 2.
plot [][0:0.4]'../c2g3l4s.txt' u 2:($1-0.5) w l t 'MooNMD', \
              '../rising2/log' u 1:2 w l t 'Basilisk'
~~~

The agreement for the bubble rise velocity with time is also good.

~~~gnuplot Rise velocity as a function of time for test case 1.
set term pop
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'../c1g3l4.txt' u 1:5 w l t 'MooNMD', \
              'out' u 1:5 w l t 'Basilisk', \
              '../rising-axi/out' u 1:5 w l t 'Basilisk (axisymmetric)'
~~~

~~~gnuplot Rise velocity as a function of time for test case 2.
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'../c2g3l4.txt' u 1:5 w l t 'MooNMD', \
              '../rising2/out' u 1:5 w l t 'Basilisk'
~~~

*/
#if MAC
event snapshot (i++)
{
  dump (file = "dump");
}
#endif

event velocity_grid_x (t += 0.1) {
  char name[80];
  sprintf (name, VELGRIDX, t);

  FILE*fp = fopen (name, "w");
    
  for (double x = 0; x <= 2.0; x += 0.02)
  {
    for (double y = 0; y <= 0.5; y += 0.02)
    {
        fprintf (fp, " %g ", interpolate (u.x, x, y));
    }
    fprintf (fp, " \n ");
  }
  fclose (fp);
}
event velocity_grid_y (t += 0.1) {
  char name[80];
  sprintf (name, VELGRIDY, t);

  FILE*fp = fopen (name, "w");
    
  for (double x = 0; x <= 2.0; x += 0.02)
  {
    for (double y = 0; y <= 0.5; y += 0.02)
    {
        fprintf (fp, " %g ", interpolate (u.y, x, y));
    }
    fprintf (fp, " \n ");
  }
  fclose (fp);
}
event pressure_grid (t += 0.1) {
  char name[80];
  sprintf (name, PGRID, t);

  FILE*fp = fopen (name, "w");
    
  for (double x = 0; x <= 2.0; x += 0.02)
  {
    for (double y = 0; y <= 0.5; y += 0.02)
    {
        fprintf (fp, " %g ", interpolate (p, x, y));
    }
    fprintf (fp, " \n ");
  }
  fclose (fp);
}
event vorticity_grid (t += 0.1) {
  char name[80];
  sprintf (name, VORTGRID, t);

  FILE*fp = fopen (name, "w");

  scalar vort[];
  vorticity (u, vort);
    
  for (double x = 0; x <= 2.0; x += 0.02)
  {
    for (double y = 0; y <= 0.5; y += 0.02)
    {
        fprintf (fp, " %g ", interpolate (vort, x, y));
    }
    fprintf (fp, " \n ");
  }
  fclose (fp);
}

#endif
