## Steklov (and mixed) Helmholtz and Laplace eigenvalue problem for piewcewise smooth domains with corners
This Steklov-Helmholtz (and Laplace) code for corners includes 14 files. The 2 main files to use are listed in items 1 and 2.

1. **corner_inputs.m** : The user selects the number of quadrature points, wave number, domain (7 available) and the polynomial mesh degree.
   - *polygon_points_clicker.m* is a graphical input script to create a polygon by clicking points on as empty figure.
   - *poly_points.csv* is a sample of my attempt at making a star.
 
2. **SH_corners.m** : This script uses inputs from **corner_inputs.m** and solves the Steklov eigenvalue problems (Helmholtz and Laplace). 
An additional input can be given to view different eigenfunctions in the *eigenfunction plots* section. Simply run this script after setting inputs in **corner_inputs.m**.

**The rest of the scripts are helper scripts and should not be edited.**

3. 6 files (3 pairs), are such that there are two files of the same type for Helmholtz and Laplace equations. 
   - **stek_helm_corners.m** and **stek_lap_corners.m**, are the solvers that are called in **SH_corners.m**. In these scripts, 
     the layer potential matrices are built and the generalized eigenvalue problems are solved. 
   - The discretized potential matrices are block matrices corresponding to interaction of the various boundary pieces. 
     The diagonal blocks are interactions of fixed points on a given boundary piece with other points on the same boundary piece. 
     These blocks are constructed in the **layer_pots2.m** for the Helmholtz equation and **laplayer_pots.m** for the Laplace equation. 
     Interactions between distinct boundary pieces are constructed in the **talking_pots.m** script for the Helmholtz equation and the 
     **laptalking_pots.m** script for the Laplace equation. 

4. **lip_curve.m** : this script constructs the boundary curve taking inputs from corner_inputs.m. New curves can be added here.
   - *slantline.m* is a helper file to fit lines between 2 points.

5. **efn_in_corners2.m** : this script reconstructs eigenfunctions, given eigendensities in the plane.
6. **save_plot.m** is a script to save plot. Provide plot name as string and integer for fontsize.

The graded mesh scheme is as described in section 3.5 of *Inverse Acoustic and Electromagnetic Scattering Theory, 3rd edition* by Colton and Kress. 
This project is part of work on a manuscript in preparation titled *The spectrum of the Steklov-Helmholtz operator* by Nilima Nigam, Kshitij Patil and Weiran Sun at Simon Fraser University, Canada.  
