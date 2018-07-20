# Esh3D, a semi-analytical-numerical approach to full and half space inclusion problems.
<p align="center">
<img src="https://github.com/Chunfang/Esh3D/blob/master/Esh3D.png">
</p>
Corresponding author: Chunfang Meng ( cmeng (at) mit.edu )  
Contributors: Pradeep Sharma ( psharma (at) uh.edu ), Will Heltsley and Tabrez Ali ( tabrez.ali (at) gmail.com )

* * *

## DESCRIPTION
This is Fortran version of the Matlab code (included) for triaxial Eshebly's solution evaluation. Numerical part of the code is derived from [Defmod](https://bitbucket.org/stali/defmod/wiki/Home) and [Defmod-SWPC](https://github.com/Chunfang/defmod-swpc). Improvements over the original code are  

* runs much quicker;  
* allows arbitrary number of inclusions;   
* (optionally) applies numerical traction cancellation to obtain half space solutions.  

## Eshelby vs Okada  
* Eshelby allows full (eigen) strain tensor.
* The semi-numerical-analytical approach to half space solution allows tomographical surface. 
* Slightly more expensive if the numerical (FEM) routine is engaged. 

## INSTALL 
* 3rd party package [PETSc](https://www.mcs.anl.gov/petsc) built with [HDF5](https://support.hdfgroup.org/HDF5), e.g. configured with `--download-hdf5`.
* At src folder, run `make all` to build.
* To obtain half space solution, [Trelis/Cubit](https://cubit.sandia.gov) mesh (hex or tet) is expected by preprocessor (python) script. 

## RUN
Analytical (full space) solution  

* At example folder edit esh3d.py file "half=False", and run `./esh3d.py` to generate esh3d.inp;  
* run `../bin/esh3d -f esh3d.inp`  

Semi-analytical-numerical (half space) solution  

* Run `trelis[cubit] -nojournal -nographics esh3d.jou` to generate FE mesh esh3d.exo.  
* Have "half=True" in esh3d.py, and run `./esh3d esh3d.exo` to produce esh3d.inp with the mesh.  
* run `mpirun ../bin/esh3d -f esh3d.inp`  

Use the Python and Cubit scripts as templates to roll custom models.

## Result visualization  
A HDF5 file can be imported to [Matlab](https://www.mathworks.com/help/matlab/high-level-functions.html), [Python](https://www.h5py.org/) and [Hdfview](https://support.hdfgroup.org/products/java/hdfview)  

* ellip, inclusion parameters, see esh3d.py for explanation; 
* ocoord, evaluation locations (x,y,z km); 
* odat, evaluation data, column 1-3 displacement (ux,uy,yz m); column 4-9 stress (s11,s22,s33,s12,s23,s33 Pa); 

Additionally, if the "half" option is on  

* scoord, surface element centroid locations;
* sdat, displacement and stress at surface element centroids for full and half space (9+9 columns).

Note, the full space solution is stored even if the "half" option is on. 

* * *

## LICENSE
MIT License, see LICENSE for details.

The authors appreciate that the users cite the following papers in any publications employing this code. For feedback or contribution, please contact the corresponding author. 

* * *

## REFERENCES
* Meng, C., W. Heltsley, and D. Pollard, 2012, "Evaluation of the Eshelby Solution for the Ellipsoidal Inclusion and Heterogeneity", Computers & Geosciences, Vols. 40, pp. 40 - 48.  
* Meng, C., F. Maerten, D. Pollard, 2012. “Modeling mixed-mode fracture propagation in isotropic elastic three dimensional solid”, International Journal of Fracture, Vols. 179, pp. 45– 57.  
* Meng, C., D. Pollard, 2014, “Eshelby’s solution for ellipsoidal inhomogeneous inclusions with applications to compaction bands”, Journal of Structural Geology, Vols. 67, pp. 1 - 19.  
* Ali, T., [Defmod](https://bitbucket.org/stali/defmod) - Parallel multiphysics finite element code for modeling crustal deformation during the earthquake/rifting cycle, Arxiv.
