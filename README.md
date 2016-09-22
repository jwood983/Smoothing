Smoothing
=========

Smoothing will smooth a set of particle positions onto a square grid, simulating a CCD image of an astronomical object. The file of particle positions are expected to be normalized such that `maxval(x) = 1`. File input and adaptive resolutions (up to 10) are accessed through the `init.nml` file, a Fortran namelist file that looks like

    &particles
     filename = 'thin_shell.txt'
     resolution = 0.5, 0.15, 0.04
    /
    
Currently, resolutions greater than 0.5 cause an error in the calculation of nearest neighbors, so keep the resolutions below that. At the same point, resolutions less than 0.004 tend to cause the program to return the original point distribution, only in a file about 400 times bigger.

Algorithm
=========

The algorithm is basically a set of loops over the points and finding the cell where the point lies within and going around the 9 boxes around it to compute where the point smooths out towards:

    do k=1,numpart
       ii = findbox(xp,xgrid)
       jj = findbox(yp,ygrid)
       do j=jj-1,jj+1
          do i=ii-1,ii+1
             ...
          enddo
       enddo
    enddo
    
The domain also extends out from -0.5 to 1.5, so the total range is a 2x2 degree image (this assumes that the input file is normalized to a maximum of 1 in either direction). The reason for this is to ensure that the lower resolution runs do not produce a segmentation fault when checking the adjacent 8 boxes in the new algorithm.

Kernels
-------

Smooth currently has four smoothing kernels:

 1. Quartic biweight
 2. Epanechnikov
 3. Gaussian
 4. Tricubic
 
Each of these are defined on the [Kernels (statistics)](http://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use) Wikipedia page. New kernels can be freely added following the format of the current selection (e.g., using the Fortran pre-processor directives). Note, though, that complex ones, such as the Gaussian, can be computationally expensive and may not necessarily be the best choice for larger datasets.

    

Compilation
===========

Smoothing uses Fortran pre-processing directives (which is pretty much the same thing as C/C++ preprocessing). In order to compile it, you need to pass the `cpp` option through it:

    gfortran -cpp gridding.f90 -O3 -xHost -o gridding
    
If you have Intel's Fortran compiler, you can ignore the need to run gcc first by using `-fpp`:

    ifort -fpp -O3 -xHost gridding.f90 -o gridding
    
I have no experience running this on Windows or Mac, but it should work on any distribution so long as there is a Fortran compiler.

Visualization
=============

The output format follows the gnuplot pm3d map input:

    x1 y1 value
    x1 y2 value
    ...
     
    x2 y1 value
    ...

Due to this, the files can get rather large (h=0.001 produces a 300 MB file). Reading this into gnuplot is simple enough:

    set xrange[-0.2:1.2]
    set yrange[-0.2:1.2]
    set pm3d map
    set palette rgbformulae -23,-28,-3
    splot <filename> u 1:2:3
    
where `<filename>` is replaced with the actual filename, which has the format `grid_<resolution>.dat`. The `rgbformulae` part heightens the contrast by forcing large values to be dark green while smaller values are white.

Disclaimer
---------

I do not take any responsibility for any damages to your computer that result from using or abusing this code.


    
    
    
