
            Welcome to the Variationally Derived Intermediates
                   Extension to GROMACS 2020 v1

                              * * * * *

This customized version is free software, distributed under GNU Lesser 
General Public License, version 2.1. Note: This is NOT the official
GROMACS software package. Instead, it is a derived version. 

Installation is identical to the one of GROMACS 2020, so for detailed 
instructions please refer to 
http://manual.gromacs.org/documentation/2020/install-guide/index.html

Note: This code package has been created in 2020. Only recently it has been moved here and was initially hosted on https://gitlab.gwdg.de/martin.reinhardt. 

The installation steps after cloning the repository are
```bash
cd gromacs-vi-extension
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
make
make check
```

Simulations can then be conducted with the 
gromacs-vi-extension/build/bin/gmx executable 

When using this code please refer to the following publication

* Martin Reinhardt and Helmut Grubmüller
  GROMACS Implementation of Non-Pairwise Variationally Derived Intermediates
  https://arxiv.org/abs/2010.14193

or, when referring to the underlying theory

Determining Free Energy Differences Through Variationally Dervied Intermediates
M.Reinhardt, H. Grubmüller
Journal of Chemical Theory and Computation, 16 (6) 2020 pp. 3504-3512
DOI: 10.1021/acs.jctc.0c00106

The structure, topology and input files used for the example cases 
provided in the GROMACS implementation publication are available in the
example_cases folder

For the underlying GROMACS package, please kindly refer to:

* GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers
  M. J. Abraham, T. Murtola, R. Schulz, S. Páll, J. C. Smith, B. Hess, E. Lindahl,
  SoftwareX, 1, (2015), 19-25
  DOI: https://doi.org/10.1016/j.softx.2015.06.001
