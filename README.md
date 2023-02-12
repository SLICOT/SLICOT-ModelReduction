# **SLICOT Model and Controller Reduction Toolbox**  

## About 

The `SLICOT Model and Controller Reduction Toolbox` (`SLICOT-ModelReduction`) includes [SLICOT](https://github.com/SLICOT/SLICOT-Reference/)-based MATLAB and Fortran tools for computing reduced-order linear models and controllers. The toolbox employs theoretically sound and numerically reliable and efficient techniques, including Balance & Truncate, singular perturbation approximation, balanced stochastic truncation, frequency-weighting balancing, Hankel-norm approximation, coprime factorization, etc.

The main functionalities of the toolbox include:

  *   order reduction for continuous-time and discrete-time multivariable models and controllers
  *   order reduction for stable or unstable models/controllers
  *   additive error model reduction
  *   relative error model and controller reduction
  *   frequency-weighted reduction with special stability/performance enforcing weights
  *   coprime factorization-based reduction of state feedback and observer-based controllers

The toolbox main features are:

  *  computational reliability using square-root and balancing-free accuracy enhancing
  *   high numerical efficiency, using latest algorithmic developments, structure exploiting algorithms, and dedicated linear algebra tools
  *   flexibility and easy-of-use
  *   enhanced functionality, e.g, for controller reduction
  *   standardized interfaces

The programs have been extensively tested on various test examples and are fully documented.

## Requirements

The codes have been tested with MATLAB 2015b through 2021b. To use the functions, the Control System Toolbox must be installed in MATLAB running under 64-bit Windows 7, 8, 8.1 or 10.  

## License

* See [`LICENSE`](https://github.com/SLICOT/SLICOT-ModelReduction/blob/master/LICENSE) for licensing information.

## References

Please cite `SLICOT-ModelReduction` using at least one of the following references: 

* A. Varga, Model reduction software in the SLICOT library, In _Applied and Computational Control, Signals, and Circuits, Ed. B. Datta, Vol. 2, pp. 239-282, Kluwer Academic Publishers, Boston, 2001._
* P. Benner, D. Kressner, V, Sima, and A. Varga, [The SLICOT Toolboxes - a Survey](https://github.com/SLICOT/SLICOT-Reports/blob/main/SLWN2009-1.pdf), _SLICOT Working Note 2009-1, August 2009._
* P. Benner, D. Kressner, V. Sima, and A. Varga, Die SLICOT-Toolboxen für Matlab - The SLICOT Toolboxes for Matlab (in German), _at – Automatisierungstechnik, 58 (2010)._


