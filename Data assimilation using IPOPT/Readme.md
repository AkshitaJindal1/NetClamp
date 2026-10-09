# Biophysical Modelling and Data Assimilation Using IPOPT

## Overview

This repository contains the computational framework, developed by Prof Alain Nogaret (University of Bath), for biophysical modelling and parameter estimation of cellular electrical activity using variational data assimilation (DA) and the Interior Point OPTimiser (IPOPT).

The framework combines conductance-based models with electrophysiological recordings to estimate ionic conductances and other biophysical parameters governing membrane voltage dynamics. By minimising the difference between model-predicted and observed voltage trajectories, while enforcing the underlying model equations, the method enables the reconstruction of cellular electrical activity and inference of model parameters (Nogaret et al., 2016).

## Framework Files

The repository contains the following Python scripts and configuration files:

1. dyn_sys.lacto : Defines the lactotroph conductance-based model, including four state variables, 21 biophysical parameters, ionic currents, and differential equations (Tabak et al., 2011).

2. Bounds.txt : Specifies the lower bounds, upper bounds, and initial estimates for model states, coupling variables, and biophysical parameters.

3. problem_info.txt : Defines the number of observations, skipped data points, integration time step, and adaptive time-stepping settings.

4. ipopt.opt : Defines IPOPT solver settings, including convergence tolerance, maximum iterations, scaling, barrier strategy, and sparse linear solver selection.

5. mss_info.txt : Configures the multiple-shooting solution (MSS) procedure, including optimisation window size, window expansion, initialisation, and restart settings.

6. methods_mss.py : Defines the dynamical system, discretises its equations using Simpson's rule, and computes the Jacobian and Hessian matrices through symbolic differentiation.

7. outputfiles_mss.py : Automatically generates the C++ source files, header files, and Makefile required to compile and execute the IPOPT optimisation.

8. compileLacto.py : Main Python script that initialises the dynamical system, loads experimental data, defines optimisation constraints, and generates the C++ files required for IPOPT.


## Computational Workflow

1. Model definition: The lactotroph model equations and parameter bounds are loaded into Python.

2. Data assimilation setup: Experimental voltage recordings and applied current inputs are incorporated into the constrained optimisation problem.

3. Symbolic differentiation: The continuous-time model is discretised, and the required Jacobian and Hessian matrices are calculated using SymPy.

4. C++ code generation: Python automatically generates the source files required by IPOPT.

5. Parameter optimisation: The compiled C++ implementation uses IPOPT to estimate the model states and biophysical parameters by minimising the mismatch between observed and predicted voltage trajectories.



## References

Nogaret, A., Meliza, C.D., Margoliash, D. and Abarbanel, H.D.I. (2016) 'Automatic construction of predictive neuron models through large scale assimilation of electrophysiological data', Scientific Reports, 6, Article 32749, pp. 1-14. Available at: https://doi.org/10.1038/srep32749.



Tabak, J., Tomaiuolo, M., Gonzalez-Iglesias, A.E., Milescu, L.S. and Bertram, R. (2011) ‘Fast-activating voltage- and calcium-dependent potassium (BK) conductance promotes bursting in pituitary cells: A dynamic clamp study’, The Journal of Neuroscience, 31(46), pp. 16855–16863. Available at: https://doi.org/10.1523/JNEUROSCI.3235-11.2011.





