"""
Example file for the generation of an executable
for the determination of the dynamical state and 
parameters of a Hodgkin Huxley system.

Using PyDSI, a set of coupled differential equations
are transformed to discrete-time mappings, which are
imposed as (strong) constraints on a cost function 
which defines a misfit between data and a model of 
the system.

Now using modified version by SAW
implementing Multiple Shooting Solution (MSS)

"""

from methods_mss import DynamicalSystem
from outputfiles_mss import BuildIpoptCpp 

# Initialize dynamical system
HH = DynamicalSystem()

# Define system via text file, indicating file path, state space dimension 
# and number of system parameters
HH.define_system('dyn_sys.txt',4,21)

# Give indices of state variables to be observed and indices of 
# state vars to be coupled with control term.
# indicies of the state variables that are being observed (in this case the voltage is being measured from the trace)
# indicies of the state variables that are being coupled with control term (uu)
meas_idx = [0] # [0,1,2,...] corresponds to V,m,h,n
ctrl_idx = [0] 
HH.set_observations(meas_idx,ctrl_idx,('Vobs','../data/V.txt','uu'))
#HH.set_initial('../data/initdata.dat')
#the initdata option has been moved to the mss_info.txt file

# Discretize continuous-time equations of motion with hermite polynomial
# interpolation as an additional constraint.
HH.simpsonify(hermite=True)

# Define path of file containing variable / parameter bounds.
HH.set_bounds(path='bounds.txt')

HH.set_forcing(('Iinj','../data/Iapp.txt'))

# Calculate the constraint Jacobian and Lagrangian Hessian 
# to be provided to IPOPT.
HH.get_ipopt_hessian()

# Build set of C++ files to be compiled.
HH_build = BuildIpoptCpp(HH, 'LactoTest')
