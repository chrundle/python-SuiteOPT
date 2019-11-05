# ----------------------------- Import libraries ----------------------------- #
import SuiteOPT
import numpy as np
# Import rosenbrock function and gradient for problem
from scipy.optimize import rosen, rosen_der

# ------------------------- Initialize problem data -------------------------- #
ncol = 5
nrow = 4 
x  = -1 * np.ones(ncol)
lo = -1 * np.ones(ncol)
hi = 2 * np.ones(ncol)
A  = np.matrix([[1, 1, 1, 1, 1],
                [1, 1, 1, 0,-1],
                [1, 0,-1,-1, 1],
                [1,-1, 0, 1, 0]])
bl = np.array([ 3,
               -1,
               -2,
               -1], dtype=np.double)
bu = np.array([ 5,
                3,
                0,
                1], dtype=np.double)
c  = np.zeros(ncol)
parm = dict(grad_tol=1e-8, bypass_checks=0)

# -------------------------- Initialize problem class ------------------------ #
problem = SuiteOPT.problem(ncol, nrow)
# Fill in problem parameters (bypass any parameters absent from your problem)
problem.x  = x   # Initial guess for solution
problem.lo = lo  # Lower bound on primal variable,   lo <= x
problem.hi = hi  # Upper bound on primal variable,   x <= hi
problem.A  = A   # Linear constraint matrix A,       bl <= Ax <= bu
problem.bl = bl  # Lower bound on linear constraint, bl <= Ax
problem.bu = bu  # Upper bound on linear constraint, Ax <= bu
problem.parm = parm  # custom parameter values (optional)
problem.objective = rosen      # objective function
problem.gradient  = rosen_der  # gradient function

# ----------------------------- Print parameters ----------------------------- #
#SuiteOPT.print_parm("all")

# -------------------------------- Call solver ------------------------------- #
newx = SuiteOPT.solve(problem)

# ------------------------------ Print solution ------------------------------ #
problem.print_stats()

if newx is not None:
    print ("SuiteOPT solution: %s" %str(newx))
