# ----------------------------- Problem Details ------------------------------ #
#  NAME          DFL001                                                            
#  An LP, contributed by Marc Meketon.
#  It "is a 'real-world' airline schedule planning
#  (fleet assignment) problem.  This LP was preprocessed by a modified
#  version of the KORBX(r) System preprocessor.  The problem reduced in
#  size (rows, columns, non-zeros) significantly.  The row and columns were
#  randomly sorted and renamed, and a fixed adjustment to the objective
#  function was eliminated.  The name of the problem is derived from the
#  initials of the person who created it."

#  Source:
#  The NETLIB collection of test problems.

# ----------------------------- Import libraries ----------------------------- #
import SuiteOPT
import numpy as np
import scipy

# ----------------------- Initialize problem parameters ---------------------- #
lo = np.load("Data/lo_dfl.npy")
hi = np.load("Data/hi_dfl.npy")
bl = np.load("Data/bl_dfl.npy")
bu = np.load("Data/bu_dfl.npy")
c = np.load("Data/c_dfl.npy")
A = scipy.sparse.load_npz("Data/A_dfl.npz")
nrow, ncol = A.shape 
x  = np.zeros(ncol)
parm = dict(PrintStatus=1)

# -------------------------- Initialize problem class ------------------------ #
problem = SuiteOPT.problem(ncol, nrow)
# Fill in problem parameters (bypass any parameters absent from your problem)
problem.x  = x   # Initial guess for solution
problem.lo = lo  # Lower bound on primal variable,   lo <= x
problem.hi = hi  # Upper bound on primal variable,   x <= hi
problem.A  = A   # Sparse matrix A
problem.bl = bl  # Lower bound on linear constraint, Bl <= Ax
problem.bu = bu  # Upper bound on linear constraint, Ax <= Bu
problem.c  = c   # Linear vector in quadratic cost
problem.parm = parm  # custom parameter values (optional)

# ----------------------------- Print parameters ----------------------------- #
#SuiteOPT.print_parm("all")

# -------------------------------- Call solver ------------------------------- #
newx = SuiteOPT.solve(problem)

# ------------------------------ Print solution ------------------------------ #
problem.print_stats()

if newx is not None:
    print ("SuiteOPT solution: %s" %str(newx))
