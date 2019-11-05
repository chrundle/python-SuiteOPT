# ----------------------------- Problem Details ------------------------------ #
#  NAME          HIER163A
#  A two-norm fitted formulation of the problem of finding the 
#  smallest perturbation of data that fits a linear model
#  arising in large-scale tabular data protection.

#  Source:
#  J. Castro, 
#  Minimum-distance controlled perturbation methods for 
#  large-scale tabular data protection, 
#  European Journal of Operational Research 171 (2006) pp 39-52.

# ----------------------------- Import libraries ----------------------------- #
import SuiteOPT
import numpy as np
import scipy

# ----------------------- Initialize problem parameters ---------------------- #
lo = np.load("Data/lo_qp.npy")
hi = np.load("Data/hi_qp.npy")
bl = np.load("Data/bl_qp.npy")
bu = np.load("Data/bu_qp.npy")
c = np.load("Data/c_qp.npy")
A = scipy.sparse.load_npz("Data/A_qp.npz")
H = scipy.sparse.load_npz("Data/H_qp.npz")
Hx = H.data    # Hessian matrix values
Hp = H.indptr  # Hessian matrix column pointers
Hi = H.indices # Hessian matrix row indices
nrow, ncol = A.shape 
x  = np.zeros(ncol)
parm = dict(PrintStatus=1)

# ----------------------- Wrappers for cutest functions ---------------------- #
def hprod(d,ifree):
  # Initialize m to length of vector d
  m = len(d)
  # Initialize n to length of vector ifree
  n = len(ifree)
  # Initialize output array to zeros
  Hd = np.zeros(m)
  # Check if ifree array is nonempty
  if len(ifree) != 0:
    for j in range(n):
      t = d[j]
      if t != 0:
        k = ifree[j]
        q = Hp[k+1]
        for p in range(Hp[k], q):
          Hd[ Hi[p] ] += t * Hx[p]
  else:
    for j in range(n):
      t = d[j]
      if t != 0:
        q = Hp[j+1]
        for p in range(Hp[j], q):
          Hd[ Hi[p] ] += t * Hx[p]
  return Hd

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
problem.hprod = hprod  # Compute product hessian and vector
problem.parm = parm  # custom parameter values (optional)

# ----------------------------- Print parameters ----------------------------- #
#SuiteOPT.print_parm("all")

# -------------------------------- Call solver ------------------------------- #
newx = SuiteOPT.solve(problem)

# ------------------------------ Print solution ------------------------------ #
problem.print_stats()

if newx is not None:
    print ("SuiteOPT solution: %s" %str(newx))
