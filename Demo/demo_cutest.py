  ############################################################################
  #                                                                          #
  #  Demo problem: Solve a problem from the CUTEst test set.                 #
  #                                                                          #
  ############################################################################

# ----------------------------- Import libraries ----------------------------- #
import SuiteOPT
import pycutest
import numpy as np
import scipy

# -------------------- Functions to preprocess CUTEst data ------------------- #
def preprocess_constraints (cute_problem, problem):
  # Get lower bound on variable from pycutest, lo <= x
  lo = cute_problem.bl
  # Set extreme negative values to -infinity
  lo[lo == -1e+20] = np.NINF
  # Check if all values are not infinite
  if np.mean(np.isinf(lo)) != 1:
    # Store lo in problem
    problem.lo = lo
  # Option to save to file
  #np.save("lo.npy", lo)

  # Get upper bound on variable from pycutest, x <= hi
  hi = cute_problem.bu
  # Set extreme positive values to infinity
  hi[hi == 1e+20] = np.Inf
  # Check if all values are not infinite
  if np.mean(np.isinf(hi)) != 1:
    # Store hi in problem
    problem.hi = hi
  # Option to save to file
  #np.save("hi.npy", hi)

  # Get constraint matrix A for polyhedral problem by evaluating at any point
  #gs, A = cute_problem.slagjac(cute_problem.x0)
  gs, A = cute_problem.slagjac(np.zeros(cute_problem.n))
  # Check if problem has linear constraint matrix A
  if A is None:
    # No linear constrints; exit program
    return
  # Convert constraint matrix to csc format for SuiteOPT
  A = A.tocsc()
  # Store A in problem
  problem.A = A
  # Option to save to file
  #scipy.sparse.save_npz("A.npz", A)

  # Get lower bound on linear constraint from pycutest, cl <= Ax
  cl = cute_problem.cl
  # Get upper bound on linear constraint from pycutest, Ax <= cu
  cu = cute_problem.cu
  # Set extreme negative values to -infinity
  cl[cl == -1e+20] = np.NINF
  # Set extreme positive values to infinity
  cu[cu == 1e+20] = np.Inf
  # Compute offset for linear constraints
  offset = cute_problem.cons(np.zeros(cute_problem.n))
  # Subtract offset vector from cl and cu and store in bl and bu
  bl = np.zeros(cute_problem.m)
  bu = np.zeros(cute_problem.m)
  for i in range(cute_problem.m):
    bl[i] = cl[i] - offset[i]
    bu[i] = cu[i] - offset[i]
  # Check if all values are not infinite
  if np.mean(np.isinf(bl)) != 1:
    # Store bl in problem
    problem.bl = bl
  # Check if all values are not infinite
  if np.mean(np.isinf(bu)) != 1:
    # Store bu in problem
    problem.bu = bu
  # Option to save to file
  #np.save("bl.npy", bl)
  #np.save("bu.npy", bu)

  # Exit program
  return

# ------------------------- Function to solve problem ------------------------ #
def solve_problem(problem_name):
  # ------------------------ Build CUTEst Test Problem ----------------------- #
  # Clear cached problem in case there is an error with it
  pycutest.clear_cache(problem_name)
  # Build requested problem
  cute_problem = pycutest.import_problem(problem_name, efirst = True, lfirst = True, nvfirst = False)

  # ---------------------- Initialize problem parameters --------------------- #
  parm = dict(PrintStat=1, PrintLevel=0, PrintParm=1)
  
  # ------------------------- Initialize problem class ----------------------- #
  problem = SuiteOPT.problem(cute_problem.n, cute_problem.m)
  
  # ------------ Fill in problem parameters based on problem type ------------ #
  problem.x  = cute_problem.x0  # Initial guess for solution
  # Polyhedral constraints, lo <= x <= hi and bl <= Ax <= bu
  preprocess_constraints (cute_problem, problem) 
  # Hessian of cost function
  H = cute_problem.ihess(np.zeros(problem.ncol))
  H = scipy.sparse.csc_matrix(H)
  # Option to save to file
  #scipy.sparse.save_npz("H.npz", H)
  # Extract sparse matrix components for hprod
  Hx = H.data    # Hessian matrix values
  Hp = H.indptr  # Hessian matrix column pointers
  Hi = H.indices # Hessian matrix row indices
  
  # ---------------------- Wrappers for cutest functions --------------------- #
  def objective(x):
    return cute_problem.obj(x)
  
  def gradient(x):
    f, g = cute_problem.obj(x, gradient=True)
    return g
  
  def valgrad(x):
    return cute_problem.obj(x, gradient=True)
  
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

  def cg_hprod(x):
    return H.dot(x)

  # ------------ Fill in problem functions based on problem type ------------- #
  # First check for constraint types that PASA cannot handle
  if problem_name in (pycutest.find_problems(constraints='Q')+
                      pycutest.find_problems(constraints='O')):
    # Problem has constraints that are nonlinear
    print("Warning: Given problem %s has nonlinear constraints\n" %problem_name)
    #print("EXITING PROGRAM\n\n") 
    #return
    # Set objective, gradient, and objgrad
    problem.objective = objective
    problem.gradient = gradient
    problem.objgrad = valgrad
  
  # Now check for specific problem types
  # Unconstrained problems
  if problem_name in pycutest.find_problems(constraints='U'):
    # Given problem is unconstrained
    problem.lo = None
    problem.hi = None
    if problem_name in pycutest.find_problems(objective='L'):
      # --- Linear programs --- #
      # Set c, lo, hi, bl, bu, and A
      f, g = cute_problem.obj(problem.x, gradient=True)
      # Linear term in cost function
      problem.c = g  
      # Option to save to file
      #np.save("c_dfl.npy", problem.c)
  
    elif problem_name in pycutest.find_problems(objective='Q'):
      # --- Quadratic programs --- #
      # Set c, hprod, lo, hi, bl, bu, and A
      f, g = cute_problem.obj(problem.x, gradient=True)
      # Linear term in cost function
      problem.c = g  
      # Hprod
      problem.cg_hprod = cg_hprod
    else:
      # Set objective, gradient, and objgrad
      problem.objective = objective
      problem.gradient = gradient
      problem.objgrad = valgrad
  elif problem_name in pycutest.find_problems(objective='L'):
    # --- Linear programs --- #
    # Set c, lo, hi, bl, bu, and A
    f, g = cute_problem.obj(problem.x, gradient=True)
    # Linear term in cost function
    problem.c = g  
    # Option to save to file
    #np.save("c_dfl.npy", problem.c)
  elif problem_name in pycutest.find_problems(objective='Q'):
    # --- Quadratic programs --- #
    # Set c and hprod
    f, g = cute_problem.obj(problem.x, gradient=True)
    # Linear term in cost function
    problem.c = g  
    #np.save("c.npy", problem.c)
    # Hprod
    problem.hprod = hprod
  else:
    # --- Nonlinear programs --- #
    # Set objective, gradient, and objgrad
    problem.objective = objective
    problem.gradient = gradient
    problem.objgrad = valgrad
  
  # Fill in problem parameters (bypass any parameters absent from your problem)
  problem.parm = parm  # custom parameter values (optional)
  
  # ---------------------------- Print Problem Name -------------------------- #
  print("Problem name: %s\n" %(cute_problem.name))
  
  # ---------------------------- Print parameters ---------------------------- #
  #SuiteOPT.print_parm("all")

  # ------------------------------- Call solver ------------------------------ #
  newx = SuiteOPT.solve(problem)
  
  # ----------------------------- Print solution ----------------------------- #
  problem.print_stats()
  
  if newx is not None:
      print ("SuiteOPT solution: %s" %str(newx))

  #g = problem.gradient(np.zeros(problem.ncol))
  #print("g(0) = ", g)
  #g = problem.gradient(np.ones(problem.ncol))
  #print("g(1) = ", g)
  #g = problem.gradient(np.arange(problem.ncol))
  #print("g([n]) = ", g)

  #print("bl = ", problem.bl)
  #print("bu = ", problem.bu)
  #print("lo = ", problem.lo)
  #print("hi = ", problem.hi)

#  # ----------------- Compare problem data to C problem data ----------------- #
#  #lo.sort()
#  lo_c = np.genfromtxt('lo_c.csv', delimiter=',')
#  #lo_c.sort()
#  lo_c[lo_c == -1e+20] = np.NINF
#  #print(lo_c)
#  #print(lo)
#  #print(lo == lo_c)
#  #print(lo - lo_c)
#  print(np.mean(problem.lo != lo_c)*len(lo_c))
#
#  #hi.sort()
#  hi_c = np.genfromtxt('hi_c.csv', delimiter=',')
#  #hi_c.sort()
#  hi_c[hi_c == -1e+20] = np.NINF
#  #print(hi_c)
#  #print(hi)
#  #print(hi == hi_c)
#  #print(hi - hi_c)
#  print(np.mean(problem.hi != hi_c)*len(hi_c))
#
#  #bl.sort()
#  bl_c = np.genfromtxt('bl_c.csv', delimiter=',')
#  #bl_c.sort()
#  bl_c[bl_c == -1e+20] = np.NINF
#  #print(bl_c)
#  #print(bl)
#  #print(bl == bl_c)
#  #print(bl - bl_c)
#  print(np.mean(problem.bl != bl_c)*len(bl_c))
#  
#  #bu.sort()
#  bu_c = np.genfromtxt('bu_c.csv', delimiter=',')
#  #bu_c.sort()
#  bu_c[bu_c == 1e+20] = np.Inf
#  #print(bu_c)
#  #print(bu)
#  #print(bu == bu_c)
#  #print(bu - bu_c)
#  print(np.mean(problem.bu != bu_c)*len(bu_c))
#  
#  #offset.sort()
#  offset_c = np.genfromtxt('offset_c.csv', delimiter=',')
#  #offset_c.sort()
#  #offset_c[offset_c == -1e+20] = np.NINF
#  #offset_c[offset_c == 1e+20] = np.Inf
#  #print(offset_c)
#  #print(offset)
#  #print(offset == offset_c)
#  #print(offset - offset_c)
##  print(np.mean(offset != offset_c)*len(offset_c))
#  
#  Ax = problem.A.data
#  #Ax.sort()
#  Ax_c = np.genfromtxt('Ax_c.csv', delimiter=',')
#  #Ax_c.sort()
#  #print(Ax_c)
#  #print(Ax)
#  #print(Ax == Ax_c)
#  #print(Ax - Ax_c)
#  print("Ax = ", problem.A.data)
#  print("Ax_c = ", Ax_c)
#  print(np.mean(problem.A.data != Ax_c)*len(Ax_c))
#
#  Ap_c = np.genfromtxt('Ap_c.csv', delimiter=',', dtype=np.int32)
#  Ai_c = np.genfromtxt('Ai_c.csv', delimiter=',', dtype=np.int32)

  # --------------- Compare problem data to python problem data -------------- #
#  #lo.sort()
#  lo_c = np.genfromtxt('lo_py.csv', delimiter=',')
#  #lo_c.sort()
#  lo_c[lo_c == -1e+20] = np.NINF
#  #print(lo_c)
#  #print(lo)
#  #print(lo == lo_c)
#  #print(lo - lo_c)
#  print(np.mean(problem.lo != lo_c)*len(lo_c))
#
#  #hi.sort()
#  hi_c = np.genfromtxt('hi_py.csv', delimiter=',')
#  #hi_c.sort()
#  hi_c[hi_c == -1e+20] = np.NINF
#  #print(hi_c)
#  #print(hi)
#  #print(hi == hi_c)
#  #print(hi - hi_c)
#  print(np.mean(problem.hi != hi_c)*len(hi_c))
#
#  #bl.sort()
#  bl_c = np.genfromtxt('bl_py.csv', delimiter=',')
#  #bl_c.sort()
#  bl_c[bl_c == -1e+20] = np.NINF
#  #print(bl_c)
#  #print(bl)
#  #print(bl == bl_c)
#  #print(bl - bl_c)
#  print(np.mean(problem.bl != bl_c)*len(bl_c))
#  for i in np.where((problem.bl != bl_c) == True)[0]:
#    print("bl[%i] = %lg\nbl_c[%i] = %lg\n" %(i, problem.bl[i], i, bl_c[i]))
#    print("bl[%i] - bl_c[%i] = %lg\n\n" %(i, i, problem.bl[i] - bl_c[i]))
#  
#  #bu.sort()
#  bu_c = np.genfromtxt('bu_py.csv', delimiter=',')
#  #bu_c.sort()
#  bu_c[bu_c == 1e+20] = np.Inf
#  #print(bu_c)
#  #print(bu)
#  #print(bu == bu_c)
#  #print(bu - bu_c)
#  print(np.mean(problem.bu != bu_c)*len(bu_c))
#  for i in np.where((problem.bu != bu_c) == True)[0]:
#    print("bu[%i] = %lg\nbu_c[%i] = %lg\n" %(i, problem.bu[i], i, bu_c[i]))
#    print("bu[%i] - bu_c[%i] = %lg\n\n" %(i, i, problem.bu[i] - bu_c[i]))
#  
#  #offset.sort()
#  #offset_c = np.genfromtxt('offset_c.csv', delimiter=',')
#  #offset_c.sort()
#  #offset_c[offset_c == -1e+20] = np.NINF
#  #offset_c[offset_c == 1e+20] = np.Inf
#  #print(offset_c)
#  #print(offset)
#  #print(offset == offset_c)
#  #print(offset - offset_c)
#  #print(np.mean(offset != offset_c)*len(offset_c))
#  
#  Ax = problem.A.data
#  #Ax.sort()
#  Ax_c = np.genfromtxt('Ax_py.csv', delimiter=',')
#  #Ax_c.sort()
#  #print(Ax_c)
#  #print(Ax)
#  #print(Ax == Ax_c)
#  #print(Ax - Ax_c)
#  #print("Ax = ", problem.A.data)
#  #print("Ax_c = ", Ax_c)
#  print(np.mean(problem.A.data != Ax_c)*len(Ax_c))
#  for i in np.where((problem.A.data != Ax_c) == True)[0]:
#    print("Ax[%i] = %lg\nAx_c[%i] = %lg\n" %(i, problem.A.data[i], i, Ax_c[i]))
#    print("Ax[%i] - Ax_c[%i] = %lg\n\n" %(i, i, problem.A.data[i] - Ax_c[i]))
#
#  #Ap_c = np.genfromtxt('Ap_c.csv', delimiter=',', dtype=np.int32)
#  #Ai_c = np.genfromtxt('Ai_c.csv', delimiter=',', dtype=np.int32)
#
#  # Check ordering of returned problem information
#  print("Equations before inequalities: ", cute_problem.eq_cons_first)
#  print("Linear before nonlinear constraints: ", cute_problem.linear_cons_first)
#  print("Nonlinear variables first: ", cute_problem.nonlinear_vars_first)

#  problem.lo = lo_c
#  problem.bl = bl_c
#  problem.bu = bu_c
#  problem.A.data = Ax_c
#  problem.A.indptr = Ap_c
#  problem.A.indices = Ai_c
#
#  # -------------------------------- Call solver ------------------------------- #
#  newx = SuiteOPT.solve(problem)
#  
#  # ------------------------------ Print solution ------------------------------ #
#  problem.print_stats()
#  
#  if newx is not None:
#      print ("SuiteOPT solution: %s" %str(newx))
  

# ------------------------ Get problem name from user ------------------------ #
pname = raw_input('Enter problem name: ')
#pname = 'AGG'
#pname = 'LUKVLE16'
solve_problem(pname)
