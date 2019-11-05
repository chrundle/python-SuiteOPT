# distutils: include_dirs = ../PASA/Include/ ../PPROJ/Include/ ../NAPHEAP/Include/ ../CGDESCENT/Include/ ../SuiteSparseX/SuiteSparseX_config/ ../SuiteSparseX/CHOLMOD/Include/ ../SuiteSparseX/CCOLAMD/Include/ ../SuiteSparseX/COLAMD/Include/ ../SuiteSparseX/CAMD/Include/ ../SuiteSparseX/AMD/Include/ ../SuiteSparseX/metis-5.1.0/include/ ../SuiteOPTconfig/

# ------------------------ SuiteOPT Module docstring ------------------------- #
"""
SuiteOPT: Optimization Software Suite
=====================================

SuiteOPT is designed to solve polyhedral constrained optimization
problems of the form

(1)        min                  f(x)  
       subject to    bl <= Ax <= bu, lo <= x <= hi

where f(x) is assumed to be a continuously differentiable function
and evaluation routines for f(x) and the gradient of f(x) are 
required. For Quadratic Programs, where f(x) = 0.5 * x^T H x + c^T x, 
the user can provide an evaluation routine to compute the product of 
the Hessian with a vector, H x,  and the linear term in the cost 
function. For Linear Programs, where f(x) = c^T x, the user need only 
provide the linear term c to evaluate f(x).

Additionally, SuiteOPT has routines designed to solve problems of
the form (1) with specific objective or constraint sets

Polyhedral projection problems 

(1a)       min            0.5 * || x - y ||^2  
       subject to    bl <= Ax <= bu, lo <= x <= hi

Napsack problems

(1b)       min             0.5 * x^T D x + c^T x  
       subject to    bl <= a^T x <= bu, lo <= x <= hi

where D is a diagonal matrix with diagonal entries d.

Unconstrained problems

(1c)   min f(x)  subject to  x in R^n

See Also
--------
problem : Class within SuiteOPT module for containing SuiteOPT 
          problem data

solve : Method within SuiteOPT for solving problems of the form
        (1), (1a), (1b), or (1c)
    
References
----------
W. W. Hager and H. Zhang, An Active Set Algorithm for Nonlinear
Optimization with Polyhedral Constraints, August, 2016, 
http://people.clas.ufl.edu/hager/papers/CG/pasa.pdf

W. W. Hager and H. Zhang, Sparse Techniques for Polyhedral Projection,
December 29, 2014, http://people.clas.ufl.edu/hager/

T. A. Davis, W. W. Hager, and J. T. Hungerford, An Efficient
Hybrid Algorithm for the Separable Convex Quadratic Knapsack Problem,
http://people.clas.ufl.edu/hager/

W. W. Hager and H. Zhang, A new conjugate gradient method
with guaranteed descent and an efficient line search,
SIAM Journal on Optimization, 16 (2005), 170-192.

W. W. Hager and H. Zhang, Algorithm 851: CG_DESCENT,
A conjugate gradient method with guaranteed descent,
ACM Transactions on Mathematical Software, 32 (2006), 113-137.

W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58.

W. W. Hager and H. Zhang, Limited memory conjugate gradients,
SIAM Journal on Optimization, 23 (2013), 2150-2168.

Example
-------
For simplicity, a short demo for setting up and solving a problem 
of the form (1) is provided below.

import SuiteOPT
import numpy as np
# Import rosenbrock function and gradient for problem
from scipy.optimize import rosen, rosen_der

# Initialize problem data
ncol = 5
nrow = 4 
x  = -1 * np.ones(ncol)
lo = -1 * np.ones(ncol)
hi = 2 * np.ones(ncol)
A  = np.matrix([[1, 1, 1, 1, 1],[1, 1, 1, 0,-1],
                [1, 0,-1,-1, 1],[1,-1, 0, 1, 0]])
bl = np.array([3,-1,-2,-1], dtype=np.double)
bu = np.array([5, 3, 0, 1], dtype=np.double)
parm = dict(grad_tol=1e-8)

# Initialize problem class
problem = SuiteOPT.pasaProblem(ncol, nrow)

# Fill in problem data (bypass any data absent from your problem)
problem.x  = x   # Initial guess for solution
problem.lo = lo  # Lower bound on primal var,     lo <= x
problem.hi = hi  # Upper bound on primal var,     x <= hi
problem.A  = A   # Linear constr. matrix A,       bl <= Ax <= bu
problem.bl = bl  # Lower bound on linear constr., bl <= Ax
problem.bu = bu  # Upper bound on linear constr., Ax <= bu
problem.objective = rosen      # objective function
problem.gradient  = rosen_der  # gradient function
# Custom parameter values (optional)
problem.parm = parm        

# Call SuiteOPT solver and store solution in newx
newx = SuiteOPT.solve(problem)

"""

# ------------------------- Import python libraries -------------------------- #
# Import libraries for compatibility with python 3
from __future__ import print_function

# Cimport definitions from cSuiteOPT.pxd file
cimport cSuiteOPT
# Cimport numpy for array operations in C
cimport numpy as np

# Import numpy for array operations in python
import numpy as np
# Import scipy csc_matrix for sparse matrix format
from scipy.sparse import csc_matrix
# Import time for getting walltime 
import time
# Import sys for printing messages 
import sys
# Import inspect for counting number of function arguments (py 2/3 handling)
try:
    from inspect import getfullargspec as _get_args
except ImportError:
    from inspect import getargspec as _get_args

# ------------------ Import C macros ------------------ #
cdef extern from "math.h":
    float INFINITY

# ----------------------- Callback function prototypes ----------------------- #
# Define C prototype for objective function
ctypedef void (*c_objective) (double *val, double *x, int n)
# Define C prototype for gradient
ctypedef void (*c_gradient) (double *grad, double *x, int n)
# Define C prototype for objective/gradient
ctypedef void (*c_objgrad) (double *val, double *grad, double *x, int n)
# Define C prototype for hprod
ctypedef void (*c_hprod) (double *Qx, double *x, int *nzIndices, int n, int nnz)
# Define C prototype for cg_hprod
ctypedef void (*c_cg_hprod) (double *Qx, double *x, int n)

# ------------------------ Callback function wrappers ------------------------ #
# ------ C wrapper for python objective function ------ #
cdef void objective_callback (double *val, double *x, int n):
    # Set global python objective function
    global py_objective
    # Initialize output for objective function
    cdef double objective_value
    # Copy x values to variable that can be passed to python 
    cdef double [:] _x   = <double [:n]>x
    # Call python objective function
    objective_value = (<object>py_objective)(_x)
    # Store value in val
    val[0] = objective_value
    return

# ------ C wrapper for python gradient ------ #
cdef void gradient_callback (double *grad, double *x, int n):
    # Set global python gradient function
    global py_gradient
    # Copy x values to variable that can be passed to python 
    cdef double [:] _x   = <double [:n]>x
    # Initialize array to store gradient values in
    cdef double [:] gradient_values
    # Call python gradient function
    gradient_values = (<object>py_gradient)(_x)
    # Copy returned gradient_values to grad
    for i in range(n):
        grad[i] = gradient_values[i]
    return

# ------ C wrapper for python objective/gradient ------ #
cdef void objgrad_callback (double *val, double *grad, double *x, int n):
    # Set global python objective and gradient functions
    global py_objgrad
    # Initialize output for objective and gradient
    cdef double objective_value
    # Initialize array to store gradient values in
    cdef double [:] gradient_values
    # Copy x values to variable that can be passed to python 
    cdef double [:] _x   = <double [:n]>x
    # Call python objective/gradient function
    objective_value, gradient_values = (<object>py_objgrad)(_x)
    # Copy values to output variables for objective and gradient
    val[0] = objective_value
    for i in range(n):
        grad[i] = gradient_values[i]
    return

# ------ C wrapper for python hprod ------ #
cdef void hprod_callback (double *Qx, double *x, int *nzIndices, int n, int nnz):
    # Set global python objective and gradient functions
    global py_hprod
    # Initialize array to store matrix product in
    cdef double [:] matrix_product
    # Copy x values to variable that can be passed to python 
    cdef double [:] _x = <double [:n]>x
    # Copy nzIndices values to variable that can be passed to python 
    cdef int [:] _nzIndices
    if nzIndices != NULL:
        _nzIndices = <int [:nnz]>nzIndices
    else:
        _nzIndices = np.arange(n, dtype=np.int32)
    # Call python objective/gradient function
    matrix_product = (<object>py_hprod)(_x, _nzIndices)
    # Copy values to output variable
    for i in range(n):
        Qx[i] = matrix_product[i]
    return

# ------ C wrapper for python cg_hprod ------ #
cdef void cg_hprod_callback (double *Qx, double *x, int n):
    # Set global python objective and gradient functions
    global py_cg_hprod
    # Initialize array to store matrix product in
    cdef double [:] matrix_product
    # Copy x values to variable that can be passed to python 
    cdef double [:] _x = <double [:n]>x
    # Call python objective/gradient function
    matrix_product = (<object>py_cg_hprod)(_x)
    # Copy values to output variable
    for i in range(n):
        Qx[i] = matrix_product[i]
    return

# ----------------------------- problem class ---------------------------- #
class problem:
    """
    Class containing all data for problem to be solved using SuiteOPT.

    The solver in SuiteOPT is designed to solve polyhedral constrained
    problems of the form
   
    (1)        min                  f(x)  
           subject to    bl <= Ax <= bu, lo <= x <= hi

    where f(x) is assumed to be a continuously differentiable function
    and evaluation routines for f(x) and the gradient of f(x) are 
    required. For Quadratic Programs, where f(x) = 0.5 * x^T H x + c^T x, 
    the user can provide an evaluation routine to compute the product of 
    the Hessian with a vector, H x,  and the linear term in the cost 
    function. For Linear Programs, where f(x) = c^T x, the user need only 
    provide the linear term c to evaluate f(x).

    Additionally, SuiteOPT has routines designed to solve problems of
    the form (1) with specific objective or constraint sets
   
    Polyhedral projection problems 

    (1a)       min            0.5 * || x - y ||^2  
           subject to    bl <= Ax <= bu, lo <= x <= hi

    Napsack problems
    
    (1b)       min             0.5 * x^T D x + c^T x  
           subject to    bl <= a^T x <= bu, lo <= x <= hi

    where D is a diagonal matrix with diagonal entries d.

    Unconstrained problems

    (1c)   min f(x)  subject to  x in R^n


    Extended Summary
    ----------------
    Setting data after initializing 'problem' class.

    For each component required by the user's problem, the user should 
    assign the corresponding value in the 'problem' class. Note that
    vectors/matrices are expected to be a numpy array/matrix. For 
    example, if the user's problem has a lower bound constraint of
    1 for each component then the user should set

    >> problem.lo = numpy.ones(ncol)

    If the user wishes to enforce the constraint lo <= x then lo is 
    required. If the user wishes to enforce the constraint x <= hi then 
    hi is required. If lo = hi for the user's problem, then lo and hi 
    must both be provided with lo = hi.

    If the user wishes to enforce the constraint bl <= Ax <= bu then A, 
    bl, and bu are all required. Note that A can be provided in dense or 
    sparse format. If the user does not wish to enforce this constraint 
    then these values should be excluded from the problem class. If 
    bl = bu in your problem, both bl and bu must be provided as the same 
    vector.

    If the user wishes to specify the number of columns to be used from  
    matrix A, then the user should provide a value for ncol. Otherwise, 
    the value of ncol will be set to the minimum length of the following 
    elements (if provided): x, lo, hi, c, y, a, d.

    Attributes
    ----------
    ncol : int 
        Dimension of primal variable (# of columns in matrix A)
    nrow : int
        Dimension of dual variable (# of rows in matrix A)
    x : numpy array
        Initial guess for primal solution
    pylambda : numpy array 
        Initial guess for dual variable (multiplier for Ax)
    A : numpy matrix or scipy.sparse.csc_matrix
        Constraint matrix (can provide in dense or sparse format)
    bl : numpy array
        Lower bounds for Ax, bl <= Ax
    bu : numpy array
        Upper bounds for Ax, Ax <= bu
    lo : numpy array
        Lower bounds for x, lo <= x
    hi : numpy array
        Upper bounds for x, x <= hi
    c : numpy array
        Linear term in objective function (if LP or QP)
    y : numpy array 
        Point to project onto polyhedron (if solving problem (1a))
    d : numpy array
        Diagonal in cost function hessian (if solving problem (1b))
    a : numpy array
        Linear constraint vector (if solving problem (1b))
    objective : function
        Function for evaluating objective function
        Should have one input argument (point to evaluate objective)
        Should have one output argument (objective value)
    gradient : function
        Function gradient of the objective function
        Should have one input argument (point to evaluate gradient)
        Should have one output argument (gradient value)
    valgrad : function
        Function for objective and its gradient
        Should have one input argument (point to evaluate obj/grad)
        Should have two output arguments (obj & grad values)
    hprod : function
        Function for evaluating hessian times vector
    cg_hprod : 
        Function for evaluating hessian times vector in 
        unconstrained problem
    parm : dict
        Provide custom solver parameters in parm dictionary
    stats : dict
        Contains problem statistics for each method after calling 
        SuiteOPT.solve(problem)
    probtype : string
        Type of problem detected by SuiteOPT
    """

    # -------- Define __init__ function -------- #
    def __init__(self, int ncol=0, int nrow=0): 
        """
        Method for initializing instance of SuiteOPT.problem class.

        Parameters
        ----------
        ncol : Dimension of primal variable (# of columns in matrix A)
        nrow : Dimension of dual variable (# of rows in matrix A)

        Note that ncol and nrow can be ommitted and they will be 
        initialized to 0.

        Returns
        -------
        An instance of the SuiteOPT.problem class
        """

        # ------ Initialize class parameters ------ #
        # ncol must be nonnegative
        if ncol > 0:
            self.ncol = ncol
        else:
            self.ncol = 0
        # nrow must be nonnegative
        if nrow > 0:
            self.nrow = nrow
        else:
            self.nrow = 0
        self.x  = None
        self.lo = None
        self.hi = None
        self.y  = None
        self.c  = None
        self.a  = None
        self.d  = None
        self.pylambda = None
        self.bl = None
        self.bu = None
        self.A  = None
        self.parm      = dict()
        self.stats     = dict()
        self.objective = None
        self.gradient  = None 
        self.objgrad   = None 
        self.hprod     = None 
        self.cg_hprod  = None 
        self.probtype  = None
        # ------ Exit init function ------ #
        return

    # -------- Define __repr__ function -------- #
    def __repr__(self): 
        """
        Method for printing instance of SuiteOPT.problem class.

        Parameters
        ----------
        None

        Returns
        -------
        String containing all attribute names and current values
        """

        # Create dict containing all attribute names and current 
        # values in instance of problem class
        attr_dict = vars(self)
        # Return all attribute names and current values in string
        prob_repr = '<SuiteOPT.problem ' 
        prob_repr += ' '.join("%s:%s" % attr for attr in attr_dict.items())
        prob_repr += '>' 
        return prob_repr

    # -------- Define __str__ function -------- #
    def __str__(self): 
        """
        Method for printing instance of SuiteOPT.problem class.

        Parameters
        ----------
        None

        Returns
        -------
        String containing all attribute names and current values
        """

        # Create dict containing all attribute names and current 
        # values in instance of problem class
        attr_dict = vars(self)
        # Return all attribute names and current values in string
        return '\n'.join("%s: %s" % attr for attr in attr_dict.items())

    # -------- Define _print_stats function -------- #
    def print_stats(self): 
      """
      Method for printing problem statistics after solving a problem.

      After solving a problem with SuiteOPT, problem statistics are 
      stored in the SuiteOPT.problem.stats member. For each solver used 
      by SuiteOPT, there is a subdictionary within stats that contains 
      statistics for that solver, either pasa, pproj, napheap, or cg. 
      This method is used to view some key statistics generated when 
      solving a problem with SuiteOPT.

      Parameters
      ----------
      None

      Returns
      -------
      None
      """

      # Check problem type
      if self.probtype == "Unconstrained":
        # Unconstrained problem (only uses CG)
        if "cg" in self.stats:
          # Set d equal to cg statistics dictionary
          d = self.stats["cg"]
          print (_colors.OKBLUE + _colors.BOLD +
               " ===================== Problem Statistics ===================="
               + "\n" + _colors.ENDC + _colors.OKBLUE +
               "  Final objective value....................: %16.6lg\n"
               "  Final error..............................: %16.6lg\n"
               "  Walltime.................................: %16.6lg\n"
               "  Objective evals in cg....................: %16i\n"
               "  Gradient evals in cg.....................: %16i\n"
               "  Iterations in cg.........................: %16i\n"
               %(d["f"], d["err"], self.stats["walltime"], d["nfunc"], 
                 d["ngrad"], d["iter"])
               +  _colors.BOLD +
               " ============================================================="
               + _colors.ENDC)
      elif self.probtype == "Constrained":
        # Check if pasa statistics generated
        if "pasa" in self.stats:
          # Set d equal to pasa statistics dictionary
          d = self.stats["pasa"]
          # Print pasa statistics
          print (_colors.OKBLUE + _colors.BOLD +
               " ===================== Problem Statistics ===================="
               + "\n" + _colors.ENDC + _colors.OKBLUE +
               "  Final objective value....................: %16.6lg\n"
               "  Final error..............................: %16.6lg\n"
               "  Walltime.................................: %16.6lg\n"
               "  Objective evals in main code.............: %16i\n"
               "  Gradient evals in main code..............: %16i\n"
               "  Iterations in gradient projection........: %16i\n"
               "  Objective evals in gradient projection...: %16i\n"
               "  Gradient evals in gradient projection....: %16i\n"
               "  Iterations in active set grad. proj......: %16i\n"
               "  Objective evals in active set grad. proj.: %16i\n"
               "  Gradient evals in active set grad. proj..: %16i\n"
               %(d["f"], d["err"], self.stats["walltime"], d["mcnf"],  
                 d["mcng"], d["gpit"], d["gpnf"], d["gpng"], d["agpit"],  
                 d["agpnf"], d["agpng"])
               +  _colors.BOLD +
               " ============================================================="
               + _colors.ENDC)
      else:
        print (_colors.OKBLUE + _colors.BOLD +
               " ===================== Problem Statistics ===================="
               + "\n" + _colors.ENDC + _colors.OKBLUE +
               "  No problem statistics generated.\n"
               +  _colors.BOLD +
               " ============================================================="
               + _colors.ENDC)
      print ("")
      # Exit functoin
      return

# ------------------------ Check user data information ----------------------- #
# Check if user provided data with given data_name. If length of provided data
# is shorter than expected, update dimension for problem (ncol or nrow) to 
# match length of shorter data and print warning message to user
def _check_data(data_name, data_memview, dim_name, exp_len_arr):
    """Function used internally by SuiteOPT to check problem data before
       passing to solver."""
    # Check if user provided data
    if data_memview is not None:
        # Set d = length of data
        d = data_memview.shape[0]
        # User provided data; check if len matches current dimension 
        if d < exp_len_arr[0]:
            # Provided data has length smaller than current dimension
            print (_colors.BOLD + _colors.WARNING + 
                   "\n\n    PASA Warning: Length mismatch\n" + 
                   _colors.ENDC + _colors.WARNING + 
                   "      User's %s has length...: %i\n"
                   "      Expected length of %s..: %i (%s)\n" 
                   %(data_name, d, data_name, exp_len_arr[0], dim_name) +
                   "    Setting %s = %i and proceeding with PASA."
                   %(dim_name, d) + _colors.ENDC)
            # Set problem dimension to data_len
            exp_len_arr[0] = d
            return 1
        elif d == exp_len_arr[0]:
            # Provided data has length equal to current dimension
            return 0
        else:
            # Provided data has length larger than current dimension
            if exp_len_arr[0] == 0:
                # Current dimension is zero, update to length of data
                print (_colors.BOLD + _colors.WARNING + 
                       "\n\n    PASA Warning: Length mismatch\n" + 
                       _colors.ENDC + _colors.WARNING + 
                       "      User's %s has length...: %i\n"
                       "      Expected length of %s..: %i (%s)\n" 
                       %(data_name, d, data_name, exp_len_arr[0], dim_name) +
                       "    Setting %s = %i and proceeding with PASA."
                       %(dim_name, d) + _colors.ENDC)
                # Set problem dimension to data_len
                exp_len_arr[0] = d
                return 1
            else:
                # Current dimension is not zero 
                print (_colors.BOLD + _colors.WARNING + 
                       "\n\n    PASA Warning: Length mismatch\n" + 
                       _colors.ENDC + _colors.WARNING + 
                       "      User's %s has length...: %i\n"
                       "      Expected length of %s..: %i (%s)\n" 
                       %(data_name, d, data_name, exp_len_arr[0], dim_name) +
                       "    Using %s = %i and proceeding with PASA."
                       %(dim_name, exp_len_arr[0]) + _colors.ENDC)
                # Set problem dimension to data_len
                return 1
    else:
        # User did not provide data; return 0 to indicate no error
        return 0

# ----------------------- Check user matrix information ---------------------- #
# Check if user provided matrix. If dimension of provided matrix is different
# than expected, return indicator that error has occurred
def _check_matrix(A, expected_nrow, expected_ncol):
    """Function used internally by SuiteOPT to check problem matrix before
       passing to solver."""
    if (type(A) is np.matrix) or (type(A) is csc_matrix): 
        # Check if A has correct dimension (nrow x ncol)
        if A.shape != (expected_nrow, expected_ncol): 
            # A has incorrect dimensions
            print (_colors.BOLD + _colors.FAIL + 
                   "\n\n    PASA Error: Constraint matrix A has incorrect "
                   "shape."
                   + _colors.ENDC)
            print (_colors.WARNING + 
                   "      User's matrix A has shape...: %s\n"
                   "      Matrix A should have shape..: (%i, %i) "
                   "(= (nrow, ncol))\n" 
                   %(str(A.shape), expected_nrow, expected_ncol) +
                   _colors.ENDC + _colors.FAIL +
                   "    Terminating PASA.\n" + _colors.ENDC)
            print (_colors.OKBLUE +
                   "    Note: If problem does not have constraint "
                   "of the form\n"
                   "                         bl <= Ax <= bu\n"
                   "          then no value for A should be provided.\n" +
                   _colors.ENDC)
            return True
        else:
            return False
    elif A is not None: 
        # Input for A of incorrect type
        # Print error message
        print (_colors.BOLD + _colors.FAIL + 
               "\n\n    PASA Error: Constraint matrix A is of incorrect type."
               + _colors.ENDC)
        print (_colors.WARNING + 
               "      User's A is of type...: %s\n"
               "      A should be a numpy or csc matrix with shape "
               "(%i, %i),\n      which is (nrow, ncol)\n" 
               %(type(A), expected_nrow, expected_ncol) +
               _colors.ENDC + _colors.FAIL +
               "    Terminating PASA.\n" + _colors.ENDC)
        return True

# ----------------------- Check user function information -------------------- #
# Check if user provided function handle. If number of input arguments is 
# different than expected or user provided non-callable function, return 
# indicator that error has occurred
def _check_func(func_handle, func_name, expected_args):
    """Function used internally by SuiteOPT to check problem functions before
       passing to solver."""
    # Check if user provided function handle 
    if func_handle is not None:
        # Check if callable function was provided for objective
        if callable(func_handle):
            # Set func_args equal to number of input arguments
            func_args = len(_get_args(func_handle).args)
            # Check if objective has expected number of input arguments
            if func_args != expected_args:
                # Print error message
                print (_colors.BOLD + _colors.FAIL + 
                       "\n\n    PASA Error: %s function has incorrect\n"
                       "    number of input arguments.\n" %(func_name) +
                       _colors.ENDC + _colors.WARNING + 
                       "      Number of input args in user %s.....: "
                       "%i\n"
                       "      Number of input args required for %s: %i\n"
                       %(func_name, func_args, func_name, expected_args) +
                       _colors.ENDC + _colors.FAIL +
                       "    Terminating PASA.\n" + _colors.ENDC)
                # Error has occurred; return True 
                return True 
            else:
                # No errors occurred; return False 
                return False 
        else: # Non callable function provided for objective
            # Print error message
            print (_colors.BOLD + _colors.FAIL + 
                   "\n\n    PASA Error: Non-callable input provided for %s\n"
                   "    function.\n" %(func_name) +
                   _colors.ENDC + _colors.WARNING + 
                   "      %s must be callable with %i input argument(s).\n" 
                   %(func_name, expected_args) +
                   _colors.ENDC + _colors.FAIL +
                   "    Terminating PASA.\n" + _colors.ENDC)
            # Error has occurred; return True 
            return True 
    else:
        # User did not provide any input for current function; return False
        return False

# ---------------------- Print parm function definition ------------------- #
def print_parm(program_name):
    """Print all parameters associated with solvers used by SuiteOPT.
    
    When using SuiteOPT to solve a problem, the user may want to provide
    custom parameter values to the solvers used by SuiteOPT. This routine
    allows the user to check the names of all parameters used within 
    SuiteOPT and their default values.

    Parameters
    ----------
    program_name : String indicating which parameters to print

    print_parm accepts a single string as an input and returns the 
    parameter names and default values to the console. There are 4 
    different solvers that may be used by SuiteOPT when solving a 
    problem: PASA, PPROJ, NAPHEAP, and CGDESCENT. The accepted strings 
    for print_parm are provided below

    'all'     - Prints parameters for all solvers used by SuiteOPT
    'pasa'    - Prints parameters for PASA
    'pproj'   - Prints parameters for PPROJ
    'napheap' - Prints parameters for NAPHEAP
    'cg'      - Prints parameters for CGDESCENT

    Returns 
    -------
    None
    """

    # Allocate memory for pasadata and initialize parameters
    pasadata = cSuiteOPT.pasa_setup()
    # Check if user requested parameters to be printed
    if program_name == "all":
        # User requested all default parameter values
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================= Default PASA parameter values "
               "================="
               + _colors.ENDC)
        # Print pasa parameter values
        cSuiteOPT.pasa_print_parm(pasadata)
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n\n================ Default PPROJ parameter values "
               "================="
               + _colors.ENDC)
        # Print pproj parameter values
        cSuiteOPT.pproj_print_parm(pasadata.ppdata)
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n\n============== Default CG_DESCENT parameter values "
               "=============="
               + _colors.ENDC)
        # Print cg parameter values
        cSuiteOPT.cg_print_parm(pasadata.cgdata)
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n\n================ Default NAPHEAP parameter values "
               "==============="
               + _colors.ENDC)
        # Print napheap parameter values
        cSuiteOPT.napheap_print_parm(pasadata.napdata)
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================================================="
               "==============="
               + _colors.ENDC)
    elif program_name == "pasa":
        # User requested default pasa parameter values
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================= Default PASA parameter values "
               "================="
               + _colors.ENDC)
        # Print pasa parameter values
        cSuiteOPT.pasa_print_parm(pasadata)
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================================================="
               "==============="
               + _colors.ENDC)
    elif program_name == "pproj":
        # User requested default pproj parameter values
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================ Default PPROJ parameter values "
               "================="
               + _colors.ENDC)
        # Print pproj parameter values
        cSuiteOPT.pproj_print_parm(pasadata.ppdata)
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================================================="
               "==============="
               + _colors.ENDC)
    elif program_name == "cg":
        # User requested default cg parameter values
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n============== Default CG_DESCENT parameter values "
               "=============="
               + _colors.ENDC)
        # Print cg parameter values
        cSuiteOPT.cg_print_parm(pasadata.cgdata)
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================================================="
               "==============="
               + _colors.ENDC)
    elif program_name == "napheap":
        # User requested default napheap parameter values
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================ Default NAPHEAP parameter values "
               "==============="
               + _colors.ENDC)
        # Print napheap parameter values
        cSuiteOPT.napheap_print_parm(pasadata.napdata)
        print (_colors.BOLD + _colors.OKBLUE + 
               "\n================================================="
               "==============="
               + _colors.ENDC)
    else: # Invalid program name
        # User requested default pproj parameter values
        print (_colors.BOLD + _colors.FAIL + 
               "\nSuiteOPT error: Invalid program name provided to print_parm."
               + _colors.ENDC)
        print (_colors.BOLD + _colors.WARNING + 
               "    Functions for which parameters can be printed: pasa, "
               "pproj, cg, and napheap."
               + _colors.ENDC)
    # -------- Free memory used by pasadata -------- #
    cSuiteOPT.pasa_terminate(&pasadata)
    # ---------------- Exit program ---------------- #
    return

# --------------------------- solve function definition ----------------------- #
# Calls the pasa function to solve the given problem using either pasa, pproj,
# napheap, or cg
def solve(problem):
    """
    Method for solving problems using SuiteOPT.

    'solve' is designed to solve polyhedral constrained optimization 
    problems of the form
   
    (1)        min                  f(x)  
           subject to    bl <= Ax <= bu, lo <= x <= hi

    where f(x) is assumed to be a continuously differentiable function
    and evaluation routines for f(x) and the gradient of f(x) are 
    required. For Quadratic Programs, where f(x) = 0.5 * x^T H x + c^T x, 
    the user can provide an evaluation routine to compute the product of 
    the Hessian with a vector, H x,  and the linear term in the cost 
    function. For Linear Programs, where f(x) = c^T x, the user need only 
    provide the linear term, c, to evaluate f(x).

    Additionally, 'solve' has routines designed to solve problems of
    the form (1) with specific objective or constraint sets
   
    Polyhedral projection problems 

    (1a)       min            0.5 * || x - y ||^2  
           subject to    bl <= Ax <= bu, lo <= x <= hi

    Napsack problems
    
    (1b)       min             0.5 * x^T D x + c^T x  
           subject to    bl <= a^T x <= bu, lo <= x <= hi

    where D is a diagonal matrix with diagonal entries d.

    Unconstrained problems

    (1c)   min f(x)  subject to  x in R^n

    Parameters
    ----------
    problem : An instance of the class SuiteOPT.problem

    Returns
    -------
    xnew : Numpy array containing solution to problem

    See Also 
    --------
    problem : Class within SuiteOPT module for containing SuiteOPT 
              problem data
    """

    # -------- Initialize variables -------- #
    # Struct to provide to pasa function
    cdef cSuiteOPT.PASAdata * pasadata
    # Memview of data to provide to pasadata struct
    cdef double[::1] x_arr_memview
    cdef double[::1] lambda_arr_memview
    cdef double[::1] lo_arr_memview
    cdef double[::1] hi_arr_memview
    cdef double[::1] bl_arr_memview
    cdef double[::1] bu_arr_memview
    cdef double[::1] Ax_arr_memview
    cdef int[::1]    Ap_arr_memview
    cdef int[::1]    Ai_arr_memview
    cdef double[::1] y_arr_memview 
    cdef double[::1] c_arr_memview 
    cdef double[::1] a_arr_memview 
    cdef double[::1] d_arr_memview 
    cdef np.ndarray[np.int32_t, ndim=1] ncol_np = np.empty(1, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] nrow_np = np.empty(1, dtype=np.int32)
    # Global variables for user defined function
    global py_objective
    global py_gradient
    global py_objgrad
    global py_hprod
    global py_cg_hprod
    # Additional variables
    cdef int num_error
    cdef double t0, wall_time

    # --------------------------- PASA running message ----------------------- #
    print (_colors.OKBLUE + _colors.BOLD +
           " =============================================================\n"
           "             SuiteOPT: Optimization Software Suite\n"
           " ============================================================="
           + _colors.ENDC)

    # ------------------------ Allocate memory for C pasa -------------------- #
    print (_colors.OKBLUE + _colors.BOLD +
           "\n  Allocating memory..."
           + _colors.ENDC, end='')
    sys.stdout.flush()

    # Allocate memory for pasadata and initialize parameters
    pasadata = cSuiteOPT.pasa_setup()

    # Finished allocating memory message
    print (_colors.OKBLUE + _colors.BOLD +
           ".................................: "
           + _colors.OKGREEN + "Done"
           + _colors.ENDC)

    # ------------------------- Import custom parameters --------------------- #
    print (_colors.OKBLUE + _colors.BOLD +
           "\n  Importing custom parameters..."
           + _colors.ENDC, end='')
    sys.stdout.flush()

    # This code copies parameters from user dict to pasadata
    # Create dictionaries of parameters structs for each solver
    cdef object pasa_parms = pasadata.Parms.pasa[0]
    cdef object pproj_parms = pasadata.Parms.pproj[0]
    cdef object napheap_parms = pasadata.Parms.napheap[0]
    cdef object cg_parms = pasadata.Parms.cg[0]
    # Check if user provided correct type for parameters
    if type(problem.parm) is dict: # Import custom parameter values
        if len(problem.parm) > 0:
            # Loop over all values in user's custom parameter dictionary
            for key in problem.parm:
                # Add custom parameter to solver dictionary
                pasa_parms[key] = problem.parm[key]
                pproj_parms[key] = problem.parm[key]
                napheap_parms[key] = problem.parm[key]
                cg_parms[key] = problem.parm[key]
            # Copy imported parameter values back into each solver's 
            # parameter struct
            pasadata.Parms.pasa[0] = pasa_parms
            pasadata.Parms.pproj[0] = pproj_parms
            pasadata.Parms.napheap[0] = napheap_parms
            pasadata.Parms.cg[0] = cg_parms
            # Finished importing parameters 
            print (_colors.OKBLUE + _colors.BOLD +
                   ".......................: "
                   + _colors.OKGREEN + "Done"
                   + _colors.ENDC)
        else: # user did not provide custom parm values
            # No parameters 
            print (_colors.OKBLUE + _colors.BOLD +
                   "..............: "
                   + _colors.OKGREEN + "None provided"
                   + _colors.ENDC)
    else: # user provided non-dict for custom parm values
        # Print warning message
        print (_colors.BOLD + _colors.WARNING + 
               "\n\n    SuiteOPT Warning: Custom parameter values (parm) "
               "provided\n    as incorrect type."
               + _colors.ENDC)
        print (_colors.WARNING + 
               "      User's custom parms provided as %s.\n" 
               "      Custom parms must be provided as <type 'dict'>."
               %str(type(problem.parm))
               + _colors.ENDC)
        print (_colors.WARNING + 
               "    Proceeding with SuiteOPT using default parameter values.\n"
               + _colors.ENDC)
        print (_colors.OKBLUE + _colors.BOLD +
               "                                .......................: "
               + _colors.WARNING + "Done"
               + _colors.ENDC)

    # ---------------------------- Check user inputs ------------------------- #
    print (_colors.OKBLUE + _colors.BOLD +
           "\n  Checking problem data..."
           + _colors.ENDC, end='')
    sys.stdout.flush()

    # Note: Passing dimension data to _check_data function requires referencing address
    # Set value of ncol_np and nrow_np
    if problem.ncol is not None:
        # Check if user value for ncol is positive
        if problem.ncol > 0:
            # Set ncol equal to user value
            ncol_np[0] = problem.ncol
        else:
            # Set ncol equal to zero
            ncol_np[0] = 0
    else:
        # Set ncol equal to zero
        ncol_np[0] = 0

    if problem.nrow is not None:
        # Check if user value for nrow is positive
        if problem.nrow > 0:
            # Set nrow equal to user value
            nrow_np[0] = problem.nrow
        else:
            # Set nrow equal to zero
            nrow_np[0] = 0
    else:
        # Set nrow equal to zero
        nrow_np[0] = 0

    # Initialize memview of input data
    x_arr_memview = problem.x
    lambda_arr_memview = problem.pylambda
    lo_arr_memview = problem.lo
    hi_arr_memview = problem.hi
    bl_arr_memview = problem.bl
    bu_arr_memview = problem.bu
    y_arr_memview  = problem.y
    c_arr_memview  = problem.c
    a_arr_memview  = problem.a
    d_arr_memview  = problem.d

    # Initialize counter for number of errors during _check_data calls
    num_error = 0 

    # Check numpy arrays of expected length ncol
    num_error += _check_data("x", x_arr_memview, "ncol", ncol_np)
    num_error += _check_data("lo", lo_arr_memview, "ncol", ncol_np)
    num_error += _check_data("hi", hi_arr_memview, "ncol", ncol_np)
    num_error += _check_data("y", y_arr_memview, "ncol", ncol_np)
    num_error += _check_data("c", c_arr_memview, "ncol", ncol_np)
    num_error += _check_data("a", a_arr_memview, "ncol", ncol_np)
    num_error += _check_data("d", d_arr_memview, "ncol", ncol_np)

    # Check numpy arrays of expected length nrow
    num_error += _check_data("pylambda", lambda_arr_memview, "nrow", nrow_np)
    num_error += _check_data("bl", bl_arr_memview, "nrow", nrow_np)
    num_error += _check_data("bu", bu_arr_memview, "nrow", nrow_np)

    # Check matrix A (for shape and type) TODO update for new memview stuff
    if _check_matrix(problem.A, nrow_np[0], ncol_np[0]):
        # User has provided matrix A resulting in error; terminate program
        return

    # Check function handles (for input args and callability)
    if (_check_func(problem.objective, "objective", 1) or
        _check_func(problem.gradient, "gradient", 1) or
        _check_func(problem.objgrad, "objgrad", 1) or
        _check_func(problem.cg_hprod, "cg_hprod", 1) or
        _check_func(problem.hprod, "hprod", 2)):
        # User has provided function resulting in error; terminate program
        return

    # -------- Import problem data to pasadata -------- #
    # Set nrow and ncol based on _check_data calls
    pasadata.nrow = nrow_np[0]
    pasadata.ncol = ncol_np[0]
    if problem.x is not None:
        pasadata.x  = &x_arr_memview[0]
    if problem.lo is not None:
        pasadata.lo = &lo_arr_memview[0]
    if problem.hi is not None:
        pasadata.hi = &hi_arr_memview[0]
    if problem.y is not None:
        pasadata.y = &y_arr_memview[0]
    if problem.c is not None:
        pasadata.c = &c_arr_memview[0]
    if problem.a is not None:
        pasadata.a = &a_arr_memview[0]
    if problem.d is not None:
        pasadata.d = &d_arr_memview[0]

    if nrow_np[0] > 0:
        # TODO: lambda has set use in python, need to figure out workaround
        #pasadata.lambda  = &lambda_arr_memview[0]
        if problem.bl is not None:
            pasadata.bl = &bl_arr_memview[0]
        if problem.bu is not None:
            pasadata.bu = &bu_arr_memview[0]
        # Check user's format of matrix data
        if type(problem.A) is csc_matrix: # Sparse matrix format
            # Extract pointers to data, column, and row pointers
            Ax_arr_memview = problem.A.data
            Ap_arr_memview = problem.A.indptr
            Ai_arr_memview = problem.A.indices
            # Set pasadata structure values for matrix A
            pasadata.Ax = &Ax_arr_memview[0]
            pasadata.Ap = &Ap_arr_memview[0]
            pasadata.Ai = &Ai_arr_memview[0]
        elif type(problem.A) is np.matrix: # Dense matrix format
            # Convert numpy matrix to sparse format
            A_csc = csc_matrix(problem.A, dtype=np.double)
            # Extract pointers to data, column, and row pointers
            Ax_arr_memview = A_csc.data
            Ap_arr_memview = A_csc.indptr
            Ai_arr_memview = A_csc.indices
            # Set pasadata structure values for matrix A
            pasadata.Ax = &Ax_arr_memview[0]
            pasadata.Ap = &Ap_arr_memview[0]
            pasadata.Ai = &Ai_arr_memview[0]

    # -------- Pass function information to pasadata -------- #
    if problem.objective is not None:
        # Set py_objective to user objective function
        py_objective = problem.objective
        # Pass objective function to pasadata struct
        pasadata.value = <c_objective> objective_callback
    if problem.gradient is not None:
        # Set py_gradient to user gradient function
        py_gradient = problem.gradient
        # Pass gradient function to pasadata struct
        pasadata.grad = <c_gradient> gradient_callback
    if problem.objgrad is not None:
        # Set py_objgrad to user objgrad function
        py_objgrad = problem.objgrad
        # Pass objgrad function to pasadata struct
        pasadata.valgrad = <c_objgrad> objgrad_callback
    if problem.hprod is not None:
        # Set py_hprod to user hprod function
        py_hprod = problem.hprod
        # Pass hprod function to pasadata struct
        pasadata.hprod = <c_hprod> hprod_callback
    if problem.cg_hprod is not None:
        # Set py_cg_hprod to user cg_hprod function
        py_cg_hprod = problem.cg_hprod
        # Pass cg_hprod function to pasadata struct
        pasadata.cg_hprod = <c_cg_hprod> cg_hprod_callback

    # Indicate data checks and imports have been finished (if no errors occured)
    if num_error > 0:
        # User received warning during data check but proceeding with PASA
        print (_colors.OKBLUE + _colors.BOLD +
               "\n                          .............................: "
               + _colors.WARNING + "Done"
               + _colors.ENDC)
    else:
        # User did not receive warning during data check
        print (_colors.OKBLUE + _colors.BOLD +
               ".............................: "
               + _colors.OKGREEN + "Done"
               + _colors.ENDC)

    # ------------------- External call to C pasa function ------------------- #
    print (_colors.OKBLUE + _colors.BOLD +
           "\n  Working on solving the problem...\n\n"
           + _colors.ENDC, end='')
    sys.stdout.flush()

    # Set print status to true to ensure user gets correct status
    pasadata.Parms.pasa.PrintStatus = True
    # Start timer for pasa
    t0 = time.time()
    # Call pasa
    cSuiteOPT.pasa(pasadata)
    # Stop timer for pasa
    wall_time = time.time() - t0
    # Add walltime to statistics
    problem.stats["walltime"] = wall_time
    # Get problem status
    status = pasadata.Stats.pasa.status
    # Print indication that problem has finished
    if status == 0:
        print (_colors.OKBLUE + _colors.BOLD +
               "                                   ....................: "
               + _colors.OKGREEN + "Done\n"
               + _colors.ENDC)
    else:
        print (_colors.OKBLUE + _colors.BOLD +
               "                                   ....................: "
               + _colors.WARNING + "Done\n"
               + _colors.ENDC)

    # ------------------ Copy results/statistics and wrapup ------------------ #
    # Copy problem solution
    newx = np.zeros(ncol_np[0], dtype=np.float64)
    for i in range(ncol_np[0]):
        newx[i] = pasadata.x[i]

    # Set problem type
    if pasadata.Stats.use_pasa == 1:
      problem.probtype = "Constrained"
    elif pasadata.Stats.use_cg == 1:
      problem.probtype = "Unconstrained"

    # Copy problem statistics
    # Create dictionaries of statistics structs for each solver
    cdef object pasa_stats = pasadata.Stats.pasa[0]
    #cdef object pproj_stats = pasadata.Stats.pproj[0]
    cdef object napheap_stats = pasadata.Stats.napheap[0]
    cdef object cg_stats = pasadata.Stats.cg[0]
    # Add each statistics struct to problem.stats dictionary
    # PASA statistics
    if pasadata.Stats.use_pasa == True:
        # Add statistics to problem.stats
        problem.stats["pasa"] = pasa_stats
    else:
        # No statistics generated for pasa
        problem.stats["pasa"] = "No statistics generated"
    # NAPHEAP statistics
    if pasadata.Stats.use_napheap == True:
        # Add statistics to problem.stats
        problem.stats["napheap"] = napheap_stats
    else:
        # No statistics generated for napheap 
        problem.stats["napheap"] = "No statistics generated"
    # CG statistics
    if pasadata.Stats.use_cg == True:
        # Add statistics to problem.stats
        problem.stats["cg"] = cg_stats
    else:
        # No statistics generated for cg
        problem.stats["cg"] = "No statistics generated"

    # --------------------- Free memory used by pasadata --------------------- #
    cSuiteOPT.pasa_terminate(&pasadata)

    # ----------------------------- Exit function ---------------------------- #
    return newx


# ----------------------------- Internal classes ----------------------------- #
# ----- _colors class ----- #
class _colors:
    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'
