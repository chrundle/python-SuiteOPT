# SuiteOPT
SuiteOPT is designed to solve optimization problems of the form

    min f(x)    subject to    bl <= Ax <= bu, lo <= x <= hi

where `x` is a vector of length `ncol`, `f(x)` is a smooth (i.e., continuously differentiable) real-valued function, `A` is a real-valued `nrow` by `ncol` matrix, `bl` and `bu` are real-valued vectors of length `nrow`, and `lo` and `hi` are real-valued vectors of length `ncol`. In problems where `f(x)` is quadratic or linear, the user can provide special evaluation routines and/or vectors for optimized performance. Additionally, SuiteOPT can solve unconstrained or partially constrained problems by omitting any of the constraints when specifying the problem data.

## Installation
SuiteOPT for Python has been tested with Python versions 2.7, 3.5, 3.6, and 3.7 on Ubuntu. SuiteOPT for python can be installed using the package manager [**pip**](https://pip.pypa.io/en/stable/) which will build all of the required dynamic libraries for SuiteOPT with Python. As part of the installation routine, SuiteOPT for Python will build the required dynamic libraries for the C/C++ packages [**SuiteSparse**](https://github.com/jluttine/suitesparse) and [**SuiteOPT**](http://users.clas.ufl.edu/hager/papers/Software/). The dynamic libraries for the C/C++ packages SuiteSparse and SuiteOPT require BLAS and LAPACK. The default installation is configured to use the [**OpenBLAS**](https://www.openblas.net/) library for BLAS and LAPACK if installed on the user's system, however, users can provide desired libraries for BLAS and LAPACK using custom installation options.

#### 1. Installing Dependencies
Several dependencies are required before calling the `pip install` command to install SuiteOPT for Python. These dependencies are outlined below.
##### 1.1 Python Dependencies: Setuptools, Numpy, and Scipy
Inside the setup.py file used to install SuiteOPT, the Python packages [**setuptools**](https://pypi.org/project/setuptools/) and [**numpy**](https://numpy.org/) are imported to complete the installation and the package and [**scipy**](https://www.scipy.org/) is required for SuiteOPT to run. As a result, users missing one or more of these packages when installing SuiteOPT will encounter an error. These packages can be installed using [**pip**](https://pip.pypa.io/en/stable/) by entering any of the following commands into the terminal to install the packages missing on the user's system:

```sh
$ pip install setuptools
$ pip install numpy
$ pip install scipy
```

Note: Numpy version >= 1.16.x should be installed as errors have been observed with previous versions. For users with Python 2.7 we have found that Scipy version = 1.2.2 works successfully while users with Python version 3.5 or greater can use Scipy 1.3.1 as it is compatible with SuiteOPT. Other versions may work as well but no errors with SuiteOPT involving Numpy or Scipy have occurred when using these versions.

##### 1.2 SuiteOPT and SuiteSparse Dependencies: BLAS, LAPACK, and CMake
Users with BLAS, LAPACK, and CMake already available on their system may proceed to step 2.

Users without BLAS and LAPACK on should choose desired BLAS and LAPACK libraries suitable for their system. SuiteOPT is configured to use [**OpenBLAS**](https://www.openblas.net/) by default but can be customized at installation to use other BLAS and LAPACK configurations. SuiteOPT has been successfully tested with [**OpenBLAS**](https://www.openblas.net/) and [**Intel's MKL BLAS and LAPACK**](https://software.intel.com/en-us/mkl/choose-download). Note: Users should ensure that all dependencies required by their version of BLAS and LAPACK are installed prior to proceeding to the following step.

Users without [**CMake**](https://cmake.org/) should install it on their system before proceeding to the following step. CMake is required to build the dynamic library for [**METIS**](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) which is used by SuiteSparse.

#### 2. Installing SuiteOPT for Python
If [**OpenBLAS**](https://www.openblas.net/) is installed on your system and the directory containing `libopenblas` and `liblapack` is located in a directory found during compilation by default, then follow the instructions for **Default Installation**. Otherwise, follow the instructions for **Custom Installation**.

- **Default Installation**
    Type the following command into the terminal window to install SuiteOPT for Python with default BLAS and LAPACK configurations:

    ```sh
    $ pip install SuiteOPT --user
    ```

    Note: If any errors involving BLAS and LAPACK are raised during the installation routine then it is likely that the installer was unable to identify information necessary for using the user's BLAS and LAPACK libraries when building the dynamic libraries for the C/C++ packages SuiteSparse and SuiteOPT. Users experiencing this problem should follow the instructions in **Custom Installation**.

- **Custom Installation**
    Users wishing to use custom BLAS and LAPACK configurations for SuiteOPT can pass them to pip during the installation procedure. There are two ways to provide custom BLAS and LAPACK configurations both of which are outlined below.

    **Option 1: Pass File Containing BLAS, LAPACK, and LDLIBS**.
    With this option, users can create a file with variables `BLAS`, `LAPACK`, and `LDLIBS` that specify their BLAS and LAPACK libraries and the location on these libraries on their system, LDLIBS. The file should have one line for each variable. An example configuration for Intel's MKL BLAS and LAPACK is provided below.

    **Example: Custom Configuration File Using Intel's MKL BLAS and LAPACK.**
    ```mk
    # Contents of user created file 'mklconfig.txt'
    BLAS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
    LAPACK = -lmkl_lapack95_lp64
    LDLIBS = -L/opt/intel/mkl/lib/intel64 -L/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64
    ```

  Alternatively, the user can omit the `LDLIBS` variable if either the location of the BLAS and LAPACK libraries is known by the compiler or the full path to each library is included. As an example, a configuration for OpenBLAS using the alternative convention is provided below.

    **Example: Custom Configuration File Using OpenBLAS.**
    ```mk
    # Contents of user created file 'openblasconfig.txt'
    BLAS = /full/path/to/libopenblas.so -lgfortran -lpthread
    LAPACK = /full/path/to/liblapack.a
    ```

  Once the user's configuration for BLAS and LAPACK is saved to a file of their choice, say `myconfig.txt`, the full path to the file can be provided as an `--install-option` during the pip installation routine for SuiteOPT by entering the following command into the terminal:

    ```sh
    $ pip install SuiteOPT --user --install-option="--blfile=/full/path/to/myconfig.txt"
    ```

  Note: In particular, the `--install-option` is used to set the flag `--blfile` equal to the full path to the user's BLAS and LAPACK configuration file.

    **Option 2: Pass BLAS, LAPACK, and LDLIBS as Individual `--install-option` Flags**.
    With this option, the user can specify their BLAS and LAPACK libraries and the location(s) of these libraries directly as part of the pip install command by using `--install-option` flags to set the `--blas`, `--lapack`, and `--ldlibs` flags equal to the user's BLAS libraries, LAPACK libraries, and the absolute path location(s) of their BLAS and LAPACK libraries, respectively. Examples using various BLAS and LAPACK configurations are provided below.

    **Example: Custom Configuration Using Intel's MKL BLAS and LAPACK.**
    ```sh
    $ pip install SuiteOPT --user \
    --install-option="--blas=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm" \
    --install-option="--lapack=-lmkl_lapack95_lp64" \
    --install-option="--ldlibs=-L/opt/intel/mkl/lib/intel64 -L/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64"
    ```
    **Example: Custom Configuration Using OpenBLAS.**
    ```sh
    $ pip install SuiteOPT --user \
    --install-option="--blas=-lopenblas -lgfortran -lpthread" \
    --install-option="--lapack=-llapack" \
    --install-option="--ldlibs=-L/full/path/to/openblas/directory"
    ```

#### Notes on Troubleshooting Issues with BLAS and LAPACK
Steps 0 -- 2 should be sufficient for installing SuiteOPT, however if errors are occurring at runtime involving BLAS and/or LAPACK please refer to this section to check for potential solutions to your problems. During runtime, SuiteOPT for Python will require the user's BLAS and LAPACK libraries and the C/C++ dynamic libraries for SuiteSparse and SuiteOPT compiled during step 2 of the installation. The installation routine in the `setup.py` file is designed to use the provided BLAS and LAPACK information to prevent runtime errors involving BLAS and LAPACK from occurring. If any of these libraries cannot be located runtime errors may occur when using SuiteOPT for Python. This step outlines how to set certain environment variables so that the necessary libraries can be found at runtime.

For users encountering runtime errors involving the libraries `libpasa`, `libpproj`, `libcg_descent`, `libnapheap`, `libamd`, `libcamd`, `libcolamd`, `libccolamd`, `libcholmod`, or `libmetis`, see 1. For users encountering runtime errors involving BLAS or LAPACK libraries see 2.

1. **Add SuiteSparse and SuiteOPT dynamic library directories to `$LD_LIBRARY_PATH` environment variable.**
    Directories containing the dynamic libraries built by SuiteOPT and SuiteSparse need to be found in the user's `$LD_LIBRARY_PATH`. During installation, the required dynamic libraries for SuiteOPT and SuiteSparse are saved in directories called `SuiteOPT-libraries` and `SuiteSparse-libraries`, respectively, and should be located where python saves installed package information. The full path to these directories should be added to the `$LD_LIBRARY_PATH` environment variable. Examples are provided below for various shell configurations.

    **Example: Bourne/bash shell**
    ```bash
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/full/path/to/SuiteOPT-libraries"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/full/path/to/SuiteSparse-libraries"
    ```

    **Example: TCSH/CSH shell**
    ```sh
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/full/path/to/SuiteOPT-libraries
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/full/path/to/SuiteSparse-libraries
    ```

2. **Add BLAS and LAPACK dynamic library directories to `$LD_LIBRARY_PATH` environment variable.**
    Directories containing the user's BLAS and LAPACK dynamic libraries also need to be found in the user's `$LD_LIBRARY_PATH`. Below we provide information and examples for two popular versions of BLAS and LAPACK that should translate to other BLAS and LAPACK configurations.

    **Users with Openblas.** When using OpenBLAS, users should be able to bypass this step as it is common for the directory containing `libopenblas` and `liblapack` to be contained in a subdirectory of `/usr/lib` (on Linux) and, hence, located in a directory that is checked by default for dynamic libraries. However, if an error is encountered indicating that OpenBLAS cannot be found, then add the absolute path to the directory containing `libopenblas` and `liblapack` to the `$LD_LIBRARY_PATH` environment variable. Examples are provided below for various shell and BLAS/LAPACK configurations.

    **Example: Bourne/bash shell using Openblas**
    ```bash
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/full/path/to/openblas/dir"
    ```

    **Example: TCSH/CSH shell using Openblas**
    ```bash
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/full/path/to/openblas/dir
    ```

  **Users with Intel's MKL BLAS and LAPACK.** It is not common for the directory containing [**Intel's MKL BLAS and LAPACK**](https://software.intel.com/en-us/mkl/choose-download) to be installed in a directory that is checked by default for dynamic libraries (at least on Linux operating systems). Hence, the absolute path to the directory containing Intel's MKL BLAS and LAPACK dynamic libraries should be added to the `$LD_LIBRARY_PATH` environment variable. Additionally, it is a common requirement to link to the library `libiomp5` included with the Intel MKL package when using Intel's MKL BLAS and LAPACK and, thus, the user should also add the full path to the directory containing the `libiomp5` dynamic library to the `$LD_LIBRARY_PATH` environment variable.

  Additionally, users should add MKL dynamic libraries to the `$LD_PRELOAD` environment variable. When the MKL dynamic libraries are not added to the `$LD_PRELOAD` environment variable it has been observed that either (a) runtime errors or (b) poor performance of SuiteOPT for Python occurs. Some information outlining the prescribed solution included here was found in [**Intel's Forum on MKL**](https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/748309). Examples are provided below for various shell configurations but note that the required dynamic libraries may vary from system to system.

  **Example: Bourne/bash shell using Intel's MKL BLAS and LAPACK**
  ```bash
  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64"
  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64"
  export LD_PRELOAD="${LD_PRELOAD}:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_def.so:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_avx2.so:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_intel_lp64.so:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_intel_thread.so:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_core.so:/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64/libiomp5.so"
  ```

  **Example: TCSH/CSH shell using Intel's MKL BLAS and LAPACK**
  ```bash
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64
  setenv LD_PRELOAD ${LD_PRELOAD}:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_def.so:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_avx2.so:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_intel_lp64.so:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_intel_thread.so:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_core.so:/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64/libiomp5.so
  ```

## Usage
We now outline how to import the SuiteOPT module, set up a problem to solve with SuiteOPT, how to solve a problem with SuiteOPT, how to view problem statistics with SuiteOPT, and conclude with references to several demo files illustrating how to solve various problem configurations using SuiteOPT. Throughout the usage section, we provide sample blocks of code illustrating how to work with each component of SuiteOPT.
#### Importing SuiteOPT module
The SuiteOPT module can be imported within python by using the following command in python:
```python
import SuiteOPT
```
#### `SuiteOPT.problem` class
Problem data is expected to be provided as an instance of the class `SuiteOPT.problem`. The attributes of the `SuiteOPT.problem` class are listed below. Note that each attribute of the `SuiteOPT.problem` class that is an array must be provided as a [**numpy**](https://numpy.org) array.
- `ncol`      : dimension of initial guess x (also, number of columns in constraint matrix A)
- `nrow`      : number of rows in matrix A
- `x`         : Initial guess for primal variable (numpy array of length ncol)
- `lo`        : lower bounds for x (numpy array of length ncol)
- `hi`        : upper bounds for x (numpy array of length ncol)
- `A`         : constraint matrix (numpy matrix or scipy.sparse.csc_matrix of dimension nrow by ncol)
- `bl`        : lower bounds for Ax (numpy array of length nrow)
- `bu`        : upper bounds for Ax (numpy array of length nrow)
- `pylambda`  : multiplier for Ax (dual variable; numpy array of length nrow)
- `y`         : Project y onto polyhedron (if solving polyhedral projection problem; numpy array of length ncol)
- `a`         : Napsack constraint vector bl <= a'x <= bu (if solving napsack problem; numpy array of length ncol)
- `d`         : Diagonal of hessian (if solving napsack problem; numpy array of length nrow)
- `c`         : Linear term in objective function (if using Quadratic or Linear mode)
- `objective` : Objective function (Python function)
- `gradient`  : Gradient of the objective function (Python function)
- `objgrad`   : Objective function and its gradient (Python function)
- `hprod`     : Product of objective hessian times vector (if using quadratic mode; Python function)
- `cg_hprod`  : Product of objective hessian times vector (unconstrained problem in quadratic mode; Python function)
- `parm`      : Custom `SuiteOPT.solve()` parameter values (Python dict; keys = parm names, values = parm values)
- `stat`      : Problem statistics after solving problem (Python dict; keys = stat names, values = stat values)

##### Initializing `SuiteOPT.problem` class
To initialize an instance of the class containing the problem data, for example, call
```python
problem = SuiteOPT.problem(ncol, nrow)
```
If the problem has no matrix `A` then the value for `nrow` can be omitted and it will be initialized to zero. Alternatively, the values for `ncol` and `nrow` can be omitted upon initialization and they will be set to zero.

##### Setting problem data
**Initial Guess.** If the user wishes to provide an initial guess for the primal (decision) variable then the attribute `x` should be updated after initializing an instance of the `SuiteOPT.problem` class. If the user wishes to provide an initial guess for the dual variable (KKT multiplier) then the attribute `pylambda` should be updated after initializing an instance of the `SuiteOPT.problem` class.

```python
# Suppose that x and pylambda have been set to desired numpy arrays
problem.x  = x   # Initial guess for solution
problem.pylambda  = pylambda   # Initial guess for dual variable
```

**Bound Constraints.** If the user wishes to enforce the constraint `lo <= x` then `lo` is required. If the user wishes to enforce the constraint `x <= hi` then `hi` is required. If `lo = hi` for the user's problem, then `lo` and `hi` must both be provided with `lo = hi`. If `lo` is not provided then it will be assumed that there is no lower bound on `x`. Similarly, if `hi` is not provided then it will be assumed that there is no upper bound on `x`.

```python
# Suppose that lo and hi have been set to desired numpy arrays
problem.lo = lo  # Lower bound on primal variable,   lo <= x
problem.hi = hi  # Upper bound on primal variable,   x <= hi
```

**Polyhedral Constraints.** If the user wishes to enforce the constraint `bl <= Ax <= bu` then `A` and at least one of the vectors `bl` and `bu` are required. If `bl` is not provided then it will be assumed that there is no lower bound on `Ax`. Similarly, if `bu` is not provided then it will be assumed that there is no upper bound on `Ax`. If `bl = bu` in your problem, both `bl` and `bu` must be provided as the same vector. If the user does not wish to enforce this constraint then `A`, `bl`, and `bu` should not be updated after initializing an instance of the `SuiteOPT.problem` class.

```python
# Suppose that bl and bl have been set to desired numpy arrays and A set to desired numpy matrix or csc_matrix
problem.A  = A   # Linear constraint matrix A,       bl <= Ax <= bu
problem.bl = bl  # Lower bound on linear constraint, bl <= Ax
problem.bu = bu  # Upper bound on linear constraint, Ax <= bu
```

**Special Attributes.** The attributes `y`, `a`, `d`, and `c` should only be provided for certain problem formulations. In particular, special implementations are used by SuiteOPT for solving polyhedral projection, separable convex quadratic knapsack, and unconstrained problems. More information on these problem formulations and how to solve them efficiently using SuiteOPT can be found in the **Special Modes and Problem Formulations in SuiteOPT** subsection.

**Objective Function.** As a general rule for determining which function attributes in `SuiteOPT.problem` are necessary for solving a problem, ensure that the attributes provided in `SuiteOPT.problem` are sufficient for computing **(a) the objective function** and **(b) the gradient of the objective function**. For straightforward use note that it is always sufficient to provide the attributes `SuiteOPT.problem.objective` and `SuiteOPT.problem.gradient`. However, alternative choices of the function attributes can be used depending on the type of problem in order to obtain improved performance. More detailed information on which attributes of `SuiteOPT.problem` are sufficient and optimized for solving problems of various types (nonlinear objective, quadratic objective, ...) can be found in the **Special Modes and Problem Formulations in SuiteOPT** subsection.

**Customizing `SuiteOPT.solve()` Parameters.** If the user wishes to provide custom parameter values for the solver `SuiteOPT.solve()`, the user can modify the values of each parameter by name by inserting it into the `SuiteOPT.problem.parm` dictionary. To see a list of all parameters and their default values, the user can call the method `SuiteOPT.print_parm()`. The following calls to this function can be used to view full or partial lists of customizable parameters within SuiteOPT:

    SuiteOPT.print_parm("all")     : List of all default parameters values
    SuiteOPT.print_parm("parm")    : List of default pasa parameter values
    SuiteOPT.print_parm("pproj")   : List of default pproj parameter values
    SuiteOPT.print_parm("cg")      : List of default cg parameter values
    SuiteOPT.print_parm("napheap") : List of default napheap parameter values

Once the desired parameters have been identified, they can be modified by inserting them into the `SuiteOPT.problem.parm` dictionary attribute using the parameter name as the key and the desired parameter value as the value.
```python
# Suppose we want to set the parameter 'PrintLevel' to be 3 (maximum printing during SuiteOPT)
# and the stopping tolerance, 'grad_tol' to be 1e-8
parm = dict(PrintLevel=3, grad_tol=1e-8)
# Provide custom parameters to instance of SuiteOPT.problem class
problem.parm = parm  # custom parameter values (optional)
```

#### Special Modes and Problem Formulations in SuiteOPT

##### Nonlinear Mode
Here we outline which functions are required and optional when `f(x)` is a nonlinear cost function.

- Required Functions : `SuiteOPT.problem.objective` and `SuiteOPT.problem.gradient`
- Optional Function  : `SuiteOPT.problem.objgrad`

Note that `objective` and `gradient` must be functions with one input (vector to evaluate function at) and one output (function/gradient value, respectively). `objgrad` must be a function with one input (vector to evaluate function/gradient at) and two outputs (function value, gradient value). `objgrad` can be provided in instances where values computed during the computation of the objective function can be reused during the computation of the gradient resulting in reduced computation cost by computing the objective and gradient values simultaneously. If the user provides `objgrad` then when SuiteOPT requires updates to both the objective and gradient the `objgrad` function will be used. Otherwise, SuiteOPT will make a call to the `objective` and `gradient` functions.

##### Quadratic Mode
If the objective `f(x)` is quadratic, say of the form `f(x) = 0.5 x^T H x + c^T x`, the user should provide a routine, `hprod`, to evaluate the product between the objective Hessian `H` and the nonzero components of a vector `x`, so as to compute `Hx`, and the linear term in the objective, `c`.

- Required Attributes : `SuiteOPT.problem.hprod` and `SuiteOPT.problem.c`

Note that `hprod` must be a function with two inputs and one output. Using the above form for `f(x)`, `hprod` should compute the vector `Hx` given a vector `x` and an array containing the nonzero indices of `x`. In particular, the first input of the function `SuiteOPT.problem.hprod` should be the vector to take the product with, the second input should be a vector containing the indices for which the input vector is nonzero, and the output should be the matrix-vector product.

For unconstrained problems with quadratic cost functions the user should instead provide `cg_hprod`. `cg_hprod` must be a function with one input, the vector `x` to compute the matrix-vector product with, and one output, the matrix-vector product `Hx`.

##### Linear Mode
If the objective `f(x)` is linear, say of the form `f(x) = c^T x`, the user should provide the linear term in the objective, `c`.

- Required Input : `SuiteOPT.problem.c`

##### Polyhedral Projection Problem
SuiteOPT implements a special routine, called [**PPROJ**](http://users.clas.ufl.edu/hager/papers/CG/pproj.pdf), for solving a polyhedral projection problem of the form

    min  0.5 * || x - y || ** 2    subject to    bl <= Ax <= bu, lo <= x <= hi

where `y` is the point to be projected onto the polyhedral set `{x : bl <= Ax <= bu, lo <= x <= hi}`. When solving a problem of this form, it is sufficient to provide the following information to compute the objective and gradient values:

- Required Input : `SuiteOPT.problem.y`

##### Separable Convex Quadratic Knapsack Problem
SuiteOPT implements a special routine, called [**NAPHEAP**](http://users.clas.ufl.edu/hager/papers/CG/knapsack.pdf), for solving a separable convex quadratic knapsack problem of the form

    min  0.5 * x^T D x - c^T x    subject to    bl <= a^T x <= bu, lo <= x <= hi

where `D` is a diagonal matrix with nonnegative diagonal equal to the vector `SuiteOPT.problem.d`, `a` is a vector of length `ncol`, and `bl` and `bu` are scalars. When solving a problem of this form, it is sufficient to provide the following information to compute the objective and gradient values:

- Required Input : `SuiteOPT.problem.d`, `SuiteOPT.problem.c`

Additionally, instead of providing a constraint matrix `A`, the user should provide a constraint vector as a numpy array in the attribute `SuiteOPT.problem.a`.

##### Unconstrained Problem
SuiteOPT implements a special routine, called [**CG_DESCENT**](http://users.clas.ufl.edu/hager/papers/CG/cg_code.pdf), for solving unconstrained optimization problems of the form

    min  f(x)    subject to    x in R^n

In addition to nonlinear mode, quadratic mode can be used to evaluate the objective and gradient for unconstrained problems. Be sure to follow the details in the subsection **Special Modes and Problem Formulations in SuiteOPT: Quadratic Mode** to use quadratic mode with unconstrained problems.

#### `SuiteOPT.solve()` method
`SuiteOPT.solve()` requires a single input which is an instance of the `SuiteOPT.problem` class and returns a single output argument which is a numpy array containing the problem solution.
```python
# Solve problem using SuiteOPT and store solution in xnew
xnew = SuiteOPT.solve(problem)
```

#### Problem Statistics
Problem statistics generated during a run of `SuiteOPT.solve()` are stored in the `SuiteOPT.problem.stats` attribute as a dictionary and can be printed using the method `SuiteOPT.problem.print_stats()`.

### SuiteOPT Examples
Example files illustrating how to setup an instance of `SuiteOPT.problem`, call `SuiteOPT.solve()`, and print `SuiteOPT.problem.stats` are available in the `Demo` directory at https://github.com/chrundle/python-SuiteOPT. The type of problem and corresponding demo file are listed below:

     Nonlinear Objective : demo.py
     Quadratic Objective : demoQP.py
     Linear Objective    : demoLP.py

Users with the [**pycutest**](https://pypi.org/project/pycutest/) package installed may also run the `demo_cutest.py` file to solve any problem from the [**CUTEst**](https://github.com/ralna/CUTEst/wiki) optimization test suite. When running `demo_cutest.py`, the user will be prompted to enter the name of the problem from the CUTEst test set that they wish to solve using SuiteOPT.

## License
[**GPLv2**](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
