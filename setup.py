import os
import sys
import site 
import numpy 
import subprocess 
from setuptools import find_packages
from setuptools import setup, Extension 
from setuptools.command.install import install as _install 

################################################################################

# Initialize package details that do not change
PACKAGE_NAME = "SuiteOPT"
VERSION = "1.0.3"
DESCRIPTION = "Python Wrapper for Optimization Software Suite SuiteOPT"
AUTHOR = "James Diffenderfer"
AUTHOR_EMAIL = "jdiffen1@ufl.edu"
URL = "https://github.com/chrundle"
KEYWORDS = "optimization"
LICENSE = "GPLv2"
CLASSIFIERS = [ 
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    "Operating System :: POSIX :: Linux",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
]
# Not using install requires at the moment; causes error on some 
# systems when passing custom user BLAS and LAPACK information
INSTALL_REQUIRES = ["numpy>=1.16.x", "scipy<=1.2.9"]
# read the contents of your README file
with open("README.md", "r") as fh:
    long_description = fh.read()

# Custom class for passing user BLAS and LAPACK information during installation
class SuiteOPTInstall(_install):
  user_options = _install.user_options + [
    ('blas=', None, None), # option that takes value
    ('lapack=', None, None), # option that takes value
    ('ldlibs=', None, None), # option that takes value
    ('blfile=', None, None), # option that takes value
  ]

  def initialize_options(self):
    _install.initialize_options(self)
    self.blas = None
    self.lapack = None
    self.ldlibs = None
    self.blfile = None

  def finalize_options(self):
    # Indicate BLFILE being used
    if self.blfile != None and os.path.exists(self.blfile):
      print("User BLFILE: %s" %(self.blfile))
      # Set blas, lapack, and ldlibs from blfile
      # Copy contents of Userconfig.mk file
      with open(self.blfile, "r") as myfile:
        for line in myfile:
          # Find line with BLAS and change to user BLAS if provided
          if line.startswith("BLAS"):
            self.blas = line[4:].replace("\n","")
            # Strip leading whitespace and = signs
            self.blas = self.blas.lstrip()
            self.blas = self.blas.lstrip("=")
            self.blas = self.blas.lstrip()
          # Find line with LAPACK and change to user LAPACK if provided
          if line.startswith("LAPACK"):
            self.lapack = line[6:].replace("\n","")
            # Strip leading whitespace and = signs
            self.lapack = self.lapack.lstrip()
            self.lapack = self.lapack.lstrip("=")
            self.lapack = self.lapack.lstrip()
          # Find line with LDLIBS and change to user LDLIBS if provided
          if line.startswith("LDLIBS"):
            self.ldlibs = line[6:].replace("\n","")
            # Strip leading whitespace and = signs
            self.ldlibs = self.ldlibs.lstrip()
            self.ldlibs = self.ldlibs.lstrip("=")
            self.ldlibs = self.ldlibs.lstrip()
    else:
      print("BLFILE not provided.")

    # Indicate BLAS being used
    if self.blas != None:
      print("User BLAS: %s" %(self.blas))
    else:
      print("BLAS not provided. Using default: -lopenblas")
    # Indicate LAPACK being used
    if self.lapack != None:
      print("User LAPACK: %s" %(self.lapack))
    else:
      print("LAPACK not provided. Using default: -llapack")
    # Indicate LDLIBS being used
    if self.ldlibs != None:
      print("User LDLIBS: %s" %(self.ldlibs))
    else:
      print("LDLIBS not provided. Using default: ")

    # Get current working directory
    cwd = os.getcwd()
    # Check if SuiteOPT C pacakge is sub or parent directory of cwd or 
    # environment variable
    if os.path.isdir(cwd + "/SuiteOPT"): 
      # SuiteOPT C package is subdirectory of cwd
      SUITEOPTDIR = cwd + "/SuiteOPT"
    elif os.path.isdir(cwd + "/../../SuiteOPT"): 
      # SuiteOPT C package is parent directory of cwd
      SUITEOPTDIR = os.path.abspath(os.path.join(cwd,os.pardir))

    # Copy contents of Userconfig.mk file
    with open(SUITEOPTDIR + "/SuiteOPTconfig/Userconfig.mk", "r") as myfile:
      # Save contents of Userconfig.mk to tmpdata
      tmpdata = myfile.readlines()
    # Set ndata equal to length of tmpdata
    ndata = len(tmpdata)
    # Find line with BLAS and change to user BLAS if provided
    if self.blas != None:
      for i in range(ndata-1,0,-1):
        if tmpdata[i].startswith("BLAS"):
          tmpdata[i] = "BLAS = " + self.blas + "\n"
          break
      # If MKL blas not user need to add DNSUPER flag to OPTFLAGS
      for i in range(ndata-1,0,-1):
        if tmpdata[i].startswith("OPTFLAGS +="):
          if "mkl" in self.blas:
            tmpdata[i] = "OPTFLAGS += \n"
          else:
            tmpdata[i] = "OPTFLAGS += -DNSUPER\n"
          break
    # Find line with LAPACK and change to user LAPACK if provided
    if self.lapack != None:
      for i in range(ndata-1,0,-1):
        if tmpdata[i].startswith("LAPACK"):
          tmpdata[i] = "LAPACK = " + self.lapack + "\n"
          break
    # Find line with LDLIBS and change to user LDLIBS if provided
    if self.ldlibs != None:
      for i in range(ndata-1,0,-1):
        if tmpdata[i].startswith("LDLIBS"):
          tmpdata[i] = "LDLIBS = " + self.ldlibs + "\n"
          break
    # Print updated data back to file
    with open(SUITEOPTDIR + "/SuiteOPTconfig/Userconfig.mk", "w") as myfile:
      # Save modified contents of tmpdata back to Userconfig.mk
      myfile.writelines(tmpdata)

    # Call standard finalize_options for install
    _install.finalize_options(self)

  def run(self):
    global blas, lapack, ldlibs, blfile
    blas = self.blas
    lapack = self.lapack
    ldlibs = self.ldlibs
    blfile = self.blfile

    # This attempts to solve problem of libpasa.so and pibpproj.so being
    # unable to locate SuiteSparse dynamic libraries at runtime from 
    # within python
    # Set environment variable for linking SuiteSparse libraries at runtime
    # Note: This environment variable is appended to LDLIBS in compilation
    # of dynamic libraries for SuiteOPT and SuiteSparse
    # ---- Set directories to search for dynamic libraries at runtime ---- #
    # Initialize list to empty list
    runtime_library_dirs = []
    # Add SuiteOPT and SuiteSparse dynamic library directories to all 
    # directories from getusersitepackages to list (typically for site-packages)
    # as these directories need to be searched for by SuiteSparse and SuiteOPT
    # dynamic libraries at runtime
    usersitepackages = site.getusersitepackages()
    if type(usersitepackages) == str:
      runtime_library_dirs += [usersitepackages + "/SuiteOPT-libraries"]
      runtime_library_dirs += [usersitepackages + "/SuiteSparse-libraries"]
    elif type(usersitepackages) == list:
      runtime_library_dirs += [x + "/SuiteOPT-libraries" for x in usersitepackages]
      runtime_library_dirs += [x + "/SuiteSparse-libraries" for x in usersitepackages]
    # Add SuiteOPT and SuiteSparse dynamic library directories to all 
    # directories from getsitepackages to list (typically for dist-packages)
    # as these directories need to be searched for by SuiteSparse and SuiteOPT
    # dynamic libraries at runtime
    userdistpackages = site.getsitepackages()
    if type(userdistpackages) == str:
      runtime_library_dirs += [userdistpackages + "/SuiteOPT-libraries"]
      runtime_library_dirs += [userdistpackages + "/SuiteSparse-libraries"]
    elif type(userdistpackages) == list:
      runtime_library_dirs += [x + "/SuiteOPT-libraries" for x in userdistpackages]
      runtime_library_dirs += [x + "/SuiteSparse-libraries" for x in userdistpackages]
    # SuiteOPT will be looking for BLAS and LAPACK libraries at runtime
    if self.ldlibs != None:
      # Append list of ldlibs provided by user after stripping -L
      runtime_library_dirs += [x[2:] for x in self.ldlibs.split()]
    # For all directories in runtime_library_dirs, prepend "-Wl,-rpath," so
    # so that compiler knows to search these directories at runtime
    runtime_library_dirs = ["-Wl,-rpath," + x for x in runtime_library_dirs]
    # Join all flags and directories into single string
    runtime_library_dirs = " ".join(runtime_library_dirs)
    os.putenv("PYSUITEOPT_RPATH", runtime_library_dirs)
    #os.putenv("SUITEOPT_RPATH", "-Wl,-rpath," + site.getusersitepackages() + "/SuiteSparse-libraries")
    # Make SuiteSparse and SuiteOPT libraries
    suiteopt_command = ["make"]
    # Check if errors occurred during make command
    if subprocess.call(suiteopt_command) != 0:
      sys.exit(-1)
    # Update LDFLAGS environment variable to include user LDLIBS 
    # This is required since the extension will be looking for these
    # at time of compilation of SuiteOPT.so
    if self.ldlibs != None:
      ldflags = self.ldlibs + " " + self.blas + " " + self.lapack + " "
      try:
        os.environ["LDFLAGS"] += ldflags
      except:
        os.environ["LDFLAGS"] = ldflags
    # Execute regular install commands
    _install.run(self)

# ---- Update information required for setup on user's system ---- #
# Initialize package details that are later updated
setup_requires = ["numpy>=1.16.x", "setuptools", "scipy"]

# Get current working directory
cwd = os.getcwd()
# Check if SuiteOPT C pacakge is sub or parent directory of cwd or 
# environment variable
if os.path.isdir(cwd + "/SuiteOPT"): 
  # SuiteOPT C package is subdirectory of cwd
  SUITEOPTDIR = "SuiteOPT/"
elif os.path.isdir(cwd + "/../../SuiteOPT"): 
  # SuiteOPT C package is parent directory of cwd
  SUITEOPTDIR = "../"

# ---- Set directories to search for dynamic libraries at compilation ---- #
# Initialize library_dirs list to be empty
library_dirs = []
# Append all required library directories to library_dirs
library_dirs.append(SUITEOPTDIR + "lib")
library_dirs.append(SUITEOPTDIR + "SuiteSparseX/lib")
#library_dirs += BLASDIRS

# ---- Set directories to search for dynamic libraries at runtime ---- #
# Initialize list to empty list
runtime_library_dirs = []
# Add SuiteOPT and SuiteSparse dynamic library directories to all 
# directories from getusersitepackages to list (typically for site-packages)
usersitepackages = site.getusersitepackages()
if type(usersitepackages) == str:
  runtime_library_dirs += [usersitepackages + "/SuiteOPT-libraries"]
  runtime_library_dirs += [usersitepackages + "/SuiteSparse-libraries"]
elif type(usersitepackages) == list:
  runtime_library_dirs += [x + "/SuiteOPT-libraries" for x in usersitepackages]
  runtime_library_dirs += [x + "/SuiteSparse-libraries" for x in usersitepackages]
# Add SuiteOPT and SuiteSparse dynamic library directories to all 
# directories from getsitepackages to list (typically for dist-packages)
userdistpackages = site.getsitepackages()
if type(userdistpackages) == str:
  runtime_library_dirs += [userdistpackages + "/SuiteOPT-libraries"]
  runtime_library_dirs += [userdistpackages + "/SuiteSparse-libraries"]
elif type(userdistpackages) == list:
  runtime_library_dirs += [x + "/SuiteOPT-libraries" for x in userdistpackages]
  runtime_library_dirs += [x + "/SuiteSparse-libraries" for x in userdistpackages]

# ---- Set extra link arguments to add directories to search at runtime ---- #
# For some reason adding the runtime_library_dirs only works for SuiteOPT
# dynamic libraries and not SuiteSparse, trying workaround
# Initialize extra_link_args
extra_link_args = []
# Append all runtime_library_dirs with -Wl and -rpath flags
for d in runtime_library_dirs:
  extra_link_args.append("-Wl,-rpath," + d) 

# ---- Set directories to search for header files at compilation ---- #
# Initialize include_dirs list to contain numpy.get_include()
include_dirs = [numpy.get_include()]
# Append all required include directories to include_dirs
include_dirs.append(SUITEOPTDIR + "include")
include_dirs.append(SUITEOPTDIR + "SuiteSparseX/include")

# ---- Set list of libraries required at compilation ---- #
# Initialize libraries list to contain required SuiteOPT and SuiteSparse libs
libraries = ["pasa", "pproj", "cg_descent", "napheap", "cholmod",
             "ccolamd", "colamd", "camd", "amd", "metis"]
# Add all additional user libraries to libraries list
#libraries += USERLIBS

# ---- Set extra compiler arguments ---- #
extra_compile_args = ["-std=c99", "-lm", "-lpthread", "-fopenmp", "-Wno-cpp", "-O3"]

# ---- Set custom command class (custom install routine so user can provide BLAS/LAPACK) ---- #
# Set cmdclass arguments
cmdclass = {}
# Customize install cmdclass to allow user to pass custom BLAS and LAPACK
cmdclass['install'] = SuiteOPTInstall

# ---- Set packages to contain dynamic libraries from SuiteSparse and SuiteOPT ---- #
# Initialize pacakges to empty list
packages = []
# Initialize package_dir as empty dicionary
package_dir = {}
# Initialize package_data as empty dicionary
package_data = {}

# ---- Add SuiteOPT dynamic library data to packages ---- #
# Add SuiteOPT-libraries to list of packages
packages.append("SuiteOPT-libraries")
# Add SuiteOPT/lib to SuiteOPT-libraries package directories in dictionary
package_dir["SuiteOPT-libraries"] = SUITEOPTDIR + "lib"
# Add SuiteOPT dynamic libraries to SuiteOPT-libraries package data in dictionary
package_data["SuiteOPT-libraries"] = ["*.so*"]

# ---- Add SuiteSparse dynamic library data to packages ---- #
# Add SuiteSparse-libraries to list of packages
packages.append("SuiteSparse-libraries")
# Add SuiteOPT/SuiteSparseX/lib to SuiteSparse-libraries package directories in dictionary
package_dir["SuiteSparse-libraries"] = SUITEOPTDIR + "SuiteSparseX/lib"
# Add SuiteSparse dynamic libraries to SuiteSparse-libraries package data in dictionary
package_data["SuiteSparse-libraries"] = ["*.so*"]

# USE_CYTHON should be True to rebuild SuiteOPT.c file using cython
# USE_CYTHON should be False to build SuiteOPT.so library from prebuilt SuiteOPT.c file (for sdist)
USE_CYTHON = False

# Set extensions, setup_requires, and run setup
if USE_CYTHON:
  # Using cython
  from Cython.Build import cythonize 
  # Setup requires
  setup_requires += ["Cython>=0.29.x"]
  # Include path specifies where the SuiteOPT.pxd file is located and
  # is a cythonize option (not an Extension option)
  include_path = ["pySuiteOPT"]
  # Set extensions
  extensions = cythonize([Extension("SuiteOPT",
                                    ["pySuiteOPT/SuiteOPT.pyx"],
                                    include_dirs = include_dirs,
                                    libraries = libraries,
                                    library_dirs = library_dirs,
                                    extra_link_args = extra_link_args,
                                    extra_compile_args = extra_compile_args)],
                         include_path = include_path)
  # Setup
  setup (
    name = PACKAGE_NAME,
    version = VERSION,
    description = DESCRIPTION,
    author = AUTHOR,
    author_email = AUTHOR_EMAIL,
    url = URL,
    keywords = KEYWORDS,
    license = LICENSE,
    classifiers = CLASSIFIERS,
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    cmdclass = cmdclass,
    packages = packages,
    package_dir = package_dir,
    package_data = package_data,
    include_package_data = True,
    ext_modules = extensions
  )
else:
  # Not using cython
  extensions = [Extension("SuiteOPT", ["pySuiteOPT/SuiteOPT.c"], 
                          include_dirs = include_dirs,
                          libraries = libraries,
                          library_dirs = library_dirs,
                          extra_link_args = extra_link_args,
                          extra_compile_args = extra_compile_args)]
  # Setup
  setup (
    name = PACKAGE_NAME,
    version = VERSION,
    description = DESCRIPTION,
    author = AUTHOR,
    author_email = AUTHOR_EMAIL,
    url = URL,
    keywords = KEYWORDS,
    license = LICENSE,
    classifiers = CLASSIFIERS,
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    cmdclass = cmdclass,
    packages = packages,
    package_dir = package_dir,
    package_data = package_data,
    include_package_data = True,
    ext_modules = extensions
  )
