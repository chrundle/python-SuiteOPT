# ------------------ Define Type Macros ------------------ #
ctypedef int PASAINT
ctypedef double PASAFLOAT
ctypedef int PPINT
ctypedef double PPFLOAT
ctypedef int CGINT
ctypedef double CGFLOAT
ctypedef int NAPINT
ctypedef double NAPFLOAT

cdef extern from "pasa.h":
    # -------- pasa functions -------- #
    void pasa (PASAdata *pasadata) 
    PASAdata *pasa_setup ()
    void pasa_terminate (PASAdata **Data)
    void pasa_print_parm (PASAdata *Data)
    void pasa_print_status (int status, PASAstats *Stats)

    # -------- pasa structs -------- #
    ctypedef struct PASAdata:
        # -------- pasa input data -------- #
        PASAFLOAT         *x # size ncol, points to a guess for the solution
                             # when pasa is run; at completion of pasa, points
                             # to the computed solution 
        PASAFLOAT   *plambda # If not NULL and pasaparm->use_lambda = T,
                             #   then lambda points to an array of size nrow
                             #   (allocated by the user) that contains a
                             #   guess for the multiplier. If not NULL,
                             #   then at completion, lambda points to
                             #   the computed multiplier.
        PASAparms     *Parms # NULL => use default parameters, otherwise
                             # contains up to 4 non-NULL structures pasa, cg,
                             # pproj, and napheap with parameters
        PASAstats     *Stats # NULL => do not return statistics structure.
                             #   Otherwise, up to 4 non-NULL structures are
                             #   returned with statistics for the codes used.
        PASAINT         nrow # number of rows in A
        PASAINT         ncol # number of components in x and number of cols in
                             # A if it exists
        # If A exists, then its nonzero elements are stored in sparse matrix
        # format corresponding to two integer arrays and one real array. When
        # A does not exist, these arrays should be NULL 
        PASAINT          *Ap # size ncol+1, column pointers
        PASAINT          *Ai # size Ap [ncol], row indices for A, increasing
                             # order in each column
        PASAFLOAT        *Ax # size Ap [ncol], numerical entries of A
        PASAFLOAT        *bl # size nrow, lower bounds for inequalities
                             # NULL => use -infinity
        PASAFLOAT        *bu # size nrow, upper bounds for inequalities
                             # NULL => use +infinity
        PASAFLOAT        *lo # size ncol, lower bounds for x
                             # NULL => use -infinity
        PASAFLOAT        *hi # size ncol, upper bounds for x
                             # NULL => use +infinity
        PASAFLOAT         *y # size ncol, project y onto polyhedron in PPROJ
                             # NULL => not polyhedral projection problem
        PASAFLOAT         *c # size ncol, linear term if objective quadratic
                             # NULL => objective not quadratic
        PASAFLOAT         *a # size ncol, napsack constraint *bl <= a'x <= *bu
                             # NULL => not napsack problem
        PASAFLOAT         *d # size ncol, diagonal of Hessain in NAPHEAP
                             # NULL => no diagonal Hessian in NAPHEAP
        PASAFLOAT     *xWork # NULL => pasa should allocate real work space.
                             # Otherwise, see allocation section of pasa to
                             # determine memory requirements.
        PASAINT       *iWork # NULL  => pasa should allocate integer work space.
                             # Otherwise, see allocation section of pasa to
                             # determine memory requirements.
        CGdata       *cgdata # input data for CG_DESCENT
        PPdata       *ppdata # input data structure for PPROJ
        NAPdata     *napdata # input data for NAPHEAP
        # for a quadratic objective, hprod evaluates Hessian times vector,
        # NULL if objective not quadratic
        void   (*hprod) (double *, double *, int *, int, int)
        # for a quadratic objective in unconstrained problem, cg_hprod evaluates 
        # Hessian times vector, NULL if problem constrained
        void   (*cg_hprod) (double *, double *, int)
        # for a general nonlinear function, value (f, x, n) is function value,
        # NULL if objective is quadratic
        void   (*value) (double *, double *, int)
        # for a general nonlinear function, grad (g, x, n) is function gradient,
        # NULL if objective is quadratic
        void    (*grad) (double *, double *, int) # grad (g, x, n)
        # for a general nonlinear function, valgrad (f, g, x, n) gives both
        # function value f and gradient g at x (size n), argument can be NULL
        void (*valgrad) (double *, double *, double *, int)
        # Hessian evaluation routine, can be NULL, not currently used
        #void    (*hess) (PASAhess *)
        # Hessian evaluation routine for cg, can be NULL, not currently used
        #void    (*cghess) (CGhess *)
        # the following are used internally when the code allocates memory
        PASAFLOAT *x_created ;      # pointer to memory created for x
        PASAFLOAT *lambda_created ; # pointer to memory created for lambda

    ctypedef struct PASAparms:
        PASAparm   *pasa # parameters for PASA, NULL => use default
        CGparm       *cg # parameters for CG_DESCENT, NULL => use default
        PPparm    *pproj # parameters for PPROJ, NULL => use default
        NAPparm *napheap # parameters for NAPHEAP, NULL => use default

    ctypedef struct PASAstats:
        PASAstat   *pasa # statistics for PASA
        CGstat       *cg # statistics for CG_DESCENT
        PPstat    *pproj # statistics for PPROJ
        NAPstat *napheap # statistics for NAPHEAP
        int     use_pasa # T => pasa was used
        int       use_cg # T => cg was used
        int    use_pproj # T => pproj was used
        int  use_napheap # T => napheap was used

    ctypedef struct PASAparm:
        # PASA parameters; descriptions in PASA/Include/pasa.h
        int UNC ;     
        int BNC ;     
        int LP  ;     
        int QP  ;     
        int NL  ;     
        int NAPSACK ; 
        int PROJ ;    
        PASAFLOAT grad_tol ;
        int PrintStatus ;
        int PrintStat ;
        int PrintParm ;
        int PrintLevel ;
        int GradProjOnly ;
        int use_activeGP ;
        int use_napheap ;
        int use_hessian ;
        int loExists ;
        int hiExists ;
        int updateorder ; 
        PASAFLOAT epsilon ; 
        PASAFLOAT cerr_decay ;
        PASAFLOAT EpsilonGrow ;
        PASAFLOAT EpsilonDecay ;
        PASAFLOAT fadjust ;
        int use_lambda ;
        int pproj_start_guess ;
        int use_penalty ;
        PASAFLOAT penalty ;
        int          debug ;
        PASAFLOAT debugtol ;
        PASAFLOAT switchfactor; 
        PASAFLOAT switchdecay ; 
        int terminate_agp ;
        PASAINT testit ; 
        PASAFLOAT GPtol ; 
        PASAINT gpmaxit ; 
        PASAFLOAT restart_fac ;
        int            L ;
        int            M ; 
        int            P ;
        PASAFLOAT gamma1 ;
        PASAFLOAT gamma2 ;
        PASAFLOAT gamma3 ;
        PASAFLOAT        lambda0 ;
        PASAFLOAT  lambda0Factor ;
        PASAFLOAT            bbk ;
        PASAFLOAT       bbexpand ;
        PASAFLOAT bbSwitchFactor ;
        int MaximumCycle ;
        int NominalCycle ;
        int approxstep ;
        PASAFLOAT ArmijoSwitchFactor ;
        int  PertRule ;
        PASAFLOAT pert_eps ;
        int neps ;
        PASAFLOAT Armijo_delta ;
        PASAFLOAT Wolfe_delta ;
        PASAFLOAT Wolfe_sigma ;
        int maxsteps ;
        PASAFLOAT stepdecay ;
        PASAFLOAT safe0 ;
        PASAFLOAT safe1 ;
        PASAFLOAT infdecay ;
        PASAFLOAT infdecay_rate ;
        int     ninf_tries ;

    ctypedef struct PASAstat:
        int         status ; #  status at termination 
        int       maxsteps ; #  max number of attempts to satisfy the grad proj
                             #  line search condition 
        PASAFLOAT    lobad ; #  invalid lower bound for a variable 
        PASAFLOAT    hibad ; #  invalid upper bound for a variable 
        PASAINT       ibad ; #  index of the bad bound 
        PASAFLOAT grad_tol ; #  KKT error tolerance 
        int     ninf_tries ; #  number of tries to find finite objective value 
        PASAINT    gpmaxit ; #  max iteration in gradient projection algorithm 
        PASAFLOAT        f ; #  function value 
        PASAFLOAT      err ; #  || P (x - g) - x || 
        PASAINT       mcnf ; #  function evaluations in main code 
        PASAINT       mcng ; #  gradient evaluations in main code 
        PASAINT       gpit ; #  number of iterations in grad_proj 
        PASAINT       gpnf ; #  function evaluations in grad_proj 
        PASAINT       gpng ; #  gradient evaluations in grad_proj 
        PASAINT      agpit ; #  number of iterations in active set grad_proj 
        PASAINT      agpnf ; #  function evaluations in active set grad_proj 
        PASAINT      agpng ; #  gradient evaluations in active set grad_proj 
        PASAINT   nproject ; #  number of projections performed 

cdef extern from "pproj.h":
    # -------- pproj functions -------- #
    void pproj_print_parm (PPdata *Data)

    # -------- pproj structs -------- #
    ctypedef struct PPdata:
        pass

    ctypedef struct PPparm:
        # PPROJ parameters; descriptions in PPROJ/Include/pproj.h
        int    PrintStatus ; 
        int      PrintStat ; 
        int     PrintLevel ; 
        int      PrintParm ; 
        int use_prior_data ; 
        int    return_data ; 
        int       loExists ; 
        int       hiExists ; 
        int      getfactor ; 
        int          debug ; 
        PPFLOAT   checktol ; 
        int    start_guess ; 
        int        permute ; 
        PPFLOAT     phase1 ; 
        int        cholmod ; 
        int     multilevel ; 
        int stop_condition ; 
        PPFLOAT      sigma ; 
        PPFLOAT     Asigma ; 
        int     ScaleSigma ; 
        PPFLOAT sigma_decay ; 
        int          nprox ; 
        PPFLOAT armijo_grow; 
        int        narmijo ; 
        int            mem ; 
        int        nsparsa ; 
        PPFLOAT      gamma ; 
        PPFLOAT        tau ; 
        PPFLOAT       beta ; 
        PPFLOAT grad_decay ; 
        PPFLOAT gamma_decay; 
        int use_coor_ascent; 
        PPFLOAT   coorcost ; 
        int      use_ssor0 ; 
        int      use_ssor1 ; 
        int     use_sparsa ; 
        int    use_startup ; 
        PPFLOAT  ssordecay ; 
        PPFLOAT   ssorcost ; 
        int        ssormem ; 
        PPINT   ssormaxits ; 
        PPFLOAT  cutfactor ; 
        PPFLOAT    tolssor ; 
        PPFLOAT    tolprox ; 
        PPFLOAT tolrefactor; 
        int badFactorCutoff; 
        int             LP ; 
        PPFLOAT  LinFactor ; 
        PPFLOAT LinGrad_tol; 

    ctypedef struct PPstat:
        int             status ; #  status of the run 
        int         parm_nprox ; #  max number of proximal updates 
        int            cholmod ; #  TRUE => cholmod is employed 
        PPFLOAT          lobad ; #  invalid lower bound for a linear inequality 
        PPFLOAT          hibad ; #  invalid upper bound for a linear inequality 
        PPINT             ibad ; #  index of the bad bound 
        PPFLOAT       grad_tol ; #  relative tolerance for dual function gradient 
        PPINT       ssormaxits ; #  upper bound on number of ssor iterations 
        PPFLOAT        errdual ; #  actual error in the dual function gradient where
                                 #  the error measure is based on the stop condition
                                 #  0 => ||grad L (lambda)||_sup / absAx_sup
                                 #       where absAx_sup = max_i sum_j |a_{ij}x_j|
                                 #       x achieves min in dual function
                                 #  1 => ||grad L (lambda)||_sup
                                 #  2 => ||grad L (lambda)||_sup /(absAx_sup + ymax)
                                 #       where ymax is the sup-norm of the
                                 #       projection point y 
        PPINT             nrow ; #  number of rows in A 
        int           *updowns ; #  number of updates and downdates of each size 
        int       size_updowns ; #  dimension of updowns array 
        int            *solves ; #  size: maxdepth+1, # of solves by level  
        int           maxdepth ; #  number of levels in the partition tree 
        int               blks ; #  number of blocks in multilevel partition of A
                                 #  each separator is also counted as a block 
        int             nchols ; #  number of Cholesky factorizations 
        int              nprox ; #  number of proximal updates 
        PPINT       phase1_its ; #  number of iterations in phase 1 
        PPINT  coor_ascent_its ; #  number of coordinate ascent iterations 
        PPINT        ssor0_its ; #  number of ssor0 iterations 
        PPINT        ssor1_its ; #  number of ssor1 iterations 
        PPINT       sparsa_its ; #  number of SpaRSA iterations 
        PPINT            coldn ; #  number of rank 1 downdates to Cholesky factor 
        PPINT            colup ; #  number of rank 1 updates to Cholesky factor 
        PPINT            rowdn ; #  number of rows dropped from Cholesky factor 
        PPINT            rowup ; #  number of rows added to Cholesky factor 
        PPINT coor_ascent_free ; #  number of variables freed in coordinate ascent
        PPINT coor_ascent_drop ; #  number of rows dropped in coordinate ascent 
        PPINT       ssor0_free ; #  number of variables freed in ssor0 
        PPINT       ssor0_drop ; #  number of rows dropped in ssor0 
        PPINT       ssor1_free ; #  number of variables freed in ssor1 
        PPINT       ssor1_drop ; #  number of rows dropped in ssor1 
        PPINT       sparsa_col ; #  number of bound constraint changes in SpaRSA 
        PPINT       sparsa_row ; #  number of changes in row constraints in SpaRSA
        PPINT sparsa_step_fail ; #  number of failures of Armijo step in SpaRSA 
        PPINT             lnnz ; #  number of nonzeros in final Cholesky factor 
    
        #  timing 
        PPFLOAT    partition ; #  compute reordering of rows of A 
        PPFLOAT   initialize ; #  initwork and initlevels, includes partition 
        PPFLOAT       phase1 ; #  phase1 
        PPFLOAT       sparsa ; #  sparsa 
        PPFLOAT  coor_ascent ; #  coor_ascent 
        PPFLOAT        ssor0 ; #  ssor0 
        PPFLOAT        ssor1 ; #  ssor1 
        PPFLOAT         dasa ; #  dasa (includes coor_ascent, ssor0, and ssor1) 
        PPFLOAT    dasa_line ; #  dasa line search 
        PPFLOAT     checkerr ; #  check_error 
        PPFLOAT  prox_update ; #  prox_update 
        PPFLOAT       invert ; #  invert permutation of rows and columns 
        PPFLOAT       modrow ; #  modrow (update L by adding or deleting rows) 
        PPFLOAT       modcol ; #  modcol (rank 1 column updates of L) 
        PPFLOAT         chol ; #  cholmod_analyze, cholmod_factorize 
        PPFLOAT      cholinc ; #  incremental cholmod_rowfac 
        PPFLOAT     dltsolve ; #  dltsolve (back solve) 
        PPFLOAT       lsolve ; #  lsolve (forward solve) 

cdef extern from "cg_descent.h":
    # -------- cg functions -------- #
    void cg_print_parm (CGdata *Data)

    # -------- cg structs -------- #
    ctypedef struct CGdata:
        pass

    ctypedef struct CGparm:
        #  T => print status of run
        #  F => do not print status of run 
        int PrintStatus ;
    
        #  T => print pasa statistics
        #  F => do not print statistics 
        int PrintStat ;
    
        #  T => print parameter values 
        #  F => do not print parameter values 
        int PrintParm ;
    
        #  Level 0  = no printing, ... , Level 2 = maximum printing 
        int PrintLevel ;
    
        #  T => objective function is a quadratic
        #  F => objective function not necessarily quadratic 
        #int QuadCost ;
    
        #  an adjustment added to a quadratic cost objective when it is evaluated 
        CGFLOAT fadjust ;
    
        #  replace the Hessian Q of a quadratic by Q + QPshift 
        CGFLOAT QPshift ;
    
        #  stopping tolerance 
        CGFLOAT grad_tol ;

        # CG_DESCENT no longer includes the stopping criterion

        #    ||proj_grad||_infty <= grad_tol*(1 + |f_k|).

        # If the optimization problem is unconstrained, then the stopping
        # criterion is

        #    ||proj_grad||_infty <= testtol

        # where testtol = max(grad_tol,initial ||grad||_infty*StopFact) and
        # the default value of StopFact is zero. If the optimization problem
        # contains constraints, testtol is the pasa stopping criterion
        # switchfactor*global_error. This value for testtol is the PASA
        # criterion to stop solving the unconstrained problem and return
        # to the gradient project algorithm. */
        CGFLOAT StopFac ;
    
        #  T => check that f_k+1 - f_k <= debugtol*C_k
        #  F => no checking of function values 
        int debug ;
        CGFLOAT debugtol ;
    
        #  if step is nonzero, it is the initial step of the initial line search 
        CGFLOAT step ;
    
        #  0 => use cg_descent
        #  1 => use L-BFGS
        #  2 => use L-BFGS when LBFGSmemory >= n, use cg_descent when memory < n
        #  3 => use L-BFGS when LBFGSmemory >= n, use limited memory CG otherwise 
        int LBFGS ;
    
        #  if LBFGS is used, then LBFGSmemory is the number of vectors in memory 
        int LBFGSmemory ;
    
        #  abort cg after maxit iterations 
        CGINT maxit ;
    
        #  conjugate gradient method restarts after (n*restart_fac) iterations 
        CGFLOAT restart_fac ;
    
        #  factor in [0, 1] used to compute average cost magnitude C_k as follows:
        #  Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k 
        CGFLOAT Qdecay ;
    
        #  terminate after nslow iterations without strict improvement in
        #  either function value or gradient 
        int nslow ;
    
        #  factor by which eps grows when line search fails during contraction 
        CGFLOAT egrow ;
    
        #  T => attempt quadratic interpolation in line search when
        #           |f_k+1 - f_k|/f_k <= QuadCutoff
        #  F => no quadratic interpolation step 
        int    QuadStep ;
        CGFLOAT QuadCutOff ;
    
        #  maximum factor by which a quad step can reduce the step size 
        CGFLOAT QuadSafe ;
    
        CGFLOAT psi_lo ; #  in performing a QuadStep, we evaluate at point
                         #  betweeen [psi_lo, psi_hi]*psi2*previous step 
        CGFLOAT psi_hi ;
        CGFLOAT   psi1 ; #  for approximate quadratic, use gradient at
                         #  psi1*psi2*previous step for initial stepsize 
    
        CGFLOAT   qeps ; #  parameter in cost error for quadratic restart
                         #  criterion 
        CGFLOAT  qrule ; #  parameter used to decide if cost is quadratic 
        int   qrestart ; #  number of iterations the function should be
                         #  nearly quadratic before a restart 
    
        #  T => when possible, use a cubic step in the line search 
        int UseCubic ;
    
        #  use cubic step when |f_k+1 - f_k|/|f_k| > CubicCutOff 
        CGFLOAT CubicCutOff ;
    
        #  |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE
        CGFLOAT SmallCost ;
    
        #  maximum factor secant step increases stepsize in expansion phase 
        CGFLOAT ExpandSafe ;
    
        #  factor by which secant step is amplified during expansion phase
        #  where minimizer is bracketed 
        CGFLOAT SecantAmp ;
    
        #  factor by which rho grows during expansion phase where minimizer is
        #  bracketed 
        CGFLOAT RhoGrow ;
    
        #  If approxstep is TRUE, use approximate Wolfe line search.
        #  If approxstep is FALSE, then use ordinary line search and
        #  switch to the approximate step when |f - fnew| < ApproxSwitchFactor*|fR|,
        #  where fR is an estimate of the function size.  In the conjugate
        #  gradient code, an averaging technique is used to estimate function
        #  size. 
        int approxstep ;
        CGFLOAT ApproxSwitchFactor ;
    
        #  As the cg_descent converges, the function values typically
        #  approach a constant value. When this happens, the cubic interpolation
        #  step in the line search loses its accuracy and it is better to
        #  use a secant step based on the derivative of the objective function.
        #  The cost has converged when the relative change in the objective
        #  function <= CostConverge 
        CGFLOAT CostConverge ;
    
        CGFLOAT FuncGradSwitchFactor ;
    
        #  When performing an approximate Wolfe line search, we require
        #  that the new function value <= perturbation of the prior
        #  function value where the perturbation tries to take into
        #  rounding errors associated with the function value.
        #  There are two different perturbations:
        #  PertRule = 1 => fpert is f + Parm->pert_eps*|f| (relative change)
        #  PertRule = 0 => fpert is f + Parm->pert_eps     (absolute change) 
        int PertRule ;
        CGFLOAT pert_eps ;
    
        #  Maximum number of contractions in cg_contract. If it cannot find a
        #  step that either satisfies the Wolfe conditions or which has derivative
        #  >= 0 within ncontract attempts, then it is felt that fpert is too
        #  small and it will be increased so that the current function value
        #  is less than fpert. To increase fpert, we increase the value of
        #  pert_eps which is used to compute fpert. 
        int ncontract ;
    
        #  When pert_eps is increased, we multiply the new value by the growth
        #  factor eps_grow to ensure a healthy growth in pert_eps. 
        CGFLOAT eps_grow ;
    
        #  Maximum number of times that pert_eps is recomputed before a line
        #  search error is declared. 
        int neps ;
    
        CGFLOAT    cgdelta ; #  Wolfe line search parameter 
        CGFLOAT    cgsigma ; #  Wolfe line search parameter 
        int       maxsteps ; #  max number of tries to find acceptable step 
        CGFLOAT  stepdecay ; #  decay factor for bracket interval width 
        CGFLOAT        rho ; #  growth factor when searching for initial
                             #  bracketing interval 
        CGFLOAT       psi0 ; #  factor used in starting guess for iteration 1 
        CGFLOAT       psi2 ; #  when starting a new cg iteration, our initial
                             #  guess for the line search stepsize is
                             #  psi2*previous step 
        CGFLOAT  BetaLower ; #  parameter connected with lower bound for beta 
        CGFLOAT      theta ; #  parameter describing the cg_descent family 
        int   AdaptiveTheta ; #  T => choose theta adaptively, F => use theta 
    
        #  When an infinite or nan objective value is encountered, the
        #  stepsize is reduced in an effort to find a finite objective
        #  value. infdecay is the initial decay factor that is used when an
        #  infinite or nan objective value is encountered, ninf_tries is
        #  the number of attempts we make to find a finite objective value,
        #  and infdecay_rate is a factor by which infdecay is multiplied
        #  after each attampt to find a finite objective value. 
        CGFLOAT cg_infdecay ;
        CGFLOAT cg_infdecay_rate ;
        int  cg_ninf_tries ;
    
    #  ============ LIMITED MEMORY CG PARAMETERS ================================ 
        #  SubCheck and SubSkip control the frequency with which the subspace
        #  condition is checked. It it checked for SubCheck*mem iterations and
        #  if it is not activated, then it is skipped for Subskip*mem iterations
        #  and Subskip is doubled. Whenever the subspace condition is satisfied,
        #  SubSkip is returned to its original value. 
        int SubCheck ;
        int SubSkip ;
    
        #  when relative distance from current gradient to subspace <= eta0,
        #  enter subspace if subspace dimension = mem (eta0 = 0 means gradient
        #  inside subspace) 
    
        CGFLOAT   eta0 ; #  corresponds to eta0*eta0 in the paper 
    
        #  when relative distance from current gradient to subspace >= eta1,
        #  leave subspace (eta1 = 1 means gradient orthogonal to subspace) 
        CGFLOAT   eta1 ; #  corresponds to eta1*eta1 in the paper 
    
        #  when relative distance from current gradient to subspace <= eta2,
        #  always enter subspace (invariant space) 
        CGFLOAT   eta2 ;

    ctypedef struct CGstat:
        int         status ; #  returned status from cg_descent 
        CGFLOAT          f ; #  function value at solution 
        CGFLOAT        err ; #  sup norm of the gradient 
        CGFLOAT   grad_tol ; #  error tolerance 
        CGFLOAT      gnorm ; #  max abs component of gradient 
        CGFLOAT        tol ; #  computing tolerance: gnorm <= tol 
        CGINT        maxit ; #  maximum number of iterations 
        int  cg_ninf_tries ; #  number of tries to find finite objective value 
        CGFLOAT       oldf ; #  old function value when debugger fails 
        CGFLOAT       newf ; #  new function value when debugger fails 
        int       maxsteps ; #  max number of attempts in the line search 
        int        NegDiag ; #  T => negative diagonal encountered in QP factor 
        CGINT         iter ; #  number of iterations in cg 
        CGINT        nfunc ; #  function evaluations in cg 
        CGINT        ngrad ; #  gradient evaluations in cg 
        int        IterSub ; #  number subspace iterations in limited memory cg 
        int         NumSub ; #  number of subspaces in limited memory cg 

cdef extern from "napheap.h":
    # -------- napheap functions -------- #
    void napheap_print_parm (NAPdata *Data)

    # -------- napheap structs -------- #
    ctypedef struct NAPdata:
        pass

    ctypedef struct NAPparm:
        #  T => print status of run 
        int PrintStatus ;
    
        #  T => print statistics of the run 
        int PrintStat ;
    
        #  T => print the parameters 
        int PrintParm ;
    
        #  T => the prior_data input argument of napheap contains
        #       problem data from a previous run 
        int use_prior_data ;
    
        #  T => return the data structure for the current solution in the
        #       napcom argument of napheap (for use in another run with a
        #       different y value) 
        int return_data ;
    
        #  F => treat input argument lo as -infinity 
        int loExists ;
    
        #  F => treat input argument hi as +infinity 
        int hiExists ;
    
        #  F => no linear constraints bl <= a'x <= bu 
        int Aexists ;
    
        #  T => diagonal of objective Hessian is all positive, and the
        #       MATLAB expression all(d>0) is assumed to be true.
        #  F => there could be zeros on the diagonal of the Hessian 
        int d_is_pos ;
    
        #  T => diagonal of objective Hessian is zero.
        #       Note: The d argument of napheap is ignored. It could be NULL
        #  F => there could be nonzeros on the diagonal of the Hessian 
        int d_is_zero ;
    
        #  T => diagonal of objective Hessian is identically 1.
        #       Note: The d argument of napheap is ignored. It could be NULL 
        int d_is_one ;
    
        #  T => use input lambda as a starting guess for multiplier
        #  F => ignore input lambda (final lambda is returned in lambda argument) 
        #int use_lambda ;
    
        #  Upper bound on the number of Newton or variable fixing iterations.
        #  After performing these K iterations, the code switches to the break
        #  point search. The Newton algorithm also switches to the breakpoint
        #  search if two iterates are generated on opposite sides of the root. 
        int K ;
    
        #  T => use the Newton to update lambda, switch to the break point
        #  search either at iteration K or at the first time that two
        #  Newton iterates lie on opposite side of the root (whichever
        #  happens first)
        #  F => use K iterations of the variable fixing algorithm before
        #  switching to a breakpoint search 
        int newton ;
    
        #  An ordinary Newton iteration often converges monotonically to the
        #  root. Hence, the Newton iterates are scaled by a factor greater than 1
        #  in an effort to produce an iterate on the opposite side of the root.
        #  When such an iterate is generated, the algorithm switches to a
        #  breakpoint search. The scaled iteration is given by
        #               
        #  lambda_new = lambda_old - newton_scale*L'(lambda_old)/L''(lambda_old) 
        #               
        NAPFLOAT newton_scale ;
    
        #  To combat the effect of rounding errors, recompute a2sum whenever
        #  a2sum <= decay*a2max, the maximum a2sum encountered 
        double decay ;
    
        #  T => perform an iterative refinement step to improve solution accurary
        #  F => no refinement 
        int refine ;
    
        #  The solution is refined if |a'x - b| > err (only if refine is true) 
        NAPFLOAT err ;
    
        #  Check the input parameters for dj < 0 or hij < loj or Bhi < Blo 
        int check ;

    ctypedef struct NAPstat:
        NAPINT     nkf ;   #  number of known free variables 
        NAPINT   nfree ;   #  number of free variables in initial heap 
        NAPINT  nbound ;   #  number of bound variables in initial heap 
        NAPINT   nbrks ;   #  number break points to reach initial solution 
        NAPINT nrefine ;   #  number of break points during refinement 
        NAPINT nvarfix ;   #  number of variable fixing iterations 
        NAPINT nnewton ;   #  number of Newton iterations performed 
        NAPINT nsecant ;   #  number of secant steps performed 
    
        #  for error reporting 
        int     status ;   #  returned status from napheap 
        NAPINT  kerror ;   #  invalid n, or index where error found 
        NAPFLOAT lobad ;   #  bad value of lo, if any 
        NAPFLOAT hibad ;   #  bad value of hi, if any 
        NAPFLOAT  dbad ;   #  bad value of d, if any 
