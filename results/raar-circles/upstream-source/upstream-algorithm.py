from numpy import exp

try:
    from tqdm import tqdm, tqdm_notebook
    load_tqdm = True
except ModuleNotFoundError:
    load_tqdm = False
    pass


class Algorithm:
    """
    Base class for ProxToolbox algorithms. Derived classes need
    to implement the evaluate() method which computes the next iterate. 
    The algorithm object is created when the Experiment object is initialized.
    """

    def __init__(self, experiment, iterateMonitor, accelerator=None):
        # instantiate prox operators
        self.proxOperators = []
        for proxClass in experiment.proxOperators:
            if proxClass != '':
                prox = proxClass(experiment)
                self.proxOperators.append(prox)
        # for convenience
        if len(self.proxOperators) >= 1:
            self.prox1 = self.proxOperators[0]
        if len(self.proxOperators) >= 2:
            self.prox2 = self.proxOperators[1]

        self.norm_data = experiment.norm_data
        self.Nx = experiment.Nx
        self.Ny = experiment.Ny
        self.Nz = experiment.Nz
        self.product_space_dimension = experiment.product_space_dimension
        self.truth = experiment.truth
        self.truth_dim = experiment.truth_dim
        self.norm_truth = experiment.norm_truth
        self.diagnostic = experiment.diagnostic

        self.maxIter = experiment.MAXIT
        self.tol = experiment.TOL
        self.lambda_0 = experiment.lambda_0
        self.lambda_max = experiment.lambda_max
        self.lambda_switch = experiment.lambda_switch

        self.iterateMonitor = iterateMonitor
        self.accelerator = accelerator

        if hasattr(experiment, 'progressbar'):
            self.progressbar = experiment.progressbar
        else:
            self.progressbar = None

        # temporaries for current iteration
        self.iter = 0
        self.u = None
        self.u_new = None

        # for debugging
        self.debug = experiment.debug

    def preprocess(self):
        """
        The default implementation calls the iterate monitor's
        preprocess() method which initializes the arrays that 
        will be used to collect statistics.      
        """
        self.iterateMonitor.preprocess(self)

    def stoppingCondition(self):
        """
        Evaluate stopping condition. The default stopping condition
        checks if the maximum number of iteration is reached or if 
        the given tolerance is attained. Can be overriden by derived
        classes.
 
        Returns
        -------
        bool
            Returns False when the algorithm must stop (the stopping
            condition is satisfied). Returns True otherwise.
        """
        return self.iter < self.maxIter and self.iterateMonitor.changes[self.iter] >= self.tol

    def evaluate(self, u):
        """
        Compute the algorithm update. Take the current iterate and returns 
        the next iterate. Must be overriden by derived class.
 
        Parameters
        ----------
        u: ndarray or a list of ndarray objects
            Current iterate.
        
        Returns
        -------
        u: ndarray or a list of ndarray objects
            Next iterate.
        """
        return u

    def updateStatistics(self):
        """
        Update statistics for the current iteration. Delegate
        the work to the iterate monitor. 
        """
        self.iterateMonitor.updateStatistics(self)

    def displayProgress(self):
        """
        Display progress. Delegate
        the work to the iterate monitor. 
        """
        self.iterateMonitor.displayProgress(self)

    def postprocess(self):
        """
        Prepare the data that will be returned by the run()
        method. Call the iterate monitor postprocess()
        method to retrieve the various statistics.
        
        Returns
        -------
        output: dictionary with the following entries:
            u: ndarray or a list of ndarray objects
                The last iterate.
            iter: int
                The number of iterations the algorithm performed

            additional entries are given by the iterate monitor

        """
        output = {'stats': {'iter': self.iter}, 'u': self.u}
        self.iterateMonitor.postprocess(self, output)
        return output

    def run(self, u):
        """
        Run the algorithm given an initial iterate u. The algorithm
        runs until the stopping condition is satisfied. Statistics are
        collected at each iteration by the iterate monitor. 

        Parameters
        ----------
        u: ndarray or a list of ndarray objects
            Initial iterate.
        
        Returns
        -------
        dictionary with the following entries:
            u: ndarray or a list of ndarray objects
                The last iterate.
            iter: natural number
                The last iteration count

            additional entries are given by the iterate monitor

        """
        self.u = u
        self.preprocess()
        self.displayProgress()

        if self.progressbar == 'tqdm_notebook':
            pbar = tqdm_notebook(total=self.maxIter)
        elif self.progressbar == 'tqdm':
            pbar = tqdm(total=self.maxIter)
        else:
            pbar = None

        while self.stoppingCondition():
            self.iter += 1
            if pbar is not None:
                pbar.update()
            self.u_new = self.evaluate(self.u)
            self.updateStatistics()
            self.displayProgress()
            self.u = self.u_new
        if pbar is not None:
            pbar.close()
        return self.postprocess()

    def getDescription(self):
        """
        Function returning a string giving the algorithm name
        and optionally a list of parameters that characterize
        this algorithm (i.e., `lambda_0`, and `lambda_max` constants).
        This string is used for graphical output.
        The default implementation only returns the name of the
        class implementing this algorithm. Derived classes may
        overrride this method to add meaningful parameters.
      
        Returns
        -------
        desc : string
            short description of the algorithm
        """
        desc = self.getDescriptionHelper()
        return desc

    def getDescriptionHelper(self, param_name=None, param_0=None, param_max=None):
        """
        Generate a string giving a short description of
        the algorithm (the name of its class name and optionally,
        the relaxation parameter used).
 
        Returns
        -------
        desc: string 
            description string genarated from the name of the 
            algorithm's class and the given parameters `param_0`
            and `param_max`.
        """
        desc = type(self).__name__
        # if param_name is not None and param_0 is not None \
        #        and param_max is not None:
            # desc += ", $" 
            # + param_name + "$ = " + str(param_0) \
            #        + " to " + str(param_max)
        return desc

    def computeRelaxation(self):
        """
        Compute the relaxation parameter based on the
        current iteration count.
 
        Returns
        -------
        lmbda: Float64 
            Relaxation parameter.
        """
        iteration = self.iter + 1  # add 1 to conform with matlab version
        # gradually changes lambda from lambda_0 into lambda_max
        a = -iteration / self.lambda_switch
        a *= a * a  # compute a^3
        s = exp(a)
        lmbda = s * self.lambda_0 + (1 - s) * self.lambda_max
        return lmbda
