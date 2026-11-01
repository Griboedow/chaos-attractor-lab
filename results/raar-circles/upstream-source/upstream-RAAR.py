from proxtoolbox.algorithms.algorithm import Algorithm
from proxtoolbox.utils.cell import Cell, isCell
from numpy import exp


class RAAR(Algorithm):
    """
    User-friendly version of the Relaxed Averaged Alternating
    Reflection algorithm. For background and analysis see:
    D.R.Luke, Inverse Problems 21:37-50(2005)
    D.R. Luke, SIAM J. Opt. 19(2): 714--739 (2008).
    D. R. Luke and A.-L. Martins, "Convergence Analysis 
    of the Relaxed Douglas Rachford Algorithm" , 
    SIAM J. on Optimization, 30(1):542--584(2020).
    https://epubs.siam.org/doi/abs/10.1137/18M1229638

    Based on Matlab code written by Russell Luke (Inst. Fuer 
    Numerische und Angewandte Mathematik, Universitaet
    Gottingen) on Aug 18, 2017.    
    """

    def evaluate(self, u):
        """
        Update for the Relaxed Averaged Alternating
        Reflection algorithm. 

        Parameters
        ----------
        u : ndarray or a list of ndarray objects
            The current iterate.
               
        Returns
        -------
        u_new : ndarray or a list of ndarray objects
            (same type as input parameter u)
            The new iterate.
        """
        lmbd = self.computeRelaxation()
        if not isCell(u):
            tmp1 = 2*self.prox2.eval(u) - u  # Reflection using prox2 on u
            tmp2 = self.prox1.eval(tmp1)     # Projection of the reflection
            # update
            u_new = (lmbd*(2*tmp2 - tmp1) + (1 - lmbd)*tmp1 + u)/2
        else:
            u_new = self.prox2.eval(u)
            tmp1 = Cell(len(u))
            for j in range(len(u)):
                tmp1[j] = 2*u_new[j] - u[j]
            tmp2 = self.prox1.eval(tmp1)
            # update
            for j in range(len(u)):
                u_new[j] = (lmbd*(2*tmp2[j]-tmp1[j]) + (1-lmbd)*tmp1[j] + u[j])/2
        return u_new

    def getDescription(self):
        return self.getDescriptionHelper("\\lambda", self.lambda_0, self.lambda_max)
