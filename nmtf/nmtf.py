from __future__ import division, print_function

import warnings
import types

import numpy as np

from scipy.linalg import eig
from sklearn.base import BaseEstimator
from sklearn.utils import check_array, check_random_state
from sklearn.utils.extmath import squared_norm

from .utils import dispersion_coefficient_rho, dispersion_coefficients_eta_v, \
                    connectivity_matrix

def _init_svd(X_sum, k):
    """C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
    start for nonnegative matrix factorization, Pattern Recognition,
    Elsevier"""
    w, v = eig(X_sum)
    pos_w = np.abs(w)
    indices = np.argsort(w)[::-1]
    sorted_eigs = pos_w[indices]
    trace = np.sum(pos_w)
    tr_i = 0
    k_new = 1
    for i in range(len(w)):
        tr_i += sorted_eigs[i]
        if tr_i/trace > 0.9:
            k_new = i+1
            break

    k = min(k, k_new)

    G = []
    for i in range(k):
        xx = v[:, indices[i]]*pos_w[indices[i]]
        xp = np.maximum(xx, 0)
        xn = xp - xx
        if np.linalg.norm(xp) > np.linalg.norm(xn):
            G.append(xp)
        else:
            G.append(xn)
    G = np.array(G)
    G[np.where(G<1e-10)] = 1e-10
    return G.T


def _get_pos_neg(X):
    abs = np.abs(X)
    return (abs + X)/2., (abs - X)/2.


def _minimize_SSNMTF(X, G, A, L_neg, L_pos, weights=[], tol=1e-2, rel_tol=1e-8,
                     max_iter= 200, verbose=0, random_state = None, epsilon=1e-10,
                     compute_ktt=False, return_n_iter=False):
    X_norm = np.sum([squared_norm(x) for x in X])
    obj = np.inf
    epsilon = np.spacing(1)
    for iter_ in range(max_iter):
        GtG = G.T.dot(G) + epsilon*np.eye(G.shape[1])
        GtG_inv = np.linalg.inv(GtG)
        # update S
        S = [w*np.linalg.multi_dot((GtG_inv, G.T, x, G, GtG_inv))
             for w, x in zip(weights,X)]
        # update G
        Gn = np.zeros_like(G)
        Gd = np.zeros_like(G)
        for i in range(len(S)):
            RiGSi = X[i].dot(G).dot(S[i])
            RiGSi_p, RiGSi_n = _get_pos_neg(RiGSi)

            SiGtGSi = S[i].dot(GtG).dot(S[i])
            SiGtGSi_p, SiGtGSi_n = _get_pos_neg(SiGtGSi)

            GSiGtGSi_p = G.dot(SiGtGSi_p)
            GSiGtGSi_n =  G.dot(SiGtGSi_n)

            Gn += weights[i] * (RiGSi_p + GSiGtGSi_n)
            Gd += weights[i] * (RiGSi_n + GSiGtGSi_p + epsilon)
           
        for i in range(len(A)):
            Gn += self.gamma[i] * L_neg[i].dot(G)
            Gd += self.gamma[i] * L_pos[i].dot(G)

        G *= np.sqrt(Gn/Gd)

        #computing RSE
        rel_error = np.sum([weights[i]*squared_norm(X[i] - G.dot(S[i]).dot(G.T))
                           for i in range(len(X))])
        penalty = np.sum([gamma[i]*np.trace(G.T.dot(L_pos[i] - L_neg[i]).dot(G))
                         for i in range(len(A))])



        rse = rel_error/X_norm
        obj_old = obj
        obj = rel_error + penalty
        obj_diff = np.abs((obj_old - obj))/np.abs(obj_old) if iter_>0 else np.inf

        #KKT error
        if compute_ktt:
            err_1 = np.zeros_lik(G)
            err_2 = np.zeros_lik(S[0])
            GtG = G.T.dot(G)
            for i in range(len(X)):
                err_1 += X[i].dot(G).dot(S[i]) - G.dot(S[i]).dot(GtG).dot(S[i])
                err_2 += G.T.dot(X[i]).dot(G) -  GtG.dot(S[i]).dot(GtG)
            err_3 = np.conj(err_1).ravel().dot(G.ravel())
            KTT_err = np.maximum(
                      np.maximum(np.linalg.norm(err_1), np.linalg.norm(err_2)),
                      err_2)

        if verbose:
            if compute_ktt:
                print("iter: %d, obj:, %.4f, obj_diff: %.4f, rse: %4f, ktt:%.4f"
                    %(iter_, obj, rse, obj_diff, KTT_err))
            else:
                print("iter: %d, obj:, %.4f, obj_diff: %.4f, rse: %4f"
                    %(iter_, obj, obj_diff, rse))

        if rse <= tol or obj_diff <= rel_tol:
            break

    return_list = [G, S, rse]
    if return_n_iter:
        return_list.append(iter_)
    return return_list


class SSNMTF(BaseEstimator):
    """
    Function for solving Graph Regularized Symmetric Non-Negative Matrix Tri-Factorization

    Params
    ------

    init: string, default='SVD'
        - 'SVD':
        - 'random':
        if G_init is set is bypassed
    """

    def __init__(self, k, adjacencies=None, weights=[], gamma=0, max_iter=500, G_init=None,
                 init='SVD', epsilon=1e-10, compute_ktt=False, tol=1e-2, rtol=1e-8,
                 verbose=0, random_state=None):
        self.k = k
        self.adjacencies = adjacencies
        self.weights = weights
        self.gamma = gamma
        self.max_iter = max_iter
        self.G_init = G_init
        self.init = init
        self.verbose = verbose
        self.random_state = random_state
        self.epsilon = epsilon
        self.tol = tol
        self.rtol = rtol
        self.compute_ktt=compute_ktt


    def fit(self, X, y=None):
        """
        Params
        ------
        X: list of r nxn symmetric relationa matrices on a set of n objects
        y: None

        Returns
        -------
        G : array_like, shape=(n, k)
            The nonnegative clustering matrix shared over all Xs.
        S: list of r matrices, shape=(r, k, k)
            The latent matrices specific for each dataset in X.
        rse: float,
            Residual Square Error
        n_iter: int,
            Number of iteration for convergence
        history: list
            History of all the iterations.
        """
        if type(X) is np.ndarray:
            X = [check_array(X, ensure_min_features=2,
                         ensure_min_samples=2, estimator=self)]
        # TODO check symmetry of matrices in X
        X = [check_array(x, ensure_min_features=2,
                         ensure_min_samples=2, estimator=self) for x in X]
        assert len(np.unique([x.shape for x in X])) == 1, \
                "All the matrices in X should have the same order"
        if self.adjacencies is not None:
            A = [check_array(x, ensure_min_features=2,
                             ensure_min_samples=2, estimator=self)
                             for x in self.adjacencies]
            assert len(np.unique([x.shape for x in A])) == 1, \
                    "All the matrices in adjacencies hould have the same order"
        else:
            A = []

        self.random_state = check_random_state(self.random_state)
        
        if len(self.weights) != 0 and len(self.weights) != len(X):
            raise ValueError("You provided %d weights but the number of input "
                             "matrices is %d. You must provide a weight for"
                             "each matrix."%(len(self.weights), len(X)))
        elif len(self.weights) == 0:
            self.weights = list(np.ones(len(X)))
        # initialization
        if self.G_init is not None:
            G = G_init
        else:
            if self.init == 'random':
                G = self.random_state.rand(X[0].shape[0], self.k)
            elif str(self.init).upper() == 'SVD':
                G = _init_svd(np.sum([x*w for x, w in zip(X, self.weights)], 0), self.k)
            else:
                warnings.warn("No initialization specified,"
                              "initializating randomly")
                G = self.random_state.rand(X[0].shape[0], self.k)
        self.G_init = G.copy()
        
        # graph Regularization
        L_pos, L_neg = list(), list()
        for i in range(len(A)):
            L_pos.append(np.diag(np.sum(A[i], 1)))
            L_neg.append(A[i])

        self.G_, self.S_, self.reconstruction_error, self.n_iter_ = _minimize_SSNMTF(
            X, G, A, L_neg, L_pos, weights=self.weights, epsilon=self.epsilon, tol=self.tol,
             rel_tol=self.rtol, max_iter=self.max_iter,
             verbose= self.verbose, random_state = self.random_state,
             compute_ktt=self.compute_ktt, return_n_iter = True)

        return self


class SSNMTF_CV(BaseEstimator):
    """
    Cross Validation to determin parameter k of SSNMTF

    Params
    ------
    ks: int or list, optional, default=10
        If int 10 numbers are selected in an appropriate interval
        If list the number in the list are used as parameters.

    number_of_repetition: int, optional default=10
        The number of times to repeat the experiments in order to get the
        dispersion coefficient of the consensus connectivity matrix
    """

    def __init__(self, ks=10, adjacencies=None, gamma=0, max_iter=500,
                 mode="kim",
                 epsilon=1e-10, compute_ktt=False, tol=1e-2, rtol=1e-8,
                 verbose=0, random_state=None, number_of_repetition=10):
        self.ks = ks
        self.adjacencies = adjacencies
        self.gamma = gamma
        self.max_iter = max_iter
        self.verbose = verbose
        self.random_state = random_state
        self.epsilon = epsilon
        self.tol = tol
        self.rtol = rtol
        self.compute_ktt=compute_ktt
        self.number_of_repetition = number_of_repetition
        self.mode = mode

    def fit(self, X, y=None):
        """
        Params
        ------
        X: list of r nxn symmetric relationa matrices on a set of n objects
        y: None

        Returns
        -------
        G : array_like, shape=(n, k)
            The nonnegative clustering matrix shared over all Xs.
        S: list of r matrices, shape=(r, k, k)
            The latent matrices specific for each dataset in X.
        rse: float,
            Residual Square Error
        n_iter: int,
            Number of iteration for convergence
        history: list
            History of all the iterations.
        """

        # parameters to test
        if  isinstance(self.ks, int):
            k_max = int(np.sqrt(X[0].shape[0]))*2
            ks = np.arange(2, max(k_max, 3), max(int(round((k_max-2)/10)),1))
        else:
            ks = self.ks
        # cross Validation
        #if self.mode == 'kim':
        best_coeff = 0
        best_k = -1
        #if self.mode == 'dognig':
        best_eta = 0
        best_eta_k = -1
        best_v = 0
        best_v_k = -1

        results = dict.fromkeys(ks)
        self.input_ = X
        for k in ks:
            consensus = np.zeros_like(X[0])
            estimators = []
            mean_re = 0
            for rep_ in range(self.number_of_repetition):
                print(rep_)
        with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                est = SSNMTF(k, self.adjacencies, self.gamma, self.max_iter,
                             init='random', epsilon=self.epsilon,
                             compute_ktt=self.compute_ktt, tol=self.tol,
                             rtol=self.rtol, verbose=max(self.verbose-1,0),
                             random_state=self.random_state)
                est.fit(X)

                C = connectivity_matrix(est.G_)
                consensus += C
                estimators.append(est)
        mean_re += est.reconstruction_error
        consensus /= self.number_of_repetition
        mean_re /= self.number_of_repetition
        coeff = dispersion_coefficient_rho(consensus)
        eta, v = dispersion_coefficients_eta_v(consensus, k)
        results[k] = [estimators, consensus, mean_re, coeff, eta, v]
        if self.verbose:
            print("k: %d, dispersion_coefficient: eta %.4f, v %.4f"
                            %(k, eta, v)) 

        self.cv_results_ = results
        return self
