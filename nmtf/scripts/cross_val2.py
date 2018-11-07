import matlab.engine
import numpy as np
eng = matlab.engine.start_matlab()
from os import listdir
from os.path import isfile, join

from nmtf.read import get_adjacency
from nmtf.utils import connectivity_matrix, dispersion_coefficient_rho, \
                        dispersion_coefficients_eta_v


folder = "/Users/veronica/Dropbox (DIBRIS)/project_UCL/BCStages-SubtypesNetworks/HER2Networks"
files = [join(folder, f) for f in listdir(folder) if isfile(join(folder, f))]
graphs = [get_adjacency(f) for f in files]


ll = [matlab.double(l.tolist()) for l in graphs]
ks = np.arange(5, 55, 10)
n_repetitions = 1

res = dict.fromkeys(ks)
for k in ks:
    consensus = np.zeros((1000,1000))
    for rep_ in range(n_repetitions):
        S, G, rse = eng.SSNMTF(ll, [], k, 0, 200, 1, nargout=3)
        conn_matrix = connectivity_matrix(np.real(np.array(G)))
        consensus += conn_matrix
    consensus /= n_repetitions
    rho = dispersion_coefficient_rho(consensus)
    eta, v = dispersion_coefficients_eta_v(consensus)
    res[k] = [consensus, rho, eta, v]
