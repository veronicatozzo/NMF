import numpy as np


# normalizzando le S
# facendo una media pesata 
def integration_SSNMTF(G, S, mode='mean'):
    if str(mode).lower() == 'mean':
        integrated = np.mean(S, axis=0)
    else:
        raise ValueError("You specified a not availale mode")

    return G.dot(integrated).dot(G.T)
