import numpy as np
from scipy.optimize import curve_fit

def langmuir(p, qmax, K):
    return (qmax * K * p) / (1 + K * p)

p_data = np.array([0,1.048e17, 2.501e17, 3.41e17,		4.64e17	])
q_data = np.array([0,1.048e17 , 2.28e17, 3.2895e17, 4.365e17])

params, cov = curve_fit(langmuir, p_data, q_data)

qmax, K = params
print(qmax, K)
