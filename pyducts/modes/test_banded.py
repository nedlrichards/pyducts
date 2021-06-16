import numpy as np
from scipy.linalg import solve, solveh_banded, solve_banded

numz = 105
d = np.full(numz, -1.9964864348874727)
C = np.diag(d) + np.diag(np.ones(numz - 1), k=1) + np.diag(np.ones(numz - 1), k=-1)

# first iteration to find eigen-vector
phi_test = np.ones(numz)
phi1 = solve(C, phi_test)

C_b = np.ones((3, numz))
C_b[1, :] = d
phi2 = solve_banded((1, 1), C_b, phi_test)

C_h = np.ones((2, numz))
C_h[1, :] = d
phi3 = solveh_banded(C_h, phi_test)
