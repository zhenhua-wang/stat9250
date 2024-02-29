import os, ctypes
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate as si
from scipy import LowLevelCallable
from hw2.q1_integrand import integrand, integrand_transformed

lib = ctypes.CDLL(os.path.abspath('./hw2/q1.so'))

lib.integrand.restype = ctypes.c_double
lib.integrand.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double))
integrand_cpp = LowLevelCallable(lib.integrand)

def F(x, delta, tau, epsabs=1.49e-14, epsrel=1.49e-15):
    integral = si.quad(integrand, a=-np.inf, b=x-1,
                       args=(x, delta, tau),
                       epsabs=epsabs, epsrel=epsrel)
    return integral[0]

F_vec = np.vectorize(F)

X_vals1 = np.linspace(-100, 100, num=2000)
X_vals2 = np.linspace(1000, 10000, num=2000)
tau = np.sqrt(10)
delta1 = 0.3
delta2 = 0.7
plt.subplot(2, 2, 1)
plt.plot(X_vals1, F_vec(X_vals1, delta1, tau))
plt.subplot(2, 2, 2)
plt.plot(X_vals1, F_vec(X_vals1, delta2, tau))
plt.subplot(2, 2, 3)
plt.plot(X_vals2, F_vec(X_vals2, delta1, tau))
plt.subplot(2, 2, 4)
plt.plot(X_vals2, F_vec(X_vals2, delta2, tau))
plt.tight_layout()
plt.show()
