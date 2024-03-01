import os, ctypes
import numpy as np
import matplotlib.pyplot as plt
from scipy import LowLevelCallable

lib_1 = ctypes.CDLL(os.path.abspath('./hw2/q1_integrand_cpp_gsl.so'))

lib_1.integrand.restype = ctypes.c_double
lib_1.integrand.argtypes = (ctypes.c_int, ctypes.c_void_p)

lib_1.f_gsl.restype = ctypes.c_double
lib_1.f_gsl.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double)

F_vec = np.vectorize(lib_1.f_gsl)

# plot
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
