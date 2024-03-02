import numpy as np

def normal_density(eps, tau):
    return np.exp(-eps*eps/(2*tau*tau))/(tau*np.sqrt(2*np.pi))

def integrand(eps, x, delta, tau):
    return (1 - (delta / (2*delta - 1)) * (x - eps)**(-(1 - delta) / delta) +
            ((1 - delta) / (2*delta - 1)) * (x - eps)**-1) * \
            normal_density(eps, tau)

def integrand_transformed(t, x, delta, tau):
    t_inv = 1/t
    return (t_inv*t_inv - (delta / (2*delta - 1)) * t_inv**(1 / delta - 3) +
            ((1 - delta) / (2*delta - 1)) * t_inv) * \
            normal_density(x - t_inv, tau)
