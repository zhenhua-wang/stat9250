import numpy as np

def normal_density(eps, tau):
    return np.exp(-eps*eps/(2*tau*tau))/(tau*np.sqrt(2*np.pi))

def integrand(eps, x, delta, tau):
    return (1 - (delta / (2*delta - 1)) * (x - eps)**(-(1 - delta) / delta) +
            ((1 - delta) / (2*delta - 1)) * (x - eps)**-1) * \
            normal_density(eps, tau)

def integrand_transformed(t, x, delta, tau):
    log_t = np.log(t)
    return (1 - (delta / (2*delta - 1)) * (1 - log_t)**(-(1 - delta) / delta) +
            ((1 - delta) / (2*delta - 1)) * (1 - log_t)**-1) * \
            normal_density(log_t + x - 1, tau) / t
