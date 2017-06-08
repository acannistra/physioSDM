import numpy as np


def sigmoid_tmin(env_var, tmin):
    return 1.0 / (1.0 + np.exp(-env_var + tmin))


def sigmoid_tmax(env_var, tmin):
    return 1.0 / (1.0 + np.exp(env_var - tmin))


def thresh_fitness_min(env_var, tmin):
    """
    0.9 if env_var > $CT_{min} , 0.1 otherwise
    """
    if isinstance(tmin, int):
        tmin = np.repeat(tmin, len(env_var))
    return [0.9 if e > c else 0.1 for e, c in zip(env_var, tmin)]


def thresh_fitness_max(env_var, tmax):
    """
    0.9 if env_var < $CT_{max}$, 0.1 otherwise
    """
    if isinstance(tmax, int):
        tmin = np.repeat(tmax, len(env_var))
    return [0.9 if e < c else 0.1 for e, c in zip(env_var, tmin)]
