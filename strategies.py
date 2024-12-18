import numpy as np

def f1(t):
    return 100*np.exp(-144*((600-t) / 3600)**2)*1000

def f2(t):
    if (600 - t) / 3600 > 0:
        return 88.42*1000
    else:
        return 0