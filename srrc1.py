from math import *
import numpy as np


def srrc1(alpha, N, Lp, Ts):
    # alpha = bandwidth excess
    # N = samples per symbol
    # Lp = truncation length
    # Ts = symbol time
    n = list(range(-Lp, Lp+1))
    for idx in range(len(n)):
        num = sin(pi * (1 - alpha) * (idx-Lp) / N) + 4 * alpha * (idx-Lp) / N * cos(pi * (1 + alpha) * (idx-Lp) / N)
        den = pi * (idx-Lp) / N * (1 - ((4 * alpha * (idx-Lp) / N) ** 2))
        
        if den == 0:
            n[idx] = sqrt(N)*(4*alpha/N-(pi*(alpha-1)/N))/(pi)
        else:
            n[idx] = ((1 / sqrt(N)) * (num / den))

    return n


def nrz(N, Lp):
    n = list(range(0, N+1))
    for idx in range(len(n)):
        shift = idx
        if shift >= 0 and shift <= N:
            n[idx] = 1/sqrt(N)
        else:
             n[idx] = 0
    return n


def manchester(N, Lp):
    
    n = list(range(0, N+1))
    for idx in range(len(n)):
        shift = idx
        if shift >= 0 and shift <= N/2:
            n[idx] = +sqrt(1/N)
        elif shift >= N/2 and shift <= N:
            n[idx] = -sqrt(1/N)
        else:
            n[idx] = 0
    return n

