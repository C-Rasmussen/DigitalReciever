import math as math

import numpy as np
import matplotlib.pyplot as plt

def find_k1(BnT, zeta, K0 = 1, Kp = 1):
    rational = BnT/(zeta + 1/(4*zeta))
    num = rational * 4 * zeta
    den = 1 + 2*zeta*rational + rational**2
    answer = (num/den)/(Kp * K0)
    return answer
    
    
def find_k2(BnT, zeta, K0 = 1, Kp = 1):
    rational = (BnT/(zeta + (1/(4*zeta))))
    num = rational**2
    num = num * 4
    den = 1 + 2*zeta*rational + rational**2
    answer = (num/den)/(Kp * K0)
    return answer

def pll(theta, k1, k2, k0):
    Omega = 2 * np.pi / 10
    theta_stored = []
    theta_e_stored = []
    theta_hat_stored = []
    
    theta_hat = .1
    theta_e_k2_old = 0
    meu_n_past = 0
    for i in range(101):
        
        #theta_e = theta - theta_hat # phase detector
        #theta_e = math.sin(theta - theta_hat)
        theta_e = math.cos(Omega*i + theta) * (-math.sin(Omega*i+theta_hat))
        
        theta_e_k1 = theta_e*k1
        theta_e_k2 = theta_e*k2 + theta_e_k2_old
        
        theta_e_k2_old = theta_e_k2
        meu_n_future = theta_e_k2 + theta_e_k1
        
        theta_hat = theta_hat + k0*meu_n_past
        
        meu_n_past = meu_n_future #update variable
        
        
        theta_stored.append(theta)
        theta_e_stored.append(theta_e)
        theta_hat_stored.append(theta_hat)
    
    plt.figure(1); plt.clf()
    plt.ylim(-1.0, 3.0)
    plt.xlim(0, 100)
    plt.plot(range(101), theta_e_stored)
    plt.xlabel('n')
    plt.ylabel('e(n)')
    plt.show()
    breakpoint()
    
    Omega = 2 * np.pi / 10  # Angular frequency for a period of 20
    N = 101                  # Number of samples
    n = np.arange(N)

# Compute the cosine values with varying phases
    cos_values_output_pll = []
    angles_output_pll = []
    angles_input_pll = []
    for idx in range(101):
        cos_val = np.cos(Omega * n[idx] + theta_hat_stored[idx])
        cos_values_output_pll.append(cos_val)
        #angles_output_pll.append(.5857*Omega * n[idx] + theta_hat_stored[idx])
        #angles_input_pll.append(.5857*Omega * n[idx] + theta_stored[idx])
        angles_output_pll.append(Omega * n[idx] + theta_hat_stored[idx])
        angles_input_pll.append(Omega * n[idx] + theta_stored[idx])
        
    
    cos_values_input_pll = np.cos(Omega * n + theta)
    plt.figure(2); plt.clf()
    plt.plot(n, cos_values_input_pll)
    plt.plot(n, cos_values_output_pll)
    plt.xlabel('n')
    plt.ylabel('real part of sinusoids')
    plt.grid(True)
    plt.show()
    breakpoint()
    
    
        
    plt.figure(3); plt.clf()
    plt.plot(n, angles_input_pll)
    plt.plot(n, angles_output_pll)
    plt.xlabel('n')
    plt.ylabel('sinusoid arguments')
    #plt.title('')
    plt.grid(True)
    plt.show()
    breakpoint()


k1 = find_k1(0.05, 1/(np.sqrt(2)), 0.5)
k2 = find_k2(0.05, 1/(np.sqrt(2)), 0.5)
ret = pll(np.pi, k1, k2, 1)

# BnT = 0.05
# zeta = 1/(np.sqrt(2))

# print (find_k1(BnT=BnT, zeta=zeta))
# print (find_k2(BnT=BnT, zeta=zeta))