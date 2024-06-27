import math
import numpy as np
from numpy.random import rand
from srrc1 import *
import matplotlib.pyplot as plt



alpha = .50 # excess bandwidth
N = 45# samples per symbol
Lp = 60 # SRRC truncation length
Ts = 1 # symbol time
T = Ts/N # sample time
srrcout = srrc1(alpha,N,Lp,Ts); # get the square-root raised cosine pulse
rcout = np.convolve(srrcout,srrcout)

Nsym = 400
#a = 1 # QPSK amplitude
a = sqrt(2) #QAM 16 amplitude

LUT_IQ = np.array([-3,-1,1,3])*a
#LUT_IQ = np.array([-1,1])*a
bits_i = np.random.randint(0, 4, Nsym)
bits_q = np.random.randint(0, 4, Nsym)
#bits_i = (rand(Nsym)> 0.5).astype(int) # generate random bits {0,1}
#bits_q = (rand(Nsym)> 0.5).astype(int)
ampa_i = LUT_IQ[bits_i] # map the bits to {+1,-1} values
ampa_q = LUT_IQ[bits_q]
#ampa = LUT_4PAM[bits]

upsampled_i = np.zeros((N*Nsym,1))
upsampled_i[range(0,N*Nsym,N)] = ampa_i.reshape(Nsym,1)

upsampled_q = np.zeros((N*Nsym,1))
upsampled_q[range(0,N*Nsym,N)] = ampa_q.reshape(Nsym,1)



s_i = np.convolve(upsampled_i.reshape((N*Nsym,)),srrcout)
s_q = np.convolve(upsampled_q.reshape((N*Nsym,)),srrcout)


#s_y = s_y.reshape(N*Nsym)
plt.figure(3); plt.clf()
plt.ylim(1.0, -1.0)
plt.xlim(-1.0, 1.0)
plt.plot(s_i, s_q)
plt.show()
test=1