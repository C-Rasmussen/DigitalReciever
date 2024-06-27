import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
from pll import find_k1, find_k2
from symbol_filters import differentiator_filter, cubic_interpolate, linear_interpolate

f_eps        = 0.0 # Carrier frequency offset percentage (0 = no offset)
Phase_Offset = 0.0 # Carrier phase offset (0 = no offset)
t_eps        = 0.0 # Clock freqency offset percentage (0 = no offset)
T_offset     = 0.0 # Clock phase offset (0 = no offset)
Ts = 1             # Symbol period
N = 4              # Samples per symbol period

fname = 'data/sim4_2024'

B = 8;            # Bits per symbol (B should be even: 8, 6, 4, 2)
# B = 4;
bits2index = 2**np.arange(B-1,-1,-1)
M = 2 ** B       # Number of symbols in the constellation
Mroot = math.floor(2**(B/2))
a = np.reshape(np.arange(-Mroot+1,Mroot,2),(2*B,1))
b = np.ones((Mroot,1))
LUT = np.hstack((np.kron(a,b), np.kron(b,a)))

Enorm = np.sum(LUT ** 2) / M
LUT = LUT/math.sqrt(Enorm)


fid= open(fname,'rb')
image = fid.read().splitlines()
for moodulated in range(len(image)):
    image[moodulated] = float(image[moodulated])
fid.close()

n = len(image)
n = np.arange(n)


cosine = np.sqrt(2)*np.cos(math.pi/2*(1+f_eps)*n) 
sine =  -np.sqrt(2)*np.sin(math.pi/2*(1+f_eps)*n)

I = image*cosine
Q = image*sine

EB = 0.7;  # Excess bandwidth
To = (1+t_eps)

Lp = 12
t = np.arange(-Lp*N,Lp*N+1) /N + 1e-8;  # +1e-8 to avoid divide by zero
tt = t + T_offset
srrc = ((np.sin(math.pi*(1-EB)*tt)+ 4*EB*tt * np.cos(math.pi*(1+EB)*tt))
    /((math.pi*tt)*(1-(4*EB*tt)**2)))

srrc = srrc/math.sqrt(N)

#match_filter
i_r = np.convolve(srrc, I)
q_r = np.convolve(srrc, Q)

i_dx_r = i_r
q_dx_r = q_r

L = 5
delay_impulse = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]

i_r = np.convolve(delay_impulse, i_r)
q_r = np.convolve(delay_impulse, q_r)

i_dx_r = differentiator_filter(L, i_dx_r, 100)
q_dx_r = differentiator_filter(L, q_dx_r, 100)

#TED init

e_arr = []
i_out = []
dx_i_out = []
q_out = []
dx_q_out = []
mu_k = []
a_hat_arr = []
i_inter_out = []
q_inter_out = []
W1 = 0
MU1 = 0
V21 = 0
NCO1 = 0
mu = 0

k1 = find_k1(0.0225, 1/(np.sqrt(2)), 0.5)
k2 = find_k2(0.0225, 1/(np.sqrt(2)), 0.5)

j = 3
c = 1


FX1 = 0; FX2 = 0; FX3 = 0; FX4 = 0; FX5 = 0;
FDX1 = 0; FDX2 = 0; FDX3 = 0; FDX4 = 0; FDX5 = 0;
FY1 = 0; FY2 = 0; FY3 = 0; FY4 = 0; FY5 = 0;
FDY1 = 0; FDY2 = 0; FDY3 = 0; FDY4 = 0; FDY5 = 0;
tempFx = 0; tempFdx = 0;
tempFy = 0; tempFdy = 0

#for n in range(0,N*(len(i_r))):
for n in range(0,(len(i_r))):
    
    temp = NCO1 - W1
    if temp < 0:
        strobe = 1
        mu = NCO1/w
        nco = 1 + temp
    else:
        strobe = 0
        mu = MU1
        nco = temp
    if strobe == 0:
        e = 0
    else:
        tempFx = -0.5*i_r[n]
        tempFdx = -0.5*i_dx_r[n]
        tempFy = -0.5*q_r[n]
        tempFdy = -0.5*q_dx_r[n]
        
        
        VF2 = -tempFx + FX1 + FX2 - FX3;
        VF1 = tempFx + (- FX1 + FX4) + FX2 + FX3;
        VF0 = FX5;
        i = (VF2*mu + VF1)*mu + VF0; # Interpolated signal
        
        VF2 = -tempFdx + FDX1 + FDX2 - FDX3;
        VF1 = tempFdx + (- FDX1 + FDX4) + FDX2 + FDX3;
        VF0 = FDX5;
        dx_i = (VF2*mu + VF1)*mu + VF0; #% Interpolated signal
        
        VF2 = -tempFy + FY1 + FY2 - FY3;
        VF1 = tempFy + (- FY1 + FY4) + FY2 + FY3;
        VF0 = FY5;
        q = (VF2*mu + VF1)*mu + VF0; # Interpolated signal
        
        VF2 = -tempFdy + FDY1 + FDY2 - FDY3;
        VF1 = tempFdy + (- FDY1 + FDY4) + FDY2 + FDY3;
        VF0 = FDY5;
        dx_q = (VF2*mu + VF1)*mu + VF0;# % Interpolated signal
        
        
        e = i*dx_i + q*dx_q
        e_arr.append(e)
        
        i_out.append(i)    
        dx_i_out.append(dx_i)
        q_out.append(q) 
        dx_q_out.append(dx_q)
        j = j + 1
        c = c + 1
        
    v1 = k1*e
    v2 = V21 + k2*e
    v = v1 + v2
    
    #NCO input
    w = v + 1/N
    
    
    FX3 = FX2;
    FX2 = FX1;
    FX1 = tempFx;
    FX5 = FX4;
    FX4 = i_r[n];
    FDX3 = FDX2;
    FDX2 = FDX1;
    FDX1 = tempFdx;
    FDX5 = FDX4;
    FDX4 = i_dx_r[n];
    
    FY3 = FY2;
    FY2 = FY1;
    FY1 = tempFy;
    FY5 = FY4;
    FY4 = q_r[n];
    FDY3 = FDY2;
    FDY2 = FDY1;
    FDY1 = tempFdy;
    FDY5 = FDY4;
    FDY4 = q_dx_r[n];

    V21 = v2
    NCO1 = nco
    MU1 = mu
    mu_k.append(mu)
    W1 = w
    
    
    
a_x = []
a_y = []
a_out = []
for sym in range(len(i_out)):
    #if sym%4 == 0: #downsample right here
        #slice right here
        for a in range(len(LUT)):
            if a == 0:
                dist = np.sqrt((LUT[a][1] - q_out[sym])**2 + (LUT[a][0]- i_out[sym])**2)
                symbol = 0
            elif (np.sqrt((LUT[a][1] - q_out[sym])**2 + (LUT[a][0] - i_out[sym])**2) < dist):
                dist = np.sqrt((LUT[a][1] - q_out[sym])**2 + (LUT[a][0] - i_out[sym])**2)
                symbol = a
        a_out.append(symbol)
        a_x.append(i_out[sym])
        a_y.append(q_out[sym])    
rows = 100+30       
columns = 0    
for p in range(2*Lp, len(a_out)):
    if p % 130 == 0:
        columns = columns + 1
        
#columns = columns - 1
          
total_symbols = rows*columns
#final_image = a_out [Lp+2:total_symbols+Lp+2]
#final_image = a_out [Lp+3:total_symbols+Lp+3]
final_image = a_out [2*Lp:total_symbols+2*Lp]
final_image = np.reshape(final_image, (int(columns), int(rows))).T

      

plt.figure(3)
plt.imshow(255-final_image,cmap=plt.get_cmap('Greys'))
plt.title('Figure 4')
plt.show()
breakpoint()    


plt.figure()    
plt.plot(a_x,a_y,'.')
plt.title('Constellation Decision Variables')
plt.show()
breakpoint()

plt.figure()
plt.plot(mu_k,'.')
plt.title('Fractional Interval')
plt.xlabel("Symbol Index")
plt.ylabel("Fractional Interval u(k)")
plt.show()
breakpoint()

plt.figure()
plt.plot(e_arr)
plt.title("Timing Error Signal")
plt.xlabel("Symbols index")
plt.ylabel("Error e(t)")
plt.ylim(-.015, .015)
plt.show()
breakpoint()

            
        
            

