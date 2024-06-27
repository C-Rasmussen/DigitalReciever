import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
from pll import find_k1, find_k2
#Fixed Phase Offset

packints = 'ii'


# Set parameters
f_eps        = 0.0 # Carrier frequency offset percentage (0 = no offset)
Phase_Offset = 0.0 # Carrier phase offset (0 = no offset)
t_eps        = 0.0 # Clock freqency offset percentage (0 = no offset)
T_offset     = 0.0 # Clock phase offset (0 = no offset)
Ts = 1             # Symbol period
N = 4              # Samples per symbol period

fname = 'data/sim3_2024'

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

i_r = np.convolve(I, srrc)
q_r = np.convolve(Q, srrc)

#correlate to find the distance between each UW
correlate = []
for i in range(len(i_r)):
    if i%4 == 0: #downsample right here
        #slice right here
        for a in range(len(LUT)):
            if a == 0:
                dist = np.sqrt((LUT[a][1] - q_r[i])**2 + (LUT[a][0]- i_r[i])**2)
                sym = 0
            elif (np.sqrt((LUT[a][1] - q_r[i])**2 + (LUT[a][0] - i_r[i])**2) < dist):
                dist = np.sqrt((LUT[a][1] - q_r[i])**2 + (LUT[a][0] - i_r[i])**2)
                sym = a
        correlate.append(sym)
        
#test = np.correlate(correlate, correlate, mode='same')
total_symbols = len(correlate) - 4*Lp
rows = 100 + 30 #uw_len = 30, data rows = 100
cols = 0
#plt.figure()
#plt.plot(test)
#plt.show()
#breakpoint()
wait = 0
for sym in range(2*Lp, len(correlate)):
    if (sym % 130 == 0):
        new_uw = correlate[sym+130:sym+160]
        wait = wait + 1
        if wait == 2:
            break
                
uw = np.array([162,29,92,47,16,112,63,234,50,7,15,211,109,124,
    239,255,243,134,119,40,134,158,182,0,101,
    62,176,152,228,36])

uw_len = uw.size
uw = uw.reshape(uw_len,)

uwsym = LUT[uw,:]
new_uw_x_y = LUT[new_uw,:]
# Build the list of four possible UW rotations
angles = 2*math.pi*np.arange(0,4)/4
uwrotsyms = np.zeros((uw_len,2,4))
    

for i in range(angles.size):
  C = math.cos(angles[i])
  S = -math.sin(angles[i])
  G = np.array([[C, -S],[S, C]])
  uwrot = uwsym @ G;  # Rotate the UW symbols
  uwrotsyms[:,:,i] = uwrot; # Save the rotated version
  
uw_rot_1_sym = []
uw_rot_2_sym = []
uw_rot_3_sym = []
uw_rot_4_sym = []

phase90 = 0
phase180 = 0
phase270 = 0
phase0 = 0
for i in range(len(uwrotsyms)):
    for j in range(len(uwrotsyms[i][1])):
        uw_rot_i = uwrotsyms[i][0][j]  # I component for all rotations
        uw_rot_q = uwrotsyms[i][1][j]  # Q component for all rotations
        
        uw_new_i = new_uw_x_y[i][0]
        uw_new_j = new_uw_x_y[i][1]

        dist = np.sqrt((uw_new_j - uw_rot_q)**2 + (uw_new_i - uw_rot_i)**2)
                     
        if j == 0:  # Append symbols based on rotation
            phase0 = phase0 + dist
        elif j == 1:
            phase90 = phase90 + dist
        elif j == 2:
            phase180 = phase180 + dist
        elif j == 3:
            phase270 = phase270 + dist
            
if phase90 < phase0 and phase90 < phase270 and phase90 < phase180:
    Phase_Offset = np.radians(90.0)
elif phase0 < phase90 and phase0 < phase270 and phase0 < phase180:
    Phase_Offset = 0.0
elif phase180 < phase0 and phase180 < phase270 and phase180 < phase90:
    Phase_Offset = np.radians(180.0)
elif phase270 < phase0 and phase270 < phase90 and phase270 < phase180:
    Phase_Offset = np.radians(270.0)

#Have now gotten initial Phase Offset to get close before pll loop filter
   
a_out = []
BnT = 0.248
zeta = 1/(np.sqrt(2))
#find k1, k2
k1 = find_k1(BnT, zeta)
k2 = find_k2(BnT, zeta)
#k1 = 2.6*10**-2
#k2 = 6.9*10**-4
theta_arr = []
theta_e_k2_old = 0
theta_e_k2 = 0
meu_n_past = 0
a_hat_arr = []
e_k_arr = []
a_x = []
a_y = []
#loop through all downsampled match filter outputs
for i in range(len(i_r)):
       #Downsample
    if i%N == 0:     
         #pll for fixed phase offset
        C = np.cos(Phase_Offset)
        S = -np.sin(Phase_Offset)
        
        xr = C*i_r[i] - S*q_r[i]
        yr = S*i_r[i] + C*q_r[i]
        
            #slice right here
        for a in range(len(LUT)):
            if a == 0:
                dist = np.sqrt((LUT[a][1] - yr)**2 + (LUT[a][0]- xr)**2)
                sym = 0
            elif (np.sqrt((LUT[a][1] - yr)**2 + (LUT[a][0] - xr)**2) < dist):
                dist = np.sqrt((LUT[a][1] - yr)**2 + (LUT[a][0] - xr)**2)
                sym = a
        a_out.append(sym)
        if i > 200:
            a_x.append(xr)
            a_y.append(yr)
            
        a_hat_0 = LUT[sym][0]
        a_hat_1 = LUT[sym][1]
        a_hat_arr.append(sym)
            
        e_k = (yr*a_hat_0) - (xr*a_hat_1)
        e_k_arr.append(e_k)
        
        theta_e_k1 = e_k*k1
        theta_e_k2 = e_k*k2 + theta_e_k2_old
            
        theta_e_k2_old = theta_e_k2
        meu_n_future = theta_e_k2 + theta_e_k1
        
        Phase_Offset = Phase_Offset + meu_n_past
        theta_arr.append(Phase_Offset)
        meu_n_past = meu_n_future #update variable
        


for i in range(2*Lp,(len(a_out)-2*Lp)):
    if (i % 130 == 0): 
            cols = cols + 1
            
cols = cols - 1 #account for the first one which is not a column
 
total_symbols = int((rows)*cols) #+5 accounts for looking at the above 5

final_image = a_out [2*Lp:total_symbols+2*Lp]

final_image = np.reshape(final_image, (int(cols), int(rows))).T

plt.figure()
plt.plot(a_x, a_y, '.')
plt.title("Constellation Decision Variables")
plt.show()
breakpoint()

plt.figure(3)
plt.imshow(255-final_image,cmap=plt.get_cmap('Greys'))
plt.title('Figure 3')
plt.show()
breakpoint()



plt.figure()
plt.plot(e_k_arr)
plt.title("Error Signal")
plt.xlabel("Symbols index")
plt.ylabel("Error")
plt.show()
breakpoint()


plt.figure()
plt.plot(theta_arr)
plt.title("Phase Estimate")
plt.xlabel("Symbols index")
plt.ylabel("Phase")
plt.figure()
plt.show()
breakpoint()