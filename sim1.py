import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
import struct   # i

packints = 'ii'

# Set parameters
f_eps        = 0.0 # Carrier frequency offset percentage (0 = no offset)
Phase_Offset = 0.0 # Carrier phase offset (0 = no offset)
t_eps        = 0.0 # Clock freqency offset percentage (0 = no offset)
T_offset     = 0.0 # Clock phase offset (0 = no offset)
Ts = 1             # Symbol period
N = 4              # Samples per symbol period

fname = 'data/sim1_2024'

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
 
a_out = []
a_x = []
a_y = []

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
        a_out.append(sym)
        a_x.append(i_r[i])
        a_y.append(q_r[i])
        
#Find the UW
uw = np.array([162,29,92,47,16,112,63,234,50,7,15,211,109,124,
    239,255,243,134,119,40,134,158,182,0,101,
    62,176,152,228,36])

uw_len = uw.size
uw = uw.reshape(uw_len,)

uwsym = LUT[uw,:]
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

for i in range(len(uwrotsyms)):
    for j in range(len(uwrotsyms[i][1])):
        uw_rot_i = uwrotsyms[i][0][j]  # I component for all rotations
        uw_rot_q = uwrotsyms[i][1][j]  # Q component for all rotations
        
        for a in range(len(LUT)):
            if a == 0:
                dist = np.sqrt((LUT[a][1] - uw_rot_q)**2 + (LUT[a][0] - uw_rot_i)**2)
                sym = 0
            elif (np.sqrt((LUT[a][1] - uw_rot_q)**2 + (LUT[a][0] - uw_rot_i)**2) < dist):
                dist = np.sqrt((LUT[a][1] - uw_rot_q)**2 + (LUT[a][0] - uw_rot_i)**2)
                sym = a
        
        if j == 0:  # Append symbols based on rotation
            uw_rot_1_sym.append(sym)
        elif j == 1:
            uw_rot_2_sym.append(sym)
        elif j == 2:
            uw_rot_3_sym.append(sym)
        elif j == 3:
            uw_rot_4_sym.append(sym)
    
columns = 0    
rows = 0 
total_rows = 0   
for x in range(len(a_out)):

    if np.array_equal(uw_rot_1_sym[0:4], a_out[x:x+4]):
        columns = columns + 1
        if total_rows == 0:
            total_rows = rows
    elif np.array_equal(uw_rot_2_sym[0:4], a_out[x:x+4]):
        columns = columns + 1
        if total_rows == 0:
            total_rows = rows
    elif np.array_equal(uw_rot_3_sym[0:4], a_out[x:x+4]):
        columns = columns + 1
        if total_rows == 0:
            total_rows = rows
    elif np.array_equal(uw_rot_4_sym[0:4], a_out[x:x+4]):
        columns = columns + 1
        if total_rows == 0:
            total_rows = rows
    else:
        rows = rows + 1
        
total_symbols = (total_rows+6)*columns #+5 accounts for looking at the above 5

final_image = a_out [2*Lp:total_symbols+2*Lp]

final_image = np.reshape(final_image, (columns, total_rows+6)).T




#create 4 30 symbol arrays....
#Loop through symbols...
#compare to first 8

#rows = number of times until I get there
#columns = number of times I get a UW
plt.figure()
plt.plot(a_x, a_y, '.')
plt.title("Constellation Decision Variables")
plt.show()

breakpoint()
plt.figure()
plt.imshow(255-final_image,cmap=plt.get_cmap('Greys'))
plt.title('Figure 1')
plt.show()
breakpoint()

#output the image
        

        

