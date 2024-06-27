# Digital communications simulation
#  Derived from Matlab code by Jake Gunther
# Date  :  April 1, 2024
# Class :  ECE 5660 (Utah State University)

import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
import struct   # i
import os

packints = 'ii'

# Set parameters
f_eps        = 0.0 # Carrier frequency offset percentage (0 = no offset)
Phase_Offset = 0.0 # Carrier phase offset (0 = no offset)
t_eps        = 0.0 # Clock freqency offset percentage (0 = no offset)
T_offset     = 0.0 # Clock phase offset (0 = no offset)
Ts = 1             # Symbol period
N = 4              # Samples per symbol period

fname = 'final_project/data/kitten1_8bit.dat'

# Select modulation type
# Use 8 bits per symbol and 256 square QAM
B = 8;            # Bits per symbol (B should be even: 8, 6, 4, 2)
# B = 4;
bits2index = 2**np.arange(B-1,-1,-1)
M = 2 ** B       # Number of symbols in the constellation
Mroot = math.floor(2**(B/2))
a = np.reshape(np.arange(-Mroot+1,Mroot,2),(2*B,1))
b = np.ones((Mroot,1))
LUT = np.hstack((np.kron(a,b), np.kron(b,a)))
# will be of the form (for example)
#  -3, -3
#  -3, -1
#  -3, 1
#  ...
#   3, 3
# of shape (B^2, 2)

# Scale the constellation to have unit energy
Enorm = np.sum(LUT ** 2) / M
LUT = LUT/math.sqrt(Enorm)

Eave = 1
Eb = Eave/B
EbN0dB = 30; # SNR in dB
N0 = Eb*10**(-EbN0dB/10)
nstd = math.sqrt(N0/2); # Noise standard deviation
# Note nstd is set to 0 below so there is no noise

if 1:
    plt.figure(1)
    # Plot the constellation
    plt.plot(LUT[:,0],LUT[:,1],'o')
    for i in range(0,M):
        plt.text(LUT[i,0]+0.02,LUT[i,1]+.02,i)

    # grid on; axis((max(axis)+0.1/B)*[-1 1 -1 1]); axis square;
    plt.xlabel('In-Phase')
    plt.ylabel('Quadrature')
    plt.title('Constellation Diagram')


# Unique word (randomly generated)
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
  
  
  
# Load and display the image
print(os.listdir('.'))
print('fname=',fname)

fid= open(fname,'rb')
twoints1 = fid.read(struct.calcsize(packints))
twoints = struct.unpack(packints,twoints1)
rows = twoints[0]
cols = twoints[1]
image = fid.read(rows*cols)
fid.close()
breakpoint()
image = np.frombuffer(image,dtype='uint8')
x = image.reshape(cols,rows).T

plt.figure(2)
plt.imshow(255-x,cmap=plt.get_cmap('Greys'))
plt.title('Original image')
plt.show()

print('rows=',rows,' cols= ',cols,' pixels=',rows*cols)



# Insert the unique word at the end of each column
x = np.vstack((x,np.matlib.repmat(uw.reshape(uw_len,1),1,cols)))
rows = rows + uw_len
plt.figure(3)
plt.imshow(255-x,cmap=plt.get_cmap('Greys'))
plt.title('With UW')
plt.show()

x = x.flatten('F')   # column scan
sym_stream = LUT[x,:]
sym_keep = sym_stream
num_syms = sym_stream.shape[0]

# Generate received signal with a clock frequency offset
print('Generating transmitted I/Q waveforms ... ')
EB = 0.7;  # Excess bandwidth
To = (1+t_eps)
if(t_eps == 0):  #  No clock skew
  Lp = 12
  t = np.arange(-Lp*N,Lp*N+1) /N + 1e-8;  # +1e-8 to avoid divide by zero
  tt = t + T_offset
  srrc = ((np.sin(math.pi*(1-EB)*tt)+ 4*EB*tt * np.cos(math.pi*(1+EB)*tt))
     /((math.pi*tt)*(1-(4*EB*tt)**2)))
  srrc = srrc/math.sqrt(N)
  Isu = np.zeros((num_syms*N,1))
  Isu[range(0,num_syms*N,N)] = sym_stream[:,0].reshape(num_syms,1)
  Qsu = np.zeros((num_syms*N,1))
  Qsu[range(0,num_syms*N,N)] = sym_stream[:,1].reshape(num_syms,1)
  I = np.convolve(srrc,Isu.reshape((N*num_syms,)))
  Q = np.convolve(srrc,Qsu.reshape((N*num_syms,)))
else: # Implement clock skew
  t = np.arange(0,num_syms*N)/N+1e-8
  I = np.zeros((1,t.length))  # In-phase pulse train
  Q = np.zeros((1,t.length))  # Quadrature pulse train
  for i in range(num_syms):
    tt = t-i*To + T_offset
    srrc = ((np.sin(math.pi*(1-EB)*tt)+4*EB*tt*np.cos(math.pi*(1+EB)*tt))
              /((math.pi*tt)*(1-(4*EB*tt)**2)))
    srrc = srrc/math.sqrt(N)
    I = I + srrc*sym_stream[i,0]
    Q = Q + srrc*sym_stream[i,1]
print('done.\n')
breakpoint()

# Modulate the pulse trains
print('Modulating I/Q waveforms ... ')
Omega0 = math.pi/2*(1+f_eps)
n = np.arange(I.size)
C =  math.sqrt(2)*np.cos(Omega0*n + Phase_Offset)
S = -math.sqrt(2)*np.sin(Omega0*n + Phase_Offset)
nstd = 0   # for this test, there is no noise
r = I * C + Q * S + nstd*np.random.normal(0,1,I.shape); # Noisy received signal
print('done.\n')

breakpoint()

