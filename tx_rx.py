import math
import numpy as np
from numpy.random import rand
from srrc1 import *
import matplotlib.pyplot as plt

def tx_rx(input_bits, Ts, Ns, alpha, fc, lut_i, lut_q):
    #Ns = samples per symbol
    #Ts = symbol period
    #input_bits = input_signal
    #alpha = excess bandwidth
    #fc = carrier frequency
    #lut_i/q = array to map input bits to , are np arrays
    Lp=60
    ###Put whatever p(t) function is desired right here....srcc default
    p_t = srrc1(alpha,Ns,Lp,Ts)
    r_t = np.convolve(p_t,p_t)
    
     #number of symbols calculated
        
    bits_2d = input_bits.reshape(-1, 2)
    print(f'input bits: {input_bits[0:20]}\n') 

    # Access each bit separately
    bit_0 = bits_2d[:, 0]  # First column
    # bit_1 = bits_2d[:, 1]  # Second column

    bits_str = np.apply_along_axis(lambda x: ''.join(map(str, x)), axis=1, arr=bits_2d)


    integer_values = np.array([int(bits, 2) for bits in bits_str])
    

    Nsym = len(bit_0)
    
    # amplitude_i = lut_i[bit_0]
    # amplitude_q = lut_q[bit_1]
    amplitude_i = lut_i[integer_values]
    amplitude_q = lut_q[integer_values]
    
    upsampled_i = np.zeros((Ns*Nsym,1))
    upsampled_i[range(0,Ns*Nsym,Ns)] = amplitude_i.reshape(Nsym,1)
    
    
    upsampled_q = np.zeros((Ns*Nsym,1))
    upsampled_q[range(0,Ns*Nsym,Ns)] = amplitude_q.reshape(Nsym,1)
    
    #pulse shape filter
    s_i = np.convolve(upsampled_i.reshape((Ns*Nsym,)),p_t)
    s_q = np.convolve(upsampled_q.reshape((Ns*Nsym,)),p_t)
    
    plt.figure(5); plt.clf()
    fft_out = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(np.correlate(s_i,s_i,'full')))))
    fft_freq = np.linspace(-Ns/2, Ns/2, (2 * len(s_i))-1, endpoint=False)
    plt.plot(fft_freq, fft_out)
    plt.show()
    
    
    #Phase trajectory Plot
    plt.figure(1); plt.clf()
    plt.ylim(1.0, -1.0)
    plt.xlim(-1.0, 1.0)
    plt.plot(s_i, s_q)
    plt.show()
    breakpoint()
    
    n = np.arange(0, len(s_i))  # Time index

# Calculate the modulated signals
    s_i = np.sqrt(2)*s_i*np.cos(2*np.pi*fc*n)
    s_q = -np.sqrt(2)*s_q*np.sin(2*np.pi*fc*n)
    

    r_t = s_i + s_q #Also s(t)
    
    plt.figure(6); plt.clf()
    
    fft_out = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(np.correlate(r_t,r_t,'full')))))
    fft_freq = np.linspace(-Ns/2, Ns/2, (2 * len(r_t))-1, endpoint=False)
    plt.plot(fft_freq, fft_out)
    #plt.plot(10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(np.correlate(r_t,r_t,'full'))))**2))
    plt.show()
    
    r_t_i = r_t*np.sqrt(2)*np.cos(2*pi*fc*n) #projection into signal space
    r_t_q = -r_t*np.sqrt(2)*np.sin(2*pi*fc*n)
    
    r_t_i = np.convolve(r_t_i, p_t)
    r_t_q = np.convolve(r_t_q, p_t)
    
    plt.figure(7); plt.clf()
    fft_out = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(np.correlate(r_t_i,r_t_i,'full'))))**2)
    fft_freq = np.linspace(-Ns/2, Ns/2, (2 * len(r_t_i))-1, endpoint=False)
    plt.plot(fft_freq, fft_out)
    plt.show()
    
    plt.figure(8); plt.clf()
    fft_out = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(np.correlate(r_t_q,r_t_q,'full'))))**2)
    fft_freq = np.linspace(-Ns/2, Ns/2, (2 * len(r_t_q))-1, endpoint=False)
    plt.plot(fft_freq, fft_out)
    plt.show()
    
    #Match filter output
    offset = (2*Lp - np.floor(Ns/2)).astype(int)
    # ˆ ˆ
    # 1st correlation |
    # peak |
    # move to center
    Nsymtoss = 2*np.ceil(Lp/Ns) # throw away symbols at the end
    nc = (np.floor((len(r_t_i) - offset - Nsymtoss*Ns)/Ns)).astype(int) # number of points of signal t
    xreshape = r_t_i[offset:offset + nc*Ns].reshape(nc,Ns)
    yreshape = r_t_q[offset:offset + nc*Ns].reshape(nc,Ns)
    
    plt.figure(2); plt.clf()
    plt.plot(np.arange(-np.floor(Ns/2), np.floor(Ns/2)+1), xreshape.T,color='b',
    linewidth=0.25)
    plt.show()
    breakpoint()
    
    plt.figure(3); plt.clf()
    plt.plot(np.arange(-np.floor(Ns/2), np.floor(Ns/2)+1), yreshape.T,color='b',
    linewidth=0.25)
    plt.show()
    breakpoint()
    
    
    plt.figure(4); plt.clf()
    plt.scatter(xreshape[:, Ns//2], yreshape[:, Ns//2])
    plt.show()
    output = []
    for sample in range(len(xreshape)):
        if xreshape[sample][Ns//2] > 0 and yreshape[sample][Ns//2] > 0:
            output.append([1, 1])
        elif xreshape[sample][Ns//2] > 0 and yreshape[sample][Ns//2] < 0:
            output.append([1, 0])
        elif xreshape[sample][Ns//2] < 0 and yreshape[sample][Ns//2] > 0:
            output.append([0, 1])
        elif xreshape[sample][Ns//2] < 0 and yreshape[sample][Ns//2] < 0:
            output.append([0, 0])

    print(f'found bits: {output[0:20]}\n')               