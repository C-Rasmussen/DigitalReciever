import numpy as np
import math as m
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

def filter_delay(data_arr, delay):
    delay_arr = np.zeros(delay)
    
    delay_arr[delay-1] = 1 #show impulse
    # plt.stem(delay_arr)
    # plt.title("Impulse Response of Delay of 100 Filter")
    # plt.show()
    
    
    
    impulse_response_padded = np.concatenate((delay_arr, np.zeros(100 * len(delay_arr))))
    frequency_response = np.abs(np.fft.fft(impulse_response_padded))
    frequency_axis = np.fft.fftfreq(len(delay_arr) + len(delay_arr)*100)
    # plt.plot(frequency_axis, frequency_response)
    # plt.xlabel('Frequency')
    # plt.ylabel('Magnitude')
    # plt.title('Magnitude Frequency Response')
    # plt.show()
    
    if delay == 0:
        return data_arr
    elif delay > 0:
        return np.concatenate((delay_arr, data_arr))
    
    
def differentiator_filter(length, data_arr, Ts):
    window = np.blackman(2*length+1)
    h_d = np.zeros(2*length+1)
    i = 0
    for idx in range (-length, length+1):
        
        
        #h_d[i] = Wc*Ts*np.cos(2*np.pi*Wc*Ts*idx)/(idx*np.pi*Ts) - 1/(np.pi*Ts)*(np.sin(Wc*Ts*idx))/(idx**2)
        if idx != 0:
            h_d[i] = (1/Ts) * (((-1)**idx)/idx)
        else:
            h_d[i] = 0
        i = i+1
       
    middle_index = len(h_d) // 2 + 1 
    new_h_d_real = np.zeros_like(h_d)
    new_h_d_real[middle_index - 2 : middle_index + 1] = h_d[middle_index - 2 : middle_index + 1]
   
    
    # plt.stem(h_d)
    # plt.title("Impulse Response of Differentiator Filter")
    # plt.show()
    
    h_real = h_d*window 
    #h_real = new_h_d_real*window
    
    
    #h_real = new_h_d_real
    impulse_response_padded = np.concatenate((h_real, np.zeros(100 * len(h_real))))
    
    frequency_response = np.abs(np.fft.fft(impulse_response_padded))
    frequency_axis = np.fft.fftfreq(len(h_real) + len(h_real)*100)
    #frequency_axis = np.fft.fftfreq(len(h_real))

    # plt.figure()
    # plt.plot(frequency_axis, frequency_response)
    # plt.xlabel('Frequency')
    # plt.ylabel('Magnitude')
    # plt.title('Magnitude Frequency Response')
    # plt.xlabel('Frequency')
    # plt.ylabel('Amplitude')
    # plt.show()
    
    frequency_axis = np.fft.fftfreq(len(np.angle(frequency_response)))
    
    plt.phase_spectrum(h_real, color ='green') 
    
    #plt.plot(frequency_axis, np.unwrap(np.angle(frequency_response)))
    plt.title('Phase of Magnitude Response, Linear Interpolator')
    plt.xlabel('Frequency [radians/sample]')
    plt.ylabel('Phase [radians]')
    plt.grid()
    plt.show()
    
    
    return np.convolve(h_real, data_arr)

def linear_interpolate(x_in, x_prev, T):
    #return T*x_in + (1-T)*x_prev
    return x_in*T + x_prev
    
    
def cubic_interpolate(x_n_2, m):
    #Following Figure 8.4.18 
   v3 = x_n_2[0]*1/6 - 1/2*x_n_2[-1] + 1/2*x_n_2[-2] -1/6*x_n_2[-3]
   v3 = v3 * m
   
   v2 = v3 + x_n_2[-1]*1/2  - x_n_2[-2] + x_n_2[-3]*1/2
   v2 = v2*m
   
   v1 = v2 - x_n_2[0]*1/6 + x_n_2[-1] -  x_n_2[-2]*1/2 - x_n_2[-3]*1/3
   v1 = v1*m
   
   return x_n_2[-2] + v1


# t = np.arange(0,10, 1/300) 
# Fc = 0.25
# input = np.sin(2*np.pi*Fc*t)
# cos = np.cos(2*np.pi*Fc*t)
# L = 6
# out = differentiator_filter(L, input, 1, np.pi)

# out = filter_delay(input, 100)



# plt.plot(input, label = "Pre Delayed Sinusoid")
# plt.plot(out, label = "Delayed Sinusoid")
# plt.title('Delay Filter')
# plt.legend()
# plt.show()

# num_coef = 2*L+1


#Differentiator
# plt.plot(t[:2700],input[:2700],label='Original Sine')
# plt.plot(t[:2700],(num_coef*3000)*out[num_coef-1:2710],label='Differentiated Sampled/Scaled Output') 
# plt.xlabel('Samples')
# plt.ylabel('Amplitude')
# plt.legend()
# plt.show()
# plt.grid()


# T = 1
# Fo = 0.15

# x = []
# #t = np.arange(0,11,1/25)
# for p in range(11):
#     x.append(np.sin(2*np.pi*Fo*p))
# #x = np.sin(2*np.pi*Fo*t)
# x_n = []
# for i in range(0,len(x)):
#     x_n.append(x[i])

# x_linear = []
# for i in range(1, len(x_n) ):
#     for j in range(10):
#         T = j/10
#         inter = linear_interpolate(x_n[i], x_n[i-1], T)
#         x_linear.append(inter)

# m = 0
# x_cubic = []
# for i in range(3,len(x_n)):
#     for j in range(10):
#         m = j/10
#         inter = cubic_interpolate(x_n[i-3:i],m)
#         x_cubic.append(inter)
    

#t1 = np.linspace(0,50,len(x_cubic))
#t2 = np.linspace(0,50,len(x_n))
# t3 = np.linspace(0,50,len(x_linear))


# t2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# t1 = np.linspace(0, 11, num=110)

# x_cubic = np.concatenate((x_cubic, np.zeros(20)))
# x_cubic = np.concatenate((np.zeros(10), x_cubic))


# x_linear = np.concatenate((x_linear, np.zeros(10) ))
# middle_values = x_linear
# middle_mask = (t1 >= 0.0) & (t1 <= 10.0) 
# y = np.ma.masked_array(middle_values, ~middle_mask)

# plt.figure()
# plt.plot(t1, y, label='Linear interpolated Function')
# plt.scatter(t2, x_n, label='Discrete Time True Function')
# plt.legend()
# plt.title('Linear Interpolation Result')
# plt.xlabel('Samples')
# plt.ylabel('Amplitude')
# plt.grid()
# plt.show()
# breakpoint()



# middle_values = x_cubic
# middle_mask = (t1 >= 1.0) & (t1 <= 9.0) 
# y = np.ma.masked_array(middle_values, ~middle_mask)

# plt.figure()
# plt.scatter(t2, x_n, label='Discrete Time True Function')
# plt.plot(t1, y,label='cubic interpolated')
# plt.legend()
# plt.title('Cubic Interpolation Result')
# plt.xlabel('Samples')
# plt.ylabel('Amplitude')
# plt.grid()
# plt.show()
# breakpoint()



