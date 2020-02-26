
################### IMPORTS ####################################
import numpy as np
# %matplotlib inline
import matplotlib.pyplot as plt
import pylab
import math
################# DECLARE GLOBAL VARIABLES ###################
# f = 500
G = 6.67*10**(-11)
c = 3*10**8
# convert 10**30 kg to natural numbers
m_one = 10**30*9.109*10**(-31)
print(m_one)
m_two = 10**30*9.109*10**(-31)
# convert 356 ly to natural numbers
M = m_one+m_two
S_0 = 9*10**(-46)
f_s = 40
f_0 = 150
t = 100
psi = 100
x = f/f_0
nu = (m_one*m_two)/M**2
alpha = []
alpha.append(1)
alpha.append(0)
SNR = 10
F_lso = 1/(6**(3/2)*np.pi*M)
# modify the variables in params to whichever ones you want to be taken into
# account for the fishermatrix (from the waveform) If a parameter is taken out
# the first line of the wavefunction has to be adjusted to not unpack that
# nonexistent parameter
params = [psi, t, M, nu]

# change this function to adjust the wavefunction to be something different
def phi_f(f, params):
    psi, t, M, nu = params
    v = (np.pi*M*f)**(-1/3)
    sum_k = 0
    for k in range(len(alpha)):
        sum_k += alpha[k]*v**k
    return 2*np.pi*f*t-psi-np.pi/4+3/(128*nu*v**5)*sum_k

###################### CODE #######################

''' '''
#first find A (integral from low to high f^-7/6)/S(f)df
def get_A(f_s, f_0, low, high, SNR, accuracy_A):
    integral = 0
    for i in range(low*accuracy_A, high*accuracy_A):
        integral+=(i/accuracy_A)**(-7/6)/s_h(i/accuracy_A)
    return SNR/np.sqrt(abs(integral))


''' Calculates h(f)'''
def get_h_f(f, params):
    phi = phi_f(f,params)
    h_f = A*f**(-7/6)*np.exp(1j * phi)
    return h_f

''' Calculates S(h)'''
def s_h(f):
    x = f/f_0
    if(f>=f_s):
        r = S_0*(((4.49*x)**(-56))+0.16*(x**(-4.52))+0.52+(0.32*(x**2)))
    else:
        r = 10000000 # if infinity leads to numerical errors
    return(r)

''' Finds derivative h'(f)'''
def get_h_f_derivative(f, i):
    delta_x = params[i]/10
    new_params =  params.copy()
    new_params[i] = new_params[i]+params[i]/10
    delta_y = get_h_f(f, new_params) - get_h_f(f, params)
    return delta_y/delta_x


''' Creates matrix using previous functions'''
def create_fisher_matrix(low, high, f, increment):
    fisher_matrix = []
    for i in range(len(params)):
        fisher_matrix.append([])
        for j in range(len(params)):
            integral = 0
            for k in range(low*increment, high*increment):
                integral += np.real(get_h_f_derivative(f, i)*get_h_f_derivative(f, j)/s_h(k/increment))
            print(integral)
            fisher_matrix[-1].append(integral)
    return(fisher_matrix)

A = get_A(f_s, f_0, 10, 1000, SNR, 100)
print(A)
matrix = create_fisher_matrix(f_s, 1000, 100, 100)


# invert fisher_matrix
inverse = np.linalg.inv(matrix)
diag = np.diagonal(inverse)
print(np.sqrt(diag))

print(matrix[0])
print(matrix[1])
print(matrix[2])
print(matrix[3])




'''The following section plots Figure 1 of the paper for the initial ligo'''
graph_array = []
for i in range(0,1000):
    graph_array.append(np.sqrt(s_h(i)))
plt.figure()
plt.plot(graph_array)
plt.yscale('log')
plt.xscale('log')
axes = plt.gca()
axes.set_ylim([10**(-24), 10**(-21)])
axes.set_xlim([10**1, 10**3])
