# Use HF to solve the molecule with 2 atoms
# Use Gaussian basis sets
# Use RHF, consider only closed-shell molecules
# Use atomic units
# We first try STO-3G basis set

# Steps (Szabo's book, p146):
# 1. define a diatomic molecule, with the distance between the atoms, the nuclear charge z1 and z2, the number of electrons num_e and basis set
# 2. calculate molecular integrals
# 3. diagonalize the overlap matrix S and obtain transformation matrix X
# 4. obtain a random guess at the density matrix P
# 5. calculate the matrix G from P
# 6. obtain F = H_core + G
# 7. calculate F' = X^{dagger} * F * X
# 8. diagonlaize F' to obatin C' and epsilon
# 9. calculate C = XC'
# 10. form a new P from C
# 11. determine whether the procedure has converged
# 12. if converged, use the resultant solution to calculate the expectation values and other quantities of interest

import numpy as np
from numpy import linalg as LA
from scipy import special

# 1. define the diatomic molecule
# Use hydrogen molecule as my first attempt
z1,z2 = 1,1 # The charge of the two nuclei
num_e = 2 # Closed-shell molecule, so the number of electrons must be even numbers
R = 1.4 # we place atom 1 at x = 0 and atom 2 at x = R (this is a 1D system)
z_list = [z1,z2] # Save the nucleus charge in this list
R_list = [0,R] # Save the coordinate of the nucleus in this list
len_CGF = 3 # The length of CGF is 3, that is to say 3 primitive gaussian functions are contracted to a CGF

# Here we want to use STO-3G to use 3 Gaussian functions to fit a Slater function
# Refer to Szabo's book p157-159, using the scaling tactics
# All the GFs in this program are the normalized gaussian functions!!!
# The basis function we will use: CGF(zeta, STO-3G) = 0.444635 * GF(zeta ** 2 * 0.109818) + 0.535328 * GF(zeta ** 2 * 0.405771) + 0.154329 * GF(zeta ** 2 * 2.22766)
# As for H atom, zeta in CGF is chosen as 1.24

zeta = 1.24
coefficient = [0.444635,0.535328,0.154329]

# Here, we represent each GF as a list with 2 numbers, represent two parameters of GF, that is alpha and R_A (R_A is the coordinate of the center nucleus)
# CGF_1 = 0.444635 * GF_1_1 + 0.535328 * GF_1_2 + 0.154329 * GF_1_3
# CGF_2 = 0.444635 * GF_2_1 + 0.535328 * GF_2_2 + 0.154329 * GF_2_3
# GF_x_x = [alpha,x_0] is the standard GF with the expression of (2 * alpha / np.pi) ** (3 / 4) * np.exp(-1 * alpha * (r - x_0) ** 2)

GF = [[None for i in range(3)] for j in range(2)]
GF[0][0] = [zeta ** 2 * 0.109818,0] # format: [alpha, x_coordinate_of_the_nucleus]
GF[0][1] = [zeta ** 2 * 0.405771,0] # format: [alpha, x_coordinate_of_the_nucleus]
GF[0][2] = [zeta ** 2 * 2.22766,0] # format: [alpha, x_coordinate_of_the_nucleus]
GF[1][0] = [zeta ** 2 * 0.109818,R] # format: [alpha, x_coordinate_of_the_nucleus]
GF[1][1] = [zeta ** 2 * 0.405771,R] # format: [alpha, x_coordinate_of_the_nucleus]
GF[1][2] = [zeta ** 2 * 2.22766,R] # format: [alpha, x_coordinate_of_the_nucleus]

# print(GF)

# Then we have do define several integrals

# First, we define the integral related to overlap matrix

def S_GF(para_1,para_2): # S_GF(para_1,para_2) is the overlap integral of GFs, para_1, para_2 is the parameters of the 2 input GFs
    alpha_1,alpha_2 = para_1[0],para_2[0]
    x_1,x_2 = para_1[1],para_2[1]
    return (4 * alpha_1 * alpha_2 / (alpha_1 + alpha_2) ** 2) ** (3 / 4) * np.exp(-1 * alpha_1 * alpha_2 * (x_1 - x_2) ** 2 / (alpha_1 + alpha_2)) # checked

# Next, we define the kinetic energy integrals

def T_GF(para_1,para_2): # T(para_1,para_2) is the kinetic energy integral of GFS, para_1, para_2 is the parameters of the 2 input GFs
    alpha_1,alpha_2 = para_1[0],para_2[0]
    x_1,x_2 = para_1[1],para_2[1]
    K = np.exp(-1 * alpha_1 * alpha_2 * (x_1 - x_2) ** 2 / (alpha_1 + alpha_2))
    p = alpha_1 + alpha_2
    return 2 ** (3 / 2) * (alpha_1 * alpha_2) ** (7 / 4) * K * p ** (-5 / 2) * (3 - 2 * alpha_1 * alpha_2 * (x_1 - x_2) ** 2 / p) # checked

# Next, we define the nuclear attraction integral

def V_one_nucl_GF(para_1,para_2,Z_C,R_C): # Z_C is the charge of the nucleus, while R_C is the coordinate (here only the x coordinate) of the nucleus
    alpha_1,alpha_2 = para_1[0],para_2[0]
    x_1,x_2 = para_1[1],para_2[1]
    K = np.exp(-1 * alpha_1 * alpha_2 * (x_1 - x_2) ** 2 / (alpha_1 + alpha_2))
    p = alpha_1 + alpha_2
    R_P = (alpha_1 * x_1 + alpha_2 * x_2) / (alpha_1 + alpha_2)
    R_Q = abs(R_P - R_C) # R_Q is always positive
    if R_Q > 10 ** (-6):
        return -1 * (4 * alpha_1 * alpha_2 / p ** 2) ** (3 / 4) * Z_C * K * special.erf(R_Q * p ** 0.5) / R_Q
    else:
        return -1 * (4 * alpha_1 * alpha_2) ** (3 / 4) * 2 * Z_C * K / (p * np.pi ** 0.5)

# Finally, we define the two-electron integral
def V_ee_GF(para_1,para_2,para_3,para_4):
    alpha_1,alpha_2,alpha_3,alpha_4 = para_1[0],para_2[0],para_3[0],para_4[0]
    x_1,x_2,x_3,x_4 = para_1[1],para_2[1],para_3[1],para_4[1]
    p_1 = alpha_1 + alpha_2
    p_2 = alpha_3 + alpha_4
    xi = (p_1 + p_2) / (4 * p_1 * p_2)
    K_1 = np.exp(-1 * alpha_1 * alpha_2 * (x_1 - x_2) ** 2 / (alpha_1 + alpha_2))
    K_2 = np.exp(-1 * alpha_3 * alpha_4 * (x_3 - x_4) ** 2 / (alpha_3 + alpha_4))
    R_P_1 = (alpha_1 * x_1 + alpha_2 * x_2) / (alpha_1 + alpha_2)
    R_P_2 = (alpha_3 * x_3 + alpha_4 * x_4) / (alpha_3 + alpha_4)
    R_L = abs(R_P_1 - R_P_2)
    if R_L > 10 ** (-6):
        return (64 / R_L) * K_1 * K_2 * (alpha_1 * alpha_2 * alpha_3 * alpha_4) ** (3 / 4) * (4 * p_1 * p_2) ** (-3 / 2) * special.erf(R_L / (2 * xi ** 0.5))
    else:
        return 64 * K_1 * K_2 * (alpha_1 * alpha_2 * alpha_3 * alpha_4) ** (3 / 4) * (4 * p_1 * p_2) ** (-3 / 2) * (np.pi * xi) ** (-1 / 2)

# Calculate matrix S
def calculate_matrix_S_element(m,n): # calculate the (m,n) element of matrix S
    output = 0
    for i in range(len_CGF):
        for j in range(len_CGF):
            output += coefficient[i] * coefficient[j] * S_GF(GF[m][i],GF[n][j])
    return output

S = np.zeros((2,2))
for m in range(2):
    for n in range(2):
        S[m][n] = calculate_matrix_S_element(m,n)
print('S = ',S) # checked
# print('type(S) = ',type(S))

# Calculate matrix T
def calculate_matrix_T_element(m,n):
    output = 0
    for i in range(len_CGF):
        for j in range(len_CGF):
            output += coefficient[i] * coefficient[j] * T_GF(GF[m][i],GF[n][j])
    return output

T = np.zeros((2,2))
for m in range(2):
    for n in range(2):
        T[m][n] = calculate_matrix_T_element(m,n)
print('T = ',T) # checked
# print('type(T) = ',type(T))

# Calculate matrix V_nucl (the energy between one nucleus and one electron, and all the nuclei are considered)
def calculate_matrix_V_nucl_element(m,n):
    output = 0
    for i in range(len(z_list)): # consider all the nuclei
        for j in range(len_CGF):
            for k in range(len_CGF):
                output += coefficient[j] * coefficient[k] * V_one_nucl_GF(GF[m][j],GF[n][k],z_list[i],R_list[i])
    return output

V_nucl = np.zeros((2,2))
for m in range(2):
    for n in range(2):
        # print(m,n)
        V_nucl[m][n] = calculate_matrix_V_nucl_element(m,n)
print('V_nucl = ',V_nucl) # checked

H_core = T + V_nucl
print('H_core = ',H_core) # checked

# Calculate matrix G
C = np.zeros((2,2)) # initial guess of the solution

# Save the double electron integral results in a 4d tensor called ee_integral_tensor
ee_integral_tensor = [[[[None for i in range(2)] for j in range(2)] for k in range(2)] for l in range(2)]

# First calculate ee_integral_tensor[0][0][0][0] = (CGF_0 CGF_0 | CGF_0 CGF_0)
output = 0
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                output += coefficient[i] * coefficient[j] * coefficient[k] * coefficient[l] * V_ee_GF(GF[0][i],GF[0][j],GF[0][k],GF[0][l])
ee_integral_tensor[0][0][0][0] = ee_integral_tensor[1][1][1][1] = output
print('(CGF_0 CGF_0 | CGF_0 CGF_0) = ',ee_integral_tensor[0][0][0][0]) # checked

# Next calculate ee_integral_tensor[0][0][1][1] = (CGF_0 CGF_0 | CGF_1 CGF_1)
output = 0
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                output += coefficient[i] * coefficient[j] * coefficient[k] * coefficient[l] * V_ee_GF(GF[0][i],GF[0][j],GF[1][k],GF[1][l])
ee_integral_tensor[0][0][1][1] = ee_integral_tensor[1][1][0][0] = output
print('(CGF_0 CGF_0 | CGF_1 CGF_1) = ',ee_integral_tensor[0][0][1][1]) # checked

# Then calculate ee_integral_tensor[1][0][1][0] = (CGF_1 CGF_0 | CGF_1 CGF_0)
output = 0
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                output += coefficient[i] * coefficient[j] * coefficient[k] * coefficient[l] * V_ee_GF(GF[1][i],GF[0][j],GF[1][k],GF[0][l])
ee_integral_tensor[1][0][1][0] = ee_integral_tensor[0][1][1][0] = ee_integral_tensor[1][0][0][1] = ee_integral_tensor[0][1][0][1] = output
print('(CGF_1 CGF_0 | CGF_1 CGF_0) = ',ee_integral_tensor[1][0][1][0]) # checked


# Then calculate ee_integral_tensor[1][0][0][0] = (CGF_1 CGF_0 | CGF_0 CGF_0)
output = 0
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                output += coefficient[i] * coefficient[j] * coefficient[k] * coefficient[l] * V_ee_GF(GF[1][i],GF[0][j],GF[0][k],GF[0][l])
ee_integral_tensor[1][0][0][0] = ee_integral_tensor[0][1][0][0] = ee_integral_tensor[0][0][1][0] = ee_integral_tensor[0][0][0][1] \
    = ee_integral_tensor[0][1][1][1] = ee_integral_tensor[1][0][1][1] = ee_integral_tensor[1][1][0][1] = ee_integral_tensor[1][1][1][0] = output
print('(CGF_1 CGF_0 | CGF_0 CGF_0) = ',ee_integral_tensor[1][0][0][0]) # checked

print(np.array(ee_integral_tensor).shape)
# Define function G(C), which outputs the matrix G as a function of variable C
# We have to first arrange the columns of matrix C according to the value of eigenvalue epsilon!!!
# Assume the first column of C corresponds to the minimal eigenvalue epsilon!!!
def calculate_matrix_F(matrix_C):
    matrix_G = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            output = 0
            for k in range(2):
                for l in range(2):
                    output += matrix_C[k][0] * matrix_C[l][0] * (2 * ee_integral_tensor[i][j][k][l] - ee_integral_tensor[i][l][k][j])
            matrix_G[i][j] = output
    matrix_F = H_core + matrix_G
    return (matrix_F)

# print(calculate_matrix_F(C))

# Calculate matrix X = Us ** (-0.5), while (U^T)SU = s and s is a diagonalized matrix
# The eigenvalues of S are all positive (exercise 3.15 in Szabo's book)
print('S = ',S)
eig_val,eig_vec = LA.eig(S) # eig_val are the eigenvalues, while eig_vec are the normalized eigenvectors
# sort the eigenvectors according to the order of eigenvalues, from smaller ones to larger ones
# print('eigenvalue = ',eig_val)
# print('eigenvector = ',eig_vec)
index = np.argsort(eig_val)
# print('index = ',index)
eig_val = eig_val[index]
eig_vec = eig_vec[:,index] # here each column of eig_vec is the eigenvector of matrix S
s = np.diag(eig_val)
U = eig_vec # (np.dot(np.dot(U.T,S),U)) is diagonalized
# print('s = ',s)
# print('S eigenvalue = ',eig_val)
# print('S eigenvector = ',eig_vec)
# print('U = ',U)
# print('np.dot(U.T,U) = ',np.dot(U.T,U))
# print('np.dot(S,U) = ',np.dot(S,U))
# print('np.dot(U,s)) = ',np.dot(U,s))
# A = np.dot(S,U)
# print(np.dot(U.T,A))
# print('np.dot(np.dot(U.T,S),U) = ', np.dot(np.dot(U.T,S),U)) # checked
X = np.dot(U,np.diag([(s[0][0]) ** (-1 / 2),(s[1][1]) ** (-1 / 2)]))
# print('np.dot(np.dot(X.T,S),X) = ',np.dot(np.dot(X.T,S),X)) # checked

F = calculate_matrix_F(C)
# print('F = ',F) # should be same as H_core, since the initial guess of C is zero matrix
F_1 = np.dot(np.dot(X.T,F),X)
eig_val,eig_vec = LA.eig(F_1)
index = np.argsort(eig_val)
eig_val = eig_val[index]
eig_vec = eig_vec[:,index]
C_1 = eig_vec

C = np.dot(X,C_1)
F = calculate_matrix_F(C)
gs_energy = eig_val[0]
excited_energy = eig_val[1]

print('\n')
calculate_times = 1
print('calculate_times = ',calculate_times)
# print('calculate_times, gs_energy, excited_energy = ',calculate_times, gs_energy, excited_energy)
print('C_1 = ',C_1)
print('C = ',C)
print('F = ',F)
print('epsilon = ',eig_val)
print('\n')

flag = True # the sign of whether the iteration should go on or not

# Now do the iteration
while flag:
    calculate_times += 1
    print('calculate_times = ',calculate_times)
    old_gs_energy = gs_energy
    old_excited_energy = excited_energy
    F_1 = np.dot(np.dot(X.T,F),X)
    eig_val,eig_vec = LA.eig(F_1)
    index = np.argsort(eig_val)
    eig_val = eig_val[index]
    eig_vec = eig_vec[:,index]
    C_1 = eig_vec
    C = np.dot(X,C_1)
    F = calculate_matrix_F(C)
    epsilon = eig_val
    gs_energy = epsilon[0]
    excited_energy = epsilon[1]
    print('C_1 = ',C_1)
    print('C = ',C)
    print('F = ',F)
    print('epsilon = ',epsilon)
    print('\n')
    # print('calculate_times, gs_energy, excited_energy = ',calculate_times, gs_energy, excited_energy)

    if abs(gs_energy - old_gs_energy) < 10 ** -6:
        flag = False

print('Final result')
print('C = ',C)
print('epsilon = ',epsilon)
