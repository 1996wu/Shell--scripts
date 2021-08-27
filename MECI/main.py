import numpy as np
import re
import os
import sys

from scipy.optimize import minimize

from GAUSSIAN import *

ang = 0.529177210903
eV = 27.211386245988
sigma = 3.5
alpha = 0.02
nstates = 1
gtol = 1.0E-05
element = []
xyz = []
conver_1 = 0.0
conver_2 = 0.0
conver_3 = 0.0
f_before = 0.0
conical_delta_ij = 0.0


def gap_function(geom: list):
    """
    This is script which can optimize Conical Intersections without Derivative Coupling Vectors
    with GAUSSIAN16 interface.
    ref: "J. Phys. Chem. B 2008, 112, 405-413"
    """
    global iteration
    global conver_1, conver_2, conver_3, f_before, conical_delta_ij

    # geom (3N,1) -> xyz(N,3) (unit Bohr)
    coord = [geom[i:i + 3] for i in range(0, len(geom), 3)]
    global xyz
    xyz = coord

    # calculate the i state
    replace_coordinate(coord)
    software_running()
    grad_i = get_grad_matrix().reshape(natom * 3)
    energy_i = get_energy()[1]
    print_traj_cicoe(iteration, flag=False)

    # calculate the j state j=i-1
    renew_calc_states(nstates - 1)
    software_running()
    read_wavefunction()
    grad_j = get_grad_matrix().reshape(natom * 3)
    energy_j = get_energy()[1]
    delete_wavefunction()
    renew_calc_states(nstates)

    delta_ij = energy_i - energy_j
    if delta_ij < 0:
        print("The geometry may be error because the excited energy of %s state is negative" % nstates)
        sys.exit()
    print_xyz(coord)

    # calculate the  value and grad of the function(F_ij, adn grad_F_ij)
    E_ij = (energy_i + energy_j) * 0.5
    G_ij = delta_ij ** 2 / (delta_ij + alpha)
    F_ij = E_ij + sigma * G_ij
    grad_E_ij = 0.5 * (grad_i + grad_j)
    grad_G_ij = (delta_ij ** 2 + 2 * alpha * delta_ij) / ((delta_ij + alpha) ** 2) * (grad_i - grad_j)
    grad_F_ij = grad_E_ij + sigma * grad_G_ij

    conical_delta_ij = 0.0
    norm_g = grad_G_ij / np.linalg.norm(grad_G_ij)
    # This equation 10, 11, 12 three criteria  convergence:
    # 1: F_ij_next - F_ij  <= 1.0E-06 Hartree
    # 2: 1 / sigma * grad_F_ij . u <= 1.0E-03 Hartree/Bohr
    # 3: ||grad_ij - (grad_F_ij . u). u || <= 1.0E-03 Hartree/Bohr
    conver_1 = F_ij - f_before
    conver_2 = 1 / sigma * np.dot(grad_F_ij, norm_g)
    conver_3 = np.linalg.norm(grad_F_ij - np.dot(grad_F_ij, norm_g) * norm_g)

    # The energy gap between of the i state and the j state is 0.001 Hartree/ 0.0273eV
    if delta_ij < 0.001:
        z_ij = 0.5 * np.dot((grad_i - grad_j), (grad_i + grad_j)) / (
                0.25 * np.dot((grad_i + grad_j), (grad_i + grad_j)))
        conical_delta_ij = alpha * (np.sqrt((sigma * z_ij) / (1 + sigma * z_ij)) - 1)

    print("The total function value {:26.6f}".format(F_ij), flush=True)
    print("The total function grad norm {:18.6}(eV)".format(np.linalg.norm(grad_F_ij) * eV), flush=True)
    print("The energy of i and j states and delta_ij  is {:18.6f} {:18.6f} {:18.5f}(eV)"
          .format(energy_i, energy_j, delta_ij * eV), flush=True)
    print("The near MECI geometry(F_ij) is{:18.6f}".format(conical_delta_ij), flush=True)
    print("convergence_1,convergence_2 and convergence_3{:18.6f}{:18.6f}{:18.6f}"
          .format(conver_1, conver_2, conver_3), flush=True)
    print("alpha = {:18.6f}   sigma = {:18.6}".format(alpha, sigma), flush=True)
    print("The the %s iteration has achieved at %s\n" % (iteration, current_time()), flush=True)

    # if iteration >= 50:
    #     print("The script iteration (%s)time too much and has been ended at %s \n" % (iteration, current_time()))
    #     # sys.exit()
    
    # gradient projection (GP): ref : J. Chem. Theory Comput. 2010, 6, 5, 1538–1545
    # g^GP = 2[E^X(Q) - E^Y(Q)] V^{DGV} + P * 0.5 [grad_E^x(Q) + grad_E^y(Q)] 
    # P = 1 -V^{DGV} * (V^{DGV})^T - V^{DCV} * (V^{DCV})^T
    # branching plane (BP)
    # v^{DGV} is an unit vector of the difference gradient vector (DGV)
    # a unit vector perpendicular to v^{DGV} insteal of V^{DCV}
    # v_DGV = (grad_i - grad_j) / np.linalg.norm(grad_i - grad_j)
    # V_DCV = 
    # P = 1 - np.dot(V_DGV, v_DGV.T) - np.dot(V_DCV, V_DCV.T)
    # g_GP = 2(energy_i - energy_j) * v_DGV + np.dot(P, (grad_i + grad_j)) * 0.5

    iteration += 1
    if iteration >= 1:
        f_before = F_ij

    return F_ij, grad_F_ij

def current_time():
    return time.asctime(time.localtime(time.time()))

def print_xyz(coordinate: list):
    with open('simulation.xyz', 'a+') as f:
        f.write(str(natom) + '\n')
        f.write('simulation time: t =' + format(iteration, '>10d') + '\n')
        for i in range(natom):
            ele = format(element[i], '<5s')
            coord = ''.join(format(x * ang, '>18.10f') for x in coordinate[i])
            f.write(ele + coord + '\n')


def read_initial_geom():
    # element x , y ,z (unit/ang)
    with open('initial_condition', 'r') as f:
        geom = []
        for value in f:
            if not value:
                break
            elem = value.split()[0].capitalize()
            data = [float(i) / ang for i in value.split()[1:4]]
            geom.append(data)
            element.append(elem)
    return np.array(geom), len(element)


def main():
    global iteration, natom, sigma, xyz
    global conver_1, conver_2, conver_3, f_before, conical_delta_ij
    max_sigma = 100
    xyz, natom = read_initial_geom()
    print("This is script which can optimize Conical Intersections without\nDerivative Coupling Vectors at TDDFT "
          "level with GAUSSIAN16 interface")
    while True:
        iteration = 0
        print("------The BFGS algorithm begin---------\n")
        res = minimize(gap_function, xyz, method='BFGS', jac=True, options={'disp': True, 'gtol': gtol, 'maxiter': 50})
        print(res)
        if abs(conver_3) <= 5.0E-03 and abs(conver_2) <= 5.0E-03 and abs(conver_1) <= 1.0E-06 and abs(
                conical_delta_ij) <= 0.001:
            print("The there criteria convergence and delta_ij  has been met{:18.6f}{:18.6f}{:18.6f}{:18.6f}"
                  .format(conver_1, conver_2, conver_3, conical_delta_ij), flush=True)
            break
        else:
            if sigma <= max_sigma:
                sigma = sigma * 4
                print("Re-adjust sigma{:18.6f} and continue to iterate".format(sigma), flush=True)
            else:
                sigma = max_sigma


if __name__ == "__main__":
    main()
