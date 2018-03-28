#Brian H Chien
#Programming Assignment 2

#global imports
import numpy as np
import matplotlib.pyplot as plt
from sympy import exp
from sympy import Symbol
from sympy import N
from sympy import Integral
from sympy import sqrt
from sympy import simplify

#problem conditions (Assume consistent units)
f = 2
A = 1
E = pow(10, 4)
P = 10

""" Preprocess Functions"""
def solve_function(global_K, force):
    #imposes BC
    before_constrain = global_K.tolist()
    for item in before_constrain:
        item[0] = 0
    force = np.delete(force, 0, 0)
    after_constrain = before_constrain

    after_constrain.pop(0)
    for item in after_constrain:
        item.pop(0)
    after_constrain = np.matrix(before_constrain)

    d = np.linalg.inv(after_constrain) * force
    d = np.insert(d, 0, [0], axis=0)

    return d

def getuppermatrix(gen_stiffness, local_dim, global_dim):
    upper_corner = []
    for i in range(0, global_dim):
        empty_row = [0] * global_dim
        if i < local_dim:
            for j in range(0, local_dim):
                empty_row[j] = gen_stiffness[i][j]
            upper_corner.append(empty_row)
        else:
            upper_corner.append(empty_row)

    return upper_corner

def getgenmatrix(gen_stiffness, local_dim, global_dim, matrix_num):
    gen_matrix = []
    gen_stiffness_address = [0, 1]
    matrix_bounds = [matrix_num, matrix_num + 1]

    for i in range(0, global_dim):
        empty_row = [0] * global_dim
        if i in matrix_bounds:
            for val in matrix_bounds:
                empty_row[val] = gen_stiffness[i - matrix_num][gen_stiffness_address[val - matrix_num]]
        gen_matrix.append(empty_row)

    return gen_matrix

def getlowermatrix(gen_stiffness, local_dim, global_dim):
    lower_corner = []
    diag_start = global_dim - 2
    for i in range(0, global_dim):
        empty_row = [0] * global_dim
        if i > global_dim - 3:
            for j in range(0, local_dim):
                empty_row[j+diag_start] = gen_stiffness[i - diag_start][j]
            lower_corner.append(empty_row)
        else:
            lower_corner.append(empty_row)
    return lower_corner

def assemble_globalstiffness(global_dimension, global_stiffness, element_stiffness):
    local_dimension = int(element_stiffness.shape[0])

    # consolidate for summation
    K_list = element_stiffness.tolist()
    global_stiffnessmatrixlist = []

    iter = 1
    matrix_num = len(range(0, global_dimension - 1))
    for i in range(0, global_dimension - 1):  # loop based on number of general matrices
        if iter == 1:
            uppercorner = getuppermatrix(K_list, local_dimension, global_dimension)
            global_stiffnessmatrixlist.append(uppercorner)
        elif iter == matrix_num:
            lowercorner = getlowermatrix(K_list, local_dimension, global_dimension)
            global_stiffnessmatrixlist.append(lowercorner)
        else:
            genmatrix = getgenmatrix(K_list, local_dimension, global_dimension, i)
            global_stiffnessmatrixlist.append(genmatrix)
        iter += 1

    jacobian_matrixsetup = []
    for item in global_stiffnessmatrixlist:
        jacobian_matrixsetup.append(np.matrix(item))

    global_K = 0
    for matrix in jacobian_matrixsetup:
        global_K += matrix

    return global_K

def preprocess(elements):
    nodes = elements + 1
    # develop general equation
    element_length = 1 / float(elements)
    a_e = A * E #serves as main coefficient

    global_stiffness = np.zeros((nodes, nodes))
    force_vector = np.zeros((nodes, 1))

    #implement internal beam force
    int_forcevector = np.zeros((nodes, 1))
    int_forcevector[0][0] = element_length
    for i in range(1, len(int_forcevector)):
        int_forcevector[i][0] = element_length * f

    #implement point load
    load_vector = np.zeros((nodes, 1))
    load_vector[nodes - 1, 0] = P

    force_vector = int_forcevector + load_vector

    #create element stiffness matrix
    K_11 = 1/float(element_length)
    K_12 = -1/float(element_length)
    K_21 = K_12
    K_22 = K_11

    element_stiffness = np.matrix([[K_11, K_12], [K_21, K_22]])

    global_K = assemble_globalstiffness(nodes, global_stiffness, element_stiffness)
    global_K = np.multiply(a_e, global_K) #AE* K

    d_params = solve_function(global_K, force_vector)
    return d_params

""" Post processing """
def postprocess(d_param, element_num):
    node_loc = []
    strain_list = []
    stress_list = []

    element_length = 1 / float(element_num)

    for i in range(0, element_num + 1):
        node_loc.append(i * element_length)

    #calculating strain
    for i in range(0, element_num):
        strain = (2/element_length) * ((d_param[i + 1] - d_param[i])/2)
        strain_list.append(strain)

    for item in strain_list:
        ind = strain_list.index(item)
        item = item.tolist()
        val = item[0][0]
        strain_list[ind] = val

    #calculating stress
    for i in range(0, element_num):
        stress = (E/element_length) * (d_param[i + 1] - d_param[i])
        stress_list.append(stress)

    stress_list.append(stress_list[-1])

    reformatted_stresslist = []

    for item in stress_list:
        item = item.tolist()
        val = item[0][0]
        reformatted_stresslist.append(val)

        stress_list = reformatted_stresslist

    x_midpt, mid_stress = calc_midptstress(element_length, node_loc)

    plot_stress(node_loc, stress_list, x_midpt, mid_stress)

def calc_midptstress(element_length, node_loc):
    half_length = element_length / 2
    updated_nodeloc = []
    midpt_stress = []

    for i in node_loc:
        if i != node_loc[-1]:
            point = i + half_length
            updated_nodeloc.append(point)

    for i in updated_nodeloc:
        stress = (2/A) * (6 - i)
        midpt_stress.append(stress)

    return updated_nodeloc, midpt_stress

def calc_qtrptstress(element_length, node_loc):
    qtr_length = element_length / 4
    updated_nodeloc = []
    qtrpt_stress = []

    for i in node_loc:
        if i != node_loc[-1]:
            point = i + qtr_length
            updated_nodeloc.append(point)

    for i in updated_nodeloc:
        stress = (2/A) * (6 - i)
        qtrpt_stress.append(stress)

    return updated_nodeloc, qtrpt_stress

def plot_stress(x1, y1, x2, y2):
    plt.figure()
    #duplicate middle nodes
    val_duplicate =
    print x1
    #plt.plot(x1, y1)
    plt.scatter(x2, y2, c = 'g')
    plt.show()
    return

def calc_error():
    return

def main():
    #Introduction to program
    elements_num = raw_input("Input the number of elements you want to implement: ")
    elements_num = int(elements_num)
    print "You have selected " + str(elements_num) + " elements."
    d_params = preprocess(elements_num)
    postprocess(d_params, elements_num)

if __name__ == "__main__": main()