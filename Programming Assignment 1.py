# Brian H Chien
# ME471 Programming Assignment 1

#global imports
import numpy as np
import matplotlib as plt
from sympy import Symbol


#Declare global boundary value problem parameters
L = float(1.0)

""" Preprocess functions"""

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
    print d

    return

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

    #consolidate for summation
    K_list = element_stiffness.tolist()
    global_stiffnessmatrixlist = []

    iter = 1
    matrix_num = len(range(0, global_dimension - 1))
    for i in range(0, global_dimension - 1): # loop based on number of general matrices
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
    #declare global stiffness and force vector matrices
    nodes = elements + 1
    global_stiffness = np.zeros((nodes, nodes))
    force_vector = np.zeros((nodes, 1))
    force_vector[nodes - 1,0] = 1

    #develop general equation
    element_length = 1/float(elements)
    K_11 = (1/float(element_length)) + (float(element_length)/3)
    K_22 = K_11
    K_21 = -(1/float(element_length)) + (float(element_length)/6)
    K_12 = K_21
    element_stiffness = np.matrix([[K_11, K_12],[K_21, K_22]])

    global_K = assemble_globalstiffness(nodes, global_stiffness, element_stiffness)

    #impose boundary conditions
    d_params = solve_function(global_K, force_vector)
    return d_params

""" Postprocess functions """

def postprocess(d_params, elements):
    #declare variables
    E = Symbol('E')
    x = Symbol('x')

    N_1 = (1/float(2)) * (1 - E)
    N_2 = (1/float(2)) * (1 + E)
    N_1e = 1/float(2)
    N_2e = 1/float(2)

    element_length = 1/float(elements)
    node_loc = []
    for i in range(0, elements + 1):
        node_loc.append(i * element_length)

    

    return

#def plotcomparisonfunction():
    #return

#def calc_error():
    #return



def main():
    #Introduction to program
    elements_num = raw_input("Input the number of elements you want to implement:   ")
    elements_num = int(elements_num)
    full_soln_flag = False
    if elements_num == 5:
        full_soln_flag = True
    print "You have selected " + str(elements_num) + " elements."
    d_params = preprocess(elements_num)

    #if full_soln_flag == True:
    postprocess(d_params, elements_num)
    #else:
        #pass
if __name__ == "__main__": main()