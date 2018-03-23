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
E = 10e4
P = 10

""" Preprocess Functions"""
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
    global_stiffness = np.zeros((nodes, nodes))
    force_vector = np.zeros((nodes, 1))

    #implement point load
    load_vector = np.zeros((nodes, 1))
    load_vector[nodes - 1, 0] = P

    # develop general equation
    element_length = 1 / float(elements)
    a_e = A * E #serves as main coefficient

    #create element stiffness matrix
    K_11 = 1/float(element_length)
    K_12 = -1/float(element_length)
    K_21 = K_12
    K_22 = K_11

    element_stiffness = np.matrix([[K_11, K_12], [K_21, K_22]])

    global_K = assemble_globalstiffness(nodes, global_stiffness, element_stiffness)
    global_K = np.multiply(a_e, global_K) #AE* K

    d_params = None
    return d_params



def main():
    #Introduction to program
    elements_num = raw_input("Input the number of elements you want to implement: ")
    elements_num = int(elements_num)
    print "You have selected " + str(elements_num) + " elements."
    d_params = preprocess(elements_num)
    #postprocess(d_params, elements_num)

if __name__ == "__main__": main()