# Brian H Chien
# ME471 Programming Assignment 1

#global imports
import numpy as np
import matplotlib.pyplot as plt
from sympy import exp
from sympy import Symbol
from sympy import N
from sympy import integrate
from sympy import sqrt

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

    #General eqs. for reference (commented for code optimization)
    #N_1 = (1/float(2)) * (1 - E)
    #N_2 = (1/float(2)) * (1 + E)
    #N_1e = 1/float(2)
    #N_2e = 1/float(2)

    element_length = 1/float(elements)
    node_loc = []
    for i in range(0, elements + 1):
        node_loc.append(i * element_length)

    E_fcn = []
    for i in range(0, len(node_loc) - 1):
        shape_fcn = ((2*x) - node_loc[i] - node_loc[i + 1])/element_length
        E_fcn.append(shape_fcn)

    u_h = {}

    for i in range(0, len(E_fcn)):
        eq_desig = 'N' + str(i + 1) + str(i+2)
        u_h[eq_desig] = {}
        u_h[eq_desig]['Constraints'] = str(node_loc[i]) + ' < ' + str(x) + ' < ' + str(node_loc[i + 1])
        u_h[eq_desig]['Equation'] = (d_params[i] * (1/float(2)) * (1 - E_fcn[i])) + (d_params[i + 1] * (1/float(2)) * (1 + E_fcn[i]))
        u_h[eq_desig]['Bounds'] = [node_loc[i], node_loc[i + 1]]
        u_h[eq_desig]['Upper Bound'] = node_loc[i + 1]
        u_h[eq_desig]['Lower Bound'] = node_loc[i]

    u_exact = ((exp(1))/(exp(2) + 1)) * (exp(x) - exp(-x))

    eval_node = create_nodeevaltable(u_exact, u_h, x, node_loc)

    #copied and pasted values for nodal value table in Excel
    #for item in eval_node.keys():
        #print item
        #for val in eval_node[item]:
            #print val

    error_L2 = calc_error(u_exact, u_h, x, node_loc, eval_node)

    if elements == 5: #Can remove this line to plot all numbers of elements
        plotcomparisonfunction(u_exact, u_h, x, node_loc)

    return

def create_nodeevaltable(u_exact, u_h, x, node_loc):
    eval_node = {}

    x_input = node_loc
    eval_node['x'] = x_input

    y_exact = []
    for i in x_input:
        y_exact.append(N(u_exact.subs(x, i)))
    eval_node['u_exact'] = y_exact

    y_uh = [0] * len(node_loc)
    for key in u_h.keys():
        for i in range(1, len(node_loc)):
            if node_loc[i] == u_h[key]['Upper Bound']:
                eq = u_h[key]['Equation']
                y_uh[i] = eq.subs(x, node_loc[i])
                break

    eval_node['u_h'] = y_uh
    return eval_node

def plotcomparisonfunction(u_exact, u_h, x, node_loc):
    #Plot exact and u_h
    #get exact soln
    x_input = np.linspace(0.0, L, num = 1000)
    x_input = x_input.tolist()
    y_exact = []

    for i in x_input:
        y_exact.append(u_exact.subs(x, i))

    y_uh = []
    x_uh = []
    for i, key in zip(range(0, len(node_loc) - 1), u_h.keys()):
        x_range = np.linspace(node_loc[i], node_loc[i + 1], 500)
        x_range = x_range.tolist()
        if node_loc[i] == u_h[key]['Lower Bound']:
            eq = u_h[key]['Equation']
            for j in x_range:
                x_uh.append(j)
                y_uh.append(eq.subs(x, j))

    plt.figure()
    plt.plot(x_input, y_exact)
    plt.plot(x_uh, y_uh, '--', linewidth = 2)
    plt.xlabel('x')
    plt.ylabel('Solution')
    plt.legend(['Exact', 'u_h'])
    plt.show()
    return

def calc_error(u_exact, u_h, x, node_loc, node_valdict):
    integral_sum = 0

    for item in u_h.keys():
        bounds = u_h[item]['Bounds']
        u_heq = u_h[item]['Equation']
        val = N(integrate(pow((u_exact - u_heq),2), (x, bounds[0], bounds[1])))
        integral_sum += val

    error_L2 = N(sqrt(integral_sum))

    return error_L2

def main():
    #Introduction to program
    elements_num = raw_input("Input the number of elements you want to implement:   ")
    elements_num = int(elements_num)
    print "You have selected " + str(elements_num) + " elements."
    d_params = preprocess(elements_num)
    postprocess(d_params, elements_num)

if __name__ == "__main__": main()