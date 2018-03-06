# Brian H Chien
# ME471 Programming Assignment 1

#global imports
import numpy as np
import scipy as sci
import matplotlib as plt
import sympy as sp


#Declare global boundary value problem parameters
L = float(1.0)

""" Preprocess functions"""
def assemble_globalstiffness(empty_matrix, element_stiffness):

    global_stiffnessmatrix = empty_matrix

    return global_stiffnessmatrix

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
    print element_stiffness

    #global_K = assemble_globalstiffness(global_stiffness, element_stiffness)


    return

""" Postprocess functions """


#def postprocess():
    #return



def main():
    #Introduction to program
    elements_num = raw_input("Input the number of elements you want to implement: ")
    elements_num = int(elements_num)
    print "You have selected " + str(elements_num) + " nodes."
    preprocess(elements_num)
    #postprocess()

if __name__ == "__main__": main()