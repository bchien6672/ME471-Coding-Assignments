# Brian H Chien
# Programming Assignment 3

""" Global Imports """
import math
import numpy as np

""" Preprocessing Functions """
def generate_elementstiffness(param_dict):
    K_list = []

    element_data = param_dict['Element Data']
    list_elements = element_data.keys()

    sorted_elements = sorted(list_elements)

    node_database = param_dict['Node Data']
    print node_database
    #assuming all coordinate basis are in the same orientation
    for element in sorted_elements:
        element_subset = element_data[element]
        first_node = element_subset['First Node']
        second_node = element_subset['Second Node']
        element_area = float(element_subset['Cross-Section Area'])
        y_mod = float(element_subset['Y-Modulus'])

        node_x1 = float(node_database[first_node]['X'])
        node_x2 = float(node_database[second_node]['X'])

        node_y1 = float(node_database[first_node]['Y'])
        node_y2 = float(node_database[second_node]['Y'])

        node_x = node_x2 - node_x1
        node_y = node_y2 - node_y1

        element_length = math.sqrt(math.pow(node_x, 2) + math.pow(node_y, 2))

        l_x = node_x/element_length
        l_y = node_y/element_length

        K_el = np.array([[pow(l_x, 2), (l_x * l_y), -(pow(l_x, 2)), -(l_x * l_y)],
                         [(l_x * l_y), (pow(l_y, 2)), -(l_x * l_y), -(pow(l_y, 2))],
                         [-(pow(l_x, 2)), -(l_x * l_y),  pow(l_x, 2), (l_x * l_y)],
                         [-(l_x * l_y), -(pow(l_y, 2)), (l_x * l_y), (pow(l_y, 2))]])

        K_el = ((element_area * y_mod)/element_length) * K_el

        K_list.append(K_el)

    return K_list

def assemble_globalstiffness(K_list):
    return

def generate_elementforce(param_dict):

    return
def assemble_globalforce():

    return

def impose_BC():

    return

def calc_disp():

    return


""" Postprocessing Functions """
def get_elementaxialforce():
    return

def calc_elementforce():

    return






""" File parsing functions """
def text_parser(txt):
    parsed_txt = []
    section = []

    for line in txt:
        section.append(line)
        if line == '':
            parsed_txt.append(section)
            section = []

    for line in parsed_txt:
        line.remove('')

    return parsed_txt


def create_inputdec(txt): #takes in text file as series of lists
    param_dict = {}

    parsed_txt = text_parser(txt)

    param_dict['Problem Param.'] = {}
    param_dict['Node Data'] = {}
    param_dict['Element Data'] = {}
    param_dict['Force Data'] = {}
    param_dict['Constraint Data'] = {}

    nodal_info = parsed_txt[0][0].split()
    param_dict['Problem Param.']['Number of Nodes'] = nodal_info[0]
    param_dict['Problem Param.']['Number of Elements'] = nodal_info[1]
    param_dict['Problem Param.']['Loaded Nodes'] = nodal_info[2]
    param_dict['Problem Param.']['Constained Nodes'] = nodal_info[3]

    element_coords = parsed_txt[1]

    basis_list = [1, 2]

    coord_basisnum = 1
    for pt in element_coords:
        pt_list = pt.split()
        node_val = pt[0]
        param_dict['Node Data'][node_val] = {}
        param_dict['Node Data'][node_val]['X'] = pt_list[1]
        param_dict['Node Data'][node_val]['Y'] = pt_list[2]
        param_dict['Node Data'][node_val]['Basis Numbers'] = []
        for num in basis_list:
            if num == 1:
                basis_set = (2 * int(node_val)) - 1
                param_dict['Node Data'][node_val]['Basis Numbers'].append(basis_set)
            elif num == 2:
                basis_set = int(node_val) * 2
                param_dict['Node Data'][node_val]['Basis Numbers'].append(basis_set)

    element_data = parsed_txt[2]
    for element in element_data:
        element_component = element.split()
        element_num = 'Element ' + str(element_component[0])
        param_dict['Element Data'][element_num] = {}
        param_dict['Element Data'][element_num]['First Node'] = element_component[1]
        param_dict['Element Data'][element_num]['Second Node'] = element_component[2]
        param_dict['Element Data'][element_num]['Cross-Section Area'] = element_component[3]
        param_dict['Element Data'][element_num]['Y-Modulus'] = element_component[4]

    force_data = parsed_txt[3]
    for force in force_data:
        force_point = force.split()
        force_at_node = str(force_point[0])
        param_dict['Force Data'][force_at_node] = {}
        param_dict['Force Data'][force_at_node]['X-Dir'] = force_point[1] #lb
        param_dict['Force Data'][force_at_node]['Y-Dir'] = force_point[2]

    constraint_data = parsed_txt[4]
    for constraint in constraint_data:
        constraint_node = constraint.split()
        node_num = str(constraint_node[0])
        param_dict['Constraint Data'][node_num] = {}

        x_constrainbool = constraint_node[1]
        y_constrainbool = constraint_node[2]

        if x_constrainbool != 1:
            param_dict['Constraint Data'][node_num]['X-disp'] = constraint_node[3]
        else:
            param_dict['Constraint Data'][node_num]['X-disp'] = 0

        if y_constrainbool != 1:
            param_dict['Constraint Data'][node_num]['Y-disp'] = constraint_node[4]
        else:
            param_dict['Constraint Data'][node_num]['Y-disp'] = 0

    return param_dict

def main():
    #Introduction to program
    print "This is Programming Assignment 3 (2D Truss FEM Solver)"
    input_file = open('class_problem.txt', 'r')
    txt = input_file.readlines()

    for line in txt:
        ind = txt.index(line)
        line = line.strip()
        txt[ind] = line

    input_param = create_inputdec(txt)
    print input_param
    #preprocess
    K_element_list = generate_elementstiffness(input_param)

    #postprocess

if __name__ == "__main__": main()