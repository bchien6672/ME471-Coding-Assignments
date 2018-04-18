# Brian H Chien
# Programming Assignment 3

""" Global Imports """
import math
import numpy as np

""" Preprocessing Functions """
def generate_elementstiffness(param_dict):
    K_dict = {}
    K_list = []

    element_data = param_dict['Element Data']
    list_elements = element_data.keys()

    sorted_elements = sorted(list_elements)

    node_database = param_dict['Node Data']

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

        basis_comp1 = node_database[first_node]['Basis Numbers']
        basis_comp2 = node_database[second_node]['Basis Numbers']

        basis = basis_comp1 + basis_comp2

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

        K_dict[element] = {}
        K_dict[element]['Element Stiffness Matrix'] = K_el
        K_dict[element]['Basis Nums'] = basis

        K_list.append(K_el)

    return K_dict

def assemble_globalstiffness(K_dict, param_dict):
    num_disp = int(param_dict['Problem Param.']['Number of Nodes']) * 2

    global_K = np.zeros((num_disp, num_disp))

    list_el = K_dict.keys()
    list_el = sorted(list_el)

    for element in list_el:
        basis_ind = K_dict[element]['Basis Nums']
        K_el = K_dict[element]['Element Stiffness Matrix']

        for basis in basis_ind:
            insert_row = basis - 1
            insertelement_row = basis_ind.index(basis)
            row_forinsertG = global_K[insert_row]
            row_forinsertE = K_el[insertelement_row]
            for basis_col in basis_ind:
                insert_col = basis_col - 1
                insert_element_col = basis_ind.index(basis_col)
                row_forinsertG[insert_col] = row_forinsertG[insert_col] + row_forinsertE[insert_element_col]

    return global_K

def generate_elementforce(param_dict): #Assume [x_dir; y_dir]
    F_list = {}
    element_data = param_dict['Element Data']
    list_elements = element_data.keys()

    sorted_elements = sorted(list_elements)

    load_list = param_dict['Force Data'].keys()

    sec_1 = [0, 1]
    sec_2 = [2, 3]

    for element in sorted_elements:
        F_el = np.zeros((4, 1))
        node_list = []

        node_list.append(element_data[element]['First Node'])
        node_list.append(element_data[element]['Second Node'])

        for node in node_list:
            if node in load_list:
                vector_ind = node_list.index(node)
                if vector_ind == 0:
                    insert_ind = sec_1
                elif vector_ind == 1:
                    insert_ind = sec_2

                F_el[insert_ind[0]] = param_dict['Force Data'][node]['X-Dir']
                F_el[insert_ind[1]] = param_dict['Force Data'][node]['Y-Dir']

        F_list[element] = F_el

    return F_list

def assemble_globalforce(param_dict):
    num_disp = int(param_dict['Problem Param.']['Number of Nodes']) * 2

    F_global = np.zeros((num_disp, 1))
    loadset = param_dict['Force Data']

    loaded_nodes = sorted(param_dict['Force Data'].keys())

    for node in loaded_nodes:
        basis_set = param_dict['Node Data'][node]['Basis Numbers']
        F_x = loadset[node]['X-Dir']
        F_y = loadset[node]['Y-Dir']

        for basis in basis_set:
            test_ind = basis_set.index(basis)
            if test_ind == 0:
                F_global[basis - 1] = F_x
            elif test_ind == 1:
                F_global[basis - 1] = F_y

    return F_global

def impose_BC(global_K, global_F, prob_params):
    constraint_set =  sorted(prob_params['Constraint Data'].keys())

    for constraint in constraint_set:
        basis_set = prob_params['Node Data'][constraint]['Basis Numbers']
        constraint_1 = prob_params['Constraint Data'][constraint]['X-disp']
        constraint_2 = prob_params['Constraint Data'][constraint]['Y-disp']

        #First constrain force vector
        for basis in basis_set:
            constrain_ind = basis_set.index(basis)
            if constrain_ind == 0:
                global_K[:, basis - 1] = global_F[basis - 1] - (global_K[:, basis - 1] * float(constraint_1))
                global_K[basis - 1] = 0
                global_K[basis - 1, basis - 1] = 1
                global_F[basis - 1] = constraint_1
            elif constrain_ind == 1:
                global_K[:, basis - 1] = global_F[basis - 1] - (global_K[:, basis - 1] * float(constraint_2))
                global_K[basis - 1] = 0
                global_K[basis - 1, basis - 1] = 1
                global_F[basis - 1] = constraint_2

    constrain_F = global_F
    constrain_K = global_K

    return constrain_K, constrain_F

def calc_disp(K_BC, F_BC):
    d_params = np.matmul(np.linalg.inv(K_BC), F_BC)

    return d_params


""" Postprocessing Functions """
def get_axialforce(forces):
    axial_forces = {}

    el_list = sorted(forces.keys())

    for element in el_list:
        F_x = forces[element][0]
        F_y = forces[element][1]
        axial_forces[element] = math.sqrt(math.pow(F_x, 2) + math.pow(F_y, 2))

    return axial_forces

def calc_elementforce(K_list, d_params):
    el_F = {}

    el_list = sorted(K_list.keys())

    for element in el_list:
        basis_set = K_list[element]['Basis Nums']
        d_set = np.zeros((len(basis_set), 1))
        K = K_list[element]['Element Stiffness Matrix']

        #assemble element disp matrix
        for basis in basis_set:
            ind = basis_set.index(basis)
            d_set[ind] = d_params[basis - 1]

        F_set = np.matmul(K, d_set)
        F_set = F_set / 1000 #convert to kip
        el_F[element] = F_set

    return el_F


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
    param_dict['Problem Param.']['Constrained Nodes'] = nodal_info[3]

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

    if nodal_info[2] != '0':
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
    else: #For Example 2
        constraint_data = parsed_txt[3]
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
    input_file = open('class_problem.txt', 'r') #change name of file to Example1.txt, Example2.txt, any new input files need two lines of blank space after the last line
    txt = input_file.readlines()

    for line in txt:
        ind = txt.index(line)
        line = line.strip()
        txt[ind] = line

    input_param = create_inputdec(txt)

    #preprocess
    #generate stiffness matrices
    K_dict = generate_elementstiffness(input_param)
    global_K = assemble_globalstiffness(K_dict, input_param)

    #generate force vectors
    #F_el = generate_elementforce(input_param)
    F_global = assemble_globalforce(input_param)

    #Solve for d
    constrain_K, constrain_F = impose_BC(global_K, F_global, input_param)
    d_params = calc_disp(constrain_K, constrain_F) #nodal displacement matrix

    #postprocess
    element_F = calc_elementforce(K_dict, d_params) #element reaction forces
    axial_Fmag = get_axialforce(element_F) #in kip

    #print out the nodal displacement (in), element reaction forces, and the magnitude of the force on each element (both in kip)
    print 'Nodal displacements (in) is: ' + str(d_params)
    print 'Element forces (kip) are: ' + str(element_F)
    print 'Axial forces (kip) are: ' + str(axial_Fmag)

if __name__ == "__main__": main()