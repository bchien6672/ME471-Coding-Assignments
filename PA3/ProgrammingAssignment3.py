# Brian H Chien
# Programming Assignment 3

""" Preprocessing Functions """






""" Postprocessing Functions """





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
    for pt in element_coords:
        pt_list = pt.split()
        node_val = 'Node ' + str(pt[0])
        param_dict['Node Data'][node_val] = {}
        param_dict['Node Data'][node_val]['X'] = pt_list[1]
        param_dict['Node Data'][node_val]['Y'] = pt_list[2]

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
        force_at_node = 'Node ' + str(force_point[0])
        param_dict['Force Data'][force_at_node] = {}
        param_dict['Force Data'][force_at_node]['X-Dir'] = force_point[1] #lb
        param_dict['Force Data'][force_at_node]['Y-Dir'] = force_point[2]

    constraint_data = parsed_txt[4]
    for constraint in constraint_data:
        constraint_node = constraint.split()
        node_num = 'Node ' + str(constraint_node[0])
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

    #preprocess

    #postprocess

if __name__ == "__main__": main()