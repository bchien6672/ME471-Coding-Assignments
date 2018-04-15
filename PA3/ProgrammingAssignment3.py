# Brian H Chien
# Programming Assignment 3

""" Preprocessing Functions """






""" Postprocessing Functions """





""" File parsing functions """
def text_parser(txt):
    parsed_txt = []
    section = []

    for line in txt:
        print line
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

    #Assume there are only four coordinates
    element_coords = parsed_txt[1]
    for pt in element_coords:
        pt_list = pt.split()
        node_val = 'Node ' + str(pt[0])
        param_dict['Node Data'][node_val] = {}
        param_dict['Node Data'][node_val]['X'] = pt_list[1]
        param_dict['Node Data'][node_val]['Y'] = pt_list[2]

    

    print param_dict

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