# Brian H Chien
# Programming Assignment 3

#global imports











""" File parsing functions """
def create_inputdec(): #takes in text file as series of lists


    return







def main():
    #Introduction to program
    print "This is Programming Assignment 3 (2D Truss FEM Solver)"
    input_file = open('class_problem.txt', 'r')
    txt = input_file.readlines()

    for line in txt:
        ind = txt.index(line)
        line = line.strip()
        txt[ind] = line
    print txt


if __name__ == "__main__": main()