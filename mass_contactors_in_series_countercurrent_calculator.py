# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:43:13 2026

@author: VEShultz75
"""

running = True


"""
First we have the functions for constructing the matrices. 
Also for future me, I have to use lambdaa as a variable because lambda is an in-built value of some kind 

***I will need to go back and make sure to make a separate case for 
"""

def construct_stage_one(phase_II_initial_conc, lambdaa, eq_const, num_stages):
    # Right now this only accounts for two or more stages. I will add another case for 
    # only one stage later. 
    eq_one_row = [-1, eq_const, 0 , 0];
    eq_two_row = [0, -1, 1/(lambdaa+eq_const), 0]
    for i in range(4, 2*num_stages):
        eq_one_row.append(0)
        eq_two_row.append(0)
    eq_one_row.append(0)
    eq_two_row.append(-(lambdaa*phase_II_initial_conc)/(lambdaa+eq_const))
    return [eq_one_row, eq_two_row];


def construct_stage_n(lambdaa, eq_const, num_stages, stage_number):
    # I can consolidate this by populating the rows with 0s then replacing 0s with the necessary values. 
    # But that will be less efficient. 
    eq_2n_1_row = [] 
    eq_2n_row = []
    for i in range(0, 2*stage_number-3):
        eq_2n_1_row.append(0)
        eq_2n_row.append(0)
    eq_2n_1_row.append(0)
    eq_2n_1_row.append(-1)
    eq_2n_1_row.append(eq_const)
    eq_2n_1_row.append(0)
    eq_2n_row.append(lambdaa/(lambdaa+eq_const))
    eq_2n_row.append(0)
    eq_2n_row.append(-1)
    eq_2n_row.append(1/(lambdaa+eq_const))
    for i in range(2*stage_number+1, 2*num_stages+1):
        eq_2n_1_row.append(0)
        eq_2n_row.append(0)
    return [eq_2n_1_row, eq_2n_row]


def construct_final_stage(phase_I_initial_conc, lambdaa, eq_const, num_stages):
    eq_2N_1_row = []
    eq_2N_row = []
    for i in range(0, 2*num_stages-3):
        eq_2N_1_row.append(0)
        eq_2N_row.append(0)
    eq_2N_1_row.append(0)
    eq_2N_1_row.append(-1)
    eq_2N_1_row.append(eq_const)
    eq_2N_1_row.append(0)
    eq_2N_row.append(lambdaa/(lambdaa+eq_const))
    eq_2N_row.append(0)
    eq_2N_row.append(-1)
    eq_2N_row.append(-phase_I_initial_conc/(lambdaa+eq_const))
    return [eq_2N_1_row, eq_2N_row];


def construct_matrix(phase_I_initial_conc, 
                     phase_II_initial_conc, lambdaa, eq_const, 
                     num_stages):
    matrix = [];
    stage_one = construct_stage_one(phase_II_initial_conc, lambdaa, eq_const, num_stages)
    matrix.append(stage_one[0]) 
    matrix.append(stage_one[1])
    for i in range(2, num_stages):
        stage_n = construct_stage_n(lambdaa, eq_const, num_stages, i)
        matrix.append(stage_n[0])
        matrix.append(stage_n[1])
    final_stage = construct_final_stage(phase_I_initial_conc, lambdaa, eq_const, num_stages)
    matrix.append(final_stage[0])
    matrix.append(final_stage[1])
    return matrix;



"""
Now we need the functions for reducing the matrix to reduced row echelon form. 
We start with basic row operations, and then make an algorithm to solve the matrix. 
""" 

def subtract_row(matrix, row_a_index, row_b_index, multiplier):
    new_row_a = []
    for i in range(0, len(matrix[0])):
        new_row_a.append(matrix[row_a_index][i] - multiplier*matrix[row_b_index][i])
    matrix[row_a_index] = new_row_a
    return matrix 

def multiply_row(matrix, row_a_index, multiplier):
    new_row_a = []
    for i in range(0, len(matrix[0])):
        new_row_a.append(multiplier*matrix[row_a_index][i]);
    matrix[row_a_index] = new_row_a
    return matrix 

def reduce_to_row_echelon_form(matrix):
    for i in range(0, len(matrix)-1):
        for j in range(i+1,len(matrix)):
            #the col index of the leading entry for i row is the same as i 
            if matrix[j][i] != 0.0:
                matrix = subtract_row(matrix, j, i, matrix[j][i]/matrix[i][i])
    return matrix

def reduce_to_rref_form(matrix):
    matrix = reduce_to_row_echelon_form(matrix)
    i = len(matrix) - 1
    for u in range(0, len(matrix)):
        matrix = multiply_row(matrix, i, 1/matrix[i][i])
        j = i - 1
        for v in range(u+1, len(matrix)):
            if matrix[j][i] != 0.0:
                matrix = subtract_row(matrix, j, i, matrix[j][i])
            j = j - 1
        i = i - 1
    return matrix 



"""
Now that the functions for creating and reducing the matrix have been specified, 
I will specify the functions for extracting the relevant information and 
displaying the results to the user. 
""" 

def extract_concentrations(rref_matrix):
    conc_list = []
    for i in range(0, len(rref_matrix)):
        conc_list.append(rref_matrix[i][len(rref_matrix)])
    return conc_list

def display_system_information(phase_I_initial_conc, phase_II_initial_conc, 
                               volumetric_flow_rate_I, volumetric_flow_rate_II, 
                               eq_const, num_stages):
    print("Multi-Stage Countercurrent System Information:")
    print("phase I initial concentration:", phase_I_initial_conc)
    print("phase II initial concentration:", phase_II_initial_conc)
    print("phase I volumetric flow rate:", volumetric_flow_rate_I)
    print("phase II volumetric flow rate:", volumetric_flow_rate_II)
    print("equilibrium constant:", eq_const)
    print("number of stages:", num_stages)
    
def display_concentrations(concentrations):
    stage_number = 1
    print("Concentrations:")
    for i in range(0, len(concentrations)):
        if i%2 == 0:
            print("phase I  stage ", stage_number, " concentration: ", round(concentrations[i], 5))
        else: 
            print("phase II stage ", stage_number, " concentration: ", round(concentrations[i], 5))
            stage_number = stage_number + 1


"""
Finally, I need the functions for interacting with the user. 

A majority of these functions are for the user to input the necessary variables 
to build the matrix and solve the system of equations. 

Note: At any time, the user can exit the program by entering “e” as an input, 
and that will be reflected in the various “get input” functions. 
"""

def string_is_float(string):
    # a valid float string must: only have one decimal, and otherwise be only 
    # comprised of numbers 
    is_string = True
    has_at_least_one_decimal = False 
    for letter in string:
        if letter.isnumeric():
            pass 
        elif letter == ".":
            if has_at_least_one_decimal != True:
                has_at_least_one_decimal = True
            else: 
                is_string = False
        else:
            is_string = False
    return is_string

def test_case(case_num):
    test_message = ""
    hand_calc_conc = []
    error_message = ""
    num_stages = 0;
    if case_num == 1:
        test_message = "Test Case 1: 2 Stages\n"
        hand_calc_conc = [129/65.0, 43/65.0, 42/65.0, 14/65.0]
        error_message = "Test Case 1 Failed. Incorrect Concentrations."
        num_stages = 2
    elif case_num == 2:
        test_message = "Test Case 2: 3 Stages\n"
        hand_calc_conc = [417/200.0, 139/200.0, 39/50.0, 13/50.0, 69/200.0, 23/200.0]
        error_message = "Test Case 2 Failed. Incorrect Concentrations."
        num_stages = 3 
    elif case_num == 3:
        test_message = "Test Case 3: 4 Stages\n"
        hand_calc_conc = [1281/605, 427/605, 498/605, 166/605, 237/605, 79/605, 30/121, 10/121]
        error_message = "Test Case 3 Failed. Incorrect Concentrations."
        num_stages = 4 
    print(test_message)
    initial_conc_I = 0.2
    initial_conc_II = 2.0
    vol_flow_rate_I = 100
    vol_flow_rate_II = 100
    eq_const = 3.0
    display_system_information(initial_conc_I, initial_conc_II, 
                               vol_flow_rate_I, vol_flow_rate_II, eq_const, 
                               num_stages)
    matrix = construct_matrix(initial_conc_I, initial_conc_II, 
                               vol_flow_rate_II/vol_flow_rate_I, eq_const, 
                               num_stages)
    rref = reduce_to_rref_form(matrix)
    program_conc = extract_concentrations(rref)
    print("\nConcentrations calculated by the program:")
    print(program_conc)
    print("\nConcentrations calculated by hand:")
    print(hand_calc_conc)
    rounded_program_conc = []
    rounded_hand_calc_conc = []
    for i in range(0, len(program_conc)):
        rounded_program_conc.append(round(program_conc[i], 10))
        rounded_hand_calc_conc.append(round(hand_calc_conc[i], 10))
    assert rounded_program_conc == rounded_hand_calc_conc, error_message 
    string = ""
    for i in range(0, 80):
        string = string + "-"
    print(string)
    
    
    

def startup_message():
    test_case(1)
    test_case(2)
    test_case(3)
    print("Enter 'e' at any time to exit the program.")
    print("All values must be greater than 0.")
    #This is an example of how the assertion works. 
    assert sum([1, 2, 3]) == 6, "Should be 6"
    

def get_num_stages():
    question_string = "How many stages?"
    invalid = True 
    while invalid:
        print(question_string)
        num_stages = input()
        if num_stages.isnumeric() or num_stages == "e":
            invalid = False
            if num_stages.isnumeric():
                num_stages = int(num_stages)
        else:
            print("Invalid Input")
    return num_stages
    
# You know what, we can consolidate these because the logic is the same. 
# The only thing that needs to change is the question string. 

float_questions = ["What is the initial concentration of 'A' in phase I?",
                   "What is the initial concentration of 'A' in phase II?",
                   "What is the volumetric flow rate of phase I?",
                   "What is the volumetric flow rate of phase II?",
                   "What is the equilibrium constant for this system?"]

def get_float(question_index):
    # Value is the index corresponding to the entry of float_questions 
    # that we are retrieving. ex. 0 = fetching equilibrium constant. 
    question_string = float_questions[question_index]
    invalid = True 
    while invalid:
        print(question_string)
        user_input = input()
        if string_is_float(user_input) or user_input == "e":
            invalid = False
            if string_is_float(user_input):
                user_input = float(user_input)
        else:
            print("Invalid Input")
    return user_input



def display_matrix(matrix):
    string = ""
    for i in range(0, 80):
        string = string + "-"
    print(string)
    print("[")
    for entry in matrix:
        print(entry)
    print("]")
    print(string)



"""
Finally, the logic of the program is contained within a while loop. 
"""
        
    
while running:
    startup_message()
    num_stages = get_num_stages()
    if num_stages == "e":
        break
    # the elements of this list are as follows: [phase I initial concentration,
    # phase II initial concentration, phase I volumetric flow rate, 
    # phase II volumetric flow rate, equilibrium constant]
    system_values = []
    exit_command = False
    for i in range(0, len(float_questions)):
        value = get_float(i)
        if value == "e":
            exit_command = True
            break 
        else:
            system_values.append(value)
    if exit_command:
        break
    initial_conc_I = system_values[0]
    initial_conc_II = system_values[1]
    vol_flow_rate_I = system_values[2]
    vol_flow_rate_II = system_values[3]
    eq_const = system_values[4]
    display_system_information(initial_conc_I, initial_conc_II, 
                               vol_flow_rate_I, vol_flow_rate_II, eq_const, 
                               num_stages)
    lambdaa = vol_flow_rate_II/vol_flow_rate_I
    matrix = construct_matrix(initial_conc_I, initial_conc_II, lambdaa, 
                              eq_const, num_stages)
    display_matrix(matrix)
    matrix = reduce_to_rref_form(matrix)
    conc_list = extract_concentrations(matrix)
    display_concentrations(conc_list)
    string = ""
    for i in range(0, 80):
        string = string + "-"
    print(string)
print("Exiting Program")
        

















