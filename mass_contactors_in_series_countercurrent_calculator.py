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

def construct_stage_one(p_II_feed_conc: float, lambdaa:float, eq_const:float, 
                        num_stages:int) -> list[list[float]]:
    """
    The first two equations corresponding to stage one will always have the 
    same configuration of variables. This function creates the first two rows 
    of the matrix that correspond to stage 1. 
    """
    # Right now this only accounts for two or more stages. I will add another case for 
    # only one stage later. 
    eq_one_row = [-1, eq_const, 0 , 0];
    eq_two_row = [0, -1, 1/(lambdaa+eq_const), 0]
    for i in range(4, 2*num_stages):
        eq_one_row.append(0)
        eq_two_row.append(0)
    eq_one_row.append(0)
    eq_two_row.append(-(lambdaa*p_II_feed_conc)/(lambdaa+eq_const))
    return [eq_one_row, eq_two_row];


def construct_stage_n(lambdaa: float, eq_const: float, num_stages: int, 
                      stage_number: int) -> list[list[float]]:
    """ 
    The corresponding equations for stages other than the first and last have 
    configurations that vary based on the stage number, but they do follow a 
    pattern that can be used to construct rows. This function creates the two 
    rows of the matrix that correspond to stage n. 
    """
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


def construct_final_stage(p_I_feed_conc: float, lambdaa: float, 
                          eq_const: float, num_stages: int) -> list[list[float]]:
    """ 
    The final two equations corresponding to the final stage will always have 
    the same configuration. This function creates the last two rows of the 
    matrix that correspond to the final stage. 
    """
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
    eq_2N_row.append(-p_I_feed_conc/(lambdaa+eq_const))
    return [eq_2N_1_row, eq_2N_row];


def construct_matrix(p_I_feed_conc: float, p_II_feed_conc: float, 
                     lambdaa: float, eq_const: float, 
                     num_stages: int) -> list[list[float]]:
    """
    Utilizes construct_stage_one, construct_stage_n, and construct_final_stage 
    to build the matrix that corresponds to the system of equations for the 
    countercurrent N-stage system. 
    """
    matrix = [];
    if num_stages == 1:
        matrix.append([-1, eq_const, 0])
        matrix.append([0, 1, (p_I_feed_conc+(lambdaa*p_II_feed_conc))/(lambdaa+eq_const)])
        return matrix
    stage_one = construct_stage_one(p_II_feed_conc, lambdaa, eq_const, num_stages)
    matrix.append(stage_one[0]) 
    matrix.append(stage_one[1])
    for i in range(2, num_stages):
        stage_n = construct_stage_n(lambdaa, eq_const, num_stages, i)
        matrix.append(stage_n[0])
        matrix.append(stage_n[1])
    final_stage = construct_final_stage(p_I_feed_conc, lambdaa, eq_const, num_stages)
    matrix.append(final_stage[0])
    matrix.append(final_stage[1])
    return matrix;



"""
Now we need the functions for reducing the matrix to reduced row echelon form. 
We start with basic row operations, and then make an algorithm to solve the 
matrix. 
""" 

def subtract_row(matrix: list[list[float]], row_a_index: int, row_b_index: int, 
                 multiplier: float) -> list[list[float]]:
    """
    Performs the row operation row_a - multiplier*row_b, and returns the 
    modified matrix. 
    """
    new_row_a = []
    for i in range(0, len(matrix[0])):
        new_row_a.append(matrix[row_a_index][i] - multiplier*matrix[row_b_index][i])
    matrix[row_a_index] = new_row_a
    return matrix 

def multiply_row(matrix: list[list[float]], row_a_index: int, 
                 multiplier: float) -> list[list[float]]:
    """ 
    Multiplies the given row_a by the multiplier: multiplier*row_a, and 
    returns the modified matrix. 
    """
    new_row_a = []
    for i in range(0, len(matrix[0])):
        new_row_a.append(multiplier*matrix[row_a_index][i]);
    matrix[row_a_index] = new_row_a
    return matrix 

def reduce_to_ref(matrix: list[list[float]]) -> list[list[float]]:
    """
    Uses the previous function subtract_row in an algorithm to reduce a matrix 
    to row echelon form. Returns the modified matrix, now in row echelon form. 
    """
    for i in range(0, len(matrix)-1):
        for j in range(i+1,len(matrix)):
            #the col index of the leading entry for i row is the same as i 
            if matrix[j][i] != 0.0:
                matrix = subtract_row(matrix, j, i, matrix[j][i]/matrix[i][i])
    return matrix

def reduce_to_rref(matrix: list[list[float]]) -> list[list[float]]:
    """ 
    Uses the previous function reduce_to_ref to get the matrix in 
    row echelon form, then uses subtract_row and multiply_row in an algorithm 
    to reduce a matrix in row_echelon form to reduced row echelon form. 
    Returns the modified matrix, now in reduced row echelon form. 
    """
    matrix = reduce_to_ref(matrix)
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

def extract_conc(rref_matrix: list[list[float]]) -> list[float]:
    """ 
    Takes the matrix in reduced row echelon form and extracts the 
    concentration values from the final row, which would be the constant side 
    of the augmented matrix, returning those concentrations in a list. 
    """
    conc_list = []
    for i in range(0, len(rref_matrix)):
        conc_list.append(rref_matrix[i][len(rref_matrix)])
    return conc_list

def display_sys_info(p_I_feed_conc: float, p_II_feed_conc: float, 
                     volumetric_flow_rate_I: float, 
                     volumetric_flow_rate_II: float, eq_const: float, 
                     num_stages: int) -> None:
    """ 
    Uses the print function to give a summary of all the 
    user-entered information
    """
    print("\nMulti-Stage Countercurrent System Information:")
    print("phase I  feed concentration:", p_I_feed_conc, " kg/mol")
    print("phase II feed concentration:", p_II_feed_conc, " kg/mol")
    print("phase I  volumetric flow rate:", volumetric_flow_rate_I, " L/min")
    print("phase II volumetric flow rate:", volumetric_flow_rate_II, "L/min")
    print("equilibrium constant:", eq_const)
    print("number of stages:", num_stages)
    
def display_concentrations(concentrations: list[float]) -> None:
    """
    Uses the print function to show the user the concentrations determined by 
    solving the system of equations. 
    """
    stage_number: int = 1
    print("\nConcentrations:")
    for i in range(0, len(concentrations)):
        if i%2 == 0:
            print("phase I  stage ", stage_number, " concentration: ", 
                  round(concentrations[i], 8), " kg/mol")
        else: 
            print("phase II stage ", stage_number, " concentration: ", 
                  round(concentrations[i], 8), " kg/mol")
            stage_number = stage_number + 1


"""
Finally, I need the functions for interacting with the user. 

A majority of these functions are for the user to input the necessary variables 
to build the matrix and solve the system of equations. 

Note: At any time, the user can exit the program by entering “e” as an input, 
and that will be reflected in the various “get input” functions. 
"""

def str_is_float(string: str) -> bool:
    """
    Loops through the given string, checking each entry in the string. If the 
    string is only composed of numbers and one period, it will be considered a 
    float, and returns True. Otherwise returns False. 
    """
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

def test_case(case_num: int) -> None:
    """
    Tests the computer-calculated solutions agree against the hand-calculated 
    solutions to the same system to ensure each concentration agrees within
    10 decimal places (to rule out negligible float errors).
    """
    test_message = "Test Case " + str(case_num) + ": " + str(case_num+1) + " Stages"
    hand_calc_conc = []
    error_message = "Test Case " + str(case_num) + " Failed. Incorrect Concentrations."
    num_stages = 0;
    if case_num == 1:
        hand_calc_conc = [129/65.0, 43/65.0, 42/65.0, 14/65.0]
        num_stages = 2
    elif case_num == 2:
        hand_calc_conc = [417/200.0, 139/200.0, 39/50.0, 13/50.0, 69/200.0, 23/200.0]
        num_stages = 3 
    elif case_num == 3:
        hand_calc_conc = [1281/605, 427/605, 498/605, 166/605, 237/605, 79/605, 30/121, 10/121]
        num_stages = 4 
    print(test_message)
    p_I_feed_conc = 0.2
    p_II_feed_conc = 2.0
    vol_flow_rate_I = 100
    vol_flow_rate_II = 100
    eq_const = 3.0
    display_sys_info(p_I_feed_conc, p_II_feed_conc, 
                               vol_flow_rate_I, vol_flow_rate_II, eq_const, 
                               num_stages)
    matrix = construct_matrix(p_I_feed_conc, p_II_feed_conc, 
                               vol_flow_rate_II/vol_flow_rate_I, eq_const, 
                               num_stages)
    rref = reduce_to_rref(matrix)
    program_conc = extract_conc(rref)
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
    print("\nThe program-calculated and hand-calculated concentrations agree to within 10 decimal places.")
    string = ""
    for i in range(0, 80):
        string = string + "-"
    print(string)
    
    
    

def startup_message() -> None:
    """
    Prints the messages seen at the start of each loop of the program.
    """
    print("Enter 'e' at any time to exit the program.")
    

def get_num_stages() -> int | str:
    """
    Retrieves either the number of stages or an exit command from the user. 
    """
    question_string = "\nHow many stages?"
    input_req = "Integers greater than 0, or 'e'"
    invalid = True 
    while invalid:
        print(question_string)
        print("Valid Inputs: " + input_req)
        num_stages = input()
        num_stages = num_stages.strip()
        if num_stages.isnumeric():
            num_stages = int(num_stages)
            if num_stages != 0:
                invalid = False
        elif num_stages == "e":
            invalid = False
        if invalid:
            print("Invalid Input")
    return num_stages
    
# You know what, we can consolidate these because the logic is the same. 
# The only thing that needs to change is the question string. 

float_questions = ["\nWhat is the feed concentration of 'A' in phase I? (kg/mol)",
                   "\nWhat is the feed concentration of 'A' in phase II? (kg/mol)",
                   "\nWhat is the volumetric flow rate of phase I? (L/min)",
                   "\nWhat is the volumetric flow rate of phase II? (L/min)",
                   "\nWhat is the equilibrium constant for this system?"]

float_reqs = ["Positive floats, 0, or 'e'",
              "Positive floats, 0, or 'e'",
              "Floats greater than 0, or 'e'",
              "Floats greater than 0, or 'e'",
              "Floats greater than 0, or 'e'"]

def get_float(question_index: int) -> float | str:
    """
    Works with a list of strings. Each element in the list is a question 
    asking the user for a piece of information, i.e. “What is the equilibrium 
    constant?” The index given corresponds to a given question. Retrieves 
    either a float or an exit command. 
    """
    # Value is the index corresponding to the entry of float_questions 
    # that we are retrieving. ex. 0 = fetching equilibrium constant. 
    question_string = float_questions[question_index]
    invalid = True 
    while invalid:
        print(question_string)
        print("Valid Inputs: " + float_reqs[question_index])
        user_input = input()
        user_input = user_input.strip()
        if str_is_float(user_input):
            user_input = float(user_input)
            if user_input != 0.0 or question_index == 0 or question_index == 1:
                invalid = False
        elif user_input == "e":
            invalid = False
        if invalid:
            print("Invalid Input")
    return user_input



"""
Finally, the logic of the program is contained within a while loop. 
"""
        
test_case(1)
test_case(2)
test_case(3)
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
    p_I_feed_conc = system_values[0]
    p_II_feed_conc = system_values[1]
    vol_flow_rate_I = system_values[2]
    vol_flow_rate_II = system_values[3]
    eq_const = system_values[4]
    try: 
        display_sys_info(p_I_feed_conc, p_II_feed_conc, 
                                   vol_flow_rate_I, vol_flow_rate_II, eq_const, 
                                   num_stages)
        lambdaa = vol_flow_rate_II/vol_flow_rate_I
        matrix = construct_matrix(p_I_feed_conc, p_II_feed_conc, lambdaa, 
                                  eq_const, num_stages)
        matrix = reduce_to_rref(matrix)
        conc_list = extract_conc(matrix)
        display_concentrations(conc_list)
    except:
        print("An unexpected error occured. Please try different inputs.")
    string = ""
    for i in range(0, 80):
        string = string + "-"
    print(string)
print("\nExiting Program")
        

















