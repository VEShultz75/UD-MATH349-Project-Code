These are the instructions for using the python program for calculating the concentrations in a multistage countercurrent extraction system. 

There are three methods documented for running the Python Program: 

1. Editor 
Open the python program mass_contactors_in_series_countercurrent_calculator.py in the python editor of your choice. Press the editor's button to run the program. 

2. Windows Powershell, Unix System Terminal 
python mass_contactors_in_series_countercurrent_calculator.py

3. Windows Command Prompt 
py mass_contactors_in_series_countercurrent_calculator.py 



The program will first run through the three test cases detailed in the accompanying paper. The final calculated concentrations will be compared to those calculated by hand, and if they agree to within 10 decimal places, the concentrations calculated by the program will be declared correct. If the concentrations are not correct, the program will exit and display an error message. The full concentrations both calculated by the computer and by hand will be displayed in their own sections. 

The program will then prompt the user for the parameters of the system. The user should enter 'e' at any time to exit the program. The program will ask for: the number of stages, the phase I feed concentration, phase II feed concentration, phase I volumetric flow rate, phase II volumetric flow rate, and the equilibrium constant, in that order. The program will indicate if the input was invalid, and allow the user to input a new value. 

Once all parameters are specified, the parameter information will be summarized and printed to the terminal. The concentrations for each stage and each phase will be calculated and displayed below the system information. 

The user will then be able to input the parameter information for a new system. 