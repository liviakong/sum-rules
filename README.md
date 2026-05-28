# Instructions for using the FlaSR package

Please follow our [guide](docs/python_help.md) to ensure your computer is properly configured to run Python code from Mathematica. This only needs to be done the first time you run the package. You can use our package by downloading the `example.nb` file and running the cells there or by following the steps below.

1. Download `FlaSR.m` and `FlaSR.py` and place both files into the same directory as your Mathematica notebook.

2. Import the Mathematica package:
   ```
   SetDirectory[NotebookDirectory[]];
   Get["FlaSR.m"];
   ```

3. This Mathematica package uses Python code from the `FlaSR.py` file to perform part of the calculations. In your Mathematica notebook, start an external session to execute Python code and load the Python file into the session:
   ```
   session = StartExternalSession["Python"];
   startPythonSession[session, "FlaSR.py"]
   ```

4. Generate the sum rules for a U-spin system. For example, the representations for a system with a singlet (u = 0) in the in state, a triplet (u = 1) Hamiltonian, and two doublets (u = 1/2, u = 1/2) in the final state are inputted as
   ```
   system = generateSRs[{0},{1},{1/2,1/2}];
   ```

5. View a detailed output of the system's representations, amplitudes, and sum rules (refer to the paper for an explanation of amplitude bases and formatting options):
   ```
   printSystem[system,ampType->{a,s},amp2Type->{Δ,Σ}];
   ```
   The program will print information about the system's representations, amplitudes, amplitude sum rules, and amplitude-squared sum rules. A sum rule displayed at order b is interpreted to be broken by corrections of the next order, O(ε<sup>b+1</sup>), where ε is the small symmetry-breaking parameter (~0.3 for U-spin).

# Help

To view the complete list of functions and variables in this package, type in
   ```
   ?FlaSR`*
   ```

For details on a function's arguments, options, and outputs, run
   ```
   FlaSRHelp[function];
   ```
This documentation is also available on our [function documentation](docs/function_docs.md) page.
