# Instructions for using the FlaSR package

Please follow our [guide](docs/python_help.md) to ensure your computer is properly configured to run Python code from Mathematica. This only needs to be done the first time you run the package.

1. Download `FlaSR.m` and `FlaSR.py` and place both files into the same directory as your Mathematica notebook.

2. Import the Mathematica package (note: the ` in the string is necessary):
   ```
   SetDirectory[NotebookDirectory[]];
   Needs["FlaSR`"];
   ```

3. This Mathematica package uses Python code from the `FlaSR.py` file to expedite calculations. In your Mathematica notebook, start an external session to execute Python code, and define the path to the `FlaSR.py` Python file:
   ```
   session = StartExternalSession["Python"];
   path = "FlaSR.py";
   ```

4. Load the Python file:
   ```
   startPythonSession[session,path];
   ```

5. Generate the sum rules for a U-spin system. For example, a system with representations u = 0 in the in state, u = 1 in the Hamiltonian, and u = 1/2 and u = 1/2 in the out state is called as
   ```
   system = generateSRs[{0},{1},{1/2,1/2}];
   ```

6. View a detailed output of the system's representations, amplitudes, and sum rules (refer to the paper for a detailed explanation of amplitude bases):
   ```
   printSystem[system,ampType->{a,s},amp2Type->{Δ,Σ}];
   ```
   The program will output information about the system's representations, amplitudes, amplitude sum rules, and squared amplitude sum rules. A sum rule holding up to order b is broken at order b+1, i.e. by corrections of O(ε<sup>b+1</sup>), where ε is the small symmetry-breaking parameter (~0.3 for U-spin).

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
