# Instructions for using the package

Please follow our [guide](docs/python_help.md) to ensure your computer is properly configured to run Python code from Mathematica. This only needs to be done the first time you run the package.

1. Download `ASRs.m` and `sum_rules.py` and place both files into the same directory as your Mathematica notebook.

2. Import the Mathematica package (note: the ` in the string is necessary):
   ```
   SetDirectory[NotebookDirectory[]];
   Needs["ASRs`"];
   ```

3. This Mathematica package uses Python code from the `sum_rules.py` file to expedite calculations. In your Mathematica notebook, start an external session to execute Python code, and define the path to the `sum_rules.py` Python file:
   ```
   session = StartExternalSession["Python"];
   path = "sum_rules.py";
   ```

4. Load the Python file:
   ```
   startPythonSession[session,path];
   ```

5. Generate the amplitude sum rules for a U-spin system. For example, a system with representations u = 0 in the in state, u = 1 in the Hamiltonian, and u = 1/2 and u = 1/2 in the out state is called as
   ```
   system = generateASRs[{0},{1},{1/2,1/2}];
   ```

6. View a detailed output of the system's representations, amplitudes, and amplitude sum rules:
   ```
   printSystem[system];
   ```

# Help

To view the complete list of functions and variables in this package, type in
   ```
   ?ASRs`*
   ```

For details on a function's arguments, options, and outputs, run
   ```
   ASRsHelp[function];
   ```
