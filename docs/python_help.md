# Instructions for running Python code from Mathematica

This Mathematica package uses Python code from the `sum_rules.py` file to expedite calculations. Please follow this guide to ensure you have correctly configured your computer to run Python code from Mathematica.

Requirements: Mathematica 12.0+, Python 3.8+  
_Note: Some slightly older versions (Mathematica 11.2+, Python 3.4+) may work, but support is limited._



## 1. Install Python
Within Mathematica, the `ExternalEvaluate` function is used to call Python code, and the `StartExternalSession` function is used to initialize a persistent Python session so that multiple calls of `ExternalEvaluate` share the same state (ex. variables, imports).

For this to be possible, you must have a Python installation on your computer. There are two routes:
* Option A (preferred): You do not already have a Python installation on your computer, or you would prefer to keep your Python installation(s) independent of Mathematica's processes.
* Option B: You would like to configure your existing Python installation(s) for Mathematica, or Option A fails to install Python.


### Option A: Wolfram-managed installation (preferred)
1. Open a Mathematica notebook and run
   ```
   StartExternalSession["Python"]
   ```
   
2. You may receive an error and be prompted to install Python. Run
   ```
   ExternalEvaluate["Python","Install"]
   ```
   to install a Wolfram-managed Miniconda Python environment.
   1. If you are using Mathematica 12.0+, the program may attempt to do this automatically.
   
3. If the Wolfram-managed installation fails, follow the steps below for Option B.


### Option B: Manual installation
1. Follow [Wolfram's guide](https://reference.wolfram.com/language/workflow/ConfigurePythonForExternalEvaluate.html) for installing Python and the required packages for using `ExternalEvaluate`. If you already have a Python installation, make sure to follow the steps to install the `pyzmq` package using the package manger `pip`. This allows Mathematica to communicate with Python.
   
2. In Mathematica, list out the Python evaluators on your computer using
   ```
   FindExternalEvaluators["Python"]
   ```
   1. You may have more than one evaluator on your system&mdash;make sure to keep track of which one you are installing packages into.
   
3. Copy the path to the Python executable (either `python.exe` or `python3`) and register it:
   ```
   RegisterExternalEvaluator[
     "Python",
     <|"SystemExecutable" -> "C:\\path\\to\\python.exe"|>
   ]
   ```
   This ensures that your computer remembers where to find the Python installation.
   1. On Windows, remember to escape the '\\' character in the string using double backslashes as shown above.
   2. If you would like to assign a name during registration, you can run the following instead (replace `environment name` with the desired name):
      ```
      RegisterExternalEvaluator[
        "Python",
        <|
          "SystemExecutable" -> "C:\\path\\to\\python.exe",
          "Name" -> "environment name"
        |>
      ]
      ```
   3. If you need to delete a registered evaluator, use the [UnregisterExternalEvaluator](https://reference.wolfram.com/language/ref/UnregisterExternalEvaluator.html) function.



## 2. Check your Python installation
1. Start an external Python session in Mathematica:
   ```
   StartExternalSession["Python"]
   ```
   
2. Use the external session to evaluate Python code:
   ```
   ExternalEvaluate["Python","import sys; sys.version"]
   ```



## 3. Install NumPy
1. The Python script requires the NumPy package to run. If you followed Option A, your installation should already include NumPy. Otherwise, you may have to manually install it. See if you can import NumPy (step 2). If it fails, enter the following on your terminal or command prompt:
   ```
   pip install numpy
   ```
   
2. In Mathematica, check that you can import NumPy:
   ```
   ExternalEvaluate["Python","import numpy as np; print(np.__version__)"]
   ```
