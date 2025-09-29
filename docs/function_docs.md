# Full documentation of functions

## List of functions

- [ASRsHelp](#ASRsHelp): ASRsHelp[function] prints extended documentation on a function's arguments, options, and outputs.

## Documentation
### ASRsHelp
ASRsHelp[function] prints extended documentation on a function's arguments, options, and outputs.

Arguments:
function (Symbol): The name of a function from the ASRs Mathematica package

### startPythonSession
startPythonSession[session,path] checks for a valid Python session and file path and loads in the Python file.

Arguments:
session (ExternalSessionObject): An active external Python session
path (String): Path to the sum_rules.py Python file

Returns:
An ExternalFunction object indicating the sum_rules.py file has been loaded into the Python session

### generateASRs
generateASRs[in,h,out] finds all amplitudes and amplitude sum rules (ASRs) for a given system.

Arguments:
in (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the incoming state
h (List): Contains U-spins (Real) OR CKM factors (List of Symbols) in the Hamiltonian
out (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the outgoing state

Options:
phys (True|False): Indicates whether function arguments contain U-spins (False) or particle multiplets/CKM factors (True). Default: phys->False.

Returns:
system (Association): All information about the system's representations, amplitudes, and amplitude sum rules. Keys and values:
- \"n doublets\" (Real): Number of would-be doublets
- \"Irreps\" (List): Inputted U-spin representations (List of Reals) in {incoming reps, H reps, outgoing reps} format
- \"n amps\" (Real): Number of amplitudes in the system
- \"Amplitudes\" (Association): Contains all amplitudes in the system. Keys and values:
	- \"Processes\" (List): Contains physical processes constructed from particle names (String) for a process and its U-spin conjugate process. Only available for physical systems.
	- \"QNs\" (List): Contains processes written using m quantum numbers (String), where m is the third component of U-spin, for a process and its U-spin conjugate process
	- \"n-tuple\" (String): Representation of an amplitude pair using a comma-separated tuple of substrings. Each substring is comprised of '-'s and '+'s and encodes u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components. By convention, we use one n-tuple (beginning with a '-' sign) to represent an amplitude pair: an amplitude and its U-spin conjugate amplitude.
	- \"Node\" (String): Representation of an amplitude pair using the coordinate notation of the lattice used to derive sum rules
	- \"Binary indices\" (List): Contains decimal indices (Real) derived from converting the n-tuple and its U-spin conjugate into binary numbers through '-' <-> 0 and '+' <-> 1
	- \"q factor\" (Real): (-1)^q_i factor for an amplitude pair accounting for movements of representations between the in state, Hamiltonian, and out state
	- \"p factor\" (Real): (-1)^p parity factor for the system
	- \"mu\" (Real): mu-factor for the node in the lattice used to derive sum rules
	- \"CG\" (Real): Clebsch-Gordan coefficient from symmetrization for systems without doublets. Equal to 1 for all amplitudes for a system with at least one doublet.
	- \"CKM\" (List): Contains weak interaction factors (Real) from the Hamiltonian. Only appears for physical systems.
- \"n ASRs\" (List): Contains the number of amplitude sum rules (Real) at each order of breaking
- \"ASRs\" (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking

### printSystem

printSystem[system] prints information about the system's representations, amplitudes, and amplitude sum rules and modifies the system to include formatted sum rules.

Arguments:
system (Association): A system association. See the documentation for generateASRs for details.

Options:
showReps (True|False): Default: showReps->True.

showAmps (True|False): Default: showAmps->True
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.

showASRs (True|False): Default: showASRs->True.
takeProd (True|False): Indicates whether to write each row of a sum rules matrix as an algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: takeProd->True.
ampFormat (String): Specified format for displaying amplitudes. Options are a/s-type amplitudes with n-tuples (\"a/s n-tuple\"), a/s-type amplitudes with numbered indices (\"a/s indices\"), a/s-type amplitudes with nodes (\"a/s nodes\"), A amplitudes with numbered indices (\"A indices\"), and A amplitudes with physical process (\"A physical\"). Default: ampFormat->\"a/s n-tuple\".
CKM (True|False): Indicates whether to include CKM factors in the sum rules (True) or not (False). Default: CKM->False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
system (Association): The inputted system association, modified to include formatted amplitude sum rules. New/modified keys and values:
- \"Formatted ASRs\" (List): Contains Lists/matrices of amplitude sum rules/amplitude sum rule coefficients formatted according to sum rule printing settings.

### numAmps
numAmps[amplitudes] returns the total number of amplitudes in the system.

Arguments:
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.

Options:
nPairs (True|False): Indicates whether to return the number of amplitudes (False) or amplitude pairs (True). Default: nPairs: False.

Returns:
Total number of amplitudes or amplitude pairs in the system (Real)

### labelAmps
labelAmps[amplitudes,colName,labels] adds a column of labels to amplitudes.

Arguments:
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.
colName (String): Name of new column to add to amplitudes association
labels (List): Contains labels (Any) for each amplitude or amplitude pair. Number of labels must equal number of amplitudes or amplitude pairs.

Options:
labeling (String): Labeling mode to indicate whether user is labeling single amplitudes (\"Amplitudes\") or amplitude pairs (\"Amplitude pairs\"). Default: labeling -> \"Amplitudes\".

Returns:
amplitudes (Association): The inputted amplitudes association, modified to include a new column of user-defined labels. New/modified keys and values:
- colName (List|Any): Either contains user-defined labels (Any) for an amplitude and its U-spin conjugate or a label (Any) for the amplitude pair

### unlabelAmps
unlabelAmps[amplitudes,colNames] removes columns from amplitudes.

Arguments:
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.
colNames (String|List): Name(s) of column(s) to remove from amplitudes association

Returns:
amplitudes (Association): The inputted amplitudes association, modified to remove specified columns

### printAmps
printAmps[amplitudes] prints a table of amplitudes.

Arguments:
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.

Options:
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.

### numASRs
numASRs[ASRs] returns the number of amplitude sum rules at each order of breaking.

Arguments:
ASRs (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking

Options:
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
Number of amplitude sum rules at each order of breaking (List of Reals)

### printASRs

printASRs[ASRs,amplitudes] prints amplitude sum rules at each order of breaking.

Arguments:
ASRs (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.

Options:
showASRs (True|False): Default: showASRs->True.
takeProd (True|False): Indicates whether to write each row of a sum rules matrix as an algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: takeProd->True.
ampFormat (String): Specified format for displaying amplitudes. Options are a/s-type amplitudes with n-tuples (\"a/s n-tuple\"), a/s-type amplitudes with numbered indices (\"a/s indices\"), a/s-type amplitudes with nodes (\"a/s nodes\"), A amplitudes with numbered indices (\"A indices\"), and A amplitudes with physical process (\"A physical\"). Default: ampFormat->\"a/s n-tuple\".
CKM (True|False): Indicates whether to include CKM factors in the sum rules (True) or not (False). Default: CKM->False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
\"Formatted ASRs\" (List): Contains Lists/matrices of amplitude sum rules/amplitude sum rule coefficients formatted according to sum rule printing settings
