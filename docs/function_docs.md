# Full documentation of functions

## List of functions

- [FlaSRHelp](#FlaSRHelp): FlaSRHelp[function] prints extended documentation on a function's arguments, options, and outputs.

## Documentation
### FlaSRHelp
FlaSRHelp[function] prints extended documentation on a function's arguments, options, and outputs.

Arguments:
- function (Symbol): The name of a function from the FlaSR Mathematica package.

### startPythonSession
startPythonSession[session,path] checks for a valid Python session and file path and loads in the Python file.

Arguments:
- session (ExternalSessionObject): An active external Python session.
- path (String): Path to the FlaSR.py Python file.

Returns:
- An ExternalFunction object indicating the FlaSR.py file has been loaded into the Python session.

### generateASRs
generateASRs[in,h,out] finds all amplitudes and amplitude sum rules (ASRs) for a given system.

Arguments:
- in (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the incoming state.
- h (List): Contains U-spins (Real) OR CKM factors (List of Symbols) in the Hamiltonian.
- out (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the outgoing state.

Options:
- phys (True|False): Indicates whether function arguments contain U-spins (False) or particle multiplets/CKM factors (True). Default: phys->False.

Returns:
- system (Association): All information about the system's representations, amplitudes, and ASRs. Keys and values:
	- \"Irreps\" (List): Inputted U-spin representations (List of Reals) in {incoming reps, H reps, outgoing reps} format.
	- \"n doublets\" (Real): Number of would-be doublets.
	- \"p factor\" (Real): (-1)^p parity factor for the system determining forms of a/s-type amplitudes.
	- \"n amps\" (Real): Number of amplitudes in the system.
	- \"Amplitudes\" (Association): Contains all amplitudes in the system. Keys and values:
		- \"Processes\" (List): Contains physical processes constructed from particle names (String) for a process and its U-spin conjugate process. Only available for physical systems.
		- \"QNs\" (List): Contains processes written using m quantum numbers (String), where m is the third component of U-spin, for a process and its U-spin conjugate process.
		- \"n-tuples\" (List): Contains representations (String) of an amplitude and its U-spin conjugate amplitude using comma-separated tuples of substrings. Each substring is comprised of '-'s and '+'s and encodes u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components.
		- \"Coords\" (String): Representation of an amplitude pair using the coordinate notation of the lattice used to derive sum rules.
		- \"Binary indices\" (List): Contains indices (Real), written in base 10, derived from converting the n-tuple and its U-spin conjugate into binary numbers through '-' <-> 0 and '+' <-> 1.
		- \"mu\" (Real): mu-factor for the coord in the lattice used to derive sum rules.
		- \"CG\" (Real): Clebsch-Gordan coefficient from symmetrization for systems without doublets. Equal to 1 for all amplitudes for a system with at least one doublet.
		- \"CKM\" (List): Contains weak interaction factors (Real) from the Hamiltonian. Only appears for physical systems.
	- \"n ASRs\" (List): Contains the number of amplitude sum rules (Real) at each order of breaking.
	- \"ASRs\" (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.

### generateSRs
generateSRs[in,h,out] finds all amplitudes, amplitude sum rules (ASRs), and squared amplitude sum rules (A2SRs) for a given system.

Arguments:
- in (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the incoming state.
- h (List): Contains U-spins (Real) OR CKM factors (List of Symbols) in the Hamiltonian.
- out (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the outgoing state.

Options:
- phys (True|False): Indicates whether function arguments contain U-spins (False) or particle multiplets/CKM factors (True). Default: phys->False.

Returns:
- system (Association): All information about the system's representations, amplitudes, ASRs, and A2SRs. Keys and values:
	- \"Irreps\" (List): Inputted U-spin representations (List of Reals) in {incoming reps, H reps, outgoing reps} format.
	- \"n doublets\" (Real): Number of would-be doublets.
	- \"p factor\" (Real): (-1)^p parity factor for the system determining forms of a/s-type amplitudes.
	- \"n amps\" (Real): Number of amplitudes in the system.
	- \"Amplitudes\" (Association): Contains all amplitudes in the system. Keys and values:
		- \"Processes\" (List): Contains physical processes constructed from particle names (String) for a process and its U-spin conjugate process. Only available for physical systems.
		- \"QNs\" (List): Contains processes written using m quantum numbers (String), where m is the third component of U-spin, for a process and its U-spin conjugate process.
		- \"n-tuples\" (List): Contains representations (String) of an amplitude and its U-spin conjugate amplitude using comma-separated tuples of substrings. Each substring is comprised of '-'s and '+'s and encodes u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components.
		- \"Coords\" (String): Representation of an amplitude pair using the coordinate notation of the lattice used to derive sum rules.
		- \"Binary indices\" (List): Contains indices (Real), written in base 10, derived from converting the n-tuple and its U-spin conjugate into binary numbers through '-' <-> 0 and '+' <-> 1.
		- \"mu\" (Real): mu-factor for the coord in the lattice used to derive sum rules.
		- \"CG\" (Real): Clebsch-Gordan coefficient from symmetrization for systems without doublets. Equal to 1 for all amplitudes for a system with at least one doublet.
		- \"CKM\" (List): Contains weak interaction factors (Real) from the Hamiltonian. Only appears for physical systems.
	- \"n ASRs\" (List): Contains the number of amplitude sum rules (Real) at each order of breaking.
	- \"ASRs\" (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.
	- \"n A2SRs\" (List): Contains the number of squared amplitude sum rules (Real) at each order of breaking.
	- \"A2SRs\" (List): Contains matrices of squared amplitude sum rule coefficients (Real) corresponding to each order of breaking.

### numAmps
numAmps[system,nPairs:False] returns the total number of amplitudes in the system.

Arguments:
- system (Association): A system association. See the documentation for generateSRs for details.

Options:
- nPairs (True|False): Indicates whether to return the number of amplitudes (False) or amplitude pairs (True). Default: nPairs: False.

Returns:
- The total number of amplitudes or amplitude pairs in the system (Real).

### labelAmps
labelAmps[system,colName,labels] modifies system to add a column of user-defined labels to system[[\"Amplitudes\"]].

Arguments:
- system (Association): A system association. See the documentation for generateSRs for details.
- colName (String): Name of new column to add to amplitudes association.
- labels (List): Contains labels (Any) for each amplitude or amplitude pair. Number of labels must equal number of amplitudes or amplitude pairs.

Options:
- labeling (String): Labeling mode to indicate whether user is labeling single amplitudes (\"Amplitudes\") or amplitude pairs (\"Amplitude pairs\"). Default: labeling -> \"Amplitudes\".

Returns:
- system[[\"Amplitudes\"]] (Association): The modified amplitudes association which includes a new column of user-defined labels. New/modified keys and values of amplitudes:
	- colName (List|Any): Either contains user-defined labels (Any) for an amplitude and its U-spin conjugate or a single label (Any) for an amplitude pair.

### unlabelAmps
unlabelAmps[system,colNames] modifies system to remove columns from system[[\"Amplitudes\"]].

Arguments:
- system (Association): A system association. See the documentation for generateSRs for details.
- colNames (String|List): Name(s) of column(s) to remove from amplitudes association.

Returns:
- system[[\"Amplitudes\"]] (Association): The modified amplitudes association from which the specified columns have been removed.

### printAmps
printAmps[system] prints the system's amplitudes and a/s, Δ/Σ amplitude definitions.

Arguments:
- system (Association): A system association. See the documentation for generateSRs for details.

Options:
- showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.

### numSRs
numSRs[system,squared:False] returns the number of amplitude or squared amplitude sum rules at each order of breaking.

Arguments:
- system (Association): A system association. See the documentation for generateSRs for details.

Options:
- squared (True|False): Indicates whether to return number of sum rules of amplitudes (False) or squared amplitudes (True). Default: squared: False.
- b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
- Number of (squared) amplitude sum rules at each order of breaking (List of Reals).

### printSRs
printSRs[system,ampType->{a,s}/{A} OR amp2Type->{Δ,Σ}/{A}] prints amplitude or squared amplitude sum rules at each order of breaking.

Arguments:
- system (Association): A system association. See the documentation for generateSRs for details.

Options:
- ampType (List): Contains 1 or 2 symbol(s) to select amplitude type. Convention is to set ampType->{A} or {a,s} for A amplitudes or a/s amplitudes. Default: ampType->None. Note: only one of ampType or amp2Type should be specified to print either ASRs or A2SRs; if both are specified, printSRs will print A2SRs by default.
- amp2Type (List): Contains 1 or 2 symbol(s) to select squared amplitude type. Convention is to set amp2Type->{A} or {Δ,Σ} for |A|^2 amplitudes or Δ/Σ. Default: amp2Type->None. Note: only one of ampType or amp2Type should be specified to print either ASRs or A2SRs; if both are specified, printSRs will print A2SRs by default.
- ampFormat (String): Specified format for displaying amplitudes. Options are physical process names (\"Processes\", only available for A amps), quantum numbers (\"QNs\"), n-tuples (\"n-tuples\"), coordinate notation (\"Coords\", only available for a/s amps), numbered indices (\"Binary indices\"), or user-defined labels for a column of the amplitudes table (name of column of amplitudes table containing custom labels). Default: ampFormat->\"n-tuples\" unless the system is a physical system, in which case ampFormat->\"Processes\".
- showSRs (True|False): Indicates whether to print sum rules. Default: showSRs->True.
- expandSRs (True|False): Indicates whether to display each row of a sum rules matrix as an algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandSRs->False.
- CKM (True|False): Indicates whether to include CKM factors in the sum rules (True) or not (False). Default: CKM->False.
- b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
- system[[\"SR extract\"]] (List): A new key in the system association containing Lists of (squared) amplitude sum rules/matrices of (squared) amplitude sum rule coefficients formatted according to the sum rule printing settings.
- Other new/modified keys and values of system:
	- \"Amp vector\" (List): Either is a vector of formatted A amplitudes (Symbols) or contains a formatted vector for a and s amplitudes (List of Symbols) for the system.

### printSystem
printSystem[system,ampType->{a,s}/{A}, amp2Type->{Δ,Σ}/{A}] prints information about the system's representations, amplitudes, and sum rules and modifies the system to include formatted sum rules.

Arguments:
- system (Association): A system association. See the documentation for generateSRs for details.

Options:
- showReps (True|False): Indicates whether to print information about the system's representations. Default: showReps->True.
---
- showAmps (True|False): Indicates whether to print information about the system's amplitudes. Default: showAmps->True.
- showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.
---
- showASRs (True|False): Indicates whether to print ASRs. Default: showASRs->True. Note: while both ampType and amp2Type can be separately specified, the other formatting options (e.g., ampFormat) will be shared for printing both ASRs and A2SRs.
- showA2SRs (True|False): Indicates whether to print A2SRs. Default: showA2SRs->False. Note: while both ampType and amp2Type can be separately specified, the other formatting options (e.g., ampFormat) will be shared for printing both ASRs and A2SRs.
- ampType (List): Contains 1 or 2 symbol(s) to select amplitude type. Convention is to set ampType->{A} or {a,s} for A amplitudes or a/s amplitudes. Default: ampType->None.
- amp2Type (List): Contains 1 or 2 symbol(s) to select squared amplitude type. Convention is to set amp2Type->{A} or {Δ,Σ} for |A|^2 amplitudes or Δ/Σ amplitudes. Default: amp2Type->None.
- ampFormat (String): Specified format for displaying amplitudes. Options are physical process names (\"Processes\", only available for A amps), quantum numbers (\"QNs\"), n-tuples (\"n-tuples\"), coordinate notation (\"Coords\", only available for a/s amps), numbered indices (\"Binary indices\"), or user-defined labels for a column of the amplitudes table (name of column of amplitudes table containing custom labels). Default: ampFormat->\"n-tuples\" unless the system is a physical system, in which case ampFormat->\"Processes\".
- showSRs (True|False): Indicates whether to print sum rules. Default: showSRs->True.
- expandSRs (True|False): Indicates whether to display each row of a sum rules matrix as an algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandSRs->False.
- CKM (True|False): Indicates whether to include CKM factors in the sum rules (True) or not (False). Default: CKM->False.
- b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= - b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
- system (Association): The inputted system association, modified to include formatted sum rule coefficient matrices and amplitude vector(s). New/modified keys and values:
	- \"SR extract\" (List): Contains Lists of (squared) amplitude sum rules/matrices of (squared) amplitude sum rule coefficients formatted according to the sum rule printing settings.
	- \"Amp vector\" (List): Either is a vector of formatted A amplitudes (Symbols) or contains a formatted vector for a and s amplitudes (List of Symbols) for the system.
