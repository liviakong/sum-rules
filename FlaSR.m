(* ::Package:: *)

BeginPackage["FlaSR`"];


FlaSRHelp::usage="FlaSRHelp[function] prints extended documentation on a function's arguments, options, and outputs.";
FlaSRHelp::details=
"Arguments:
function (Symbol): The name of a function from the FlaSR Mathematica package."

startPythonSession::usage="startPythonSession[session,path] checks for a valid Python session and file path and loads in the Python file.";
startPythonSession::details=
"Arguments:
session (ExternalSessionObject): An active external Python session.
path (String): Path to the FlaSR.py Python file.

Returns:
An ExternalFunction object indicating the FlaSR.py file has been loaded into the Python session."

generateASRs::usage="generateASRs[in,h,out] finds all amplitudes and amplitude sum rules (ASRs) for a given system.";
generateASRs::details=
"Arguments:
in (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the incoming state.
h (List): Contains U-spins (Real) OR coefficients (List of Symbols) in the Hamiltonian.
out (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the outgoing state.

Options:
phys (True|False): Indicates whether function arguments contain U-spins (False) or particle multiplets/CKM factors (True). Default: phys->False.

Returns:
system (Association): All information about the system's representations, amplitudes, ASRs, and A2SRs. Keys and values:
- \"Irreps\" (List): Inputted U-spin representations (List of Reals) in {{in reps}, {H rep}, {out reps}} format.
- \"Multiplets\" (List): Inputted multiplets (List of Strings) and factors (List of Symbols) in {{in multiplets}, {H factors}, {out multiplets}} format for physical systems. Empty for group-theoretic systems.
- \"n doublets\" (Real): Number of would-be doublets.
- \"p factor\" (Real): (-1)^p factor for defining a/s-type amplitudes.
- \"n amps\" (Real): Number of amplitudes in the system.
- \"Amplitudes\" (List): Contains an Association for each amplitude pair in the system. Keys and values:
	- \"Processes\" (List): Contains physical processes (String) for an amplitude and its U-spin conjugate. Only available for physical systems.
	- \"QNs\" (List): Contains m quantum number labels (String), where m is the third component of U-spin, for an amplitude and its U-spin conjugate.
	- \"n-tuples\" (List): Contains n-tuple labels (String) for an amplitude and its U-spin conjugate. n-tuples represent amplitudes as comma-separated tuples of substrings, where each substring is comprised of '-'s and '+'s and encodes the u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components.
	- \"Coords\" (String): Coordinate (String) for an amplitude pair in the lattice used to derive sum rules.
	- \"Binary indices\" (List): Contains indices (Real), written in base 10, for an amplitude and its U-spin conjugate. Indices are derived by converting the n-tuples into binary numbers through '-' <-> 0 and '+' <-> 1 and removing commas.
	- \"mu\" (Real): mu-factor for the coordinate in the lattice used to derive sum rules.
	- \"CG\" (Real): Clebsch-Gordan coefficient from symmetrization for systems without doublets. Equal to 1 for all amplitudes for a system with at least one doublet.
	- \"CKM\" (List): Contains weak interaction factors (Real) from the Hamiltonian. Only appears for physical systems.
- \"n ASRs\" (List): Contains the number of amplitude sum rules (Real) at each order of breaking.
- \"ASRs\" (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.";

findA2SRMat::usage="findA2SRMat[ASRMat] finds the A2SR matrix for a given ASR matrix.";
findA2SRMat::usage=
"Arguments:
ASRMat (List): Matrix of ASR coefficients.

Returns:
A2SRMat (List): Matrix of A2SR coefficients derived from the ASR matrix. Note: this assumes differential observables.";

generateSRs::usage="generateSRs[in,h,out] finds all amplitudes, amplitude sum rules (ASRs), and amplitude-squared sum rules (A2SRs) for a given system.";
generateSRs::details=
"Arguments:
in (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the incoming state.
h (List): Contains U-spins (Real) OR coefficients (List of Symbols) in the Hamiltonian.
out (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the outgoing state.

Options:
phys (True|False): Indicates whether function arguments contain U-spins (False) or particle multiplets/CKM factors (True). Default: phys->False.
obs (\"Diff\"|\"Int\"): Indicates whether A2SRs should be treated as the sum rules between differential (\"Diff\") or integrated (\"Int\") observables. Only in effect when phys->True. Default: obs->\"Diff\".

Returns:
system (Association): All information about the system's representations, amplitudes, ASRs, and A2SRs. Keys and values:
- \"Irreps\" (List): Inputted U-spin representations (List of Reals) in {{in reps}, {H rep}, {out reps}} format.
- \"Multiplets\" (List): Inputted multiplets (List of Strings) and factors (List of Symbols) in {{in multiplets}, {H factors}, {out multiplets}} format for physical systems. Empty for group-theoretic systems.
- \"n doublets\" (Real): Number of would-be doublets.
- \"p factor\" (Real): (-1)^p factor for defining a/s-type amplitudes.
- \"n amps\" (Real): Number of amplitudes in the system.
- \"Amplitudes\" (List): Contains an Association for each amplitude pair in the system. Keys and values:
	- \"Processes\" (List): Contains physical processes (String) for an amplitude and its U-spin conjugate. Only available for physical systems.
	- \"QNs\" (List): Contains m quantum number labels (String), where m is the third component of U-spin, for an amplitude and its U-spin conjugate.
	- \"n-tuples\" (List): Contains n-tuple labels (String) for an amplitude and its U-spin conjugate. n-tuples represent amplitudes as comma-separated tuples of substrings, where each substring is comprised of '-'s and '+'s and encodes the u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components.
	- \"Coords\" (String): Coordinate (String) for an amplitude pair in the lattice used to derive sum rules.
	- \"Binary indices\" (List): Contains indices (Real), written in base 10, for an amplitude and its U-spin conjugate. Indices are derived by converting the n-tuples into binary numbers through '-' <-> 0 and '+' <-> 1 and removing commas.
	- \"mu\" (Real): mu-factor for the coordinate in the lattice used to derive sum rules.
	- \"CG\" (Real): Clebsch-Gordan coefficient from symmetrization for systems without doublets. Equal to 1 for all amplitudes for a system with at least one doublet.
	- \"CKM\" (List): Contains weak interaction factors (Real) from the Hamiltonian. Only appears for physical systems.
	- \"Multiplet indices\" 
- \"n ASRs\" (List): Contains the number of amplitude sum rules (Real) at each order of breaking.
- \"ASRs\" (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.
- \"Unique amp pairs\" (List): Contains indices (Real) of amplitude pairs corresponding to unique channels after integration. Only appears for physical systems with obs->\"Int\".
- \"n A2SRs\" (List): Contains the number of amplitude-squared sum rules (Real) at each order of breaking.
- \"A2SRs\" (List): Contains matrices of amplitude-squared sum rule coefficients (Real) corresponding to each order of breaking.
- \"SR extract\" (List): Contains sum rule coefficient matrices at the selected b. If only ampType (amp2Type) is specified, contains only ASR (A2SR) coefficients. If both ampType and amp2Type are specified, contains A2SR coefficients. Initialized to None by generateSRs and redefined after running printSystem.
- \"Amp vector\" (List): Either is a vector of formatted A-type amplitudes (Symbols) or contains formatted vectors for a- and s-type amplitudes (List of Symbols) for the system. Initialized to None by generateSRs and redefined after running printSystem.";

numAmps::usage="numAmps[system,nPairs:False] returns the total number of amplitudes in the system.";
numAmps::details=
"Arguments:
system (Association): A system association. See the documentation for generateSRs for details.

Options:
nPairs (True|False): Indicates whether to return the number of amplitudes (False) or amplitude pairs (True). Default: nPairs: False.

Returns:
The total number of amplitudes or amplitude pairs in the system (Real).";

labelAmps::usage="labelAmps[system,colName,labels] modifies system to add a column of user-defined labels to system[[\"Amplitudes\"]]."
labelAmps::details=
"Arguments:
system (Association): A system association. See the documentation for generateSRs for details.
colName (String): Name of new column to add to amplitudes table.
labels (List): Contains labels (String) for each amplitude or amplitude pair. Number of labels must equal number of amplitudes or amplitude pairs.

Options:
labeling (String): Labeling mode to indicate whether user is labeling single amplitudes (\"Amplitudes\") or amplitude pairs (\"Amplitude pairs\"). Default: labeling -> \"Amplitudes\".

Returns:
system[[\"Amplitudes\"]] (List): The modified amplitudes table which includes a new column of user-defined labels. New/modified keys and values of amplitudes associations:
- colName (List|Any): Either contains user-defined labels (Any) for an amplitude and its U-spin conjugate or a single label (Any) for an amplitude pair.";

unlabelAmps::usage="unlabelAmps[system,colNames] modifies system to remove columns from system[[\"Amplitudes\"]]."
unlabelAmps::details=
"Arguments:
system (Association): A system association. See the documentation for generateSRs for details.
colNames (String|List): Name(s) of column(s) to remove from amplitudes table.

Returns:
system[[\"Amplitudes\"]] (List): The modified amplitudes table from which the specified column(s) have been removed.";

printAmps::usage="printAmps[system] prints the system's amplitudes and definitions for a/s-type amplitudes and \[CapitalDelta]/\[CapitalSigma]-type amplitudes-squared.";
printAmps::details=
"Arguments:
system (Association): A system association. See the documentation for generateSRs for details.

Options:
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Default: showFactors->False.";

numSRs::usage="numSRs[system,squared:False] returns the number of amplitude or amplitude-squared sum rules at each order of breaking.";
numSRs::details=
"Arguments:
system (Association): A system association. See the documentation for generateSRs for details.

Options:
squared (True|False): Indicates whether to return counts of amplitude (False) or amplitude-squared (True) sum rules. Default: squared: False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real, 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
Number of amplitude (or amplitude-squared) sum rules at each order of breaking (List of Reals).";

printSRs::usage="printSRs[system,ampType->{a,s}/{A} OR amp2Type->{\[CapitalDelta],\[CapitalSigma]}/{A}] prints amplitude or amplitude-squared sum rules at each order of breaking, returns the printed matrices, and modifies the amplitude vector(s) in the system association.";
printSRs::details=
"Arguments:
system (Association): A system association. See the documentation for generateSRs for details.

Options:
ampType (List): Contains 1 or 2 symbol(s) to select amplitude type. Convention is to set ampType->{A} for A-type amplitudes and ampType->{a,s} for a/s-type amplitudes. Default: ampType->None. Note: only one of ampType or amp2Type should be specified to print either ASRs or A2SRs; if both are specified, printSRs will print A2SRs by default.
amp2Type (List): Contains 1 or 2 symbol(s) to select squared amplitude type. Convention is to set amp2Type->{A} for |A\!\(\*SuperscriptBox[\(|\), \(2\)]\) amplitudes-squared and amp2Type->{\[CapitalDelta],\[CapitalSigma]} for \[CapitalDelta]/\[CapitalSigma]-type amplitudes-squared. Default: amp2Type->None. Note: only one of ampType or amp2Type should be specified to print either ASRs or A2SRs; if both are specified, printSRs will print A2SRs by default.
ampFormat (String): Labeling convention for displaying amplitudes. Options are physical process names (\"Processes\", only available for A-type amps), m quantum numbers (\"QNs\"), n-tuples (\"n-tuples\"), coordinate notation (\"Coords\", only available for a/s-type amps), numbered indices (\"Binary indices\"), or user-defined labels for a column of the amplitudes table (name of custom column). Default: ampFormat->\"n-tuples\" unless the system is a physical system, in which case ampFormat->\"Processes\".
showSRs (True|False): Indicates whether to print sum rules. Default: showSRs->True.
expandSRs (True|False): Indicates whether to display each row of a sum rules matrix as an expanded algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandSRs->False.
CKM (True|False): Indicates whether to include CKM factors in the sum rules. Default: CKM->False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real, 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.
amp2Quad (True|False): Indicates whether the symbol in amp2Type is quadratically (True) or linearly (False) dependent on A. Default: amp2Quad->False.

Returns:
system[[\"SR extract\"]] (List): Modified value of the \"SR extract\" key to the system association. Contains sum rule coefficient matrices at the selected b. If only ampType (amp2Type) is specified, contains only ASR (A2SR) coefficients. If both ampType and amp2Type are specified, contains A2SR coefficients.
Other modified keys and values of the system assocation:
- \"Amp vector\" (List): Either is a vector of formatted A-type amplitudes (Symbols) or contains formatted vectors for a- and s-type amplitudes (List of Symbols) for the system.";

printSystem::usage="printSystem[system,ampType->{a,s}/{A}, amp2Type->{\[CapitalDelta],\[CapitalSigma]}/{A}] prints information about the system's representations, amplitudes, and sum rules and adds formatted sum rules and amplitude vectors to the system association.";
printSystem::details=
"Arguments:
system (Association): A system association. See the documentation for generateSRs for details.

Options:
showReps (True|False): Indicates whether to print the system summary. Default: showReps->True.
-----
showAmps (True|False): Indicates whether to print the amplitudes table. Default: showAmps->True.
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Default: showFactors->False.
-----
showSRs (True|False): Indicates whether to print sum rules. Default: showSRs->True.
ampType (List): Contains 1 or 2 symbol(s) to select amplitude type. Convention is to set ampType->{A} for A-type amplitudes and ampType->{a,s} for a/s-type amplitudes. Default: ampType->None. Note: only one of ampType or amp2Type should be specified to print either ASRs or A2SRs; if both are specified, printSRs will print A2SRs by default.
amp2Type (List): Contains 1 or 2 symbol(s) to select squared amplitude type. Convention is to set amp2Type->{A} for |A\!\(\*SuperscriptBox[\(|\), \(2\)]\) amplitudes-squared and amp2Type->{\[CapitalDelta],\[CapitalSigma]} for \[CapitalDelta]/\[CapitalSigma]-type amplitudes-squared. Default: amp2Type->None. Note: only one of ampType or amp2Type should be specified to print either ASRs or A2SRs; if both are specified, printSRs will print A2SRs by default.
ampFormat (String): Labeling convention for displaying amplitudes. Options are physical process names (\"Processes\", only available for A-type amps), m quantum numbers (\"QNs\"), n-tuples (\"n-tuples\"), coordinate notation (\"Coords\", only available for a/s-type amps), numbered indices (\"Binary indices\"), or user-defined labels for a column of the amplitudes table (name of custom column). Default: ampFormat->\"n-tuples\" unless the system is a physical system, in which case ampFormat->\"Processes\".
expandSRs (True|False): Indicates whether to display each row of a sum rules matrix as an expanded algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandSRs->False.
CKM (True|False): Indicates whether to include CKM factors in the sum rules. Default: CKM->False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real, 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.
amp2Quad (True|False): Indicates whether the symbol in amp2Type is quadratically (True) or linearly (False) dependent on A. Default: amp2Quad->False.

Returns:
system (Association): The inputted system association, modified to include formatted sum rule coefficient matrices and amplitude vector(s). New/modified keys and values:
- \"SR extract\" (List): Contains sum rule coefficient matrices at the selected b. If only ampType (amp2Type) is specified, contains only ASR (A2SR) coefficients. If both ampType and amp2Type are specified, contains A2SR coefficients. Initialized to None by generateSRs and redefined after running printSystem.
- \"Amp vector\" (List): Either is a vector of formatted A-type amplitudes (Symbols) or contains formatted vectors for a- and s-type amplitudes (List of Symbols) for the system. Initialized to None by generateSRs and redefined after running printSystem.";


Begin["`Private`"];


FlaSRHelp[function_Symbol]:=Module[{usage, details},
usage=Quiet[MessageName[function,"usage"]];
details=Quiet[MessageName[function,"details"]];

If[StringQ[usage],
(Print[usage];
If[StringQ[details],Print[details]];
),
Print["No documentation found for ",function,"."];
];
];


$FlaSRSession=.;


startPythonSession[session_,path_String]:=Module[{},
$FlaSRSession=If[MatchQ[session,_ExternalSessionObject],session,Message[startPythonSession::nosession];Return[$Failed]];
If[FileExistsQ[path],Null,Message[startPythonSession::nofile,path];Return[$Failed]];

ExternalEvaluate[$FlaSRSession,File[path]]
];

startPythonSession::nosession="No active Python session provided.";
startPythonSession::nofile="File `1` does not exist.";


(* Wrapper function for using ExternalEvaluate in the active Python session *)
pyEval[expr_,args_:<||>]:=ExternalEvaluate[$FlaSRSession,<|"Command"->expr,"Arguments"->args|>];


(* Extracts amplitudes from system (a Python System object) *)
Options[extractAmps]={partVal->{}};
extractAmps[system_,OptionsPattern[]]:=Module[{amplitudes,colNames,extractParticles,partVal=OptionValue[partVal]},
amplitudes=pyEval["System.extract_amps",system];
colNames={"Processes","QNs","n-tuples","Coords","Binary indices","q factor","mu","CG"};
amplitudes=Map[AssociationThread[colNames,#]&]@amplitudes;
amplitudes[[All,"Multiplet indices"]]=Map[{#[[1]],#[[3]]}&,amplitudes[[All,"Processes"]],{2}];

(* Processes from particle names *)
extractParticles[particles_,indices_]:=MapThread[MapThread[Part,{#1,#2}]&,{particles,indices}];
If[Length[partVal]>0,
(amplitudes[[All,"Processes"]]=Map[{extractParticles[partVal,#[[1]]],extractParticles[partVal,#[[2]]]}&,amplitudes[[All,"Processes"]]];
amplitudes[[All,"CKM"]]=Map[#[[2]]&,amplitudes[[All,"Processes"]],{2}];
amplitudes[[All,"Processes"]]=Map[{#[[1]],#[[3]]}&,amplitudes[[All,"Processes"]],{2}];
amplitudes[[All,"Processes"]]=Map[StringRiffle,amplitudes[[All,"Processes"]],{-2}];
amplitudes[[All,"Processes"]]=Map[{StringJoin[#[[1]]," \[Rule] ",#[[2]]]}&,amplitudes[[All,"Processes"]],{2}];
amplitudes[[All,"Processes"]]=Flatten/@amplitudes[[All,"Processes"]];
),
amplitudes=KeyDrop[#,"Processes"]&/@amplitudes;
];

amplitudes[[All,"QNs"]]=Map[If[#>0,"+"<>ToString[Rationalize[#],StandardForm],ToString[Rationalize[#],StandardForm]]&,amplitudes[[All,"QNs"]],{-1}];
amplitudes[[All,"QNs"]]=Map[StringRiffle,amplitudes[[All,"QNs"]],{-2}];
amplitudes[[All,"QNs"]]=Map[{StringJoin[#[[1]],ToString[Overscript[" \[Rule] ",#[[2]]],StandardForm],#[[3]]]}&,amplitudes[[All,"QNs"]],{2}];
amplitudes[[All,"QNs"]]=Flatten/@amplitudes[[All,"QNs"]];

amplitudes[[All,"mu"]]=Map[Sqrt[#[[1]]]*#[[2]]&,amplitudes[[All,"mu"]]]; (* mu factors *)
amplitudes[[All,"CG"]]=Map[ClebschGordan@@Rationalize[#]&,amplitudes[[All,"CG"]]]; (* CG coeffs from symmetrization *)
amplitudes
];


(* Constructs a Python System object and lists the multiplet contents for physical systems *)
defineSystem[in_,h_,out_,phys_]:=Module[{repSpins,reps,system,particles,n},
repSpins[list_]:=Table[(Length[list[[i]]]-1)/2,{i,Length[list]}];

If[phys,
(particles={in,h,out};
reps=Map[repSpins,particles];
),
(reps={in,h,out};
particles={};
)
];

n=Total[2*Flatten@reps];
If[EvenQ[n],Null,Message[defineSystem::argx,n];Return[$Failed]];

system=pyEval["define_system",{reps,phys}];
{system,particles}
];

defineSystem::argx="Invalid input. Number of would-be doublets is `1`. Please enter a system with an even number of doublets.";


Options[generateASRs]={phys->False};
generateASRs[in_,h_,out_,OptionsPattern[]]:=Module[{system,phys=OptionValue[phys],particles,aux,ASRs,amplitudes,dupPairs,factorsMat,simplifyFactors,irreps,n,p,nAmps,nASRs},
{system,particles}=defineSystem[in,h,out,phys]; (* constructs Python System object *)

aux=pyEval["System.extract_sys",system][[1]];

{system,ASRs,dupPairs}=pyEval["generate_srs",system]; (* M values for symmetrized system, contains duplicate amplitudes *)
amplitudes=extractAmps[system];

factorsMat=DiagonalMatrix[amplitudes[[All,"q factor"]]*amplitudes[[All,"mu"]]*amplitudes[[All,"CG"]]]; (* ASR matrices include q factors *)
ASRs=Map[# . factorsMat&,ASRs]; (* SRs for symmetrized system *)

(* Correct for duplicate amplitudes from symmetrization *)
Do[ASRs[[b]][[All,dupPairs[[All,1]]]]+=ASRs[[b]][[All,dupPairs[[All,2]]]],{b,Length[ASRs]}]; (* adds together cols of duplicate amp pairs *)
ASRs=Transpose[Delete[Transpose[#],List/@dupPairs[[All,2]]]]&/@ASRs; (* deletes duplicate cols *)
system=pyEval["System.remove_dups",{system,dupPairs}];
amplitudes=extractAmps[system,partVal->particles];

amplitudes=KeyDrop[#,"q factor"]&/@amplitudes;
simplifyFactors[mat_]:=(#/(If[Positive[DeleteCases[#,0][[1]]],1,-1]*Sqrt[Apply[GCD,DeleteCases[#,0]^2]]))&/@mat;
ASRs=If[aux,Map[RowReduce,ASRs],Map[simplifyFactors,ASRs]]; (* simplifies factors/row reduces (for systems requiring symmetrization) *)
ASRs=Map[Cases[Except@{0..}],ASRs];

{irreps,n,p}=pyEval["System.extract_sys",{system,True}][[2;;]];

system=<|"Amplitudes"->amplitudes,"ASRs"->ASRs|>; (* from here on, system is an association *)
nAmps=numAmps[system];
nASRs=numSRs[system];
system=<|"Irreps"->irreps,"Multiplets"->particles,"n doublets"->n,"p factor"->p,"n amps"->nAmps,"Amplitudes"->amplitudes,"n ASRs"->nASRs,"ASRs"->ASRs|>;

system
];


findA2SRMat[ASRMat_]:=Module[{ASRMatRR,pivotCols,freeCols,freeMat,groupedIndices,indices,colPairs,xMat,xMatNullSpace,groupedA2SR,A2SRMat},
ASRMatRR=RowReduce[ASRMat];

pivotCols=Map[Position[#,1][[1,1]]&,ASRMatRR];
freeCols=Complement[Range[Length[Transpose[ASRMatRR]]],pivotCols];
freeMat=ASRMatRR[[All,freeCols]]; (* extracts the free cols of an ASR matrix *)

groupedIndices=Join[pivotCols,freeCols];
indices=Table[Position[groupedIndices,i][[1,1]],{i,Length[groupedIndices]}]; (* value at the ith position gives the col of the grouped matrix corresponding to the original ith col *)

If[Dimensions[freeMat][[2]]>=2,
(* case with cross terms: find null space of cross terms matrix *)
(colPairs=Subsets[Transpose[freeMat],{2}]; (* list of all unique column pairs *)
xMat=Transpose@Map[Times@@#&,colPairs]; (* cross terms matrix: cols are pairwise multiples of free cols *)
xMatNullSpace=NullSpace[Transpose@xMat]; (* xMatNullSpace.xMat = 0; i.e., left mult of cross terms mat by null space mat sets cross terms to 0 *)

groupedA2SR=If[xMatNullSpace=={},{},Join[xMatNullSpace,-xMatNullSpace . (freeMat^2),2]]; (* multiplies free cols by xMatNullSpace and brings back to LHS *)
),
(* case with no cross terms: A^2SR has the trivial form s^2-s^2 = 0 *)
(xMat={{0}};
groupedA2SR=Join[IdentityMatrix[Length[ASRMatRR]],-(freeMat^2),2];
)
];

A2SRMat=If[groupedA2SR=={},{},Transpose@Table[Transpose[groupedA2SR][[indices[[i]]]],{i,Length[indices]}]]; (* rearranges A^2SR mat into original order *)

A2SRMat
];


Options[generateSRs]={phys->False,obs->"Diff"};
generateSRs[in_,h_,out_,opts:OptionsPattern[]]:=Module[{system,phys=OptionValue[phys],obs=OptionValue[obs],ASRs,A2SRs,integrateA2SRs,nA2SRs},
system=generateASRs[in,h,out,Sequence@@FilterRules[{opts},Options[generateASRs]]];

ASRs=system[["ASRs"]];
A2SRs=Map[findA2SRMat,ASRs];

(* diff and int observables *)
integrateA2SRs[]:=Module[{inMulti,outMulti,ampIndices,uniqueInMulti,uniqueOutMulti,canonInMulti,canonOutMulti,formIndexPairs,ampIndexPairs,uniqueAmps,convertToAmpPairs,negCols,identicalColGroups,negMatCols,integrateA2SRMat},
inMulti=system[["Multiplets"]][[1]];
outMulti=system[["Multiplets"]][[3]];
ampIndices=Flatten[system[["Amplitudes"]][[All,"Multiplet indices"]],1]; (* list out all individual amplitudes {{m in},{m out}} *)

uniqueInMulti=PositionIndex[inMulti]; (* assoc of all unique in state multiplets, {multi} -> {indices in amp table} *)
uniqueOutMulti=PositionIndex[outMulti];

canonInMulti=Map[Position[Keys[uniqueInMulti],#][[1,1]]&,inMulti]; (* canonical indices of in state multiplets in unique list *)
canonOutMulti=Map[Position[Keys[uniqueOutMulti],#][[1,1]]&,outMulti];

formIndexPairs[{ampIn_,ampOut_}]:={MapIndexed[{canonInMulti[[#2]][[1]],#1}&,ampIn],MapIndexed[{canonOutMulti[[#2]][[1]],#1}&,ampOut]};
ampIndexPairs=formIndexPairs/@ampIndices; (* amp indices list becomes {multiplet index, component} list *)

uniqueAmps=PositionIndex[Map[({#[[1]],Sort[#[[2]]]})&,ampIndexPairs]]; (* assoc of unique amps -> col indices in amp list. initial states must be exact matches, final states up to permutations *)
convertToAmpPairs[uniqueAmps_]:=Module[{ampAssoc=uniqueAmps,ampPairIndices,colList},
ampAssoc=Select[ampAssoc,OddQ[First[#]]&]; (* keep only cases with lower amps as first elements *)
ampPairIndices[indices_]:=Module[{pairs=Ceiling[indices/2],keepQ},
keepQ=MapIndexed[
Function[{i,idx},
OddQ[i]||idx[[1]]==1||pairs[[idx[[1]]]]!=pairs[[idx[[1]]-1]] (* keep lower amps, first elem, conj amps of unique amp pairs *)
],
indices
];
{Pick[pairs,keepQ],Pick[pairs,MapThread[EvenQ[#1]&&#2&,{indices,keepQ}]]} (* {unique amp pairs, neg cols} *)
];
ampAssoc=ampPairIndices/@ampAssoc;
colList=Flatten@Values[ampAssoc][[All,2]];
ampAssoc=#[[1]]&/@ampAssoc;
{ampAssoc,colList}
];
{uniqueAmps,negCols}=convertToAmpPairs[uniqueAmps]; (* uniqueAmps indices now refer to pairs, negCols indicate which cols need to be negated for delta matrices *)

identicalColGroups=Select[Values[uniqueAmps],Length[#]>1&]; (* list of groups of identical cols, cols refer to amp pairs now *)

uniqueAmps=#[[1]]&/@uniqueAmps;
uniqueAmps=AssociationMap[Reverse,uniqueAmps];
uniqueAmps=(Times@@Factorial[Values[Counts[#[[2]]]]])&/@uniqueAmps; (* assoc with unique col -> k! *)

negMatCols[A2SRMat_]:=Module[{negMatCol,newA2SRMat},
negMatCol[matTranspose_,col_]:=Module[{m=matTranspose},m[[col]]=-1*m[[col]];m];
newA2SRMat=If[A2SRMat=={},{},Fold[negMatCol,Transpose@A2SRMat,negCols]]; (* multiply neg cols by -1 *)

Transpose@newA2SRMat
];

integrateA2SRMat[A2SRMat_]:=Module[{addMatCols,scaleMatCol,newA2SRMat},
addMatCols[mat_,cols_]:=Module[{m=mat},m[[All,cols[[1]]]]=Total[m[[All,cols]],{2}];m];
scaleMatCol[matTranspose_,col_]:=Module[{m=matTranspose},m[[col]]=uniqueAmps[col]*m[[col]];m];

If[A2SRMat=={},
newA2SRMat={},
(newA2SRMat=A2SRMat;
newA2SRMat=Fold[addMatCols,newA2SRMat,identicalColGroups]; (* combine identical cols *)
newA2SRMat=Transpose@Fold[scaleMatCol,Transpose@newA2SRMat,Keys[uniqueAmps]]; (* multiply unique cols by k! *)
newA2SRMat=newA2SRMat[[All,Keys[uniqueAmps]]]; (* keep only unique cols *)
)
];

newA2SRMat
];

A2SRs=MapAt[negMatCols,A2SRs,List/@Range[1,Length[A2SRs],2]]; (* multiply neg cols in delta matrices by -1 *)
A2SRs=integrateA2SRMat/@A2SRs;
A2SRs=Map[If[#=={},{},Cases[RowReduce@#,Except@{0..}]]&,A2SRs];
system[["Unique amp pairs"]]=Keys[uniqueAmps]; (* add unique amp pairs key to system assoc *)
];

If[phys&&(obs==="Int"),integrateA2SRs[],Null];
system[["Amplitudes"]]=KeyDrop["Multiplet indices"]/@system[["Amplitudes"]]; (* removes this key only when using generateSRs, but not generateASRs *)
nA2SRs=Table[Length[A2SRs[[i]]],{i,Length[A2SRs]}];

AssociateTo[system,<|"n A2SRs"->nA2SRs,"A2SRs"->A2SRs,"Amp vector"->None,"SR extract"->None|>]
];


(* Returns the total number of amplitudes in the system *)
numAmps[system_,nPairs_:False]:=Module[{amplitudes=system[["Amplitudes"]]},
If[!nPairs,
Length@DeleteDuplicates@Flatten@amplitudes[[All,"Binary indices"]],
Length@Flatten@amplitudes[[All,"Coords"]]
]
];


defaultAmpKeys={"Processes","QNs","n-tuples","Coords","Binary indices","mu","CG","CKM","Multiplet indices"};


(* Adds a column to amplitudes *)
SetAttributes[labelAmps,HoldFirst];
Options[labelAmps]={labeling->"Amplitudes"};
labelAmps[system_,colName_String,labels_List,OptionsPattern[]]:=Module[{sysVal=Evaluate[system],amplitudes,pairs,nAmps,labelVals,labeling=OptionValue[labeling],indices,labelIndices},
If[MemberQ[defaultAmpKeys,colName],Message[labelAmps::argval,colName];Return[$Failed]];

amplitudes=sysVal[["Amplitudes"]];
pairs=If[labeling=="Amplitudes",False,True];
nAmps=numAmps[sysVal,pairs];
If[nAmps==Length@Flatten[labels],Null,Message[labelAmps::arglen,labels,nAmps,Length@Flatten[labels]];Return[$Failed]];

labelVals=Which[
labeling=="Amplitude pairs",
labels,
labeling=="Amplitudes",
(indices=amplitudes[[All,"Binary indices"]];
labelIndices=Sort[Join[Range[nAmps],Table[If[indices[[i,1]]==indices[[i,2]],2*i-1,Nothing],{i,Length[indices]}]]];
Partition[Table[Flatten[labels][[i]],{i,labelIndices}],2]
),
True,Message[labelAmps::badmode,labeling];Return[$Failed]
];

amplitudes=MapThread[Prepend,{amplitudes,colName->#&/@labelVals}];
system[["Amplitudes"]]=amplitudes;
amplitudes
];

labelAmps::argval="Invalid column name `1`. Use a column name that does not conflict with the built-in amplitude association keys.";
labelAmps::arglen="The labels list `1` has an incorrect number of labels. Expected `2` labels, got `3`. Check that the labeling option is correct.";
labelAmps::badmode="Unknown labeling mode `1`. Use \"Amplitudes\" or \"Amplitude pairs\".";


(* Removes columns from amplitudes *)
SetAttributes[unlabelAmps,HoldFirst];
unlabelAmps[system_,colNames_]:=Module[{sysVal=Evaluate[system],amplitudes},
If[ContainsAny[defaultAmpKeys,Flatten@List@colNames],Message[unlabelAmps::argval,colNames];Return[$Failed]];

amplitudes=sysVal[["Amplitudes"]];
amplitudes=KeyDrop[colNames]/@amplitudes;
system[["Amplitudes"]]=amplitudes;
amplitudes
];

unlabelAmps::argval="Invalid column name(s) `1`. Cannot remove any of the built-in amplitude association keys.";


(* Prints a table of amplitudes *)
Options[printAmps]={showFactors->False};
printAmps[system_,OptionsPattern[]]:=Module[{amplitudes=system[["Amplitudes"]],indices,selfConj,nAmps,showFactors=OptionValue[showFactors],signs},
(* Delete self-conjugate duplicate values from display *)
indices=amplitudes[[All,"Binary indices"]];
selfConj=Table[If[indices[[i,1]]==indices[[i,2]],i,Nothing],{i,Length[indices]}];
amplitudes[[selfConj]]=Map[If[ListQ[#],#[[1]],#]&,amplitudes[[selfConj]],{2}];

nAmps=numAmps[system];
If[!showFactors,amplitudes=KeyDrop[#,{"Binary indices","mu","CG","Multiplet indices"}]&/@amplitudes,Null]; (* show/hide internal factors from display *)
signs=If[system[["p factor"]]==1,{"-","+"},{"+","-"}];

Print["Amplitude table","\n",
"Number of amplitudes: ",nAmps,"\n",
"a/s definitions: \!\(\*SubscriptBox[\(a\), \(i\)]\) = \!\(\*SubscriptBox[\(A\), \(i\)]\) ",signs[[1]]," \!\(\*SubscriptBox[OverscriptBox[\(A\), \(_\)], \(i\)]\), \!\(\*SubscriptBox[\(s\), \(i\)]\) = \!\(\*SubscriptBox[\(A\), \(i\)]\) ",signs[[2]]," \!\(\*SubscriptBox[OverscriptBox[\(A\), \(_\)], \(i\)]\)","\n",
"\[CapitalDelta]/\[CapitalSigma] definitions: \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(i\)]\) = |\!\(\*SubscriptBox[\(A\), \(i\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\) - |\!\(\*SubscriptBox[OverscriptBox[\(A\), \(_\)], \(i\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\), \!\(\*SubscriptBox[\(\[CapitalSigma]\), \(i\)]\) = |\!\(\*SubscriptBox[\(A\), \(i\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\) + |\!\(\*SubscriptBox[OverscriptBox[\(A\), \(_\)], \(i\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)","\n",
TableForm[Prepend[Values/@amplitudes,Keys[amplitudes[[1]]]]]
];
];


(* Returns a list of the orders of breaking *)
listbOrders[b_,nOrders_]:=Module[{bVal=b,bList},
bVal=b/.{
{i_,j_,k_:1}:>Span[
Replace[i,{Null->1,x_Integer:>x+1}],
Replace[j,{Null->All,x_Integer:>x+1}],
k
],
{i_}:>Span[Replace[i,{Null->1,x_Integer:>x+1}],All],
i_Integer:>Span[i+1,i+1],
All:>All
};
bList=Flatten@{(Range[nOrders]-1)[[bVal]]};
bList
];


(* Returns the number of amplitude or amplitude-squared sum rules at each order of breaking *)
Options[numSRs]={b->All};
numSRs[system_,squared_:False,OptionsPattern[]]:=Module[{b=OptionValue[b],SRs,indices,sublist},
SRs=If[squared,
If[KeyExistsQ[system,"A2SRs"],system[["A2SRs"]],(Message[numSRs::missingkey];Return[$Failed])],
system[["ASRs"]]
];
If[b===All,
Table[Length[SRs[[i]]],{i,Length[SRs]}],
(indices=listbOrders[b,Length[SRs]]+1;
sublist=SRs[[indices]];
Table[Length[sublist[[i]]],{i,Length[sublist]}]
)
]
];

numSRs::missingkey="This system does not contain the A2SRs key.";


(* Prints sum rules at each order of breaking *)
SetAttributes[printSRs,HoldFirst];
Options[printSRs]={ampType->None,amp2Type->None,ampFormat->None,showSRs->True,expandSRs->False,CKM->False,b->All,amp2Quad->False};
printSRs[system_,opts:OptionsPattern[]]:=Module[{sysVal=Evaluate[system],ampType=OptionValue[ampType],amp2Type=OptionValue[amp2Type],ampFormat=OptionValue[ampFormat],showSRs=OptionValue[showSRs],expandSRs=OptionValue[expandSRs],CKM=OptionValue[CKM],b=OptionValue[b],amp2Quad=OptionValue[amp2Quad],formatSRMats,ampsToVectors,printWrittenSRs,syms,squared,p,SRs,amplitudes,indices,selfConj,physAmps,ampVectors,bList},
(* Function definitions *)
(* Double widths of SRs matrices, add CKM and p factors *)
formatSRMats[amps_]:=Module[{factorsMats,delSelfConj},
factorsMats=If[squared,
{DiagonalMatrix@Flatten@Table[{1,-1},Length[amps]],IdentityMatrix[2*Length[amps]]}, (* {\[CapitalDelta] mat,\[CapitalSigma] mat} *)
{DiagonalMatrix@Flatten@Table[{1,-p},Length[amps]],DiagonalMatrix@Flatten@Table[{1,p},Length[amps]]} (* {a mat,s mat} *)
];
SRs=Map[Transpose@Flatten[{#,#}&/@Transpose@#,1]&,SRs]; (* duplicates cols of SRs matrix *)
SRs=MapIndexed[If[#1=={},{},If[OddQ[First[#2]],#1 . factorsMats[[1]],#1 . factorsMats[[2]]]]&,SRs]; (* adds factors to SRs matrices *)

delSelfConj[mat_]:=Module[{matVal=mat},
matVal=If[matVal=={},{},
(matVal[[All,(2*selfConj[[1]]-1)]]+=matVal[[All,2*selfConj[[1]]]]; (* adds together self-conj cols *)
Transpose[Delete[Transpose[matVal],2*selfConj[[1]]]]) (* deletes extra self-conj col *)
]
];

SRs=If[Length[selfConj]>0,Map[delSelfConj[#]&,SRs],SRs]; (* corrects for self-conj amps if necessary *)
SRs
];

(* Return a list of two amplitude vectors (either a,s or two A) according to formatting options *)
ampsToVectors[amps_]:=Module[{ampsToVector,vec1,vecList},
If[ampFormat=="Processes"&&!KeyExistsQ[amps[[1]],"Processes"],Message[printSRs::invalidformat];Return[$Failed],Null];

ampsToVector[ampSym_Symbol]:=Module[{vector,rule},
vector=Map[ampSym[#]&,
If[!physAmps,
(Switch[ampFormat, (* a/s amps *)
"Processes",(Message[printSRs::invalidformat];Return[$Failed]),
"QNs",amps[[All,"QNs",1]], (* a/s-type amps with QNs *)
"n-tuples",amps[[All,"n-tuples",1]], (* a/s-type amps with n-tuples *)
"Coords",amps[[All,"Coords"]], (* a/s-type amps with coords *)
"Binary indices",amps[[All,"Binary indices",1]], (* a/s-type amps with numbered subscripts *)
_String/;KeyExistsQ[amps[[1]],ampFormat],amps[[All,ampFormat]], (* custom amplitude format *)
_,(Message[printSRs::invalidformat];Return[$Failed])
]
),
(Switch[ampFormat, (* A amps *)
"Processes",Flatten@amps[[All,"Processes"]], (* A amps with physical processes *)
"QNs",Flatten@amps[[All,"QNs"]], (* A amps with QNs *)
"n-tuples",Flatten@amps[[All,"n-tuples"]], (* A amps with n-tuples *)
"Coords",(Message[printSRs::invalidformat];Return[$Failed]),
"Binary indices",Flatten@amps[[All,"Binary indices"]], (* A amps with numbered subscripts *)
_String/;KeyExistsQ[amps[[1]],ampFormat],Flatten@amps[[All,ampFormat]], (* custom amplitude format *)
_,(Message[printSRs::invalidformat];Return[$Failed])
]
)
]
];

If[physAmps&&CKM,vector=vector/Flatten@amps[[All,"CKM"]],Null]; (* divide by CKM *)

rule=If[ampFormat=="Processes"||ampFormat=="QNs",
expr_Symbol[i_]/;expr=!=List:>StringJoin[ToString[expr],"(",i,")"],
expr_Symbol[i_]/;expr=!=List:>Subscript[expr,i]
];

vector=vector/.rule; (* print formatting of amplitudes *)
If[physAmps&&squared&&(!amp2Quad),vector=Abs[#]^2&/@vector,Null]; (* square the vector *)

vector
];

vec1=ampsToVector[syms[[1]]];
vecList=If[!physAmps,{vec1,ampsToVector[syms[[2]]]},{vec1,vec1}];
vecList=If[Length[selfConj]>0,
(If[physAmps,
Map[Delete[#,2*selfConj[[1]]]&,vecList],
Delete[vecList,{If[EvenQ[p],1,2],selfConj[[1]]}]
]
),
vecList
]; (* for A/A^2, delete conj of self-conj amp; for a/s/\[CapitalDelta]/\[CapitalSigma], delete 0 cols *)
vecList
];

printWrittenSRs[]:=Module[{numSRsList,matForm,writtenSRs},
numSRsList=numSRs[sysVal,squared,Sequence@@FilterRules[{opts},Options[numSRs]]];
matForm:=If[expandSRs,MatrixForm[#]&,MatrixForm[#[[2;;]],TableHeadings->{None,#[[1]]}]&];
writtenSRs=If[expandSRs,
MapIndexed[If[#1=={},{},If[OddQ[First[#2]],#1 . ampVectors[[1]],#1 . ampVectors[[2]]]]&,SRs],
MapIndexed[If[#1=={},{},If[OddQ[First[#2]],Prepend[#1,ampVectors[[1]]],Prepend[#1,ampVectors[[2]]]]]&,SRs]
]; (* combine amp vecs with SRs matrices *)

If[squared,Print["Amplitude-squared sum rules"],Print["Amplitude sum rules"]];
MapIndexed[
Print["b = ",bList[[#2[[1]]]],"\n",
"Number of SRs: ",numSRsList[[#2[[1]]]],"\n",
If[#1=={},"No sum rules at this order.",matForm[#1]//TraditionalForm]
]&,
writtenSRs[[bList+1]]
];
];



(* Function calls *)
If[(ampType==None)&&(amp2Type==None),(Message[printSRs::missingtype];Return[$Failed]),None]; (* error if no amp type is specified *)
{syms,squared}=If[(amp2Type=!=None),{amp2Type,True},{ampType,False}]; (* only one opt should be specified, but if both are specified, keep amp2Type by default *)
p=sysVal[["p factor"]];
SRs=If[squared,
If[KeyExistsQ[sysVal,"A2SRs"],sysVal[["A2SRs"]],(Message[printSRs::missingkey];Return[$Failed])],
sysVal[["ASRs"]]
]; (* set SRs list to either ASRs or A2SRs *)
amplitudes=sysVal[["Amplitudes"]];
If[(KeyExistsQ[sysVal,"Unique amp pairs"])&&squared,amplitudes=amplitudes[[sysVal[["Unique amp pairs"]]]],Null];(* if obs -> True and printing A2SRs, reduce amps assoc to unique amps *)
indices=amplitudes[[All,"Binary indices"]];
selfConj=Table[If[indices[[i,1]]==indices[[i,2]],i,Nothing],{i,Length[indices]}];

physAmps=Switch[Length[syms],
1,True, (* A/A^2 *)
2,False, (* a/s/\[CapitalDelta]/\[CapitalSigma] *)
_,(Message[printSRs::invalidformat];Return[$Failed])
]; (* indicate whether amplitudes are phys (True) or group-theoretic (False) *)
If[ampFormat==None,
If[physAmps&&(KeyExistsQ[amplitudes[[1]],"Processes"]),ampFormat="Processes",ampFormat="n-tuples"],
Null
]; (* if no ampFormat was specified, set to n-tuples by default unless it's a physical system with A amplitudes *)
If[physAmps,
formatSRMats[amplitudes],
If[Length[selfConj]>0,
SRs=MapIndexed[If[#1=={},{},If[EvenQ[p]==OddQ[First[#2]],Transpose[Delete[Transpose[#1],selfConj[[1]]]],#1]]&,SRs],
Null
]
]; (* for phys amps, double widths of SRs matrices, add CKM and p factors. for group-theoretic, correct for self-conj amps *)
ampVectors=ampsToVectors[amplitudes]; (* return list of two vectors of amplitudes; the two vecs are the same in phys case *)

bList=listbOrders[b,Length[SRs]];
If[showSRs,printWrittenSRs[],Null]; (* print SRs *)

system[["Amp vector"]]=If[physAmps,ampVectors[[1]],ampVectors];
SRs=If[Length[bList]==1,SRs[[bList+1]][[1]],SRs[[bList+1]]];
system[["SR extract"]]=SRs;
SRs
];

printSRs::missingtype="Missing amplitude type. Please specify an amplitude type for either ampType or amp2Type.";
printSRs::invalidformat="Invalid or missing amplitude format.";
printSRs::missingkey="This system does not contain the A2SRs key.";


SetAttributes[printSystem,HoldFirst];
Options[printSystem]=Join[{showReps->True,showAmps->True,showASRs->True,showA2SRs->True},Options[printAmps],Options[printSRs]];
printSystem[system_,opts:OptionsPattern[]]:=Module[{sysVal=Evaluate[system],showReps=OptionValue[showReps],showAmps=OptionValue[showAmps],showASRs=OptionValue[showASRs],showA2SRs=OptionValue[showA2SRs],printLine,irreps},
printLine[]:=Print["-------------------------"];
irreps=sysVal[["Irreps"]];

If[showReps,
(Print["System: ",Sort@Flatten@irreps];
printLine[];
Print["Number of would-be doublets: ",sysVal[["n doublets"]],"\n",
"In: ",irreps[[1]],"\n",
"H: ",irreps[[2]],"\n",
"Out: ",irreps[[3]]
];
),
Null
];

If[showAmps,
If[showReps,printLine[]];
printAmps[sysVal,Sequence@@FilterRules[{opts},Options[printAmps]]];
];

If[showASRs,
If[showReps||showAmps,printLine[]];
printSRs[sysVal,Sequence@@FilterRules[FilterRules[{opts},Except[amp2Type]],Options[printSRs]],showSRs->showASRs,amp2Type->None];
];

If[showA2SRs,
If[showReps||showAmps||showASRs,printLine[]];
printSRs[sysVal,Sequence@@FilterRules[FilterRules[{opts},Except[ampType]],Options[printSRs]],showSRs->showA2SRs,ampType->None];
];

system=sysVal;
sysVal
];


End[]


EndPackage[]
