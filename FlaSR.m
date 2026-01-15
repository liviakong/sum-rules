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
path (String): Path to the sum_rules.py Python file.

Returns:
An ExternalFunction object indicating the sum_rules.py file has been loaded into the Python session."

generateASRs::usage="generateASRs[in,h,out] finds all amplitudes and amplitude sum rules (ASRs) for a given system.";
generateASRs::details=
"Arguments:
in (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the incoming state.
h (List): Contains U-spins (Real) OR CKM factors (List of Symbols) in the Hamiltonian.
out (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the outgoing state.

Options:
phys (True|False): Indicates whether function arguments contain U-spins (False) or particle multiplets/CKM factors (True). Default: phys->False.

Returns:
system (Association): All information about the system's representations, amplitudes, and amplitude sum rules. Keys and values:
- \"Irreps\" (List): Inputted U-spin representations (List of Reals) in {incoming reps, H reps, outgoing reps} format.
- \"n doublets\" (Real): Number of would-be doublets.
- \"p factor\" (Real): (-1)^p parity factor for the system determining forms of a/s-type amplitudes.
- \"n amps\" (Real): Number of amplitudes in the system.
- \"Amplitudes\" (Association): Contains all amplitudes in the system. Keys and values:
	- \"Processes\" (List): Contains physical processes constructed from particle names (String) for a process and its U-spin conjugate process. Only available for physical systems.
	- \"QNs\" (List): Contains processes written using m quantum numbers (String), where m is the third component of U-spin, for a process and its U-spin conjugate process.
	- \"n-tuples\" (List): Contains representations (String) of an amplitude and its U-spin conjugate amplitude using comma-separated tuples of substrings. Each substring is comprised of '-'s and '+'s and encodes u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components.
	- \"Nodes\" (String): Representation of an amplitude pair using the coordinate notation of the lattice used to derive sum rules.
	- \"Binary indices\" (List): Contains indices (Real), written in base 10, derived from converting the n-tuple and its U-spin conjugate into binary numbers through '-' <-> 0 and '+' <-> 1.
	- \"mu\" (Real): mu-factor for the node in the lattice used to derive sum rules.
	- \"CG\" (Real): Clebsch-Gordan coefficient from symmetrization for systems without doublets. Equal to 1 for all amplitudes for a system with at least one doublet.
	- \"CKM\" (List): Contains weak interaction factors (Real) from the Hamiltonian. Only appears for physical systems.
- \"n ASRs\" (List): Contains the number of amplitude sum rules (Real) at each order of breaking.
- \"ASRs\" (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.";

numAmps::usage="numAmps[system] returns the total number of amplitudes in the system.";
numAmps::details=
"Arguments:
system (Association): A system association. See the documentation for generateASRs for details.

Options:
nPairs (True|False): Indicates whether to return the number of amplitudes (False) or amplitude pairs (True). Default: nPairs: False.

Returns:
The total number of amplitudes or amplitude pairs in the system (Real).";

labelAmps::usage="labelAmps[system,colName,labels] modifies system to add a column of user-defined labels to system[[\"Amplitudes\"]]."
labelAmps::details=
"Arguments:
system (Association): A system association. See the documentation for generateASRs for details.
colName (String): Name of new column to add to amplitudes association.
labels (List): Contains labels (Any) for each amplitude or amplitude pair. Number of labels must equal number of amplitudes or amplitude pairs.

Options:
labeling (String): Labeling mode to indicate whether user is labeling single amplitudes (\"Amplitudes\") or amplitude pairs (\"Amplitude pairs\"). Default: labeling -> \"Amplitudes\".

Returns:
system[[\"Amplitudes\"]] (Association): The modified amplitudes association which includes a new column of user-defined labels. New/modified keys and values:
- colName (List|Any): Either contains user-defined labels (Any) for an amplitude and its U-spin conjugate or a single label (Any) for an amplitude pair.";

unlabelAmps::usage="unlabelAmps[system,colNames] modifies system to remove columns from system[[\"Amplitudes\"]]."
unlabelAmps::details=
"Arguments:
system (Association): A system association. See the documentation for generateASRs for details.
colNames (String|List): Name(s) of column(s) to remove from amplitudes association.

Returns:
system[[\"Amplitudes\"]] (Association): The modified amplitudes association from which the specified columns have been removed.";

printAmps::usage="printAmps[system] prints a table of amplitudes.";
printAmps::details=
"Arguments:
system (Association): A system association. See the documentation for generateASRs for details.

Options:
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.";

numSRs::usage="numSRs[system] returns the number of amplitude or squared amplitude sum rules at each order of breaking.";
numSRs::details=
"Arguments:
system (Association): A system association. See the documentation for generateASRs for details.

Options:
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.
squared (True|False): Indicates whether to return number of sum rules of amplitudes (False) or squared amplitudes (True). Default: squared->False.

Returns:
Number of (squared) amplitude sum rules at each order of breaking (List of Reals).";

printSRs::usage="printSRs[system] prints amplitude or squared amplitude sum rules at each order of breaking.";
printSRs::details=
"Arguments:
system (Association): A system association. See the documentation for generateASRs for details.

Options:
showSRs (True|False): Default: Indicates whether to print sum rules. showSRs->True.
expandSRs (True|False): Indicates whether to write each row of a sum rules matrix as an algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandSRs->True.
ampFormat (String): Specified format for displaying amplitudes. Options are a/s-type amplitudes with n-tuples (\"a/s n-tuples\"), a/s-type amplitudes with numbered indices (\"a/s indices\"), a/s-type amplitudes with nodes (\"a/s nodes\"), A amplitudes with numbered indices (\"A indices\"), and A amplitudes with physical process (\"A physical\"). User-defined labels for a column of amplitudes (or amplitude pairs) can also be used using the syntax \"A col\" (or \"a/s col\" for pairs), where \"col\" is the name of the column containing custom labels. Default: ampFormat->\"a/s n-tuples\".
CKM (True|False): Indicates whether to include CKM factors in the sum rules (True) or not (False). Default: CKM->False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.
squared (True|False): Indicates whether to print amplitude sum rules (False) or squared amplitude sum rules (True). Default: squared->False.

Returns:
system[[\"Formatted ASRs\"]] (or system[[\"Formatted A2SRs\"]]) (List): A new key in the system association containing Lists of (squared) amplitude sum rules/matrices of (squared) amplitude sum rule coefficients formatted according to the sum rule printing settings.";

printSystem::usage="printSystem[system] prints information about the system's representations, amplitudes, and sum rules and modifies the system to include formatted sum rules.";
printSystem::details=
"Arguments:
system (Association): A system association. See the documentation for generateASRs for details.

Options:
showReps (True|False): Indicates whether to print information about the system's representations. Default: showReps->True.
-----
showAmps (True|False): Indicates whether to print information about the system's amplitudes. Default: showAmps->True.
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.
-----
showASRs (True|False): Default: showASRs->True. Note: ASRs and A2SRs use the same formatting options.
showA2SRs (True|False): Default: showA2SRs->False. Note: ASRs and A2SRs use the same formatting options.
expandSRs (True|False): Indicates whether to write each row of a sum rules matrix as an algebraic expression of (squared) amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandSRs->True.
ampFormat (String): Specified format for displaying (squared) amplitudes. Options are a/s-type amplitudes with n-tuples (\"a/s n-tuples\"), a/s-type amplitudes with numbered indices (\"a/s indices\"), a/s-type amplitudes with nodes (\"a/s nodes\"), A amplitudes with numbered indices (\"A indices\"), and A amplitudes with physical process (\"A physical\"). Default: ampFormat->\"a/s n-tuples\".
CKM (True|False): Indicates whether to include CKM factors in the sum rules (True) or not (False). Default: CKM->False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
system (Association): The inputted system association, modified to include formatted amplitude sum rules. New/modified keys and values:
- \"Formatted ASRs\" (or \"Formatted A2SRs\") (List): Contains Lists of (squared) amplitude sum rules/matrices of (squared) amplitude sum rule coefficients formatted according to the sum rule printing settings.";

generateA2SRs::usage="generateA2SRs[in,h,out] finds all amplitudes, amplitude sum rules (ASRs), and squared amplitude sum rules (A2SRs) for a given system.";
generateA2SRs::details=
"Arguments:
in (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the incoming state.
h (List): Contains U-spins (Real) OR CKM factors (List of Symbols) in the Hamiltonian.
out (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the outgoing state.

Options:
phys (True|False): Indicates whether function arguments contain U-spins (False) or particle multiplets/CKM factors (True). Default: phys->False.

Returns:
system (Association): All information about the system's representations, amplitudes, and amplitude sum rules. Keys and values:
- \"Irreps\" (List): Inputted U-spin representations (List of Reals) in {incoming reps, H reps, outgoing reps} format.
- \"n doublets\" (Real): Number of would-be doublets.
- \"p factor\" (Real): (-1)^p parity factor for the system determining forms of a/s-type amplitudes.
- \"n amps\" (Real): Number of amplitudes in the system.
- \"Amplitudes\" (Association): Contains all amplitudes in the system. Keys and values:
	- \"Processes\" (List): Contains physical processes constructed from particle names (String) for a process and its U-spin conjugate process. Only available for physical systems.
	- \"QNs\" (List): Contains processes written using m quantum numbers (String), where m is the third component of U-spin, for a process and its U-spin conjugate process.
	- \"n-tuples\" (List): Contains representations (String) of an amplitude and its U-spin conjugate amplitude using comma-separated tuples of substrings. Each substring is comprised of '-'s and '+'s and encodes u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components.
	- \"Nodes\" (String): Representation of an amplitude pair using the coordinate notation of the lattice used to derive sum rules.
	- \"Binary indices\" (List): Contains indices (Real), written in base 10, derived from converting the n-tuple and its U-spin conjugate into binary numbers through '-' <-> 0 and '+' <-> 1.
	- \"mu\" (Real): mu-factor for the node in the lattice used to derive sum rules.
	- \"CG\" (Real): Clebsch-Gordan coefficient from symmetrization for systems without doublets. Equal to 1 for all amplitudes for a system with at least one doublet.
	- \"CKM\" (List): Contains weak interaction factors (Real) from the Hamiltonian. Only appears for physical systems.
- \"n ASRs\" (List): Contains the number of amplitude sum rules (Real) at each order of breaking.
- \"ASRs\" (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.
- \"n A2SRs\" (List): Contains the number of squared amplitude sum rules (Real) at each order of breaking.
- \"A2SRs\" (List): Contains matrices of squared amplitude sum rule coefficients (Real) corresponding to each order of breaking.
- \"n thetas\" (List): Contains the number of columns (Real) in the cross terms matrix at each order of breaking.
- \"Theta ranks\" (List): Contains the rank of the cross terms matrix (Real) at each order of breaking.";


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


pyEval[expr_,args_:<||>]:=ExternalEvaluate[$FlaSRSession,<|"Command"->expr,"Arguments"->args|>];


(* Extracts amplitudes from system (a Python System object) *)
Options[extractAmps]={partVal->{}};
extractAmps[system_,OptionsPattern[]]:=Module[{amplitudes,colNames,extractParticles,partVal=OptionValue[partVal]},
amplitudes=pyEval["System.extract_amps",system];
colNames={"Processes","QNs","n-tuples","Nodes","Binary indices","q factor","mu","CG"};
amplitudes=Map[AssociationThread[colNames,#]&]@amplitudes;

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
system=<|"Irreps"->irreps,"n doublets"->n,"p factor"->p,"n amps"->nAmps,"Amplitudes"->amplitudes,"n ASRs"->nASRs,"ASRs"->ASRs|>;

system
];


(* Returns the total number of amplitudes in the system *)
numAmps[system_,nPairs_:False]:=Module[{amplitudes=system[["Amplitudes"]]},
If[!nPairs,
Length@DeleteDuplicates@Flatten@amplitudes[[All,"Binary indices"]],
Length@Flatten@amplitudes[[All,"Nodes"]]
]
];


(* Adds a column to amplitudes *)
SetAttributes[labelAmps,HoldFirst];
Options[labelAmps]={labeling->"Amplitudes"};
labelAmps[system_,colName_String,labels_List,OptionsPattern[]]:=Module[{sysVal=Evaluate[system],amplitudes,pairs,nAmps,labelVals,labeling=OptionValue[labeling],indices,labelIndices},
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

labelAmps::arglen="The argument supplied to `1` has an incorrect number of labels. Expected `2` labels, got `3`. Check that the labeling option is correct.";
labelAmps::badmode="Unknown labeling mode `1`. Use \"Amplitudes\" or \"Amplitude pairs\".";


(* Removes columns from amplitudes *)
SetAttributes[unlabelAmps,HoldFirst];
unlabelAmps[system_,colNames_]:=Module[{sysVal=Evaluate[system],amplitudes},
amplitudes=sysVal[["Amplitudes"]];
amplitudes=KeyDrop[colNames]/@amplitudes;
system[["Amplitudes"]]=amplitudes;
amplitudes
];


(* Prints a table of amplitudes *)
Options[printAmps]={showFactors->False};
printAmps[system_,OptionsPattern[]]:=Module[{amplitudes=system[["Amplitudes"]],indices,selfConj,nAmps,showFactors=OptionValue[showFactors],signs},
(* Delete self-conjugate duplicate values from display *)
indices=amplitudes[[All,"Binary indices"]];
selfConj=Table[If[indices[[i,1]]==indices[[i,2]],i,Nothing],{i,Length[indices]}];
amplitudes[[selfConj]]=Map[If[ListQ[#],#[[1]],#]&,amplitudes[[selfConj]],{2}];

nAmps=numAmps[system];
If[!showFactors,amplitudes=KeyDrop[#,{"Binary indices","mu","CG"}]&/@amplitudes,Null]; (* show/hide internal factors from display *)
signs=If[system[["p factor"]]==1,{"-","+"},{"+","-"}];

Print["Amplitude table","\n",
"Number of amplitudes: ",nAmps,"\n",
"a/s definitions: \!\(\*SubscriptBox[\(a\), \(i\)]\) = \!\(\*SubscriptBox[\(A\), \(i\)]\) ",signs[[1]]," \!\(\*OverscriptBox[SubscriptBox[\(A\), \(i\)], \(_\)]\), \!\(\*SubscriptBox[\(s\), \(i\)]\) = \!\(\*SubscriptBox[\(A\), \(i\)]\) ",signs[[2]]," \!\(\*OverscriptBox[SubscriptBox[\(A\), \(i\)], \(_\)]\)","\n",
"\[CapitalDelta]/\[CapitalSigma] definitions: \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(i\)]\) = |\!\(\*SubscriptBox[\(A\), \(i\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\) ",signs[[1]]," |\!\(\*OverscriptBox[SubscriptBox[\(A\), \(i\)], \(_\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\), \!\(\*SubscriptBox[\(\[CapitalSigma]\), \(i\)]\) = |\!\(\*SubscriptBox[\(A\), \(i\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\) ",signs[[2]]," |\!\(\*OverscriptBox[SubscriptBox[\(A\), \(i\)], \(_\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)","\n",
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


(* Returns the number of amplitude or squared amplitude sum rules at each order of breaking *)
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
Options[printSRs]={showSRs->True,expandSRs->True,CKM->False,b->All}; 
printSRs[system_,{pwr_,label_,symbols__Symbol},opts:OptionsPattern[]]:=Module[{sysVal=Evaluate[system],showSRs=OptionValue[showSRs],expandSRs=OptionValue[expandSRs],CKM=OptionValue[CKM],b=OptionValue[b],formatSRMats,ampsToVectors,printWrittenSRs,syms={symbols},squared,p,SRs,amplitudes,indices,selfConj,physAmps,ampVectors,writtenSRs},
(* Function definitions *)
(* Double widths of SRs matrices, add CKM and p factors *)
formatSRMats[amps_]:=Module[{factorsMats,delSelfConj},
factorsMats={DiagonalMatrix@Flatten@Table[{1,-p},Length[amps]],DiagonalMatrix@Flatten@Table[{1,p},Length[amps]]}; (* {a mat,s mat} *)
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
If[label=="Processes"&&!KeyExistsQ[amps[[1]],"Processes"],Message[printSRs::invalidformat];Return[$Failed],Null];

ampsToVector[ampSym_Symbol]:=Module[{vector},
vector=Map[ampSym[#]&,
If[!physAmps,
(Switch[label,
"Processes",(Message[printSRs::invalidformat];Return[$Failed]),
"QNs",amps[[All,"QNs",1]], (* a/s-type amps with QNs *)
"n-tuples",amps[[All,"n-tuples",1]], (* a/s-type amps with n-tuples *)
"Nodes",amps[[All,"Nodes"]], (* a/s-type amps with nodes *)
"Binary indices",amps[[All,"Binary indices",1]], (* a/s-type amps with numbered subscripts *)
_String/;KeyExistsQ[amps[[1]],label],amps[[All,label]], (* custom amplitude format *)
_,(Message[printSRs::invalidformat];Return[$Failed])
]
),
(Switch[label,
"Processes",Flatten@amps[[All,"Processes"]], (* A amps with physical processes *)
"QNs",Flatten@amps[[All,"QNs"]], (* A amps with QNs *)
"n-tuples",Flatten@amps[[All,"n-tuples"]], (* A amps with n-tuples *)
"Nodes",(Message[printSRs::invalidformat];Return[$Failed]),
"Binary indices",Flatten@amps[[All,"Binary indices"]], (* A amps with numbered subscripts *)
_String/;KeyExistsQ[amps[[1]],label],Flatten@amps[[All,label]], (* custom amplitude format *)
_,(Message[printSRs::invalidformat];Return[$Failed])
]
)
]
];

If[physAmps&&CKM,vector=vector/Flatten@amps[[All,"CKM"]],Null];
If[physAmps&&squared,vector=#^2&/@vector,Null]; (* note that this is actually |A|^2, not AA *)

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

(* Print SRs with formatting and SR info *)
printWrittenSRs[]:=Module[{bList,numSRsList,matForm,rule},
bList=listbOrders[b,Length[SRs]];
numSRsList=numSRs[sysVal,squared,Sequence@@FilterRules[{opts},Options[numSRs]]][[bList+1]];
matForm:=If[expandSRs,MatrixForm[#]&,MatrixForm[#[[2;;]],TableHeadings->{None,#[[1]]}]&];
rule=If[label=="Processes"||label=="QNs",
expr_Symbol[i_]/;expr=!=List:>StringJoin[ToString[expr],"(",i,")"],
expr_Symbol[i_]/;expr=!=List:>Subscript[expr,i]
];

If[squared,Print["Squared amplitude sum rules"],Print["Amplitude sum rules"]];
MapIndexed[
Print["b = ",bList[[#2[[1]]]],"\n",
"Number of SRs: ",numSRsList[[#2[[1]]]],"\n",
If[#1=={},"No sum rules at this order.",matForm@If[expandSRs,#1/.rule,ReplacePart[#1,1->(#1[[1]]/.rule)]]]
]&,
writtenSRs[[bList+1]]
];
];



(* Function calls *)
If[label==None,Message[printSRs::invalidformat];Return[$Failed],Null];

squared=(pwr==2);
p=sysVal[["p factor"]];
SRs=If[squared,
If[KeyExistsQ[sysVal,"A2SRs"],sysVal[["A2SRs"]],(Message[printSRs::missingkey];Return[$Failed])],
sysVal[["ASRs"]]
]; (* set SRs list to either ASRs or A2SRs *)
amplitudes=sysVal[["Amplitudes"]];
indices=amplitudes[[All,"Binary indices"]];
selfConj=Table[If[indices[[i,1]]==indices[[i,2]],i,Nothing],{i,Length[indices]}];

physAmps=Switch[Length[syms],
1,True, (* A/A^2 *)
2,False, (* a/s/\[CapitalDelta]/\[CapitalSigma] *)
_,(Message[printSRs::invalidformat];Return[$Failed])
]; (* indicate whether amplitudes are phys (True) or group-theoretic (False) *)
If[physAmps,
formatSRMats[amplitudes],
If[Length[selfConj]>0,
SRs=MapIndexed[If[#1=={},{},If[EvenQ[p]==OddQ[First[#2]],Transpose[Delete[Transpose[#1],selfConj[[1]]]],#1]]&,SRs],
Null
]
]; (* for phys amps, double widths of SRs matrices, add CKM and p factors. for group-theoretic, correct for self-conj amps *)
ampVectors=ampsToVectors[amplitudes]; (* return list of two vectors of amplitudes; the two vecs are the same in phys case *)

writtenSRs=If[expandSRs,
MapIndexed[If[#1=={},{},If[OddQ[First[#2]],#1 . ampVectors[[1]],#1 . ampVectors[[2]]]]&,SRs],
MapIndexed[If[#1=={},{},If[OddQ[First[#2]],Prepend[#1,ampVectors[[1]]],Prepend[#1,ampVectors[[2]]]]]&,SRs]
]; (* combine amp vecs with SRs matrices *)

If[showSRs,printWrittenSRs[],Null]; (* print SRs *)
If[squared,system[["Formatted A2SRs"]]=writtenSRs,system[["Formatted ASRs"]]=writtenSRs];
writtenSRs
];

printSRs::invalidformat="Invalid or missing amplitude format.";
printSRs::missingkey="This system does not contain the A2SRs key.";


SetAttributes[printSystem,HoldFirst];
Options[printSystem]=Join[{showReps->True,showAmps->True,showASRs->True,showA2SRs->True},Options[printAmps],Options[printSRs]];
printSystem[system_,ASRFormat_:None,A2SRFormat_:None,opts:OptionsPattern[]]:=Module[{sysVal=Evaluate[system],showReps=OptionValue[showReps],showAmps=OptionValue[showAmps],showASRs=OptionValue[showASRs],showA2SRs=OptionValue[showA2SRs],irreps},
irreps=sysVal[["Irreps"]];
Print["System: ",Sort@Flatten@irreps];

If[showReps,
(Print["-------------------------"];
Print["Number of would-be doublets: ",sysVal[["n doublets"]],"\n",
"In: ",irreps[[1]],"\n",
"H: ",irreps[[2]],"\n",
"Out: ",irreps[[3]]
];
),
Null
];

If[showAmps,
(Print["-------------------------"];
printAmps[sysVal,Sequence@@FilterRules[{opts},Options[printAmps]]];
),
Null
];

If[showASRs,Print["-------------------------"],Null];
printSRs[sysVal,ASRFormat,Sequence@@FilterRules[{opts},Options[printSRs]],showSRs->showASRs];

If[showA2SRs,Print["-------------------------"],Null];
printSRs[sysVal,A2SRFormat,Sequence@@FilterRules[{opts},Options[printSRs]],showSRs->showA2SRs];

system=sysVal;
sysVal
];


Options[generateA2SRs]={phys->False};
generateA2SRs[in_,h_,out_,opts:OptionsPattern[]]:=Module[{system,phys=OptionValue[phys],ASRs,nAmps,nASRs,nThetas,findA2SRMat,A2SRs,nA2SRs,thetaRanks,xMats},
system=generateASRs[in,h,out,Sequence@@FilterRules[{opts},Options[generateASRs]]];
ASRs=system[["ASRs"]];
nAmps=numAmps[system,True];
nASRs=system[["n ASRs"]];
nThetas=((nAmps-nASRs)(nAmps-nASRs-1))/2;

(* Find the A^2SR matrix for a given ASR matrix *)
findA2SRMat[ASRMat_]:=Module[{ASRMatRR,pivotCols,freeCols,freeMat,groupedIndices,indices,colPairs,xMat,xMatNullSpace,groupedA2SR,thetaRank,A2SRMat},
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

thetaRank=MatrixRank[xMat];
A2SRMat=If[groupedA2SR=={},{},Transpose@Table[Transpose[groupedA2SR][[indices[[i]]]],{i,Length[indices]}]]; (* rearranges A^2SR mat into original order *)

{A2SRMat,thetaRank,xMat}
];

{A2SRs,thetaRanks,xMats}=Transpose@Map[findA2SRMat,ASRs];
nA2SRs=Table[Length[A2SRs[[i]]],{i,Length[A2SRs]}];

AssociateTo[system,<|"n A2SRs"->nA2SRs,"A2SRs"->A2SRs,"n thetas"->nThetas,"Theta ranks"->thetaRanks,"Thetas"->xMats|>]
];


End[]


EndPackage[]
