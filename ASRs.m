(* ::Package:: *)

BeginPackage["ASRs`"];


ASRsHelp::usage="ASRsHelp[function] prints extended documentation on a function's arguments, options, and outputs.";
ASRsHelp::details=
"Arguments:
function (Symbol): The name of a function from the ASRs Mathematica package."

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
	- \"n-tuple\" (String): Representation of an amplitude pair using a comma-separated tuple of substrings. Each substring is comprised of '-'s and '+'s and encodes u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components. By convention, one n-tuple (beginning with a '-' sign) is used to represent an amplitude pair: an amplitude and its U-spin conjugate amplitude.
	- \"Node\" (String): Representation of an amplitude pair using the coordinate notation of the lattice used to derive sum rules.
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
ampFormat (String): Specified format for displaying amplitudes. Options are a/s-type amplitudes with n-tuples (\"a/s n-tuple\"), a/s-type amplitudes with numbered indices (\"a/s indices\"), a/s-type amplitudes with nodes (\"a/s nodes\"), A amplitudes with numbered indices (\"A indices\"), and A amplitudes with physical process (\"A physical\"). User-defined labels for a column of amplitudes (or amplitude pairs) can also be used using the syntax \"A col\" (or \"a/s col\" for pairs), where \"col\" is the name of the column containing custom labels. Default: ampFormat->\"a/s n-tuple\".
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
showReps (True|False): Default: showReps->True.
-----
showAmps (True|False): Default: showAmps->True.
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.
-----
showASRs (True|False): Default: showASRs->True. Note: ASRs and A2SRs use the same formatting options.
showA2SRs (True|False): Default: showA2SRs->False. Note: ASRs and A2SRs use the same formatting options.
expandSRs (True|False): Indicates whether to write each row of a sum rules matrix as an algebraic expression of (squared) amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandSRs->True.
ampFormat (String): Specified format for displaying (squared) amplitudes. Options are a/s-type amplitudes with n-tuples (\"a/s n-tuple\"), a/s-type amplitudes with numbered indices (\"a/s indices\"), a/s-type amplitudes with nodes (\"a/s nodes\"), A amplitudes with numbered indices (\"A indices\"), and A amplitudes with physical process (\"A physical\"). Default: ampFormat->\"a/s n-tuple\".
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
	- \"n-tuple\" (String): Representation of an amplitude pair using a comma-separated tuple of substrings. Each substring is comprised of '-'s and '+'s and encodes u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components. By convention, one n-tuple (beginning with a '-' sign) is used to represent an amplitude pair: an amplitude and its U-spin conjugate amplitude.
	- \"Node\" (String): Representation of an amplitude pair using the coordinate notation of the lattice used to derive sum rules.
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


A=.;
a=.;
s=.;


Begin["`Private`"];


ASRsHelp[function_Symbol]:=Module[{usage, details},
usage=Quiet[MessageName[function,"usage"]];
details=Quiet[MessageName[function,"details"]];

If[StringQ[usage],
(Print[usage];
If[StringQ[details],Print[details]];
),
Print["No documentation found for ",function,"."];
];
];


$ASRsSession=.;


startPythonSession[session_,path_String]:=Module[{},
$ASRsSession=If[MatchQ[session,_ExternalSessionObject],session,Message[startPythonSession::nosession];Return[$Failed]];
If[FileExistsQ[path],Null,Message[startPythonSession::nofile,path];Return[$Failed]];

ExternalEvaluate[$ASRsSession,File[path]]
];

startPythonSession::nosession="No active Python session provided.";
startPythonSession::nofile="File `1` does not exist.";


pyEval[expr_,args_:<||>]:=ExternalEvaluate[$ASRsSession,<|"Command"->expr,"Arguments"->args|>];


(* Extracts amplitudes from system (a Python System object) *)
Options[extractAmps]={partVal->{}};
extractAmps[system_,OptionsPattern[]]:=Module[{amplitudes,colNames,extractParticles,partVal=OptionValue[partVal]},
amplitudes=pyEval["System.extract_amps",system];
colNames={"Processes","QNs","n-tuple","Node","Binary indices","q factor","mu","CG"};
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


defineSystem[in_,h_,out_,phys_]:=Module[{repSpins,reps,system,particles},
repSpins[list_]:=Table[(Length[list[[i]]]-1)/2,{i,Length[list]}];

If[phys,
(particles={in,h,out};
reps=Map[repSpins,particles];
),
(reps={in,h,out};
particles={};
)
];

system=pyEval["define_system",{reps,phys}];
{system,particles}
];


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
simplifyFactors[mat_]:=(#/Sqrt[Apply[GCD,DeleteCases[#,0]^2]])&/@mat;
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
Length@Flatten@amplitudes[[All,"Node"]]
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
"a/s definitions: \!\(\*SubscriptBox[\(a\), \(i\)]\) = \!\(\*SubscriptBox[\(A\), \(i\)]\) ",signs[[1]]," \!\(\*SubscriptBox[\(A\), \(l\)]\), \!\(\*SubscriptBox[\(s\), \(i\)]\) = \!\(\*SubscriptBox[\(A\), \(i\)]\) ",signs[[2]]," \!\(\*SubscriptBox[\(A\), \(l\)]\)","\n",
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
Options[numSRs]={b->All,squared->False};
numSRs[system_,OptionsPattern[]]:=Module[{squared=OptionValue[squared],b=OptionValue[b],SRs,indices,sublist},
SRs=If[squared,
If[KeyExistsQ[system,"A2SRs"],system[["A2SRs"]],(Message[numSRs::missingkey];Return[$Failed])],
system[["ASRs"]]
];
If[b===All,
Table[Length[SRs[[i]]],{i,Length[SRs]}],
(indices=listbOrders[b,Length[SRs]]+1;
sublist=SRs[[indices]];
Table[Length[sublist[[i]]],{i,Length[sublist]}])
]
];

numSRs::missingkey="This system does not contain the A2SRs key.";


(* Prints sum rules at each order of breaking *)
SetAttributes[printSRs,HoldFirst];
Options[printSRs]={showSRs->True,expandSRs->True,ampFormat->"a/s n-tuple",CKM->False,b->All,squared->False};
printSRs[system_,opts:OptionsPattern[]]:=Module[{sysVal=Evaluate[system],SRs,amplitudes,showSRs=OptionValue[showSRs],expandSRs=OptionValue[expandSRs],ampFormat=OptionValue[ampFormat],CKM=OptionValue[CKM],b=OptionValue[b],squared=OptionValue[squared],pwr,physicalAmps,ampsToVector,ampVector,writtenSRs,bList,numSRsList,matForm,rule},
SRs=If[squared,
If[KeyExistsQ[sysVal,"A2SRs"],sysVal[["A2SRs"]],(Message[printSRs::missingkey];Return[$Failed])],
sysVal[["ASRs"]]
];
amplitudes=sysVal[["Amplitudes"]];
pwr=If[squared,2,1];

physicalAmps[amp_,asType_,col_,CKM_]:=Module[{asSign,CKMfactors},
asSign=If[SymbolName[asType]=="a",-1,+1];
CKMfactors=If[CKM,{amp[["CKM",1]],amp[["CKM",2]]},{1,1}];
A[amp[[col,1]]]/CKMfactors[[1]]+asSign*sysVal[["p factor"]]*A[amp[[col,2]]]/CKMfactors[[2]]
];

ampsToVector[amps_,asType_,ampFormat_,CKM_,pwr_]:=Module[{vector,x,label},
If[ampFormat=="A physical"&&!KeyExistsQ[amps[[1]],"Processes"],Message[printSRs::invalidformat];Return[$Failed],Null];
{x,label}=StringSplit[ampFormat,WhitespaceCharacter,2];

vector=Switch[ampFormat,
"a/s n-tuple",Map[asType[#]&,amps[[All,"n-tuple"]]], (* a/s-type amps with n-tuples; default option *)
"a/s nodes",Map[asType[#]&,amps[[All,"Node"]]], (*a/s-type amps with nodes*)
"a/s indices",Map[asType[#]&,amps[[All,"Binary indices",1]]], (* a/s-type amps with numbered subscripts *)
"A physical",Map[physicalAmps[#,asType,"Processes",CKM]&,amps], (* A amps with physical processes *)
"A QNs",Map[physicalAmps[#,asType,"QNs",CKM]&,amps], (* A amps with QNs *)
"A indices",Map[physicalAmps[#,asType,"Binary indices",CKM]&,amps], (* A amps with numbered subscripts *)
_String/;KeyExistsQ[amps[[1]],label],If[x=="a/s",Map[asType[#]&,amps[[All,label]]],Map[physicalAmps[#,asType,label,CKM]&,amps]], (* custom amplitude format *)
_,(Message[printSRs::invalidformat];Return[$Failed])
];
vector^pwr (* pwr works for A and A^2 *)
];

ampVector={ampsToVector[amplitudes,a,ampFormat,CKM,pwr],ampsToVector[amplitudes,s,ampFormat,CKM,pwr]};

writtenSRs=If[expandSRs,
MapIndexed[If[#1=={},{},If[OddQ[First[#2]],#1 . ampVector[[1]],#1 . ampVector[[2]]]]&,SRs],
MapIndexed[If[#1=={},{},If[OddQ[First[#2]],Prepend[#1,ampVector[[1]]],Prepend[#1,ampVector[[2]]]]]&,SRs]
];

(* Cosmetics for printing *)
bList=listbOrders[b,Length[SRs]];
numSRsList=numSRs[sysVal,Sequence@@FilterRules[{opts},Options[numSRs]]][[bList+1]];
matForm:=If[expandSRs,MatrixForm[#]&,MatrixForm[#[[2;;]],TableHeadings->{None,#[[1]]}]&];
rule=If[ampFormat=="A physical"||ampFormat=="A QNs",
x_Symbol[i_]/;x=!=List:>StringJoin[ToString[x],"(",i,")"],
x_Symbol[i_]/;x=!=List:>Subscript[x,i]
];

If[showSRs,
(If[squared,Print["Squared amplitude sum rules"],Print["Amplitude sum rules"]];
MapIndexed[
Print["b = ",bList[[#2[[1]]]],"\n",
"Number of SRs: ",numSRsList[[#2[[1]]]],"\n",
If[#1=={},"No sum rules at this order.",matForm@If[expandSRs,#1/.rule,ReplacePart[#1,1->(#1[[1]]/.rule)]]]
]&,
writtenSRs[[bList+1]]
];
),
Null
];

If[squared,system[["Formatted A2SRs"]]=writtenSRs,system[["Formatted ASRs"]]=writtenSRs];
writtenSRs
];

printSRs::invalidformat="Invalid amplitude format for this system.";
printSRs::missingkey="This system does not contain the A2SRs key.";


SetAttributes[printSystem,HoldFirst];
Options[printSystem]=Join[{showReps->True,showAmps->True,showASRs->True,showA2SRs->False},Options[printAmps],Options[printSRs]];
printSystem[system_,opts:OptionsPattern[]]:=Module[{sysVal=Evaluate[system],showReps=OptionValue[showReps],showAmps=OptionValue[showAmps],showASRs=OptionValue[showASRs],showA2SRs=OptionValue[showA2SRs],irreps,writtenSRs},
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
printSRs[sysVal,Sequence@@FilterRules[{opts},Options[printSRs]],showSRs->showASRs];

If[showA2SRs,Print["-------------------------"],Null];
printSRs[sysVal,Sequence@@FilterRules[{opts},Options[printSRs]],showSRs->showA2SRs,squared->True];

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

AssociateTo[system,<|"n A2SRs"->nA2SRs,"A2SRs"->A2SRs,"n thetas"->nThetas,"Theta ranks"->thetaRanks,"Thetas"->xMats|>] (************************************** do we still need these last two? *)
];


End[]


EndPackage[]
