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
- \"n doublets\" (Real): Number of would-be doublets.
- \"Irreps\" (List): Inputted U-spin representations (List of Reals) in {incoming reps, H reps, outgoing reps} format.
- \"n amps\" (Real): Number of amplitudes in the system.
- \"Amplitudes\" (Association): Contains all amplitudes in the system. Keys and values:
	- \"Processes\" (List): Contains physical processes constructed from particle names (String) for a process and its U-spin conjugate process. Only available for physical systems.
	- \"QNs\" (List): Contains processes written using m quantum numbers (String), where m is the third component of U-spin, for a process and its U-spin conjugate process.
	- \"n-tuple\" (String): Representation of an amplitude pair using a comma-separated tuple of substrings. Each substring is comprised of '-'s and '+'s and encodes u and m QNs of a component of a participating multiplet. Signs are inverted for initial state and Hamiltonian components. By convention, one n-tuple (beginning with a '-' sign) is used to represent an amplitude pair: an amplitude and its U-spin conjugate amplitude.
	- \"Node\" (String): Representation of an amplitude pair using the coordinate notation of the lattice used to derive sum rules.
	- \"Binary indices\" (List): Contains decimal indices (Real) derived from converting the n-tuple and its U-spin conjugate into binary numbers through '-' <-> 0 and '+' <-> 1.
	- \"q factor\" (Real): (-1)^q_i factor for an amplitude pair accounting for movements of representations between the in state, Hamiltonian, and out state.
	- \"p factor\" (Real): (-1)^p parity factor for the system.
	- \"mu\" (Real): mu-factor for the node in the lattice used to derive sum rules.
	- \"CG\" (Real): Clebsch-Gordan coefficient from symmetrization for systems without doublets. Equal to 1 for all amplitudes for a system with at least one doublet.
	- \"CKM\" (List): Contains weak interaction factors (Real) from the Hamiltonian. Only appears for physical systems.
- \"n ASRs\" (List): Contains the number of amplitude sum rules (Real) at each order of breaking.
- \"ASRs\" (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.";

numAmps::usage="numAmps[amplitudes] returns the total number of amplitudes in the system.";
numAmps::details=
"Arguments:
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.

Options:
nPairs (True|False): Indicates whether to return the number of amplitudes (False) or amplitude pairs (True). Default: nPairs: False.

Returns:
The total number of amplitudes or amplitude pairs in the system (Real).";

labelAmps::usage="labelAmps[amplitudes,colName,labels] adds a column of user-defined labels to amplitudes."
labelAmps::details=
"Arguments:
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.
colName (String): Name of new column to add to amplitudes association.
labels (List): Contains labels (Any) for each amplitude or amplitude pair. Number of labels must equal number of amplitudes or amplitude pairs.

Options:
labeling (String): Labeling mode to indicate whether user is labeling single amplitudes (\"Amplitudes\") or amplitude pairs (\"Amplitude pairs\"). Default: labeling -> \"Amplitudes\".

Returns:
amplitudes (Association): The inputted amplitudes association, modified to include a new column of user-defined labels. New/modified keys and values:
- colName (List|Any): Either contains user-defined labels (Any) for an amplitude and its U-spin conjugate or a single label (Any) for an amplitude pair.";

unlabelAmps::usage="unlabelAmps[amplitudes,colNames] removes columns from amplitudes."
unlabelAmps::details=
"Arguments:
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.
colNames (String|List): Name(s) of column(s) to remove from amplitudes association.

Returns:
amplitudes (Association): The inputted amplitudes association, modified to remove specified columns.";

printAmps::usage="printAmps[amplitudes] prints a table of amplitudes.";
printAmps::details=
"Arguments:
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.

Options:
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.";

numASRs::usage="numASRs[ASRs] returns the number of amplitude sum rules at each order of breaking.";
numASRs::details=
"Arguments:
ASRs (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.

Options:
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
Number of amplitude sum rules at each order of breaking (List of Reals).";

printASRs::usage="printASRs[ASRs,amplitudes] prints amplitude sum rules at each order of breaking.";
printASRs::details=
"Arguments:
ASRs (List): Contains matrices of amplitude sum rule coefficients (Real) corresponding to each order of breaking.
amplitudes (Association): An amplitude association. See the documentation for generateASRs for details.

Options:
showASRs (True|False): Default: showASRs->True.
expandASRs (True|False): Indicates whether to write each row of a sum rules matrix as an algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandASRs->True.
ampFormat (String): Specified format for displaying amplitudes. Options are a/s-type amplitudes with n-tuples (\"a/s n-tuple\"), a/s-type amplitudes with numbered indices (\"a/s indices\"), a/s-type amplitudes with nodes (\"a/s nodes\"), A amplitudes with numbered indices (\"A indices\"), and A amplitudes with physical process (\"A physical\"). Default: ampFormat->\"a/s n-tuple\".
CKM (True|False): Indicates whether to include CKM factors in the sum rules (True) or not (False). Default: CKM->False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
\"Formatted ASRs\" (List): Contains Lists/matrices of amplitude sum rules/amplitude sum rule coefficients formatted according to sum rule printing settings";

printSystem::usage="printSystem[system] prints information about the system's representations, amplitudes, and amplitude sum rules and modifies the system to include formatted sum rules.";
printSystem::details=
"Arguments:
system (Association): A system association. See the documentation for generateASRs for details.

Options:
showReps (True|False): Default: showReps->True.
-----
showAmps (True|False): Default: showAmps->True.
showFactors (True|False): Indicates whether to print internal calculation factors from the sum rule algorithm in the amplitudes table (True) or not (False). Not to be confused with the amplitude sum rule matrices. Default: showFactors->False.
-----
showASRs (True|False): Default: showASRs->True.
expandASRs (True|False): Indicates whether to write each row of a sum rules matrix as an algebraic expression of amplitudes (True) or to keep each row as a list of coefficients (False). Default: expandASRs->True.
ampFormat (String): Specified format for displaying amplitudes. Options are a/s-type amplitudes with n-tuples (\"a/s n-tuple\"), a/s-type amplitudes with numbered indices (\"a/s indices\"), a/s-type amplitudes with nodes (\"a/s nodes\"), A amplitudes with numbered indices (\"A indices\"), and A amplitudes with physical process (\"A physical\"). Default: ampFormat->\"a/s n-tuple\".
CKM (True|False): Indicates whether to include CKM factors in the sum rules (True) or not (False). Default: CKM->False.
b (All|Real|List): Breaking order(s) at which to print sum rules. User can print sum rules to all possible orders of breaking (All), at a particular order (Real s.t. 0 <= b <= highest order of breaking), or over a range of orders of breaking ({start b (min: 0), end b (max: highest order of breaking, or All), increment}). Default: b->All.

Returns:
system (Association): The inputted system association, modified to include formatted amplitude sum rules. New/modified keys and values:
- \"Formatted ASRs\" (List): Contains Lists/matrices of amplitude sum rules/amplitude sum rule coefficients formatted according to sum rule printing settings.";


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


(* Extracts amplitudes from Python System object *)
Options[extractAmps]={partVal->{}};
extractAmps[system_,OptionsPattern[]]:=Module[{amplitudes,colNames,extractParticles,partVal=OptionValue[partVal]},
amplitudes=pyEval["extract_amps",system];
colNames={"Processes","QNs","n-tuple","Node","Binary indices","q factor","p factor","mu","CG"};
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
generateASRs[in_,h_,out_,OptionsPattern[]]:=Module[{system,phys=OptionValue[phys],particles,ASRs,amplitudes,dupPairs,factorsMat,n,irreps,nAmps,nASRs},
{system,particles}=defineSystem[in,h,out,phys]; (* constructs Python System object *)

{system,ASRs,dupPairs}=pyEval["generate_srs",system]; (* M values for symmetrized system, contains duplicate amplitudes *)
amplitudes=extractAmps[system];

factorsMat=DiagonalMatrix[amplitudes[[All,"mu"]]*amplitudes[[All,"CG"]]];
ASRs=Map[# . factorsMat&,ASRs]; (* SRs for symmetrized system *)

(* Correct for duplicate amplitudes from symmetrization *)
Do[ASRs[[b]][[All,dupPairs[[All,1]]]]+=ASRs[[b]][[All,dupPairs[[All,2]]]],{b,Length[ASRs]}]; (* adds together cols of duplicate amp pairs *)
ASRs=Transpose[Delete[Transpose[#],List/@dupPairs[[All,2]]]]&/@ASRs; (* deletes duplicate cols *)
system=pyEval["remove_dups",{system,dupPairs}];
amplitudes=extractAmps[system,partVal->particles];

ASRs=Map[Cases[Except@{0..}],Map[RowReduce,ASRs]];

{n,irreps}=pyEval["extract_sys",system];
nAmps=numAmps[amplitudes];
nASRs=numASRs[ASRs];
system=<|"n doublets"->n,"Irreps"->irreps,"n amps"->nAmps,"Amplitudes"->amplitudes,"n ASRs"->nASRs,"ASRs"->ASRs|>;

system
];


(* Returns the total number of amplitudes in the system *)
numAmps[amplitudes_,nPairs_:False]:=If[!nPairs,Length@DeleteDuplicates@Flatten@amplitudes[[All,"Binary indices"]],Length@Flatten@amplitudes[[All,"Node"]]];


(* Adds a column to amplitudes *)
SetAttributes[labelAmps,HoldFirst];
Options[labelAmps]={labeling->"Amplitudes"};
labelAmps[amplitudes_,colName_String,labels_List,OptionsPattern[]]:=Module[{pairs,nAmps,labelVals,labeling=OptionValue[labeling],indices,labelIndices},
pairs=If[labeling=="Amplitudes",False,True];
nAmps=numAmps[amplitudes,pairs];
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
amplitudes
];

labelAmps::arglen="The argument supplied to `1` has an incorrect number of labels. Expected `2` labels, got `3`. Check that the labeling option is correct.";
labelAmps::badmode="Unknown labeling mode `1`. Use \"Amplitudes\" or \"Amplitude pairs\".";


(* Removes columns from amplitudes *)
SetAttributes[unlabelAmps,HoldFirst];
unlabelAmps[amplitudes_,colNames_]:=(amplitudes=KeyDrop[colNames]/@amplitudes);


(* Prints a table of amplitudes *)
Options[printAmps]={showFactors->False};
printAmps[amplitudes_,OptionsPattern[]]:=Module[{amplitudesVal=amplitudes,indices,selfConj,nAmps,showFactors=OptionValue[showFactors]},
(* Delete self-conjugate duplicate values from display *)
indices=amplitudesVal[[All,"Binary indices"]];
selfConj=Table[If[indices[[i,1]]==indices[[i,2]],i,Nothing],{i,Length[indices]}];
amplitudesVal[[selfConj]]=Map[If[ListQ[#],#[[1]],#]&,amplitudesVal[[selfConj]],{2}];

nAmps=numAmps[amplitudesVal];

If[!showFactors,amplitudesVal=KeyDrop[#,{"Binary indices","q factor","p factor","mu","CG"}]&/@amplitudesVal,Null]; (* show/hide internal factors from display *)

Print["Amplitude table","\n",
"Number of amplitudes: ",nAmps,"\n",
TableForm[Prepend[Values/@amplitudesVal,Keys[amplitudesVal[[1]]]]]
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
]


(* Returns the number of amplitude sum rules at each order of breaking *)
Options[numASRs]={b->All};
numASRs[ASRs_,OptionsPattern[]]:=Module[{b=OptionValue[b],indices,sublist},
If[b===All,
Table[Length[ASRs[[i]]],{i,Length[ASRs]}],
(indices=listbOrders[b,Length[ASRs]]+1;
sublist=ASRs[[indices]];
Table[Length[sublist[[i]]],{i,Length[sublist]}])
]
];


(* Prints amplitude sum rules at each order of breaking *)
Options[printASRs]={showASRs->True,expandASRs->True,ampFormat->"a/s n-tuple",CKM->False,b->All};
printASRs[ASRs_,amplitudes_,OptionsPattern[]]:=Module[{showASRs=OptionValue[showASRs],expandASRs=OptionValue[expandASRs],ampFormat=OptionValue[ampFormat],CKM=OptionValue[CKM],b=OptionValue[b],physicalAmps,ampsToVector,ampVector,writtenASRs,bList,numASRsList,matForm,rule},
physicalAmps[amp_,asType_,col_,CKM_]:=Module[{asSign,CKMfactors},
asSign=If[SymbolName[asType]=="a",-1,+1];
CKMfactors=If[CKM,{amp[["CKM",1]],amp[["CKM",2]]},{1,1}];
amp[["q factor"]](A[amp[[col,1]]]/CKMfactors[[1]]+asSign*amp[["p factor"]]*A[amp[[col,2]]]/CKMfactors[[2]])
];

ampsToVector[amps_,asType_,ampFormat_,CKM_]:=Module[{vector,x,label},
If[ampFormat=="A physical"&&!KeyExistsQ[amps[[1]],"Processes"],Message[printASRs::invalidformat];Return[$Failed],Null];
{x,label}=StringSplit[ampFormat,WhitespaceCharacter,2];

vector=Switch[ampFormat,
"a/s n-tuple",Map[asType[#]&,amps[[All,"n-tuple"]]], (* a/s-type amps with n-tuples; default option *)
"a/s nodes",Map[asType[#]&,amps[[All,"Node"]]], (*a/s-type amps with nodes*)
"a/s indices",Map[asType[#]&,amps[[All,"Binary indices",1]]], (* a/s-type amps with numbered subscripts *)
"A physical",Map[physicalAmps[#,asType,"Processes",CKM]&,amps], (* A amps with physical processes *)
"A QNs",Map[physicalAmps[#,asType,"QNs",CKM]&,amps], (* A amps with QNs *)
"A indices",Map[physicalAmps[#,asType,"Binary indices",CKM]&,amps], (* A amps with numbered subscripts *)
_String/;KeyExistsQ[amplitudes[[1]],label],If[x=="a/s",Map[asType[#]&,amps[[All,label]]],Map[physicalAmps[#,asType,label,CKM]&,amps]], (* custom amplitude format *)
_,(Message[printASRs::invalidformat];Return[$Failed])
];
vector
];

ampVector={ampsToVector[amplitudes,a,ampFormat,CKM],ampsToVector[amplitudes,s,ampFormat,CKM]};

writtenASRs=If[expandASRs,
MapIndexed[If[OddQ[First[#2]],#1 . ampVector[[1]],#1 . ampVector[[2]]]&,ASRs],
MapIndexed[If[OddQ[First[#2]],Prepend[#1,ampVector[[1]]],Prepend[#1,ampVector[[2]]]]&,ASRs]
];

(* Cosmetics for printing *)
bList=listbOrders[b,Length[ASRs]];
numASRsList=numASRs[ASRs][[bList+1]];
matForm:=If[expandASRs,MatrixForm[#]&,MatrixForm[#[[2;;]],TableHeadings->{None,#[[1]]}]&];
rule=If[ampFormat=="A physical"||ampFormat=="A QNs",
x_Symbol[i_]/;x=!=List:>StringJoin[ToString[x],"(",i,")"],
x_Symbol[i_]/;x=!=List:>Subscript[x,i]
];

If[showASRs,
(Print["Sum rules"];
MapIndexed[
Print["b = ",bList[[#2[[1]]]],"\n",
"Number of ASRs: ",numASRsList[[#2[[1]]]],"\n",
matForm@If[expandASRs,#1/.rule,ReplacePart[#1,1->(#1[[1]]/.rule)]]
]&,
writtenASRs[[bList+1]]
];
),
Null
];

writtenASRs
];

printASRs::invalidformat="Invalid amplitude format for this system.";


SetAttributes[printSystem,HoldFirst];
Options[printSystem]=Join[{showReps->True,showAmps->True},Options[printAmps],Options[printASRs]];
printSystem[system_,opts:OptionsPattern[]]:=Module[{showReps=OptionValue[showReps],showAmps=OptionValue[showAmps],showASRs=OptionValue[showASRs],irreps,writtenASRs},
irreps=system[["Irreps"]];

Print["System: ",Sort@Flatten@irreps];

If[showReps,
(Print["-------------------------"];
Print["Number of would-be doublets: ",system[["n doublets"]],"\n",
"In: ",irreps[[1]],"\n",
"H: ",irreps[[2]],"\n",
"Out: ",irreps[[3]]
];
),
Null];

If[showAmps,
(Print["-------------------------"];
printAmps[system[["Amplitudes"]],Sequence@@FilterRules[{opts},Options[printAmps]]];
),
Null];

If[showASRs,Print["-------------------------"],Null];
writtenASRs=printASRs[system[["ASRs"]],system[["Amplitudes"]],Sequence@@FilterRules[{opts},Options[printASRs]]];
system[["Formatted ASRs"]]=Evaluate[writtenASRs];

system
];


End[]


EndPackage[]
