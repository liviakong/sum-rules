(* ::Package:: *)

BeginPackage["ASRs`"];


ASRsHelp::usage="ASRsHelp[function] prints extended documentation on how to use a function."

initializePythonSession::usage="initializePythonSession[session,path] checks for a valid Python session and file path and loads in the Python file.";

defineSystem::usage="defineSystem[in,h,out] constructs a U-Spin System object in Python.";
defineSystem::details=
"Arguments:
in (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the incoming state
h (List): Contains U-spins (Real) OR CKM factors (List of Symbols) in the Hamiltonian
out (List): Contains U-spins (Real) OR particle multiplets (List of Strings) in the outgoing state

Options:
phys (True|False): Set to False (default) to input U-spins and True to input particle multiplets/CKM factors

Returns:
system (Python System object): U-spin system represented as a Python System object
particles (List): Contains particle multiplets (List of Strings/Symbols) in the in state, Hamiltonian, and out state. Only returned if phys->True.";

generateASRs::usage="generateASRs[system,particles->{}] takes in a Python U-spin System object and optionally a list of particle multiplets and returns the updated Python U-spin System object, a list of amplitude sum rule coefficients at each order of breaking, and an association containing all amplitudes in the system.";

numAmps::usage="numAmps[amplitudes] returns the total number of amplitudes in the system.";

printAmps::usage="printAmps[amplitudes] prints a table of amplitudes and returns the number of amplitudes and the amplitude association.";

numASRs::usage="numASRs[ASRs] returns the number of amplitude sum rules at each order of breaking.";

printASRs::usage="printASRs[ASRs,amplitudes,takeProd->False,ampFormat->\"a/s n-tuple\",CKM->False,b->All] prints amplitude sum rules at each order of breaking and returns the number of sum rules at each order of breaking and the expanded or coefficient form of the sum rules, depending on user choice. Valid options for takeProd: True, False; ampFormat: \"a/s n-tuple\", \"a/s indices\", \"a/s nodes\", \"A indices\", \"A physical\"; CKM: True, False; b: an integer 0 <= b <= highest order of breaking, All, {start b (min: 0), end b (max: highest order of breaking, or All), increment}.";

printSystem::usage="printSystem[system,ASRs,amplitudes] prints information about the system's representations, amplitudes, and amplitude sum rules.";
printSystem::details=
"Arguments:
system
ASRs
amplitudes

Options:
Uses the same options as printASRs to format the amplitude sum rules

Returns:
Association: <|\"n doublets\"->n,\"Irreps\"->irreps,\"n amps\"->nAmps,\"Amplitudes\"->amplitudes,\"n ASRs\"->nASRs,\"ASRs\"->writtenASRs|>";


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


initializePythonSession[session_,path_String]:=Module[{},
$ASRsSession=If[MatchQ[session,_ExternalSessionObject],
session,
Message[initializePythonSession::nosession];Return[$Failed]];
If[FileExistsQ[path],
Null,
Message[initializePythonSession::nofile,path];Return[$Failed]];

ExternalEvaluate[$ASRsSession,File[path]]
];
initializePythonSession::nosession="No active Python session provided.";
initializePythonSession::nofile="File `1` does not exist.";


pyEval[expr_,args_:<||>]:=ExternalEvaluate[$ASRsSession,<|"Command"->expr,"Arguments"->args|>];


Options[defineSystem]={phys->False};
defineSystem[in_,h_,out_,OptionsPattern[]]:=Module[{phys=OptionValue[phys],repSpins,uIn,uH,uOut,pIn,pH,pOut,system,particles},
repSpins[list_]:=Table[(Length[list[[i]]]-1)/2,{i,Length[list]}];
If[phys,
({pIn,pH,pOut}={in,h,out};
{uIn,uH,uOut}=Map[repSpins,{pIn,pH,pOut}];
particles={pIn,pH,pOut};
),
{uIn,uH,uOut}={in,h,out};
];
system=pyEval["define_system",{{uIn,uH,uOut},phys}];
If[phys,{system,particles},system]
];


Options[extractAmps]={particles->{}};
extractAmps[system_,OptionsPattern[]]:=Module[{amplitudes,colNames,partVal=OptionValue[particles]},
amplitudes=pyEval["System.extract_amplitudes",system];
colNames={"Processes","Indices","n-tuple","Node","q factor","p factor","mu","CG"};
amplitudes=Map[AssociationThread[colNames,#]&]@amplitudes;

extractParticles[particles_,indices_]:=MapThread[MapThread[Part,{#1,#2}]&,{particles,indices}];
If[Length[partVal]>0,
(amplitudes[[All,"Processes"]]=Map[{extractParticles[partVal,#[[1]]],extractParticles[partVal,#[[2]]]}&,amplitudes[[All,"Processes"]]];
amplitudes[[All,"CKM"]]=Map[#[[2]]&,amplitudes[[All,"Processes"]],{2}];
amplitudes[[All,"Processes"]]=Map[{#[[1]],#[[3]]}&,amplitudes[[All,"Processes"]],{2}];
amplitudes[[All,"Processes"]]=Map[{StringRiffle[#[[1]]," "],StringRiffle[#[[2]]," "]}&,amplitudes[[All,"Processes"]],{2}];
amplitudes[[All,"Processes"]]=Map[{StringJoin[#[[1]]," \[Rule] ",#[[2]]]}&,amplitudes[[All,"Processes"]],{2}];
amplitudes[[All,"Processes"]]=Flatten/@amplitudes[[All,"Processes"]];
),
Null
]; (* particles involved in specific processes *)

amplitudes[[All,"mu"]]=Map[Sqrt[#[[1]]]*#[[2]]&,amplitudes[[All,"mu"]]]; (* mu factors *)
amplitudes[[All,"CG"]]=Map[ClebschGordan@@Rationalize[#]&,amplitudes[[All,"CG"]]]; (* CG coeffs from symmetrization *)
amplitudes
];


Options[generateASRs]={particles->{}};
generateASRs[system_,OptionsPattern[]]:=Module[{systemVal=system,partVal=OptionValue[particles],ASRs,amplitudes,dupPairs,factorsMat},
{systemVal,ASRs,dupPairs}=pyEval["generate_srs",systemVal]; (* M values for symmetrized system, contains duplicate amplitudes *)
amplitudes=extractAmps[systemVal];

factorsMat=DiagonalMatrix[amplitudes[[All,"mu"]]*amplitudes[[All,"CG"]]];
ASRs=Map[# . factorsMat&,ASRs]; (* SRs for symmetrized system *)

(* Correct for duplicate amplitudes from symmetrization *)
Do[ASRs[[b]][[All,dupPairs[[All,1]]]]+=ASRs[[b]][[All,dupPairs[[All,2]]]],{b,Length[ASRs]}]; (* adds together cols of duplicate amp pairs *)
ASRs=Transpose[Delete[Transpose[#],List/@dupPairs[[All,2]]]]&/@ASRs; (* deletes duplicate cols *)
systemVal=pyEval["remove_dups",{systemVal,dupPairs}];
amplitudes=extractAmps[systemVal,particles->partVal];

ASRs=Map[Cases[Except@{0..}],Map[RowReduce,ASRs]];

{systemVal,ASRs,amplitudes}
];


numAmps[amplitudes_]:=Length@DeleteDuplicates@Flatten@amplitudes[[All,"Indices"]];


printAmps[amplitudes_]:=Module[{amplitudesVal=amplitudes,indices,selfConjs,nAmps},
If[amplitudesVal[[1,"Processes"]]=!="n/a",
(indices=amplitudesVal[[All,"Indices"]];
selfConjs=Table[If[indices[[i,1]]==indices[[i,2]],i,Nothing],{i,Length[indices]}];

amplitudesVal[[selfConjs,"Processes"]]=DeleteDuplicates/@amplitudesVal[[selfConjs,"Processes"]];
amplitudesVal[[selfConjs,"CKM"]]=DeleteDuplicates/@amplitudesVal[[selfConjs,"CKM"]];
),
Null
];

amplitudesVal[[All,"Indices"]]=DeleteDuplicates/@Table[amplitudesVal[[i,"Indices"]],{i,Length[amplitudesVal]}];
nAmps=numAmps[amplitudesVal];

Print["Amplitude table","\n",
"Number of amplitudes: ",nAmps,"\n",
TableForm[Prepend[Values/@amplitudesVal,Keys[amplitudesVal[[1]]]]]
];

{nAmps,amplitudesVal}
];


numASRs[ASRs_]:=Table[Length[ASRs[[i]]],{i,Length[ASRs]}];


Options[printASRs]={takeProd->False,ampFormat->"a/s n-tuple",CKM->False,b->All};
printASRs[ASRs_,amplitudes_,OptionsPattern[]]:=Module[{takeProd=OptionValue[takeProd],ampFormat=OptionValue[ampFormat],CKM=OptionValue[CKM],b=OptionValue[b],physicalAmps,ampsToVector,ampVector,writtenASRs,blist,numASRsList,matForm,rule},
physicalAmps[amp_,asType_,col_,CKM_]:=Module[{asSign,CKMfactors},
asSign=If[SymbolName[asType]=="a",-1,+1];
CKMfactors=If[CKM,{amp[["CKM",1]],amp[["CKM",2]]},{1,1}];
amp[["q factor"]](A[amp[[col,1]]]/CKMfactors[[1]]+asSign*amp[["p factor"]]*A[amp[[col,2]]]/CKMfactors[[2]])
];

ampsToVector[amps_,asType_,ampFormat_,CKM_]:=Module[{vector},
vector=Switch[ampFormat,
"a/s indices",Map[asType[#]&,amps[[All,"Indices",1]]], (* a/s-type amps with numbered subscripts *)
"a/s nodes",Map[asType[#]&,amps[[All,"Node"]]], (* a/s-type amps with nodes *)
"A indices",Map[physicalAmps[#,asType,"Indices",CKM]&,amps], (* A amps with numbered subscripts *)
"A physical",Map[physicalAmps[#,asType,"Processes",CKM]&,amps], (* A amps with physical processes *)
_,Map[asType[#]&,amps[[All,"n-tuple"]]] (* a/s-type amps with n-tuples; default option *)
];
vector
];

ampVector={ampsToVector[amplitudes,a,ampFormat,CKM],ampsToVector[amplitudes,s,ampFormat,CKM]};

writtenASRs=If[takeProd,
MapIndexed[If[OddQ[First[#2]],#1 . ampVector[[1]],#1 . ampVector[[2]]]&,ASRs],
MapIndexed[If[OddQ[First[#2]],Prepend[#1,ampVector[[1]]],Prepend[#1,ampVector[[2]]]]&,ASRs]
];

b=b/.{
{i_,j_,k_:1}:>Span[
Replace[i,{Null->1,x_Integer:>x+1}],
Replace[j,{Null->All,x_Integer:>x+1}],
k
],
{i_}:>Span[Replace[i,{Null->1,x_Integer:>x+1}],All],
i_Integer:>Span[i+1,i+1],
All:>All
};
blist=Flatten@{(Range[Length[ASRs]]-1)[[b]]};
numASRsList=numASRs[ASRs];
matForm:=If[takeProd,MatrixForm[#]&,MatrixForm[#[[2;;]],TableHeadings->{None,#[[1]]}]&];
rule=x_Symbol[i_]/;x=!=List:>Subscript[x,i];

Print["Sum rules"];
MapIndexed[
Print["b = ",blist[[#2[[1]]]],"\n",
"Number of ASRs: ",numASRsList[[blist[[#2[[1]]]]+1]],"\n",
matForm@If[takeProd,#1/.rule,ReplacePart[#1,1->(#1[[1]]/.rule)]]
]&,
writtenASRs[[b]]
];

{numASRsList,If[takeProd,writtenASRs,ASRs]}
];


Options[printSystem]=Options[printASRs];
printSystem[system_,ASRs_,amplitudes_,opts:OptionsPattern[]]:=Module[{n,irreps,nAmps,nASRs,writtenASRs},
{n,irreps}=pyEval["extract_sys",system];
Print["System","\n",
"Number of doublets: ",n,"\n",
"System irreps: ",Sort@Flatten@irreps,"\n",
"In: ",irreps[[1]],"\n",
"H: ",irreps[[2]],"\n",
"Out: ",irreps[[3]]
];
Print["-------------------------"];
nAmps=printAmps[amplitudes][[1]];
Print["-------------------------"];
{nASRs,writtenASRs}=printASRs[ASRs,amplitudes,Sequence@@FilterRules[{opts},Options[printASRs]]];

<|"n doublets"->n,"Irreps"->irreps,"n amps"->nAmps,"Amplitudes"->amplitudes,"n ASRs"->nASRs,"ASRs"->writtenASRs|>
];


End[]


EndPackage[]
