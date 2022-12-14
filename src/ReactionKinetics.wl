(* ::Package:: *)

(* :Title:   ReactionKinetics  *)
(* :Context: ReactionKinetics` *)
(* :Version: 1.0               *)

(* :Authors:    Janos TOTH,     Attila Laszlo NAGY,         David PAPP          *)
(* :E-mails: jtoth@math.bme.hu, nagyal@math.bme.hu, dpapp@iems.northwestern.edu *)

(* :References: Janos TOTH, Attila Laszlo NAGY and David PAPP. 
Reaction Kinetics: Exercises, Programs and Theorems. Mathematica for deterministic and stochastic kinetics, Springer, 2018. *)

(* :Summary:  This is a Wolfram Mathematica package intended to facilitate the work of researchers in reaction kinetics. *)

(* :Keywords: chemical reactions, stoichiometry, deterministic models, induced kinetic differential equation, 
              mass action type kinetics, general kinetics, stochastic models, induced kinetic Markov process,
			  stochastic simulation, graphs of reactions, decompositions, Arrhenius equation, CHEMKIN standards. *)

(* :Acknowledgements: Daniel Lichtbau, Balazs Boros, Bela Vizvari, Yifan Hu, Benedek Kovacs, Tamas Ladics, Ilona Nagy and others. *)

(* :History: 
    Version 0.1 by David PAPP, 2005. 
    Version 0.3.2 by David PAPP, 2006. *)


(* ::Subsection::Closed:: *)
(*license*)


(*
	ReactionKinetics.m is a Mathematica package aimed to be a tool for  
	analyzing (bio)chemical reactions and models.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*)


(*
	''You acknowledge that Software is not designed, licensed or intended for
	use in the design, construction, operation or maintenance of any nuclear
	facility.'' (Java SDK demo applet)
*)


(* ::Subsection::Closed:: *)
(*begin*)


BeginPackage["ReactionKinetics`"]
Unprotect["ReactionKinetics`*"];
ClearAll["ReactionKinetics`*"];
ClearAll["ReactionKinetics`Private`*"];


If[$VersionNumber<6, 
	Off[DiscreteMath`Combinatorica`Star]; 
	Needs["Graphics`"];
	Needs["Graphics`Legend`"];
	Needs["LinearAlgebra`MatrixManipulation`"];
	Needs["Utilities`"];
	On[DiscreteMath`Combinatorica`Star];
]; 

(* If[$VersionNumber >= 6,
	Needs["XML`"];
	Needs["JLink`"];
	Needs["GUIKit`"];
]; *)
(* XML` is automatically loaded in versions 4 & 5 *)

(*
If[$VersionNumber>5,
	Needs["DifferentialEquations`NDSolveProblems`"];
	Needs["DifferentialEquations`NDSolveUtilities`"];
	Needs["FunctionApproximations`"];
	Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];
	Needs["PlotLegends`"];
	Needs["MathSBML`"];
	Needs["EquationTrekker`"];
];
*)

If[$VersionNumber < 8,
	Off[General::compat];
	Needs["Combinatorica`"];
	Needs["GraphUtilities`"];
    On[General::compat];
]; 


Unprotect[
	$ReactionKinetics,
	$ReactionKineticsPackageLoaded,
	$ReactionKineticsVersionNumber
];


(* ================================================================ *)
(* usage messages for the exported functions and the context itself *)
(* ================================================================ *)


$ReactionKinetics::usage = "ReactionKinetics.m is a Mathematica package that is no doubt the nicest package ever created to help you do fine reaction kinetics.";

$ReactionKineticsVersionNumber::usage = "$ReactionKineticsVersion returns the version number and release date of the package which is currently being used.";

$ReactionKineticsPackageLoaded::usage = "$ReactionKineticsPackageLoaded returns true if the current version of ReactionKinetics is loaded.";


$ReactionKineticsVersionNumber = "1.0 [March 25, 2018]";


(* ::Subsection::Closed:: *)
(*main*)


(* =============================================== *)
(* begin the private context (implementation part) *)
(* =============================================== *)

Begin["`Private`"]


(* ::Subsubsection::Closed:: *)
(*Options*)


Bipartite::usage = "An option for function ShowVolpertGraph. Bipartite \[Rule] True arranges the vertices of the (bipartite) Volpert graph \
into two disjoint sets (of reaction steps and species) higlighting the edges connecting them.";
CombinatorialMoments::usage = "An option for function MomentEquations. If CombinatorialMoments \[Rule] True, then combinatorial moments are taken into account.";
ComplexColors::usage = "An option for function ShowFHJGraph. E.g. ComplexColors \[Rule] listcolors colors the complexes (i.e. edges) of the FHJ graph according to listcolors.";
Conditions::usage = "An option for function StationaryPoints. E.g. Conditions \[Rule] list adds the list to the constraints stationary points should obey.";
ContejeanDevie::usage = "Method \[Rule] ContejeanDevie is an option for Decompositions for decomposing an overall reaction.";
EdgeLabels::usage = "An option for function ShowVolpertGraph. If EdgeLabels \[Rule] True, all the stoichiometric coefficients are put \
on the directed edges of the Volpert graph.";
ExternalSpecies::usage = "An option for several functions including e.g. ReactionsData, ShowFHJGraph, Concentrations, Simulation, Decompositions etc. \
The list of external species of a reaction can be given by ExternalSpecies \[Rule] listexternals.";
FastSelection::usage = "ObjectiveFunction \[Rule] FastSelection is an option for \
CoveringDecompositionSet to specify how the consecutive decompositions should be searched for.";
Filter::usage = "An option for function Decompositions. E.g. Filter \[Rule] MinimalDecompositions filters out minimal decompositions.";
FixedStepSize::usage = "An option for function Simulation (see approximation methods).";
FormattedOutput::usage = "An option for function ToAtomMatrix. FormattedOutput \[Rule] True pretty prints a table for atomic matrix.";
GeneratedRateCoefficient::usage = "An option for function DetailedBalanced. GeneratedRateCoefficient \[Rule] rratecoeffs uses rratecoeffs as reaction rate coefficients.";
GreedySelection::usage = "ObjectiveFunction \[Rule] GreedySelection is an option for \
CoveringDecompositionSet to specify how the consecutive decompositions should be searched for.";
Highlight::usage = "An option for function ShowVolpertGraph. Highlight \[Rule] listvertices highlists the vertices or the subgraph induced by vertices listed in listvertices \
depending on SubgraphHighlight \[Rule] False or True.";
Highlighted::usage = "An option for function ShowVolpertGraph. Highlighted \[Rule] listvertices highlists the vertices or the subgraph induced by vertices listed in listvertices.";
Indexed::usage = "An option for function ShowVolpertGraph. Indexed \[Rule] listspecies displays the Volpert graph with Volpert indexes, \
where listspecies is considered as the initial set of species.";
LPBased::usage = "Method \[Rule] LPBased is an option for Decompositions using linear programming to obtain decompositions for an overall reaction.";
MassAction::usage = "An option for functions RightHandSide, DeterministicModel, Concentrations, StationaryPoints and EigensystemJacobian. When \
MassAction \[Rule] True, mass action type kinetics are taken into account with reaction rate coefficients.";
MaxIteration::usage = "An option for function Simulation. It specifies the maximum number of iterations that should be used during the simulation.";
Memo::usage = "An option for function CHEMKINImport. Memo \[Rule] True forces CHEMKINImport to use memorization i.e. to remember outputs of already given inputs.";
Method::usage = "Method \[Rule] method is an option for several functions letting method be used during the calculation.";
MinimalDecompositions::usage = "An option for function Decompositions, e.g. Filter \[Rule] MinimalDecompositions.";
Numbered::usage = "An option for functions ShowFHJGraph and ShowVolpertGraph. \
Numbered \[Rule] True only shows the graphs with numbers assigned to each reaction step.";
ObjectiveFunction::usage = "An option for function Decompositions, e.g. ObjectiveFunction \[Rule] OriginalSelection.";
OriginalSelection::usage = "ObjectiveFunction \[Rule] OriginalSelection is an option for \
CoveringDecompositionSet to specify how the consecutive decompositions should be searched for. \
(The naming is purely historical and has no mathematical justification.)";
PlotFunction::usage = "An option for functions ShowFHJGraph, ShowVolpertGraph, SimulationPlot, SimulationPlot2D and SimulationPlot3D. \
E.g. PlotFunction \[Rule] \"GraphPlot\" uses \"GraphPlot\" to visualize the results.";
Positivity::usage = "An option for function StationaryPoints. Positivity \[Rule] True tries to find positive stationary points.";
Preprocess::usage = "An option for function Decompositions. Preprocess \[Rule] True does preprocessing.";
Side::usage = "An option for function FilterReactions. Side \[Rule] \"Product\" (or \"Reactant\" or \"All\") filters all the reaction steps, \
where one of the reactant or product or (reactant or product) species is contained in a given list of species.";
Species::usage = "An option for functions SimulationPlot, SimulationPlot2D and SimulationPlot3D. \
Species \[Rule] listspecies only plots those species given by listspecies.";
StronglyConnectedComponentsColors::usage = "An option for function ShowFHJGraph. StronglyConnectedComponentsColors \[Rule] listcolors colors the strongly connected \
components of the FHJ graph according to listcolors.";
TimeLimit::usage = "An option for function AbsoluteConcentrationRobustness. \
TimeLimit \[Rule] timelimit suppresses calculations if its duration exceeds timelimit.";
Tolerance::usage = "An option for function Simulation (see approximation methods).";
Verbose::usage = "An option for several functions. Verbose \[Rule] True (default) gives more detailed information about the calculation.";
Volume::usage = "An option for function Simulation. Volume \[Rule] vol simulates the process in volume vol.";


Global`aExplodator::usage = "Coefficient in the Explodator reaction.";
Global`bExplodator::usage = "Coefficient in the Explodator reaction.";
Global`fOregonator::usage = "Coefficient in the Oregonator reaction.";
Global`t::usage = "Default symbol for time.";
Global`k::usage = "Default symbol for the reaction rate coefficients.";
Global`c::usage = "Default symbol for concentrations in the induced kinetic differential equation.";
Global`\[DoubleStruckCapitalE]::usage = "Default symbol for the expected value.";
Global`g::usage = "Default symbol for the probability generating function.";
Global`\[CapitalPi]::usage = "Default symbol for stationary distribution.";
Global`P::usage = "Default symbol in the master equation for probability.";
Global`z::usage = "Default symbol for the variable of the probability generating function.";


(* ::Subsubsection::Closed:: *)
(*Constants*)


AvogadrosNumber::usage = "AvogadrosNumber (1/Mole) returns the Avogadro constant that is the number of constituent particles contained \
in the amount of substance given by one mole.";

AvogadrosNumber = QuantityMagnitude[UnitConvert[Quantity["AvogadroNumber"]]]; (*1 / Mole*)

MolarGasConstant::usage = "MolarGasConstant (Joule/(Kelvin*Mole)) returns the ideal gas constant that is the constant of proportionality \
expressing the relation between the energy scale and the temperature scale in physics.";

MolarGasConstant = QuantityMagnitude[UnitConvert[Quantity["MolarGasConstant"]]]; (*Joule / Kelvin / Mole*)


(* ::Subsubsection::Closed:: *)
(*Auxiliary functions*)


ZeroVectorQ::usage = "ZeroVectorQ[list] checks whether list has only zero elements.";
ZeroVectorQ[v_List] := Max[Abs[v]] === 0;


VariablesQ[x_] := And@@((Head[#]===Symbol)&/@x);

MyToRules[list_] := Flatten[If[Head[#]===Equal,ToRules[#],#]&/@list];

MyFilterOptions[f_, opts___]:= 
(*	If[$VersionNumber<6,
			Utilities`FilterOptions`FilterOptions[f, opts],*)
			Sequence @@ FilterRules[Flatten[{opts}], Options[f]];

SetAttributes[{PositivePart,NegativePart}, Listable];
PositivePart[x_] := Max[x,0];
NegativePart[x_] := -Min[x,0];

If[$VersionNumber < 10,
	(*SubsetQ[x_,y_] := (Intersection[x,y] == x);*)
	SubsetQ[x_,y_] := And @@ (MemberQ[x,#] &/@ y);
];

MinimalQ[v_?VectorQ, vlist_List] :=
	If[Scan[If[VectorLessEqual[#, v], Return[False]] &, vlist] === False, False, True];

SMaux[minimals_List, {}] := minimals;
SMaux[minimals_List, {next_List, rest___List}] :=
	SMaux[Join[minimals, Select[next, MinimalQ[#, minimals] &]], {rest}];

MyFloor[x_, prec_] := Floor[x/prec]*prec;
MyCeiling[x_, prec_] := Ceiling[x/prec]*prec;

MyMinimize[var_, conlist_, vars_] :=
	Module[{c, b, A, lp},
		(* EZ LENNE A NORMALIS MEGOLDAS *)
		(* First[Minimize[var, Thread[conlist >= 0], vars]] *)
		(* c = Normal[Last[CoefficientArrays[var, vars]]]; *)
		c = Coefficient[var, #] & /@ vars;
		{b, A} = Normal[CoefficientArrays[conlist, vars]];
		lp = LinearProgramming[c, A, -b, -Infinity, Method -> "Simplex"]; (* Automatic, Simplex, RevisedSimplex, InteriorPoint *)
		If[Head[lp] === LinearProgramming, "no solution", c . lp]
	];

MyMaximize[var_, conlist_, vars_] := -MyMinimize[-var, conlist, vars];


(* ::Subsection::Closed:: *)
(*models*)


Models::usage = "Models contains a list of names of some of the widespread reactions in reaction kinetics.";

Models := 
{
	"Autocatalator",
	"Belousov-Zhabotinsky",
	"Briggs-Rauscher",
	"Brusselator",
	"Chapman cycle",
	"Clarke",
	"Colquhoun",
	"Consecutive",
	"Decomposition",
	"DeYoung-Keizer",
	"Dimerization",
	"Edelstein",
	"Eigen",
	"\[CapitalEAcute]rdi-Ropolyi",
	"Explodator",
	"Horn-Jackson",
	"Feinberg-Horn I",
	"Feinberg-Horn II",
	"FKN mechanism",
	"Frank",
	"Goldbeter",
	"Hudson-R\[ODoubleDot]ssler",
	"Huxel",
	"Hyver",
	"Inflow",
	"Ivanova",
	"Kliemann",
	"Leonard-Reichl",
	"Lotka",
	"Lotka-Volterra",
	"Michaelis-Menten",
	"Mole",
	"Ogg",
	"Oregonator",
	"Outflow",
	"Petri",
	"Robertson",
	"Schl\[ODoubleDot]gl",
	"Simple birth-death",
	"Triangle",
	"Tur\[AAcute]nyi-Gy\[ODoubleDot]rgyi-Field",
	"Vaiman",
	"Verhulst",
	"Volpert",
	"Wegscheider",
	"Willamowski-R\[ODoubleDot]ssler"
};


Reactions::usage = "Reactions is a list of some of the widespread reactions in reaction kinetics.";

Reactions :=
Thread[
	Rule[Models,
		{
	(*Autocatalator*) 
{"A"->"X","X"->"Y","X"+2"Y"->3"Y","Y"->"P"}, 

	(*Belousov-Zhabotinsky*) 
{"X"+"Y"+"H"->2"V",
"Y"+"A"+2"H"->"X"+"V",
2"X"->"V",
1/2"X"+"A"+"H"->"X"+"Z",
"X"+"Z"->1/2"X",
"Z"+"M"->"Q",
"Z"+"V"->"Y",
"V"->"Y",
"X"->"0",
"Y"->"0",
"Z"->"0",
"V"->"0"}, 

	(*Briggs-Rauscher*) 
{"IO3 -"+2"H2O2"+"CH2(CO2H)2"+"H +"->"ICH(CO2H)2"+2"O2"+3"H2O",
"IO3 -"+2"H2O2"+"H +"->"HIO"+2"O2"+2"H2O",
"HIO"+"CH2(CO2H)2"->"ICH(CO2H)2"+"H2O",
"IO3 -"+"I -"+2"H +"->"HIO2"+"HIO",
"HIO2"+"I -"+"H +"->2"HIO",
"HIO"+"H2O2"->"I -"+"O2"+"H +"+"H2O",
"IO3 -"+"HIO2"+"H +"->2"IO2"+"H2O",
"IO2"+"Mn 2+"+"H2O"->"HIO2"+"Mn(OH) 2+",
"Mn 2+"+"H2O2"->"Mn 2+"+"H2O"+"HO2",
2"HO2"->"H2O2"+"O2",
2"HIO2"->"IO3 -"+"HIO"+"H +",
"I -"+"HIO"+"H +"->"I2"+"H2O",
"I2"+"CH2(CO2H)2"->"ICH(CO2H)2"+"H +"+"I -"},

	(*Brusselator*) 
{"X"\[LeftArrow]"A","B"+"X"->"Y"+"D","Y"+2"X"->3"X","X"->"E"}, 
(*externals: "A","B","D","E"*)

	(*Chapman cycle*)
{"O2"->"O"+"O", "O2"+"O"+"M"->"O3"+"M", "O3"->"O2"+"O", "O"+"O3"->"O2"+"O2"},
(*externals: "M"*)

	(*Clarke*)
{"A"+"T"->"M"+"S",
"T"+"Y"+"H"->"V"+"S",
"X"+"Y"+"H"->2"T",
"B"+"Y"+2"H"->"X"+"T",
2"X"->"B"+"T"+"H",
"B"+"X"+"H"->2"U"+"S",
"U"+"W"+"H"->"X"+"Z",
"U"+"Z"+"S"->"B"+"W"+2"H",
"A"+"V"->"M"+"Y"+"H",
"A"+6"Z"+2"S"->"P"+2"Q"+6"W"+6"H",
"M"+4"Z"+2"S"->"P"+2"Q"+4"W"+"Y"+5"H",
"P"+"T"->"Q"+"Y"+"H"+"S"},

	(*Colquhoun*)
{"AR"\[Equilibrium]"A"+"R"\[Equilibrium]"RA","A"+"AR"\[Equilibrium]"ARA"\[Equilibrium]"A"+"RA"},

	(*Consecutive*)
{"A"->"B"->"C"},

	(*Decomposition*)
{"C" -> "A"+"B"},

	(*DeYoung-Keizer*)
{"S(000)"\[Equilibrium]"S(001)"\[Equilibrium]"S(100)"\[Equilibrium]"S(010)","S(100)"\[Equilibrium]"S(101)"\[Equilibrium]"S(110)","S(001)"\[Equilibrium]"S(101)"\[Equilibrium]"S(011)",
"S(011)"\[Equilibrium]"S(010)"\[Equilibrium]"S(111)","S(110)"\[Equilibrium]"S(010)"\[Equilibrium]"S(111)","S(111)"\[Equilibrium]"S(101)"},

	(*Dimerization*)
{"A"+"B" -> "C"},

	(*Edelstein*)
{"X"\[DoubleLeftRightArrow]2"X","X"+"Y"\[DoubleLeftRightArrow]"Z"\[DoubleLeftRightArrow]"Y"},

	(*Eigen*) 
{"X"->0,"X"->2"X"},

	(*\[CapitalEAcute]rdi-Ropolyi*)
{"R"+4"T"\[Equilibrium]"T4R1"\[Equilibrium]"T4R2"\[Equilibrium]"T4R3"},
	
	(*Explodator*) 
{"A"+"X"->(1+Global`aExplodator)"X","X"+"Y"->"Z","Z"->(1+Global`bExplodator)"Y","Y"->"P"},

	(*Horn-Jackson*)
{3"X"->"X"+2"Y"->3"Y"->2"X"+"Y"->3"X"},

	(*Feinberg-Horn I*)
{2"J" -> "G", 2"J" \[LeftRightArrow] "H", "G"->"H", "A"+"B" -> "G", "A"+"B" \[LeftRightArrow] "C", "C" \[RightArrow] "D"+"E", "D"+"E" \[LeftRightArrow] "F"},

	(*Feinberg-Horn II*) 
{"X"->"Y","X"+"Z"->"T"->"Y"+"U"->"X"+"Z"},

	(*FKN mechanism*)
{"HBrO2"+"Br -"->2"HOBr","Br -"+"BrO3 -"->"HBrO2"+"HOBr",2"HBrO2"->"BrO3 -"+"HOBr","HOBr"+"Br -"->"Br2"+"H2O","BrO3 -"+"HBrO2"->2"BrO2"+"H2O",
"HBrO2"+"H2O"->"BrO3 -","Br2"+"MA"->"BrMA"+"Br -","MA"+2"H2O"->C+2"CO2","BrMA"+2"H2O"->"Br -"+"HCO2H"+2"CO2"},

	(*Frank*)
{"A" -> "R", "A" -> "S", "A"+"R" -> 2"R", "A"+"S" -> 2"S"},

	(*Goldbeter*)
{"M"\[Equilibrium]"0"\[Equilibrium]"C","0"\[Equilibrium]"X","C"\[RightArrow]"0"},

	(*Hudson-R\[ODoubleDot]ssler*) 
{"P"->"A","Q"->"B","A"+"B"->2"B","B"->"R","A"\[LeftRightArrow]"C"}, 
(*externals: "P","Q","R"*)
	
	(*Huxel (special Lotka-Volterra)*)
{"X"\[Equilibrium]2"X", "X"+"Y"->2"Y","Y"->"0"},

	(*Hyver*) 
{"K1"->"A1"->"A2"->"A3"->"A4"->"A5"->"A6"->"A7", "A7"+"B"->"C"+"E",
"K"->"B"->"0", "A1"+"E"->"0", "A1"+"C"->"0", "\[Mu]"->"D", "D"+"C"->"0", "D"+"E"->"0"},

	(*Inflow*)
{"0"->"X"},

	(*Ivanova*)
{"X"+"Y"->2"Y","Y"+"Z"->2"Z","X"+"Z"->2"X"},

	(*Kliemann*)
{"X"\[Equilibrium]"Y",2"X"\[Equilibrium]"X"+"Y"\[Equilibrium]2"Y"},

	(*Leonard-Reichl*) 
{2"X"\[LeftArrowRightArrow]"Y"},

	(*Lotka*) 
{"A"->"X","X"+"Y"->2"Y","Y"->"P"},

	(*Lotka-Volterra*)
{"A"+"X"->2"X","X"+"Y"->2"Y","B"\[LeftArrow]"Y"},

	(*Michaelis-Menten*)
{"E"+"S"\[Equilibrium]"C"->"E"+"P"},

	(*Mole*)
{"X" + "Y" \[RightArrow] 2 "X" + 2 "Y", "X" \[LeftRightArrow] "0" \[LeftRightArrow] "Y"},

	(*Ogg*)
{"\!\(\*SubscriptBox[\(N\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \(5\)]\)" \[LeftRightArrow] "\!\(\*SubscriptBox[\(NO\), \(2\)]\)"+"\!\(\*SubscriptBox[\(NO\), \(3\)]\)", "\!\(\*SubscriptBox[\(NO\), \(2\)]\)"+"\!\(\*SubscriptBox[\(NO\), \(3\)]\)"->"\!\(\*SubscriptBox[\(NO\), \(2\)]\)"+"NO"+"\!\(\*SubscriptBox[\(O\), \(2\)]\)", "\!\(\*SubscriptBox[\(NO\), \(3\)]\)"+"NO"->2"\!\(\*SubscriptBox[\(NO\), \(2\)]\)"},

	(*Oregonator*) 
{"A"+"Y"->"X","X"+"Y"->"P","B"+"X"->2"X"+"Z",2"X"->"Q","Z"->Global`fOregonator*"Y"},

	(*Outflow*)
{"X"->"0"},

	(*Petri*)
{"C"+"\!\(\*SubscriptBox[\(O\), \(2\)]\)" -> "\!\(\*SubscriptBox[\(CO\), \(2\)]\)", "\!\(\*SubscriptBox[\(CO\), \(2\)]\)"+"NaOH" -> "\!\(\*SubscriptBox[\(NaHCO\), \(3\)]\)", "\!\(\*SubscriptBox[\(NaHCO\), \(3\)]\)"+"HCl" -> "\!\(\*SubscriptBox[\(H\), \(2\)]\)O"+"NaCl"+"\!\(\*SubscriptBox[\(CO\), \(2\)]\)"},

	(*Robertson*)
{"A"->"B", 2"B"->"C"+"B", "B"+"C"->"A"+"C"},

	(*Schl\[ODoubleDot]gl*) 
{0\[Equilibrium]"X",2"X"\[Equilibrium]3"X"},

	(* Simple birth-death *)
{"X"->2"X", "X"->"A"},

	(*Triangle*)
{"A"->"B"\[RightArrow]"C"->"A"},

	(*Tur\[AAcute]nyi-Gy\[ODoubleDot]rgyi-Field*)
{"X"+"Y" -> 2"P", "Y"+"A" -> "X"+"P", 2"X" -> "P"+"A", "X"+"A" -> 2"X"+2"Z", "X"+"Z" -> 0.5"X"+"A", "Z"+"M"->"Y"-"Z"},

	(*Vaiman*)
{"O2" + 2 "K" -> 2 "KO", "H2S" + "K" -> "KH2S", "O2" + 2 "KH2S" -> 2 "H2O" + "K" + "KS2", 
    "H2S" + "KO" -> "H2O" + "KS", "KO" + "KH2S" -> "H2O" + "KS" + "K", 2 "KS" -> "K" + "KS2", "KS2" -> "S2" + "K"},

	(*Verhulst*)
{"X"\[Equilibrium]2"X"},

	(*Volpert*)
{"X"+"Y"->"T","Y"+"Z"->"U"},

	(*Wegscheider*)
{"A"\[Equilibrium]"B",2"A"\[Equilibrium]"A"+"B"},

	(*Willamowski-R\[ODoubleDot]ssler*)
{"A"+"X"\[ReverseEquilibrium]2"X","X"+"Y"\[ReverseEquilibrium]2"Y","Q"+"Y"\[ReverseEquilibrium]"B","X"+"Z"\[ReverseEquilibrium]"C","P"+"Z"\[ReverseEquilibrium]2"Z"}
		}
		]
];


BioModels::usage = "BioModels contains a list of names of some of the widespread biochemical reactions.";

BioModels := 
{
	"Calvin Cycle (Photosynthesis)",
	"Glycolysis",
	"Glyoxylate Cycle",
	"SzentGy\[ODoubleDot]rgyi-Krebs Cycle"
};


BioReactions::usage = "BioReactions is a list of some of the widespread biochemical reactions.";

BioReactions=
Thread[
	Rule[BioModels,
		{
	(* Calvin Cycle (Photosynthesis) *)
{"DG3P"->"GP","DF1,6BP"->"GP"+"DG3P","DF1,6BP"+"H2O"->"DF6P"+"Pi","S7P"+"DG3P"->"DR5P"+"DX5P","DR5P"->"DRL5P","DRL5P"->"DX5P",
"ATP"+"DRL5P"->"ADP"+"DRL1,5BP","3PDG"+2"H+"->"DRL1,5BP"+"CO2"+"H2O" ,"ATP"+"3PDG"->"ADP"+"3PDGP","DG3P"+"Pi"+"NADP+"->"3PDGP"+ "NADPH" + "H+"}/.
	{"DG3P"->Tooltip["DG3P","D-glyceraldehyde 3-phosphate"],
	"GP"->Tooltip["GP","glycerone phosphate"],
	"DRL5P"->Tooltip["DRL5P","D-ribulose 5-phosphate"],
	"DR5P"->Tooltip["DR5P","D-ribose 5-phosphate"],
	"DX5P"->Tooltip["DX5P","D-xylulose 5-phosphate"],
	"S7P"->Tooltip["S7P","sedoheptulose 7-phosphate"],
	"DF6P"->Tooltip["DF6P","D-fructose 1,6-bisphosphate"],
	"DG3P"->Tooltip["DG3P","D-glyceraldehyde 3-phosphate"],
	"DF1,6BP"->Tooltip["DF1,6BP","D-fructose 1,6-bisphosphate"],
	"3PDG"->Tooltip["3PDG","3-phospho-D-glycerate"],
	"Pi"->Tooltip["Pi","phosphate"],
	"DRL1,5BP"->Tooltip["DRL1,5BP","D-ribulose 1,5-bisphosphate"],
	"3PDGP"->Tooltip["3PDGP","3-phospho-D-glyceroyl phosphate"],
	"H2O"->Tooltip["H2O","dihydrogen monoxide"],
	"H+"->Tooltip["H+","hydrogen ion"],
	"CO2"->Tooltip["CO2","carbon dioxide"],
	"ATP"->Tooltip["ATP","adenosine-5'-triphosphate"],
	"ADP"->Tooltip["ADP","adenosine diphosphate"],
	"NADP+"->Tooltip["NADP+","nicotinamide adenine dinucleotide phosphate"],
	"NADPH"->Tooltip["NADPH","NADPH"](*???*)},

	(* Glycolysis *)
{"Glc"+"ATP"-> "G6P"+"ADP"+"H+","G6P"\[Equilibrium]"F6P","F6P"+"ATP"->"ADP"+"F1,6BP","F1,6BP"\[RightArrowLeftArrow]"GADP"+"DHAP","DHAP"\[Equilibrium]"GADP",
"GADP"+"NAD+"+"Pi"-> "NADH"+"H+"+"1,3BPG","1,3BPG"->"GADP","1,3BPG"+"ADP"->"ATP"+"3PG","3PG"\[ReverseEquilibrium]"2PG",
"2PG"->"PEP"+"H2O","PEP"->"2PG","PEP"+"ADP"-> "ATP"+"Pyr"}/.
	{"Glc"->Tooltip["Glc","D-glucose"],
	"G6P"->Tooltip["G6P","\[Alpha]-D-glucose-6-phosphate"],
	"F6P"->Tooltip["F6P","\[Beta]-D-fructose 6-phosphate"],
	"Pi"->Tooltip["Pi","phosphate"],
	"F1,6BP"->Tooltip["F1,6BP","\[Beta]-D-fructose 1,6-bisphosphate"],
	"GADP"->Tooltip["GADP","D-glyceraldehyde 3-phosphate"],
	"DHAP"->Tooltip["DHAP","dihydroxyacetone phosphate"],
	"1,3BPG"->Tooltip["1,3BPG","D-1,3-bisphosphoglycerate"],
	"3PG"->Tooltip["3PG","3-phosphoglycerate"],
	"2PG"->Tooltip["2PG","2-phosphoglycerate"],
	"PEP"->Tooltip["PEP","phosphoenolpyruvate"],
	"Pyr"->Tooltip["Pyr","pyruvate"],
	"H+"->Tooltip["H+","hydrogen ion"],
	"H2O"->Tooltip["H2O","dihydrogen monoxide"],
	"ATP"->Tooltip["ATP","adenosine-5'-triphosphate"],
	"ADP"->Tooltip["ADP","adenosine diphosphate"],
	"NAD+"->Tooltip["NAD+","nicotinamide adenine dinucleotide"],
	"NADH"->Tooltip["NADH","NADH"](*???*)},

	(* Glyoxylate Cycle *)
{"SM"+"NAD+"->"Oxc" + "NADH" + "H+","ACoA"+"H2O"+"Oxc"->"C"+"HS-CoA"+"H+","ACoA"+"H2O"+"G"->"SM"+"HS-CoA"+"H+","IsoC"->"S"+"G","C"->"cA"+"H2O","cA"+"H2O"->"IsoC"}/.
	{"Oxc"-> Tooltip["Oxc","oxalocetate"],
	"C"->Tooltip["C","citrate"],
	"ACoA"-> Tooltip["ACoA","acetyl CoA"],
	"SM"->Tooltip["SM","S-malate"],
	"S"-> Tooltip["S","succinate"],
	"IsoC"->Tooltip["IsoC","isocitrate"],
	"cA"->Tooltip["cA","cis-aconitate"],
	"G"->Tooltip["G","glyoxylate"],
	"H+"->Tooltip["H+","hydrogen ion"],
	"H2O"->Tooltip["H2O","dihydrogen monoxide"],
	"NAD+"->Tooltip["NAD+","nicotinamide adenine dinucleotide"],
	"HS-CoA"->Tooltip["HS-CoA","coenzyme A"],
	"NADH"->Tooltip["NADH","NADH"](*???*)},

	(* SzentGy\[ODoubleDot]rgyi-Krebs Cycle *)
{"Oc"+"ACoA"+"H2O"->"C"+"CSH","C"->"cA"+"H2O"->"IsoC","IsoC"+"NAD+"-> "Os"+"NADH"+"H+","Os"-> "\[Alpha]K"+"CO2","\[Alpha]K"+"NAD+"+"CSH"->"SCoA"+"NADH"+"H+"+"CO2",
"SCoA"+"GDP"+"Pi"->"S"+"CSH"+"GTP","S"+"Q"->"F"+"QH2","F"+"H2O"->"LM","LM"+"NAD+"->"Oc"+"NADH"+"H+"}/.
	{"Oc"->Tooltip["Oc","oxalocetate"],
	"ACoA"->Tooltip["ACoA","acetyl CoA"],
	"C"->Tooltip["C","citrate"],
	"CSH"->Tooltip["CSH","CoA-SH"],
	"cA"->Tooltip["cA","cis-aconitate"],
	"IsoC"->Tooltip["IsoC","isocitrate"],
	"Os"->Tooltip["Os","oxalosuccinate"],
	"\[Alpha]K"->Tooltip["\[Alpha]K","\[Alpha]-ketoglutarate"],
	"SCoA"->Tooltip["SCoA","succinyl-CoA"],
	"S"->Tooltip["S","succinate"],
	"Q"->Tooltip["Q","ubiquinone"],
	"F"->Tooltip["F","fumarate"],
	"QH2"-> Tooltip["QH2","ubiquinol"],
	"LM"->Tooltip["LM","L-malate"],
	"Pi"->Tooltip["Pi","phosphate"],
	"H2O"->Tooltip["H2O","dihydrogen monoxide"],
	"H+"->Tooltip["H+","hydrogen ion"],
	"CO2"->Tooltip["CO2","carbon dioxide"],
	"NAD+"->Tooltip["NAD+","nicotinamide adenine dinucleotide"],
	"GDP"->Tooltip["GDP","guanosine diphosphate"],
	"GTP"->Tooltip["GTP","guanosine-5'-triphosphate"],
	"NADH"->Tooltip["NADH","NADH"](*???*)}
		}
		]
];


ModelsSumPrivate := Join[BioModels, Models];

ReactionsSumPrivate := Join[BioReactions, Reactions];


GetReaction::usage = "GetReaction[model] returns a reaction if model is one of the built-in models. \
GetReaction[\"Models\"] returns all the available built-in models. See also Models.";


GetReaction::badarg = "Illegal argument of function GetReaction.";
GetReaction::nvmod = "Argument '`1`' is not a valid built-in (bio)model. \
See GetReaction[\"Models\"] or Models.";

SyntaxInformation[GetReaction]={"ArgumentsPattern"->{__}};

GetReaction["Models"] := ModelsSumPrivate; (*Models*)

GetReaction[x_?StringQ] := 
	If[
		MemberQ[ModelsSumPrivate,x], (*Models*)
			x /. ReactionsSumPrivate, (*Reactions*)
			Message[GetReaction::"nvmod",x];
			$Failed
	  ];

SetAttributes[GetReaction,Listable];

GetReaction[x__?StringQ] := GetReaction[{x}];

GetReaction[___] := (Message[GetReaction::"badarg"]; $Failed)


(* ::Subsection::Closed:: *)
(*fitting to standards*)


MyTotal[x_] := 
	Module[{ pos },

		Total[
			If[ DigitQ[StringTake[#,1]],

				pos = Position[DigitQ /@ Characters[#], False, 1, 1][[1,1]]-1;

				ToExpression[StringTake[#,pos]]*StringDrop[#,pos],

				#
			] & /@ x
		]
	];

EqualDependentQ := StringCount[#,"="] === 1&;

SS[x_] := 
	Module[{ r, read },
		read = (If[
					MemberQ[
						r=Quiet[
							ReadList[
								StringToStream[
									StringReplace[#, "."~~Longest[a__?DigitQ]~~x1_?(#=!="E"&&DigitQ[#]&)~~"+"~~b_?DigitQ :> "."~~a~~x1~~"E+"~~b]
								],
								Number]
						   ], $Failed],
						#,
					    Sequence@@r
				]) &/@ x;

				read /. {xx___,Longest[yy__?StringQ],Longest[zz__?NumberQ],uu___} :> 
					{
						StringSplit[
							StringReplace[
								StringJoin[yy],
								{"(+"~~Shortest[___]~~")":>""}
							],
							{"=>"->0,"="->1,"<=>"->1,"+"}
						] /. {{a__, 0, b__} :> {RightArrow[MyTotal[{a}],MyTotal[{b}]]}, {a__, 1, b__} :> {Equilibrium@@(MyTotal/@{{a},{b}})}}, {zz}
					}
];


(* MEMO IMPORT FUNCTIONS *)

MyMemoImport[filename_] := MyMemoImport[filename] = 
		Check[Import[filename,"Text"], Return[$Failed];, {Import::nffil,Import::chtype}];


MyImport[filename_, data_] := MyImport[filename, data] =
	Module[{ m, mm = filename, x, y, z, arr },

		m = MyMemoImport[mm]; (*text-kent beolvasva*)

		If[m===$Failed, Message[CHEMKINImport::"bimp"]; Return[$Failed];];

		m = StringReplace[m,"!"~~Shortest[___]~~"\n"->"\n"]; (*szemet-kommentek kitorolve*)

		If[ StringCases[m, "ELEMENTS"] =!= {}, x = "ELEMENTS"; , 
				If[StringCases[m, "ELEMS"] =!= {}, x = "ELEMS";, x = ""; ]];

		If[ StringCases[m, "SPECIES"] =!= {}, y = "SPECIES"; , y = ""; ];

		If[ StringCases[m, "REACTIONS"] =!= {}, z = "REACTIONS"; , z = ""; ];

		If[ z === "", 
				Message[CHEMKINImport::misdat,"REACTIONS"];		
				Return[{}],

				arr = Cases[
						StringSplit[
							First[StringCases[m,z~~Shortest[a___]~~"END":>a]],
								"\n"],
						x_?EqualDependentQ :> (SS[StringSplit[x]] /. Thread[{"M","HE","AR"}:>0] /. {"BASE",___}:>Sequence[])
					  ];
		];

		Switch[ {x =!= "",y =!= ""},

			{True, True},   Switch[data,
								"chemkinelements", Flatten[StringSplit[StringCases[m, x~~Shortest[a___]~~"END":>a], {Whitespace,"\n"}]],

								"chemkinspecies", Flatten[StringSplit[StringCases[m, y~~Shortest[a___]~~"END":>a], {Whitespace,"\n"}]],

								"chemkinarrhenius" | "chemkinreactions",
									If[data === "chemkinarrhenius", {#[[1,1]], #[[2]]}& /@ arr, DeleteDuplicates[Flatten[First /@ arr]]]
								],
				
			{True, False},  
							Switch[data,
								"chemkinelements", Flatten[StringSplit[StringCases[m, x~~Shortest[a___]~~"END":>a], {Whitespace,"\n"}]],

								"chemkinarrhenius" | "chemkinreactions",
									If[data === "chemkinarrhenius", {#[[1,1]], #[[2]]}& /@ arr, DeleteDuplicates[Flatten[First /@ arr]]],
								_,
									Message[CHEMKINImport::misdat,"SPECIES"];  
									Return[{}]
							],

			{False, True},   
							Switch[data,
								"chemkinspecies", Flatten[StringSplit[StringCases[m, y~~Shortest[a___]~~"END":>a], {Whitespace,"\n"}]],

								"chemkinarrhenius" | "chemkinreactions",
									If[data === "chemkinarrhenius", {#[[1,1]], #[[2]]}& /@ arr, DeleteDuplicates[Flatten[First /@ arr]]],
								_, 
									Message[CHEMKINImport::misdat,"ELEMENTS"]; 
									Return[{}]
							],

			{False, False}, 
							Switch[data,
								"chemkinarrhenius" | "chemkinreactions",
									If[data === "chemkinarrhenius", {#[[1,1]], #[[2]]}& /@ arr, DeleteDuplicates[Flatten[First /@ arr]]],
								_,
									Message[CHEMKINImport::misdat,"ELEMENTS"];  
									Message[CHEMKINImport::misdat,"SPECIES"];
									Return[{}]
							]
		]
	];


PropertiesCHEMKINImport = 
	{
		"chemkinelements","chemkinspecies","chemkinreactions","chemkinarrhenius"
	};


Options[CHEMKINImport] := {Memo -> False};


CHEMKINImport::usage = "CHEMKINImport[file] imports all the relevant data from a CHEMKIN file.";


CHEMKINImport::illarg = "The argument of CHEMKINImport needs a string with a file name the extension of which is .dat or .inp.";
CHEMKINImport::badname = "At least one of the arguments '`1`' is a non-identified property. Try CHEMKINImport[\"Properties\"].";
CHEMKINImport::bimp = "Some problem occurred in the import procedure. For possible further issues, see the function Import."; 
CHEMKINImport::misdat = "Missing DATA: the CHEMKIN file to be imported does not contain data for '`1`'.";
CHEMKINImport::badarg = "Illegal argument of function CHEMKINImport.";

(*SyntaxInformation[CHEMKINImport]={"ArgumentsPattern"->{__}};*)

CHEMKINImport[filename_?StringQ]["Properties"] := PropertiesCHEMKINImport;

CHEMKINImport[filename_?StringQ][] := Thread[PropertiesCHEMKINImport->CHEMKINImport[filename]["All"]];

CHEMKINImport[filename_?StringQ]["All"] := CHEMKINImport[filename][PropertiesCHEMKINImport];

CHEMKINImport[filename_?StringQ][data__?StringQ] := 
	(
	If[ Last[StringSplit[StringTrim[filename],"."]] =!= "dat" && Last[StringSplit[StringTrim[filename],"."]] =!= "inp",
		Message[CHEMKINImport::"illarg"];
		Return[$Failed];
	];
	If[Nand @@ Map[MemberQ[PropertiesCHEMKINImport,#]&,{data}],
		Message[CHEMKINImport::"badname",data];
		Return[$Failed];
	];
	(*Catch[*)
		Part[Map[(*Check[*)MyImport[filename, #](*,Throw[] ...]*)&, {data}], Length[{data}] /. x_/;x>1 :> (1;;x)
		]
	(*]*)
	);

CHEMKINImport[filename_?StringQ][data__] := 
	CHEMKINImport[filename][Sequence @@ (ToString /@ Flatten[{data}])];

(*CHEMKINImport[files__?StringQ] := CHEMKINImport /@ Flatten[{files}];*)

Format[CHEMKINImport[filename_?StringQ]] := 
	DynamicModule[{x, chi},
					x = Dynamic[Refresh[Round[Clock[Infinity]],UpdateInterval->1]];
					Monitor[chi = CHEMKINImport[filename][];,
									Column[{Row[{"CHEMKINImport is now importing and calculating ",Dynamic[Mod[First[x],5]/.{0->"",1->".",2->"..",3->"...",4->"...."}]}],
											ProgressIndicator[Dynamic[Clock[Infinity]],Indeterminate,ImageMargins->1,BaselinePosition->Center],
											Row[{x,". ","seconds passed"}]
									}]	
					];
					FinishDynamic[];
					chi(*FormatString[{reactions},{externals}]*)
	];

CHEMKINImport[___][___] := (Message[CHEMKINImport::"badarg"]; $Failed)


CHEMKINExport::usage = "CHEMKINExport[file,reaction] exports the reaction into a CHEMKIN file.";


CHEMKINExport::illarg = "The first argument of CHEMKINExport needs a string with a file name the extension of which is .dat.";
CHEMKINExport::fexist = "The given file already exists and will be overwritten.";
CHEMKINExport::wrgformat = "The given list of reactions has wrong format.";
CHEMKINExport::badarg = "Illegal argument of function CHEMKINExport.";

CHEMKINExport[filename_?StringQ, reactions_List] := 
	Module[{filestr, path, i, r, mmax, str, n, sq, step, rstep, arrh, w, nform, narrh, b = 0},
		
		If[ Last[StringSplit[StringTrim[filename],"."]] =!= "dat", 
			Message[CHEMKINExport::"illarg"];
			Return[$Failed];
		];

		If[FileNameDepth[filename] <= 1 || StringMatchQ[filename, "*:/*"], 
			path = filename;, path = FileNameJoin[{Directory[],filename}];
		];

		If[FileExistsQ[path], Message[CHEMKINExport::"fexist"]; ]; (*Do you want to delete or rename DeleteFile[path]*)
		
		filestr = Check[OpenWrite[path], Return[$Failed];, {General::"aofil", OpenWrite::"noopen"}];

		If[ filestr === $Failed, Return[$Failed]; ];

		r = Length[reactions];
		mmax = Max[StringLength[StringReplace[ToString[#]," "->""]]& /@ Flatten[{reactions}]];

		Do[
			step = reactions[[i]];			
			rstep = {};

			If[VectorQ[{step}],
				rstep = step;
				arrh = {};
			];

			If[Head[step]===List,
				If[Length[step]===2 && VectorQ[Last[step]],
					rstep = First[step];
					arrh = Last[step];
				];
			];
	
			If[rstep === {},
				Message[CHEMKINExport::wrgformat];
				Return[$Failed];
			];

			If[
				Length[
					Join@@(Position[rstep,#]&/@
								{Rule,RightArrow,LongRightArrow,Equilibrium,ReverseEquilibrium,LeftRightArrow,LongLeftRightArrow,LongLeftArrow,LeftArrow}
							)] =!= 1,
					Message[CHEMKINExport::wrgformat];
					Return[$Failed];,

					If[
						MemberQ[{LongLeftArrow, LeftArrow}, Head[rstep]],
							w = ToString[rstep[[2]]->rstep[[1]]];,
							w = ToString[rstep];
					];
					(*ns = If[StringQ[#], NumericQ[ToExpression[#,InputForm]], NumericQ[#]]& /@ arrh;*)
					w = StringReplace[w, " "->""];
					(*ns = NumericQ /@ arrh;
					If[And @@ ns,*)
						(*
						sq = MapIndexed[{#1, First[#2]}&, StringQ /@ arrh];
						str = StringToStream[""];
						n = Map[
							If[First[#], 
									Close[str]; 
									str = StringToStream[arrh[[Last[#]]]]; 
									N[Read[str, Number]], 
									N[arrh[[Last[#]]]]
							]&, sq
						];
						Close[str];
						*)
						nform = StringReplace[ToString[FortranForm[#]],"e"~~x_:>If[x==="-","E"~~x,"E+"~~x]]& /@ arrh;
						If[nform === {},
							narrh = {};,
							narrh = StringJoin[Append[(#<>StringJoin[ConstantArray[" ", Max[mmax+4-StringLength[#],2]]])& /@ Most[nform], Last[nform]]];
									(*nform[[1]]<>StringJoin[ConstantArray[" ",mmax+4-StringLength[nform[[1]]]]]<>
									nform[[2]]<>StringJoin[ConstantArray[" ",mmax+4-StringLength[nform[[2]]]]]<>nform[[3]];*)
						];
						(*nform = ScientificForm[#,{Infinity,10},NumberFormat->(StringJoin[{#1,"E",If[#3==="","+0",If[StringTake[#3,1]==="-",#3,"+"~~#3]],"\t"}]&)] &/@ n;*)
						(*nform = NumberForm[#,NumberFormat->(StringJoin[{#1,"E+",If[#3==="","0",#3],"\t"}]&)] &/@ n;*)
						WriteString[filestr, ToString[i]<>"\n"];
						WriteString[filestr,
							StringReplace[w, {" ":>"","->"->"=>","\[RightArrow]"->"=>","\[LongRightArrow]"->"=>","\[Equilibrium]":>"<=>","\[ReverseEquilibrium]":>"<=>","\[LeftRightArrow]":>"<=>","\[LongLeftRightArrow]":>"<=>"}]<>
												StringJoin[ConstantArray[" ", Max[mmax+10-StringLength[w],2]]]<>narrh<>"\n"]; 
						(*StringJoin[ToString /@ nform]*)
						b = i;(*,

						Message[CHEMKINExport::wrgformat];
						Return[$Failed];
					];*)
					
			];,
			{i, 1, r}
		];

		Close[filestr];

		If[ b === r, 
			Print["All given data have been written successfully in the file(path): " <> filename];,
			Return[$Failed];
		];

	];

(* {FileNameSetter[Dynamic[f]],Dynamic[f]} *)

CHEMKINExport[___][___] := (Message[CHEMKINExport::"badarg"]; $Failed)


(* ::Subsection::Closed:: *)
(*stoichiometry*)


(* ::Subsubsection::Closed:: *)
(*ReactionsData*)


SpeciesQ := UpperCaseQ[StringTake[ToString[#],1]]&;

ExternalsQ := OptionQ[{#}]&&Union[First/@Flatten[{#}]]==={ExternalSpecies}&;

NotExternalsQ := Not[ExternalsQ[#]]&;

AllExternals := Flatten[#/.Rule->List/.ExternalSpecies->Sequence[]]&;


SymbolQ := Head[#]===Symbol&;

PlusQ := Head[#]===Plus&;

TimesQ := Head[#]===Times&;

SubscriptQ := Head[#]===Subscript&;

SuperscriptQ := Head[#]===Superscript&;

RightArrowQ := (Head[#]===RightArrow&&Length[#]===2&&
				And@@Map[StringQ[#]||SymbolQ[#]||PlusQ[#]||TimesQ[#]||SubscriptQ[#]||SuperscriptQ[#]&,#/.RightArrow->List])&;

LeftArrowQ := (Head[#]===LeftArrow&&Length[#]===2&&
				And@@Map[StringQ[#]||SymbolQ[#]||PlusQ[#]||TimesQ[#]||SubscriptQ[#]||SuperscriptQ[#]&,#/.LeftArrow->List])&;

EquilibriumQ := (Head[#]===Equilibrium&&Length[#]===2&&
				And@@Map[StringQ[#]||SymbolQ[#]||PlusQ[#]||TimesQ[#]||SubscriptQ[#]||SuperscriptQ[#]&,#/.Equilibrium->List])&;

ReverseEquilibriumQ := (Head[#]===ReverseEquilibrium && Length[#]===2 && 
						And@@Map[StringQ[#]||SymbolQ[#]||PlusQ[#]||TimesQ[#]||SubscriptQ[#]||SuperscriptQ[#]&,#/.ReverseEquilibrium->List])&;


ReactionsToList[{reactions__?RightArrowQ}]:= DeleteDuplicates[{reactions}];

ReactionsToList[{reactions__?LeftArrowQ}] := DeleteDuplicates[Reverse /@ {reactions} /. LeftArrow->RightArrow];

ReactionsToList[{reactions__?EquilibriumQ}]:=
		DeleteDuplicates[Sequence@@({#,Reverse[#]}/.Equilibrium->RightArrow)&/@{reactions}];

ReactionsToList[{reactions__?ReverseEquilibriumQ}]:=
		DeleteDuplicates[Sequence@@({Reverse[#],#}/.ReverseEquilibrium->RightArrow)&/@{reactions}];

ReactionsToList[{reactions__}]:=
	Module[{form, sp, a},
			form = DeleteDuplicates[ToCanonicalForm[reactions]];
			DeleteDuplicates[
			Flatten[
				Map[
					Switch[#,
						_?RightArrowQ,
							#,
						_?LeftArrowQ,
							Reverse[#] /. LeftArrow->RightArrow,
						_?EquilibriumQ,
							{#, Reverse[#]} /. Equilibrium->RightArrow,
						_?ReverseEquilibrium,
							{Reverse[#], #} /. ReverseEquilibrium->RightArrow,
						_,
							(sp = Flatten[ToBoxes[#,StandardForm]/.RowBox->List/.Thread[{"(",")","{","}",","}:>Sequence[]]];
							(*MakeBoxes*)
							a = Sort[Flatten[Position[sp,#]&/@{"\[Equilibrium]","\[RightArrow]","\[ReverseEquilibrium]","\[LeftArrow]"}]];
							(Partition[
									Riffle[
										ReleaseHold[
											MakeExpression[RowBox[Take[sp,#+{1,-1}]],StandardForm]
										] &/@ Partition[Join[{0},a,{Length[sp]+1}],2,1],
									Part[sp, a]],
							3,2] /. {
										{x_,"\[Equilibrium]",y_}->Sequence[RightArrow[x,y],RightArrow[y,x]],
										{x_,"\[ReverseEquilibrium]",y_}->Sequence[RightArrow[y,x],RightArrow[x,y]],
										{x_,"\[LeftArrow]",y_}->RightArrow[y,x],
										{x_,"\[RightArrow]",y_}->RightArrow[x,y]
									 }
							)
							)
					]&,
					form
				] /. RightArrow[A_,A_]->Sequence[]
			]
			]
];


complextospecies[cplxs_] := DeleteDuplicates[
								((Flatten[{cplxs} /. Plus:>List] 
										/. r_*A_?StringQ :> A /. A_?StringQ*r_ :> A)
										/. r_*A_?SpeciesQ :> A /. A_?SpeciesQ*r_ :> A)
										/. A_?StringQ :> A /. A_?SpeciesQ :> A /. r_*"0":>"0"
							];
									(*collect, symbolic coefficient, multip coeff*)
									(*Union*)

complextocoefflist[cplxs_] := Block[{spec},
								spec = complextospecies[cplxs];
								Thread[List[spec, Coefficient[cplxs, spec]]]
							 ];

subgraphedges[graph_,vtceslist_] := 
			If[ Length[#] === 1, {{}, #}, Join[{EdgeList[Subgraph[graph, #]] /. DirectedEdge->Rule}, {#}]] & /@ vtceslist;



Rdata[{reactions__},externals__] := Rdata[{reactions},externals] = 
	Module[
			{ reacs,
			  m, builtin, exs, exsrule, 
			  steps, steprarrow, steprule, fhjgraph, complexes, vgggraph, vgraph, vggraphwzerocomplex, sp, 
			  indexed, a, b, e, nloc, lloc, sloc,
			  lsp, lrs, alpha, beta, gamma, fhjweakvertices, fhjweakcomponents, 
			  fhjstrongvertices, fhjstrongcomponents,
			  fhjterminalstrongvertices, fhjterminalstrongcomponents
			},

			reacs = Flatten[{reactions}];

			If[(m = Cases[reacs,_String])=!={},

				builtin = Check[(#->GetReaction[#])& /@ m, Return[$Failed];, GetReaction::"nvmod"];(**)

				Rdata[DeleteDuplicates[Join[Flatten[reacs /. builtin], reacs /. Thread[m->Sequence[]]]], externals]
				(*DeleteDuplicates*),

				exs = DeleteDuplicates[Prepend[AllExternals[{externals}],"0"]/.(0->"0")]; (*internalspecies megadas*)

				exsrule = Thread[exs -> 0];
				(**)
				(* MAIN *)

				steps = ReactionsToList[reacs] /. (0->"0"); (*"0" - zero complex*)

				steprarrow = DeleteDuplicates[steps /. exsrule /. (0->"0")];

				steprule = steprarrow /. RightArrow -> Rule;

				(*Print[steprule];*)
				fhjgraph = Graph[steprule];

				complexes = DeleteDuplicates[Flatten[steprarrow /. RightArrow -> List]]; (*Sort*)
				(*VertexList[fhjgraph];*)

				vgggraph = Flatten[(Thread /@ {#->complextocoefflist[Last[#]], complextocoefflist[First[#]]->#}) &/@ steprarrow];
				vgraph = vgggraph /. Rule[A_RightArrow,B_]:>{Rule[A, First[B]], Last[B]} /. Rule[A_,B_RightArrow]:>{Rule[First[A], B], Last[A]};
				(*Print[vgraph];*)

				sp = complextospecies[complexes] /. "0" -> Sequence[]; (*komplexek csak belso anyagfajta erejeig egyertelmuek*)

				indexed = Flatten[MapIndexed[#1->First[#2]&, #]& /@ {sp, steprarrow}];
				(*Print[indexed];*)

				{lsp, lrs} = Length /@ {sp, steprule};

				vggraphwzerocomplex = vgraph /. {Rule["0",y_RightArrow],z_} :> Sequence[] /. {Rule[x_RightArrow,"0"],z_} :> Sequence[];

				a = Cases[vggraphwzerocomplex, {Rule[x_,y_RightArrow],z_}:> Rule[{x,y} /. indexed, z]];
				b = Cases[vggraphwzerocomplex, {Rule[x_RightArrow,y_],z_}:> Rule[{y,x} /. indexed, z]];
				(*Print[{Join[{a,b}],lsp,lrs}];*)
				alpha = SparseArray[a, {lsp, lrs}];
				beta = SparseArray[b, {lsp, lrs}];
				gamma = beta-alpha;
				(*SparseArray[Join[a,b] //. {x___,Rule[{y_,z_},k_],u___,Rule[{y_,z_},l_],v___}:>{x,Rule[{y,z},k+l],u,v}, {lsp,lrs}];*)
				(*Print[gamma];*)

				fhjweakvertices = ConnectedComponents[UndirectedGraph[fhjgraph]];
				fhjweakcomponents = subgraphedges[fhjgraph, fhjweakvertices];

				fhjstrongvertices = ConnectedComponents[fhjgraph];
				fhjstrongcomponents = subgraphedges[fhjgraph, fhjstrongvertices];

(*				volpertgraph = Graph[First /@ vgraph];
				volpertbasis = Flatten[Map[If[
											VertexInDegree[volpertgraph,#]===0, 
											#, {}
										    ]&, 
											Append[sp,"0"]
									      ] /. "0" :> {}
							         ];
*)
(*				fhjstronglength = Length[fhjstrongvertices];
				fhjstronggraph = Graph[
									DeleteDuplicates[
										steprarrow /. 
											Sort[Flatten[MapThread[Thread[#1 -> #2]&, {fhjstrongvertices, Range[fhjstronglength]}]], 
													Length[First[#1]]>Length[First[#2]]&
												]
									] /. (RightArrow[y_,y_] :> RightArrow[y,"\[Infinity]"]) /. RightArrow -> Rule];
				(*Print[fhjstronggraph];*)
				volpertbasis = Map[If[
									VertexInDegree[fhjstronggraph,#]===0, 
									Append[fhjstrongcomponents[[#]],Flatten[(First /@ complextocoefflist[#])& /@ fhjstrongcomponents[[#,2]]]], {{}}
									]&,
									(VertexList[fhjstronggraph] /. "\[Infinity]":>Sequence[])
								] /. {{}} :> Sequence[];
*)		

				fhjterminalstrongvertices = Map[
												Complement[
													Union @@ Map[Cases[steprarrow, RightArrow[#,y_]:>y]&, #], #
												]==={} /. True-># /. False->Sequence[]&, fhjstrongvertices
											];
				fhjterminalstrongcomponents = subgraphedges[fhjgraph, fhjterminalstrongvertices];

				Dispatch[{

					"species" -> sp,

					"M" -> lsp,

					"externalspecies" -> (e = exs /. "0" -> Sequence[]),
					(*e = If[MemberQ[complexes,"0"], exs, Rest[exs]]*)

					"E" -> Length[e],

					"reactionsteps" -> steps, (*steprule*)

					"R" -> lrs,

					"complexes" -> complexes,

					"fhjgraphedges" -> steprule, 
					(*FromOrderedPairs[steprule /. Rule[c1_,c2_] :> {Position[complexes, c1][[1, 1]], Position[complexes, c2][[1, 1]]}]*)

					"fhjweaklyconnectedcomponents" -> fhjweakcomponents,
					(*(fhjweak = WeakComponents[steprule])*)

					"fhjstronglyconnectedcomponents" -> fhjstrongcomponents,
					(*(fhjstrong = StrongComponents[steprule])*)

					"fhjterminalstronglyconnectedcomponents" -> fhjterminalstrongcomponents,

					(*"volpertbasis" -> {},*)

					"N" -> (nloc = Length[complexes]),

					"L" -> (lloc = Length[fhjweakvertices]),

					"S" -> (sloc = MatrixRank[gamma]),
					(*MatrixRank[gamma = First[Differences[ab = (Normal/@{a,b})]]]*)

					"deficiency" -> "\[Delta]=N-L-S="<>ToString[nloc]<>"-"<>ToString[lloc]<>"-"<>ToString[sloc]<>"="<>ToString[nloc-lloc-sloc],

					"\[Alpha]" -> alpha,

					"\[Beta]" -> beta,

					"\[Gamma]" -> gamma,

					"reactionsteporders" -> (Total /@ Transpose[alpha]),

					"variables" -> (Subscript[Global`c, #] & /@ sp),

					"volpertgraph" -> {vggraphwzerocomplex, Drop[indexed, lsp]}

				}]
			]
	];


PropertiesReactionsData = 
	{
		"species","M","externalspecies","E","reactionsteps","R","complexes",
		"fhjgraphedges","fhjweaklyconnectedcomponents","fhjstronglyconnectedcomponents","fhjterminalstronglyconnectedcomponents",
		"N","L","S","deficiency","\[Alpha]","\[Beta]","\[Gamma]","reactionsteporders","variables","volpertgraph"
	};


ReactionsForm[reactionlist_] := 
	If[Length[Flatten[reactionlist]]<=2,ToString[reactionlist,StandardForm],"{"<>ToString[First[reactionlist],StandardForm]<>",...,"<>ToString[Last[reactionlist],StandardForm]<>"}"];


ReactionsData::usage = "ReactionsData[{reactions},externalspecies] \
returns the most important characteristics of the set of reactions such as the species, number of species, \
(possible) external species, complexes, number of complexes, (genuine) reaction steps, the edges of the Feinberg-Horn-Jackson (FHJ) graph, \
weakly and strongly connected (terminal) components of the FHJ graph, Volpert graph, number of linkage classes, deficiency, \[Alpha], \[Beta] and \[Gamma], the order of the reaction steps. \  
With the optional argument externalspecies one can treat external species. \n\
ReactionsData[{reactions},externalspecies][properties] gives only those characteristics that are listed in \
properties. Check ReactionsData[\"Properties\"].";


ReactionsData::badarg = "Illegal argument of function ReactionsData.";
ReactionsData::badname = "At least one of the arguments '`1`' is a non-identified property. Try \
ReactionsData[\"Properties\"].";
ReactionsData::wrreac = "Argument `1` may not be in a correct form. Check OpenReactionKineticsPalette[].";


SyntaxInformation[ReactionsData]={"ArgumentsPattern"->{__,OptionsPattern[]}};

Options[ReactionsData] = {ExternalSpecies->{}};

ReactionsData["Properties"] := PropertiesReactionsData;

ReactionsData[{},externals:(_?VectorQ|OptionsPattern[])][] := {};

ReactionsData[{reactions__},externals:(_?VectorQ|OptionsPattern[])][] := 
	Check[Rdata[{reactions},{externals}],
					Message[ReactionsData::"wrreac",ReactionsForm[{reactions}]];
					Return[$Failed]
	];

ReactionsData[{},externals:(_?VectorQ|OptionsPattern[])]["Properties"] := PropertiesReactionsData;

ReactionsData[{reactions__},externals:(_?VectorQ|OptionsPattern[])]["Properties"] := PropertiesReactionsData;
(*ReactionsData[Longest[reactions__?NotExternalsQ],Shortest[externals:(___?ExternalsQ)]]["Properties"] := PropertiesReactionsData;*)

ReactionsData[{},externals:(_?VectorQ|OptionsPattern[])]["All"] := {};

ReactionsData[{reactions__},externals:(_?VectorQ|OptionsPattern[])]["All"] := ReactionsData[{reactions},externals][PropertiesReactionsData];
(*ReactionsData[Longest[reactions__?NotExternalsQ],Shortest[externals:(___?ExternalsQ)]]["All"] := ReactionsData[{reactions},externals][PropertiesReactionsData];*)

ReactionsData[{},externals:(_?VectorQ|OptionsPattern[])][data__?StringQ] := 
	(
	If[Nand @@ Map[MemberQ[PropertiesReactionsData,#]&,{data}],
		Message[ReactionsData::"badname",data];
	];
	Return[{}]
	);

ReactionsData[{reactions__},externals:(_?VectorQ|OptionsPattern[])][data__?StringQ] := 
	(
	If[Nand @@ Map[MemberQ[PropertiesReactionsData,#]&,{data}],
		Message[ReactionsData::"badname",data];
	];
	Catch[
		Part[
			ReplaceAll[{data},
					Check[Rdata[{reactions},{externals}],
							Throw[
								Message[ReactionsData::"wrreac",ReactionsForm[{reactions}]];
								Return[$Failed]
							]
					]
			], Length[{data}] /. x_/;x>1 :> (1;;x)
		]
	]
	);

ReactionsData[{},externals:(_?VectorQ|OptionsPattern[])][data__] := ReactionsData[{}][Sequence @@ ToString /@ Flatten[{data}]];

ReactionsData[{reactions__},externals:(_?VectorQ|OptionsPattern[])][data__] := 
		Block[{d, pos, rdata, sol = {data}},

				d = Flatten[sol];

				rdata = ReactionsData[{reactions},externals][Sequence @@ ToString /@ d];
				(*Check[, Return[$Failed], {ReactionsData::wrreac,ReactionsData::badarg,ReactionsData::badname}]];*)

				If[Length[d]===1, rdata = {rdata}];				

				pos = Flatten[Position[sol, #]& /@ d, 1];

				Do[Part[sol, Sequence@@pos[[n]]] = rdata[[n]], {n, Length[d]}];
				
				If[ Length[sol] === 1,
					Flatten[sol, 1],
					sol
				]
		];

Format[ReactionsData[{reactions__},externals:(_?VectorQ|OptionsPattern[])]] := 
	DynamicModule[{x, rd},
					x = Dynamic[Refresh[Round[Clock[Infinity]],UpdateInterval->1]];
					Monitor[rd = ReactionsData[{reactions},externals][];,
									Column[{Row[{"ReactionsData is now calculating ",Dynamic[Mod[First[x],5]/.{0->"",1->".",2->"..",3->"...",4->"...."}]}],
											ProgressIndicator[Dynamic[Clock[Infinity]],Indeterminate,ImageMargins->1,BaselinePosition->Center],
											Row[{x,". ","seconds passed"}]
									}]	
					];
					FinishDynamic[];
					rd(*FormatString[{reactions},{externals}]*)
	];

ReactionsData[___][___] := (Message[ReactionsData::"badarg"]; Return[$Failed])


(* ::Subsubsection::Closed:: *)
(*ToCanonicalForm, ToReversible, ReversibleFHJRepresentation, ReversibleQ, WeaklyReversibleQ*)


ToCanonicalForm::usage = "ToCanonicalForm[reactions] returns the canonical form of reactions \
using only right, left and left-right arrows.";

ToCanonicalForm::badarg = "Illegal argument of function ToCanonicalForm.";


SyntaxInformation[ToCanonicalForm]={"ArgumentsPattern"->{__}};

ToCanonicalForm[reactions__]:=
	Flatten[{reactions} /. Reactions] /. 
					  Join[Thread[{LeftRightArrow,DoubleLeftRightArrow,LongLeftRightArrow,DoubleLongLeftRightArrow,RightArrowLeftArrow,TwoWayRule,UndirectedEdge}->Equilibrium],
									{LeftArrowRightArrow->ReverseEquilibrium},
									Thread[{Rule,ShortRightArrow,DoubleRightArrow,LongRightArrow,DoubleLongRightArrow,DirectedEdge}->RightArrow],
									Thread[{ShortLeftArrow,DoubleLeftArrow,LongLeftArrow,DoubleLongLeftArrow}->LeftArrow]];

ToCanonicalForm[___] := (Message[ToCanonicalForm::"badarg"]; $Failed)


ToReversible::usage = "ToReversible[reactions] make all the reaction steps of reactions reversible.";

ToReversible::badarg = "Illegal argument of function ToReversible.";


SyntaxInformation[ToReversible]={"ArgumentsPattern"->{__}};

ToReversible[reactions__] := ToCanonicalForm[reactions]/.RightArrow->Equilibrium/.LeftArrow->Equilibrium;

ToReversible[___] := (Message[ToReversible::"badarg"]; $Failed)


ReversibleFHJRepresentation::usage = "ReversibleFHJRepresentation[{reactions},options] results in the Feinberg-Horn-Jackson \
graph making all its edges reversible.";

ReversibleFHJRepresentation::badarg = "Illegal argument of function ReversibleFHJRepresentation.";


SyntaxInformation[ReversibleFHJRepresentation]={"ArgumentsPattern"->{__,OptionsPattern[]}};

Options[ReversibleFHJRepresentation] = {ExternalSpecies->{}};

ReversibleFHJRepresentation[{reactions__},externals:(_?VectorQ|OptionsPattern[])] :=
	ReplaceRepeated[
		Check[ReactionsData[{reactions},externals]["fhjgraphedges"],Abort[];,{ReactionsData::wrreac,ReactionsData::badarg}],
		{x___, Rule[y_,z_], t___, Rule[z_,y_], u___} :> {x, Equilibrium[y,z], t, u}
	];

(*avagy: Complement[EdgeList[g], EdgeList[ReverseGraph[g]]], where g = FHJ Graph object*)

ReversibleFHJRepresentation[___] := (Message[ReversibleFHJRepresentation::"badarg"]; $Failed)


ReversibleQ::usage = "ReversibleQ[{reactions},options] yields True if and only if the reaction is reversible.";

ReversibleQ::badarg = "Illegal argument of function ReversibleQ.";


SyntaxInformation[ReversibleQ]={"ArgumentsPattern"->{__}};

Options[ReversibleQ] = {ExternalSpecies->{}};

ReversibleQ[{reactions__},externals:(_?VectorQ|OptionsPattern[])] := 
		FreeQ[ReversibleFHJRepresentation[{reactions},externals],_Rule];

ReversibleQ[___] := (Message[ReversibleQ::"badarg"]; $Failed)


WeaklyReversibleQ::usage = "WeaklyReversibleQ[{reactions},options] yields True if and only if the reaction is weakly reversible.";

WeaklyReversibleQ::badarg = "Illegal argument of function WeaklyReversibleQ.";


SyntaxInformation[WeaklyReversibleQ]={"ArgumentsPattern"->{__}};

Options[WeaklyReversibleQ] = {ExternalSpecies->{}};

WeaklyReversibleQ[{reactions__},externals:(_?VectorQ|OptionsPattern[])] := 
		(Length[Check[ReactionsData[{reactions},externals]["fhjstronglyconnectedcomponents"],Abort[];,{ReactionsData::wrreac,ReactionsData::badarg}]] === 1);

WeaklyReversibleQ[___] := (Message[WeaklyReversibleQ::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*FilterReactions, FromStoichiometry, DeleteAutocatalysis*)


FilterReactions::usage = "FilterReactions[{reactions},species,options] filters all the reaction steps in which \
the given set of species shows up as either reactant or product or (reactant or product) species. \
This can be controlled by option Side \[Rule] \"Reactant\", \"Product\" or \"All\".";

FilterReactions::badsp = "Some of the species may be not contained in the list of internal species of the reaction.";
FilterReactions::badopt = "Option Side only takes three values: \"Reactant\", \"Product\" or \"All\".";
FilterReactions::badarg = "Illegal argument of function FilterReactions.";

Options[FilterReactions]:={ExternalSpecies -> {}, Side -> "All"};

FilterReactions[{reactions__},specs_?VectorQ,opts : OptionsPattern[]]:=
	Module[{ options, rdata, fhj, sp, spp, side },

		options = Flatten[{opts, Options[FilterReactions]}];
		rdata = Check[ReactionsData[{reactions},ExternalSpecies /. options]["fhjgraphedges","species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		If[rdata === $Failed, Return[$Failed];];
		{fhj, sp} = rdata;
			
		If[(spp = Intersection[specs,sp]) =!= {},

			If[ MemberQ[{"Reactant","Product","All"}, (side = Side /. options)],	
				Switch[side,
					"Reactant", 
						Flatten[If[Intersection[complextospecies[First[#]],specs]=!={}, #, {}] &/@ fhj],
					"Product", 
						Flatten[If[Intersection[complextospecies[Last[#]],specs]=!={}, #, {}] &/@ fhj],
					"All", 
						Flatten[If[Intersection[complextospecies[# /. Rule -> Plus],specs]=!={}, #, {}] &/@ fhj]
				],
				Message[FilterReactions::"badopt"]; 
				Return[$Failed]
			],
			Message[FilterReactions::"badsp"];
			Return[$Failed]
		]
	];

FilterReactions[___][___] := (Message[FilterReactions::"badarg"]; $Failed)


FromStoichiometry::usage = "FromStoichiometry[\[Alpha],\[Beta],species] builds up the set of reaction steps from the matrixes of stoichiometric coefficients \[Alpha] and \[Beta] \
using species. The species argument is optional.";

FromStoichiometry::dimmx = "Matrixes \[Alpha] and \[Beta] must have the same size.";
FromStoichiometry::vars = "The number of species does not match with the dimensions of \[Alpha] and \[Beta].";
FromStoichiometry::badarg = "Illegal argument of function FromStoichiometry.";

FromStoichiometry[alpha_?MatrixQ,beta_?MatrixQ] := 
	FromStoichiometry[alpha, beta, Array[Global`X[#]&,Length[alpha]]];

FromStoichiometry[alpha_?MatrixQ,beta_?MatrixQ,specs_?VectorQ]:=
	Module[{ m, r },
		{m, r} = Dimensions[alpha];

		If[ {m, r} === Dimensions[beta],

			If[ m === Length[specs],

				Thread[specs . alpha -> specs . beta], (*Rightarrow?*)

				 Message[FromStoichiometry::"vars"];
				 Return[$Failed];
			],

		 Message[FromStoichiometry::"dimmx"];
		 Return[$Failed];
		]
	];

FromStoichiometry[___][___] := (Message[FromStoichiometry::"badarg"]; $Failed)


DeleteAutocatalysis::usage = "DeleteAutocatalysis[{reactions}] and DeleteAutocatalysis[\[Alpha],\[Beta],species] deletes the autocatalytic steps of the reaction.";

DeleteAutocatalysis::dimmx = "Matrixes \[Alpha] and \[Beta] must have the same size.";
DeleteAutocatalysis::vars = "The number of species does not match with the dimensions of \[Alpha] and \[Beta].";
DeleteAutocatalysis::badarg = "Illegal argument of function DeleteAutocatalysis.";

Options[DeleteAutocatalysis] := {ExternalSpecies->{}};

DeleteAutocatalysis[{reactions__}, OptionsPattern[]] := 
	DeleteAutocatalysis @@ Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["\[Alpha]","\[Beta]","species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

DeleteAutocatalysis[{reactions__}, variables_?VectorQ, OptionsPattern[]] := 
	DeleteAutocatalysis[Sequence @@ Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["\[Alpha]","\[Beta]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],variables];

DeleteAutocatalysis[alpha_?MatrixQ,beta_?MatrixQ] := 
	DeleteAutocatalysis[alpha, beta, Array[Global`X[#]&,Length[alpha]]];

DeleteAutocatalysis[alpha_?MatrixQ,beta_?MatrixQ,specs_?VectorQ]:=
	Module[{ gamma, m, r, reacs, revreacs },
		{m, r} = Dimensions[alpha];
		If[ {m, r} === Dimensions[beta],

			If[ m === Length[specs],
				gamma = beta-alpha;

				reacs = DeleteDuplicates[Thread[specs . NegativePart[gamma] -> specs . PositivePart[gamma]]];
				(*Thread[Rule@@(species.#[gamma]&/@{NegativePart,PositivePart})]]*)
				revreacs = Reverse /@ reacs;

				Join[Complement[reacs,revreacs] /. Rule -> RightArrow, Union[Sort /@ Intersection[reacs, revreacs]] /. Rule -> Equilibrium],
		
				Message[DeleteAutocatalysis::vars];
				Return[$Failed];
			],

		 Message[DeleteAutocatalysis::dimmx];
		 Return[$Failed];
		]
	];

DeleteAutocatalysis[___][___] := (Message[DeleteAutocatalysis::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*MinFHJWeaklyConnectedComponents, MinFHJStronglyConnectedComponents, MaxFHJWeaklyConnectedComponents, MaxFHJStronglyConnectedComponents*)


MinFHJWeaklyConnectedComponents::usage = "MinFHJWeaklyConnectedComponents[{reactions},options] returns the smallest weakly connected components of the given reaction.";

Options[MinFHJWeaklyConnectedComponents] := {ExternalSpecies->{}};

MinFHJWeaklyConnectedComponents[{reactions__}, OptionsPattern[]] := 
	Module[{fhjwcc, lengthvtxes, minvtxes},

		fhjwcc = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["fhjweaklyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		lengthvtxes = Length[Last[#]] &/@ fhjwcc;

		minvtxes = Min[lengthvtxes];

		(*fhjwcc[[Flatten[Position[lengthvtxes, minvtxes]]]]*)
		Select[fhjwcc, Length[Last[#]]===minvtxes &]

	];

MinFHJWeaklyConnectedComponents[___][___] := (Message[MinFHJWeaklyConnectedComponents::"badarg"]; $Failed)


MinFHJStronglyConnectedComponents::usage = "MinFHJStronglyConnectedComponents[{reactions},options] returns the smallest strongly connected components of the given reaction.";

Options[MinFHJStronglyConnectedComponents] := {ExternalSpecies->{}};

MinFHJStronglyConnectedComponents[{reactions__}, OptionsPattern[]] := 
	Module[{fhjscc, lengthvtxes, minvtxes},

		fhjscc = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["fhjstronglyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		lengthvtxes = Length[Last[#]] &/@ fhjscc;

		minvtxes = Min[lengthvtxes];

		(*fhjscc[[Flatten[Position[lengthvtxes, minvtxes]]]]*)
		Select[fhjscc, Length[Last[#]]===minvtxes &]

	];

MinFHJStronglyConnectedComponents[___][___] := (Message[MinFHJStronglyConnectedComponents::"badarg"]; $Failed)


MaxFHJWeaklyConnectedComponents::usage = "MaxFHJWeaklyConnectedComponents[{reactions},options] returns the largest weakly connected components of the given reaction.";

Options[MaxFHJWeaklyConnectedComponents] := {ExternalSpecies->{}};

MaxFHJWeaklyConnectedComponents[{reactions__}, OptionsPattern[]] := 
	Module[{fhjwcc, lengthvtxes, maxvtxes},

		fhjwcc = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["fhjweaklyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		lengthvtxes = Length[Last[#]] &/@ fhjwcc;

		maxvtxes = Max[lengthvtxes];

		(*fhjwcc[[Flatten[Position[lengthvtxes, maxvtxes]]]]*)
		Select[fhjwcc, Length[Last[#]]===maxvtxes &]

	];

MaxFHJWeaklyConnectedComponents[___][___] := (Message[MaxFHJWeaklyConnectedComponents::"badarg"]; $Failed)


MaxFHJStronglyConnectedComponents::usage = "MaxFHJStronglyConnectedComponents[{reactions},options] returns the largest strongly connected components of the given reaction.";

Options[MaxFHJStronglyConnectedComponents] := {ExternalSpecies->{}};

MaxFHJStronglyConnectedComponents[{reactions__}, OptionsPattern[]] := 
	Module[{fhjscc, lengthvtxes, maxvtxes},

		fhjscc = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["fhjstronglyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		lengthvtxes = Length[Last[#]] &/@ fhjscc;

		maxvtxes = Max[lengthvtxes];

		(*fhjscc[[Flatten[Position[lengthvtxes, maxvtxes]]]]*)
		Select[fhjscc, Length[Last[#]]===maxvtxes &]

	];

MaxFHJStronglyConnectedComponents[___][___] := (Message[MaxFHJStronglyConnectedComponents::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*ShowFHJGraph*)


GraphPlotFunctionQ := MemberQ[{"Graph","GraphPlot","LayeredGraphPlot","GraphPlot3D","TreePlot"},#]&;


ShowFHJGraph::usage = "ShowFHJGraph[{reactions},rratecoeffs,options] displays the Feinberg-Horn-Jackson graph of the reaction \
with rratecoeffs as reaction rate coeffcients.";

ShowFHJGraph::badarg = "Illegal argument of function ShowFHJGraph.";
ShowFHJGraph::funcarg = "Argument `1` is not a valid graph plot function. Try Graph, GraphPlot, GraphPlot3D, LayeredGraphPlot or TreePlot.";
ShowFHJGraph::ccols = "The list of colors has wrong shape.";
ShowFHJGraph::srates = "The list of reaction rate coefficients has wrong shape.";
ShowFHJGraph::sccomps = "The list of colors (for strongly connected components) may have wrong shape or complex colors are also given.";

Options[ShowFHJGraph]:={ExternalSpecies -> {}, ComplexColors -> {}, PlotFunction -> "GraphPlot", StronglyConnectedComponentsColors -> {}, Numbered -> False};

(*SyntaxInformation[ShowFHJGraph]={"ArgumentsPattern"->{{__},___}};*)

ShowFHJGraph[{reactions__}, rates_List:{}, opts___?OptionQ] := 
	Module[{ cc, externals, rdata, exs, 
			 n, r, fhj, cxs, fhjscc, graphfunc, nmd,
			 nopt, c, sc, rscc, rs 
		   },

			cc = FilterRules[{opts},Options[ShowFHJGraph]];
			externals = FilterRules[{opts},ExternalSpecies];

			rdata = Check[ReactionsData[{reactions},Flatten[externals]]["N","R","fhjgraphedges","complexes","fhjstronglyconnectedcomponents"], Return[$Failed], {ReactionsData::wrreac,ReactionsData::badarg}];
			If[rdata === $Failed, Return[$Failed]];

			{n, r, fhj, cxs, fhjscc} = rdata;

			graphfunc = (PlotFunction /. cc) /. PlotFunction->"GraphPlot";			
			If[GraphPlotFunctionQ[graphfunc],
				graphfunc = ToExpression[graphfunc],
				Message[ShowFHJGraph::"funcarg",graphfunc];
				Return[$Failed]				
			];

			fhjscc = Last /@ fhjscc; 
			If[((Numbered /. cc) /. Numbered -> False)===True, 
				cxs = Thread[Rule[cxs,Range[n]]];
				fhj = fhj /. cxs;
				fhjscc = fhjscc /. cxs;
				cxs = Last /@ cxs;
			];

			nopt = FilterRules[{opts},Options[graphfunc]];
			c = ComplexColors /. cc;
			sc = StronglyConnectedComponentsColors /. cc;

			If[sc =!= {} && sc =!= StronglyConnectedComponentsColors,

				If[(VectorQ[sc] && Length[fhjscc] === Length[sc]) && (c === {} || c === ComplexColors),
					rscc = fhj /. Flatten[Thread[First[#]->Thread[Style[First[#], Plain, Last[#]]]] &/@ Thread[{fhjscc,sc}]];
						If[rates==={},
							graphfunc[rscc,nopt],
							If[r===Length[rates],
								graphfunc[Thread[{rscc,rates}],nopt],
								Message[ShowFHJGraph::"srates"];
								Return[$Failed]
							]
						],
					Message[ShowFHJGraph::"sccomps"];
					$Failed
				],

				If[c =!= {} && c =!= ComplexColors,

					If[VectorQ[c] && n===Length[c],
						rs = fhj /. (First[#]->Style[First[#], Plain, Last[#]] &/@ Thread[{cxs,c}]);
						If[rates==={},
							graphfunc[rs,nopt],
							If[r===Length[rates],
								graphfunc[Thread[{rs,rates}],nopt],
								Message[ShowFHJGraph::"srates"];
								Return[$Failed]
							]
						],
						Message[ShowFHJGraph::"ccols"];
						$Failed
					],

					If[rates === {},
						GraphPlot[fhj,nopt],
						If[r===Length[rates],
								graphfunc[Thread[{fhj,rates}],nopt],
								Message[ShowFHJGraph::"srates"];
								Return[$Failed]
							]
					]
				]
			]
	];

ShowFHJGraph[___] := (Message[ShowFHJGraph::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*Volpertindexing, ShowVolpertGraph, AcyclicVolpertGraphQ, VolpertSpecies*)


VolpertIndexing::usage = "VolpertIndexing[{reactions},initspecies,options] gives the Volpert indexes \
of species and reactions of the given reaction with respect to the initial species given by initspecies.";

VolpertIndexing::badarg = "Illegal argument of function VolpertIndexing.";
VolpertIndexing::specerr = "None of the species in `1` is in the list of internal species of the reaction.";
VolpertIndexing::verbtrumis = "The species and/or reaction steps may not appear in the given reaction.";

VolpertIndexing[alpha_?MatrixQ, beta_?MatrixQ, m0_?VectorQ] :=
	Module[
		{iterationlist = FixedPointList[
				{#[[1]](1 - Sign[beta . (1 - #[[2]])]),#[[2]]} & @ {#[[1]], #[[2]] Sign[#[[1]] . alpha]} &,
						{
						 MapAt[# - 1 &, Array[1 &, Length[alpha]], Transpose[{m0}]],
						 Array[1 &, Length[First[alpha]]]
						}
						]
		},
		MapAt[# - 1 &, (Last[iterationlist] /. (1 -> Infinity)) + Total[iterationlist], 2]
	];

Options[VolpertIndexing] := {ExternalSpecies -> {}, Verbose -> False};

VolpertIndexing[{reactions__}, initspecies_?VectorQ, specsreacs_List : {}, OptionsPattern[]]:=
	Module[{
			rdata, m, r, sp, reacs, vggraph, vgraph, w, w0, indexed, complx, initcomplx,
			indexedr, indexedsp, a, b, v, req, z1, z2, t1, t2, vindeces
		   },

		rdata = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["M","R","species","fhjgraphedges","volpertgraph"],Return[$Failed],
							     {ReactionsData::wrreac,ReactionsData::badarg}];

		If[rdata === $Failed, Return[$Failed], {m, r, sp, reacs, vggraph} = rdata; ];
 
		complx = Transpose[reacs /. Rule->List];
		initcomplx = First[complx];

		If[MemberQ[Flatten[complx],"0"], sp = Append[sp,"0"]; m = m+1; ];
		If[MemberQ[initcomplx,"0"], w0 = {m};, w0 = {};];

		w = Flatten[Position[sp, #] &/@ DeleteDuplicates[initspecies]];
		If[w =!= {} || initspecies === {},
	
			reacs = reacs /. Rule->RightArrow;
			indexedr = Last[vggraph];
			indexedsp = MapIndexed[#1->First[#2]&, sp];
			indexed = Join[indexedr,indexedsp];
			vgraph = First[vggraph];

			a = Cases[vgraph, {Rule[x_,y_RightArrow],z_}:> Rule[{x,y} /. indexed, z]];
			b = Cases[vgraph, {Rule[x_RightArrow,y_],z_}:> Rule[{y,x} /. indexed, z]];

			v = VolpertIndexing[SparseArray[a, {m, r}], SparseArray[b, {m, r}], DeleteDuplicates[Join[w,w0]]];
			vindeces = MapThread[Thread[Rule[#1,#2]]&, {{sp, reacs}, v}];,

			Message[VolpertIndexing::"specerr",initspecies];
			Return[$Failed]

		];

		If[OptionValue[Verbose],

			req = Flatten[specsreacs];		
			If[req === {}, 
				z1 = sp; 
				z2 = reacs;,
				z1 = Intersection[sp, req];
				z2 = Intersection[reacs, req];
				If[ z1 === {}, 
					If[ z2 === {},
						Message[VolpertIndexing::"verbtrumis"]; 
						Return[$Failed], 
						z1 = sp;
					];,
					If[ z2 === {}, z2 = reacs; ];
				];
			];
			z2 = z2 /. Rule->RightArrow;
			t1 = Sort[Transpose[{z1, z1 /. First[vindeces]}], #1[[2]]<#2[[2]]&];
			t2 = Sort[Transpose[{z2, z2 /. Last[vindeces]}], #1[[2]]<#2[[2]]&];
			{
				Grid[Join[{{"Species","Indexes"}},t1],Dividers->{{False,True},{False,True}}],
				Grid[Join[{{"Reaction steps","Indexes"}},t2],Dividers->{{False,True},{False,True}}]
			},

			vindeces
		]
	];


GraphPlotFunction2Q := MemberQ[{"Graph","GraphPlot","LayeredGraphPlot","GraphPlot3D","TreePlot"},#]&;
(*"ShowGraph","ShowLabeledGraph"*)


ShowVolpertGraph::usage = "ShowVolpertGraph[{reactions},options] visualizes the Volpert graph of the reaction.";

Options[ShowVolpertGraph] := {ExternalSpecies->{}, PlotFunction->"GraphPlot", EdgeLabels -> False, Highlighted -> {}, Indexed -> False, Numbered -> False};

ShowVolpertGraph::badarg = "Illegal argument of function ShowVolpertGraph.";
ShowVolpertGraph::funcarg = "Argument `1` is not a valid graph plot function. Try Graph, \
GraphPlot, GraphPlot3D, LayeredGraphPlot or TreePlot.";

ShowVolpertGraph[{reactions__},opts___?OptionQ] := 
	Module[{
			fullopts, graphfunc, exps, vgp, vindices, vtempind, 
			vg, vind, complx, numbd, ind, initsp, highlighted, edgelabels,
			nopts, g, bipd, fhj, sp, m, r
		   },

		fullopts = Join[{opts},Options[ShowVolpertGraph]];		
		graphfunc = PlotFunction /. FilterRules[fullopts, PlotFunction];

		If[GraphPlotFunction2Q[graphfunc],

			exps = FilterRules[{opts}, ExternalSpecies];
			vgp = Check[ReactionsData[{reactions},exps]["volpertgraph"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
			
			If[ vgp === $Failed, Return[$Failed];, {vg, vind} = vgp;];

			highlighted = Flatten[{Highlighted /. FilterRules[fullopts, Highlighted]}]; (*ToCanonicalForm*)

			ind = Indexed /. FilterRules[fullopts, Indexed];
			sp = ReactionsData[{reactions},exps]["species"];
			initsp = Intersection[Flatten[{ind}],sp];
			If[ ind =!= False, 
				vindices = VolpertIndexing[{reactions},initsp] /. Rule[x_,y_] :> Rule[x,{x,y}];

				If[ind =!= True && initsp === {} && Flatten[{ind}] =!= {}, Message[VolpertIndexing::"specerr",Flatten[{ind}]];];
				vg = Cases[vg, {Rule[x_,y_],z_} :> If[Head[x]===RightArrow,
														{Rule[x /. Last[vindices], y /. First[vindices]],z},
														{Rule[x /. First[vindices], y /. Last[vindices]],z}] 
					 ];
				(*
				vg = vg /. x_RightArrow :> {x /. vtempind, Last[x /. Last[vindices]]};
				vg = vg /. First[vindices] /. (Reverse /@ vtempind);
				*)
				vtempind = vind /. Rule[x_,y_] :> Rule[x,Subscript[xyz,y]];
				highlighted = highlighted /. x_RightArrow :> {x /. vtempind, Last[x /. Last[vindices]]};
				highlighted = highlighted /. First[vindices] /. (Reverse /@ vtempind);
				(*If[Head[#]===RightArrow, {#, # /. Last[vindices]}, {#, # /. First[vindices]}] &/@ highlight;*)
			];

			numbd = (Numbered /. FilterRules[fullopts, Numbered]) === True; 
			If[ numbd, 
					vg = vg /. vind; 
					highlighted = highlighted /. vind;
			];

			edgelabels = (EdgeLabels /. FilterRules[fullopts, EdgeLabels]) === True;

			If[ highlighted =!= {},
				
				If[ edgelabels, vg = Labeled@@# &/@ vg;, vg = First /@ vg; ];
				(*ez itt erzeketlen a plotfunction-ra es a bipartite-re*)
				g = Graph[vg];
				highlighted = Subgraph[g,highlighted];
				HighlightGraph[g, highlighted, FilterRules[{opts},Options[HighlightGraph]]],
				
				nopts = FilterRules[{opts},Options[ToExpression[graphfunc]]];		
				Switch[graphfunc,
						"Graph",
							If[ edgelabels, vg = Labeled@@# &/@ vg;, vg = First /@ vg; ];
							(*ez meg itt erzeketlen a bipartite-re*)
							Graph[vg, nopts],
						_, 
							bipd = (GraphLayout /. FilterRules[fullopts, GraphLayout]) === "Bipartite";
							nopts = Join[(FilterRules[nopts, EdgeLabeling] /. {}->{EdgeLabeling->False}), nopts];
							If[ edgelabels, nopts = Join[{EdgeLabeling -> True}, nopts]; ];
							If[ bipd,
								fhj = ReactionsData[{reactions},exps]["fhjgraphedges"] /. Rule->RightArrow;
								{m, r} = ReactionsData[{reactions},exps]["M","R"];
								complx = ReactionsData[{reactions},exps]["complexes"];
								If[MemberQ[complx,"0"], sp = Append[sp, "0"]; m = m+1;];
								If[ ind =!= False, 
									sp = sp /. Flatten[vindices];
									fhj = fhj /. Flatten[vindices];
								];
								If[ numbd, fhj = fhj /. vind;];
								(*Print[{vg,sp,fhj,m,r,nopts}];*)

								GraphPlot[vg,
									EdgeRenderingFunction -> (
										If[
											MemberQ[{RightArrow,Integer},Head[If[VectorQ[#2[[1]]],#2[[1,1]],#2[[1]]]]], 
											{Orange, Arrow[#1,0.15], If[edgelabels, Inset[Cases[vg,{Rule@@#2,x_}:>x][[1]], Mean[#1], Background->White], Sequence[]]},
											{Blue, Arrow[#1,0.15], If[edgelabels, Inset[Cases[vg,{Rule@@#2,x_}:>x][[1]], Mean[#1], Background->White], Sequence[]]}
										]&),
									VertexCoordinateRules -> 
										Join[
											Table[(sp[[i]])->{0,i},{i,1,m,1}],Table[(fhj[[i]])->{Ceiling[Min[{m,r}]/2],i},{i,1,r,1}]
										], nopts
								],
								ToExpression[graphfunc][vg, nopts]
							]
					]
			  ],
			Message[ShowVolpertGraph::"funcarg",graphfunc];
			Return[$Failed]
		]
	];

ShowVolpertGraph[___] := (Message[ShowVolpertGraph::"badarg"]; $Failed)


AcyclicVolpertGraphQ::usage = "AcyclicVolpertGraphQ[reactions,options] checks whether the volpert graph of the reaction is acyclic.";

AcyclicVolpertGraphQ::badarg = "Illegal argument of function AcyclicVolpertGraphQ.";

Options[AcyclicVolpertGraphQ] := {ExternalSpecies -> {}};

AcyclicVolpertGraphQ[{reactions__}, externals:(_?VectorQ|OptionsPattern[])] := 
	AcyclicGraphQ[Graph[First/@First[Check[ReactionsData[{reactions},externals]["volpertgraph"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]]]];

AcyclicVolpertGraphQ[___] := (Message[AcyclicVolpertGraphQ::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*SCLGraph, ShowSCLGraph*)


ShowSCLGraph::usage = "Given a reactionData returns the SCL graph of the mechanism. The linkage classes are denoted by \!\(\*SubscriptBox[\(L\), \(i\)]\) where i is theindex of the linkage class in\
 WeaklyConnectedGraphComponents[reactionData[\"fhjgraphedges\"]]";
SCLGraph[reactionData_,opts___?OptionQ] :=
 Block[{linkageClasses = WeaklyConnectedGraphComponents[Graph[reactionData["fhjgraphedges"]]],species=reactionData["species"],edges, edgeLabels},
 {edges, edgeLabels} = Flatten[
         Function[complex, Function[species,{#[[1]]->species,complex}] /@ complextospecies[complex]]
         /@ DeleteDuplicates[Flatten[ReactionsToList[ToCanonicalForm[EdgeList[#[[2]]]]]/.RightArrow->List]] &
         /@ ({ Subscript[Global`L, #]& /@ Range[Length[linkageClasses]],linkageClasses}//Transpose),
         2
 ]// Transpose;
 UndirectedGraph[edges,EdgeLabels->Table[edges[[i]]->edgeLabels[[i]],{i, Length[edges]}]]
]


ShowSCLGraph::usage = "Given a ReactionData[reactions] draws the SCL graph of the mechanism. The linkage classes are denoted by numbers. The number of a vertex represents\
 the index of its linkage class in WeaklyConnectedGraphComponents[reactionData[\"fhjgraphedges\"]]";

ShowSCLGraph::badarg = "Illegal argument of function ShowSCLGraph.";

ShowSCLGraph[reactionData_,opts___?OptionQ]:=Block[{g=SCLGraph[reactionData,opts]},
 GraphPlot[g,opts,GraphLayout->"BipartiteEmbedding"]
]


(* ::Subsubsection::Closed:: *)
(*Atoms, ToAtomMatrix, AtomConservingQ, FromAtomMatrix*)


Atoms::usage = "Atoms returns the list of atoms from Ac through Uuh to Y.";

Atoms = 
	{
	  "Ubn","Ubu","Uub","Uue","Uuh","Uuo","Uup","Uuq","Uus","Uut",
	  "Ac","Ag","Al","Am","Ar","As","At","Au","Ba","Be","Bh","Bi","Bk","Br","Ca","Cd","Ce","Cf","Cl","Cm","Co","Cr","Cs","Cu",
	  "Db","Ds","Dy","Er","Es","Eu","Fe","Fm","Fr","Ga","Gd","Ge","He","Hf","Hg","Ho","Hs","In","Ir","Kr","La","Li","Lr","Lu",
	  "Md","Mg","Mn","Mo","Mt","Na","Nb","Nd","Ne","Ni","No","Np","Os","Pa","Pb","Pd","Pm","Po","Pr","Pt","Pu","Ra","Rb","Re","Rf","Rg","Rh","Rn","Ru",
	  "Sb","Sc","Se","Sg","Si","Sm","Sn","Sr","Ta","Tb","Tc","Te","Th","Ti","Tl","Tm","Xe","Yb","Zn","Zr",
	  "B","C","F","H","I","K","N","O","P","S","U","V","W","Y"
	};
(* sorrend! *)


NestedMoleculeQ[s_?StringQ] := StringFreeQ[NestWhile[StringReplace[#,("("~~x___~~")"/;StringCount[x,{"(",")"}]===0)->""]&,s,UnsameQ[##]&,2],{"(",")"}];

MoleculeToAtoms[s_?StringQ,ats_?VectorQ]:=
	StringCases[s, 
			Join[" "~~r1___?DigitQ~~#:>ToExpression[#<>If[r1=="","1",r1]]*"charge"&/@{"+","-"},
			#~~r2___?DigitQ:>ToExpression[If[r2 == "","1",r2]]*#&/@ats,
			{"("~~x___?NestedMoleculeQ~~")"~~r3__?DigitQ:>Sequence@@Times[ToExpression[If[r3=="","1",r3]],MoleculeToAtoms[x,ats]]}],
				Overlaps->False];


ToAtomMatrix::usage = "ToAtomMatrix[molecules,options] constructs the atomic matrix using the molecular structure of the species. \
It prints a table if the option FormattedOutput is set to be True.";
(*Checks if the atoms are conserved in the reaction steps and list those steps where they are not*)

ToAtomMatrix::badarg = "Illegal argument of function ToAtomMatrix.";

Options[ToAtomMatrix] := {FormattedOutput->False};

ToAtomMatrix[molecules_?VectorQ,opts___?OptionQ]:=
	Module[{ats, l, ind, atomm, fout},

		ats = Append[Reverse[Sort[Union@@StringCases[molecules,Atoms],StringLength[#1]<StringLength[#2]&]], "charge"];
		l = Length[ats];
		atomm = Transpose[Map[Total[MoleculeToAtoms[#,ats]]&,molecules]/.Thread[ats->(UnitVector[l,#]&/@Range[l])]];
		fout = (FormattedOutput /. Flatten[{opts}]) /. FormattedOutput -> False;

		If[fout,
			Grid[MapThread[Prepend,{Prepend@@({atomm, molecules}),Prepend[ats,"Atoms\\Molecules"]},1], FilterRules[{opts}, Options[Grid]]],

			ind = Sort[MapIndexed[{#1,First[#2]}&, Most[ats]]];
			Join[{Join[First/@ind,{"charge"}]}, {Join[atomm[[#]]& /@ (Last/@ind), {Last[atomm]}]}]
		]
	];

ToAtomMatrix[___] := (Message[ToAtomMatrix::"badarg"]; $Failed)


AtomConservingQ::usage = "AtomConservingQ[{reactions},options] returns True if and only if the \
number of atoms together with the charges is conserved in each reaction step.";

AtomConservingQ::badarg = "Illegal argument of function AtomConservingQ.";

Options[AtomConservingQ]:={ExternalSpecies->{}};

SyntaxInformation[AtomConservingQ]={"ArgumentsPattern"->{{__},___}};

AtomConservingQ[{reactions__},externals:(_?VectorQ|OptionsPattern[])] :=
	Module[{species, gamma, atommatrix, nofatch},
		{species, gamma} = Check[ReactionsData[{reactions},externals]["species","\[Gamma]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		atommatrix = Last[ToAtomMatrix[species]];
		nofatch = Length[atommatrix];

		(#===ConstantArray[0,nofatch]) &/@ Transpose[atommatrix . gamma]
	];

AtomConservingQ[___] := (Message[AtomConservingQ::"badarg"]; $Failed)


AtomsQ::usage = "AtomsQ[list] returns True if and only if each element of list is an atom \
listed in Atoms.";

AtomsQ::badarg = "Illegal argument of function AtomsQ.";

SyntaxInformation[AtomsQ] = {"ArgumentsPattern"->{_}};

AtomsQ[v_?VectorQ] := And@@(MemberQ[Atoms,#]&/@v);

AtomsQ[___] := (Message[AtomsQ::"badarg"]; $Failed)


FromAtomMatrix::usage = "FromAtomMatrix[atommatrix,atoms] constructs the molecular structure from an atomic matrix using the atoms.";

FromAtomMatrix::badarg = "Illegal argument of function FromAtomMatrix.";

SyntaxInformation[FromAtomMatrix]={"ArgumentsPattern"->{_,_}};

FromAtomMatrix[atommatrix_?MatrixQ,atoms_?AtomsQ] := 
	Inner[
		Switch[{#1,#2},{0,_},"",{_,"charge"},StringJoin[" ",If[Abs[#1]===1,"",ToString[Abs[#1]]],If[#1<0,"-","+"]],{1,_},#2,_,StringJoin[#2,ToString[#1]]]&,
		Transpose[atommatrix],
		Append[atoms,"charge"],
		StringJoin
	];

FromAtomMatrix[___] := (Message[FromAtomMatrix::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*DetailedBalanced*)


USet[a_Integer,b_Integer,s_List] :=
	Module[{sa = FSet[a,s], sb = FSet[b,s], set = s},
	      If[ set[[sa,2]] < set[[sb,2]], {sa,sb} = {sb,sa} ];
	        set[[sa]] = {sa, Max[ set[[sa,2]], set[[sb,2]]+1 ]};
	        set[[sb]] = {sa, set[[sb,2]]};
	        set
	];

FSet[n_Integer,s_List] := 
        Block[{$RecursionLimit = Infinity}, 
              If [n == s[[n,1]], n, FSet[s[[n,1]],s]]
        ];

MST[e_List,n_Integer]:=
	Module[{sorte, i, s},
			sorte = Sort[e, (#1[[2]] <= #2[[2]])&];
			s = Table[{i,1},{i,n}];
			Select[sorte,(If[FSet[#[[1]],s]!=FSet[#[[2]], s], s = USet[#[[1]],#[[2]], s]; True,False])&]
	];


DetailedBalanced::usage = "DetailedBalanced[{reactions},rrcoeffs,options] gives the necessary and sufficient conditions \
for the given reaction to be detailed balanced using rrcoeffs as reaction rate coefficients.";

DetailedBalanced::badarg = "Illegal argument of function DetailedBalanced.";
DetailedBalanced::genparerr = "The option GeneratedRateCoefficient must be assigned a symbol.";
DetailedBalanced::notrev = "The reaction is not reversible.";
DetailedBalanced::defeqz = "The reaction has deficiency zero so only the circuit condition(s) are taken into account.";
DetailedBalanced::nocycle = "The FHJ-graph of the reaction has no cycle so only the spanning forest condition(s) are taken into account.";
DetailedBalanced::nocycdefz = "The reaction has deficiency zero and its FHJ-graph has no cycle so the reaction is detailed balanced.";
DetailedBalanced::rates = "The number of reaction steps does not match with that of reaction rate coefficients: `1`.";
DetailedBalanced::misscirc = "We may have missing circuit condition(s).";
DetailedBalanced::missspf = "We may have missing spanning forest condition(s).";

Options[DetailedBalanced] := {ExternalSpecies -> {}, GeneratedRateCoefficient -> Global`k};

SyntaxInformation[DetailedBalanced]={"ArgumentsPattern"->{{__},___}};

DetailedBalanced[{reactions__}, opts : OptionsPattern[]] :=
	Block[ { genpar, rr },

		genpar = GeneratedRateCoefficient /. Flatten[{opts, Options[DetailedBalanced]}];

		If[ Head[genpar] === Symbol,

			rr = Check[ReactionsData[{reactions}, ExternalSpecies /. Flatten[{opts, Options[DetailedBalanced]}]]["R"], Return[$Failed],
							{ReactionsData::wrreac,ReactionsData::badarg}];

			(*If[rr === $Failed, Return[$Failed]];*)

			DetailedBalanced[{reactions}, Array[Subscript[genpar,#]&,rr], opts],

			Message[DetailedBalanced::genparerr];
			Return[$Failed];

		]
	];

DetailedBalanced[{reactions__}, rates_?VectorQ, opts : OptionsPattern[]] := 
	Module[{ exs, cxes, m, r, nn, ll, fhj, gamma, deficiency, fhjunlist, fhjun, rsteplistind, rrates, shadow, shw,
			 fhjun0, mintree, mintreeun, c0, circs, spanreacs, mm, ccc, circconds, spanforconds, rev
		   },

			exs = ExternalSpecies /. Flatten[{opts, Options[DetailedBalanced]}];

			If[Not[ReversibleQ[{reactions},exs]],
				Message[DetailedBalanced::"notrev"];
				Return[$Failed];
			];

			{cxes, m, r, nn, ll, fhj, gamma, deficiency} = 
					Check[ReactionsData[{reactions},exs]["complexes","M","R","N","L","fhjgraphedges","\[Gamma]","deficiency"], Return[$Failed], 
						 {ReactionsData::wrreac,ReactionsData::badarg}];

			If[Length[rates] =!= r,
				Message[DetailedBalanced::"rates",rates];
				Return[$Failed];
			];
			
			deficiency = ToExpression[StringReplace[deficiency,Longest[__]~~"="->""]];
							(*ToExpression[Last[StringSplit[deficiency,"="]]];*)

			fhj = fhj /. Sort[Thread[cxes->Range[nn]], Length[#1[[1]]/.Plus->List]>Length[#2[[1]]/.Plus->List]&];

			fhjunlist = Union[Sort /@ (fhj /. Rule->UndirectedEdge)];
			fhjun = Graph[fhjunlist];
			rsteplistind = fhj /. Rule->List;
			shadow = Array[Subscript[shw,#]&,r];
			rrates = Thread[rsteplistind -> shadow];
			(**)
			fhjun0 = fhjunlist /. UndirectedEdge->List;
			mintreeun = MST[fhjun0,nn];
			mintree = Graph[UndirectedEdge@@# &/@ mintreeun];
			(**)
			rev = Thread[shadow -> rates];
			(**)
			circconds := (
							If[Length[c0 = Complement[fhjun0, mintreeun]] =!= r/2-nn+ll,
								Message[DetailedBalanced::"misscirc"];
							];
							circs = Map[Partition[Append[FindShortestPath[mintree,#[[1]],#[[2]]],#[[1]]],2,1]&, c0];

							Map[
								(Times @@ (# /. rrates) == Times @@ (Map[Reverse,#] /. rrates)) /. rev
								&, circs
							]
						 );
			(**)
			spanforconds := (
							spanreacs = Flatten[Position[rsteplistind, #]& /@ mintreeun];
							If[Length[mm = NullSpace[Part[gamma, All, spanreacs]]] =!= deficiency,
								Message[DetailedBalanced::"missspf"];
							];

							Map[
								Numerator[Together[
											Times @@ ((mintreeun /. rrates)^#) + Times @@ ((Map[Reverse,mintreeun] /. rrates)^#)
										 ]
								] /. Plus -> Equal /. rev				
								&, mm
							]
							);

			If[ deficiency === 0,

				If[AcyclicGraphQ[fhjun],

					Message[DetailedBalanced::"nocycdefz"];
					Return[True];,

					Message[DetailedBalanced::"defeqz"];
					circconds
				],

				If[AcyclicGraphQ[fhjun],

					Message[DetailedBalanced::"nocycle"];
					spanforconds,

					Join[circconds,spanforconds]	
				]
			]
	];

DetailedBalanced[___] := (Message[DetailedBalanced::"badarg"]; Return[$Failed];)


(* ::Subsection::Closed:: *)
(*deterministic approach*)


(* ::Subsubsection::Closed:: *)
(*MassActionKinetics, RightHandSide, DeterministicModel*)


MassActionKinetics::usage = "MassActionKinetics[{reactions},rratecoeffs,vars] builds up the mass action type kinetics attached to the reaction using \
rratecoeffs as reaction rate coefficients and vars as the names of independent variables of the kinetic function.";

MassActionKinetics::args = "The number of variables or that of reaction rate coefficients do not match with the dimensions of \[Alpha].";
MassActionKinetics::badarg = "Illegal argument of function MassActionKinetics.";

Options[MassActionKinetics] := {ExternalSpecies->{}};

MassActionKinetics[{reactions__}, rratecoeffs_?VectorQ, OptionsPattern[]] := 
	MassActionKinetics[
			Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["\[Alpha]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],
			rratecoeffs,
			Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]];

MassActionKinetics[{reactions__}, rratecoeffs_?VectorQ, variables_?VectorQ, OptionsPattern[]] := 
	MassActionKinetics[
			Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["\[Alpha]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],
			rratecoeffs,
			variables];

MassActionKinetics[alpha_?MatrixQ,rratecoeffs_?VectorQ,variables_?VectorQ] := 
	Module[{s},
		s = Dimensions[alpha];
		If[Length[rratecoeffs]===Last[s]&&Length[variables]===First[s],
			rratecoeffs*((Times @@ Power[variables,#]) &/@ Transpose[alpha]),
			Message[MassActionKinetics::args];
			$Failed
		]
	];

MassActionKinetics[___] := (Message[MassActionKinetics::"badarg"]; $Failed)


RightHandSide::usage = "RightHandSide[{reactions},rates,vars,options] returns the right-hand side of the \
induced kinetic differential equation assigned to the reaction using vars as the names of independent variables and rates as the kinetic function. \
MassAction \[Rule] True returns the right-hand side of the induced kinetic differential equation endowed with mass action type kinetics.";

RightHandSide::args = "The number of variables or that of rates do not match with the dimensions of the \[Alpha] matrix.";
RightHandSide::badarg = "Illegal argument of function RightHandSide.";

Options[RightHandSide] := {ExternalSpecies->{},MassAction->True};

RightHandSide[{reactions__}, rates_?VectorQ, opts : OptionsPattern[]] := 
	RightHandSide[{reactions},rates,Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],opts];

RightHandSide[{reactions__}, rates_?VectorQ, vars_?VectorQ, OptionsPattern[]] := 
	Module[{alpha, gamma, species, nofr, nofs},
		{alpha, gamma, species, nofr, nofs} = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["\[Alpha]","\[Gamma]","species","R","M"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		If[Length[vars]===nofs && Length[rates]===nofr,
			If[OptionValue[MassAction],
				Dot[gamma, MassActionKinetics[alpha,rates,vars]],
				Dot[gamma, rates /. Thread[species -> vars]]
			],
			Message[RightHandSide::args];
			$Failed
		]
	];

RightHandSide[___] := (Message[RightHandSide::"badarg"]; $Failed)


DeterministicModel::usage = "DeterministicModel[{reactions},rates,vars,options] returns the induced kinetic differential equation \
assigned to the reaction using rates as the kinetic function and vars as the names of the independent variables (both are optional arguments). \
MassAction \[Rule] True yields the induced kinetic differential equation endowed with mass action type kinetics. \n\
DeterministicModel[{reactions},arrheniuscoeffs,enthalpies,c0,V,Ta,tres,p1,p2,options] returns a general evolution equation for reactions involving temperature effects. \
The kinetics is assumed to be of the mass action type, the temperature dependence is allowed to be of the generalized Arrhenius type. \
The Arrhenius triples are given by arrheniuscoeffs.";

DeterministicModel::badarg = "Illegal argument of function DeterministicModel.";
DeterministicModel::args = "The number of rates or that of enthalpiesis are unequal to the number of reaction steps.";
DeterministicModel::steady = "The given set of steady state concentrations do not match with the number of variables.";

Options[DeterministicModel] := {ExternalSpecies->{}, MassAction->True};

DeterministicModel[{reactions__},s_Symbol:Global`t, opts : OptionsPattern[]] :=
	DeterministicModel[{reactions},Array[Subscript[Global`k,#]&,
		Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]],s,opts];

DeterministicModel[{reactions__},rates_?VectorQ,s_Symbol:Global`t, opts : OptionsPattern[]] :=
	DeterministicModel[{reactions},rates,
		Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],s,opts];	

DeterministicModel[{reactions__},rates_?VectorQ,vars_?VectorQ,s_Symbol:Global`t, opts : OptionsPattern[]] :=
	Module[{v},
		v = Through[vars[s]];		
		{Thread[D[v,s] == Check[RightHandSide[{reactions},rates,v,opts],Return[$Failed];,{RightHandSide::args,RightHandSide::badarg}]], v}
	];

(*Arrhenius reaction rate coefficients*)
DeterministicModel[{reactions__},arrhenius_?MatrixQ,enthalpies_?VectorQ,C0_?VectorQ,V_,Ta_,tres_,p1_,p2_,s_Symbol:Global`t,OptionsPattern[]] :=
	Module[{kkk, maf, alpha, beta, vars, reacs},
		reacs = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		If[reacs===Length[arrhenius] && reacs===Length[enthalpies],
			{alpha, beta, vars} = ReactionsData[{reactions},OptionValue[ExternalSpecies]]["\[Alpha]","\[Beta]","variables"];
			vars = Through[vars[s]];
			If[Length[C0]===Length[vars],
				kkk[{a_,n_,e_}] := a Global`T[s]^n Exp[- e /( Global`T[s] MolarGasConstant )];
				maf = MassActionKinetics[alpha,kkk/@arrhenius,vars];
				{Join[Thread[D[vars,s] == (Dot[beta-alpha,maf] + (C0-vars)/tres)], {D[Global`T[s],s] == - Dot[maf,enthalpies]/p1 - (1/tres + p2/p1)(Global`T[s]-Ta)/V}], Join[vars,{Global`T[s]}]},
				Message[DeterministicModel::"steady"];
				$Failed
			],
			Message[DeterministicModel::"args"];
			$Failed
		]
	];

DeterministicModel[___] := (Message[DeterministicModel::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*Concentrations*)


Concentrations::usage = "Concentrations[{reactions},rates,initvalues,vars,options] attempts to solve the \
induced kinetic differential equation endowed with some kinetics described by the rates. Initial concentrations are given by initvalues, \
vars contain the names of the independent variables. This way Concentrations uses DSolve, accepts symbolic parameters, and \
tries to return the symbolic solution to the equation (i.e. concentrations versus time curves). \n\
Concentrations[{reactions},rates,initvalues,{t0,t1},vars,options] uses NDSolve to give the numerical solution from t0 to t1. \n\
Concentrations[{reactions},arrheniuscoeffs,enthalpies,c0,V,Ta,tres,p1,p2,initvalues,{t0,t1},options] tries to solve the general evolution equation for reactions \
involving temperature effects. In this case the kinetics is assumed to be of the mass action type, the temperature dependence is allowed to be of the \
generalized Arrhenius type. The Arrhenius triples are given by arrheniuscoeffs.";

Options[Concentrations] := {ExternalSpecies->{},MassAction->True};

Concentrations::badarg = "Illegal argument of function Concentrations.";
Concentrations::init = "The initial vector has wrong size.";

(* Solving the ODE symbolically *)
Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,s_Symbol:Global`t,opts___?OptionQ] := 
	Concentrations[{reactions},rates,init,0,s,opts];

Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,vars_?VariablesQ,s_Symbol:Global`t,opts___?OptionQ] := 
	Concentrations[{reactions},rates,init,0,vars,s,opts];

Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,t0_?NumericQ(*Symbol*),s_Symbol:Global`t,opts___?OptionQ] :=
	Module[{kineq, vars, exs},
		exs = Sequence@@FilterRules[{opts},Options[DeterministicModel]];
		{kineq, vars} = Check[DeterministicModel[{reactions},rates,s,exs],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===ReactionsData[{reactions},exs]["M"],
			{
				vars,
				First[DSolve[Join[kineq, Thread[(vars /. s -> t0) == init]], vars, s, FilterRules[{opts},Options[DSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,t0_?NumericQ(*Symbol*),vars_?VectorQ,s_Symbol:Global`t,opts___?OptionQ] :=
	Module[{kineq, newvars, exs},
		exs = Sequence@@FilterRules[{opts},Options[DeterministicModel]];
		{kineq, newvars} = Check[DeterministicModel[{reactions},rates,vars,s,exs],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===ReactionsData[{reactions},exs]["M"],
			{
				newvars,
				First[DSolve[Join[kineq, Thread[(newvars /. s -> t0) == init]], newvars, s, FilterRules[{opts},Options[DSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

(* Solving the ODE numerically *)
Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,{t0_?NumericQ,t1_?NumericQ},s_Symbol:Global`t,opts___?OptionQ] :=
	Module[{kineq, vars, exs},
		exs = Sequence@@FilterRules[{opts},Options[DeterministicModel]];
		{kineq, vars} = Check[DeterministicModel[{reactions},rates,s,exs],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===ReactionsData[{reactions},exs]["M"],
			{ 
				vars,
				First[NDSolve[Join[kineq, Thread[(vars /. s -> t0) == init]], vars, {s, t0, t1}, FilterRules[{opts},Options[NDSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,{t0_?NumericQ,t1_?NumericQ},vars_?VectorQ,s_Symbol:Global`t,opts___?OptionQ] :=
	Module[{kineq, newvars, exs},
		exs = Sequence@@FilterRules[{opts},Options[DeterministicModel]];
		{kineq, newvars} = Check[DeterministicModel[{reactions},rates,vars,s,exs],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===ReactionsData[{reactions},exs]["M"],
			{ 
				newvars,
				First[NDSolve[Join[kineq, Thread[(newvars /. s -> t0) == init]], newvars, {s, t0, t1}, FilterRules[{opts},Options[NDSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

(*Arrhenius reaction rate coefficients*)
Concentrations[{reactions__},arrhenius_?MatrixQ,enthalpies_?VectorQ,C0_?VectorQ,V_,Ta_,tres_,p1_,p2_,init_?VectorQ,{t0_?NumericQ,t1_?NumericQ},s_Symbol:Global`t,opts___?OptionQ] :=
	Module[{kineq, vars, exs},
		exs = Sequence@@FilterRules[{opts},Options[DeterministicModel]];
		{kineq, vars} = Check[DeterministicModel[{reactions},arrhenius,enthalpies,C0,V,Ta,tres,p1,p2,s,exs],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===(ReactionsData[{reactions},exs]["M"]+1),
			{ 
				vars,
				First[NDSolve[Join[kineq, Thread[(vars /. s -> t0) == init]], vars, {s, t0, t1}, FilterRules[{opts},Options[NDSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

Concentrations[___] := (Message[Concentrations::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*StationaryPoints, GammaLeftNullSpace, MassConservationRelations, EigensystemJacobian*)


StationaryPoints::usage = "StationaryPoints[{reactions},rates,initvalues,vars] returns the stationary points, found in a reaction simplex, \
of the induced kinetic differential equation endowed with some kinetics described by the rates.";

Options[StationaryPoints] := {ExternalSpecies->{}, Positivity->False, Conditions->{}, Method->"Automatic", MassAction->True};

StationaryPoints::badarg = "Illegal argument of function StationaryPoints.";
StationaryPoints::args = "One of the arguments `1`, `2` or `3` may have wrong shape.";
StationaryPoints::noeql = "The reaction simplex contains no stationary point.";
StationaryPoints::unknowm = "The system could not recognize the method `1`. Try Automatic or GammaLeftNullSpace.";

StationaryPoints[{reactions__},opts___?OptionQ] := 
	Module[{r, vars},
		{r, vars} = Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["R","variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		StationaryPoints[{reactions},Array[Subscript[Global`k,#]&,r],Superscript[#,0]&/@vars,Superscript[#,"*"]& /@ vars,opts]
	];

StationaryPoints[{reactions__},rates_?VectorQ,opts___?OptionQ] := 
	Module[{vars},
		vars = Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		StationaryPoints[{reactions},rates,Superscript[#,0]&/@vars,Superscript[#,"*"]& /@ vars,opts]
	];

StationaryPoints[{reactions__},rates_?VectorQ,init_?VectorQ,opts___?OptionQ] := 
	StationaryPoints[{reactions},rates,init,Check[Superscript[#,"*"]& /@ ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],opts];

StationaryPoints[{reactions__},rates_?VectorQ,init_?VectorQ,vars_?VectorQ,opts___?OptionQ]:=
	Module[{equilibriums, ropts, gamma, m, r, species, method, y, z, statp, balances, rateconst, inimass, nonneg, sol, fless, null, nullspacecond},

			ropts = Sequence@@FilterRules[{opts},Options[Reduce]];
			{gamma, m, r, species} = Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["\[Gamma]","M","R","species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
			fless = If[(Positivity/.{opts})/.Positivity->False,Less,LessEqual];

			If[r===Length[rates] && m===Length[init] && m===Length[vars],
			
				Switch[(method = (Method/.{opts}) /. Method->"Automatic"),

					"Automatic",(*default method*)

					z = Subscript[y,ToString[#]]& /@ Range[r];

					{statp, balances, rateconst, inimass, nonneg} = 
						Thread/@{Check[RightHandSide[{reactions},rates,vars,FilterRules[{opts},Options[RightHandSide]]],Return[$Failed];,{RightHandSide::args,RightHandSide::badarg}] == 0, 
							gamma . z == vars-init, (rates /. Thread[species -> vars]) > 0, Cases[init,Except[0]]>0, fless[0,vars]};

					If[(sol = Reduce[Join[statp, balances, rateconst, inimass, nonneg, (Conditions/.{opts})/.Conditions->{}], Join[vars,z], Reals, Backsubstitution->True, ropts])===False,
						Message[StationaryPoints::noeql];
						Return[$Failed],

						Join[{vars},{Sort[Sort[MyToRules[Select[#/.And->List,Not[MemberQ[z,First[#]]]&]]]&/@Flatten[{sol/.Thread[rateconst->True]/.Thread[inimass->True]/.Thread[nonneg->True]/.Or->List}]]}]
					],

					"GammaLeftNullSpace",(*using the nullspace of gamma*)

					null = NullSpace[Transpose[gamma]] /. {}->{ConstantArray[0,m]};

					{statp, nullspacecond, rateconst, inimass, nonneg} = 
						Thread/@{Check[RightHandSide[{reactions},rates,vars,FilterRules[{opts},Options[RightHandSide]]],Return[$Failed];,{RightHandSide::args,RightHandSide::badarg}] == 0, 
							null . (vars-init) == 0, (rates /. Thread[species -> vars]) > 0, Cases[init,Except[0]]>0,fless[0,vars]};

					sol := Reduce[Join[statp, nullspacecond, rateconst, inimass, nonneg, (Conditions/.{opts})/.Conditions->{}], vars, Reals, Backsubstitution->True, ropts];

					If[sol === False,
						Message[StationaryPoints::noeql];
						Return[$Failed],

						Join[{vars},{Sort[Sort[MyToRules[# /. And->List]]& /@ Flatten[{sol/.Thread[rateconst->True]/.Thread[inimass->True]/.Thread[nonneg->True]/.Or->List}]]}]
					],

					_,

					Message[StationaryPoints::unknowm,method];
					Return[$Failed]
				],

				Message[StationaryPoints::args,rates,init,vars];
				Return[$Failed]
			]
	];

StationaryPoints[___] := (Message[StationaryPoints::"badarg"]; $Failed)


GammaLeftNullSpace::usage = "GammaLeftNullSpace[{reactions},init,vars,opts] gives the left null space of the stoichiometric matrix \[Gamma].";

Options[GammaLeftNullSpace] := {ExternalSpecies->{}};

GammaLeftNullSpace::badarg = "Illegal argument of function GammaLeftNullSpace.";
GammaLeftNullSpace::args = "One of the arguments `1` or `2` may have wrong shape.";
GammaLeftNullSpace::null = "Left null space of \[Gamma] is null-dimensional.";

GammaLeftNullSpace[{reactions__},opts___?OptionQ] :=
	Module[{lowerspecies},
		lowerspecies = ToLowerCase[Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]];
		GammaLeftNullSpace[{reactions}, ToExpression/@(#<>"0" &/@ lowerspecies), ToExpression/@lowerspecies, opts]
	];

GammaLeftNullSpace[{reactions__},init_?VectorQ,opts___?OptionQ] := 
	GammaLeftNullSpace[{reactions},init,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["variables"],Return[$Failed];,
													{ReactionsData::wrreac,ReactionsData::badarg}],opts];

GammaLeftNullSpace[{reactions__},init_?VectorQ,vars_?VectorQ,opts___?OptionQ] :=
	Module[{gamma, m, ns, conds, sol},
			{gamma, m} = Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["\[Gamma]","M"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

			If[m===Length[init] && m===Length[vars],

				ns = NullSpace[Transpose[gamma]];
				If[ns === {},
					Message[GammaLeftNullSpace::null];
					$Failed,

					Thread[ns . vars == ns . init]
				],

				Message[GammaLeftNullSpace::args,init,vars];
				$Failed
			]
	];

GammaLeftNullSpace[___] := (Message[GammaLeftNullSpace::badarg]; $Failed)


MassConservationRelations::usage = "MassConservationRelations[{reactions},vars] gives a set of conditions for a (strictly positive) mass vector \
showing whether the given reaction is mass conserving or not.";

Options[MassConservationRelations] := {ExternalSpecies->{}};

MassConservationRelations::badarg = "Illegal argument of function MassConservationRelations.";

MassConservationRelations[{reactions__},opts___?OptionQ] := 
	MassConservationRelations[{reactions}, Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["variables"],Return[$Failed];,
													{ReactionsData::wrreac,ReactionsData::badarg}] /. (Global`c -> Global`\[Rho]), opts];

MassConservationRelations[{reactions__},vars_?VectorQ,opts___?OptionQ] := 
	Module[{gamma, ropts},
			gamma = Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["\[Gamma]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
			ropts = FilterRules[{opts},Options[Reduce]];

			Reduce[Flatten[Thread/@{Transpose[gamma] . vars==0, vars>0}], vars, Reals, Backsubstitution->True, ropts] /. And -> List /. Or -> List
	];

MassConservationRelations[___] := (Message[MassConservationRelations::"badarg"]; $Failed)


EigensystemJacobian::usage = "EigensystemJacobian[{reactions},rates,vars] returns the eigensystem of the Jacobian matrix of the right-hand side of the \
induced kinetic differential equation endowed with some kinetics given by the rates.";

EigensystemJacobian::badarg = "Illegal argument of function EigensystemJacobian.";
EigensystemJacobian::args = "One of the arguments `1` or `2` may have wrong shape.";

Options[EigensystemJacobian] := {ExternalSpecies->{},MassAction->True};

EigensystemJacobian[{reactions__}, opts : OptionsPattern[]] := 
	Module[{r, vars},
		{r, vars} = Check[ReactionsData[{reactions},opts]["R","variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		EigensystemJacobian[{reactions},Array[Subscript[Global`k,#]&,r],vars,opts]
	];

EigensystemJacobian[{reactions__},rates_?VectorQ,opts : OptionsPattern[]] := 
	EigensystemJacobian[{reactions},rates,
			Check[ReactionsData[{reactions},opts]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],opts];

EigensystemJacobian[{reactions__},rates_?VectorQ,vars_?VectorQ, opts : OptionsPattern[]] := 
	Module[{m, r},
		{m, r} = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["M","R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
	
		If[r===Length[rates] && m===Length[vars],

			Eigensystem[D[Check[RightHandSide[{reactions},rates,vars,opts],Return[$Failed];,{RightHandSide::args,RightHandSide::badarg}],{vars}]],

			Message[EigensystemJacobian::args,rates,vars];
			$Failed
		]
	];

EigensystemJacobian[___] := (Message[EigensystemJacobian::badarg]; $Failed)


(* ::Subsubsection::Closed:: *)
(*AbsoluteConcentrationsRobustness*)


AbsoluteConcentrationRobustness::usage = "AbsoluteConcentrationRobustness[{reactions},options] checks the Shinar-Feinberg conditions for the reaction.";

AbsoluteConcentrationRobustness::badarg = "Illegal argument of function AbsoluteConcentrationRobustness.";
AbsoluteConcentrationRobustness::maybe = "`1`.";
AbsoluteConcentrationRobustness::timeexc = "Time spent on the calculations of a positive stationary point exceeded `1` seconds \
without returning any result (the deficiency is one and there are nonterminal complexes differing only in one species). Try also StationaryPoints.";

Options[AbsoluteConcentrationRobustness] := { TimeLimit -> 120, ExternalSpecies -> {} };

AbsoluteConcentrationRobustness[{reactions__},opts___?OptionQ] := 
	Module[{exs, time, delta, fhjstrong, fhjterminal, nonterms, acr},

		exs = Sequence @@ FilterRules[{opts},ExternalSpecies];
		time = (TimeLimit /. FilterRules[{opts},TimeLimit]) /. TimeLimit->120;
		{delta, fhjstrong, fhjterminal} = 
			Check[ReactionsData[{reactions},exs]["deficiency","fhjstronglyconnectedcomponents","fhjterminalstronglyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		If[ Last[StringSplit[delta,"="]]=!="1",

			Message[AbsoluteConcentrationRobustness::"maybe","Deficiency does not equal to 1"];
			Return[{reactions}],

			nonterms = Flatten[fhjstrong/.Thread[fhjterminal->Sequence[]]] /. "0"->0 /. Thread[exs->0];
			acr = Flatten[(nonterms-# &/@ nonterms) /. Plus[A_,B_]->0 /. 0->{}]; (*/. r_*A_:>A;*)

			If[ acr === {},

				Message[AbsoluteConcentrationRobustness::"maybe","Deficiency is one but there do not exist nonterminal complexes differing only in one species"];
				Return[{reactions}],

				Switch[TimeConstrained[Check[Length[StationaryPoints[{reactions},exs,Positivity->True,Sequence@@FilterRules[{opts},Conditions],Sequence@@FilterRules[{opts},Method]]]>=2,False,
								{StationaryPoints::unknowm,StationaryPoints::args,StationaryPoints::badarg,StationaryPoints::noeql}],time,"Abort"],
						"Abort",
							Message[AbsoluteConcentrationRobustness::"timeexc",time];
							Return[{reactions}],
						False,
							Message[AbsoluteConcentrationRobustness::"maybe","There are no positive stationary points."]; 
							Return[{reactions}],
						True,
							Union[acr /. r_*A_?StringQ :> A /. A_?StringQ*r_ :> A /. r_*A_?SpeciesQ :> A /. A_?SpeciesQ*r_ :> {A}](*!*)
				]	
			]
		]
	];

AbsoluteConcentrationRobustness[___][___] := (Message[AbsoluteConcentrationRobustness::"badarg"]; $Failed)


(* ::Subsection::Closed:: *)
(*stochastic approach*)


Intensity[rates_,X_,alphatr_] := (*Intensity[alphatr,X,rates] =*) rates * (Times @@ FactorialPower[X,#] &/@ alphatr);

INTENSITY = Compile[{{x,_Integer,1},{alphatr,_Integer,2}},

				Times @@ FactorialPower[x,#] &/@ alphatr,
				(*RuntimeAttributes->{Listable},*) 
				Parallelization->Automatic,
				RuntimeOptions->{"WarningMessages"->False,"CatchMachineIntegerOverflow"->False} (*"EvaluateSymbolically"*)
			];

ChangeUnits[V_,rates_,alphatr_] := rates*(AvogadrosNumber*V)^(1 - Total /@ alphatr);


(* ::Subsubsection::Closed:: *)
(*MasterEquation*)


MasterEquation::usage = "MasterEquation[{reactions},rratecoeffs,vars,options] releases the master equation of the induced kinetic Markov process \
endowed with stochastic mass action type kinetics with stochastic reaction rate coefficients given by rratecoeffs. \
The optional argument vars contains the names of the independent variables.";

Options[MasterEquation] := {ExternalSpecies->{}};

MasterEquation::badarg = "Illegal argument of function MasterEquation.";
MasterEquation::args = "Argument `1` or `2` has wrong shape.";

MasterEquation[{reactions__}, h_Symbol:Global`P, opts : OptionsPattern[]] := 
	MasterEquation[{reactions},
		Array[Subscript[Global`k,#]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["R"],Return[$Failed];,
												{ReactionsData::wrreac,ReactionsData::badarg}]], h, opts];

MasterEquation[{reactions__}, rates_?VectorQ, h_Symbol:Global`P, opts : OptionsPattern[]] := 
	MasterEquation[{reactions},rates,
		Prepend[Array[Subscript[Global`x, #]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["M"],Return[$Failed];,
														{ReactionsData::wrreac,ReactionsData::badarg}]],Global`t], h, opts];

MasterEquation[{reactions__}, rates_?VectorQ, vars_?VectorQ, h_Symbol:Global`P, OptionsPattern[]] := 
	Module[{alpha, gamma, m, r, alphatr, gammatrind, n},
			
			{alpha, gamma, m, r} = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["\[Alpha]","\[Gamma]","M","R"],Return[$Failed];,
												{ReactionsData::wrreac,ReactionsData::badarg}];	
			alphatr = Transpose[alpha];
			gammatrind = MapIndexed[{#1, First[#2]}&, Transpose[gamma]];

			If[ Length[rates] === r && Length[vars] === m+1,

				{
					D[h @@ vars, First[vars]] == 
					- (h @@ vars)*Total[FunctionExpand[Intensity[rates, Rest[vars], alphatr]]]
					+ Total[Map[
						h@@(vars-Prepend[First[#],0])*
							FunctionExpand[rates[[Last[#]]]*(Times@@FactorialPower[Rest[vars]-First[#],alphatr[[Last[#]]]])] &, gammatrind
							]
					  ]
				},
				
				Message[MasterEquation::args, rates, vars];
				Return[$Failed];
			  ]
	];

MasterEquation[___] := (Message[MasterEquation::badarg]; $Failed)


(* ::Subsubsection::Closed:: *)
(*StationaryProbabilityDistributionEquation*)


StationaryProbabilityDistributionEquation::usage = "StationaryProbabilityDistributionEquation[{reactions},rratecoeffs,vars,options] releases \
the system of equations for the stationary distribution of the induced kinetic Markov process endowed with stochastic mass action type kinetics \
with stochastic reaction rate coefficients given by rratecoeffs.";

Options[StationaryProbabilityDistributionEquation] := {ExternalSpecies->{}};

StationaryProbabilityDistributionEquation::badarg = "Illegal argument of function StationaryProbabilityDistributionEquation.";
StationaryProbabilityDistributionEquation::args = "Argument `1` or `2` has wrong shape.";

StationaryProbabilityDistributionEquation[{reactions__}, h_Symbol:Global`\[CapitalPi], opts : OptionsPattern[]] := 
	StationaryProbabilityDistributionEquation[{reactions},
		Array[Subscript[Global`k,#]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["R"],Return[$Failed];,
										  {ReactionsData::wrreac,ReactionsData::badarg}]], h, opts];

StationaryProbabilityDistributionEquation[{reactions__}, rates_?VectorQ, h_Symbol:Global`\[CapitalPi], opts : OptionsPattern[]] := 
	StationaryProbabilityDistributionEquation[{reactions},rates,
		Array[Subscript[Global`x, #]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["M"],Return[$Failed];,
										   {ReactionsData::wrreac,ReactionsData::badarg}]], h, opts];

StationaryProbabilityDistributionEquation[{reactions__}, rates_?VectorQ, vars_?VectorQ, h_Symbol:Global`\[CapitalPi], OptionsPattern[]] := 
	Module[{alpha, gamma, m, r, exs = OptionValue[ExternalSpecies]},
			
			{m, r} = Check[ReactionsData[{reactions},exs]["M","R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
			
			If[ Length[rates] === r && Length[vars] === m,

				(Last[#]==0/.(h[Global`t,q___]:>h[q])) &/@ MasterEquation[{reactions}, rates, Prepend[vars,Global`t], h, ExternalSpecies->exs],

				Message[StationaryProbabilityDistributionEquation::args,rates,vars];
				$Failed
			  ]
	];

StationaryProbabilityDistributionEquation[___] := (Message[StationaryProbabilityDistributionEquation::badarg]; $Failed)


(* ::Subsubsection::Closed:: *)
(*ProbabilityGeneratingFunctionEquation, SolveProbabilityGeneratingFunctionEquation*)


ProbabilityGeneratingFunctionEquation::usage = "ProbabilityGeneratingFunctionEquation[{reactions},rrcoeffs,vars,options] \
displays the partial differential equation for the probability generating function of \
the induced kinetic Markov process endowed with stochastic mass action type kinetics with stochastic reaction rate coefficients given by rrcoeffs.";

ProbabilityGeneratingFunctionEquation::badarg = "Illegal argument of function ProbabilityGeneratingFunctionEquation.";
ProbabilityGeneratingFunctionEquation::args = "Argument `1` or `2` has wrong shape.";

Options[ProbabilityGeneratingFunctionEquation] := {ExternalSpecies->{}};

ProbabilityGeneratingFunctionEquation[{reactions__}, h_Symbol:Global`g, opts : OptionsPattern[]] := 
	ProbabilityGeneratingFunctionEquation[{reactions},
		Array[Subscript[Global`k,#]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["R"],Return[$Failed];,
										  {ReactionsData::wrreac,ReactionsData::badarg}]], h, opts];

ProbabilityGeneratingFunctionEquation[{reactions__}, rates_?VectorQ, h_Symbol:Global`g, opts : OptionsPattern[]] := 
	ProbabilityGeneratingFunctionEquation[{reactions},rates,
		Prepend[Array[Subscript[Global`z, #]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["M"],Return[$Failed];,
												   {ReactionsData::wrreac,ReactionsData::badarg}]
					 ], Global`t], h, opts];

ProbabilityGeneratingFunctionEquation[{reactions__}, rates_?VectorQ, vars_?VectorQ, h_Symbol:Global`g, OptionsPattern[]] := 
	Module[{alpha, beta, m, r, tv = First[vars], zv = Rest[vars]},
			{alpha, beta, m, r} = Check[ReactionsData[{reactions},OptionValue[ExternalSpecies]]["\[Alpha]","\[Beta]","M","R"],Return[$Failed];,
									   {ReactionsData::wrreac,ReactionsData::badarg}];

			If[Length[rates] === r && Length[vars] === m+1,
				{alpha, beta} = Transpose /@ {alpha, beta};
				D[h @@ vars,tv]==Total[rates*MapThread[(Times@@(zv^#1)-Times@@(zv^#2))*Derivative[0,Sequence@@#2][h]@@vars&,{beta,alpha}]],

				Message[ProbabilityGeneratingFunctionEquation::args,rates,vars]; 
				Return[$Failed];
			  ]
	];

ProbabilityGeneratingFunctionEquation[___] := (Message[ProbabilityGeneratingFunctionEquation::badarg]; $Failed)


SolveProbabilityGeneratingFunctionEquation::usage = "SolveProbabilityGeneratingFunctionEquation[{reactions},rratecoeffs,init,vars,options] \
attempts to solve the partial differential equation with initial and boundary conditions for the probability generating function of the induced kinetic Markov process \
endowed with stochastic mass action type kinetics with stochastic reaction rate coefficients given by rratecoeffs. \n\
The default method (\"Built-in\") uses DSolve to obtain the symbolic solution for the initial value problem with boundary condition. \
The methods \"MatrixExponential\" and \"Characteristics\" can only be used for first order reactions where \
the characteristic equations are being solved. The method \"MatrixExponential\" only works for compartmental systems.";

SolveProbabilityGeneratingFunctionEquation::badarg = "Illegal argument of function SolveProbabilityGeneratingFunctionEquation.";
SolveProbabilityGeneratingFunctionEquation::args = "One of the arguments `1`, `2` or `3` has wrong shape.";
SolveProbabilityGeneratingFunctionEquation::methoderr = "Argument `1` is not a valid method. Try \"Built-in\", \"Characteristics\" or \"MatrixExponential\".";
SolveProbabilityGeneratingFunctionEquation::orderlone = "The order of the reaction is at most one, thus the method \"Characteristics\" can also be applied.";
SolveProbabilityGeneratingFunctionEquation::ordererr = "The given method `1` cannot be applied since the order of the reaction is greater than one. Try the method \"Built-in\".";
SolveProbabilityGeneratingFunctionEquation::betamaxlone = "All the entries of \[Beta] are at most one, thus the method \"MatrixExponential\" can also be applied.";
SolveProbabilityGeneratingFunctionEquation::betamaxerr = "There exists an entry in \[Beta] which is greater than one. Thus the method `1` cannot be applied. \
Try \"Characteristics\" or \"Built-in\".";

Options[SolveProbabilityGeneratingFunctionEquation] := {ExternalSpecies->{}, Method->"Built-in"};

SolveProbabilityGeneratingFunctionEquation[{reactions__}, h_Symbol:Global`g, opts___?OptionQ] := 
	SolveProbabilityGeneratingFunctionEquation[{reactions}, 
			Array[Subscript[Global`k,#]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["R"],Return[$Failed];,
													{ReactionsData::"wrreac",ReactionsData::"badarg"}]], 
			Array[Superscript[Subscript[Global`x, #],0]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["M"],Return[$Failed];,
															{ReactionsData::wrreac,ReactionsData::badarg}]], h, opts];

SolveProbabilityGeneratingFunctionEquation[{reactions__}, rates_?VectorQ, h_Symbol:Global`g, opts___?OptionQ] := 
	SolveProbabilityGeneratingFunctionEquation[{reactions}, rates, 
			Array[Superscript[Subscript[Global`x, #],0]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["M"],Return[$Failed];,
															{ReactionsData::wrreac,ReactionsData::badarg}]], h, opts];

SolveProbabilityGeneratingFunctionEquation[{reactions__}, rates_?VectorQ, init_?VectorQ, h_Symbol:Global`g, opts___?OptionQ]:=
	SolveProbabilityGeneratingFunctionEquation[{reactions}, rates, init, 
		Prepend[Array[Subscript[Global`z, #]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["M"],Return[$Failed];,
														{ReactionsData::wrreac,ReactionsData::badarg}]], Global`t], h, opts];

SolveProbabilityGeneratingFunctionEquation[{reactions__}, rates_?VectorQ, init_?VectorQ, vars_?VectorQ, h_Symbol:Global`g, opts___?OptionQ] := 
	Module[{ 
			 exs = FilterRules[{opts},ExternalSpecies], tv = First[vars], zv = Rest[vars], method,
			 rd, m, alpha, beta, alphatr, betatr, rsteporders, order, maxbeta, 
			 zrb, coeffarray, bb, Ab, u, ub, ubs, v, vb, cb, s, eqsol, homsol, inhomsol, pgfe
		   },
			rd = Check[ReactionsData[{reactions},exs],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

			If[(m = rd["M"]) === Length[init] && Length[rates] === rd["R"] && Length[vars] === m+1,

				method = Method /. Flatten[{opts, Options[SolveProbabilityGeneratingFunctionEquation]}];
				{alpha, beta} = rd["\[Alpha]","\[Beta]"];
				{alphatr, betatr} = Transpose /@ {alpha, beta};
				rsteporders = rd["reactionsteporders"];
				order = Max[rsteporders];
				maxbeta = Max[betatr];

				ub = Array[Subscript[u, #]&, m];
				zrb = rates * MapThread[(Times@@(ub^#1)-Times@@(ub^#2))&, {betatr, alphatr}];

				Switch[method,

					"Built-in",
						If[order <= 1,
							Message[SolveProbabilityGeneratingFunctionEquation::orderlone];
						];

						pgfe = (D[h @@ vars, tv] == Total[rates*MapThread[(Times@@(zv^#1)-Times@@(zv^#2))*Derivative[0,Sequence@@#2][h]@@vars&, {betatr, alphatr}]]);

						{h @@ vars, Flatten[DSolve[{pgfe, ((h @@ vars)/.{tv->0})==Times@@(Rest[vars]^init)},h @@ vars,vars,FilterRules[{opts},Options[DSolve]]]]},

					"Characteristics",
						If[order > 1,
							Message[SolveProbabilityGeneratingFunctionEquation::ordererr,method];
							Return[$Failed];,
							If[maxbeta <= 1,
								Message[SolveProbabilityGeneratingFunctionEquation::betamaxlone];
							]
						];

						vb = Array[Subscript[v, #]&, m];
						ubs = Through[ub[s]];
						eqsol = ubs /. 
										Flatten[
											DSolve[Join[Thread /@ {D[ubs, s] == -Total[(zrb /. Thread[ub -> ubs]) * alphatr], (ubs /. s->0)==zv}], ubs, s]
										];

						homsol = (zv /. Flatten[Quiet[Solve[Thread[eqsol==vb], zv]]]) /. Thread[vb -> zv];
						inhomsol = Integrate[((1 - Sign[rsteporders]) . zrb) /. Thread[ub -> homsol], {s, 0, tv}];

						{h @@ vars, h @@ vars -> (Times@@((homsol /. s->tv)^init)) * Exp[inhomsol]},

					"MatrixExponential",
						If[order > 1,
							Message[SolveProbabilityGeneratingFunctionEquation::ordererr,method];
							Return[$Failed];,
							If[maxbeta>1,
								Message[SolveProbabilityGeneratingFunctionEquation::betamaxerr,method];
								Return[$Failed];,
								Message[SolveProbabilityGeneratingFunctionEquation::orderlone];
							]
						];

						coeffarray = CoefficientArrays[alpha . zrb, ub];

						If[ZeroVectorQ[coeffarray],
							bb = ConstantArray[0, m];
							Ab = ConstantArray[0, {m,m}];,
							{bb, Ab} = coeffarray;
						];
						cb = Flatten[ub /. Quiet[Solve[Thread[Ab . ub == -bb], ub]]] /. Thread[ub -> 0];

						homsol = MatrixExp[Ab s] . (zv - cb) + cb;
						inhomsol = Integrate[((1 - Sign[rsteporders]) . zrb) /. Thread[ub -> homsol], {s, 0, tv}];

						{h @@ vars, h @@ vars -> (Times@@((homsol /. s->tv)^init)) * Exp[inhomsol]},

					_,
						Message[SolveProbabilityGeneratingFunctionEquation::methoderr,method];
						Return[$Failed];						
				],
				Message[SolveProbabilityGeneratingFunctionEquation::args,rates,init,vars]; 
				Return[$Failed];
			]
	];

SolveProbabilityGeneratingFunctionEquation[___] := (Message[SolveProbabilityGeneratingFunctionEquation::badarg]; $Failed)


(* ::Subsubsection::Closed:: *)
(*MomentEquations*)


NonNegativeIntegerNotZeroVectorQ := Max[#]>0 && MatchQ[#, {__Integer?NonNegative}]&;

MomentswithSpeciesQ := MatchQ[#,{_,_Integer?NonNegative}]&;

ExpectedSubscript[e_,newvars_,ylist_,vars_,var_] := 
	Block[{ coeffrules },
		coeffrules = CoefficientRules[Expand[FunctionExpand[Times@@FactorialPower[vars,ylist]]],vars];
		(coeffrules/.Rule[a_,b_]:>b) . Map[Subscript[e,ToString[#,StandardForm]][var]&,coeffrules/.Rule[a_,b_]:>Inner[Power,newvars,a,Times]]
	];


MomentEquations::usage = "MomentEquations[{reactions},a,rratecoeffs,options] displays the a^th moment equations of the induced kinetic Markov process endowed with \
stochastic mass action type kinetics with stochastic reaction rate coefficients given by rratecoeffs. CombinatorialMoments \[Rule] True calculates with combinatorial moments.";

MomentEquations::badarg = "Illegal argument of function MomentEquations.";
MomentEquations::args = "Argument `1` has wrong shape.";
MomentEquations::nonnegarg = "Argument `1` must be a non-zero vector having nonnegative integer entries.";
MomentEquations::speciesarg = "Argument `1` contains invalid internal species.";

Options[MomentEquations] := {CombinatorialMoments->False, ExternalSpecies->{}};

MomentEquations[{reactions__}, a_, opts : OptionsPattern[]] :=
	MomentEquations[{reactions},a,
		Array[Subscript[Global`k,#]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["R"],Return[$Failed];,
												{ReactionsData::"wrreac",ReactionsData::"badarg"}]], Global`\[DoubleStruckCapitalE], Global`t, opts];

MomentEquations[{reactions__}, a_, rates_?VectorQ, opts : OptionsPattern[]] := MomentEquations[{reactions}, a, rates, Global`\[DoubleStruckCapitalE], Global`t, opts];

(*{e_Symbol, var_Symbol} or e_Symbol:Global`E, optional*)
MomentEquations[{reactions__}, a_, e_Symbol, var_Symbol, opts : OptionsPattern[]] :=
	MomentEquations[{reactions},a,
		Array[Subscript[Global`k,#]&,Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["R"],Return[$Failed];,
												{ReactionsData::wrreac,ReactionsData::badarg}]], e, var, opts];

MomentEquations[{reactions__}, a_, rates_?VectorQ, e_Symbol, var_Symbol, OptionsPattern[]] := 
	Module[{ 
			species, m, l, sppos, a1, pgfe, vars, newvars, exp, es,
			exs = OptionValue[ExternalSpecies], combmom = OptionValue[CombinatorialMoments] 
		   },

			{species, m} = Check[ReactionsData[{reactions},exs]["species","M"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

			Switch[a, 
				{__?MomentswithSpeciesQ},
					If[(sppos = Flatten[Position[species,#]& /@ First[a1 = Transpose[a]]]) =!= {},
						MomentEquations[{reactions}, #, rates, e, var, ExternalSpecies->exs, CombinatorialMoments->combmom]& /@ ((UnitVector[m,#]& /@ sppos)*Last[a1]),
						Message[MomentEquations::speciesarg,a];
						Return[$Failed];
					],
				_?NonNegativeIntegerNotZeroVectorQ,
					l = Length[a];
					If[l === m,
						If[combmom,
							pgfe = Check[ProbabilityGeneratingFunctionEquation[{reactions},rates,ExternalSpecies->exs],Return[$Failed];,
											{ProbabilityGeneratingFunctionEquation::args,ProbabilityGeneratingFunctionEquation::badarg}];
							vars = First[pgfe] /. Derivative[x__][y_][z__]:>{z};
							newvars = ToString /@ species;
							D[pgfe,Sequence@@Thread[List[Rest[vars],a]]]/.Thread[Rest[vars]->1]
												/. Derivative[x_,y__][f_][z__]:>D[ExpectedSubscript[e,newvars,{y},Rest[vars],var],{var,x}]
												/. Global`g[Global`t,Sequence@@Array[1&,l]]->1 /. func_[Global`t]:>func[var],
							D[Subscript[e,ToString[Inner[Power,species,a,Times],StandardForm]][var],var] == 
								Total[(Times@@StirlingS2[a,#])*Last[MomentEquations[{reactions},#,rates,e,var,ExternalSpecies->exs,CombinatorialMoments->True]]& /@ Rest[Tuples[Prepend[Range[Max[a]],0],l]]]
						],
						Message[MomentEquations::args,a];
						Return[$Failed];
					],
				_,
					Message[MomentEquations::nonnegarg,a];
					Return[$Failed];
			]
	];

MomentEquations[___] := (Message[MomentEquations::badarg]; $Failed)


(* ::Subsection::Closed:: *)
(*stochastic simulation*)


MonitorMSG[prop_,pill_,species_] := 
		Column[{
				Row[{ProgressIndicator[prop]," " <> ToString[Ceiling[100*prop]]<>"%"}],
				"Current time: " <> ToString[First[pill],TraditionalForm] <> " units.",
				Grid[Prepend[Thread[List[ToString /@ species, Last[pill]]], {"Species","Current state"}],
					 Dividers->{{False,True}, {False,True}}]
			   }
		];


(* ::Subsubsection::Closed:: *)
(*DependencyGraph*)


DependencyGraph::usage = "DependencyGraph[{reactions},options] returns the dependency graph of the given reaction.";

DependencyGraph::badarg = "Illegal argument of function DependencyGraph.";
DependencyGraph::args = "Argument `1` has wrong shape.";

Options[DependencyGraph] := {ExternalSpecies->{}};

DependencyGraph[{reactions__},opts___?OptionQ] :=
	DependencyGraph[{reactions}, Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["\[Alpha]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}], opts];

DependencyGraph[{reactions__},rates_?VectorQ,opts___?OptionQ] :=
	Module[{species, r},
		{species, r} = Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["species","R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		If[Length[rates] === r,
			DependencyGraph[{reactions}, MapThread[(1 - Boole[SameQ[#1,#2]])&,{rates /. # -> 0, rates}] &/@ species, opts],

			Message[DependencyGraph::args, rates];
			Return[$Failed];
		]
	];

DependencyGraph[{reactions__},depend_?MatrixQ,opts___?OptionQ] :=
	Module[{m, r, gamma},
		{m, r, gamma} = Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["M","R","\[Gamma]"],Return[$Failed];,
											{ReactionsData::wrreac,ReactionsData::badarg}];
		If[Dimensions[depend] === {m,r},
			Sign[Transpose[Abs[Sign[gamma]]] . depend], (*this is the adjacency matrix of the dependency graph*)

			Message[DependencyGraph::args, depend];
			Return[$Failed];
		]
	];

DependencyGraph[__] := (Message[DependencyGraph::badarg]; $Failed);


ShowDependencyGraph::usage = "ShowDependencyGraph[{reactions},options] displays the dependency graph of the reaction.";

ShowDependencyGraph::badarg = "Illegal argument of function ShowDependencyGraph.";

Options[ShowDependencyGraph] := {ExternalSpecies->{}};

ShowDependencyGraph[{reactions__},opts___?OptionQ] :=
	ShowDependencyGraph[{reactions}, Check[DependencyGraph[{reactions},opts],Return[$Failed];,{DependencyGraph::badarg,DependencyGraph::args}], opts];

ShowDependencyGraph[{reactions__},rates_?VectorQ,opts___?OptionQ] :=
	ShowDependencyGraph[{reactions}, Check[DependencyGraph[{reactions},rates,opts],Return[$Failed];,{DependencyGraph::badarg,DependencyGraph::args}], opts];

ShowDependencyGraph[{reactions__},adjmxdepend_?MatrixQ,opts___?OptionQ] := 
	Module[{rsteps, graphrules},
		rsteps = Check[ReactionsData[{reactions},FilterRules[{opts},ExternalSpecies]]["reactionsteps"],Return[$Failed];,
								{ReactionsData::wrreac,ReactionsData::badarg}];
		graphrules = Cases[ArrayRules[adjmxdepend], Rule[{x_,y_},z_]/;(z=!=0) :> Rule[rsteps[[x]],rsteps[[y]]]];

		GraphPlot[graphrules, FilterRules[{opts},Options[GraphPlot]]]
	];

ShowDependencyGraph[__] := (Message[ShowDependencyGraph::badarg]; $Failed);


(* ::Subsubsection::Closed:: *)
(*Direct and first reaction methods*)


MyExponentialDistribution[0.`] := Infinity;
MyExponentialDistribution[0] := Infinity;
MyExponentialDistribution[x_] := RandomReal[ExponentialDistribution[x]];

SetAttributes[MyExponentialDistribution,Listable];


(*$IterationLimit, $RecursionLimit*)
(* t0 *)
DirectMethod[species_,alphatr_,gammatr_,rates_,X0_,maxtime_,maxiteration_,verbose_] :=
	Module[{it, itotal, nwl, pill, cfint},
		nwl := NestWhileList[
				(pill:=#;
				# + {RandomReal[ExponentialDistribution[itotal]],First[RandomChoice[it->gammatr,1]]})&,
				{0,X0}, (First[#]<=maxtime && (itotal = Total[it = Intensity[rates,Last[#],alphatr]])>0)&, 1, maxiteration
			   ];
		If[verbose,
			Monitor[nwl,
				MonitorMSG[First[pill]/maxtime,pill,species]
			],
			nwl
		]
	];

FirstReactionStepMethod[species_,alphatr_,gammatr_,rates_,X0_,maxtime_,maxiteration_,verbose_] := 
	Module[{itd, it, min, pill, nwl},
		nwl := NestWhileList[
				(pill := #;
				itd = MyExponentialDistribution[it]; (*Parallelize*)
				min = Min[itd];
				# + {min, min/.Thread[itd->gammatr]})&,
				{0,X0}, (First[#]<=maxtime && Total[it = Intensity[rates,Last[#],alphatr]]>0)&, 1, maxiteration
			   ];
		If[verbose,
			Monitor[nwl,
				MonitorMSG[First[pill]/maxtime,pill,species]
			],
			nwl
		]
	];


(* ::Subsubsection::Closed:: *)
(*Next reaction method*)


NextReactionStepMethod[species_,alphatr_,gammatr_,rates_,X0_,maxtime_,maxiteration_,verbose_] := 
	Module[{r, rr, ss, ittime, nexttime, nextstep, rmin, itold, itnew, xi, txi, transt, pill, nwl},
		r = Length[rates];
		rr = Range[r];
		itnew = N[Intensity[rates,X0,alphatr]]; (**)
		ss = MyExponentialDistribution[itnew];

		ittime = Transpose[{ss, ss, itnew, ConstantArray[0,r]}];
		ittime = If[First[ittime[[#]]]===Infinity, {Infinity, Infinity, 0.`, -Infinity}, ittime[[#]]] &/@ rr;

		nwl := NestWhileList[
				(pill := #;
				(*Print[ittime];*)
				itold = itnew;
				ss = First /@ ittime;
				nexttime = Min[ss];
				rmin = nexttime /. Thread[ss->rr];
				nextstep = Last[#] + gammatr[[rmin]];
				itnew = N[Intensity[rates,nextstep,alphatr]]; (**)

				xi = MyExponentialDistribution[itnew[[rmin]]];

				ittime = If[itnew[[#]]==0.`, Join[{Infinity}, Rest[ittime[[#]]]],
							If[#===rmin, {nexttime + xi, nexttime + xi, itnew[[#]], nexttime}, 
								If[Last[ittime[[#]]]===-Infinity,
									txi = MyExponentialDistribution[itnew[[#]]];
									{nexttime + txi, nexttime + txi, itnew[[#]], nexttime},
									transt = (ittime[[#]][[2]] / itnew[[#]]) * (ittime[[#]][[2]] - Last[ittime[[#]]]);
									{nexttime + transt, nexttime + transt, itnew[[#]], nexttime}
								]
							]
						] &/@ rr;

				{nexttime, nextstep})&,
				{0,X0}, (First[#]<=maxtime && Total[itnew]>0)&, 1, maxiteration
			   ];
		If[verbose,
			Monitor[nwl,
				MonitorMSG[First[pill]/maxtime,pill,species]
			],
			nwl
		]
	];


(* ::Subsubsection::Closed:: *)
(*Tau-leaping methods*)


MyPoissonDistribution[0.`] := 0;
MyPoissonDistribution[0] := 0;
MyPoissonDistribution[x_] := RandomInteger[PoissonDistribution[x]];

SetAttributes[MyPoissonDistribution,Listable];


HOR[{1, 1}, x_] := 1;
HOR[{2, 1}, x_] := 2;
HOR[{2, 2}, x_] := 2 + 1/(x - 1);
HOR[{3, 1}, x_] := 3;
HOR[{3, 2}, x_] := 3 + 3/(2*x-2);
HOR[{3, 3}, x_] := 3 + 1/(x - 1) + 2/(x - 2);

HOR[alphatr_?MatrixQ] := 
	Module[{xm, maxs, rsteporders, m},
		m = Length[First[alphatr]]; (**)
		rsteporders = Total /@ alphatr; (**)
		xm = Array[Subscript[Global`x, #]&, m];
		maxs = First[MaximalBy[Transpose[{rsteporders, #}], Times, 1]] &/@ Transpose[alphatr]; (*Requires at least Mma 10*)
		MapThread[HOR[#1,#2]&, {maxs, xm}, 1]
	];

HOR[___] := 0;


TauLeap[gamma_,intensity_,tolX_,hor_] := 
	Module[{qt, rr},
		qt = Max[#,1] &/@ (Quiet[tolX / hor]/.{ComplexInfinity->Infinity,Indeterminate->Infinity}); (**)
		Min[Join[Quiet[qt / Abs[(# &/@ gamma) . intensity]], Quiet[qt^2 / ((#^2 &/@ gamma) . intensity)]]/.{ComplexInfinity->Infinity,Indeterminate->Infinity}] (**)
	];


ExplicitTauLeapingMethod[species_,alphatr_,gammatr_,rates_,X0_,maxtime_,maxiteration_,stepsize_,tolerance_,verbose_] :=
	Module[{ nwl, pill, it, tau, xm, hor },
		xm = Array[Subscript[Global`x, #]&, Length[species]]; (**)
		hor = HOR[alphatr];

		nwl := NestWhileList[
				(pill := #;
				tau = N[Min[stepsize, TauLeap[Transpose[gammatr], it, tolerance * Last[#], N[Quiet[hor /. Thread[xm -> Last[#]]]]]]]; (**)
				# + {tau, Total[MyPoissonDistribution[it * tau] * gammatr]})&,
				{0,X0}, (First[#]<=maxtime && Min[Last[#]]>=0 && Total[it = Intensity[rates,Last[#],alphatr]]>0)&, 1, maxiteration
			   ];
		If[verbose,
			Monitor[nwl,
				MonitorMSG[First[pill]/maxtime,pill,species]
			],
			nwl
		]
	];


ImplicitTauLeapingMethod[species_,alphatr_,gammatr_,rates_,vars_,X0_,maxtime_,maxiteration_,stepsize_,tolerance_,verbose_] :=
	Module[{nwl, pill, it, tau, xm, hor, poi, rhs, sol, nonneg, intvars},
		nonneg = Thread[vars>=0];
		intvars = Intensity[rates,vars,alphatr];
		xm = Array[Subscript[Global`x, #]&, Length[species]]; (**)
		hor = HOR[alphatr];

		nwl := NestWhileList[
				(pill := #;
				tau = N[Min[stepsize, TauLeap[Transpose[gammatr], it, tolerance * Last[#], N[Quiet[hor /. Thread[xm -> Last[#]]]]]]];
				poi = MyPoissonDistribution[it*tau];
				rhs = poi-it*tau+intvars*tau;
				sol = Flatten[FindInstance[Join[nonneg,Thread[vars-Last[#]==Total[rhs*gammatr]]],vars,Reals]];(*what if there is no solution*)
				# + {tau, Total[Round[rhs/.sol]*gammatr]})&,
				{0,X0}, (First[#]<=maxtime && Min[Last[#]]>=0 && Total[it = Intensity[rates,Last[#],alphatr]]>0)&, 1, maxiteration
			   ];
		If[verbose,
			Monitor[nwl,
				MonitorMSG[First[pill]/maxtime,pill,species]
			],
			nwl
		]
	];


TrapezoidalTauLeapingMethod[species_,alphatr_,gammatr_,rates_,vars_,X0_,maxtime_,maxiteration_,stepsize_,tolerance_,verbose_] :=
	Module[{nwl, pill, it, tau, xm, hor, poi, rhs, sol, nonneg, intvars},
		nonneg = Thread[vars>=0];
		intvars = Intensity[rates,vars,alphatr];
		xm = Array[Subscript[Global`x, #]&, Length[species]]; (**)
		hor = HOR[alphatr];

		nwl := NestWhileList[
				(pill := #;
				tau = N[Min[stepsize, TauLeap[Transpose[gammatr], it, tolerance * Last[#], N[Quiet[hor /. Thread[xm -> Last[#]]]]]]];
				poi = MyPoissonDistribution[it*tau];
				rhs = poi-it*tau/2+intvars*tau/2;
				(*sol = Flatten[FindInstance[Thread[vars-Last[#]==Total[rhs*gammatr]],vars,Reals]];*)
				sol = Flatten[FindInstance[Join[nonneg,Thread[vars-Last[#]==Total[rhs*gammatr]]],vars,Reals]];
				(*gamma.rhs, what if there is no solution*)
				# + {tau, Total[Round[rhs/.sol]*gammatr]})&,
				{0,X0}, (First[#]<=maxtime && Min[Last[#]]>=0 && Total[it = Intensity[rates,Last[#],alphatr]]>0)&, 1, maxiteration
			   ];
		If[verbose,
			Monitor[nwl,
				MonitorMSG[First[pill]/maxtime,pill,species]
			],
			nwl
		]
	];


(* ::Subsubsection::Closed:: *)
(*Simulation*)


Simulation::usage = "Simulation[{reactions},rratecoeffs,init,maxtime,options] simulates the \
induced kinetic Markov process of the reaction endowed with stochastic mass action type kinetics with stochastic reaction rate coefficients given by rratecoeffs. \
The initial state are given by init, while the stopping time is specified by maxtime. \n\
Options include ExternalSpecies, FixedStepSize, Tolerance (for tau-leaping), MaxIteration, Method, Verbose and Volume.
The following methods can be used: \"Direct\" (default), \"FirstReaction\", \
\"NextReaction\", \"ExplicitTau-Leaping\", \"ImplicitTau-Leaping\" and \"TrapezoidalTau-Leaping\".";


Options[Simulation] = {ExternalSpecies->{}, FixedStepSize -> 1000, Tolerance -> 0.01, MaxIteration -> Infinity, Method -> "Direct", Verbose -> True, Volume -> Infinity}; 

Simulation::badarg = "Illegal argument of function Simulation.";
Simulation::method = "The given method '`1`' is not available.";
Simulation::time = "Argument `1` is not positive.";
Simulation::rsarg = "Argument '`1`' or '`2`' may have wrong shape or have inappropriate entries.";
Simulation::intarg = "The initial number of species must be nonnegative integers.";
Simulation::volarg = "Argument `1` is not positive.";

SyntaxInformation[Simulation]={"ArgumentsPattern"->{{__},{__},{__},_,OptionsPattern[]}};

Simulation[{reactions__}, rates_?VectorQ, X0_?VectorQ, maxtime_?NumericQ, opts : OptionsPattern[]] := 
	Module[{
			exs, fss, tol, maxit, method, verbose, v, cou, coi, rd, res, alpha, beta, gamma, 
			m, r, sp, vars, alphatr, betatr, gammatr, x, error
			},
		{exs, fss, tol, maxit, method, verbose, v} = {ExternalSpecies, FixedStepSize, Tolerance, MaxIteration, Method, Verbose, Volume} /. Flatten[{opts, Options[Simulation]}];
		If[v===Infinity && (Or @@ (Not[IntegerQ[#]]& /@ X0)), Message[Simulation::intarg]; Return[$Failed];];
		If[Positive[v],
			cou[vol_,atr_] := If[vol===Infinity, rates, ChangeUnits[vol,rates,atr]];
			coi[vol_,atr_] := If[vol===Infinity, atr, Round[atr*AvogadrosNumber*vol]];,
			Message[Simulation::volarg,v];
			Return[$Failed];
		];
		rd := Check[ReactionsData[{reactions},Flatten[{exs}]]["\[Alpha]","\[Beta]","\[Gamma]","M","R","species","variables"], Return[$Failed];,
						{ReactionsData::wrreac,ReactionsData::badarg}]; 
		If[verbose,
			x = Dynamic[Refresh[Round[Clock[Infinity]],UpdateInterval->1]];
			Monitor[res = rd,
						Column[{Row[{"ReactionsData is now calculating ",Dynamic[Mod[First[x],5]/.{0->"",1->".",2->"..",3->"...",4->"...."}]}],
								ProgressIndicator[Dynamic[Clock[Infinity]],Indeterminate,ImageMargins->1,BaselinePosition->Center],
								Row[{x,". ","seconds passed"}],
								Row[{"Simulation process will be started soon."}]
						}]	
			];
			FinishDynamic[];,
			res = rd;
		];
		{alpha, beta, gamma, m, r, sp, vars} = res;
		{alphatr, betatr, gammatr} = Normal /@ (Transpose /@ {alpha, beta, gamma});

		If[r===Length[rates] && m===Length[X0] && Min[X0]>=0 && Min[rates]>=0,
			If[maxtime>0,
				Switch[method,
					"Direct",
						DirectMethod[sp,alphatr,gammatr,cou[v,alphatr],coi[v,X0],maxtime,maxit,verbose],
					"FirstReaction",
						FirstReactionStepMethod[sp,alphatr,gammatr,cou[v,alphatr],coi[v,X0],maxtime,maxit,verbose],
					"NextReaction", 
						NextReactionStepMethod[sp,alphatr,gammatr,cou[v,alphatr],coi[v,X0],maxtime,maxit,verbose],
					"ExplicitTau-Leaping",
						ExplicitTauLeapingMethod[sp,alphatr,gammatr,cou[v,alphatr],coi[v,X0],maxtime,maxit,fss,tol,verbose],
					"ImplicitTau-Leaping",
						ImplicitTauLeapingMethod[sp,alphatr,gammatr,cou[v,alphatr],vars,coi[v,X0],maxtime,maxit,fss,tol,verbose],
					"TrapezoidalTau-Leaping",
						TrapezoidalTauLeapingMethod[sp,alphatr,gammatr,cou[v,alphatr],vars,coi[v,X0],maxtime,maxit,fss,tol,verbose],
					_, 
						Message[Simulation::method,method]; 
						Return[$Failed];
				],
				Message[Simulation::time,maxtime]; 
				Return[$Failed];
			],
			Message[Simulation::rsarg,rates,X0]; 
			Return[$Failed];
		]
	];

Simulation[___] := (Message[Simulation::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*SimulationPlot, SimulationPlot2D, SimulationPlot3D*)


SimulationPlot::usage = "SimulationPlot[{reactions},rratecoeffs,init,maxtime,options] simulates and plots the trajectories of species, where \
the induced kinetic Markov endowed with stochastic mass action type kinetics with stochastic reaction rate coefficients given by rratecoeffs.\n\
Options include ExternalSpecies, PlotFunction (default is \"ListLinePlot\") and Species (default is All).";


listplotfunctions = {"ListPointPlot3D", "ListLinePlot", "ListPlot", "ListLogPlot"};


Options[SimulationPlot] = {ExternalSpecies->{}, PlotFunction -> "ListLinePlot", Species -> All};

SimulationPlot::badarg = "Illegal argument of function SimulationPlot.";
SimulationPlot::species = "Wrong format or none of the species '`1`' is in the list of internal species of reaction {`2`,...}.";
SimulationPlot::plotfunc = "The symbol '`1`' is not a list plotting function. Try \"ListPointPlot3D\", \"ListLinePlot\", \"ListPlot\" or \"ListLogPlot\".";

SimulationPlot[{reactions__}, rates_?VectorQ, X0_?VectorQ, maxtime_?NumericQ, opts___?OptionQ] := 
	Module[{ops, plot, sp, givensp, proc, w},
			ops = Flatten[{opts, Options[SimulationPlot]}];

			If[MemberQ[listplotfunctions,(plot = (PlotFunction /. ops))],
				(*Monitor*)
				plot = ToExpression[plot]; (*Flatten[{ExternalSpecies/.ops}]*)
				sp = Check[ReactionsData[{reactions},ExternalSpecies/.ops]["species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
				givensp = (Species /. ops) /. All->sp;

				If[(w = Flatten[Position[sp,#] &/@ Union[Flatten[{givensp}]]]) =!= {},
					proc = Check[Simulation[{reactions},rates,X0,maxtime, FilterRules[{opts},Options[Simulation]]],Return[$Failed];,
								{Simulation::volarg,Simulation::rsarg,Simulation::intarg,Simulation::time,Simulation::method,Simulation::badarg}];
					(*Transpose[Map[Thread[List[#]]&,proc]][[w]], maskent kiszedve az anyagokat*) (* Outer[Flatten[#2][[{1,#1+1}]]&,w,proc,1] *)
					plot[Transpose[Map[Thread[List@@#]&,proc]][[w]], FilterRules[{opts},Options[plot]]/.(Method->x_):>Sequence[]], 
					(*Show, PlotLegend*)
					Message[SimulationPlot::"species",givensp,First[Flatten[{reactions}]]];
					Return[$Failed];
				],
				Message[SimulationPlot::"plotfunc",plot];
				Return[$Failed];
			]
	];

SimulationPlot[___] := (Message[SimulationPlot::"badarg"]; $Failed)


SimulationPlot2D::usage = "SimulationPlot2D[{reactions},rratecoeffs,init,maxtime,options] simulates and plots the trajectories of pairs of species, where \
the induced kinetic Markov process endowed with stochastic mass action type kinetics with stochastic reaction rate coefficients given by rratecoeffs.\n\
Options include ExternalSpecies, PlotFunction (default is \"ListPlot\") and Species (default is All).";


SpecSubsets[l_,w_] := Part[#,w] &/@ l;

listplot2dfunctions = {"ListPlot","ListLinePlot"}; (*ParametricPlot*)


Options[SimulationPlot2D] = {ExternalSpecies->{}, PlotFunction->"ListPlot", Species -> All};

SimulationPlot2D::badarg = "Illegal argument of function SimulationPlot2D.";
SimulationPlot2D::species = "Wrong format or none of the given list of pairs of species '`1`' is in the list of internal species of reaction {`2`,...}.";
SimulationPlot2D::plotfunc = "The symbol '`1`' is not a valid 2D plotting function. Try \"ListPlot\" or \"ListLinePlot\".";

SimulationPlot2D[{reactions__}, rates_?VectorQ, X0_?VectorQ, maxtime_?NumericQ, opts___?OptionQ] := 
	Module[{ ops, plot, plotopts, sp, givensp, proc, w },
			ops = Flatten[{opts, Options[SimulationPlot2D]}];
			If[MemberQ[listplot2dfunctions,(plot = PlotFunction /. ops)],
				(*Monitor*)
				plot = ToExpression[plot];
				sp = Check[ReactionsData[{reactions},ExternalSpecies/.ops]["species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
				givensp = (Species /. ops) /. All->Subsets[sp,{2}];
				If[MatrixQ[givensp] && Last[Dimensions[givensp]]===2 &&(w = Flatten[Position[sp,#]&/@#] &/@ givensp) =!= {},
					proc = Check[Last /@ Simulation[{reactions},rates,X0,maxtime,FilterRules[{opts},Options[Simulation]]],Return[$Failed];,
								{Simulation::volarg,Simulation::rsarg,Simulation::intarg,Simulation::time,Simulation::method,Simulation::badarg}];
					(*Print[w];*)
					Show[
						plot[SpecSubsets[proc,#]& /@ w,FilterRules[{opts},Options[plot]]/.(Method->x_):>Sequence[]]
					],
					(*Show, PlotLegend*)
					Message[SimulationPlot2D::species,givensp,First[Flatten[{reactions}]]];
					Return[$Failed];
				],
				Message[SimulationPlot2D::plotfunc,plot];
				Return[$Failed];
			]
	];

SimulationPlot2D[___] := (Message[SimulationPlot2D::badarg]; $Failed)


SimulationPlot3D::usage = "SimulationPlot3D[{reactions},rratecoeffs,init,maxtime,options] simulates and plots the trajectories of triples of species, where \
the induced kinetic Markov process endowed with stochastic mass action type kinetics with stochastic reaction rate coefficients given by rratecoeffs.\n\
Options include ExternalSpecies, PlotFunction (default is \"ListPointPlot3D\") and Species (default is All).";


listplot3dfunctions = {"ListPlot3D","ListPointPlot3D"}; (*ParametricPlot*)


Options[SimulationPlot3D] = {ExternalSpecies->{}, PlotFunction->"ListPointPlot3D", Species -> All};

SimulationPlot3D::badarg = "Illegal argument of function SimulationPlot3D.";
SimulationPlot3D::species = "Wrong format or none of the given list of species '`1`' is in the list of internal species of reaction {`2`,...}.";
SimulationPlot3D::plotfunc = "The given symbol '`1`' is not a valid 3D plotting function. Try \"ListPlot3D\" or \"ListPointPlot3D\".";

SimulationPlot3D[{reactions__}, rates_?VectorQ, X0_?VectorQ, maxtime_?NumericQ, opts___?OptionQ]:=
	Module[{ ops, plot, plotopts, sp, givensp, proc, w },
			ops = Flatten[{opts, Options[SimulationPlot3D]}];
			If[MemberQ[listplot3dfunctions,(plot = PlotFunction /. ops)],
				plot = ToExpression[plot];
				(*Monitor*)
				sp = Check[ReactionsData[{reactions},ExternalSpecies/.ops]["species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
				givensp = (Species /. ops) /. All->Subsets[sp,{3}];
				If[MatrixQ[givensp] && Last[Dimensions[givensp]]===3 &&(w = Flatten[Position[sp,#]&/@#] &/@ givensp) =!= {},
					proc = Check[Last /@ Simulation[{reactions},rates,X0,maxtime,FilterRules[{opts},Options[Simulation]]],Return[$Failed];,
								{Simulation::volarg,Simulation::rsarg,Simulation::intarg,Simulation::time,Simulation::method,Simulation::badarg}];
					Show[
						plot[SpecSubsets[proc,#]& /@ w,FilterRules[{opts},Options[plot]]/.(Method->x_):>Sequence[]]
					],
					Message[SimulationPlot3D::species,givensp,First[Flatten[{reactions}]]];
					Return[$Failed];
				],
				Message[SimulationPlot3D::plotfunc,plot];
				Return[$Failed];
			]
	];

SimulationPlot3D[___] := (Message[SimulationPlot3D::badarg]; $Failed)


(* ::Subsection::Closed:: *)
(*decompositions*)


solvedLPs = 0;
currentbest = Infinity;


VectorFromPoints::badarg = "Illegal argument of function VectorFromPoints.";

VectorFromPoints[len_Integer, pts_List, vals_List, pad_:0] :=
	Module[ {u = Table[pad, {len}]}, 
			(u[[ First[#] ]] = Last[#]) & /@ Transpose[{pts, vals}];
			 u
	];

VectorFromPoints[___] := (Message[VectorFromPoints::"badarg"]; Return[$Failed];)

FromStep::badarg = "Illegal argument of function FromStep.";

FromStep[step_, species_] :=
	Subtract @@ (VectorFromPoints[Length[species],
					Sequence @@ Transpose[#[Flatten[List[#]] & /@ (step /. Plus -> List) /.
				  	coeff_. spec_?(MemberQ[species, #]&) :> {Position[species, spec][[1,1]], coeff}]]] & /@ {Last, First}
	);

FromStep[___] := (Message[FromStep::badarg]; Return[$Failed];)

ToStep::badarg = "Illegal argument of function ToStep.";

ToStep[vector_?ReactionQ, species_List] :=
	First[
    		Inner[	Rule,
			{Abs[#2[#1, #3@Positive]] . #4},
			{#2[#1, #3@Negative] . #4},
			List
		] & [vector, ReplaceAll, _?# -> 0 &, species]
	];

ToStep[___] := (Message[ToStep::"badarg"]; Return[$Failed];)


ReactionQ = MatchQ[#, {__Integer}] &;

DecompositionQ = MatchQ[#, {__Integer?NonNegative}] &;

ElementaryReactionQ = ReactionQ[#] && Total[Select[#, Negative]] >= -2 &;

GeneralizedDecompositionQ = MatchQ[#, { (_Integer?NonNegative | _Rational?NonNegative).. }] &;


Circular[overall_?ReactionQ, steps : {__?ReactionQ}] :=
	Block[{$Messages},
    	Check[
			LinearProgramming[
				Array[-1 &, Length[steps]],
				Transpose[steps],
				{#, 0} & /@ overall
			],
      		Return[True],
      		LinearProgramming::lpsub
      	]; False
    ];


CoveringDecompositionSet::usage = "CoveringDecompositionSet[overallreaction,reactionsteps,options] gives some of the decompositions \
of the overall reaction into the reaction steps denoted by reactionsteps. Note that CoveringDecompositionSet deletes autocatalytic steps of the reaction.\n\
CoveringDecompositionSet[overall,gamma,options] gives the same where overall is the transpose of the stoichiometrix matrix of \
the overall reaction, while gamma is the stoichiometric matrix of the reaction steps.";

CoveringDecompositionSet::badarg = "Illegal argument of function CoveringDecompositionSet.";
CoveringDecompositionSet::args = "One of the arguments has wrong shape.";
CoveringDecompositionSet::badmeth = "The given method cannot be recognized. Try GreedySelection, FastSelection or OriginalSelection.";
CoveringDecompositionSet::nulldecomp = "Overall reaction contains a species which is not among those of the given reaction steps.";
CoveringDecompositionSet::isnotoverall = "The overal reaction has to be a single reaction step.";

Options[CoveringDecompositionSet] = {ObjectiveFunction -> GreedySelection, Verbose -> True}

CoveringDecompositionSet[{overallreaction_?OverallQ},{reactions__},opts___?OptionQ] := 
	Module[{elemgamma, elemspecs, overallgamma, overallspecs, orderelem, ordercompl, steps, overall},

		{overallgamma, overallspecs} = Check[ReactionsData[Flatten[{overallreaction}],MyFilterOptions[ReactionsData,opts]]["\[Gamma]","species"],Return[$Failed],
												{ReactionsData::wrreac,ReactionsData::badarg}];
		If[Length[First[overallgamma]]=!=1,Message[CoveringDecompositionSet::isnotoverall];Return[$Failed];];

		{elemgamma, elemspecs} = Check[ReactionsData[Flatten[{reactions}],MyFilterOptions[ReactionsData,opts]]["\[Gamma]","species"],Return[$Failed],
											{ReactionsData::wrreac,ReactionsData::badarg}];

		If[Complement[overallspecs,elemspecs]=!={},Message[CoveringDecompositionSet::nulldecomp];Return[{}]];

		orderelem = Flatten[Position[elemspecs,#] &/@ overallspecs];
		ordercompl = Flatten[Position[elemspecs,#] &/@ Complement[elemspecs,overallspecs]];

		steps = Normal[Join[elemgamma[[orderelem]],elemgamma[[ordercompl]]]]; (*SparseArray*)
		overall = Join[Flatten[Normal[overallgamma]],ConstantArray[0,Length[ordercompl]]]; (*SparseArray*)

		CoveringDecompositionSet[overall,steps,opts]

	];

CoveringDecompositionSet[overall_?ReactionQ, steps_?MatrixQ, opts___?OptionQ] :=
	Module[{$Messages,
			verbose = Verbose /. Flatten[{opts, Options[CoveringDecompositionSet]}],
			objective = ObjectiveFunction /. Flatten[{opts, Options[CoveringDecompositionSet]}],
			reactions,
			markedspecies = Array[1 &, Length[First[steps]]],
			decompositions = {},
			nextdecomposition,
			coefficientvector,
			mainf,
			monitormsg = "Working..."
	      },

	If[Not[VectorQ[overall]]||Length[steps]=!=Length[overall],Message[CoveringDecompositionSet::args];Return[$Failed];];
	If[Not[MemberQ[{ReactionKinetics`GreedySelection,ReactionKinetics`FastSelection,ReactionKinetics`OriginalSelection},objective]],
				Message[CoveringDecompositionSet::badmeth];Return[$Failed];];

	coefficientvector :=
		Switch[objective,
				ReactionKinetics`GreedySelection, 1 - markedspecies,
				ReactionKinetics`FastSelection, Array[0 &, Length[First[steps]]],
				ReactionKinetics`OriginalSelection, markedspecies
			   ];

	mainf := While[True,
		If[Max[markedspecies] === 0, Break[]]; (**)
		nextdecomposition =
			Check[
				LinearProgramming[
					coefficientvector,
					Prepend[steps, markedspecies],
					Prepend[{#,0} & /@ overall, {1,1}]
				],
          			Break[],
          			LinearProgramming::"lpsnf"
          		];
		PrependTo[decompositions, nextdecomposition];
		markedspecies *= (1 - Sign[nextdecomposition]);
		monitormsg = "Found " <> ToString[Length[decompositions]] <>
					  " decomposition(s) covering " <> ToString[Total[1 - markedspecies]] <>
					  " of the " <> ToString[Length[First[steps]]] <> " reaction steps."
	];

	If[verbose, Monitor[mainf, monitormsg], mainf];

	If[verbose, Print["Found " <> ToString[Length[decompositions]] <> " decomposition(s) covering " <>
		ToString[Total[1 - markedspecies]] <> " of the " <> ToString[Length[First[steps]]] <> " reaction steps."
	]];

	decompositions

];

ElementaryReactions::usage = "ElementaryReactions[species,maxproduct,options] computes every possible \
elementary steps among the species having at most maxproduct products. \
The second argument is optional, its default value is Infinity.\n\
ElementaryReactions[atommatrix,maxproduct,options] does the same from the atomic matrix given by atommatrix. \
This latter returns the transpose of the stoichiometric matrix \[Gamma] of the possible elementary reaction steps.";

ElementaryReactions::badarg = "Illegal argument of function ElementaryReactions.";

Options[ElementaryReactions] = {Verbose -> True}

ElementaryReactions[species_?VectorQ, maxlen : (_Integer?Positive | Infinity) : Infinity, opts___?OptionQ] :=
	Module[{headers, atomm, elemgamma},
		{headers,atomm} = Check[ToAtomMatrix[species,FormattedOutput->False],Return[$Failed];,{ToAtomMatrix::badarg}];
		elemgamma = Check[Transpose[ElementaryReactions[atomm,maxlen,opts]],Return[$Failed];,{ElementaryReactions::badarg}];

		FromStoichiometry[NegativePart[elemgamma],PositivePart[elemgamma],species]
];

ElementaryReactions[atomm_?MatrixQ, maxlen : (_Integer?Positive | Infinity) : Infinity, opts___?OptionQ] :=
	Module[{$Messages, verbose = Verbose /. Flatten[{opts, Options[ElementaryReactions]}],
		molec, others, reac, i, j, atommatrix, fda, monitorf, monitormsg = "Working...", casesdone},

		atommatrix = Transpose[atomm];
		fda = First[Dimensions[atommatrix]];

		monitorf := Join[
			Flatten[
				Table[
					molec = atommatrix[[i]];
					casesdone = 2*i;
	    			others = Delete[atommatrix, {{i}}];
	    			{
						Insert[#, -1, {{i}}] & /@ Decompositions[  molec, Transpose[others], maxlen, Flatten[{Verbose -> False, MyFilterOptions[Decompositions, opts]}]],
						Insert[#, -2, {{i}}] & /@ Decompositions[2 molec, Transpose[others], maxlen, Flatten[{Verbose -> False, MyFilterOptions[Decompositions, opts]}]]
    				},
					{i, fda}
				],
				2
			],
			Flatten[
				Table[
    				molec = atommatrix[[i]] + atommatrix[[j]];
					casesdone = -1 + fda - i/2 + fda*i - i^2/2 + j; 
	    			others = Delete[atommatrix, {{i}, {j}}];
					Insert[#, -1, {{i}, {j-1}}] & /@ Decompositions[molec, Transpose[others], maxlen, Flatten[{Verbose -> False, MyFilterOptions[Decompositions, opts]}]]
    				,
	    			{i, fda-1}, {j, i+1, fda}
    			],
    			2
    		]
    	];
		If[verbose, Monitor[monitorf, ToString[IntegerPart[100.*casesdone/((fda*(3 + fda))/2)]]<>"% done"], monitorf]
    ];

ElementaryReactions[___] := (Message[ElementaryReactions::badarg]; $Failed)


Obligatory::usage = "Obligatory[overall,gamma] selects those reaction steps from gamma \
which are provably contained in each decomposition of the overall reaction into reaction steps. \
The function Obligatory returns a list of pairs that contain the mandatory reaction steps \
and their weights in a decomposition.";

Obligatory::badarg = "Illegal argument of function Obligatory.";

Obligatory[overall_?ReactionQ, reactions : {__?ReactionQ}] :=
	Transpose[{
		Flatten[Position[#, _?Positive]], Select[#, Positive]
			  }] & [Array[
				LinearProgramming[
					Table[If[i===#, 1, 0], {i,Length[First[reactions]]}],
					reactions,
					{#, 0} & /@ overall
				][[#]] &,
				Length[First[reactions]]
			]];

Obligatory[___] := (Message[Obligatory::badarg]; $Failed)

Omittable::usage = "Omittable[overall,gamma,options] selects those reaction steps from gamma \
which can be omitted from the decompisitions of the overall reaction into reaction steps.";

Omittable::badarg = "Illegal argument of function Omittable.";

Omittable[overall_?ReactionQ, reactions_?MatrixQ, opts___?OptionQ] :=
	Complement[
		Range[Length[First[reactions]]],
		Union[Last /@ Position[CoveringDecompositionSet[overall, reactions, MyFilterOptions[CoveringDecompositionSet,opts]], _?Positive]]
	];

Omittable[___] := (Message[Omittable::badarg]; $Failed)

SelectMinimalDecompositions::usage = "SelectMinimalDecompositions[list] selects those vectors \
from the given list which are minimal with respect to the \"component-wise less or equal\" partial order. \
Note that using this function to select minimal decompositions will not necessarily give truly \
minimal decompositions, unless all decompositions up to a certain length are given in the list.";

SelectMinimalDecompositions::badarg = "Illegal argument of function SelectMinimalDecompositions.";

SelectMinimalDecompositions[vlist : {__?DecompositionQ}] :=
	SMaux[First[#], Rest[#]] &[
		Map[Last,
			Split[Sort[{Total[#], #} & /@ vlist, First[#1] < First[#2] &], First[#1] === First[#2] &],
			{2}]
		];

SelectMinimalDecompositions[___] := (Message[SelectMinimalDecompositions::badarg]; $Failed)


Options[CDPartitions] = {Verbose -> True, Filter -> MinimalDecompositions}

CDPartitions[overall_?ReactionQ, steps : {__?ReactionQ}, stepno_, opts___?OptionQ] :=
	Module[
		{verbose, filter,
		 vecs = Append[steps, -overall],
		 im = IdentityMatrix[Length[steps] + 1],
		 solutions = {},
		 stack, product, tmpsol, stepstogo,
		 chooseable, properoriented,
		 x, verboseStepnum = 0, verboseMaxdepth = 0,
		 monitormsg = "Working...", monitorf},

		{verbose, filter} = {Verbose, Filter} /. Flatten[{opts, Options[Decompositions]}];
		If[filter === All, Message[Decompositions::cdwall]; Return[$Failed]];

		stack = Array[
					{steps[[#]] - overall,
					 im[[#]] + im[[Length[vecs]]],
					 stepno,
					 Range[#]} &,
					Length[steps]
				];

		monitorf := While[Length[stack] =!= 0,

			{product, tmpsol, stepstogo, chooseable} = First[stack];
			++verboseStepnum;
			monitormsg = "Running... " <> ToString[Length[solutions]] <>
								  " decomposition(s) found in " <> ToString[verboseStepnum] <>
								  " reaction steps.";
			verboseMaxdepth = Max[verboseMaxdepth, Total[tmpsol]];

			Which[
				(*Total[tmpsol] > stepno ||*) stepstogo === 0 || ! MinimalQ[tmpsol, solutions], stack = Rest[stack],
				ZeroVectorQ[product], PrependTo[solutions, tmpsol];
									  stack = Rest[stack],
				True, properoriented = Select[chooseable, product . vecs[[#]] <= -(product . product)/stepstogo && (* this is the improvement *)
														  product . vecs[[#]] < 0 &&
														  MinimalQ[tmpsol + im[[#]], solutions] &];
					  stack = Join[
								Map[Function[x, {product + vecs[[x]],
            									 tmpsol + im[[x]],
												 stepstogo - 1,
												 Complement[chooseable, Select[properoriented, # > x &]]
												}],
        							properoriented
								],
								Rest[stack]];
			];
		];
		If[verbose, Monitor[monitorf, monitormsg], monitorf];

		If[verbose,
			Print[Length[solutions], " decomposition(s) found in ", verboseStepnum,
			 " steps.\nMaximum search depth was ", verboseMaxdepth, ", the longest decomposition consists of ",
        	 Max[Total /@ solutions] - 1, " reaction step(s)."
			]
		];
		Most /@ solutions
	];

(* ----------------------------------------------------------------------------------------------------- *)

UnrestrictedSolutions[m_?MatrixQ, b_List] :=
	Module[{ls, ns},
		Check[
			ls = LinearSolve[m, b],
			Return[{}],
			LinearSolve::nosol
		];
    	If [(ns = NullSpace[m]) === {}, ls, ls + Transpose[ns] . (x /@ Range[Length[ns]])]
    ];

Options[LPPartitions] = {Filter -> All, Preprocess -> False, Verbose -> True}

LPPartitions[overall_?ReactionQ, steps : {__?ReactionQ}, maxlen_, opts___?OptionQ] :=
	Block[{verbose, filter, preprocess, generalsol, vars, solvedLPs = 0, constraints},

		{verbose, filter, preprocess} = {Verbose, Filter, Preprocess} /. Flatten[{opts, Options[LPPartitions]}];

		If[maxlen === Infinity,
			If[Circular[overall, steps] && filter === All, Message[Decompositions::unbnd]; Return[$Failed]];
			If[Circular[overall, steps], Message[Decompositions::wrmet]; Return[$Failed]]
		];

		If[ (generalsol = UnrestrictedSolutions[Transpose[steps], overall]) === {}, Return[{}] ];

		vars = Union[Cases[generalsol, x[_], Infinity]];
		If[vars === {}, Return[
							If[And @@ (IntegerQ[#] && NonNegative[#] & /@ generalsol) && Total[generalsol] <= maxlen,
          					   {generalsol},
							   {}
							]]
    	];

    	constraints = If[maxlen =!= Infinity, Append[generalsol, maxlen - Total[generalsol]], generalsol];

    	(* Check whether the polyhedron of the relaxed problem is empty *)
    	(* Saves a lot of breath and time if it is, good for nothing if not *)
    	If[MyMinimize[x[1], constraints, vars] === "no solution",
    		If[verbose,	Print["0 solution(s) found with 1 solved LPs."]];
			Return[{}]
		];
(**********constraintsols, nem constraints? *)
		constraintsols = (* FindAllSolutionsWithPreprocessing[constraints, vars, verbose]; *)
			If[preprocess,
				FindAllSolutionsWithPreprocessing[constraints, vars, verbose],
				FindAllSolutions[constraints, vars]
			];

		If[verbose,
			Print[ToString[Length[constraintsols]] <> " solution(s) found with " <>
				  ToString[solvedLPs + If[preprocess, 2 Length[vars], 0]] <> " solved LPs."]
		];
    	generalsol /. Thread[vars -> #] & /@ constraintsols

    ];

(* Az uj valtozathoz *)
FindAllSolutionsWithPreprocessing[conlist_List, vars_List, verbose_] :=
	Module[{bounds, bnd, bnd2, filtered = {}, f1 = 0, f2 = 0, newf2 = True, unfiltered, newconlist = conlist,
	filteredsolutions, i, lvar = Length[vars], steprec, monitorf, monitormsg = "Working..."},

		(* If[verbose, Print["Preprocessing..."]]; *)

	monitorf = While[newf2 && Length[filtered] < lvar,
		newf2 = False;
		For[i = 1, i <= lvar, i++,
			monitormsg = ToString[i-1] <> "/" <> ToString[lvar] <> " variables examined, " <>
											ToString[f1] <> " + " <> ToString[f2] <> " variables filtered by Rule #1 and #2.";
			If[! And @@ (IntegerQ /@ Select[newconlist, NumberQ]), Return[{}]];
			If[MemberQ[filtered, i], Continue[]];

			(* steprec = GCD @@ (Cases[newconlist, a_. + b_. x[i] /; NumberQ[a] :> b, {1}]);
			If[steprec == 0, Print[{i, conlist, newconlist}]]; *)

			solvedLPs += 2;
			Check[
				bnd = #[x[i], newconlist, vars] & /@ {MyMinimize, MyMaximize},
				Return[{}],
				LinearProgramming::lpsnf
				(* ?? Miert nem lep ki soha? Akkor a kovetkezo mar nem is kene *)
			];
			If[ !MatchQ[bnd, {_?NumberQ, _?NumberQ}], Return[{}] ];

			If[ SameQ @@ bnd,

				(* Print["Variable ", i, " is FILTERED by Rule no. 1."]; *)
				AppendTo[filtered, i];
				(* AppendTo[bounds, bnd]; *)
				bounds[i] = bnd;
				newconlist = newconlist /. x[i] -> First[bnd];
				++f1;
				If[! And @@ (IntegerQ /@ Select[newconlist, NumberQ]), Return[{}]],

				If[ Length[
						bnd2 = Intersection @@
								(Cases[newconlist, a_. + b_. x[i] /; NumberQ[a] :> Range[
																					NestWhile[# + 1/b &, -a/b, If[Sign[b]>0, # < First[bnd], # > Last[bnd]] &],
																					If[Sign[b]>0, Last[bnd], First[bnd]],
																					1/b
																				   ],
								 {1}])
					] === 1,

					(* Print["Variable ", i, " is FILTERED by Rule no. 2."]; *)
					(*Print[newconlist];
					Print[bnd];
					Print[bnd2];*)

				 	(* Ceiling[First[bnd] steprec] + 1 > steprec Last[bnd] *)

					bnd2 = First[bnd2];
					AppendTo[filtered, i];
					(* bnd2 = Ceiling[First[bnd] steprec] / steprec; *)
					(* AppendTo[bounds, {bnd2, bnd2}]; *)
					bounds[i] = {bnd2, bnd2};
					newconlist = newconlist /. x[i] -> bnd2;
					(* Print[newconlist]; *)
					++f2;
					If[! And @@ (IntegerQ /@ Select[newconlist, NumberQ]), Return[{}]];
					newf2 = True,

						(*Print["Variable ", i, " is NOT filtered."]; *)
						(* AppendTo[bounds, bnd]; *)
						bounds[i] = bnd;
				]
			]
		]
	];

	If[verbose, Monitor[monitorf, monitormsg], monitorf];

		unfiltered = Complement[Range[lvar], filtered];

		If[verbose, Print[Length[filtered], " out of ", lvar, " variable(s) filtered out in preprocessing."]
		];

		If[unfiltered === {},
			Return[ {Array[First[bounds[#]] &, {lvar}]} ],

			(* Print[newconlist /. ReactionKinetics`Private`x -> Global`x];
			Print[First /@ bounds /@ unfiltered];
			Print[Last /@ bounds /@ unfiltered]; *)
			filteredsolutions = FindAllSolutions[
									newconlist,
									x /@ unfiltered,
									(* bounds[[unfiltered, 2]] *)
									First /@ bounds /@ unfiltered,
									Last /@ bounds /@ unfiltered
								];
			(* Print["fsols: ", filteredsolutions]; *)
		];

		(* Return *)
		Array[If[MemberQ[unfiltered, #], 0, (* bounds[[#, 1]] *) First[bounds[#]]] &, {lvar}] +
		VectorFromPoints[lvar, unfiltered, #] & /@ filteredsolutions
	];

(* UJ VALTOZAT, ALSO-FELSO KORLATOKKAL, LESZAMLALASSAL *)
FindAllSolutions[conlist_List, vars : {var_, restvars___}, lowerbounds_List, upperbounds_List] :=
	Module[{lower, upper, newcons, a, b, step, offset},

		{step, offset} = First[Sort[
			Cases[conlist, a_. + b_. var /; NumberQ[a] :> {Abs[1/b], -a/b}, {1}],
			First[#1] > First[#2] &
		]];

		(*	step   = 1 / GCD @@ (Cases[conlist, a_. + b_. var /; NumberQ[a] :> b, {1}]);
			offset = First[Cases[conlist, a_. + b_. var /; NumberQ[a] :> -a/b, {1}, 1]]; *)

		(* Print[{var, step, offset}]; *)
		(* Print[conlist /. ReactionKinetics`Private`x -> Global`x]; *)
		(* If[vars =!= Union[Cases[conlist,x[_],Infinity]], Print[vars]]; *)
		(* Print[solvedLPs]; *)
		(* Print[{Length[vars], step, offset}]; *)

		(* solvedLPs += 2;
		Check[
			lower = MyCeiling[MyMinimize[var, conlist, vars] - offset, step] + offset;
			upper = MyFloor[MyMaximize[var, conlist, vars] - offset, step] + offset,
			Return[{}],
			LinearProgramming::lpsnf
		];*)

		lower = MyCeiling[First[lowerbounds] - offset, step] + offset;
		upper = MyFloor[First[upperbounds] - offset, step] + offset;
		While[Or @@ Negative[conlist /. var -> lower] && (lower <= upper), lower += step];
		While[Or @@ Negative[conlist /. var -> upper] && (lower <= upper), upper -= step];

		(* Return *)
		If[lower > upper,
			{},

			If[Length[vars] === 1,

				(* Print[Length[Select[Range[lower, upper, step], And @@ IntegerQ /@ (conlist /. var -> #) &]]]; *)
				Transpose[{Select[Range[lower, upper, step], And @@ IntegerQ /@ (conlist /. var -> #) &]}],

				Join @@ Table[
					newcons = conlist /. var -> i;
					If[ Not[And @@ IntegerQ /@ Select[newcons, NumberQ]],
					    {},
						Prepend[#, i] & /@ FindAllSolutions[
												DeleteCases[newcons, _Integer],
												{restvars},
												Rest[lowerbounds],
												Rest[upperbounds]
										   ]
					], {i, lower, upper, step}
				] (* Join@@Table*)
			] (* If length *)
		] (* If lower *)
	]; (* Module *)

(* REGEBBI, NINCS ELOFElDOLGOZAS *)
FindAllSolutions[conlist_List, vars : {var_, restvars___}] :=
	Module[{lower, upper, newcons, a, b, step, offset},

		{step, offset} = First[Sort[
			Cases[conlist, a_. + b_. var /; NumberQ[a] :> {Abs[1/b], -a/b}, {1}],
			First[#1] > First[#2] &
		]];

		(* Print[{var, step, offset}]; *)
		(* Print[conlist /. ReactionKinetics`Private`x -> Global`x]; *)
		(* If[vars =!= Union[Cases[conlist,x[_],Infinity]], Print[vars]]; *)
		(* Print[solvedLPs]; *)
		(* Print[{Length[vars], step, offset}]; *)

		solvedLPs += 2;
		Check[
			lower = MyCeiling[MyMinimize[var, conlist, vars] - offset, step] + offset;
			upper = MyFloor[MyMaximize[var, conlist, vars] - offset, step] + offset,
			Return[{}],
			LinearProgramming::lpsnf
		];
		(*
		lower = MyCeiling[First[lowerbounds] - offset, step] + offset;
		upper = MyFloor[First[upperbounds] - offset, step] + offset;
		While[Or @@ Negative[conlist /. var -> lower] && (lower <= upper), lower += step];
		While[Or @@ Negative[conlist /. var -> upper] && (lower <= upper), upper -= step];
		*)
		(* Return *)
		If[lower > upper,
			{},
			If[Length[vars] === 1,
				(* Print[Length[Select[Range[lower, upper, step], And @@ IntegerQ /@ (conlist /. var -> #) &]]]; *)
				Transpose[{Select[Range[lower, upper, step], And @@ IntegerQ /@ (conlist /. var -> #) &]}],
				Join @@ Table[
					newcons = conlist /. var -> i;
					If[ Not[And @@ IntegerQ /@ Select[newcons, NumberQ]],
					    {},
						Prepend[#, i] & /@ FindAllSolutions[
												DeleteCases[newcons, _Integer],
												{restvars}
										   ]
					], {i, lower, upper, step}
				] (* Join@@Table*)
			] (* If length *)
		] (* If lower *)
	]; (* Module *)

(* ----------------------------------------------------------------------------------------------------- *)

Decompositions::usage = "Decompositions[overallreaction,reactionsteps,maxlen,options] gives the decompositions of the overall \
reaction into the reaction steps. Note that Decompositions deletes autocatalytic steps of the reaction.\n\
Decompositions[overall,gamma,maxlen,options] gives the same where overall is the transpose of the stoichiometrix matrix of \
the overall reaction, while gamma is the stoichiometric matrix of the reaction steps.";

Decompositions::badarg = "Illegal argument of function Decompositions.";
Decompositions::cdwall = "Method -> ContejeanDevie can yield the minimal \
solutions only. Choose another method or use the option Filter -> MinimalDecompositions.";
Decompositions::unbnd = "The number of the decompositions is unbounded. Give \
upper bound for solution length or use the option Filter -> MinimalDecompositions.";
Decompositions::wrmet = "Incompatible method and filter option (perhaps due to the circularity of the network of reaction steps).";
(*(Ez fel\[EAcute]r egy bad command or file name-mel, egyelore...)*)

Decompositions::badmeth = "The given method cannot be recognized. Try one of the following: LPBased or ContejeanDevie.";
Decompositions::args = "One of the arguments has wrong shape.";
Decompositions::nulldecomp = "Overall reaction contains a species which is not among those of the given reaction steps.";
Decompositions::isnotoverall = "The overal reaction has to be a single reaction step.";

Options[Decompositions] = {Filter->MinimalDecompositions,Method->ContejeanDevie,Preprocess->True,Verbose->True}
(*Filter\[Rule]All*)

OverallQ := MemberQ[{Rule,
		ShortRightArrow,RightArrow,LongRightArrow,
		ShortLeftArrow,LeftArrow,LongLeftArrow,
		DoubleRightArrow,DoubleLongRightArrow,
		DoubleLeftArrow,DoubleLongLeftArrow}, Head[#]]&;

Decompositions[{overallreaction_?OverallQ},{reactions__},
			   maxlen : (_Integer?Positive | Infinity) : Infinity,opts___?OptionQ] := 
	Module[{elemgamma, elemspecs, overallgamma, overallspecs, orderelem, ordercompl, steps, overall},

		{overallgamma, overallspecs} = Check[ReactionsData[Flatten[{overallreaction}],MyFilterOptions[ReactionsData,opts]]["\[Gamma]","species"],Return[$Failed],
												{ReactionsData::wrreac,ReactionsData::badarg}];
		If[Length[First[overallgamma]]=!=1,Message[Decompositions::isnotoverall];Return[$Failed];];

		{elemgamma, elemspecs} = Check[ReactionsData[Flatten[{reactions}],MyFilterOptions[ReactionsData,opts]]["\[Gamma]","species"],Return[$Failed],
											{ReactionsData::wrreac,ReactionsData::badarg}];

		If[Complement[overallspecs,elemspecs]=!={},Message[Decompositions::nulldecomp];Return[{}]];

		orderelem = Flatten[Position[elemspecs,#] &/@ overallspecs];
		ordercompl = Flatten[Position[elemspecs,#] &/@ Complement[elemspecs,overallspecs]];

		steps = Normal[Join[elemgamma[[orderelem]],elemgamma[[ordercompl]]]]; (*SparseArray*)
		overall = Join[Flatten[Normal[overallgamma]],ConstantArray[0,Length[ordercompl]]]; (*SparseArray*)

		Decompositions[overall,steps,maxlen,opts]

	];

Decompositions[overall_?ReactionQ, steps : {__?ReactionQ},
			   maxlen : (_Integer?Positive | Infinity) : Infinity, opts___?OptionQ] :=
	Module[{method, $RecursionLimit = Infinity, stepstr},
		If[Not[VectorQ[overall]]||Length[overall]=!=Length[steps]||Not[MatrixQ[steps]],Message[Decompositions::args]; Return[$Failed];];

		{method} = {Method} /. Flatten[{opts, Options[Decompositions]}];
		stepstr = Transpose[steps];

		Switch[method,
				LPBased, LPPartitions[overall, stepstr, maxlen, MyFilterOptions[LPPartitions,opts]],
				ContejeanDevie, CDPartitions[overall, stepstr, maxlen, MyFilterOptions[CDPartitions,opts]],
				_, Message[Decompositions::badmeth]; Return[$Failed];
		]
    ];

Decompositions[___] := (Message[Decompositions::badarg]; Return[$Failed];)


(* ::Subsection::Closed:: *)
(*palettes, notebooks*)


OpenReactionKineticsPalette::usage = "OpenReactionKinetics[] provides a palette including the basic symbols \
in ReactionKinetics.";

OpenReactionKineticsNamesPalette::usage = "OpenReactionKineticsNamesPalette[] provides a palette including the basic names and notions \
in ReactionKinetics.";


OpenReactionKineticsPalette::badarg = "Illegal argument of function OpenReactionKineticsPalette.";

OpenReactionKineticsNamesPalette::badarg = "Illegal argument of function OpenReactionKineticsNamesPalette.";


SyntaxInformation[OpenReactionKineticsPalette] = { "ArgumentsPattern"->{} };

OpenReactionKineticsPalette[] :=
    NotebookPut[Notebook[{
        Cell[BoxData[GridBox[
		{
		 {ButtonBox["\[Alpha]"], ButtonBox["\[Beta]"], ButtonBox["\[Gamma]"], ButtonBox["\[Delta]"], ButtonBox["\[Theta]"], ButtonBox["\[Lambda]"]},
		 {ButtonBox["\"M\""], ButtonBox["\"E\""], ButtonBox["\"R\""], ButtonBox["(\[Placeholder])"], ButtonBox["{\[Placeholder]}"], ButtonBox["[\[Placeholder]]"]},
		 {ButtonBox["\[ShortRightArrow]"], ButtonBox["\[RightArrow]"], ButtonBox["\[LongRightArrow]"], ButtonBox["\[ShortLeftArrow]"],ButtonBox["\[LeftArrow]"], ButtonBox["\[LongLeftArrow]"]},
		 {ButtonBox["\[LeftRightArrow]"], ButtonBox["\[LongLeftRightArrow]"], ButtonBox["\[Equilibrium]"], ButtonBox["\[ReverseEquilibrium]"], ButtonBox["\[RightArrowLeftArrow]"], ButtonBox["\[LeftArrowRightArrow]"]}, 
		 {ButtonBox["\[DoubleRightArrow]"], ButtonBox["\[DoubleLongRightArrow]"], ButtonBox["\[DoubleLeftArrow]"], ButtonBox["\[DoubleLongLeftArrow]"], ButtonBox["\[DoubleLeftRightArrow]"], ButtonBox["\[DoubleLongLeftRightArrow]"]},
			(*ButtonBox["\[Placeholder]\!\(\*OverscriptBox[\"\[RightArrow]\", \"\[Placeholder]\"]\)\[Placeholder]"]},*)
		 {ButtonBox[SubscriptBox["\[Placeholder]", "\[Placeholder]"]], ButtonBox[SuperscriptBox["\[Placeholder]", "\[Placeholder]"]], ButtonBox[SubsuperscriptBox["\[Placeholder]", "\[Placeholder]", "\[Placeholder]"]],
		  ButtonBox["\!\(\*SubscriptBox[\"\[InvisiblePrefixScriptBase]\", \"\[Placeholder]\"]\)\[Placeholder]"], ButtonBox["\!\(\*SuperscriptBox[\"\[InvisiblePrefixScriptBase]\", \"\[Placeholder]\"]\)\[Placeholder]"], ButtonBox["\!\(\*SubsuperscriptBox[\"\", \"\[Placeholder]\", \"\[Placeholder]\"]\)\[Placeholder]"]}
		}, 
		RowSpacings -> 1.5, ColumnSpacings -> 1.5]],
        FontFamily -> "TimesNewRoman", FontWeight -> "Normal", FontSize -> 24]}, (*Bold*)
		Enabled->True,
        CellMargins -> {{0, 0}, {0, 0}},
		WindowMargins->{{Automatic,0},{Automatic,0}},
        Editable -> True, (*False*)
        ShowCellBracket -> False,
        WindowElements -> {},
        WindowFrame -> "Palette",
        WindowFrameElements -> "CloseBox",
        WindowSize -> {Fit, Fit},
        WindowToolbars -> {},
        WindowClickSelect -> False,
        WindowTitle -> "ReactionKinetics Palette",
		WindowFloating->True]
];

OpenReactionKineticsPalette[__] := (Message[OpenReactionKineticsPalette::"badarg"]; $Failed);

SyntaxInformation[OpenReactionKineticsNamesPalette] = {"ArgumentsPattern"->{}};

OpenReactionKineticsNamesPalette[] := 
	NotebookPut[Notebook[{
        Cell[BoxData[GridBox[
		{
			{ButtonBox["\"atommatrix\""]},
			{ButtonBox["\"atoms\""]},
			{ButtonBox["\"complexes\""]},
			{ButtonBox["\"deficiency\""]},
			{ButtonBox["\"externalspecies\""]},
			{ButtonBox["\"fhjgraphedges\""]},
			{ButtonBox["\"fhjweaklyconnectedcomponents\""]},
			{ButtonBox["\"fhjstronglyconnectedcomponents\""]},
			{ButtonBox["\"fhjterminalstronglyconnectedcomponents\""]},
			{ButtonBox["\"reactionsteporders\""]},
			{ButtonBox["\"reactionsteps\""]},
			{ButtonBox["\"species\""]},
			{ButtonBox["\"variables\""]},
			{ButtonBox["\"volpertgraph\""]},
		    {ButtonBox["Models"]},
			{ButtonBox["Reactions"]}
		},RowSpacings ->.9, ColumnSpacings ->1.]],
        FontFamily -> "TimesNewRoman", FontWeight -> "Normal", FontSize -> 22]}, (*Bold*)
		Enabled -> True,
        CellMargins -> {{0, 0}, {0, 0}},
		WindowMargins->{{Automatic,0},{0,Automatic}},
        Editable -> False,
        ShowCellBracket -> False,
        WindowElements -> {},
        WindowFrame -> "Palette",
        WindowFrameElements -> "CloseBox",
        WindowSize -> {Fit, Fit},
        WindowToolbars -> {},
        WindowClickSelect -> False,
        WindowTitle -> "ReactionKinetics Names",
		WindowFloating->True,WindowElements->{"VerticalScrollBar"}]]

OpenReactionKineticsNamesPalette[__] := (Message[OpenReactionKineticsNamesPalette::"badarg"]; $Failed);


ReactionRatesNotebook::usage = "ReactionRatesNotebook[{reactions},options] opens a new notebook where reaction rates can be assigned to the given reaction.";


ReactionRatesNotebook::badarg = "Illegal argument of function ReactionRatesNotebook.";


ReactionKinetics`Private`myReactionRatesCounter = 0;

Options[ReactionRatesNotebook] := {ExternalSpecies->{}};

SyntaxInformation[ReactionRatesNotebook]={"ArgumentsPattern"->{__,OptionsPattern[]}};

ReactionRatesNotebook[{reactions__}] := ReactionRatesNotebook[{reactions},{}];

ReactionRatesNotebook[{reactions__},externals:(_?VectorQ|OptionsPattern[])] :=
	(++ReactionKinetics`Private`myReactionRatesCounter;
	NotebookPut[Notebook[
		{
		Cell[TextData[{
			"Fill the following table in every placeholder with reaction rate coeficients or functions.\nEvaluate the cell when finished."}], 
			"Text", FontSize -> 15, FontColor -> RGBColor[0,0, 0], FontWeight -> Bold
			],
		Cell[BoxData[RowBox[{RowBox[{"myReactionRates" <> ToString[ReactionKinetics`Private`myReactionRatesCounter], "=",
				RowBox[{"Last", "/@",GridBox[
					Thread[List @@ Through[List[ToBoxes/@First[#]&,ToBoxes/@Array[Placeholder[{##}]&,Last[#]]&][
												Check[ReactionsData[{reactions},externals]["reactionsteps","R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]]]],
						RowLines -> True, RowSpacings -> 2.5, ColumnSpacings -> 3, ColumnAlignments -> {Center, Right}]
				}]}]}]], "Input"
			]
		},
		WindowSize -> {720, 500},
		WindowTitle -> "Reaction Rates Notebook " <> ToString[ReactionKinetics`Private`myReactionRatesCounter]
		]
	];)

ReactionRatesNotebook[___] := (Message[ReactionRatesNotebook::badarg]; $Failed)


(* ::Subsection::Closed:: *)
(*welcoming words*)


If[Length[Names["FLAGS`ECHOLOAD"]]>0,
	ReactionKinetics`Private`msgflag=ToExpression["FLAGS`ECHOLOAD"],
	ReactionKinetics`Private`msgflag=True;
];

If[ReactionKinetics`Private`msgflag,
	Print[Style["\!\(\*
		StyleBox[\"ReactionKinetics\",\nFontFamily->\"Courier New\"]\) Version "<> $ReactionKineticsVersionNumber <>" using Mathematica Version "<>
		$Version<>" (Version "<>ToString[$VersionNumber]<>", Release "<>ToString[$ReleaseNumber]<>")"<>
		" loaded "<>DateString[{"Day"," ","MonthName", " ", "Year", " at ", "Hour24",":", "Minute", " ", "TimeZone"}]
		<>"\nGNU General Public License (GPLv3) Terms Apply. ",FontWeight-> Bold,FontColor->Black]];
];

If[ReactionKinetics`Private`msgflag,
	Print[Style["Please report any issues, comments, complaint related to ReactionKinetics at\n",FontWeight-> Bold,FontColor->Black], 
			Hyperlink["jtoth@math.bme.hu", "mailto:jtoth@math.bme.hu"],Style[", "],
			Hyperlink["nagyal@math.bme.hu","mailto:nagyal@math.bme.hu"], Style[" or "],
			Hyperlink["dpapp@iems.northwestern.edu","mailto:dpapp@iems.northwestern.edu"]]
];

FLAGS`ECHOLOAD = ReactionKinetics`Private`msgflag;
Remove[ReactionKinetics`Private`msgflag];


(* ::Subsection::Closed:: *)
(*end*)


End[ ];


If[ $ReactionKineticsPackageLoaded =!= True, 
	If[ Length[Notebooks["ReactionKinetics Palette"]] === 0,
			OpenReactionKineticsPalette[];
	];
	If[ Length[Notebooks["ReactionKinetics Names"]] === 0, 
			OpenReactionKineticsNamesPalette[];
	];
];

$ReactionKineticsPackageLoaded = True;


SetAttributes[
	{
		$ReactionKinetics,
		$ReactionKineticsPackageLoaded,
		$ReactionKineticsVersionNumber
	},
	{Protected}
];


Protect[
	AvogadrosNumber,
	Global`aExplodator,
	Global`bExplodator,
	Global`fOregonator,
	Global`c,
	Global`t,
	Global`k,
	Global`\[DoubleStruckCapitalE],
	Global`g,
	Global`\[CapitalPi],
	Global`P,
	Global`\[Rho],
	Global`T,
	Global`X,
	Global`z,
	Bipartite,
	CombinatorialMoments,
	ComplexColors,
	Conditions,
	EdgeLabels,
	ExternalSpecies,
	FixedStepSize,
	FormattedOutput,
	GeneratedRateCoefficient,
	MassAction,
	Highlight,
	Indexed, 
	MaxIteration,
	Memo,
	Method,
	MolarGasConstant,
	Numbered,
	PlotFunction,
	Positivity,
	Side,
	Species,
	StronglyConnectedComponentsColors,
	SubgraphHighlight,
	TimeLimit,
	Verbose,
	Volume
];


SetAttributes[
	{
		Models,
		Reactions,
		BioModels,
		BioReactions,
		GetReaction,
		(**)
		ToCanonicalForm,
		ToReversible,
		ReactionsData,
		ReversibleQ,
		WeaklyReversibleQ,
		ShowFHJGraph,
		ShowVolpertGraph,
		AcyclicVolpertGraphQ,
		Atoms,
		ToAtomMatrix,
		AtomConservingQ,
		AtomsQ,
		FromAtomMatrix,
		DetailedBalanced,
		(**)
		MassActionKinetics,
		RightHandSide,
		DeterministicModel,
		Concentrations,
		StationaryPoints,
		GammaLeftNullSpace,
		EigensystemJacobian,
		AbsoluteConcentrationRobustness,
		(**)
		VolpertIndexing,
		(**)
		MassConservationRelations,
		MasterEquation,
		StationaryProbabilityDistributionEquation,
		ProbabilityGeneratingFunctionEquation,
		SolveProbabilityGeneratingFunctionEquation,
		MomentEquations,
		DependencyGraph,
		ShowDependencyGraph,
		Simulation,
		SimulationPlot,
		SimulationPlot2D,
		SimulationPlot3D,
		(**)
		OpenReactionKineticsPalette,
		OpenReactionKineticsNamesPalette,
		ReactionRatesNotebook,
		(**)
		CHEMKINExport,
		CHEMKINImport, 
		FilterReactions,
		FromStoichiometry, 
		DeleteAutocatalysis,
		ReversibleFHJRepresentation,
		MinFHJWeaklyConnectedComponents,
		MinFHJStronglyConnectedComponents,
		MaxFHJWeaklyConnectedComponents,
		MaxFHJStronglyConnectedComponents
		(*VolpertSpecies*)

	},
	{Protected,ReadProtected}
];


SetAttributes[
	{
	ContejeanDevie,
	CoveringDecompositionSet,
	(*DecompositionQ,*)
	Decompositions,
	(*ElementaryReactionQ,*)
	ElementaryReactions,
	FastSelection,
	Filter,
	(*GeneralizedDecompositionQ,*)
	GreedySelection,
	LPBased,
	MinimalDecompositions,
	ObjectiveFunction,
	Obligatory,
	Omittable,
	OriginalSelection,
	Preprocess,
	(*ReactionQ,*)
	SelectMinimalDecompositions,
(**)
	ZeroVectorQ 
    },
	{Protected,ReadProtected}
];


EndPackage[ ];
