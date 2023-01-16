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
(*Models*)


Get[FindFile["ReactionKineticsModels.wl"]]


(* ::Subsection::Closed:: *)
(*CHEMKIN*)


Get[FindFile["CHEMKIN.wl"]]


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


(* ::Subsubsection:: *)
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
(*stoichiometry*)


Get[FindFile["Stoichiometry.wl"]];


(* ::Subsection::Closed:: *)
(*deterministic approach*)


(* ::Subsubsection::Closed:: *)
(*MassActionKinetics, RightHandSide, DeterministicModel*)


MassActionKinetics::usage = "MassActionKinetics[{reactions},rratecoeffs,vars] builds up the mass action type kinetics attached to the reaction using \
rratecoeffs as reaction rate coefficients and vars as the names of independent variables of the kinetic function.";

MassActionKinetics::args = "The number of variables or that of reaction rate coefficients do not match with the dimensions of \[Alpha].";
MassActionKinetics::badarg = "Illegal argument of function MassActionKinetics.";

Options[MassActionKinetics] := Options[ReactionsData];

MassActionKinetics[{reactions__}, rratecoeffs_?VectorQ, opts:OptionsPattern[]] :=
	MassActionKinetics[
			Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["\[Alpha]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],
			rratecoeffs,
			Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]];

MassActionKinetics[{reactions__}, rratecoeffs_?VectorQ, variables_?VectorQ, opts:OptionsPattern[]] :=
	MassActionKinetics[
			Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["\[Alpha]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],
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

Options[RightHandSide] := Join[{MassAction->True}, Options[ReactionsData]];

RightHandSide[{reactions__}, rates_?VectorQ, opts : OptionsPattern[]] :=
	RightHandSide[{reactions},rates,Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],opts];

RightHandSide[{reactions__}, rates_?VectorQ, vars_?VectorQ, opts:OptionsPattern[]] :=
	Module[{alpha, gamma, species, nofr, nofs},
		{alpha, gamma, species, nofr, nofs} = Check[ReactionsData[{reactions},FilterRules[opts, ReactionsData]]["\[Alpha]","\[Gamma]","species","R","M"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
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

Options[DeterministicModel] := Join[{MassAction->True},Options[ReactionsData]];

DeterministicModel[{reactions__},s_Symbol:Global`t, opts : OptionsPattern[]] :=
	DeterministicModel[{reactions},Array[Subscript[Global`k,#]&,
		Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]],s,opts];

DeterministicModel[{reactions__},rates_?VectorQ,s_Symbol:Global`t, opts : OptionsPattern[]] :=
	DeterministicModel[{reactions},rates,
		Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],s,opts];

DeterministicModel[{reactions__},rates_?VectorQ,vars_?VectorQ,s_Symbol:Global`t, opts : OptionsPattern[]] :=
	Module[{v},
		v = Through[vars[s]];
		{Thread[D[v,s] == Check[RightHandSide[{reactions},rates,v,opts],Return[$Failed];,{RightHandSide::args,RightHandSide::badarg}]], v}
	];

(*Arrhenius reaction rate coefficients*)
DeterministicModel[{reactions__},arrhenius_?MatrixQ,enthalpies_?VectorQ,C0_?VectorQ,V_,Ta_,tres_,p1_,p2_,s_Symbol:Global`t,opts:OptionsPattern[]] :=
	Module[{kkk, maf, alpha, beta, vars, reacs},
		reacs = Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		If[reacs===Length[arrhenius] && reacs===Length[enthalpies],
			{alpha, beta, vars} = ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["\[Alpha]","\[Beta]","variables"];
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

Options[Concentrations] := Join[Options[DeterministicModel], Options[ReactionsData], Options[DSolve], Options[NDSolve]];

Concentrations::badarg = "Illegal argument of function Concentrations.";
Concentrations::init = "The initial vector has wrong size.";

(* Solving the ODE symbolically *)
Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,s_Symbol:Global`t,opts:OptionsPattern[]] :=
	Concentrations[{reactions},rates,init,0,s,opts];

Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,vars_?VariablesQ,s_Symbol:Global`t,opts:OptionsPattern[]] :=
	Concentrations[{reactions},rates,init,0,vars,s,opts];

Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,t0_?NumericQ(*Symbol*),s_Symbol:Global`t,opts:OptionsPattern[]] :=
	Module[{kineq, vars},
		{kineq, vars} = Check[DeterministicModel[{reactions},rates,s,FilterRules[opts,Options[DeterministicModel]]],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["M"],
			{
				vars,
				First[DSolve[Join[kineq, Thread[(vars /. s -> t0) == init]], vars, s, FilterRules[opts,Options[DSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,t0_?NumericQ(*Symbol*),vars_?VectorQ,s_Symbol:Global`t,opts:OptionsPattern[]] :=
	Module[{kineq, newvars},
		exs = Sequence@@FilterRules[{opts},Options[DeterministicModel]];
		{kineq, newvars} = Check[DeterministicModel[{reactions},rates,vars,s,FilterRules[opts,Options[DeterministicModel]]],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===ReactionsData[{reactions},exs]["M"],
			{
				newvars,
				First[DSolve[Join[kineq, Thread[(newvars /. s -> t0) == init]], newvars, s, FilterRules[opts,Options[DSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

(* Solving the ODE numerically *)
Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,{t0_?NumericQ,t1_?NumericQ},s_Symbol:Global`t,opts:OptionsPattern[]] :=
	Module[{kineq, vars},
		{kineq, vars} = Check[DeterministicModel[{reactions},rates,s,FilterRules[opts,Options[DeterministicModel]]],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["M"],
			{
				vars,
				First[NDSolve[Join[kineq, Thread[(vars /. s -> t0) == init]], vars, {s, t0, t1}, FilterRules[opts,Options[NDSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

Concentrations[{reactions__},rates_?VectorQ,init_?VectorQ,{t0_?NumericQ,t1_?NumericQ},vars_?VectorQ,s_Symbol:Global`t,opts:OptionsPattern[]] :=
	Module[{kineq, newvars},
		{kineq, newvars} = Check[DeterministicModel[{reactions},rates,vars,s,FilterRules[opts,Options[DeterministicModel]]],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["M"],
			{
				newvars,
				First[NDSolve[Join[kineq, Thread[(newvars /. s -> t0) == init]], newvars, {s, t0, t1}, FilterRules[opts,Options[NDSolve]]]]
			},
			Message[Concentrations::"init"];
			$Failed
		]
	];

(*Arrhenius reaction rate coefficients*)
Concentrations[{reactions__},arrhenius_?MatrixQ,enthalpies_?VectorQ,C0_?VectorQ,V_,Ta_,tres_,p1_,p2_,init_?VectorQ,{t0_?NumericQ,t1_?NumericQ},s_Symbol:Global`t,opts:OptionsPattern[]] :=
	Module[{kineq, vars},
		{kineq, vars} = Check[DeterministicModel[{reactions},arrhenius,enthalpies,C0,V,Ta,tres,p1,p2,s,FilterRules[opts,Options[DeterministicModel]]],Return[$Failed];,
									{DeterministicModel::badarg,DeterministicModel::args,DeterministicModel::steady,RightHandSide::args,RightHandSide::badarg}];
		If[Length[init]===(ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["M"]+1),
			{
				vars,
				First[NDSolve[Join[kineq, Thread[(vars /. s -> t0) == init]], vars, {s, t0, t1}, FilterRules[opts,Options[NDSolve]]]]
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

Options[StationaryPoints] := Join[{Positivity->False, Conditions->{}, Method->"Automatic", MassAction->True},Options[ReactionsData],Options[Reduce],Options[RightHandSide]];

StationaryPoints::badarg = "Illegal argument of function StationaryPoints.";
StationaryPoints::args = "One of the arguments `1`, `2` or `3` may have wrong shape.";
StationaryPoints::noeql = "The reaction simplex contains no stationary point.";
StationaryPoints::unknowm = "The system could not recognize the method `1`. Try Automatic or GammaLeftNullSpace.";

StationaryPoints[{reactions__},opts:OptionsPattern[]] :=
	Module[{r, vars},
		{r, vars} = Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["R","variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		StationaryPoints[{reactions},Array[Subscript[Global`k,#]&,r],Superscript[#,0]&/@vars,Superscript[#,"*"]& /@ vars,opts]
	];

StationaryPoints[{reactions__},rates_?VectorQ,opts:OptionsPattern[]] :=
	Module[{vars},
		vars = Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		StationaryPoints[{reactions},rates,Superscript[#,0]&/@vars,Superscript[#,"*"]& /@ vars,opts]
	];

StationaryPoints[{reactions__},rates_?VectorQ,init_?VectorQ,opts:OptionsPattern[]] :=
	StationaryPoints[
		{reactions},
		rates,
		init,
		Check[Superscript[#,"*"]& /@ ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["variables"],
			Return[$Failed];,
			{ReactionsData::wrreac,ReactionsData::badarg}
		],
		opts
	];

StationaryPoints[{reactions__},rates_?VectorQ,init_?VectorQ,vars_?VectorQ,opts:OptionsPattern[]]:=
	Module[{equilibriums, gamma, m, r, species, method, y, z, statp, balances, rateconst, inimass, nonneg, sol, fless, null, nullspacecond},

			{gamma, m, r, species} = Check[
				ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["\[Gamma]","M","R","species"],
				Return[$Failed];,
				{ReactionsData::wrreac,ReactionsData::badarg}
			];
			fless = If[OptionValue[Positivity],Less,LessEqual];

			If[r===Length[rates] && m===Length[init] && m===Length[vars],

				Switch[OptionValue[Method],

					"Automatic",(*default method*)

					z = Subscript[y,ToString[#]]& /@ Range[r];

					{statp, balances, rateconst, inimass, nonneg} =
						Thread/@{Check[RightHandSide[{reactions},rates,vars,FilterRules[opts,Options[RightHandSide]]],Return[$Failed];,{RightHandSide::args,RightHandSide::badarg}] == 0,
							gamma . z == vars-init, (rates /. Thread[species -> vars]) > 0, Cases[init,Except[0]]>0, fless[0,vars]};

					If[(sol = Reduce[Join[statp, balances, rateconst, inimass, nonneg, OptionValue[Conditions]], Join[vars,z], Reals, Backsubstitution->True, FilterRules[opts,Options[Reduce]]])===False,
						Message[StationaryPoints::noeql];
						Return[$Failed],
						Join[{vars},{Sort[Sort[MyToRules[Select[#/.And->List,Not[MemberQ[z,First[#]]]&]]]&/@Flatten[{sol/.Thread[rateconst->True]/.Thread[inimass->True]/.Thread[nonneg->True]/.Or->List}]]}]
					],

					"GammaLeftNullSpace",(*using the nullspace of gamma*)

					null = NullSpace[Transpose[gamma]] /. {}->{ConstantArray[0,m]};

					{statp, nullspacecond, rateconst, inimass, nonneg} =
						Thread/@{Check[RightHandSide[{reactions},rates,vars,FilterRules[opts,Options[RightHandSide]]],Return[$Failed];,{RightHandSide::args,RightHandSide::badarg}] == 0,
							null . (vars-init) == 0, (rates /. Thread[species -> vars]) > 0, Cases[init,Except[0]]>0,fless[0,vars]};

					sol := Reduce[Join[statp, nullspacecond, rateconst, inimass, nonneg, OptionValue[Conditions]], vars, Reals, Backsubstitution->True, FilterRules[opts,Options[Reduce]]];

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

Options[GammaLeftNullSpace] := Options[ReactionsData];

GammaLeftNullSpace::badarg = "Illegal argument of function GammaLeftNullSpace.";
GammaLeftNullSpace::args = "One of the arguments `1` or `2` may have wrong shape.";
GammaLeftNullSpace::null = "Left null space of \[Gamma] is null-dimensional.";

GammaLeftNullSpace[{reactions__},opts:OptionsPattern[]] :=
	Module[{lowerspecies},
		lowerspecies = ToLowerCase[Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]];
		GammaLeftNullSpace[{reactions}, ToExpression/@(#<>"0" &/@ lowerspecies), ToExpression/@lowerspecies, opts]
	];

GammaLeftNullSpace[{reactions__},init_?VectorQ,opts:OptionsPattern[]] :=
	GammaLeftNullSpace[{reactions},init,Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["variables"],Return[$Failed];,
													{ReactionsData::wrreac,ReactionsData::badarg}],opts];

GammaLeftNullSpace[{reactions__},init_?VectorQ,vars_?VectorQ,opts:OptionsPattern[]] :=
	Module[{gamma, m, ns, conds, sol},
			{gamma, m} = Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["\[Gamma]","M"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

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

Options[MassConservationRelations] := Options[ReactionsData];

MassConservationRelations::badarg = "Illegal argument of function MassConservationRelations.";

MassConservationRelations[{reactions__},opts:OptionsPattern[]] :=
	MassConservationRelations[{reactions}, Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["variables"],Return[$Failed];,
													{ReactionsData::wrreac,ReactionsData::badarg}] /. (Global`c -> Global`\[Rho]), opts];

MassConservationRelations[{reactions__},vars_?VectorQ,opts:OptionsPattern[]] :=
	Module[{gamma, ropts},
			gamma = Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["\[Gamma]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
			ropts = FilterRules[{opts},Options[Reduce]];

			Reduce[Flatten[Thread/@{Transpose[gamma] . vars==0, vars>0}], vars, Reals, Backsubstitution->True, ropts] /. And -> List /. Or -> List
	];

MassConservationRelations[___] := (Message[MassConservationRelations::"badarg"]; $Failed)


EigensystemJacobian::usage = "EigensystemJacobian[{reactions},rates,vars] returns the eigensystem of the Jacobian matrix of the right-hand side of the \
induced kinetic differential equation endowed with some kinetics given by the rates.";

EigensystemJacobian::badarg = "Illegal argument of function EigensystemJacobian.";
EigensystemJacobian::args = "One of the arguments `1` or `2` may have wrong shape.";

Options[EigensystemJacobian] := Join[{MassAction->True},OptionsPattern[ReactionsData],Options[RightHandSide]];

EigensystemJacobian[{reactions__}, opts : OptionsPattern[]] :=
	Module[{r, vars},
		{r, vars} = Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["R","variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		EigensystemJacobian[{reactions},Array[Subscript[Global`k,#]&,r],vars,opts]
	];

EigensystemJacobian[{reactions__},rates_?VectorQ,opts : OptionsPattern[]] :=
	EigensystemJacobian[{reactions},rates,
			Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["variables"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],opts];

EigensystemJacobian[{reactions__},rates_?VectorQ,vars_?VectorQ, opts : OptionsPattern[]] :=
	Module[{m, r},
		{m, r} = Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["M","R"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		If[r===Length[rates] && m===Length[vars],

			Eigensystem[D[Check[RightHandSide[{reactions},rates,vars,FilterRules[opts,Options[RightHandSide]]],Return[$Failed];,{RightHandSide::args,RightHandSide::badarg}],{vars}]],

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

Options[AbsoluteConcentrationRobustness] := Join[{TimeConstraint -> 120}, Options[ReactionsData],Options[StationaryPoints]];

AbsoluteConcentrationRobustness[{reactions__},opts:OptionsPattern[]] :=
	Module[{delta, fhjstrong, fhjterminal, nonterms, acr},

		{delta, fhjstrong, fhjterminal} =
			Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["deficiency","fhjstronglyconnectedcomponents","fhjterminalstronglyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		If[ Last[StringSplit[delta,"="]]=!="1",

			Message[AbsoluteConcentrationRobustness::"maybe","Deficiency does not equal to 1"];
			Return[{reactions}],

			nonterms = Flatten[fhjstrong/.Thread[fhjterminal->Sequence[]]] /. "0"->0 /. Thread[OptionValue[ExternalSpecies]->0];
			acr = Flatten[(nonterms-# &/@ nonterms) /. Plus[A_,B_]->0 /. 0->{}]; (*/. r_*A_:>A;*)

			If[ acr === {},

				Message[AbsoluteConcentrationRobustness::"maybe","Deficiency is one but there do not exist nonterminal complexes differing only in one species"];
				Return[{reactions}],

				Switch[TimeConstrained[Check[Length[StationaryPoints[{reactions},Positivity->True,FilterRules[opts, StationaryPoints]]]>=2,False,
								{StationaryPoints::unknowm,StationaryPoints::args,StationaryPoints::badarg,StationaryPoints::noeql}],OptionValue[TimeConstraint],"Abort"],
						"Abort",
							Message[AbsoluteConcentrationRobustness::"timeexc",OptionValue[TimeConstraint]];
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


End[];


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
		(**)
		ToCanonicalForm,
		ToReversible,
		ReactionsData,
		ReversibleQ,
		WeaklyReversibleQ,
		FHJGraph,
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
