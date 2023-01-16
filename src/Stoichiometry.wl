(* ::Package:: *)

(* ::Subsubsection::Closed:: *)
(*ToCanonicalForm*)


ToCanonicalForm::usage = "ToCanonicalForm[reactions] returns the canonical form of reactions \
using only right, left and left-right arrows.";
ToCanonicalForm::badarg = "Illegal argument of function ToCanonicalForm.";


SyntaxInformation[ToCanonicalForm]={"ArgumentsPattern"->{__}};
(* Reactions is from ReactionKineticsModels *)
ToCanonicalForm[reactions__]:=
	Flatten[{reactions} /. ReactionKineticsModels`Reactions] /.
					  Join[Thread[{LeftRightArrow,DoubleLeftRightArrow,LongLeftRightArrow,DoubleLongLeftRightArrow,RightArrowLeftArrow,TwoWayRule,UndirectedEdge}->Equilibrium],
									{LeftArrowRightArrow->ReverseEquilibrium},
									Thread[{Rule,ShortRightArrow,DoubleRightArrow,LongRightArrow,DoubleLongRightArrow,DirectedEdge}->RightArrow],
									Thread[{ShortLeftArrow,DoubleLeftArrow,LongLeftArrow,DoubleLongLeftArrow}->LeftArrow]];
ToCanonicalForm[___] := (Message[ToCanonicalForm::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*ReactionsData*)


SpeciesQ := UpperCaseQ[StringTake[ToString[#],1]]&;


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


RdataOptions = {ExternalSpecies->{},InternalSpecies->{}};
Rdata[{reactions__},opts:OptionsPattern[RdataOptions]] := Rdata[{reactions},opts] =
	Module[
			{ reacs,
			  m, builtin, externals, allExternals, exsrule,
			  reactionsteps, internalReactionSteps, fhjgraph, complexes, vgggraph, vgraph, vggraphwzerocomplex, species,
			  indexed, a, b, e, nloc, lloc, sloc,
			  lsp, lrs, alpha, beta, gamma, fhjweakvertices, fhjweakcomponents,
			  fhjstrongvertices, fhjstrongcomponents,
			  fhjterminalstrongvertices, fhjterminalstrongcomponents
			},

			reacs = Flatten[{reactions}];

			If[(m = Cases[reacs,_String])=!={},

				builtin = Check[(#->GetReaction[#])& /@ m, Return[$Failed];, GetReaction::"nvmod"];(**)
				Rdata[DeleteDuplicates[Join[Flatten[reacs /. builtin], reacs /. Thread[m->Sequence[]]]], opts],

				reactionsteps = ReactionsToList[reacs] /. (0->"0"); (*"0" - zero complex*)
				If[OptionValue[ExternalSpecies]=!={} && OptionValue[InternalSpecies]=!={},
				    (Message[ReactionsData::"internalexternal"]; Return[$Failed]),
					If[OptionValue[InternalSpecies]=!={},
						Block[{allSpecies,allComplexes},
							allComplexes = DeleteDuplicates[Flatten[reactionsteps /. RightArrow -> List]];
							allSpecies=complextospecies[allComplexes] /. "0" -> Sequence[];
							externals=Complement[allSpecies,OptionValue[InternalSpecies]]
						],
						externals=OptionValue[ExternalSpecies];
					]
				];
				allExternals = DeleteDuplicates[Prepend[externals,"0"]/.(0->"0")]; (*internalspecies megadas*)
				exsrule = Thread[allExternals -> 0];

				internalReactionSteps = DeleteDuplicates[reactionsteps /. exsrule /. (0->"0")];

				fhjgraph = Graph[Rule@@@internalReactionSteps];

				complexes = DeleteDuplicates[Flatten[internalReactionSteps /. RightArrow -> List]]; (*Sort*)

				vgggraph = Flatten[(Thread /@ {#->complextocoefflist[Last[#]], complextocoefflist[First[#]]->#}) &/@ internalReactionSteps];
				vgraph = vgggraph /. Rule[A_RightArrow,B_]:>{Rule[A, First[B]], Last[B]} /. Rule[A_,B_RightArrow]:>{Rule[First[A], B], Last[A]};

				species = complextospecies[complexes] /. "0" -> Sequence[]; (*komplexek csak belso anyagfajta erejeig egyertelmuek*)

				indexed = Flatten[MapIndexed[#1->First[#2]&, #]& /@ {species, internalReactionSteps}];

				{lsp, lrs} = Length /@ {species, internalReactionSteps};

				vggraphwzerocomplex = vgraph /. {Rule["0",y_RightArrow],z_} :> Sequence[] /. {Rule[x_RightArrow,"0"],z_} :> Sequence[];

				a = Cases[vggraphwzerocomplex, {Rule[x_,y_RightArrow],z_}:> Rule[{x,y} /. indexed, z]];
				b = Cases[vggraphwzerocomplex, {Rule[x_RightArrow,y_],z_}:> Rule[{y,x} /. indexed, z]];
				(*Print[{Join[{a,b}],lsp,lrs}];*)
				alpha = SparseArray[a, {lsp, lrs}];
				beta = SparseArray[b, {lsp, lrs}];
				gamma = beta-alpha;
				(*SparseArray[Join[a,b] //. {x___,Rule[{y_,z_},k_],u___,Rule[{y_,z_},l_],v___}:>{x,Rule[{y,z},k+l],u,v}, {lsp,lrs}];*)

				fhjweakvertices = ConnectedComponents[UndirectedGraph[fhjgraph]];
				fhjweakcomponents = subgraphedges[fhjgraph, fhjweakvertices];

				fhjstrongvertices = ConnectedComponents[fhjgraph];
				fhjstrongcomponents = subgraphedges[fhjgraph, fhjstrongvertices];
				fhjterminalstrongvertices = Map[
												Complement[
													Union @@ Map[Cases[internalReactionSteps, RightArrow[#,y_]:>y]&, #], #
												]==={} /. True-># /. False->Sequence[]&, fhjstrongvertices
											];
				fhjterminalstrongcomponents = subgraphedges[fhjgraph, fhjterminalstrongvertices];

				Dispatch[{
					"species" -> species,
					"M" -> lsp,
					"externalspecies" -> (e = allExternals /. "0" -> Sequence[]),
					(*e = If[MemberQ[complexes,"0"], allExternals, Rest[allExternals]]*)
					"E" -> Length[e],
					"reactionsteps" -> reactionsteps, (*steprule*)
					"R" -> lrs,
					"complexes" -> complexes,
					"fhjgraphedges" -> Rule@@@internalReactionSteps,
					"fhjweaklyconnectedcomponents" -> fhjweakcomponents,
					"fhjstronglyconnectedcomponents" -> fhjstrongcomponents,
					"fhjterminalstronglyconnectedcomponents" -> fhjterminalstrongcomponents,
					"N" -> (nloc = Length[complexes]),
					"L" -> (lloc = Length[fhjweakvertices]),
					"S" -> (sloc = MatrixRank[gamma]),
					(*MatrixRank[gamma = First[Differences[ab = (Normal/@{a,b})]]]*)
					"deficiency" -> "\[Delta]=N-L-S="<>ToString[nloc]<>"-"<>ToString[lloc]<>"-"<>ToString[sloc]<>"="<>ToString[nloc-lloc-sloc],
					"\[Alpha]" -> alpha,
					"\[Beta]" -> beta,
					"\[Gamma]" -> gamma,
					"reactionsteporders" -> (Total /@ Transpose[alpha]),
					"variables" -> (Subscript["c", #] & /@ species),
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
ReactionsData::internalexternal = "Only one of InternalSpecies and ExternalSpecies can be present.";
ReactionsData::badname = "At least one of the arguments '`1`' is a non-identified property. Try \
ReactionsData[\"Properties\"].";
ReactionsData::wrreac = "Argument `1` may not be in a correct form. Check OpenReactionKineticsPalette[].";


SyntaxInformation[ReactionsData]={"ArgumentsPattern"->{__,OptionsPattern[]}};
Options[ReactionsData] = RdataOptions;
ReactionsData["Properties"] := PropertiesReactionsData;
ReactionsData[{},OptionsPattern[]][] := {};
ReactionsData[{reactions__},opts:OptionsPattern[]][] :=
	Check[Rdata[{reactions},opts],
					Message[ReactionsData::"wrreac",ReactionsForm[{reactions}]];
					Return[$Failed]
	];
ReactionsData[{},OptionsPattern[]]["Properties"] := PropertiesReactionsData;
ReactionsData[{reactions__},OptionsPattern[]]["Properties"] := PropertiesReactionsData;
ReactionsData[{},OptionsPattern[]]["All"] := {};
ReactionsData[{reactions__},opts:OptionsPattern[]]["All"] := ReactionsData[{reactions},opts][PropertiesReactionsData];
ReactionsData[{},opts:OptionsPattern[]][data__?StringQ] := (
	If[Nand @@ Map[MemberQ[PropertiesReactionsData,#]&,{data}],
		Message[ReactionsData::"badname",data];
	];
	Return[{}]
);
ReactionsData[{reactions__},opts:OptionsPattern[]][data__?StringQ] := (
	If[Nand @@ Map[MemberQ[PropertiesReactionsData,#]&,{data}],
		Message[ReactionsData::"badname",data];
	];
	Catch[
		Part[
			ReplaceAll[{data},
					Check[Rdata[{reactions},opts],
							Throw[
								Message[ReactionsData::"wrreac",ReactionsForm[{reactions}]];
								Return[$Failed]
							]
					]
			], Length[{data}] /. x_/;x>1 :> (1;;x)
		]
	]
);
ReactionsData[{},OptionsPattern[]][data__] := ReactionsData[{}][Sequence @@ ToString /@ Flatten[{data}]];
ReactionsData[{reactions__},opts:OptionsPattern[]][data__] :=
		Block[{d, pos, rdata, sol = {data}},

				d = Flatten[sol];
				rdata = ReactionsData[{reactions},opts][Sequence @@ ToString /@ d];
				If[Length[d]===1, rdata = {rdata}];
				pos = Flatten[Position[sol, #]& /@ d, 1];
				Do[Part[sol, Sequence@@pos[[n]]] = rdata[[n]], {n, Length[d]}];

				If[ Length[sol] === 1,
					Flatten[sol, 1],
					sol
				]
		];
Format[ReactionsData[{reactions__},opts:OptionsPattern[]]] :=
	DynamicModule[{x, rd},
					x = Dynamic[Refresh[Round[Clock[Infinity]],UpdateInterval->1]];
					Monitor[rd = ReactionsData[{reactions},opts][];,
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
(*ToReversible, ReversibleFHJRepresentation, ReversibleQ, WeaklyReversibleQ*)


ToReversible::usage = "ToReversible[reactions] make all the reaction steps of reactions reversible.";
ToReversible::badarg = "Illegal argument of function ToReversible.";


SyntaxInformation[ToReversible]={"ArgumentsPattern"->{__}};
ToReversible[reactions__] := ToCanonicalForm[reactions]/.RightArrow->Equilibrium/.LeftArrow->Equilibrium;
ToReversible[___] := (Message[ToReversible::"badarg"]; $Failed)


ReversibleFHJRepresentation::usage = "ReversibleFHJRepresentation[{reactions},options] results in the Feinberg-Horn-Jackson \
graph making all its edges reversible.";
ReversibleFHJRepresentation::badarg = "Illegal argument of function ReversibleFHJRepresentation.";


Options[ReversibleFHJRepresentation] = Options[ReactionsData];
SyntaxInformation[ReversibleFHJRepresentation]={"ArgumentsPattern"->{__,OptionsPattern[]}};
ReversibleFHJRepresentation[{reactions__},opts:OptionsPattern[]] :=
	ReplaceRepeated[
		Check[ReactionsData[{reactions},opts]["fhjgraphedges"],Abort[];,{ReactionsData::wrreac,ReactionsData::badarg}],
		{x___, Rule[y_,z_], t___, Rule[z_,y_], u___} :> {x, Equilibrium[y,z], t, u}
	];

(*avagy: Complement[EdgeList[g], EdgeList[ReverseGraph[g]]], where g = FHJ Graph object*)
ReversibleFHJRepresentation[___] := (Message[ReversibleFHJRepresentation::"badarg"]; $Failed)


ReversibleQ::usage = "ReversibleQ[{reactions},options] yields True if and only if the reaction is reversible.";
ReversibleQ::badarg = "Illegal argument of function ReversibleQ.";


SyntaxInformation[ReversibleQ]={"ArgumentsPattern"->{__}};
Options[ReversibleQ] = Options[ReversibleFHJRepresentation];
ReversibleQ[{reactions__},opts:OptionsPattern[]] :=
		FreeQ[ReversibleFHJRepresentation[{reactions},opts],_Rule];
ReversibleQ[___] := (Message[ReversibleQ::"badarg"]; $Failed)


WeaklyReversibleQ::usage = "WeaklyReversibleQ[{reactions},options] yields True if and only if the reaction is weakly reversible.";
WeaklyReversibleQ::badarg = "Illegal argument of function WeaklyReversibleQ.";


SyntaxInformation[WeaklyReversibleQ]={"ArgumentsPattern"->{__}};
Options[WeaklyReversibleQ] = Options[ReactionsData];
WeaklyReversibleQ[{reactions__},opts:OptionsPattern[]] :=
		(Length[Check[ReactionsData[{reactions},opts]["fhjstronglyconnectedcomponents"],Abort[];,{ReactionsData::wrreac,ReactionsData::badarg}]] === 1);
WeaklyReversibleQ[___] := (Message[WeaklyReversibleQ::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*FilterReactions, FromStoichiometry, DeleteAutocatalysis*)


FilterReactions::usage = "FilterReactions[{reactions},species,options] filters all the reaction steps in which \
the given set of species shows up as either reactant or product or (reactant or product) species. \
This can be controlled by option Side \[Rule] \"Reactant\", \"Product\" or \"All\".";

FilterReactions::badsp = "Some of the species may be not contained in the list of internal species of the reaction.";
FilterReactions::badopt = "Option Side only takes three values: \"Reactant\", \"Product\" or \"All\".";
FilterReactions::badarg = "Illegal argument of function FilterReactions.";

Options[FilterReactions]:=Join[{Side -> "All"}, Options[ReactionsData]];

FilterReactions[{reactions__},specs_?VectorQ,opts : OptionsPattern[]]:=
	Module[{ rdata, fhj, sp, spp },

		rdata = Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["fhjgraphedges","species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];
		If[rdata === $Failed, Return[$Failed];];
		{fhj, sp} = rdata;

		If[(spp = Intersection[specs,sp]) =!= {},

			If[ MemberQ[{"Reactant","Product","All"}, OptionValue[Side]],
				Switch[OptionValue[Side],
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
	FromStoichiometry[alpha, beta, Array["X"[#]&,Length[alpha]]];

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

Options[DeleteAutocatalysis] := Options[ReactionsData];

DeleteAutocatalysis[{reactions__}, opts:OptionsPattern[]] :=
	DeleteAutocatalysis @@ Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["\[Alpha]","\[Beta]","species"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

DeleteAutocatalysis[{reactions__}, variables_?VectorQ, opts:OptionsPattern[]] :=
	DeleteAutocatalysis[Sequence @@ Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["\[Alpha]","\[Beta]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}],variables];

DeleteAutocatalysis[alpha_?MatrixQ,beta_?MatrixQ] :=
	DeleteAutocatalysis[alpha, beta, Array["X"[#]&,Length[alpha]]];

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

Options[MinFHJWeaklyConnectedComponents] := Options[ReactionsData];

MinFHJWeaklyConnectedComponents[{reactions__}, opts:OptionsPattern[]] :=
	Module[{fhjwcc, lengthvtxes, minvtxes},

		fhjwcc = Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["fhjweaklyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		lengthvtxes = Length[Last[#]] &/@ fhjwcc;

		minvtxes = Min[lengthvtxes];

		(*fhjwcc[[Flatten[Position[lengthvtxes, minvtxes]]]]*)
		Select[fhjwcc, Length[Last[#]]===minvtxes &]

	];

MinFHJWeaklyConnectedComponents[___][___] := (Message[MinFHJWeaklyConnectedComponents::"badarg"]; $Failed)


MinFHJStronglyConnectedComponents::usage = "MinFHJStronglyConnectedComponents[{reactions},options] returns the smallest strongly connected components of the given reaction.";

Options[MinFHJStronglyConnectedComponents] := Options[ReactionsData];

MinFHJStronglyConnectedComponents[{reactions__}, opts:OptionsPattern[]] :=
	Module[{fhjscc, lengthvtxes, minvtxes},

		fhjscc = Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["fhjstronglyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		lengthvtxes = Length[Last[#]] &/@ fhjscc;

		minvtxes = Min[lengthvtxes];

		(*fhjscc[[Flatten[Position[lengthvtxes, minvtxes]]]]*)
		Select[fhjscc, Length[Last[#]]===minvtxes &]

	];

MinFHJStronglyConnectedComponents[___][___] := (Message[MinFHJStronglyConnectedComponents::"badarg"]; $Failed)


MaxFHJWeaklyConnectedComponents::usage = "MaxFHJWeaklyConnectedComponents[{reactions},options] returns the largest weakly connected components of the given reaction.";

Options[MaxFHJWeaklyConnectedComponents] := Options[ReactionsData];

MaxFHJWeaklyConnectedComponents[{reactions__}, opts:OptionsPattern[]] :=
	Module[{fhjwcc, lengthvtxes, maxvtxes},

		fhjwcc = Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["fhjweaklyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		lengthvtxes = Length[Last[#]] &/@ fhjwcc;

		maxvtxes = Max[lengthvtxes];

		(*fhjwcc[[Flatten[Position[lengthvtxes, maxvtxes]]]]*)
		Select[fhjwcc, Length[Last[#]]===maxvtxes &]

	];

MaxFHJWeaklyConnectedComponents[___][___] := (Message[MaxFHJWeaklyConnectedComponents::"badarg"]; $Failed)


MaxFHJStronglyConnectedComponents::usage = "MaxFHJStronglyConnectedComponents[{reactions},options] returns the largest strongly connected components of the given reaction.";

Options[MaxFHJStronglyConnectedComponents] := Options[ReactionsData];

MaxFHJStronglyConnectedComponents[{reactions__}, opts:OptionsPattern[]] :=
	Module[{fhjscc, lengthvtxes, maxvtxes},

		fhjscc = Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["fhjstronglyconnectedcomponents"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

		lengthvtxes = Length[Last[#]] &/@ fhjscc;

		maxvtxes = Max[lengthvtxes];

		(*fhjscc[[Flatten[Position[lengthvtxes, maxvtxes]]]]*)
		Select[fhjscc, Length[Last[#]]===maxvtxes &]

	];

MaxFHJStronglyConnectedComponents[___][___] := (Message[MaxFHJStronglyConnectedComponents::"badarg"]; $Failed)


(* ::Subsubsection::Closed:: *)
(*FHJGraph, ShowFHJGraph*)


FHJGraph::usage = "FHJGraph[{reactions},options] returns the Feinberg-Horn-Jackson graph of the reaction \
of the reaction with options from ReactionsData.";
Options[FHJGraph] = Options[ReactionsData];

FHJGraph[{reactions__}, opts:OptionsPattern[]] := Graph[ReactionsData[{reactions}, FilterRules[opts, Options[ReactionsData]]]["fhjgraphedges"]]


GraphPlotFunctionQ := MemberQ[{"Graph","GraphPlot","LayeredGraphPlot","GraphPlot3D","TreePlot"},#]&;


ShowFHJGraph::usage = "ShowFHJGraph[{reactions},rratecoeffs,options] displays the Feinberg-Horn-Jackson graph of the reaction \
with rratecoeffs as reaction rate coeffcients.";

ShowFHJGraph::badarg = "Illegal argument of function ShowFHJGraph.";
ShowFHJGraph::funcarg = "Argument `1` is not a valid graph plot function. Try Graph, GraphPlot, GraphPlot3D, LayeredGraphPlot or TreePlot.";
ShowFHJGraph::ccols = "The list of colors has wrong shape.";
ShowFHJGraph::srates = "The list of reaction rate coefficients has wrong shape.";
ShowFHJGraph::sccomps = "The list of colors (for strongly connected components) may have wrong shape or complex colors are also given.";

Options[ShowFHJGraph]:=Join[
  {ComplexColors -> {}, PlotFunction -> "GraphPlot", StronglyConnectedComponentsColors -> {}, Numbered -> False},
   Options[ReactionsData],
   Options[GraphPlot],
   Options[LayeredGraphPlot],
   Options[GraphPlot3D],
   Options[TreePlot]
];

(*SyntaxInformation[ShowFHJGraph]={"ArgumentsPattern"->{{__},___}};*)

ShowFHJGraph[{reactions__}, rates_List:{}, opts:OptionsPattern[]] :=
	Module[{ cc, externals, rdata, exs,
			 n, r, fhj, cxs, fhjscc, graphfunc, nmd,
			 nopt, c, sc, rscc, rs
		   },

			cc = FilterRules[opts,Options[ShowFHJGraph]];

			rdata = Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["N","R","fhjgraphedges","complexes","fhjstronglyconnectedcomponents"], Return[$Failed], {ReactionsData::wrreac,ReactionsData::badarg}];
			If[rdata === $Failed, Return[$Failed]];

			{n, r, fhj, cxs, fhjscc} = rdata;

			graphfunc = OptionValue[PlotFunction];
			If[GraphPlotFunctionQ[graphfunc],
				graphfunc = ToExpression[graphfunc],
				Message[ShowFHJGraph::"funcarg",graphfunc];
				Return[$Failed]
			];

			fhjscc = Last /@ fhjscc;
			If[OptionValue[Numbered]===True,
				cxs = Thread[Rule[cxs,Range[n]]];
				fhj = fhj /. cxs;
				fhjscc = fhjscc /. cxs;
				cxs = Last /@ cxs;
			];

			nopt = FilterRules[opts,Options[graphfunc]];
			c = OptionValue[ComplexColors];
			sc = OptionValue[StronglyConnectedComponentsColors];

			If[sc =!= {},

				If[(VectorQ[sc] && Length[fhjscc] === Length[sc]) && c === {},
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

				If[c =!= {},

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

Options[VolpertIndexing] := Join[{Verbose -> False}, Options[ReactionsData]];

VolpertIndexing[{reactions__}, initspecies_?VectorQ, specsreacs_List : {}, opts:OptionsPattern[]]:=
	Module[{
			rdata, m, r, sp, reacs, vggraph, vgraph, w, w0, indexed, complx, initcomplx,
			indexedr, indexedsp, a, b, v, req, z1, z2, t1, t2, vindeces
		   },

		rdata = Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["M","R","species","fhjgraphedges","volpertgraph"],Return[$Failed],
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

Options[ShowVolpertGraph] := Join[
  {PlotFunction->"GraphPlot", EdgeLabels -> False, Highlighted -> {}, Indexed -> False, Numbered -> False, GraphLayout->"BipartiteEmbedding"},
  Options[ReactionsData],
  Options[GraphPlot],
  Options[LayeredGraphPlot],
  Options[GraphPlot3d],
  Options[TreePlot]
 ];

ShowVolpertGraph::badarg = "Illegal argument of function ShowVolpertGraph.";
ShowVolpertGraph::funcarg = "Argument `1` is not a valid graph plot function. Try Graph, \
GraphPlot, GraphPlot3D, LayeredGraphPlot or TreePlot.";

ShowVolpertGraph[{reactions__},opts:OptionsPattern[]] :=
	Module[{
			graphfunc, vgp, vindices, vtempind,
			vg, vind, complx, numbd, ind, initsp, highlighted, edgelabels,
			nopts, g, bipd, fhj, species, m, r
		   },

		graphfunc = OptionValue[PlotFunction];

		If[GraphPlotFunction2Q[graphfunc],

			vgp = Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["volpertgraph"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

			If[ vgp === $Failed, Return[$Failed];, {vg, vind} = vgp;];

			highlighted = OptionValue[Highlighted]; (*ToCanonicalForm*)

			ind = OptionValue[Indexed];
			species = ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["species"];
			initsp = Intersection[Flatten[{ind}],species];
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

			numbd = OptionValue[Numbered] === True;
			If[ numbd,
					vg = vg /. vind;
					highlighted = highlighted /. vind;
			];

			edgelabels = OptionValue[EdgeLabels] === True;

			If[ highlighted =!= {},

				If[ edgelabels, vg = Labeled@@# &/@ vg;, vg = First /@ vg; ];
				(*ez itt erzeketlen a plotfunction-ra es a bipartite-re*)
				g = Graph[vg];
				highlighted = Subgraph[g,highlighted];
				HighlightGraph[g, highlighted, FilterRules[opts,Options[HighlightGraph]]],

				nopts = FilterRules[{opts},Options[ToExpression[graphfunc]]];
				Switch[graphfunc,
						"Graph",
							If[ edgelabels, vg = Labeled@@# &/@ vg;, vg = First /@ vg; ];
							(*ez meg itt erzeketlen a bipartite-re*)
							Graph[vg, nopts],
						_,
							bipd = OptionValue[GraphLayout] === "BipartiteEmbedding";
							nopts = Join[(FilterRules[nopts, EdgeLabeling] /. {}->{EdgeLabeling->False}), nopts];
							If[ edgelabels, nopts = Join[{EdgeLabeling -> True}, nopts]; ];
							If[ bipd,
								fhj = ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["fhjgraphedges"] /. Rule->RightArrow;
								{m, r} = ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["M","R"];
								complx = ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["complexes"];
								If[MemberQ[complx,"0"], species = Append[species, "0"]; m = m+1;];
								If[ ind =!= False,
									species = species /. Flatten[vindices];
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
											Table[(species[[i]])->{0,i},{i,1,m,1}],Table[(fhj[[i]])->{Ceiling[Min[{m,r}]/2],i},{i,1,r,1}]
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

Options[AcyclicVolpertGraphQ] := Options[ReactionsData];

AcyclicVolpertGraphQ[{reactions__}, opts:OptionsPattern[]] :=
	AcyclicGraphQ[Graph[First/@First[Check[ReactionsData[{reactions},FilterRules[opts,Options[ReactionsData]]]["volpertgraph"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}]]]];

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
         /@ ({ Subscript["L", #]& /@ Range[Length[linkageClasses]],linkageClasses}//Transpose),
         2
 ]// Transpose;
 UndirectedGraph[edges,EdgeLabels->Table[edges[[i]]->edgeLabels[[i]],{i, Length[edges]}]]
]


ShowSCLGraph::usage = "Given a ReactionData[reactions] draws the SCL graph of the mechanism. The linkage classes are denoted by numbers. The number of a vertex represents\
 the index of its linkage class in WeaklyConnectedGraphComponents[reactionData[\"fhjgraphedges\"]]";

ShowSCLGraph::badarg = "Illegal argument of function ShowSCLGraph.";
Options[ShowSCLGraph] = Options[GraphPlot];

ShowSCLGraph[reactionData_,opts:OptionsPattern[]]:= GraphPlot[SCLGraph[reactionData,opts], FilterRules[opts, Options[GraphPlot]], GraphLayout->"BipartiteEmbedding"];


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

Options[AtomConservingQ]:=Options[ReactionsData];

SyntaxInformation[AtomConservingQ]={"ArgumentsPattern"->{{__},___}};

AtomConservingQ[{reactions__},opts:OptionsPattern[]] :=
	Module[{species, gamma, atommatrix, nofatch},
		{species, gamma} = Check[ReactionsData[{reactions},opts]["species","\[Gamma]"],Return[$Failed];,{ReactionsData::wrreac,ReactionsData::badarg}];

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

Options[DetailedBalanced] := Join[Options[ReactionsData], {GeneratedRateCoefficient -> "k"}];

SyntaxInformation[DetailedBalanced]={"ArgumentsPattern"->{{__},___}};

DetailedBalanced[{reactions__}, opts : OptionsPattern[]] :=
	Block[ { genpar, rr },

		genpar = OptionValue[GeneratedRateCoefficient];

		If[ Head[genpar] === Symbol || Head[genpar] === String,

			rr = Check[ReactionsData[{reactions}, FilterRules[opts,Options[ReactionsData]]]["R"], Return[$Failed],
							{ReactionsData::wrreac,ReactionsData::badarg}];

			(*If[rr === $Failed, Return[$Failed]];*)

			DetailedBalanced[{reactions}, Array[Subscript[genpar,#]&,rr], opts],

			Message[DetailedBalanced::genparerr];
			Return[$Failed];

		]
	];

DetailedBalanced[{reactions__}, rates_?VectorQ, opts : OptionsPattern[]] :=
	Module[{ cxes, m, r, nn, ll, fhj, gamma, deficiency, fhjunlist, fhjun, rsteplistind, rrates, shadow, shw,
			 fhjun0, mintree, mintreeun, c0, circs, spanreacs, mm, ccc, circconds, spanforconds, rev
		   },

			If[Not[ReversibleQ[{reactions},FilterRules[opts, Options[ReversibleQ]]]],
				Message[DetailedBalanced::"notrev"];
				Return[$Failed];
			];

			{cxes, m, r, nn, ll, fhj, gamma, deficiency} =
					Check[ReactionsData[{reactions},FilterRules[opts, Options[ReactionsData]]]["complexes","M","R","N","L","fhjgraphedges","\[Gamma]","deficiency"], Return[$Failed],
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
