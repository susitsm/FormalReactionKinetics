(* ::Package:: *)

BeginPackage["ReactionKineticsModels`"]
Unprotect["ReactionKineticsModels`*"];
ClearAll["ReactionKineticsModels`*"];
ClearAll["ReactionKineticsModels`Private`*"];


Reactions::usage = "Reactions is a list of some of the widespread reactions in reaction kinetics.";
Reactions = Association[
{
	"Autocatalator" -> {"A"->"X","X"->"Y","X"+2"Y"->3"Y","Y"->"P"},
	"Belousov-Zhabotinsky"-> {
            "X"+"Y"+"H"->2"V",
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
            "V"->"0"
        },
	"Briggs-Rauscher" -> {
            "IO3 -"+2"H2O2"+"CH2(CO2H)2"+"H +"->"ICH(CO2H)2"+2"O2"+3"H2O",
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
            "I2"+"CH2(CO2H)2"->"ICH(CO2H)2"+"H +"+"I -"
        },
	"Brusselator" -> {"X"\[LeftArrow]"A","B"+"X"->"Y"+"D","Y"+2"X"->3"X","X"->"E"},
	"Chapman cycle" -> {"O2"->"O"+"O", "O2"+"O"+"M"->"O3"+"M", "O3"->"O2"+"O", "O"+"O3"->"O2"+"O2"},
	"Clarke" -> {
            "A"+"T"->"M"+"S",
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
            "P"+"T"->"Q"+"Y"+"H"+"S"
        },
	"Colquhoun" -> {"AR"\[Equilibrium]"A"+"R"\[Equilibrium]"RA","A"+"AR"\[Equilibrium]"ARA"\[Equilibrium]"A"+"RA"},
	"Consecutive" -> {"A"->"B"->"C"},
	"Decomposition" -> {"C" -> "A"+"B"},
	"DeYoung-Keizer" -> {
            "S(000)"\[Equilibrium]"S(001)"\[Equilibrium]"S(100)"\[Equilibrium]"S(010)","S(100)"\[Equilibrium]"S(101)"\[Equilibrium]"S(110)","S(001)"\[Equilibrium]"S(101)"\[Equilibrium]"S(011)",
            "S(011)"\[Equilibrium]"S(010)"\[Equilibrium]"S(111)","S(110)"\[Equilibrium]"S(010)"\[Equilibrium]"S(111)","S(111)"\[Equilibrium]"S(101)"
        },
	"Dimerization" -> {"A"+"B" -> "C"},
	"Edelstein" -> {"X"\[DoubleLeftRightArrow]2"X","X"+"Y"\[DoubleLeftRightArrow]"Z"\[DoubleLeftRightArrow]"Y"},
	"Eigen" -> {"X"->0,"X"->2"X"},
	"\[CapitalEAcute]rdi-Ropolyi" -> {"R"+4"T"\[Equilibrium]"T4R1"\[Equilibrium]"T4R2"\[Equilibrium]"T4R3"},
	"Explodator" -> {"A"+"X"->(1+Global`aExplodator)"X","X"+"Y"->"Z","Z"->(1+Global`bExplodator)"Y","Y"->"P"},
	"Horn-Jackson" -> {3"X"->"X"+2"Y"->3"Y"->2"X"+"Y"->3"X"},
	"Feinberg-Horn I" -> {2"J" -> "G", 2"J" \[LeftRightArrow] "H", "G"->"H", "A"+"B" -> "G", "A"+"B" \[LeftRightArrow] "C", "C" \[RightArrow] "D"+"E", "D"+"E" \[LeftRightArrow] "F"},
	"Feinberg-Horn II" -> {"X"->"Y","X"+"Z"->"T"->"Y"+"U"->"X"+"Z"},
	"FKN mechanism" -> {
            "HBrO2"+"Br -"->2"HOBr","Br -"+"BrO3 -"->"HBrO2"+"HOBr",2"HBrO2"->"BrO3 -"+"HOBr","HOBr"+"Br -"->"Br2"+"H2O","BrO3 -"+"HBrO2"->2"BrO2"+"H2O",
            "HBrO2"+"H2O"->"BrO3 -","Br2"+"MA"->"BrMA"+"Br -","MA"+2"H2O"->C+2"CO2","BrMA"+2"H2O"->"Br -"+"HCO2H"+2"CO2"
        },
	"Frank" -> {"A" -> "R", "A" -> "S", "A"+"R" -> 2"R", "A"+"S" -> 2"S"},
	"Goldbeter" -> {"M"\[Equilibrium]"0"\[Equilibrium]"C","0"\[Equilibrium]"X","C"\[RightArrow]"0"},
	"Hudson-R\[ODoubleDot]ssler" -> {"P"->"A","Q"->"B","A"+"B"->2"B","B"->"R","A"\[LeftRightArrow]"C"},
	"Huxel" -> {"X"\[Equilibrium]2"X", "X"+"Y"->2"Y","Y"->"0"},
	"Hyver" -> {
            "K1"->"A1"->"A2"->"A3"->"A4"->"A5"->"A6"->"A7", "A7"+"B"->"C"+"E",
            "K"->"B"->"0", "A1"+"E"->"0", "A1"+"C"->"0", "\[Mu]"->"D", "D"+"C"->"0", "D"+"E"->"0"
        },
	"Inflow" -> {"0"->"X"},
	"Ivanova" -> {"X"+"Y"->2"Y","Y"+"Z"->2"Z","X"+"Z"->2"X"},
	"Kliemann" -> {"X"\[Equilibrium]"Y",2"X"\[Equilibrium]"X"+"Y"\[Equilibrium]2"Y"},
	"Leonard-Reichl" -> {2"X"\[LeftArrowRightArrow]"Y"},
	"Lotka" -> {"A"->"X","X"+"Y"->2"Y","Y"->"P"},
	"Lotka-Volterra" -> {"A"+"X"->2"X","X"+"Y"->2"Y","B"\[LeftArrow]"Y"},
	"Michaelis-Menten" -> {"E"+"S"\[Equilibrium]"C"->"E"+"P"},
	"Mole" -> {"X" + "Y" \[RightArrow] 2 "X" + 2 "Y", "X" \[LeftRightArrow] "0" \[LeftRightArrow] "Y"},
	"Ogg" -> {
            "\!\(\*SubscriptBox[\(N\), \(2\)]\)\!\(\*SubscriptBox[\(O\), \(5\)]\)" \[LeftRightArrow] "\!\(\*SubscriptBox[\(NO\), \(2\)]\)"+"\!\(\*SubscriptBox[\(NO\), \(3\)]\)",
            "\!\(\*SubscriptBox[\(NO\), \(2\)]\)"+"\!\(\*SubscriptBox[\(NO\), \(3\)]\)"->"\!\(\*SubscriptBox[\(NO\), \(2\)]\)"+"NO"+"\!\(\*SubscriptBox[\(O\), \(2\)]\)",
            "\!\(\*SubscriptBox[\(NO\), \(3\)]\)"+"NO"->2"\!\(\*SubscriptBox[\(NO\),
            \(2\)]\)"
        },
	"Oregonator" -> {"A"+"Y"->"X","X"+"Y"->"P","B"+"X"->2"X"+"Z",2"X"->"Q","Z"->Global`fOregonator*"Y"},
	"Outflow" -> {"X"->"0"},
	"Petri" -> {
            "C"+"\!\(\*SubscriptBox[\(O\), \(2\)]\)" -> "\!\(\*SubscriptBox[\(CO\), \(2\)]\)",
            "\!\(\*SubscriptBox[\(CO\), \(2\)]\)"+"NaOH" -> "\!\(\*SubscriptBox[\(NaHCO\), \(3\)]\)",
            "\!\(\*SubscriptBox[\(NaHCO\), \(3\)]\)"+"HCl" -> "\!\(\*SubscriptBox[\(H\), \(2\)]\)O"+"NaCl"+"\!\(\*SubscriptBox[\(CO\), \(2\)]\)"
        },
	"Robertson" -> {"A"->"B", 2"B"->"C"+"B", "B"+"C"->"A"+"C"},
	"Schl\[ODoubleDot]gl" -> {0\[Equilibrium]"X",2"X"\[Equilibrium]3"X"},
	"Simple birth-death" -> {"X"->2"X", "X"->"A"},
	"Triangle" -> {"A"->"B"\[RightArrow]"C"->"A"},
	"Tur\[AAcute]nyi-Gy\[ODoubleDot]rgyi-Field" -> {"X"+"Y" -> 2"P", "Y"+"A" -> "X"+"P", 2"X" -> "P"+"A", "X"+"A" -> 2"X"+2"Z", "X"+"Z" -> 0.5"X"+"A", "Z"+"M"->"Y"-"Z"},
	"Vaiman" -> {"O2" + 2 "K" -> 2 "KO", "H2S" + "K" -> "KH2S", "O2" + 2 "KH2S" -> 2 "H2O" + "K" + "KS2", "H2S" + "KO" -> "H2O" + "KS", "KO" + "KH2S" -> "H2O" + "KS" + "K", 2 "KS" -> "K" + "KS2", "KS2" -> "S2" + "K"},
	"Verhulst" -> {"X"\[Equilibrium]2"X"},
	"Volpert" -> {"X"+"Y"->"T","Y"+"Z"->"U"},
	"Wegscheider" -> {"A"\[Equilibrium]"B",2"A"\[Equilibrium]"A"+"B"},
	"Willamowski-R\[ODoubleDot]ssler" -> {"A"+"X"\[ReverseEquilibrium]2"X","X"+"Y"\[ReverseEquilibrium]2"Y","Q"+"Y"\[ReverseEquilibrium]"B","X"+"Z"\[ReverseEquilibrium]"C","P"+"Z"\[ReverseEquilibrium]2"Z"}
}];


Models::usage = "Models contains a list of names of some of the widespread reactions in reaction kinetics.";
Models = Sort[Keys[Reactions]];


BioReactions::usage = "BioReactions is a list of some of the widespread biochemical reactions.";
BioReactions=Association[{
        "Calvin Cycle (Photosynthesis)" -> {
            "DG3P"->"GP","DF1,6BP"->"GP"+"DG3P","DF1,6BP"+"H2O"->"DF6P"+"Pi","S7P"+"DG3P"->"DR5P"+"DX5P","DR5P"->"DRL5P","DRL5P"->"DX5P",
            "ATP"+"DRL5P"->"ADP"+"DRL1,5BP","3PDG"+2"H+"->"DRL1,5BP"+"CO2"+"H2O" ,"ATP"+"3PDG"->"ADP"+"3PDGP","DG3P"+"Pi"+"NADP+"->"3PDGP"+ "NADPH" + "H+"
        } /. {
            "DG3P"->Tooltip["DG3P","D-glyceraldehyde 3-phosphate"],
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
            "NADPH"->Tooltip["NADPH","NADPH"](*???*)
        },

	"Glycolysis" -> {
            "Glc"+"ATP"-> "G6P"+"ADP"+"H+","G6P"\[Equilibrium]"F6P","F6P"+"ATP"->"ADP"+"F1,6BP","F1,6BP"\[RightArrowLeftArrow]"GADP"+"DHAP","DHAP"\[Equilibrium]"GADP",
            "GADP"+"NAD+"+"Pi"-> "NADH"+"H+"+"1,3BPG","1,3BPG"->"GADP","1,3BPG"+"ADP"->"ATP"+"3PG","3PG"\[ReverseEquilibrium]"2PG",
            "2PG"->"PEP"+"H2O","PEP"->"2PG","PEP"+"ADP"-> "ATP"+"Pyr"
        } /. {
            "Glc"->Tooltip["Glc","D-glucose"],
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
            "NADH"->Tooltip["NADH","NADH"](*???*)
        },

	"Glyoxylate Cycle" -> {
            "SM"+"NAD+"->"Oxc" + "NADH" + "H+",
            "ACoA"+"H2O"+"Oxc"->"C"+"HS-CoA"+"H+",
            "ACoA"+"H2O"+"G"->"SM"+"HS-CoA"+"H+",
            "IsoC"->"S"+"G",
            "C"->"cA"+"H2O",
            "cA"+"H2O"->"IsoC"
        } /. {
            "Oxc"-> Tooltip["Oxc","oxalocetate"],
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
            "NADH"->Tooltip["NADH","NADH"](*???*)
        },

	"SzentGy\[ODoubleDot]rgyi-Krebs Cycle" -> {
            "Oc"+"ACoA"+"H2O"->"C"+"CSH",
            "C"->"cA"+"H2O"->"IsoC",
            "IsoC"+"NAD+"-> "Os"+"NADH"+"H+",
            "Os"-> "\[Alpha]K"+"CO2",
            "\[Alpha]K"+"NAD+"+"CSH"->"SCoA"+"NADH"+"H+"+"CO2",
            "SCoA"+"GDP"+"Pi"->"S"+"CSH"+"GTP",
            "S"+"Q"->"F"+"QH2",
            "F"+"H2O"->"LM",
            "LM"+"NAD+"->"Oc"+"NADH"+"H+"
        } /. {
            "Oc"->Tooltip["Oc","oxalocetate"],
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
            "NADH"->Tooltip["NADH","NADH"](*???*)
        }
    }
];


BioModels::usage = "BioModels contains a list of names of some of the widespread biochemical reactions.";
BioModels = Sort[Keys[BioReactions]];


GetReaction::usage = "GetReaction[model] returns a reaction if model is one of the built-in models. \
GetReaction[\"Models\"] returns all the available built-in models. See also Models.";
Begin["`Private`"];
ModelsSumPrivate = Join[BioModels, Models];
ReactionsSumPrivate = Join[BioReactions, Reactions];


GetReaction::badarg = "Illegal argument of function GetReaction.";
GetReaction::nvmod = "Argument '`1`' is not a valid built-in (bio)model. \
See GetReaction[\"Models\"] or Models.";

SyntaxInformation[GetReaction]={"ArgumentsPattern"->{__}};

GetReaction["Models"] := ModelsSumPrivate; (*Models*)
GetReaction[x_?StringQ] :=
	If[{Print["asdf"],
		KeyExistsQ[ReactionsSumPrivate,x]}[[2]], (*Models*)
		    ReactionsSumPrivate[x],
			Message[GetReaction::"nvmod",x];
			$Failed
	  ];

SetAttributes[GetReaction,Listable]
GetReaction[x__?StringQ] := GetReaction[{x}];
GetReaction[___] := (Message[GetReaction::"badarg"]; $Failed)


End[]
EndPackage[];
Protect["ReactionKineticsModels`*"];
