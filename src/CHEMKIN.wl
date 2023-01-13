(* ::Package:: *)

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
CHEMKINExport[___][___] := (Message[CHEMKINExport::"badarg"]; $Failed)
