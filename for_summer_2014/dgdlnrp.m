(* ::Package:: *)

BeginPackage["dgdlnrp`"]

d\[Gamma]dlnrp::usage = "Returns number per yr as a function of";
\[Gamma]listintegrate;
\[Gamma]direct;
Doratecalcsingle;

Begin["`Private`"]




<<dgdlnrp_helpers;


(* ::Input:: *)
(*G=6.67259 10^-8;*)
(*msun=1.989 10^33;*)
(*rsun=6.9599 10^10;*)
(*km=10^5;*)
(*yr=365 24 3600;*)


yr=365 24 3600;


d\[Gamma]dlnrp[u_,GD_List,fcns_List,tolin_,Sersic_,Emin_:.01,Emax_:100]:=Module[{MBHomsun,MBHnorm,rborT,tdynb,qmin,prefactor,sumtop,q,lnRlcm,f,counter=0,ansratio=1,ansnew,ans=0,ansratlim=tolin},


	Print["


__________________________	"];
	Print["u = ", u];

	Print["Remember to set Sersic to True or False"];
(* If Sersic is just a number, "Sersic" will be indeterminate, so the line below will evaluate the last expression *)
	If[ Sersic, {MBHomsun,MBHnorm,rborT,tdynb}=GD[[{3,4,6,8}]],
	{MBHomsun,MBHnorm,rborT,tdynb}=GD[[{5,6,7,9}]], {MBHomsun,MBHnorm,rborT,tdynb}=GD[[{5,6,7,9}]] ];

	prefactor = 8 \[Pi]^2 MBHomsun rborT^-1 tdynb^-1;
	Print["prefactor (per yr) = ", prefactor yr];
	Print["MBHomsun = ", MBHomsun, ", rborT = ", rborT, ", tdynb(yr) = ", tdynb/yr];

If[u==0, ans=0,	

	Off[InverseFunction::"ifun"];
	Off[Solve::"tdep"];
	Off[NIntegrate::"slwcon"];
	Off[NIntegrate::"ncvb"];

	{q,lnRlcm,f}=fcns[[{1,2,3}]];

(*	Emin= 0.01; *)
(*	Emin=0.001; *)
(*	Emax= MBHnorm rborT; *)
(*	Emax=10^3; *)
(*Print["counter, ansratio ", counter, ansratio];*)
	If[tolin == Null, ansratlim=0.01];
	
	While[(counter==0),
(*	While[(counter<=4 && ansratio>ansratlim), *)
		
		Print["Emin ", Emin, " Emax ", Emax];

		qmin=q[Emax];
		Print["qmin = ",qmin,", qmax = ", q[Emin]];
		Print["sumtop = ",2qmin^(-1/2)];

(* qmin=300;*)
(* sumtop=Max[200,2qmin^(-1/2)]; *)
(* Print[sumtop];
	Print[alphamfastfull[7 10^6]];
ans = Sum[E^(-(alphamfastfull[m]^2)/4) ,{m,1,9.9 10^5}]; *)
	(* ans= Sum[k^2,{k,1,20}]; *)
 	(*  ans = Sum[(E^(-(alphamfastfull[m]^2)/4) Besselfin[u,m])/mpiecegood[m],{m,1,20}];   *)
	 (* ans = Sum[E^(-(alphamfastfull[m]^2)/4),{m,1,20}]; *)
		ansnew= prefactor u^2 NIntegrate[f[EE]/(1+q[EE]^-1 xigood[q[EE]]lnRlcm[EE]) 
			(1-2Sum[(E^(-(alphamfastfull[m]^2 q[EE])/4) Besselfin[u,m])/mpiecegood[m],
			{m,1,Max[200,2qmin^(-1/2)]}]),{EE,Emin,Emax},AccuracyGoal->4];
		Print["temp ans: rate per lnrp per yr = ",ansnew yr];
		ansratio=Abs[ans-ansnew]/ansnew; (* first time through loop, ansratio = 1, so we run it again *)
		Print["ansratio = ", ansratio];
	(*	Emin = Emin / 2;
		Emax = Emax * 2; *)
		counter++;
		ans = ansnew;
		]


	If[counter>4,Print["**** Warning: did not converge in Emax! counter = ", counter]];
]
	Print["

rate per lnrp per yr = ",ans yr];

	On[InverseFunction::"ifun"];
	On[Solve::"tdep"];
	On[NIntegrate::"slwcon"];
	On[NIntegrate::"ncvb"];


	ans
]


\[Gamma]direct[GD_List,fcns_List,tolin_,Sersic_]:=Module[{MBHomsun,MBHnorm,rborT,tdynb,Emin=0.1,Emax=10,xicomp,intfcn,intfcncomp,qmin,prefactor,sumtop,q,lnRlcm,f,theintegral,theintegralcomp,counter=0,ansratio=1,ansnew,ans=0,ansratlim=tolin},

	Off[InverseFunction::"ifun"];
	Off[Solve::"tdep"];
	Off[NIntegrate::"slwcon"];
	Off[NIntegrate::"ncvb"];


	Print["Remember to set Sersic to True or False"];
	If[ Sersic, {MBHomsun,MBHnorm,rborT,tdynb}=GD[[{3,4,6,8}]],
	{MBHomsun,MBHnorm,rborT,tdynb}=GD[[{5,6,7,9}]], {MBHomsun,MBHnorm,rborT,tdynb}=GD[[{5,6,7,9}]] ];
	
	Print[MBHomsun];
	
	prefactor = 8 \[Pi]^2 MBHomsun rborT^-1 tdynb^-1;
	Print["prefactor (per yr) = ", prefactor yr];

	{q,lnRlcm,f}=fcns[[{1,2,3}]];

(*	If[Emin==0,Emin= 0.001];
(*	If[Emax==Null,Emax=MBHnorm rborT]; *)
	If[Emax==0,Emax=10^4]; *)


(*	Print[Emin];
	Print[Emax]; *)
	Print[q[1]];

	intfcn=f[#] q[#] / (q[#]/xigood[q[#]] + lnRlcm[#]) &;
	Print["compare plot:  at E = 42, we have ", 42 intfcn[42]];


	If[tolin == Null, ansratlim=0.01];
	qmin=q[Emax];
	
	While[(counter<=4 && ansratio>ansratlim),
		
		Print["qmin = ",q[Emax],", qmax = ", q[Emin]];
	
	
(* 	theintegral = NIntegrate[f[E] q[E] / (q[E]/xigood[q[E]] + lnRlcm[E]),{E,Emin,Emax}]; *)
		theintegral = NIntegrate[intfcn[EE],{EE,Emin,Emax}]; (* don't use E as int variable -- Mathematica thinks of it as the number e!! *)

		Print["tempans: the integral = ", theintegral];
		ansnew= prefactor theintegral;
		Print["tempans: rate per yr = ",ansnew yr];
		Print["tempans: LOG10 rate per yr = ", Log10[ansnew yr]];

		ansratio=Abs[ans-ansnew]/ansnew; (* first time through loop, ansratio = 1, so we run it again *)

		Print["ansratio = ", ansratio];
		Print[ans, " ", ansnew, " ", ansratio];

(*		qmin = qmin / 100;  *)
		Emin = Emin / 10;
		Emax = Emax * 10;
		counter++;
		ans = ansnew;
		]

	If[counter>4,Print["**** Warning: did not converge in Emax! counter = ", counter]];
	
	Print["
		"];
	Print["FINAL ANS: the integral = ", theintegral];
	Print["FINAL ANS: rate per yr = ",ans yr];
	Print["FINAL ANS: LOG10 rate per yr = ", Log10[ans yr]];
	Print["
		"];

	xicomp=Function[Piecewise[{
		{#,#>1},
		{0.186#+0.824Sqrt[#],# <=1}
		}]];

	intfcncomp=f[#] q[#] / (xicomp[q[#]] + lnRlcm[#]) &;
	theintegralcomp = NIntegrate[intfcncomp[EE],{EE,Emin,Emax}]; (* don't use E as int variable -- Mathematica thinks of it as the number e!! *)
	Print["with approx xi: the integral = ", theintegralcomp];
	Print["with approx xi:  rate per yr = ",prefactor theintegralcomp yr];
	Print["with approx xi:  LOG10 rate per yr = ", Log10[prefactor theintegralcomp yr]];


	Print["intfcn(emin), intfcn(emax), ", intfcn[Emin], ", ", intfcn[Emax]];
	Print["intfcncomp(emin), intfcncomp(emax), ", intfcncomp[Emin], ", ", intfcncomp[Emax]];
	Print["intfcn(emin*2), intfcn(emax/2), ", intfcncomp[Emin*2], ", ", intfcncomp[Emax/2]];
	Print["intfcncomp(emin*2), intfcncomp(emax/2), ", intfcncomp[Emin*2], ", ", intfcncomp[Emax/2]];


(*	Show[Plot[{intfcncomp[EE]},{EE,Emin,Emax}]]; *)


	On[InverseFunction::"ifun"];
	On[Solve::"tdep"];
	On[NIntegrate::"slwcon"];
	On[NIntegrate::"ncvb"];

(*	{ans, intfcn, intfcncomp} *)
	ans
]


\[Gamma]listintegrate[rportd\[Gamma]dlnrptable_List]:=Module[{fcntable,interp,bottom,dgdlnrpcont,ans},
	
	If[rportd\[Gamma]dlnrptable[[1]]=={0,0},
		fcntable=rportd\[Gamma]dlnrptable[[2;;Length[rportd\[Gamma]dlnrptable]]],
		fcntable=rportd\[Gamma]dlnrptable
	];
		(* interp = Interpolation[Log10[rportd\[Gamma]dlnrptable[[2;;Length[rportd\[Gamma]dlnrptable]]]]], *)
	interp = Interpolation[Log10[fcntable]];
	bottom = fcntable[[1,1]];
	dgdlnrpcont = Function[Piecewise[{
		{10^interp[Log10[#]],#>=bottom},
		{fcntable[[1,2]] #/bottom,#<bottom}}]];
	ans = NIntegrate[dgdlnrpcont[rport]/rport,{rport,0,1}];
(* The problem is that I can't calculate d\[Gamma]dlnrp at rp = rT, and so I have to extrapolate that point in the integral, 
and it can affect the answer somewhat because it's a large number (even if only over a small range.) *)

	ans
]



Doratecalcsingle[GD_List,fcns_List,stepsizein_,tol_,Sersic_]:=Module[{filestem,stepsize=stepsizein,utable,d\[Gamma]dlnrptable,result,intresult,directresult},

filestem="~/Documents/research/tidal_dynamics/dgdlnrp_tables_Faber/" ;
(* Get[filestem <>ToString[GDall[[i,1]]]<>".fcns"]; *)


If[stepsizein==Null, stepsize=0.01];

(*
utable95 = Table[u,{u,0,0.95,0.01}]; (* u = Sqrt[rp/rT] *) 
d\[Gamma]dlnrptable95 = Table[d\[Gamma]dlnrp[u,10^3,GD,fcns],{u,utable95}];
(* Emax can be lower and we've already converged *)

utable5 = {0.96,0.97,0.98,0.99}; (* u = Sqrt[rp/rT] *) 
d\[Gamma]dlnrptable5 = Table[d\[Gamma]dlnrp[u,10^5,GD,fcns],{u,utable5}];
(* close to rT, it's better converged with higher Emax.  see my convergence testing further down this notebook. *)


utable=Join[utable95,utable5];
d\[Gamma]dlnrptable=Join[d\[Gamma]dlnrptable95,d\[Gamma]dlnrptable5];
*)

utable = Table[u,{u,0,1-stepsize,stepsize}]; (* u = Sqrt[rp/rT] *) 
Print["utable", utable];
d\[Gamma]dlnrptable = Table[d\[Gamma]dlnrp[u,GD,fcns,tol,Sersic],{u,utable}];

result = Partition[Riffle[utable^2,d\[Gamma]dlnrptable  yr],2];
(* rport = u^2 *)


Print[ListLogLogPlot[result,Mesh->All,Joined->True,AxesLabel->{"\!\(\*FractionBox[SubscriptBox[\(r\), \(p\)], SubscriptBox[\(r\), \(T\)]]\)","\!\(\*FractionBox[\(d\[Gamma]\), \(d\\\ ln\\\ \*SubscriptBox[\(r\), \(p\)]\)]\)(\!\(\*SuperscriptBox[\(yr\), \(-1\)]\))"}]];

intresult=\[Gamma]listintegrate[result];
Print["integrating result: total rate per yr = ",intresult];

directresult=\[Gamma]direct[GD,fcns,tol,Sersic];
Print["direct result: total rate per yr = ",directresult  yr];

(* This DOES NOT WORK!! IT DOESN'T SAVE THE LOCAL VARIABLE RESULT; IT TRIES TO SAVE THE GLOBAL VARIABLE RESULT!!
Save[filestem <> ToString[GD[[1]]]<>".dgdlnrp","result"]; (* note the quotes around result *)
*)


result

]


(* nupts = 10;
qmin = 10^-3;
(* qmin[u_]:=Piecewise[{{10^-3,u<=0.9},{10^-4,u>0.9}}]; (* don't go above 0.99 *) *)
qmax=100;
preamble={{MBH,0},{sigmah,0}, {\[Gamma],0},{qmax,0}};(* For putting in file header *)
(* utable = Table[u2,{u2,0,1,nupts^-1}]; (* u2 = u^2 = rp/rT *) *) *)


(*
utable = Table[u,{u,0,1,nupts^-1}]; (* u = Sqrt[rp/rT] *) 
d\[Gamma]dlnrptable = Table[d\[Gamma]dlnrp[u,MBH,sigmah,\[Gamma],qmin,qmax],{u,utable}]
result = Partition[Riffle[utable^2,d\[Gamma]dlnrptable yr],2];
totalvec = Join[preamble, result];(* To make exporting easier *)
filename="~/Documents/research/tidal_dynamics/dgdlnrp_tables/dgdlnrp" <>filenum<>".dat"; *)


(* ***************************** End ***************************** *)

End[] (* `Private` context *)

EndPackage[](* dgdlnrp` *)
