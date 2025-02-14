(*
::Package::

 :Context: PlaneCurvePlot`

 :Title: PlaneCurvePlot

 :Author: Xah Lee

 : Copyright 1994-2024 by Xah Lee.

 :Summary:

This package exports the function PlaneCurvePlot. PlaneCurvePlot is like ParametricPlot with extra options to plot tangent, secant, normal lines, parallel curves, caustics, evolutes, radial, concoid, and inversion of a curve. These features are useful for studing plane curves.

 :Keywords: graphics, geometry, calculus, curve, plot, tangent, secant, normal, parallel, caustics, involute, evolute, radial, concoid, inversion, pedal

:Homepage: http://xahlee.info/M/plane_curve_plot.html

 :Package Version: 4.0.20250214021935

 :History:
Version 4. 2024-03-01. all rewritten for wolfram lang 13.
Version 3.4. 199803xx. Converted the package to 3.0 format. Minor editing here and there.
Version 3.33 199609xx. Added DeleteCases in the main body of PlaneCurvePlot so the output grahpics does not contain empty lists or lists containing graphics directives only.
Also fixed a precaution so that single item fed to Options can be
specified whether it has the head List. e.g. CurveColorFunction -> {Hue[#]&} or CurveColorFunction -> (Hue[#]&) are both okay.
Version 3.32 199510xx. Loading package do not switch off messages.
Version 3.3 199508xx. Added option to plot inversion and involute.
Version 3.2 19950520. Added option to plot Conchoid, Pedal, Evolute,
Radial, and Osculating circle.
Version 3.1 19950500. Added two new options Caustic and Parallel.
Version 3.0 19950400. Complete redesign and rewrite.
Version 2.0 19940700. First published on CompuServe.

 :Discussion:

Load the package.
Type
?PlaneCurvePlot
for documentation.
all symbols has inline documentation.

Xah Lee
http://xahlee.org/

*)

BeginPackage["PlaneCurvePlot`"];

Clear[CausticLineLength,
CausticOrigin,
CausticLineStyle,
ConchoidSetting,
CurveColorFunction,
CurveStyle,
DotStyle,
InversionCircle,
LineStyle,
NormalLineLength,
NormalLineStyle,
ParallelCurveOffset,
PedalPoint,
PlotCurve,
PlotDot,
PlotOsculatingCircle,
PlotInvolute,
PlotEvolute,
PlotParameterValue,
PlotRadial,
SecantConnection,
SecantLineStyle,
TangentLineLength,
TangentLineStyle
];

CausticLineLength::usage=
"CausticLineLength is an options for PlaneCurvePlot.
Specifies the length of catacaustic rays (lines) to be plotted.
That is, lines is drawn starting from CausticOrigin -> {0,0} (the default) toward the curve, and reflected off from it.
If the value is a negative number, draw diacaustic rays.
The default CausticLineLength -> 0, plot no rays.
Example:
(* Parabola and its optical property *)
PlaneCurvePlot[ {(1 -2 x + x^2)/4+1, x},{x, -2, 4, .4}, CausticLineLength -> 5, CausticOrigin -> {2,1}, PlotDot -> False]
";

CausticOrigin::usage=
"see ?CausticLineLength";

CausticLineStyle::usage=
"CausticLineStyle is an option for PlaneCurvePlot.
List of graphics directives for CausticLineLength.
Elements are used cyclically.
Example:
PlaneCurvePlot[{2 Cos[t], 1 Sin[t]}, {t, 0, 2 Pi}, Range[0.3, 2.5, .18], CausticLineLength -> 3.9, CausticOrigin -> {1.737, 0}, PlotDot -> False, CausticLineStyle -> {Yellow, Green}, Epilog -> {Red, PointSize[.02], Point[{-1.737, 0}]}]
";

ConchoidSetting::usage=
"ConchoidSetting is an option for PlaneCurvePlot.
ConchoidSetting -> {{x,y},{n1,n2,...}} specifies the plotting of conchoid with respect to the point {x,y}, and offsets {n1, n2,...}.
Default is ConchoidSetting -> {{0,0},{0}}.
Use the option DotStyle and LineStlye to control rendering.
(*Conchoid is a method to generate new curves.
Given any curve C and a point O, a conchoid of that curve with respect to O and offset k is defined.
Let P be any point on the curve.
Draw a line passing O and P.
Mark off points p1 and p2 on the line with distance +k,-k from P.
The locus of p1 and p2 is the conchoid of the curve with respect to O and offset k.*)
Example:
(* Conchoid of Nicomedes *)
PlaneCurvePlot[{x,1},{x,-3,3,.2}, ConchoidSetting -> {{0,0},{-2,2,1}}, DotStyle -> {Red,Blue,Magenta}, PlotDot -> False]
";

CurveColorFunction::usage=
"CurveColorFunction is an option for PlaneCurvePlot that specifies how the curve should be colored varying with the parameter t.
CurveColorFunction -> myColorFunc, where myColorFunc must be a function that can be mapped into a real number list between 0 and 1 and return a color directive.
Default is CurveColorFunction -> Function[Hue[.65,#,.8]].
Example:
(* Hypotrochoid *)
PlaneCurvePlot[ {(5 Cos[t])/8 + 0.1 Cos[(5 t)/3], (5 Sin[t])/8 - 0.1 Sin[(5 t)/3]}, {t, 0, 3 2 Pi}, CurveColorFunction -> Hue, CurveStyle -> {Thickness[.015]}, PlotDot -> False]
";

CurveStyle::usage=
"CurveStyle is an option for PlaneCurvePlot.
List of graphics directives for curve style except color.
To specify color, use CurveColorFunction.
Default is CurveStyle -> {Thickness[0.005]}.";

DotStyle::usage=
"DotStyle is an option for PlaneCurvePlot.
List of graphics directives for dots on the curve.
Elements are used cyclically.
Example:
PlaneCurvePlot[{Sin[t 2 Pi/12],Cos[t 2 Pi/12]}, {t, 0, 11}, Range[0,12], PlotParameterValue -> True, DotStyle -> { {PointSize[.04],Red}, {PointSize[.02],Hue[.4,.8,.7]}, {PointSize[.02],Hue[.4,.8,.7]} }, Ticks -> None]
";

InversionCircle::usage=
"InversionCircle is an option for PlaneCurvePlot.
InversionCircle -> {{x,y},r} plots the inversion of the curve with respect to the circle with center {x,y} and radius r.
Default is InversionCircle -> None.
Use DotStyle to control the rendering of dots.
Example:
(* Rose and Epi spiral *)
PlaneCurvePlot[{Cos[3 t] Cos[t], Cos[3 t] Sin[t]}, {t, 0, Pi, Pi/120}, InversionCircle -> {{0, 0}, 1}, PlotDot -> False, PlotRange -> {{-4, 4}, {-4, 4}}]
";

LineStyle::usage=
"LineStyle is an option for PlaneCurvePlot.
List of graphics directives for lines or circles on the curve.
Elements are used cyclically.
Example:
PlaneCurvePlot[{t - Sin[t], 1 - Cos[t]}, {t, 0, 2 2 Pi, (2 2 Pi)/60}, PlotRadial -> True, LineStyle -> {Red, {Yellow, Thickness[.005]}}, AspectRatio -> Automatic]
";

NormalLineLength::usage=
"NormalLineLength is an option for PlaneCurvePlot.
Specifies the length of normal lines to be plotted.
The default NormalLineLength -> 0, plot no normal lines.
Example:
(* Parabola and its Evolute *)
PlaneCurvePlot[ {t, t^2}, {t, -3, 3, .1}, NormalLineLength -> 15, NormalLineStyle -> {Red}, PlotDot -> False, PlotRange -> {{-6,6},{-2,5}}]
";

NormalLineStyle::usage=
"NormalLineStyle is an option for PlaneCurvePlot.
List of graphics directives for NormalLineLength.
Elements are used cyclically.
Example:
PlaneCurvePlot[{t, t^2}, {t, -3, 3, .1}, NormalLineLength -> 15, NormalLineStyle -> {Black, White}, PlotDot -> False, PlotRange -> {{-6, 6}, {-2, 5}}, Background -> Gray]
";

ParallelCurveOffset::usage=
"ParallelCurveOffset is an option for PlaneCurvePlot.
Specifies the offset distance from the curve to its parallel curve.
The default ParallelCurveOffset -> {0}, plot no parallel curves.
You can specify more than one parallel curve by giving a list {offset1, offset2,...}.
Example:
(* Lemniscate of Bernoulli *)
PlaneCurvePlot[{Cos[t]/(1+Sin[t]^2), Sin[t] Cos[t]/(1+Sin[t]^2)}, {t, 0, 2 Pi,.05}, ParallelCurveOffset -> {.3,-.3}, PlotDot -> False, DotStyle -> {Red, Blue}]
";

PedalPoint::usage=
"PedalPoint is an option for PlaneCurvePlot.
PedalPoint -> {x,y} plots the pedal of a curve with respect to {x,y} together with tangent lines showing how pedals are generated.
Default is PedalPoint -> None.
Use LineStyle and DotStyle to control the rendering of lines and dots.
What is the Pedal of a curve?
Pedal is a way of defining a new curve from a given curve and a point.
Given a curve C (e.g. a circle), and a point M.
Draw a tangent at any point P on the curve.
Mark a point Q on the tangent line such that MQ is perpendicular to the tangent line.
The locus of Q is the pedal of the given curve with respect to the point M.
Example:
(* Parabola and Line *)
PlaneCurvePlot[{t,1/2 t^2},{t,-2,2,.2}, PedalPoint -> {0,1/2}, PlotDot -> False, DotStyle -> {{Red,PointSize[.02]}}, LineStyle -> {Red}, CurveStyle -> {Thickness[.008]}, PlotRange -> {{-1.2,1.2},{-.2,.8}}]
";

PlotCurve::usage=
"PlotCurve is an option for PlaneCurvePlot.
Example:
PlaneCurvePlot[{Cos[t], Sin[t]}, {t, 0, 6, 6/30}, PlotCurve -> False, PlotDot -> True]
";

PlotDot::usage=
"PlotDot is an option for PlaneCurvePlot.
Specifies whether the dots on the curve generated by the parameter t list should be plotted.
Example:
PlaneCurvePlot[{Cos[t], Sin[t]}, {t, 0, 6}, PlotDot -> False]
";

PlotOsculatingCircle::usage=
"PlotOsculatingCircle is an option for PlaneCurvePlot.
The default is PlotOsculatingCircle -> False.
Use LineStyle to control the rendering.
Example:
(* Equiangular spiral *)
PlaneCurvePlot[ {1/4 E^(t Cot[1.4]) Cos[t], 1/4 E^(t Cot[1.4]) Sin[t] }, {t, 0, 4 Pi, 4 Pi/40}, PlotOsculatingCircle -> True, LineStyle -> {Red}, PlotRadial -> True, PlotCurve -> False, PlotDot -> False]
";

PlotInvolute::usage=
"PlotInvolute is an option for PlaneCurvePlot.
PlotInvolute -> True plot the involute of the given curve.
Default is PlotInvolute -> False.
Use the option DotStyle and LineStyle to control rendering.
Example:
PlaneCurvePlot[{t,Sin[t]}, {t, 0, 2 Pi, 2 Pi/80}, PlotInvolute -> True, PlotDot -> False]
";

PlotEvolute::usage=
"PlotEvolute is an option for PlaneCurvePlot.
The default is PlotEvolute -> False.
Evolute of a curve is the locus of its center of curvature.
Example:
(* Tracktrix and Catenary *)
PlaneCurvePlot[{ Log[ Sec[t] + Tan[t]] -Sin[t], Cos[t]}, {t, -1.53, 1.53,.1}, PlotEvolute -> True, PlotRange -> {{-3,3},{0,3}}]
";

PlotParameterValue::usage=
"PlotParameterValue is an option for PlaneCurvePlot.
Specifies whether the the value of t should be superscripted in the plot.
The default is PlotParameterValue -> False.
Example:
PlaneCurvePlot[{t, Sin[t]}, {t, 0, 6, 1}, PlotParameterValue -> True]
";

PlotRadial::usage=
"PlotRadial is an option for PlaneCurvePlot.
The default is PlotRadial -> False.
Radial of a curve is the locus of vectors defined by P2-P1, where P1 is any point on the curve and P2 is the center of curvature at P1.
Example:
(* Astroid and Rose *)
PlaneCurvePlot[ {(3 Cos[t])/4 + Cos[3 t]/4, (3 Sin[t])/4 - Sin[3 t]/4}, {t,0, 2 Pi, 2 Pi/100}, PlotRadial -> True, PlotDot -> False, LineStyle -> {Red},DotStyle -> {{Hue[.9], PointSize[.01]}}, Axes -> False]
";

SecantConnection::usage=
"SecantConnection is an option for PlaneCurvePlot, specifying which dots on the curve to be connected to which.
SecantConnection -> {m,n} connect every (1+i*n)th point with (1+i*n+m)th point, where i starts from 0.
The default SecantConnection -> {1,1} do not plot any secant line.
Example:
(* Circle *)
PlaneCurvePlot[{Cos[t/5],Sin[t/5]}, {t, 0, 30,1}, SecantConnection -> {6,2}, SecantLineStyle -> {Hue[.4,.8,.7]}, PlotParameterValue -> True, DotStyle -> {{PointSize[.02],Hue[.4,.8,.7]},{PointSize[.02],Red}}, Ticks -> None]
";

SecantLineStyle::usage=
"SecantLineStyle is an option for PlaneCurvePlot.
Specifies how secant lines should be rendered.
Example:
PlaneCurvePlot[{t Cos[t], t Sin[t]}, {t, 0, 4 2 Pi, 8 Pi/125}, SecantConnection -> {10, 1}, PlotDot -> False, SecantLineStyle -> Table[Hue[i], {i, 0, 1, 1/125}]]
";

TangentLineLength::usage=
"TangentLineLength is an option for PlaneCurvePlot.
Specifies the length of tangent lines to be plotted.
The default TangentLineLength -> 0, plot no tangent lines.
Use TangentLineStyle to control how you want these lines to be rendered.
Example:
(* TangentCrawl *)
PlaneCurvePlot[ {t, Cos[t] Sin[t] +t}, {t, -5, 10 }, Range[3, 5 , 0.05], TangentLineLength -> 6, PlotDot -> False, PlotRange -> {{-6,11},{-6,11}}]
";

TangentLineStyle::usage=
"TangentLineStyle is an option for PlaneCurvePlot.
List of graphics directives for TangentLineLength.
Elements are used cyclically.
Example:
(* colorful tangents on epitrochoid *)
PlaneCurvePlot[ {Cos[t]+Cos[t/5] 4, Sin[t]+Sin[t/5] 4}, {t, 0, 9 Pi,.1}, TangentLineLength -> 6, PlotDot -> False, PlotCurve -> False, TangentLineStyle -> Table[{Hue[i]},{i,0,1,.05}]]
";

Clear[ PlaneCurvePlot ];

PlaneCurvePlot::usage = "PlaneCurvePlot[{x[t],y[t]},{t,tmin,tmax}] plot the parametric curve {x[t],y[t]} where the curve color vary with t.
PlaneCurvePlot[{x[t],y[t]},{t,tmin,tmax,dt}] show dots on curve with parameter values in increment of dt.
PlaneCurvePlot[{x[t],y[t]},{t,tmin,tmax}, {t1,t2,...,tn}] show dots with explicit parameter values t1, t2, etc.
Options allow one to plot tangents, normals, caustics, parallel curves, conchoid, pedal curves, inversion, osculating circles, evolute, radial, involute.";

Begin[ "Private`" ]

Clear[ normalLineGP ];

normalLineGP::usage="normalLineGP[tPts, tangentVecs, lnlen]
tPts is a list of points {x,y} representing the points of a parametric curve.
tangentVecs is a list of tangent vectors.
lnlen is a number representing a length.
normalLineGP return a list of Line[...].
These are lines normal to the curve at tPts with length lnlen.";

normalLineGP = Function[{tPts, tangentVecs, lnlen},
Function[{norms},  Line /@ Transpose[{tPts - norms, tPts + norms}] ][Map[ Function[ Function[{-Last[#], First[#]}][ Normalize[ # ] lnlen/2. ] ] , tangentVecs] ]
];

Clear[ tangentLineGP ];

tangentLineGP::usage="tangentLineGP[tPts, tangentVecs, lnlen]
tangentLineGP return a list of Line[...].
These are lines tangent to the curve at tPts with length lnlen.";

tangentLineGP = Function[{tPts, tangentVecs, lnlen},
Function[{tans},
Line /@ Transpose[{tPts - tans, tPts + tans}]]
[ Map[ Function[ Normalize[ # ] lnlen/2. ] , tangentVecs] ]
];

parallelDotGP::usage="parallelDotGP[tPts, tangentVecs, distList]
is similar to normalLineGP. distList is a list of number each representing a distance.
return a list of the form
{
{Point[{xa1,ya1}], Point[{xa2,ya2}],...},
{Point[{xb1,yb1}], Point[{xb2,yb2}], ...},
...
}
The ith list of points represents the curve parallel to the curve represented by tPts, with distance distList[[i]]";

parallelDotGP = Function[{cuvPts, velociVectors, distList},
With[{xnomals = Map[ Function[ Function[{-Last[#], First[#]}][Normalize[ # ]] ], velociVectors]},
Map[Point, Map[ Function[ cuvPts + xnomals*#1 ] , distList], {2}]
] ];

Clear[ reflDirVector ];

reflDirVector::usage="reflDirVector[rayV, wallV]
return a unit vector that is the ray vector rayV reflected off a wall vector wallV.
The ray and wall vectors do not have to be unit vectors, and if the wall is the zero vector, then zero vector is returned.";

reflDirVector = Function[{rayV, wallV},
If[Norm[ wallV ] == 0 , {0, 0},
Normalize[- rayV + (2*( wallV . rayV * wallV))/(wallV . wallV) ]]
] ;

Clear[ causticLineGP ];

causticLineGP::usage="causticLineGP[tPts, tangentVecs, rayOrigin, lnlen]
return a list of Line[...].
They represent rays from rayOrigin shooting towards the points on the parametric curve and reflected.
Each line length is lnlen, including the reflected portion.";

causticLineGP = Function[{tPts, tangentVecs, rayOrigin, lnlen},
With[ { xsign = Sign[lnlen], lnlen2 = Abs[lnlen]},
If[ Norm[ rayOrigin ] == 0,
MapThread[
If[Norm[#1] >= lnlen2,
Line[{{0,0}, Normalize[#1]*lnlen2}],
Line[{{0,0}, #1, #1 + reflDirVector[#1, #2]* (lnlen2 - Norm[#1])* xsign}]
] &, {tPts, tangentVecs}]
,
MapThread[
If[Norm[#1-rayOrigin] >= lnlen2,
Line[{rayOrigin , rayOrigin + Normalize[#1-rayOrigin]*lnlen2}],
Line[{rayOrigin , #1, #1 + reflDirVector[#1-rayOrigin, #2]* (lnlen2 - Norm[#1-rayOrigin])* xsign}]
] &, {tPts, tangentVecs}]
] ] ];

Clear[ secantLineGP ];

secantLineGP::usage="secantLineGP[ {{x1,y1},{x2,y2},..}, integer, integer]
return a list of the form
{Line[{{a,b},{c,d}}], Line[],...}
These lines connect points given by the first argument.
It connects every (1+i*n)th point with (1+i*n+m)th point, where i starts from 0 (until there are no more points to connect).";

secantLineGP = Function[{tPts, m, n}, Map[ Function[(Line[{First[#1], Last[#1]}] )], Partition[tPts, m + 1, n]]];

Clear[ conchoidLineAndDotGP ];

conchoidLineAndDotGP::usage="conchoidLineAndDotGP[point, tPts, distList]
return a list of the form
{
{Line[..],Line[..],..},
{{Point[..],Point[..],..},
..}
}
To see what are these, try
Show[Graphics[ conchoidLineAndDotGP[{0,0},{{1,1},{2,1}},{-.5,.5,1}] ],Axes->True]";

conchoidLineAndDotGP = Function[{pointP, tPts, distList},
With[{ dirVList = Map[ Function[Normalize[#1 - pointP]], N[tPts]] },
{Line /@ Transpose[{
tPts + Min[distList]*dirVList,
tPts + Max[distList]*dirVList}],
Map[Point, Function[tPts + dirVList*#1 ] /@ distList, {2}]}] ];

Clear[ getPedalPoint ];

getPedalPoint::usage="getPedalPoint[ {t1,t2}, {tp1,tp2}, {p1,p2}]
return a point of the form {q1,q2}
such that {q1,q2} == t {t1,t2} + {tp1,tp2} for some t,
and
Dot[ {q1,q2}-{p1,p2}, {t1,t2}] == 0.

Explanation: {t1,t2} represent a tangent vector of a curve at {tp1,tp2}.
{q1,q2} is a point on this tangent line such that {q1,q2} to {p1,p2} is perpendicular to it.";

getPedalPoint = Function[{cT, tp, p},
Function[If[ SameQ[# == {0, 0}, True ] , {}, tp + ((p - tp) . #) * #]][Normalize[cT]]
 ];

Clear[ pedalLineAndDotGP ];

pedalLineAndDotGP::usage="pedalLineAndDotGP[pointP, tPts, tangentVecs]
return a list of the form
{
 {Line[ { tPts[[1]], pedalpoint1, pointP} ],
  Line[ { tPts[[2]], pedalpoint3, pointP} ],
  Line[ { tPts[[3]], pedalpoint3, pointP} ],
  ...},
 {Point[pedalpoint1], Point[pedalpoint2], Point[pedalpoint3], ...}
}.
the three points inside Line[] makes a right angle, and (ti-pedalpointi) have the same direction as the ith part of tangentVecs.";

pedalLineAndDotGP = Function[{pointP, tPts, tangentVecs},
Function[{pedalPoints},
{MapThread[If[#2 === {}, {}, Line[{#1, #2, pointP}]] & , {tPts, pedalPoints}],
(If[#1 === {}, {}, Point[#1]] & ) /@ pedalPoints}
][ MapThread[getPedalPoint[#1, #2, pointP] & , {tangentVecs, tPts}]]
];

Clear[ inversionDotGP ];

inversionDotGP::usage="inversionDotGP[{c1,c2}, r, {p2, p2,...}]
return a list of Point[...].
These points are the inverse of points {{tx1,ty1},{tx2,ty2},...} with respect to the circle centered at {c1,c2} and radius r.";

inversionDotGP = Function[{center, radius, tPts},
If[ Norm[center] == 0 ,
Point /@
Map[ Function[ If[Norm[#1] == 0, Nothing, (radius^2*Normalize[#1])/ Norm[#1]]] , tPts] ,
Point /@
Map[ Function[ If[Norm[#1 - center] == 0, Nothing, (radius^2*Normalize[#1 - center])/ Norm[#1 - center] + center  ]] , tPts]
]
];

Clear[involuteLineAndDotGP];

involuteLineAndDotGP::usage="involuteLineAndDotGP[ fx, fy, tList, curvePoints, velociVectors]
return a list of the form
{ {Line[..],Line[..],...}, {Point[..],Point[..],...} }
fx and fy are the curve functions.
tList is a list of the cuvrve's parameters.
curvePoints is a List of {x,y}, numerical values, that are points on the curve.
velociVectors are the tangent vectors of the curve at curvePoints.

The list of graphic primitives returned represents the involute of the curve.
The Point primitives is the involute of the curve at curvePoints, starting at
{x1, y1} going in the direction of the curve.
The Line primitives connect these dots to curvePoints.";

involuteLineAndDotGP = Function[{fx, fy, tList, curvePoints, velociVectors},
 Module[{a, b, t, nArcLenF, f, arcLengths, involutePoints},
 {a, b} = N[{First[tList], Last[tList]}];
nArcLenF = NDSolve[{Derivative[1][f][t] == Sqrt[Derivative[1][fx][t]^2 + Derivative[1][fy][t]^2], f[a] == 0}, f, {t, a, b}][[1,1,2]];
 arcLengths = nArcLenF /@ tList;
involutePoints = arcLengths*Normalize /@ (-velociVectors) + curvePoints;
{MapThread[Line[{#1, #2}] & , {curvePoints, involutePoints}], Point /@ involutePoints}]
];

Clear[ nCurvatureK ];

nCurvatureK::usage="nCurvatureK[ velociVectors, acclVectors]
return the list of real numbers that represent the curvature of a curve.
velociVectors is a list of velocity vectors,
acclVectors is a list of accelation vectors.";

nCurvatureK = Function[
MapThread[
Function[ With[{xnorm=Norm[ #1 ]}, If[ xnorm == 0 , 0, (#2 . {-Last[#1], First[#1]})/ xnorm^3, 0] ] ] ,
{#1, #2}] ];

Clear[ curvatureGPs ];

curvatureGPs::usage="curvatureGPs[tPts, tangentVecs, acclVecs, circleEvoluteRadial]
return a list of graphics primitives of the form
{
{Circles[],..},
{{Line[{{a,b},{c,d}}],..}, {Point[],..}},
{{Line[{{e,f},{e,h}}],..}, {Point[],..}}
}
any of the list may be empty such as returning
{ {}, { {}, {} }, { {}, {} } },
but it always have the same structure.
The first list are osculating circles, second are lines and points for evolutes, third list are lines and points for radials.
Whether one of the list is empty depends whether True/False value given in circleEvoluteRadial.
circleEvoluteRadial is a list of boolean {osculatingCircle, evolute, radial}.";

curvatureGPs = Function[{tPts, tangentVs, acceVs, circleEvoluteRadial },
Function[{pts2, tanV2, accV2},
With[ {
normalUnitVList = Map[ Function[ Function[{-Last[#], First[#]}][Normalize[ # ]] ], tanV2],
radiusList = 1/nCurvatureK[ tanV2, accV2]
},
Block[ { centerVectors, centerPoints },
centerVectors = normalUnitVList*radiusList;
centerPoints = centerVectors + pts2;
{
If[circleEvoluteRadial[[1]], MapThread[Circle[#1, #2] & , {centerPoints, Abs[radiusList]}], {}],
If[circleEvoluteRadial[[2]], {MapThread[Line[{#1, #2}] & , {pts2, centerPoints}], Point /@ centerPoints}, {{}, {}}],
If[circleEvoluteRadial[[3]], {(Line[{{0, 0}, #1}] & ) /@ centerVectors, Point /@ centerVectors} , {{}, {}}]}
 ] ] ] @@
(Function[ (Delete@Position[ ((Norm@#) &) /@ #2 , 0.|0]) /@ {#1,#2,#3} ] [tPts, tangentVs, acceVs])
];

Clear[ curvatureNewGP ];

curvatureNewGP[fx_,fy_, tlist_, tPts_, tVelocities_, {cir_, evo_, rad_}] := Module[
{normalUnitVList, radiusList, centerVectors, centerPoints},
 normalUnitVList = Map[ Function[ Function[{-Last[#], First[#]}][Normalize[ # ]] ], tVelocities];
    radiusList = 1/ Table[ Evaluate@ FrenetSerretSystem[ {fx[t833],fy[t833]}, t833 ][[1,1]], {t833, tlist} ] ;
    centerVectors = normalUnitVList*radiusList;
    centerPoints = centerVectors + tPts;
    {If[cir, MapThread[Circle[#1, #2] & , {centerPoints, Abs[radiusList]}], {}],
     If[evo, {MapThread[Line[{#1, #2}] & , {tPts, centerPoints}], Point /@ centerPoints}, {{}, {}}],
     If[rad, {(Line[{{0, 0}, #1}] & ) /@ centerVectors, Point /@ centerVectors}, {{}, {}}]}];

Clear[ injectStyle ];

injectStyle::usage="injectStyle[graphicsPrims, xstyles]
pair the elements of these two lists.
xstyles may be a non-list, or a list.
if xstyles is non-list, simply return {xstyles, graphicsPrims}.
if xstyles is list of single item, simply return {First[xstyles], graphicsPrims}.

If xstyles is longer than graphicsPrims, then the rest of xstyles are dropped.
If xstyles is shorter, then xstyles is reused cyclically to match the length of graphicsPrims.

Each member of xstyles is flattened before being paired with elements in graphicsPrims.
For example,
graphicsPrims = {1,2,3,..},
xstyles = {{a1,a2},b,c,d,...},
 then the result is
{{a1,a2,1},{b,2},{c,3},...},
Warning: not
{{ {a1,a2} ,1},{b,2},{c,3},...}
When the length of xstyles is 0 or 1, special methods are used to save time.

you can test how style apply to graphics primitives by:

Graphics[ { Red, Circle[ {0,0},1 ]} ] (*work*)
Graphics[ { Red, {Circle[ {0,0},1 ]}} ] (*work*)
Graphics[ { {Red}, Circle[ {0,0},1 ]} ] (*no*)";

injectStyle[graPrims_, xstyle_] :=
Which[
Head[xstyle] =!= List, {xstyle, graPrims},
xstyle === {}, graPrims,
And[ Head[xstyle] === List, Length[xstyle] === 1 ], {Sequence @@ Flatten[{xstyle}] , graPrims},
True, Flatten /@ Transpose[
     {Take[Join @@ Table[xstyle, {Ceiling[Length[graPrims]/Length[xstyle]]}],
       Length[graPrims]], graPrims}]];

Options[PlaneCurvePlot] = Join[{
"PlotCurve" -> True,
"PlotDot" -> True,
"PlotParameterValue" -> False,
"NormalLineLength" -> 0,
"TangentLineLength" -> 0,
"CausticLineLength" -> 0,
"CausticOrigin" -> {0, 0},
"SecantConnection" -> {1, 1},
"ParallelCurveOffset" -> {0},
"ConchoidSetting" -> {{0, 0}, {0}},
"PedalPoint" -> None,
"InversionCircle" -> None,
"PlotOsculatingCircle" -> False,
"PlotInvolute" -> False,
"PlotEvolute" -> False,
"PlotRadial" -> False,
"CurveColorFunction" -> Function[Hue[0.65, #1, 0.8]],
"CurveStyle" -> {Thickness[0.005]},
"DotStyle" -> {{Red, PointSize[0.01]}},
"LineStyle" -> {{Green, Thickness[0.004]}},
"NormalLineStyle" -> {{Hue[0.6], Thickness[0.004]}},
"TangentLineStyle" -> {{Red, Thickness[0.004]}},
"SecantLineStyle" -> {{Red, Thickness[0.004]}},
"CausticLineStyle" -> {{Hue[0.5], Thickness[0.004]}}
}, Options[ParametricPlot]];

PlaneCurvePlot[{xf_Function, yf_Function}, {tmin_, tmax_, dt_:Automatic}, tlist_List:{Automatic}, xOpts:OptionsPattern[]] :=

With[{
xpedalPt = OptionValue["PedalPoint"],
invCircle = OptionValue["InversionCircle"],
xdotstyl = OptionValue["DotStyle"],
xlineStyle = OptionValue["LineStyle"],
tminN = N[tmin],
tmaxN = N[tmax],
dtN = If[dt === Automatic, N[tmax - tmin]/4., N[dt]]
},

Module[{
tVals,
tPoints,
tangentVecs,
acclVecs,
tDotGP,
tTextGP,
xnormalGP,
xtanlineGP,
xcausticLineGP,
xsecantLineGP,
xparallPointGP,
xconcoidLineGP,
xconcoidPointGP,
xpedalLineGP,
xpedalPointGP,
xinversPointGP,
xinvoluteLineGP,
xinvolutePointGP,
xoscCircleGP,
xevoluteLineGP,
xevolutePointGP,
xradialLineGP,
xradialPointGP,
xcurveGP},

tVals = If[tlist === {Automatic}, N@ Range[tmin, tmax, dtN], N[tlist]];
tPoints = Map[ Function[{xf[#1], yf[#1]}], tVals];

tangentVecs = If[
Or[
OptionValue["NormalLineLength"] != 0,
OptionValue["TangentLineLength"] != 0,
OptionValue["CausticLineLength"] != 0,
OptionValue["ParallelCurveOffset"] != {0},
OptionValue["PedalPoint"] =!= None,
OptionValue["PlotInvolute"],
OptionValue["PlotOsculatingCircle"],
OptionValue["PlotEvolute"],
OptionValue["PlotRadial"]
] ,
Map[Function[{Derivative[1][xf][#1], Derivative[1][yf][#1]} ], tVals],
 {} ];

acclVecs =
If[ Or[
OptionValue["PlotOsculatingCircle"],
OptionValue["PlotEvolute"],
OptionValue["PlotRadial"]] ,
Map[Function[{Derivative[2][xf][#1], Derivative[2][yf][#1]}], tVals],
{}
 ];

tTextGP = If[OptionValue@ "PlotParameterValue" === True, Table[Text[NumberForm[tVals[[t]], 4], tPoints[[t]], {-2, -1}], {t, 1, Length[tVals]}], {}];

tDotGP = If[OptionValue@ "PlotDot" === True, Point /@ tPoints, {}];
xnormalGP = If[OptionValue["NormalLineLength"] === 0, {}, normalLineGP[tPoints, tangentVecs, OptionValue["NormalLineLength"]]];

xtanlineGP = If[OptionValue["TangentLineLength"] === 0, {},
 tangentLineGP[tPoints, tangentVecs, OptionValue["TangentLineLength"]]];

xcausticLineGP = If[OptionValue["CausticLineLength"] === 0, {}, causticLineGP[tPoints, tangentVecs, N[OptionValue["CausticOrigin"]], N[OptionValue["CausticLineLength"]]]];

xsecantLineGP = If[OptionValue["SecantConnection"] === {1, 1}, {}, secantLineGP[tPoints, First[OptionValue["SecantConnection"]], Last[OptionValue["SecantConnection"]]]];

xparallPointGP = If[OptionValue["ParallelCurveOffset"] === {0}, {}, parallelDotGP[tPoints, tangentVecs, OptionValue["ParallelCurveOffset"]]];

 {xconcoidLineGP, xconcoidPointGP} = If[Last[OptionValue["ConchoidSetting"]] === {0}, {{}, {}}, conchoidLineAndDotGP[ First[OptionValue["ConchoidSetting"]], tPoints, Last[OptionValue["ConchoidSetting"]]]];

{xpedalLineGP, xpedalPointGP} = If[xpedalPt === None, {{}, {}}, pedalLineAndDotGP[xpedalPt, tPoints, tangentVecs]];

xinversPointGP = If[invCircle === None, {}, inversionDotGP[First@invCircle, Last@ invCircle, tPoints]];

{xinvoluteLineGP, xinvolutePointGP} = If[OptionValue["PlotInvolute"] === False, {{}, {}}, involuteLineAndDotGP[xf, yf, tVals, tPoints, tangentVecs]];

{xoscCircleGP, {xevoluteLineGP, xevolutePointGP}, {xradialLineGP, xradialPointGP}} =
If[acclVecs === {},
{{}, {{}, {}}, {{}, {}}},
curvatureGPs[tPoints, tangentVecs, acclVecs, {OptionValue["PlotOsculatingCircle"], OptionValue["PlotEvolute"], OptionValue["PlotRadial"]}]
];

xcurveGP = If[OptionValue[ "PlotCurve" ] === True,
First@ ParametricPlot[ {xf[t], yf[t]}, {t, tminN, tmaxN},
DisplayFunction -> Identity,
ColorFunction :> Function[ First[ Flatten[ {OptionValue["CurveColorFunction"]} ] ][#3] ],
PlotStyle :> OptionValue["CurveStyle"],
Evaluate[FilterRules[{xOpts}, Options[ ParametricPlot]]]
], {}];

Show[Graphics[ {
injectStyle[xradialLineGP, xlineStyle],
injectStyle[xoscCircleGP, xlineStyle],
injectStyle[xinvoluteLineGP, xlineStyle],
injectStyle[xevoluteLineGP, xlineStyle],
injectStyle[xtanlineGP , OptionValue["TangentLineStyle"]],
injectStyle[xnormalGP, OptionValue["NormalLineStyle"]] ,
injectStyle[xsecantLineGP, OptionValue["SecantLineStyle"]],
injectStyle[xcausticLineGP, OptionValue["CausticLineStyle"]],
injectStyle[xconcoidLineGP, xlineStyle],
injectStyle[xpedalLineGP, xlineStyle],
xcurveGP,
injectStyle[xinvolutePointGP, xdotstyl],
injectStyle[xevolutePointGP, xdotstyl],
injectStyle[xradialPointGP, xdotstyl],
injectStyle[xparallPointGP, xdotstyl],
injectStyle[xconcoidPointGP, xdotstyl],
injectStyle[xpedalPointGP, xdotstyl],

If[invCircle === None, {}, {Yellow, Circle @@ invCircle}],
injectStyle[xinversPointGP, xdotstyl],
injectStyle[tDotGP, xdotstyl],
tTextGP,
If[OptionValue["CausticLineLength"] === 0, {}, {Red, PointSize[0.02], Point[OptionValue["CausticOrigin"]]}],
If[Last[OptionValue["ConchoidSetting"]] === {0}, {}, {Blue, PointSize[0.02], Point[First[OptionValue["ConchoidSetting"]]]}],
If[xpedalPt === None, {}, {Blue, PointSize[0.02], Point[xpedalPt]}]
} ,
FilterRules[{xOpts}, Options[Graphics]], Axes -> True]]
] ] /;
And[NumericQ[tmin], NumericQ[tmax], Or[NumericQ[dt], SameQ[dt, Automatic]]];

PlaneCurvePlot[{fx_,fy_}, {t838_Symbol, tmin_, tmax_, dt_:Automatic}, tlist_List:{Automatic}, xOpts:OptionsPattern[]] :=
PlaneCurvePlot[ {Function[t838, fx ], Function[t838, fy]}, {tmin, tmax, dt}, tlist, xOpts] /;
And[NumericQ[tmin], NumericQ[tmax], Or[NumericQ[dt], SameQ[dt, Automatic]]];

End[];

EndPackage[];

