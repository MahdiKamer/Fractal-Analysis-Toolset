# Fractal analysis toolset
```mathematica
cm=72/2.54 ;(*centimetres*)
asr=1./GoldenRatio;(**aspect Ratio**)
ttick[min_,max_]:=Table[{N[i*10^-3],i},{i,Ceiling[min*10^3],Floor[max*10^3],1}];(*time plot ticks definition*)
ftick[min_,max_]:=Table[{i*1000.,i},{i,0,Floor[max/1000.],Floor[max/10000.]}];(*frequency plot ticks definition*)

tlabel=Style[" time (ms)",FontSize->10,FontFamily->"Times New Roman",Italic,Black];
flabel=Style[" frequency (kHz)",FontSize->10,FontFamily->"Times New Roman",Italic,Black];

SetOptions[ListLinePlot, PlotRange->All,AspectRatio->asr,ImageSize->18 cm,Axes->False,Frame->True,FrameStyle->Directive[Thickness[0.0015],FontFamily->"Times New Roman",FontSize->10,Black]];
SetOptions[ListPlot,PlotRange->All,AspectRatio->asr,ImageSize->18cm,Axes->False,Frame->True,FrameStyle->Directive[Thickness[0.0015],FontFamily->"Times New Roman",FontSize->10,Black]];
SetOptions[Plot,PlotRange->All,AspectRatio->asr,ImageSize->7.5 cm,Axes->False,Frame->True,FrameStyle->Directive[Thickness[0.0015],FontFamily->"Times New Roman",FontSize->10,Black]];
SetOptions[Histogram3D,PlotRange->All,AspectRatio->1,ImageSize->7.5 cm];
SetOptions[MatrixPlot,PlotRange->All,AspectRatio->1,ImageSize->7.5 cm,Axes->False,Frame->True,FrameStyle->Directive[Thickness[0.0015],FontFamily->"Times New Roman",FontSize->10,Black]];
SetOptions[ArrayPlot,PlotRange->All,AspectRatio->1,ImageSize->7.5 cm,Axes->False,Frame->True,FrameStyle->Directive[Thickness[0.0015],FontFamily->"Times New Roman",FontSize->10,Black]];
SetOptions[ReliefPlot,PlotRange->All,AspectRatio->1,ImageSize->7.5cm,Axes->False,Frame->True,FrameStyle->Directive[Thickness[0.0015],FontFamily->"Times New Roman",FontSize->10,Black]];
SetOptions[ListDensityPlot,PlotRange->All,AspectRatio->1,ImageSize->7.5 cm,Axes->False,Frame->True,FrameStyle->Directive[Thickness[0.0015],FontFamily->"Times New Roman",FontSize->10,Black]];
```


## Box Counting method
```mathematica
Dbox[dlis_,vmin_,vmax_]:=Module[{\[Epsilon],tmp,logbox,fitbox,minbox,maxbox,tbox},
\[Epsilon]=(vmax-vmin)/(2.^12.);
logbox=Table[Log[{1./((2.^i )\[Epsilon]),N@Total@Map[Count[Except[0]],BinCounts[dlis,{vmin,vmax,(2.^i )\[Epsilon]},{vmin,vmax,(2.^i )\[Epsilon]}]]}],{i,0.,12.,1.}];
logbox=Select[logbox,Abs[#[[2]]]!= Infinity&];
logbox]
```

## Correlation Method
##### Please only uncomment one of the correlation methods that suits your application. Alternatively you can rename modules if you want to utilise all three correlation methods.
### method 1
```mathematica
(*Dcr[dlis_,vmin_,vmax_]:=Module[{lognearest,\[Epsilon]},
\[Epsilon]=(vmax-vmin)/((2.^24.) Sqrt[2.]);
lognearest=N@Table[{Log[2.,(2.^i )\[Epsilon]],Log[2.,Length@Nearest[dlis,{0.,0.},{All,(2.^i )\[Epsilon]}]]},{i,0.,24.,1.}];
lognearest=Select[lognearest,Abs[#[[2]]]\[NotEqual] Infinity&];
lognearest]*)
```
### method 2
```mathematica
(*Dcr[dlis_,vmin_,vmax_]:=Module[{lognearest,\[Epsilon],lis},
\[Epsilon]=(vmax-vmin)/((2.^22.) Sqrt[2.]);
lis=Flatten[Partition[dlis,1,1000],1];
lognearest=Table[{Log[2.,(2.^i) \[Epsilon]],Log[2.,Total@Table[N@Length@Nearest[Drop[lis,k],lis[[k]],{All,(2.^i )\[Epsilon]}],{k,1,-1+Length@lis}]]},{i,0.,22.,1.}];
lognearest=Select[lognearest,Abs[#[[2]]]\[NotEqual] Infinity&];
lognearest]*)
```
### method 3
```mathematica
Dcr[dlis_,vmin_,vmax_]:=Module[{logcr,lis,\[Epsilon],dist,eliminate,pos},
\[Epsilon]=(vmax-vmin)/(2.^22.);
lis=Flatten[Partition[dlis,1,1000],1];
dist=DistanceMatrix[lis,DistanceFunction->SquaredEuclideanDistance];
eliminate[i_]:=UnitStep[ConstantArray[(2.^(2. i)) \[Epsilon]^2.,{Length@lis,Length@lis}]-dist]-IdentityMatrix[Length@lis];
pos[i_]:=SparseArray[eliminate[i]]["NonzeroPositions"];
logcr=Table[{Log[(2.^i) \[Epsilon]],Log[0.5 Length@pos[i]]},{i,0.,22.,1.}];
Select[logcr,Abs[#[[2]]]!= Infinity&]
]
```
## Fourier Analysis
```mathematica
ffans[lis_,fband_,fm_]:=Module[{lfcf,dmfft,freqlist},
lfcf=Fourier[lis,FourierParameters->{-1,-1}];
dmfft=Min[Floor[Length@lfcf/2.],Floor[(fband/fm)+1]];
freqlist=N@Table[(s-1) fm,{s,dmfft}];
Transpose@{freqlist,Abs@Take[lfcf,dmfft]}]
```
## Carier Producer
```mathematica
Fvc[fm_,df_,cm_,fband_]:= Module[{res=1,fsw,dt,nf,nd,lfc,dlfc,fv,fc,vc,xl,s,txt,gfc,gbhfc,gffc,gtbfc,gtcfc},
fsw=2000.;
nf=2. 1/fm/dt;
nd=1/fm/dt;
```
## selecting the desired carrier wave model by the variable "cm"
```mathematica
dt=.005/fsw;
lfc=Table[{i dt,SawtoothWave[{-1.,1.},fsw i dt]},{i,0.,nf}];
dlfc=Table[{SawtoothWave[{-1.,1.},fsw i dt],SawtoothWave[{-1,1},fsw ((i-nd)( dt))]},{i,nd,nf}];
txt={"Conventional Periodic Carrier","  "}
,
cm==2.,
fv[t_]:=fsw+df*Sin[2\[Pi] fm t];
fc[0.]=fv[0.];
vc[0.]=-1.;
dt=.0001/fsw;
Do[
If[vc[i]<=1,
vc[(i+1)]=vc[i]+2 fc[i] dt;fc[i+1]=fc[i],
fc[(i+1)]=fv[(i)*dt];vc[(i+1)]=-1],
{i,0.,nf}];
lfc=Table[{i dt,vc[(i)]},{i,0.,nf}];
dlfc=Table[{vc[(i)],vc[(i-nd)]},{i,nd,nf}];
txt={"Frequecy modulated with sine wave Periodic Carrier","Modulation Frequency="<>ToString[fm]<>"Hz, Deviation Frequency="<>ToString[df]<>"Hz"}
,
cm==3.,
xl[0]=N[.011];
fc[0.]=fsw;
vc[0.]=-1.;
dt=.0005/fsw;
s=Table[xl[i+1]=4.xl[i](1-xl[i]),{i,0.,nf}];
Do[
If[vc[i]<=1.,
vc[(i+1)]=vc[i]+2. fc[i] dt;fc[i+1]=fc[i],
fc[(i+1)]=fsw+df*N[s[[i+1]]];vc[(i+1)]=-1.],
{i,0.,nf}];
lfc=Table[{i dt,vc[(i)]},{i,0.,nf}];
dlfc=Table[{vc[(i)],vc[(i-nd)]},{i,nd,nf}];
txt={"Frequecy modulated with logestic map Chaotic Carrier","Modulation Frequency="<>ToString[fm]<>"Hz, Deviation Frequency="<>ToString[df]<>"Hz"}
,
cm==4.,
xl[0]=N[1./GoldenRatio,10000];
fc[0.]=fsw;
vc[0.]=-1.;
dt=.0005/fsw;
s=Table[xl[i+1]=N[Mod[2xl[i],1],1300],{i,0.,nf}];
Do[
If[vc[i]<=1.,
vc[(i+1)]=vc[i]+2. fc[i] dt;fc[i+1]=fc[i],
fc[(i+1)]=fsw+df*N[s[[i+1]]];vc[(i+1)]=-1.],
{i,0.,nf}];
lfc=Table[{i dt,vc[(i)]},{i,0.,nf}];
dlfc=Table[{vc[(i)],vc[(i-nd)]},{i,nd,nf}];
txt={"Frequecy modulated with bernoli map Chaotic Carrier","Modulation Frequency="<>ToString[fm]<>"Hz, Deviation Frequency="<>ToString[df]<>"Hz"}
];

gffc=ffans[lfc[[1;;Floor@nd,2]],fband,fm];
gtbfc=Dbox[dlfc,-1.,1.];
gtcfc=Dcr[dlfc,-1.,1.];
```
## return structure of the module
```mathematica
{Text[Style[txt[[1]],FontSize->14,Red]],Text[Style[txt[[2]],FontSize->14,Red]],Text[Style["RMS Value="NumberForm[RootMeanSquare[lfc[[All,2]]],{5,4}],FontSize->14]],0,0,lfc,Flatten[Partition[dlfc,1,res],1],gffc,gtbfc,gtcfc,dt}
]
```
## Common mode circuit modeler
```mathematica
Fcmv[lfc_,dt_,fm_,mi_,fbandv_,fbandi_]:= Module[{res=1,lisvc,lt,nf,nd,ndd,lrefA,lrefB,lrefC,lv0inv,dlv0inv,vno,lvno,dlvno,vbo,lvbo,dlvbo,ino,lino,dlino,ibo,libo,dlibo,slic,vcm,vs,p,vcwf0,vcwf0p,vcrf0,vcrf0p,vcww0p,vcwr0p,il0,vcb0,slb,gv0inv,gbhv0inv,gfv0inv,gtbv0inv,gtcv0inv,gvno,gbhvno,gfvno,gtbvno,gtcvno,gino,gbhino,gfino,gtbino,gtcino,
gvbo,gbhvbo,gfvbo,gtbvbo,gtcvbo,gibo,gbhibo,gfibo,gtbibo,gtcibo,il,vcwf,vcrf,vcb,tpl,
(**************************************Lumped Model Parameters**************************************************************)
\[Epsilon]0=1./(36.0*\[Pi]*10^9),Cb,Cwr1=7.5*^-12,Cwr2=5*^-12,Crf=830*^-12,Rb=0.24,Rr=0.01,Rww=425.,Cww=3.027999712106278*^-10, Cwf=1.1758462604578913*^-9,Lww=0.004165082331641172},
Cb=(2.*(2.5*\[Epsilon]0*(2*\[Pi]*0.004*0.02)/0.002)+3*\[Epsilon]0*(0.06*0.2-2*\[Pi]*(0.006)^2)/0.0012);
lisvc=lfc[[All,2]];
lt=lfc[[All,1]];
nf=Length@lisvc-1;
ndd=Round[1/fm/dt];
nd=1/fm/dt;
lrefA=Table[mi Sin[2.\[Pi] fm i dt],{i,0.,nf}];
lrefB=Table[mi Sin[2.\[Pi] fm i dt-(2.\[Pi])/3.],{i,0.,nf}];
lrefC=Table[mi Sin[2.\[Pi] fm i dt-(4.\[Pi])/3.],{i,0.,nf}];
lv0inv=1./3. (Sign[lrefA-lisvc]+Sign[lrefB-lisvc]+Sign[lrefC-lisvc]);

vcm=Interpolation[Transpose@{lt,lv0inv}];
slic=First@Solve[{Cww p vs == (Cww+Cwf+Cwr2)p vcwf-Cwr2 p vcrf, Cwr1 p (vs-vcrf)+Cwr2 p( vcwf-vcrf)==Crf p vcrf+2. Cb  p vcb,0.5*^-7==(vcrf-vcb)/(Rr+Rb/2.)},{vcwf,vcrf,vcb}];
vcwf0=vcwf/.slic/.vs-> vcm[0.];
vcwf0p=vcwf/.slic/.vs-> vcm'[0.];
vcrf0=vcrf/.slic/.vs-> vcm[0.];
vcrf0p=vcrf/.slic/.vs-> vcm'[0.];
vcww0p=vcm'[0.]-vcwf0p;
vcwr0p=vcwf0p-vcrf0p;
il0=-Cww vcww0p+Cwf vcwf0p+Cwr2 vcwr0p;
vcb0=vcb/.slic/.vs-> vcm[0.];
slb=First@NDSolve[{-vcm[t]+Rww il[t]+ Lww il'[t]+vcwf[t]== 0.,
2.Cb vcb'[t]==(vcrf[t]-vcb[t])/(Rr+Rb/2.),
Cwr2 (vcwf'[t]-vcrf'[t])+Cwr1 (vcm'[t]-vcrf'[t])==(vcrf[t]-vcb[t])/(Rr+Rb/2.)+Crf vcrf'[t],
il[t]+Cww (vcm'[t]-vcwf'[t])==Cwf vcwf'[t]+Cwr2 (vcwf'[t]-vcrf'[t]),
il[0.]== il0, vcwf[0.]== vcwf0, vcrf[0.]==vcrf0,vcb[0.]==vcb0},{il,vcwf,vcrf,vcb},{t,0.,vcm[[1,1,2]]},MaxSteps->\[Infinity](*,Method->{"EquationSimplification"->"Residual"}*)];

vbo[t_]:=Evaluate[(vcb[t]+Rb/2. (vcrf[t]-vcb[t])/(Rr+Rb/2.))/.slb];
ibo[t_]:=Evaluate[0.5((vcrf[t]-vcb[t])/(Rr+Rb/2.))/.slb];
vno[t_]:=Evaluate[vcwf[t]/.slb];
ino[t_]:=Evaluate[(il[t]+Cww(vcm[t]-vcwf[t])+Cwr1(vcm[t]-vcrf[t]))/.slb];

lvno=Table[vno[i dt],{i,0.,nf}];
lino=Table[ino[i dt],{i,0.,nf}];
libo=Table[ibo[i dt],{i,0.,nf}];
lvbo=Table[vbo[i dt],{i,0.,nf}];

dlv0inv=Table[{lv0inv[[i]],lv0inv[[i-ndd]]},{i,ndd+1,2 ndd +1}];
dlvno=Table[{lvno[[i]],lvno[[i-ndd]]},{i,ndd+1,2 ndd +1}];
dlino=Table[{lino[[i]],lino[[i-ndd]]},{i,ndd+1,2 ndd +1}];
dlvbo=Table[{lvbo[[i]],lvbo[[i-ndd]]},{i,ndd+1,2 ndd +1}];
dlibo=Table[{libo[[i]],libo[[i-ndd]]},{i,ndd+1,2 ndd+1}];

tpl=Floor[nf/5.];

gfv0inv=ffans[lv0inv[[1;;Floor@nd]],fbandv,fm];
gtbv0inv=Dbox[dlv0inv,Min@lv0inv,Max@lv0inv];
gtcv0inv=Dcr[dlv0inv,Min@lv0inv,Max@lv0inv];

gfvno=ffans[lvno[[1;;Floor@nd]],fbandv,fm];
gtbvno=Dbox[dlvno,Min@lvno,Max@lvno];
gtcvno=Dcr[dlvno,Min@lvno,Max@lvno];

gfino=ffans[lino[[1;;Floor@nd]],fbandi,fm];
gtbino=Dbox[dlino,Min@lino,Max@lino];
gtcino=Dcr[dlino,Min@lino,Max@lino];

gfvbo=ffans[lvbo[[1;;Floor@nd]],fbandv,fm];
gtbvbo=Dbox[dlvbo,Min@lvbo,Max@lvbo];
gtcvbo=Dcr[dlvbo,Min@lvbo,Max@lvbo];

gfibo=ffans[libo[[1;;Floor@nd]],fbandi,fm];
gtbibo=Dbox[dlibo,Min@libo,Max@libo];
gtcibo=Dcr[dlibo,Min@libo,Max@libo];
```
## return structure of the module
```mathematica
{Text[Style["Line-End Voltage, \!\(\*SubscriptBox[\(v\), \(0  inv\)]\)",FontSize->14,Bold,Brown]],Text[Style["Modulation Index="NumberForm[mi,{3,2}],FontSize->14,Bold,Brown]],Text[Style["RMS Value="NumberForm[N@RootMeanSquare[lv0inv],{5,4}],FontSize->14]],0,0,Transpose@{lt[[1;;tpl]],lv0inv[[1;;tpl]]},Flatten[Partition[dlv0inv,1,res],1],gfv0inv,gtbv0inv,gtcv0inv,
Text[Style["Open-End Voltage, \!\(\*SubscriptBox[\(v\), \(no\)]\)",FontSize->14,Bold,Brown]],Null,Text[Style["RMS Value="NumberForm[N@RootMeanSquare[lvno],{5,4}],FontSize->14]],0,0,Transpose@{lt[[1;;tpl]],lvno[[1;;tpl]]},Flatten[Partition[dlvno,1,res],1],gfvno,gtbvno,gtcvno,
Text[Style["Coil current, \!\(\*SubscriptBox[\(i\), \(no\)]\)",FontSize->14,Bold,Brown]],Null,Text[Style["RMS Value="EngineeringForm@N@RootMeanSquare[lino],FontSize->14]],0,0,Transpose@{lt[[1;;tpl]],lino[[1;;tpl]]},Flatten[Partition[dlino,1,res],1],gfino,gtbino,gtcino,
Text[Style["Bearing voltage, \!\(\*SubscriptBox[\(v\), \(bo\)]\)",FontSize->14,Bold,Brown]],Null,Text[Style["RMS Value="EngineeringForm@N@RootMeanSquare[lvbo],FontSize->14]],0,0,Transpose@{lt[[1;;tpl]],lvbo[[1;;tpl]]},Flatten[Partition[dlvbo,1,res],1],gfvbo,gtbvbo,gtcvbo,
Text[Style["Bearing Current, \!\(\*SubscriptBox[\(i\), \(bo\)]\)",FontSize->14,Bold,Brown]],Null,Text[Style["RMS Value="EngineeringForm@N@RootMeanSquare[libo],FontSize->14]],0,0,Transpose@{lt[[1;;tpl]],libo[[1;;tpl]]},Flatten[Partition[dlibo,1,res],1],gfibo,gtbibo,gtcibo}
]
```
## Simulation Module
```mathematica
Fsim[fm_,df_,cm_,fband_,mi_]:= Module[{ CarWav,CMVsol},
CarWav=Fvc[fm,df,cm,fband];
CMVsol=Fcmv[CarWav[[6]],CarWav[[11]],fm,mi,fband,fband];
Flatten[{CarWav[[1;;10]],CMVsol},1]
]
```
## graphic functions
```mathematica
gfft=(ListPlot[#,PlotLabel->"Fourier Spectrum",Filling->Axis,FrameTicks->{{All,None},{ftick,None}},FrameLabel->{{,},{flabel, }}])&;
gden=(MatrixPlot[BinCounts[#,{Min[#[[All,1]]],Max[#[[All,1]]],(Max[#[[All,1]]]-Min[#[[All,1]]])/10.},{Min[#[[All,2]]],Max[#[[All,2]]],(Max[#[[All,2]]]-Min[#[[All,2]]])/10.}], PlotLabel->"Matrix Plot of the Attractor Matrix"])&;
(*gclus=(ListPlot[FindClusters[#], PlotLabel\[Rule]"Clusterized Plot of the Attractor Matrix"])&;*)
gtim=ListLinePlot[#,PlotLabel->"time Wave",FrameTicks->{{All,None},{ttick,None}},FrameLabel->{{,},{tlabel, }}]&;
gatt=ListLinePlot[#,PlotLabel->"Reconstructed Delayed Trajectories"]&;
gcr=Show[ListPlot[#1,PlotStyle->{Black,PointSize[Medium]}],Plot[Normal@LinearModelFit[Drop[Drop[#1,#2],#3],y,y]/.y-> x,{x,Min[#1[[All,1]]],Max[#1[[All,1]]]}],PlotLabel->"Correlation Dimension",FrameLabel->{"log[\!\(\*SuperscriptBox[\(2\), \(i\)]\)\[Epsilon]]","log[N(\!\(\*SuperscriptBox[\(2\), \(i\)]\)\[Epsilon])]"}]&;
tcr=Text[Style["Correlation Dimension="NumberForm[Normal[LinearModelFit[Drop[Drop[#1,#2],#3],x,x]][[2,1]],{5,4}],FontSize->14]]&;
gbox=Show[ListPlot[#1,PlotStyle->{Black,PointSize[Medium]}],Plot[Normal@LinearModelFit[Drop[Drop[#1,#2],#3],y,y]/.y-> x,{x,Min[#1[[All,1]]],Max[#1[[All,1]]]}],PlotLabel->"Box Counting Dimension",FrameLabel->{"log[1/(\!\(\*SuperscriptBox[\(2\), \(i\)]\)\[Epsilon])]","log[N(\!\(\*SuperscriptBox[\(2\), \(i\)]\)\[Epsilon])]"}]&;
tbox=Text[Style["BoxCounting Dimension="NumberForm[Normal[LinearModelFit[Drop[Drop[#1,#2],#3],x,x]][[2,1]],{5,4}],FontSize->14]]&;
Fplot=Grid[{#1[[1;;3]],{gtim[(#1[[6]])[[1;;Floor[Length[#1[[6]]]/10.]]]],gatt[#1[[7]]],gfft[#1[[8]]]},{tbox[#1[[9]],#2,#3],tcr[#1[[10]],#4,#5],""},{gbox[#1[[9]],#2,#3],gcr[#1[[10]],#4,#5],gden[#1[[7]]]},
#[[11;;13]],{gtim[(#1[[16]])],gatt[#1[[17]]],gfft[#1[[18]]]},{tbox[#1[[19]],#6,#7],tcr[#1[[20]],#8,#9],""},{gbox[#1[[19]],#6,#7],gcr[#1[[20]],#8,#9],gden[#1[[17]]]},
#[[21;;23]],{gtim[(#1[[26]])],gatt[#1[[27]]],gfft[#1[[28]]]},{tbox[#1[[29]],#10,#11],tcr[#1[[30]],#12,#13],""},{gbox[#1[[29]],#10,#11],gcr[#1[[30]],#12,#13],gden[#1[[27]]]},
#[[31;;33]],{gtim[(#1[[36]])],gatt[#1[[37]]],gfft[#1[[38]]]},{tbox[#1[[39]],#14,#15],tcr[#1[[40]],#16,#17],""},{gbox[#1[[39]],#14,#15],gcr[#1[[40]],#16,#17],gden[#1[[37]]]},
#[[41;;43]],{gtim[(#1[[46]])],gatt[#1[[47]]],gfft[#1[[48]]]},{tbox[#1[[49]],#18,#19],tcr[#1[[50]],#20,#21],""},{gbox[#1[[49]],#18,#19],gcr[#1[[50]],#20,#21],gden[#1[[47]]]},
#[[51;;53]],{gtim[(#1[[56]])],gatt[#1[[57]]],gfft[#1[[58]]]},{tbox[#1[[59]],#22,#23],tcr[#1[[60]],#24,#25],""},{gbox[#1[[59]],#22,#23],gcr[#1[[60]],#24,#25],gden[#1[[57]]]}},Background->{None,{Lighter[Yellow,.9],{White,Lighter[Blend[{Blue,Green}],.99]}}},Dividers->{{Darker[Gray,.6],{Lighter[Gray,.5]},Darker[Gray,.6]},{Darker[Gray,.6],Darker[Gray,.6],{False},Darker[Gray,.6]}},Alignment->{{Center,Center}},Frame->Darker[Gray,.6],ItemStyle->14,Spacings->{Automatic,.8}]&;
```

```mathematica
SetDirectory[NotebookDirectory[]];
Put[1,"ComModFracDimFuncs.nbd"]
Save["ComModFracDimFuncs.nbd",{Fsim,Fplot,ttick,ftick,tlabel,flabel}]
```