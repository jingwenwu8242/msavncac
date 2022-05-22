clear all; 

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=256; ny =256;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);


ss=sprintf('./data/accuracy/datac1.m'); accuracyphi = load(ss);
ss=sprintf('./data/accuracy/datac2.m'); accuracyphi2 = load(ss);
ss=sprintf('./data/accuracy/datac3.m'); accuracyphi3 = load(ss);
ss=sprintf('./data/accuracy/datau.m'); accuracyu = load(ss);
ss=sprintf('./data/accuracy/datav.m'); accuracyv = load(ss);
ss=sprintf('./data/accuracy/datap.m'); accuracyp = load(ss);

ss=sprintf('./data/accuracy2dt/datac1.m'); phidt2 = load(ss);
ss=sprintf('./data/accuracy2dt/datac2.m'); phi2dt2= load(ss);
ss=sprintf('./data/accuracy2dt/datac3.m'); phi3dt2 = load(ss);
ss=sprintf('./data/accuracy2dt/datau.m'); accuracyudt2 = load(ss);
ss=sprintf('./data/accuracy2dt/datav.m'); accuracyvdt2 = load(ss);
ss=sprintf('./data/accuracy2dt/datap.m'); accuracypdt2= load(ss);

ss=sprintf('./data/accuracy4dt/datac1.m'); phidt4 = load(ss);
ss=sprintf('./data/accuracy4dt/datac2.m'); phi2dt4 = load(ss);
ss=sprintf('./data/accuracy4dt/datac3.m'); phi3dt4 = load(ss);
ss=sprintf('./data/accuracy4dt/datau.m'); accuracyudt4 = load(ss);
ss=sprintf('./data/accuracy4dt/datav.m'); accuracyvdt4 = load(ss);
ss=sprintf('./data/accuracy4dt/datap.m'); accuracypdt4 = load(ss);



ss=sprintf('./data/accuracy8dt/datac1.m'); phidt8 = load(ss);
ss=sprintf('./data/accuracy8dt/datac2.m'); phi2dt8 = load(ss);
ss=sprintf('./data/accuracy8dt/datac3.m'); phi3dt8 = load(ss);
ss=sprintf('./data/accuracy8dt/datau.m'); accuracyudt8 = load(ss);
ss=sprintf('./data/accuracy8dt/datav.m'); accuracyvdt8 = load(ss);
ss=sprintf('./data/accuracy8dt/datap.m'); accuracypdt8 = load(ss);

ss=sprintf('./data/accuracy16dt/datac1.m'); phidt16 = load(ss);
ss=sprintf('./data/accuracy16dt/datac2.m'); phi2dt16 = load(ss);
ss=sprintf('./data/accuracy16dt/datac3.m'); phi3dt16 = load(ss);
ss=sprintf('./data/accuracy16dt/datau.m'); accuracyudt16 = load(ss);
ss=sprintf('./data/accuracy16dt/datav.m'); accuracyvdt16 = load(ss);
ss=sprintf('./data/accuracy16dt/datap.m'); accuracypdt16 = load(ss);


accuracyA = accuracyphi;
accuracyB = accuracyphi2;
accuracyC = 1-accuracyA-accuracyB;





Adt2=phidt2;
Bdt2=phi2dt2;
Cdt2=1-Adt2-Bdt2;


Adt4=phidt4;
Bdt4=phi2dt4;
Cdt4=1-Adt4-Bdt4;

Adt8=phidt8;
Bdt8=phi2dt8;
Cdt8=1-Adt8-Bdt8;

Adt16=phidt16;
Bdt16=phi2dt16;
Cdt16=1-Adt16-Bdt16;

Adt2error=sqrt(sum(sum((accuracyA-Adt2).^2))*h^2)
Bdt2error=sqrt(sum(sum((accuracyB-Bdt2).^2))*h^2)
Cdt2error = sqrt(sum(sum((accuracyC-Cdt2).^2))*h^2)
Udt2error = sqrt(sum(sum((accuracyu-accuracyudt2).^2))*h^2)
Vdt2error = sqrt(sum(sum((accuracyv-accuracyvdt2).^2))*h^2)
Pdt2error = sqrt(sum(sum((accuracyp-accuracypdt2).^2))*h^2)



Adt4error=sqrt(sum(sum((accuracyA-Adt4).^2))*h^2)
Bdt4error=sqrt(sum(sum((accuracyB-Bdt4).^2))*h^2)
Cdt4error = sqrt(sum(sum((accuracyC-Cdt4).^2))*h^2)
Udt4error = sqrt(sum(sum((accuracyu-accuracyudt4).^2))*h^2)
Vdt4error = sqrt(sum(sum((accuracyv-accuracyvdt4).^2))*h^2)
Pdt4error = sqrt(sum(sum((accuracyp-accuracypdt4).^2))*h^2)

Adt8error=sqrt(sum(sum((accuracyA-Adt8).^2))*h^2)
Bdt8error=sqrt(sum(sum((accuracyB-Bdt8).^2))*h^2)
Cdt8error = sqrt(sum(sum((accuracyC-Cdt8).^2))*h^2)
Udt8error = sqrt(sum(sum((accuracyu-accuracyudt8).^2))*h^2)
Vdt8error = sqrt(sum(sum((accuracyv-accuracyvdt8).^2))*h^2)
Pdt8error = sqrt(sum(sum((accuracyp-accuracypdt8).^2))*h^2)


Adt16error=sqrt(sum(sum((accuracyA-Adt16).^2))*h^2)
Bdt16error=sqrt(sum(sum((accuracyB-Bdt16).^2))*h^2)
Cdt16error = sqrt(sum(sum((accuracyC-Cdt16).^2))*h^2)
Udt16error = sqrt(sum(sum((accuracyu-accuracyudt16).^2))*h^2)
Vdt16error = sqrt(sum(sum((accuracyv-accuracyvdt16).^2))*h^2)
Pdt16error = sqrt(sum(sum((accuracyp-accuracypdt16).^2))*h^2)


Arate1=log2(Adt4error/Adt2error)
Brate1=log2(Bdt4error/Bdt2error)
Crate1=log2(Cdt4error/Cdt2error)
Urate1=log2(Udt4error/Udt2error)
Vrate1=log2(Vdt4error/Vdt2error)
Prate1=log2(Pdt4error/Pdt2error)

Arate2=log2(Adt8error/Adt4error)
Brate2=log2(Bdt8error/Bdt4error)
Crate2=log2(Cdt8error/Cdt4error)
Urate2=log2(Udt8error/Udt4error)
Vrate2=log2(Vdt8error/Vdt4error)
Prate2=log2(Pdt8error/Pdt4error)

Arate3=log2(Adt16error/Adt8error)
Brate3=log2(Bdt16error/Bdt8error)
Crate3=log2(Cdt16error/Cdt8error)
Urate3=log2(Udt16error/Udt8error)
Vrate3=log2(Vdt16error/Vdt8error)
Prate3=log2(Pdt16error/Pdt8error)












