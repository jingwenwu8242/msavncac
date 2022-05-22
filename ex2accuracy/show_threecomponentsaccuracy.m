clear all; 

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=256; ny =256;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);


ss=sprintf('./dataout/3componentaccuracy/datac.m'); accuracyphi = load(ss);
ss=sprintf('./dataout/3componentaccuracy/datac2.m'); accuracyphi2 = load(ss);
ss=sprintf('./dataout/3componentaccuracy/datac3.m'); accuracyphi3 = load(ss);

ss=sprintf('./dataout/3component2dt/datac.m'); phidt2 = load(ss);
ss=sprintf('./dataout/3component2dt/datac2.m'); phi2dt2= load(ss);
ss=sprintf('./dataout/3component2dt/datac3.m'); phi3dt2 = load(ss);

ss=sprintf('./dataout/3component4dt/datac.m'); phidt4 = load(ss);
ss=sprintf('./dataout/3component4dt/datac2.m'); phi2dt4 = load(ss);
ss=sprintf('./dataout/3component4dt/datac3.m'); phi3dt4 = load(ss);

ss=sprintf('./dataout/3component8dt/datac.m'); phidt8 = load(ss);
ss=sprintf('./dataout/3component8dt/datac2.m'); phi2dt8 = load(ss);
ss=sprintf('./dataout/3component8dt/datac3.m'); phi3dt8 = load(ss);

ss=sprintf('./dataout/3component16dt/datac.m'); phidt16 = load(ss);
ss=sprintf('./dataout/3component16dt/datac2.m'); phi2dt16 = load(ss);
ss=sprintf('./dataout/3component16dt/datac3.m'); phi3dt16 = load(ss);


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

Adt4error=sqrt(sum(sum((accuracyA-Adt4).^2))*h^2)
Bdt4error=sqrt(sum(sum((accuracyB-Bdt4).^2))*h^2)
Cdt4error = sqrt(sum(sum((accuracyC-Cdt4).^2))*h^2)

Adt8error=sqrt(sum(sum((accuracyA-Adt8).^2))*h^2)
Bdt8error=sqrt(sum(sum((accuracyB-Bdt8).^2))*h^2)
Cdt8error = sqrt(sum(sum((accuracyC-Cdt8).^2))*h^2)


Adt16error=sqrt(sum(sum((accuracyA-Adt16).^2))*h^2)
Bdt16error=sqrt(sum(sum((accuracyB-Bdt16).^2))*h^2)
Cdt16error = sqrt(sum(sum((accuracyC-Cdt16).^2))*h^2)


Arate1=log2(Adt4error/Adt2error)
Brate1=log2(Bdt4error/Bdt2error)
Crate1=log2(Cdt4error/Cdt2error)


Arate2=log2(Adt8error/Adt4error)
Brate2=log2(Bdt8error/Bdt4error)
Crate2=log2(Cdt8error/Cdt4error)

Arate3=log2(Adt16error/Adt8error)
Brate3=log2(Bdt16error/Bdt8error)
Crate3=log2(Cdt16error/Cdt8error)













