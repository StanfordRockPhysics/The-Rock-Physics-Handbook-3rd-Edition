function [vp1,vs]=biothfgs(vpdry,vsdry,k0,mu0,ro0,rofl,kfl,por,alfa)
%function [vp1,vs]=biothf(vpdry,vsdry,k0,mu0,ro0,rofl,kfl,por,alfa)

rodry=(1-por)*ro0; ro=(1-por)*ro0+por*rofl;
mudry=rodry*vsdry^2; kdry=rodry*vpdry^2-(4/3)*mudry; b=kdry/k0;
robiot=ro0*(1-por)+por*rofl*(1-1/alfa);
ro12=(1-alfa)*por*rofl; ro11=(1-por)*ro0-ro12; ro22=por*rofl*alfa;
ro=ro11+2*ro12+ro22; rol=(ro12+ro22)/por; roc=ro22/por^2;
den=((1-por-b)/k0 + por/kfl);
h=(1-b)^2/den + kdry + (4/3)*mudry;
k=(1-b)/den; l=1/den;
vp1=sqrt((l*ro + h*roc - 2*rol*k)/(ro*roc-rol^2));
vs=sqrt(mudry/robiot);
