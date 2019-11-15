function logax(ylmt)
%function logax(ylmt)
%LOGAX changes the depth (y) axis limits of all subplots plotted by LPLOT.
%ylmt=[ymin ymax]

%Written by T. Mukerji

h=findobj(gcf,'type','axes');

for k=1:length(h) set(h(k),'ylim',ylmt); end;
