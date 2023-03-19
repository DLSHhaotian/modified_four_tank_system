function [xmean,s,xmeanp2s,xmeanm2s]=ScalarSampleMeanStdVar(x)

xmean = mean(x,1);
s = std(x,0,1);
xmeanp2s = xmean + 2*s;
xmeanm2s = xmean - 2*s;