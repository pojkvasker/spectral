clear all
close all

load('lynxdata.mat')
load('sunspotdata.mat')

n1 = length(loglynx);
n2 = length(lynx);
n3 = length(sunspot);

N1 = n1/10;
[mu1,sig1] = lsar(loglynx,N1);