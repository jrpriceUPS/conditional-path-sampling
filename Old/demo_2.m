%a periodic cosine potential example

addpath('../Functions')

clear all;close all;

kB = 1.38064852*10^-23;

SDE.drift      =  @(x) 3*sin(x);
SDE.potential  =  @(x) 3*cos(x);
SDE.noise      =  1/2;
SDE.initial    =  -pi;


meta.weight    =  0.3;
meta.width     =  0.3;
meta.tempered  =  1;
meta.freq      =  10;
meta.temp      =  1000/kB;


domain.dt        =  0.1;
domain.endtime   =  500;
domain.periodic  =  1;
domain.edges     =  [-2*pi,2*pi];


plots.show        =  1;
plots.axes        =  [-2*pi,2*pi,-4,10];
plots.resolution  =  1000;


output = metadynamics(SDE,meta,domain,plots);