%a double gaussian well example

addpath('../Functions')

clear all;close all;

kB = 1.38064852*10^-23;

SDE.drift      =  @(x) 5*gaussian_deriv(x,2/0.75,1/sqrt(2),1)...
                       + 10*gaussian_deriv(x,-2/0.75,1/sqrt(2),1)...
                       - (4*(x+4).^3).*(x<-4.3193) - (4*(x-4).^3).*(x>4.25882);
SDE.potential  =  @(x) -5*gaussian(x,2/0.75,1/sqrt(2),1)...
                       - 10*gaussian(x,-2/0.75,1/sqrt(2),1)...
                       + ((x+4).^4-1.690133).*(x<-4.3193)...
                       + ((x-4).^4-0.845067).*(x>4.25882);
SDE.noise      =  @(x) 1/2;
SDE.initial    =  -2/0.75;


meta.weight    =  0.25;
meta.width     =  0.25;
meta.tempered  =  0;
meta.freq      =  10;
meta.temp      =  1000/kB;


domain.dt        =  0.01;
domain.endtime   =  10;
domain.periodic  =  0;


plots.show        =  1;
plots.axes        =  [-4,4,-10,5];
plots.resolution  =  1000;


output = metadynamics(SDE,meta,domain,plots);