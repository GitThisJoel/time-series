clear
close all

addpath("../CourseMaterial/Code/data");
addpath("../CourseMaterial/Code/functions");

rng(0);

%%

n = 500; % Number o f s am ple s
A3 = [ 1 .5 ] ;
C3 = [ 1 -.3 .2 ] ;
w = sqrt ( 2 ) * randn( n + 100 , 1 ) ;
x = filter(C3 , A3,w ) ; % Crea te t h e i n p u t
A1 = [ 1 -.65];
A2 = [ 1 .90 .78 ] ;
C = 1 ;
B = [ 0 0 0 0 .4 ] ;
e = sqrt ( 1.5 ) * randn( n + 100 , 1 ) ;
y = filter(C, A1, e ) + filter(B, A2, x ) ; % Crea te t h e o u t p u t
x = x ( 101 : end); 
y = y ( 101 : end); % Omit i n i t i a l s am ple s
clear A1 A2 C B e w A3 C3 