clc; clear; close all;

a=[2 3;
    1 4];
inv(a)

h=[1 0 0 0;
    0 0 1 0];
h'*inv(a)
p=[1 1 1 11;
    2 2 2 2;
    1 1 1 1;
    2 2 2 2];
w= p*h'*inv(a)
p-w*inv(a)*w'

test=[2.7991e+00   2.4395e+00
   2.4395e+00   6.2094e+02];
inv(test)