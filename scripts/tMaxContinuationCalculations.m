
clear; close all; clc

% Continuation steps from Zheng and Topputo
tMaxs    = [10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0,0.6]';
tfs      = [8.6404,9.5548,10.6174,11.8946,15.1533,16.8524,20.8663,27.4773,37.8965,84.3688,140.2678]';

% Fit curve to continuation steps
f   = fit(tMaxs,tfs,'exp2');

% Plot fit
plot(f,tMaxs,tfs)