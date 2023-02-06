clear; close all; clc

% Read in solution
sol     = readmatrix("./data/solution.txt");
halo    = readmatrix("./data/halo.txt");

% Plot trajectory
figure()
plot3(sol(:,2), sol(:,3), sol(:,4))

hold on
plot3(halo(:,1), halo(:,2), halo(:,3), "k")

grid on 
axis equal