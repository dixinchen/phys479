function [x,y,z] = plot_helmholtz_coils(radius,xCenter,yCenter,zCenter)

% Make an array for all the angles:
theta = linspace(0, 2 * pi, 2000);
% Create the x and y locations at each angle:
x = radius * cos(theta) + xCenter;
y = radius * sin(theta) + yCenter;
% Need to make a z value for every (x,y) pair:
z = zeros(1, numel(x)) + zCenter;