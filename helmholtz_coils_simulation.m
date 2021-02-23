clc;clear;close all

%% Define constants and parameters

% Define global variables
global u0

u0=4*pi*1e-7; % Permeability of free space

% Define coil parameters
I0=4; % Coil current in Amps
a=.4; % Coil radius in m

x_p1=0; y_p1=0; z_p1=0; % Define coordinates of coil 1 center point
x_p2=0; y_p2=0; z_p2=a; % Define coordinates of coil 2 center point

%% Plot the helmholtz coils

[xx1,yy1,zz1] = plot_helmholtz_coils(a,x_p1,y_p1,z_p1); % Calculate coordinates of coil 1
[xx2,yy2,zz2] = plot_helmholtz_coils(a,x_p2,y_p2,z_p2); % Calculate coordinates of coil 2

f1 = figure;
plot3(y_p1, z_p1, x_p1, 'r.', 'LineWidth', 2, 'MarkerSize', 10); % Plot the center point of coil 1
hold on
plot3(y_p2, z_p2, x_p2, 'r.', 'LineWidth', 2, 'MarkerSize', 10); % Plot the center point of coil 2
plot3(yy1, zz1, xx1, 'b-', 'LineWidth', 2); % Plot coil 1
plot3(yy2, zz2, xx2, 'b-', 'LineWidth', 2); % Plot coil 2
grid on
axis('square');
xlabel('y', 'FontSize', 20);
ylabel('z', 'FontSize', 20);
zlabel('x', 'FontSize', 20);

%% Calculate magnetic field

% Input mesh of points in 2D plane
x=0;
% This is a 2d plane over the x=0 plane that extends away from the coils in the yz plane. 
[y,z]=meshgrid(linspace(-a/2,a/2,25),linspace(0,a,25)); 

[Bx1,By1,Bz1] = magnetic_field_current_loop(x,y,z,x_p1,y_p1,z_p1,a,I0); % Calculate field from first coil
[Bx2,By2,Bz2] = magnetic_field_current_loop(x,y,z,x_p2,y_p2,z_p2,a,I0); % Calculate field from second coil

% Add the components
Bx=Bx1+Bx2;
By=By1+By2;
Bz=Bz1+Bz2;

% Calculate the magnitude of the vector
B_mag=sqrt(Bx.^2+By.^2+Bz.^2);

%% Plot helmholtz coils and the corresponding magnetic field

f2 = figure;
% surf(y,z,B_mag, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');
surf(y,z,B_mag);
hold on
% hold on;
% plot(y_p1, z_p1, 'r.', 'LineWidth', 2, 'MarkerSize', 10);
% plot(y_p2, z_p2, 'r.', 'LineWidth', 2, 'MarkerSize', 10);
% plot(yy1, zz1, 'b-', 'LineWidth', 2);
% plot(yy2, zz2, 'b-', 'LineWidth', 2);
xlabel('y [m]')
ylabel('z [m]')
zlabel('B [T]')
title('2D magnetic field')
colorbar %add colorbar
% shading flat %Removes black lines from the mesh

movegui(f1,[100 600]);
movegui(f2,[700 600]);