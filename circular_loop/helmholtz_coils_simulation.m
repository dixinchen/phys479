clc;clear;close all

%% Define constants and parameters

% Define global variables
global u0

u0=4*pi*1e-7; % Permeability of free space

% Define coil parameters
I0=40; % Coil current in Amps
a=.4; % Coil radius in m

x_p1=0; y_p1=0; z_p1=-a/2; % Define coordinates of coil 1 center point
x_p2=0; y_p2=0; z_p2=a/2; % Define coordinates of coil 2 center point

%% Plot 3D helmholtz coils

[x1,y1,z1] = plot_helmholtz_coils(a,x_p1,y_p1,z_p1); % Calculate coordinates of coil 1
[x2,y2,z2] = plot_helmholtz_coils(a,x_p2,y_p2,z_p2); % Calculate coordinates of coil 2

f1 = figure;
% plot3(y_p1, z_p1, x_p1, 'r.', 'LineWidth', 2, 'MarkerSize', 10); % Plot the center point of coil 1
% hold on
% plot3(y_p2, z_p2, x_p2, 'r.', 'LineWidth', 2, 'MarkerSize', 10); % Plot the center point of coil 2
plot3(y1, z1, x1, 'k-', 'LineWidth', 2); % Plot coil 1
hold on
plot3(y2, z2, x2, 'k-', 'LineWidth', 2); % Plot coil 2 on x direction
% plot3(x1, y1, z1, 'k--', 'LineWidth', 1);
% plot3(x2, y2, z2, 'k--', 'LineWidth', 1);
% plot3(z1, x1, y1, 'k--', 'LineWidth', 1);
% plot3(z2, x2, y2, 'k--', 'LineWidth', 1);
grid on
axis('square');
xlabel('x [m]', 'FontSize', 20);
ylabel('y [m]', 'FontSize', 20);
zlabel('z [m]', 'FontSize', 20);

%% Calculate magnetic field

% Input mesh of points in 2D plane
x=0;
% This is a 2d plane over the x=0 plane that extends away from the coils in the yz plane. 
[y,z]=meshgrid(linspace(-a/2,a/2,25),linspace(-a/2,a/2,25)); 

[Bx1,By1,Bz1] = magnetic_field_current_loop(x,y,z,x_p1,y_p1,z_p1,a,I0); % Calculate field from first coil
[Bx2,By2,Bz2] = magnetic_field_current_loop(x,y,z,x_p2,y_p2,z_p2,a,I0); % Calculate field from second coil

% Add the components
Bx=Bx1+Bx2;
By=By1+By2;
Bz=Bz1+Bz2;
% Bx=Bx1+Bx2+By1+By2+Bz1+Bz2;
% By=By1+By2+Bz1+Bz2+Bx1+Bx2;
% Bz=Bx1+Bx2+By1+By2+Bz1+Bz2;

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
xlabel('x [m]','FontSize', 20)
ylabel('y [m]','FontSize', 20)
zlabel('B [T]','FontSize', 20)
title('2D magnetic field','FontSize', 20)
colorbar %add colorbar
% shading flat %Removes black lines from the mesh

movegui(f1,[100 600]);
movegui(f2,[700 600]);


%%
% This is a 2d plane over the x=0 plane that extends away from the coils in the yz plane. 
[x,y,z]=meshgrid(linspace(-a,a,50),linspace(-a,a,50),linspace(-a,a,50)); 

[Bx1,By1,Bz1] = magnetic_field_current_loop(x,y,z,x_p1,y_p1,z_p1,a,I0); % Calculate field from first coil
[Bx2,By2,Bz2] = magnetic_field_current_loop(x,y,z,x_p2,y_p2,z_p2,a,I0); % Calculate field from second coil

% Add the components
Bx=Bx1+Bx2+By1+By2+Bz1+Bz2;
By=By1+By2+Bz1+Bz2+Bx1+Bx2;
Bz=Bx1+Bx2+By1+By2+Bz1+Bz2;

f3 = figure;
quiver3(x,y,z,Bx,By,Bz);
hold on
plot3(y1, z1, x1, 'k-', 'LineWidth', 1);
plot3(y2, z2, x2, 'k-', 'LineWidth', 1);
plot3(x1, y1, z1, 'k-', 'LineWidth', 1);
plot3(x2, y2, z2, 'k-', 'LineWidth', 1);
plot3(z1, x1, y1, 'k-', 'LineWidth', 1);
plot3(z2, x2, y2, 'k-', 'LineWidth', 1);
grid on
axis("square");
xlabel('x [m]', 'FontSize', 20);
ylabel('y [m]', 'FontSize', 20);
zlabel('z [m]', 'FontSize', 20);