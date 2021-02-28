clc;clear;close all

global u0
u0=4*pi*1e-7; % Permeability of free space

position = [6 6 6];
length = 50;
current = 3;

[X,Y,Z,BX,BY,BZ] = straightwiremagneticfield(position,length,current);

figure(1)
quiver3(X,Y,Z,BX,BY,BZ,2);
hold on
line([-5 5],[0 0], [0 0],'linewidth',3,'color','r');
 axis([-5 5 -5 5 -5 5])
xlabel('X-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
zlabel('Z-axis','fontsize',14)
title('B-field of a current wire along X-axis','fontsize',14)
h=gca; 
set(h,'FontSize',14)
fh = figure(1); 
set(fh, 'color', 'white'); 

figure(2)
quiver(Y((position(1)-1)/2,:,:),Z((position(1)-1)/2,:,:),BY((position(1)-1)/2,:,:),BZ((position(1)-1)/2,:,:),2);
hold on
G1=plot(0,0,'.','markersize',6);
set(G1,'MarkerEdgeColor','r')
axis([ -5 5 -5 5])
xlabel('Y-axis','fontsize',14)
ylabel('Z-axis','fontsize',14)
title('B-field YZ plane','fontsize',14)
h=gca; 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(2); 
set(fh, 'color', 'white'); 