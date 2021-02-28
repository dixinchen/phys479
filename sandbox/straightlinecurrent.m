% function [BX,BY,BZ] = straightlinecurrent(x1,y1,z1,x2,y2,z2,current)
clc;clear;close all
x1 = 0;y1 = 0;z1 = 0;
x2 = 3;y2 = 3;z2 = 4;

current = -7;
u0=4*pi*1e-7;

Nx=11;
Nz=11;  
Ny=11;  
N=50;

Xw=linspace(x1,x2,N+1); % X-coordinates of the wire
Yw=linspace(y1,y2,N+1).'; % Y-coordinates of the wire
Zw=linspace(z1,z2,N+1).'; % Z-coordinates of the wire

I=current;    % current 

xp(1:Nx)=-(Nx-1)/2:(Nx-1)/2;  
yp(1:Ny)=-(Ny-1)/2:(Ny-1)/2;    
zp(1:Nz)=-(Nz-1)/2:(Nz-1)/2;  

X(1:Nx,1:Ny,1:Nz)=0;
Y(1:Nx,1:Ny,1:Nz)=0; 
Z(1:Nx,1:Ny,1:Nz)=0;

for i=1:Nx
    X(i,:,:)=xp(i); 
end
for i=1:Ny
    Y(:,i,:)=yp(i); 
end
for i=1:Nz
    Z(:,:,i)=zp(i);
end

for a=1:Nx  
for b=1:Ny  
for c=1:Nz

for i=1:N-1
Rx(i)=X(a,b,c)-(0.5*(Xw(i)+Xw(i+1)));
Ry(i)=(Y(a,b,c)-(0.5*(Yw(i)+Yw(i+1))));
Rz(i)=(Z(a,b,c)-(0.5*(Zw(i)+Zw(i+1))));
dlx(i)=Xw(i+1)-Xw(i);
dly(i)=Yw(i+1)-Yw(i);
dlz(i)=Zw(i+1)-Zw(i);
end
Rx(N)=(X(a,b,c)-0.5*(Xw(N)+1));
Ry(N)=(Y(a,b,c)-(0.5*(Yw(N)+1)));
Rz(N)=(Z(a,b,c)-(0.5*(Zw(N)+1)));
dlx(N)=-Xw(N)+1;
dly(N)=-Yw(N)+1;
dlz(N)=-Zw(N)+1;

for i=1:N
Xcross(i)=(dly(i).*Rz(i))-(dlz(i).*Ry(i));
Ycross(i)=(dlz(i).*Rx(i))-(dlx(i).*Rz(i));
Zcross(i)=(dlx(i).*Ry(i))-(dly(i).*Rx(i));
R(i)=sqrt(Rx(i).^2+Ry(i).^2+Rz(i).^2);
end

Bx1=(I*u0./(4*pi*(R.^3))).*Xcross;
By1=(I*u0./(4*pi*(R.^3))).*Ycross;
Bz1=(I*u0./(4*pi*(R.^3))).*Zcross;

BX(a,b,c)=0;       % Initialize sum magnetic field to be zero first
BY(a,b,c)=0;
BZ(a,b,c)=0;

for i=1:N   % loop over all current elements along coil    
    BX(a,b,c)=BX(a,b,c)+Bx1(i);
    BY(a,b,c)=BY(a,b,c)+By1(i);
    BZ(a,b,c)=BZ(a,b,c)+Bz1(i);
end

end
end

end


figure(1)

quiver3(X,Y,Z,BX,BY,BZ,2);
hold on
line([x1 x2],[y1 y2], [z1 z2],'linewidth',3,'color','r');
xlabel('X-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
zlabel('Z-axis','fontsize',14)
title('B-field of a current wire along X-axis','fontsize',14)
h=gca; 
set(h,'FontSize',14)
fh = figure(1); 
set(fh, 'color', 'white'); 

% figure(2)
% quiver(Y((Nx-1)/2,:,:),Z((Nx-1)/2,:,:),BY((Nx-1)/2,:,:),BZ((Nx-1)/2,:,:),2);
% hold on
% G1=plot(y1,z1,'.','markersize',6);
% set(G1,'MarkerEdgeColor','r')
% xlabel('Y-axis','fontsize',14)
% ylabel('Z-axis','fontsize',14)
% title('B-field YZ plane','fontsize',14)
% h=gca; 
% set(h,'FontSize',14)
% h = get(gca, 'ylabel');
% fh = figure(2); 
% set(fh, 'color', 'white'); 