function [X,Y,Z,BX,BY,BZ] = straightwiremagneticfield(x1,y1,z1,x2,y2,z2,current)

global u0 %permeability of free space is a global variable

Nx=11;
Nz=11;  
Ny=11;  
N=((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^0.5;   % length of the wire 

Xw=(x1:x2); % X-coordinates of the wire
Yw=zeros(N+1,1); % Y-coordinates of the wire
Zw=zeros(N+1,1); % Z-coordinates of the wire

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
