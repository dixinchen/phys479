function [Bx,By,Bz] = magnetic_field_current_loop(x,y,z,xc,yc,zc,n,a,I0)

global u0 %permeability of free space is a global variable

% (x,y,z) is the coordinate of an arbitrary point in the space
% (xc,yc,zc) is the coordinate of the center of the current loop
% n is the normal vector of the center of the current loop
% a is the radius of the loop
% I0 is the current flowing through the loop

% Rotate the loop so that it aligns with z-axis
r = vrrotvec(n,[0,0,1]); % the rotation needed to transform the normal to z
[nx,ny,nz] = r.*[xn,yn,zn];
[cx,cy,cz] = r.*[xc,yc,zc];
% Change from cartesian to cylindrical coord.
[phi,rc,zz] = cart2pol((x-cx),(y-cy),(z-cz));

% Calculate k_c^2
kc2 = (4.*a.*rc)/((rc+a).^2+(z-cz).^2)

% Calculate the elliptic integral of the first and second kind
[kofk,eofk] = ellipke()

kofkc=(pi/2)+(pi/8).*m+(9*pi/128).*m.^2; %K(k) elliptical function, this is a taylor expansion of the K elliptical integral.
eofkc=(pi/2)+(-pi/8).*m+(-3*pi/128).*m.^2;%E(k) elliptical function this is a taylor expansion of the E elliptical integral.
eofkc=(pi/2)+(3pi/8).*m+(45*pi/128).*m.^2;%E(k) elliptical function this is a taylor expansion of the E elliptical integral.

%Note for improved accuracy, Matlab has built in elliptical integral
%calculation but these expressions here are still very accurate when rc < a

Brc=(u0.*I0./(2.*pi.*rc)).*(z-z_p).*((((rc+a).^2)+((z-z_p).^2)).^(-.5)).*(-kofkc+eofkc.*((rc.^2+a.^2+(z-z_p).^2)./(((rc-a).^2)+((z-z_p).^2)))); %radial component of B%
Bz=(u0.*I0./(2.*pi)).*((((rc+a).^2)+((z-z_p).^2)).^(-.5)).*(kofkc-eofkc.*((rc.^2-a.^2+(z-z_p).^2)./(((rc-a).^2)+((z-z_p).^2)))); %axial component of B
Bx=Brc.*(x-x_p)./rc; %This converts the polar component into cartesian form.
By=Brc.*(y-y_p)./rc;

%The following sets any terms that result in Inf to zero, this occurs at
%the points near the coil itself.
Bx(isnan(Bx)) = 0 ;
By(isnan(By)) = 0 ;
Bz(isnan(Bz)) = 0 ;