function [Bx,By,Bz] = magnetic_field_current_loop(x,y,z,x_p,y_p,z_p,a,I0)

global u0 %permeability of free space is a global variable

rc=((x-x_p).^2+(y-y_p).^2).^.5; %Radial component is required for cylindrical coordinate system.
m=(4.*a.*rc).*(((rc+a).^2)+((z-z_p).^2)).^(-1); %This is a parameter for calculating the Elliptical integrals
kofkc=(pi/2)+(pi/8).*m+(9*pi/128).*m.^2; %K(k) elliptical function, this is a taylor expansion of the K elliptical integral.
eofkc=(pi/2)+(-pi/8).*m+(-3*pi/128).*m.^2;%E(k) elliptical function this is a taylor expansion of the E elliptical integral.

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