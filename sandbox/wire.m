clc;clear;close all

% Permeability of free space
mu0=4*pi*1e-7;

% Define the coordinates of the two ends of the wire
v1 = [1,2,3];
v2 = [4,5,9];
length = ((v1(1)-v2(1))^2+(v1(2)-v2(2))^2+(v1(3)-v2(3))^2)^0.5;
% Define current
current = 7;

% Define a volume of interests
xc = (v1(1)+v2(1))/2; % xyz center point
yc = (v1(2)+v2(2))/2;
zc = (v1(3)+v2(3))/2;
n = 10; % number of grids
[xi,yi,zi] = meshgrid(linspace((xc-length/2),(xc+length/2),n),linspace((yc-length/2),(yc+length/2),n),linspace((zc-length/2),(zc+length/2),n));
[bx,by,bz] = meshgrid(linspace((xc-length/2),(xc+length/2),n),linspace((yc-length/2),(yc+length/2),n),linspace((zc-length/2),(zc+length/2),n));

% Calculate spherical coord. for each point in the subspace
% for i=1:n
%     rho(i) = norm(cross(v1-v2,[xi(i),yi(i),zi(i)]-v2)) / norm(v1-v2);
%     cosTheta1(i)=(z1);
% end
% for i=1:n
%     for j=1:n
%         for k=1:n
%             rho = norm(cross(v1-v2,[xi(i,j,k),yi(i,j,k),zi(i,j,k)]-v2)) / norm(v1-v2);
%             cosTheta1 = (zi(i,j,k)-v1(3))/((rho^2+(zi(i,j,k)-v1(3))^2)^0.5);
%             cosTheta2 = (zi(i,j,k)-v2(3))/((rho^2+(zi(i,j,k)-v2(3))^2)^0.5);
% %             [bx(i,j,k),by(i,j,k),bz(i,j,k)] = sph2cart((mu0*current/(4*pi*rho))*(cosTheta2-cosTheta1),0,rho);
%             bx(i,j,k) = -sin((mu0*current/(4*pi*rho))*(cosTheta2-cosTheta1));
%             bz(i,j,k) = cos((mu0*current/(4*pi*rho))*(cosTheta2-cosTheta1));
%             by(i,j,k) = rho;
%         end
%     end
% end

c = mu0*current/(4*pi);
xA = v1(1);
yA = v1(2);
zA = v1(3);
xB = v2(1);
yB = v2(2);
zB = v2(3);
b = zeros(n,n,n);
for i=1:n
    for j=1:n
        for k=1:n
            x = xi(i,j,k);
            y = yi(i,j,k);
            z = zi(i,j,k);
            r1 = ((x-xA)^2+(y-yA)^2+(z-zA)^2)^0.5;
            r2 = ((x-xB)^2+(y-yB)^2+(z-zB)^2)^0.5;
            cosTheta1 = (r2^2-r1^2-length^2)/(2*length*r1);
            cosTheta2 = (r2^2-r1^2-length^2)/(2*length*r2);
            distance = ((2*r1^2*r2^2+2*r1^2*length^2+2*r2^2*length^2-r1^4-r2^4-length^4)^0.5)/(2*length);
            b(i,j,k) = c*(cosTheta2-cosTheta1)/distance;
        end
    end
end


% quiver3(xi,yi,zi,bx,by,bz,2)
surf(squeeze(yi(:,1,:)),squeeze(zi(:,1,:)),squeeze(b(1,:,:)));
% hold on
% line([v1(1) v2(1)],[v1(2) v2(2)],[v1(3) v2(3)],'linewidth',3,'color','r');
axis equal



