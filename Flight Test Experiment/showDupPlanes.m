clear all; close all; clc;

%my points
%x = [6 5 4];
%y = [8 4 2];
%z = [4 10 6];

%make a bunch of random planes
A = 4;
B = 6;
C = -3;
D = 5;
%Ax + By + Cz + D = 0
[x1,y1] = meshgrid(linspace(-5,5,5), linspace(-5,5,5));
z = (D - A .* x1 - B .* y1) ./ C
z1 = z + rand(size(z)) * 2;

A = 10;
B = -40;
C = -10;
D = 1;
%Ax + By + Cz + D = 0
[x2,y2] = meshgrid(linspace(-15, -5, 5), linspace(-15, -5, 5));
z = (D - A .* x2 - B .* y2) ./ C;
z2 = z + rand(size(z)) * 2;
surf(x1,y1,z1)
hold on
surf(x2,y2,z2)

%x = [reshape(x1,1,length(x1) * length(x1(1,:))) reshape(x2,1,length(x2)*length(x2(1,:)))];
%y = [reshape(y1,1,length(y1) * length(y1(1,:))) reshape(y2,1,length(y2)*length(y2(1,:)))];
%z = [reshape(z1,1,length(z1) * length(z1(1,:))) reshape(z2,1,length(z2)*length(z2(1,:)))];

accumSize = 50; 
minPercentVotes = 50;  
dupRhoInterval = 1;
dupThetaInterval = 1;
dupPhiInterval = 1;

x = reshape(x1,1,length(x1) * length(x1(1,:)));
y = reshape(y1,1,length(y1) * length(y1(1,:)));
z = reshape(z1,1,length(z1) * length(z1(1,:)));

[theta, phi, rho] = hough_3D(x, y, z, accumSize, minPercentVotes, dupRhoInterval, dupThetaInterval, dupPhiInterval)
%%
plotx = [-1000 -1000 1000 1000];
ploty = [-1000 1000 1000 -1000];
axis([min(x) max(x) min(y) max(y) min(z) max(z)])
for j = 1:length(rho) %for all detected planes
    points = zeros(3, length(plotx)); %store the x,y,z corner points of plane
    for k = 1:length(plotx) %for each of the four corner points of the plane
        %find corresponding z value
        plotz = (rho(j) - plotx(k).*cosd(theta(j)).*sind(phi(j)) - ploty(k).*sind(phi(j)).*sind(theta(j))) ./ cosd(phi(j));
        %pitfall is that if phi = 90 or -90, z goes to infinity
        %need to figure out a more robust version of this
        points(1,k) = plotx(k); %append x coord
        points(2,k) = ploty(k); %append y coord
        points(3,k) = plotz; %append z coord
    end
    %plot four corner points, filling the area between them
    fill3(points(1,:),points(2,:),points(3,:),'r')
    grid on
    alpha(0.3)
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    hold on

end