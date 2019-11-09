clear all; close all; clc;

%make a bunch of random planes
A = 4;
B = 6;
C = -3;
D = 5;
%Ax + By + Cz + D = 0
[x1,y1] = meshgrid(-5:2:5, -5:2:5);
z = (D - A .* x1 - B .* y1) ./ C;
z1 = z + rand(size(z)) * 2;

A = 10;
B = -40;
C = -10;
D = 1;
%Ax + By + Cz + D = 0
[x2,y2] = meshgrid(-15:2:-5, -15:2:-5);
z = (D - A .* x2 - B .* y2) ./ C;
z2 = z + rand(size(z)) * 2;
surf(x1,y1,z1)
hold on
surf(x2,y2,z2)

x = [reshape(x1,1,length(x1) * length(x1(1,:))) reshape(x2,1,length(x2)*length(x2(1,:)))];
y = [reshape(y1,1,length(y1) * length(y1(1,:))) reshape(y2,1,length(y2)*length(y2(1,:)))];
z = [reshape(z1,1,length(z1) * length(z1(1,:))) reshape(z2,1,length(z2)*length(z2(1,:)))];

% x = reshape(x1,1,length(x1) * length(x1(1,:)));
% y = reshape(y1,1,length(y1) * length(y1(1,:)));
% z = reshape(z1,1,length(z1) * length(z1(1,:)));

%
% 1: while still enough points in point set P do
% 2:    randomly pick three points p1, p2, p3 from the set of
%       points P
% 3:    if p1, p2 and p3 fulfill the distance criterion then
% 4:        calculate plane (?, ?, ?) spanned by p1, p2, p3
% 5:        increment cell A(?, ?, ?) in the accumulator space
% 6:            if the counter |A(?, ?, ?)| equals threshold t then
% 7:                (?, ?, ?) parameterize the detected plane
% 8:                delete all points close to (?, ?, ?) from P
% 9:                reset the accumulator
% 10:           end if
% 11:   else
% 12:       continue
% 13:   end if
% 14: end while