clear all; close all; clc;

theta = linspace(0, 360, 10);
phi = linspace(-90, 90, 10);
rho = linspace(-10, 10, 10);
neg10 = ones(10) * -10;
surf(theta, phi, neg10)
xlabel('\theta', 'Fontsize', 20)
ylabel('\phi', 'Fontsize', 20)
zlabel('\rho', 'Fontsize', 20)
pos10 = ones(10) * 10;
hold on
surf(theta, phi, pos10)
title('Sample Accumulator Matrix', 'Fontsize', 20)
