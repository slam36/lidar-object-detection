clear all; close all; 

pointA = [-50, -50, -79.3];
pointB = [50, 50, 58.76];
pointC = [0, 0, -10.26];

points=[pointA' pointB' pointC']; % using the data given in the question
fill3(points(1,:),points(2,:),points(3,:),'r-o')
grid on
alpha(0.3)
xlabel('x')
ylabel('y')
zlabel('z')