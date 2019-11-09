function [theta, phi, rho] = hough_3D(x, y, z, accumSize, minPercentVotes, dupRhoInterval, dupThetaInterval, dupPhiInterval)

% %my points
% %x = [6 5 4];
% %y = [8 4 2];
% %z = [4 10 6];
% 
% %make a bunch of random planes
% A = 4;
% B = 6;
% C = -3;
% D = 5;
% %Ax + By + Cz + D = 0
% [x1,y1] = meshgrid(-5:2:5, -5:2:5);
% z = (D - A .* x1 - B .* y1) ./ C;
% z1 = z + rand(size(z)) * 2;
% 
% A = 10;
% B = -40;
% C = -10;
% D = 1;
% %Ax + By + Cz + D = 0
% [x2,y2] = meshgrid(-15:2:-5, -15:2:-5);
% z = (D - A .* x2 - B .* y2) ./ C;
% z2 = z + rand(size(z)) * 2;
% surf(x1,y1,z1)
% hold on
% surf(x2,y2,z2)
% 
% x = [reshape(x1,1,length(x1) * length(x1(1,:))) reshape(x2,1,length(x2)*length(x2(1,:)))];
% y = [reshape(y1,1,length(y1) * length(y1(1,:))) reshape(y2,1,length(y2)*length(y2(1,:)))];
% z = [reshape(z1,1,length(z1) * length(z1(1,:))) reshape(z2,1,length(z2)*length(z2(1,:)))];
% 
% % x = reshape(x1,1,length(x1) * length(x1(1,:)));
% % y = reshape(y1,1,length(y1) * length(y1(1,:)));
% % z = reshape(z1,1,length(z1) * length(z1(1,:)));
% fileName = 'Frame5807_tree.csv';
%  
% % Load data
% data = load(fileName);
% x = data(:,1);
% y=data(:,2);
% z=data(:,3);

%%
%try to make the points as close as (0,0,0) as possible to decrease rho
%size
% xSubtract = mean(x);
% ySubtract = mean(y);
% zSubtract = mean(z);
% x = x - xSubtract;
% y = y - ySubtract;
% z = z - zSubtract;

allRhos = sqrt(x.^2+y.^2+z.^2); %calculate vector magnitudes to determine the maximum values of rho
%so if it was 36, your theta values would be 0, 10, 20,... 350, 360
maxRho = max(allRhos);
rhoRange = linspace(-maxRho, maxRho, accumSize);
thetaRange = linspace(0, 360, accumSize);
%thetaRange = linspace(0, 180, accumSize);
phiRange = linspace(0, 90, accumSize);
%phiRange = linspace(-90, 90, accumSize);
%this is a three dimensional matrix with each dimension of size accumSize
accumulator = zeros(length(rhoRange), length(rhoRange), length(rhoRange)); %(theta, phi, rho)

%the accumulator matrix contains a lot of possible values of theta,phi, and
%rho. for each point, we will iterate through all of the possible values of
%theta and phi, find the corresponding rho value for that theta and phi,
%and then find the nearest rho value in the accumulator array and add a
%vote to that cell
disp('Filling accumulator array...')
for j = 1:length(x) %iterate through points
    for k = 1:length(thetaRange) %iterate through values of theta
       for l = 1:length(phiRange) %iterate through values of phi
           %calculate the rho value at this x,y,z,phi,theta
           currentRho = x(j) * cosd(thetaRange(k)) * sind(phiRange(l)) + y(j) * sind(thetaRange(k)) * sind(phiRange(l)) + z(j) * cosd(phiRange(l));
           %find the nearest value of the current rho to the rho values on
           %the accumulator array
           %for example, if your accumulator array rho values are 0:0.1:2,
           %and your rho value is 1.67, it will round it to 1.7
           interpRho = interp1(rhoRange, rhoRange, currentRho, 'nearest');
           %find the index that has that rho value in the accumulator
           index = find(rhoRange==interpRho);
           %increment that index (add a vote)
           accumulator(k, l, index) = accumulator(k, l, index) + 1;
       end
    end
end


%set the threshold number of votes for a certain theta,phi,rho value
%currently i just have it set to whatever the max number of votes for a
%certain accumulator cell was
%thres = max(max(max(accumulator))) * minPercentVotes / 100
%let's try making it a percentage of the total number of points 

thres = minPercentVotes * length(x) / 100;
%now we need to store the theta,rho,phi values that represent "planes"
%we will use vector append for this, even tho it is not efficient 
%maybe there is a way to do it using logical indexing
rho = zeros(size(x))';
theta = zeros(size(x))';
phi = zeros(size(x))';

counter = 1;
for j = 1:length(accumulator) %iterate theta values
    for k = 1:length(accumulator) %iterate phi values
        for l = 1:length(accumulator) %iterate rho values
            %check the number of votes for current theta,phi,rho are
            %greater than or equal to the threshold
            if accumulator(j,k,l) >= thres
                %append theta, phi, and rho value 
                theta(counter) = thetaRange(j);
                phi(counter) = phiRange(k);
                rho(counter) = rhoRange(l);
                counter = counter + 1;
            end
        end
    end
end
theta(counter:end) = [];
phi(counter:end) = [];
rho(counter:end) = [];

%%%%%%%
%%
% try to capture duplicates
%planes that have similar rho, theta, phi values are probably repeats
%set thresholds for eliminating repeat planes
rhoThres = maxRho / (accumSize / dupRhoInterval);
thetaThres = 180 / (accumSize / dupRhoInterval);
phiThres = 90 / (accumSize / dupRhoInterval);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%daniel's solution
%if two out of the three plane parameters are close to each other, consider
%it a duplicate
%do this a couple times...
    for repeat = 1:5
        for j = 1:length(theta) %iterate through all possible planes
            %check all of the remaining planes
            deleteIndex = [];
            for k = j+1:length(theta)
                counter = 0;
                if abs(theta(j) - theta(k)) <= thetaThres
                    counter = counter + 1;
                end
                if abs(phi(j) - phi(k)) <= phiThres
                    counter = counter + 1;
                end
                if abs(rho(j) - rho(k)) <= rhoThres
                    counter = counter + 1;
                end
                if counter > 2
                   deleteIndex = [deleteIndex k];
                end
            end

            %take average of all indexes you are going to delete
            if length(deleteIndex) > 0
               theta(j) = mean([theta(j) theta(deleteIndex)]);
               rho(j) = mean([rho(j) rho(deleteIndex)]);
               phi(j) = mean([phi(j) phi(deleteIndex)]);
            end
            theta(deleteIndex) = [];
            phi(deleteIndex) = [];
            rho(deleteIndex) = [];
        end
    end
    
    


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ryan and i's solution
% %initialize empty arrays of theta,phi,rho values that we will append true
% %planes to. this is horribly inefficient, will edit later if have time
% %this is kind of like a fencepost problem, so we start the loop at
% %iteration 2 and append the first plane by itself
% trueTheta = zeros(size(theta));
% truePhi = zeros(size(phi));
% trueRho = zeros(size(rho));
% trueTheta(1) = theta(1);
% truePhi(1) = phi(1);
% trueRho(1) = rho(1);
% tic
% counter = 2;
% for j = 2:length(theta) %for all detected planes
%     %note that theta is sorted!!! this is very convenient
%     %this means trueTheta will also be sorted
%     %iterate through all trueTheta values that are within threshold
%     k = length(trueTheta);
%     tic
%     while (k > 0 && theta(j) - trueTheta(k) <= thetaThres)
%         %check rho and phi values
%         if abs(rho(j) - trueRho(k)) < rhoThres
%             break
%         end
%         if abs(phi(j) - truePhi(k)) < phiThres
%             break
%         end
%         %this is a match! append to the true planes arrays
%         trueTheta(counter) = theta(j);
%         truePhi(counter) = phi(j);
%         trueRho(counter) = rho(j);
%         counter = counter + 1;
%         k = k - 1;
%     end
%     toc
%     
% end
% toc

%trueTheta(counter + 1 : end) = [];
%truePhi(counter + 1 : end) = [];
%trueRho(counter + 1 : end) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end






