%Post process data from the LiDAR.
%
%Instructions
%   1.  Upload your data to the K drive (see CL or HR)
%   2.  Modify the robocopy script located at
%       \\Mapping\Software\UW\LiDARPostProcessing\robocopy\CopyKDriveToCDrive_Data.cmd
%   3.  Run the robocopy script to copy files down to your local machine.
%   4.  Run this Matlab script
%
%Created by Christopher Lum
%lum@uw.edu
 
%Version History
%05/23/18: Created
 
clear
clc
close all
 
t0 = cputime;
 
%% Constants
%dataPath = '/Users/slammer/Desktop/PleaseWork/Mapping/Software/UW/LiDARPostProcessing';   %path to the data directory
%addpath(dataPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% RUN THIS SCRIPT IN THE LIDARPostProcessing FOLDER BECAUSE IT %%%%%%
%%%%%% CONTAINS ALL OF THE LOCAL FILES THAT YOU MIGHT NEED!!!       %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHOOSE YOUR FILE!!!
%fileName = '20180305_155956_fountain_night_dual_2.csv';
%fileName = '20180305_161614_redsquare_night_dual.csv';

%Is the file from the 5/12/18 flight test or the sample data taken on
%campus? Comment out the other one.
%fileType = 'sample';
fileType = 'flighttest';
fileName = '2018-05-12_141619 (Frame 5807).csv';
%set threshold distance for two points within the same object (meters)
threshold = 3;
minPoints = 10;
%STANDARD HOUGH TRANSFORM PARAMETERS
accumSize = 100; %This is the number of spacings of the accumulator array
minPercentVotes = 22; %This is the percentage of the number of points, 
%number of votes in accumulator cell must be greater than this to be recognized as a plane

%These next three variables are to try to remove duplicate planes (planes
%that essentially capture the same points. These variables are the number
%of linspace intervals you want to check for theta, phi, and rho values.
%For example, if your accumulator array for theta is 0,3,6,9, 3n..., and your
%dupThetaInterval is 4, then any plane that has theta values within 3*4 =
%12 degrees of each other will be classified as duplicates
dupRhoInterval = 20;
dupThetaInterval = 25;
dupPhiInterval = 25;
%Note that all of these parameters that are passed into the hough_3D
%function are pretty dependent upon the dataset... this is what sucks about
%this. Hopefully once the randomized hough transform is coded then there
%will be less parameters and they won't matter as much.
numClustersToHough = 2; %number of clusters to perform hough transform on... to not blow up memory


%% Read .csv files



if strcmp(fileType, 'flighttest')
    %this is flight test data. we need to do a coordinate rotation
    %this function also eliminates some datapoints coming from the drone
    %itself, and removes all data points under 1.5 meters, because those
    %points are mostly the ground
    disp(['Scanning flight test data...', char(10)])
    [x,y,z] = rotatePCAP(fileName);
else 
    %this is sample data given to us by AFSL, we simply read the .csv file
    disp(['Scanning AFSL sample data...', char(10)])
    data = load(fileName);
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
end
 
3%% RAW 3D DATA
figure(1)
plot3(x, y, z, '.')
xlabel('x (m)', 'Fontsize', 20)
ylabel('y (m)', 'Fontsize', 20)
zlabel('z (m)', 'Fontsize', 20)
grid on
title('Raw Data of Red Square', 'Fontsize', 20)
 
%% TAKE EVERY NTH POINT SO YOU DONT BLOW UP YOUR MEMORY

if strcmp(fileType, 'sample')
    figure(2)
    n = 3;
    deleteIndex = [];
    x = data(:, 1);
    x = x(1:n:end);
    y = data(:, 2);
    y = y(1:n:end);
    z = data(:, 3);
    z = z(1:n:end);

    %delete points near Z = 0 (these are those weird circle things around the
    %sensor that arent actually objects on the map)
    for j = 1:length(x)
        if (z(j) > -2 && z(j) < 2) %choose cutoff altitude
            deleteIndex = [deleteIndex j];
        end
    end
    x(deleteIndex) = [];
    y(deleteIndex) = [];
    z(deleteIndex) = [];
    %PLOT TOP DOWN VIEW OF DATA AFTER REMOVING Z = 0 POINTS AND TAKING EVERY
    %NTH POINT
    plot3(x, y, z, '.')
    xlabel('x (m)', 'Fontsize', 20)
    ylabel('y (m)', 'Fontsize', 20)
    zlabel('z (m)', 'Fontsize', 20)
    grid on
    %title('Top Down view after removal of unwanted points')
    title('Ground-Subtracted Data of Red Square', 'Fontsize', 20)
end

%% THIS IS THE CODE TO MAKE CLUSTERS
%set threshold length between two objects (in meters)
disp(['Clustering ', num2str(length(x)), ' points into objects based on threshold distance...'])
timeStart = cputime;
clusters = tryCluster2(x, y, z, threshold, minPoints);
disp(['Found ', num2str(length(clusters)), ' cluster(s)'])
disp(['Time required for object clustering = ',num2str(cputime - timeStart),' seconds', char(10)])
 

%% THE SAME TOP DOWN VIEW AS BEFORE, BUT NOW CLUSTERS ARE SEPARATED BY COLOR
%you can change the code to view it in 3d if u want
figure(3)
colors = jet(length(clusters));
colors = [0 0 1; 0 1 0; 1 0 0; 0.5 0 0.5; 0 0.5 0.5; 0.5 0.5 0;0 0 1; 0 1 0; 1 0 0; 0.5 0 0.5; 0 0.5 0.5; 0.5 0.5 0;0 0 1; 0 1 0; 1 0 0; 0.5 0 0.5; 0 0.5 0.5; 0.5 0.5 0; 0 0 1; 0 1 0;]
for j = 1:length(clusters)
   plot3(clusters{j}(:,1), clusters{j}(:,2), clusters{j}(:,3),'.','color', colors(j,:))
   hold on
end
xlabel('x (m)', 'Fontsize', 20)
ylabel('y (m)', 'Fontsize', 20)
zlabel('z (m)', 'Fontsize', 20)
title('Clustered Objects in Flight Test Data', 'Fontsize', 20)

%% IF YOU WANT TO SEE EACH INDIVIDUAL 3D CLUSTER. THIS WILL GENERATE HELLA PLOTS
for j = 1:length(clusters)
   figure(j + 3)
   plot3(clusters{j}(:,1), clusters{j}(:,2), clusters{j}(:, 3), '*')
   xlabel('x')
   ylabel('y')
   zlabel('z')
end


%% Perform Standard Hough Transform
%these combinations work
%accsumSize, min%Vote, dupRho, dupTheta, dupPhi
%50, 60, 20, 20, 20
%100 60 15 20 20

disp('Perfoming Hough Transform... this will take a while...')
for j = 2
    disp([char(10), 'Performing Hough Transform of cluster ', num2str(j), '...'])
    startTime = cputime;
    x = clusters{j}(:,1);
    y = clusters{j}(:,2);
    z = clusters{j}(:,3);
    xSubtract = mean(x);
    ySubtract = mean(y);
    zSubtract = mean(z);
    x = x - xSubtract;
    y = y - ySubtract;
    z = z - zSubtract;
    [theta, phi, rho] = hough_3D(x, y, z, accumSize, minPercentVotes, dupRhoInterval, dupThetaInterval, dupPhiInterval);
    
    figure
    %plot hough transform
    plot3(x,y,z,'o', 'Linewidth', 2)
    hold on
    
    %try to remove the planes that are parallel to the lidar lines 
    %this is currently hardcoded, but will eventually need to make it a
    %parameter
    deleteIndex = [];
    for j = 1:length(rho)
%         if abs(phi(j)) > 75 && ((abs(theta(j) - 85) < 10) || (abs(theta(j) - 265) < 10))
%             deleteIndex = [deleteIndex j];
%         end
    end
    rho(deleteIndex) = [];
    theta(deleteIndex) = [];
    phi(deleteIndex) = [];
    %%
    %this is all code to plot planes
    %we will make the four corner points of a plane with huge x and y bounds
    %this will give the appearance of an infinite plane MOST OF THE TIME
    plot3(x,y,z,'bo', 'Linewidth', 2)
    axis([min(x) max(x) min(y) max(y) min(z) max(z)])
    plotx = [-1000 -1000 1000 1000];
    ploty = [-1000 1000 1000 -1000];
    for j = 2:1:3 %for all detected planes
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
        axis([min(x) max(x) min(y) max(y) min(z) max(z)])
        pause(1)
        grid on
        alpha(0.3)
        xlabel('x (m)', 'Fontsize', 20)
        ylabel('y (m)', 'Fontsize', 20)
        zlabel('z (m)', 'Fontsize', 20)
        hold on
      
    end
    %title(['Cluster ', num2str(j), ': ', num2str(length(rho)), ' plane(s) found'])
    title(['Cluster ', num2str(j), ': ', num2str(2), ' plane(s) found'], 'Fontsize', 20)
    disp(['Time required to process cluster ', num2str(j), ' = ' ,num2str(cputime - startTime),' seconds'])

end

%% Clean up
%rmpath(dataPath)
 
disp(['Total time required for processing = ',num2str(cputime - t0),' seconds'])
disp('DONE!')