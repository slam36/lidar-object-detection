# LiDAR Man-Made Object Detection in Nature
This mini-project was the first stepping stone towards developing a real-time aerial LiDAR flagging software that detects the presence of man-made objects (as opposed to objects of nature) in LiDAR point cloud data. This could be used towards multiple applications such as alerting pilots of nearby poachers in wildlife environments. Note that the purpose of the project is not to identify these objects, but to simply flag them to be verified by a human operator. This project was done primarily in Matlab, as it is a good prototying scripting language for the research stage of this project. We were required to implement all code from scratch so no external matlab libraries were used.

To gather data, a LiDAR sensor was mounted on a remotely piloted octocopter. We performed flight tests at the University of Washington Carnation Farms UAS Test Site, which is mostly an empty patch of farmland. We placed cars and tents and large boxes on the ground amongst trees and bushes and flew around the farm and collected LiDAR point cloud data.

Our flagging algorithm is based on the assumption that man-made objects contain flat plane surfaces as opposed to nature. For example, the windshield of a car is flat. It is rare for objects in nature to have this characteristic. Trees and bushes (unless they are trimmed) have very complex shapes. 

### Concise Explanation of End to End Algorithm 

First, the objects in the point cloud frame need to be clustered so the detection algorithm can be run on each individual object. Clusters are defined as a group of points where there are always two points in the cluster that are less than a threshold distance from each other. 

Once the point clouds for the individual objects are isolated, we can start the point cloud analysis. The 2D Hough transform is a feature extraction used to detect lines and circles in images or two dimensional point clouds. This algorithm can be extended to the third dimension to detect planes in three dimensional point clouds. 

We faced a few challenges with this project. First, we had to define a good threshold for the minimum number of points that lie near a plane in order to be considered a valid plane. Without a threshold, any group of three points could be identified as a plane. Even with the threshold, the algorithm would still produce duplicate planes, which were two or more planes that are very similar to each other in terms of location and orientation. 

The flagging criteria has not been fully developed, but it is proposed that if more than 75% of the points in the dataset are within a certain distance of the detected planes, the object can be flagged as man-made.


### Code and File Structure
All of the code is in the /Flight Test Experiment folder. Run MAIN_LiDAR_Post_Processing.m. The code is split into sections for each step of the processing and algorithm. A sample csv file containing a frame of point cloud data from the test flight is located in this directory too which the code currently reads from. 

The /Flight Planning directory contains the Mission Planner waypoints for our test flight. 

The report that I wrote for this project is also in the root folder. This describes the project and all of its steps in detail.





