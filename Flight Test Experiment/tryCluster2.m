function clusters = tryCluster(x, y, z, threshold, minPoints)
    %Attempt to cluster points within a certain threshold

    %data structure of all nearest neighbors
    %data structure: [pointIndex1 pointIndex2 dist]
    %calculate number of lines from number of vertices
    numLines = sum(1:1:length(x)-1);
    fullGraph = zeros(numLines, 3);
    ii = 1;
    %figure(1)
    %hold on
    %axis([0 10 0 10])
    for j = 1:length(x) %for each point
        for k = j+1:length(x)
            dist = norm([(y(k) - y(j)), (x(k) - x(j))]);
            fullGraph(ii, :) = [j k dist];
            %plot([x(j) x(k)], [y(j) y(k)])
            ii = ii + 1;
        end
    end
    %title('All Possible Lines')

    %sort graph by distances so we can pop the shortest line each iteration
    fullGraph = sortrows(fullGraph, 3);
    %need a variable to keep track of all the different trees, which will
    %eventually be combined into one big tree
    trackTree = zeros(1,length(x));
    numFinalLines = length(x) - 1;
    usedIndex = zeros(1, numFinalLines); %keep track of indexes of lines that are kept
    lineCounter = 0; %keep track of number of current lines
    j = 1; %iterator for fullGraph
    numTreesMade = 0;
    numClusters = 1; %there will be n+1 clusters, where n is number of lines past threshold
    while (lineCounter < numFinalLines) %&& (j <= length(fullGraph(:, 1))) 
        %while there is less than n - 1 lines and there is still more of
        %fullGraph to iterate, might not need the second boolean check

        %check what tree the two points of the current line are in
        %pointTree = 0 means the point is isolated
        %pointTree variables which the same number are in the same tree
        %we only want to connect two points that are in different trees
        %sometimes we will have to merge trees
        %we also need to keep track of the number of total trees that were made
        %doesnt matter how many actual trees are actually there
        pointTree1 = fullGraph(j, 1);
        pointTree2 = fullGraph(j, 2);
        trackTree1 = trackTree(pointTree1);
        trackTree2 = trackTree(pointTree2);

        if fullGraph(j, 3) > threshold
            %if the distance is too big, skip it
            %we will now have one less line in the end
            numClusters = numClusters + 1;
            numFinalLines = numFinalLines - 1;
        elseif trackTree1 == trackTree2 && trackTree1 ~= 0
           %if they are already in the same tree (nonzero)
        elseif trackTree1 == 0 && trackTree2 == 0 
           %if both numbers are 0
           %update number of trees made
           numTreesMade = numTreesMade + 1;
           %update the same tree number to both points
           trackTree(pointTree1) = numTreesMade;
           trackTree(pointTree2) = numTreesMade;  
           lineCounter = lineCounter + 1;
           usedIndex(lineCounter) = j;
        elseif trackTree1 == 0 && trackTree2 ~= 0
            %if point 1 is by itself and point 2 is already in a tree
            %add point 1 to point 2 tree
            trackTree(pointTree1) = trackTree(pointTree2);
            lineCounter = lineCounter + 1;
            usedIndex(lineCounter) = j;
        elseif trackTree1 ~= 0 && trackTree2 == 0
            %if point 2 is by itself and point 1 is already in a tree
            %add point 2 to point 1 tree
            trackTree(pointTree2) = trackTree(pointTree1);
            lineCounter = lineCounter + 1;
            usedIndex(lineCounter) = j;
        else
           %both points are already in a tree
           %renumber all points in tree 1 to tree 2
           %find all indices that are in tree 1
           k = find(trackTree == trackTree1);
           %make them equal to tree 2's number
           trackTree(k) = trackTree2;
           lineCounter = lineCounter + 1;
           usedIndex(lineCounter) = j;
        end
        j = j + 1;
    end
    %get rid of the array indexes for the lines that were past the
    %threshold
    usedIndex = usedIndex(usedIndex > 0);
    %delete lines
    minSpanGraph = fullGraph(usedIndex, :);
%     figure(2)
%     hold on
    %axis([0 10 0 10])
%     for j = 1:length(minSpanGraph(:, 1))
%        plot(x(minSpanGraph(j, 1:2)), y(minSpanGraph(j, 1:2)), '-o')
%     end
% 
%     title('Minimum Span Graph')
    
    %NOTE: SOME VALUES IN TRACKTREE WILL BE ZERO, THIS MEANS THAT IT IS BY
    %THEMSELVES, THERE ARE NO NEIGHBORS WITHIN THE TRESHOLD DISTANCE
    
    %THIS IS ALL CODE TO ORGANIZE THE DATA
    clusters = cell(1, numClusters);
    clusterNumbers = -1 * ones(1, numClusters);
    k = 1;
    %get all the tree numbers
    for j = 1:length(trackTree)
        %find a new tree number that has not been used
        current = trackTree(j);
        if ~ismember(current, clusterNumbers) %first time seeing this number
            clusterNumbers(k) = current;
            k = k + 1;
        end
    end

    %take out the -1 indices on clusterNumbers
    [~, index] = min(clusterNumbers);
    clusterNumbers = clusterNumbers(1:index - 1);
    numClusters = length(clusterNumbers);
    %find number of points in each tree
    %get rid of trees less than a # of points threshold
    threshold = minPoints; %points
    counter = 0;
    keepClusters = -1 * ones(1, numClusters);
    for j = 1:numClusters
       current = length(trackTree(trackTree == clusterNumbers(j)));
       if current >= threshold && clusterNumbers(j) ~= 0
           counter = counter + 1;
           clusters{counter} = zeros(current, 3); %3 columns: x, y, z
           keepClusters(counter) = clusterNumbers(j);
       end

    end

    % %delete all empty cells
    newClusters = cell(1, counter);
    for j = 1:length(newClusters)
       newClusters{j} = clusters{j}; 
    end
    keepClusters = keepClusters(1:counter);
    clusters = newClusters;

    %dump all the points into clusters
    %need an array to keep track of all the indexing...
    storeIndices = ones(1, length(clusters));
    for j = 1:length(trackTree)
        currentNum = trackTree(j);
        for k = 1:length(clusters)
            if keepClusters(k) == currentNum
                clusters{k}(storeIndices(k), 1) = x(j);
                clusters{k}(storeIndices(k), 2) = y(j);
                clusters{k}(storeIndices(k), 3) = z(j);
                storeIndices(k) = storeIndices(k) + 1;
            end
        end
    end

    %DATA IS FINALLY ORGANIZED.....
    %ALL DATA IS STORED INTO THE CELL ARRAY 'clusters'
    %To call a certain cluster, just do 'clusters{j}'
    %this will give you the x,y,z vectors of each point which you can then do
    %the hough transform and other analysis on
end