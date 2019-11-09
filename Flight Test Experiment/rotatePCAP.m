function [X,Y,Z] = rotatePCAP(fileName)

    f = readtable(fileName);
    X = f{:,1};
    Y = f{:,2};
    Z = f{:,3};

    xregress = [ones(length(X),1) X];
    b = xregress \ Y;

    theta = -atand(b(2));

    x = X.*cosd(theta) - Y.*sind(theta);
    y = X.*sind(theta) + Y.*cosd(theta);
    z = Z;



    theta = -90;
    X = x;
    Y = y.*cosd(theta) - z.*sind(theta);
    Z = y.*sind(theta) + z.*cosd(theta);


    minZ = min(Z);
    Z = Z - minZ;

    maxZ = max(Z);
    deleteIndex = [];
    for j = 1:length(Z)
        if Z(j) > maxZ - 10 %there were some sparse points from drone itself
            deleteIndex = [deleteIndex j];
        end
        if Z(j) < 1.5 %only capture points above altitude of 1.5 m 
            deleteIndex = [deleteIndex j];
        end
    end

    X(deleteIndex) = [];
    Y(deleteIndex) = [];
    Z(deleteIndex) = [];

    % plot3(X,Y,Z,'.')
    % xlabel('x')
    % ylabel('y')
    % zlabel('z')

    %csvwrite('Frame5807Reduced.csv', [X Y Z])

    %axis([-50 -0 0 10 1.5 3])
end





