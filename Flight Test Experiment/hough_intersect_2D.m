function [theta, rho] = hough_intersect_2D(hough1, hough2, tol)
    %tests if two 2-D hough transform curves intersect each other
    %tol = difference between two point to classify it as an intersection
    %hough1 and hough2 must have same theta values per indices
    %returns theta value of intersection and averaged rho

    diff = [hough1(:,1) abs(hough2(:,2) - hough1(:,2))];
    [minDiff, minIndex] = min(diff(:,2));
    theta = NaN;
    rho = NaN;
    if minDiff < tol
        theta = hough1(minIndex, 1);
        rho = (hough1(minIndex, 2) + hough2(minIndex, 2)) / 2;
    end
end