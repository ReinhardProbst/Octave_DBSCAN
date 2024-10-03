##
## Originally based on ChatGPT code creation
##

# Sample points
points = [ 1.0,  1.1;
           1.2,  1.3;
           0.9,  0.8;
          10.0, 10.1;
          10.2, 10.3;
           9.9,  9.8;
           5.0,  5.1];

# The DBSCAN parameter
eps = 1.5;
minPts = 3;

fig_no = 0;

function [clusterIds, visited] = expand_cluster(points, clusterIds, pointIdx, neighbors, clusterId, eps, minPts, visited)
    clusterIds(pointIdx) = clusterId;
    seeds = neighbors;

    idx = 1;
    while idx <= length(seeds)
        currentP = seeds(idx);
        if ~visited(currentP)
            visited(currentP) = 1;
            result = region_query(points, points(currentP, :), eps);
            if length(result) >= minPts
                #printf("Found core point with idx: %d x: %.2f y: %.2f\n", currentP, points(currentP,1), points(currentP,2));
                seeds = [seeds; result];       # Add new neighbors to list
            end
        end
        if clusterIds(currentP) == -1
            clusterIds(currentP) = clusterId;  # Set to cluster if not yet explored
        end
        idx += 1;
    end
end

function neighbors = region_query(points, p, eps)
    numPoints = size(points, 1);
    neighbors = [];
    for i = 1:numPoints
        if distance_man(points(i, :), p) <= eps  # Euclidean or Manhattan distance
            neighbors = [neighbors; i];
        end
    end
end

function d = distance(a, b)
    d = sqrt((a(1) - b(1))^2 + (a(2) - b(2))^2);  # Euclidean distance
end

function d = distance_man(a, b)
    d = abs(a(1) - b(1)) + abs(a(2) - b(2));      # Manhattan distance
end

function clusterIds = dbscan_fct(points, eps, minPts)
    clusterId = 0;
    numPoints = size(points, 1);
    clusterIds = -ones(numPoints, 1);  # -1 meeans unexplored, 0 means noice, 1.. means clusterID
    visited = zeros(numPoints, 1);     # 0 means not yet visited, 1 means already visited

    for i = 1:numPoints
        if ~visited(i)
            visited(i) = 1;
            neighbors = region_query(points, points(i, :), eps);
            if length(neighbors) < minPts
                clusterIds(i) = 0;  # Mark as noice
            else
                clusterId += 1;
                [clusterIds, visited] = expand_cluster(points, clusterIds, i, neighbors, clusterId, eps, minPts, visited);
            end
        end
    end

#    for i = 1:numPoints
#        fprintf('Point (%.2f, %.2f) - Cluster ID: %d Visited: %d\n', points(i, 1), points(i, 2),
#                clusterIds(i), visited(i));
#    end
end

points = dlmread("./smile_face.csv", ";");
eps = 2;
minPts = 5;

figure(++fig_no);
plot(points(:,1), points(:,2), "kx", "markersize", 15);
grid on;
title("Raw points");

# Run DBSCAN algorithm
clusterIds = dbscan_fct(points, eps, minPts);

# Show detailed results
for i = 0:max(clusterIds)
    c = find(clusterIds==i);

    printf("Points in cluster %d\n", i);
    pts = points(c,:);
    display(pts);

    figure(++fig_no);
    plot(pts(:,1), pts(:,2), "rx", "markersize", 15);
    grid on;
    title(strcat("Points of cluster:", num2str(i)));
end
