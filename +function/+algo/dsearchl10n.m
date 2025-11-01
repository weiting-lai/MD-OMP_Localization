function [k,distance,distance_xy,distance_z] = dsearchl10n(source,estimated_pt)
% Function  : dsearchl10n
% Creation  : 21-07-24
% Author    : Weiting Lai
% Version   : 21-07-24
%
% Description
%
%   Find the nearest data point to each query point, and compute the 
%   corresponding distances.
%
%   Function won't re-select either the same data point or query point.
%
% Inputs
%
%   source          True source points.
%   estimated_pt    Estimated source points.
%
% Outputs
%
%   k               Indices of the data points closest to the query points.
%   distance        Distances between each query point and the closest 
%                   input point.
%
dist = zeros(size(source,1));
for i = 1:size(source,1)
    [~, d] = dsearchn(source(i,:),estimated_pt);
    dist(i,:) = d; 
end
idx_pt = zeros(size(dist,2),2);
distance = zeros(size(dist,2),1);
distance_xy = zeros(size(dist,2),1);
distance_z = zeros(size(dist,2),1);
for i = 1:size(source,1)
    [a,idx] = min(dist(:));
    distance(i) = a;
    idx_pt(i) = idx;
    [row, col] = ind2sub(size(dist), idx);
    dist(row,:) = Inf;
    dist(:,col) = Inf;
    idx_pt(i,1) = row;
    idx_pt(i,2) = col;
    distance_xy(i) = norm(source(row,1:2)-estimated_pt(col,1:2));
    distance_z(i) = norm(source(row,3)-estimated_pt(col,3));
end
[~,pt_sort] = sort(idx_pt(:,1));
k = idx_pt(pt_sort, 2);
distance = distance(pt_sort);
distance_xy = distance_xy(pt_sort);
distance_z = distance_z(pt_sort);
end