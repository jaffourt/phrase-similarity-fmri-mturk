function sphereCoords = searchlight_sphereInds(relativeCoords, gmCoords, centerCoords, volSize)

%%% Finds the *absolute* voxel indices that constitute a sphere around a
%%% given center location

%%% INPUT:
%%% relativeCoords = output from searchlight_makeSphere.m
%%% gmCoords = can be one of two types:
%%%            (1) 3D binary matrix, with 1 for gray-matter voxels
%%%            (2) Nx3 array of x-y-z coordinates of gray-matter voxels
%%% centerCoords = Mx3 array of x-y-z coordinates for center of searchlight
%%%                sphere
%%% volSize = 1x3 array of the size of the brain volume (in voxels)

%%% OUTPUT:
%%% sphereCoords = Mx1 cell, each entry contains a Kx3 array with x-y-z coordinates of K gray
%%%                matter voxels around the corresponding center voxel in centerCoords

%%% Idan Blank, Feb 20, 2017; copied from Kriegeskorte's RSA toolbox
%%%                           (with changes to variable names)

%%% SANITY CHECKED %%%

if ndims(gmCoords) == 2
    gmInds = sub2ind(volSize, gmCoords(:,1), gmCoords(:,2), gmCoords(:,3));
else
    gmInds = find(gmCoords > 0);
end

sphereCoords = cell(size(centerCoords,1),1);
for ii = 1:size(centerCoords,1)
    currVoxelCoords = repmat(centerCoords(ii,:),[size(relativeCoords,1), 1])+relativeCoords; 
        % x,y,z coordinates of voxels around the center that are illuminated by sphere

    badVoxels=(currVoxelCoords(:,1)<1 | currVoxelCoords(:,1) > volSize(1)|...
                currVoxelCoords(:,2)<1 | currVoxelCoords(:,2) > volSize(2)|...
                currVoxelCoords(:,3)<1 | currVoxelCoords(:,3) > volSize(3));   % exclude out-of-volume voxels
    currVoxelCoords = currVoxelCoords(~badVoxels,:);
    
    currVoxelInds = sub2ind(volSize, currVoxelCoords(:,1), currVoxelCoords(:,2), currVoxelCoords(:,3)); 
        % indices (not x,y,z coordinates) of voxels around the center that are illuminated by sphere
    currVoxelInds = intersect(gmInds, currVoxelInds);
    [x,y,z] = ind2sub(volSize, currVoxelInds);
    sphereCoords{ii} = [x,y,z];
end