function relativeCoords = searchlight_makeSphere(radius_mm, voxelSize_mm)
    
%%% Finds the *relative* voxel indices that constitute a sphere of a given size
%%% (the central voxel is denoted as 0)

%%% INPUT:
%%% radius_mm = 1x3 vector of numbers, radius of sphere in each direction (mm)
%%% voxelSize_mm = 1x3 vector of numbers, size of voxel dimensions (mm)

%%% OUTPUT:
%%% relativeSubs = Nx3 matrix of integers, containing x, y, z coordinates
%%%                included in the sphere relative to its center (0,0,0)

%%% Idan Blank, Feb 20, 2017; copied from Kriegeskorte's RSA toolbox
%%%                           (with changes to variable names)

%%% SANITY CHECKED %%%

radius_voxels = radius_mm./voxelSize_mm;        % from mm to number of voxels
minMargin_vox = floor(radius_voxels);           % take only full voxels
[x,y,z] = meshgrid(-minMargin_vox(1):minMargin_vox(1),...
    -minMargin_vox(2):minMargin_vox(2),...
    -minMargin_vox(3):minMargin_vox(3));

sphere=((x*voxelSize_mm(1)).^2 + ...
    (y*voxelSize_mm(2)).^2 + ...
    (z*voxelSize_mm(3)).^2) <= (radius_mm^2);   % volume with sphere voxels marked 1 and the outside 0

nVoxelsSphere = [size(sphere),ones(1,3-ndims(sphere))];   % enforce 3D, because matlab stupidly autosqueezes trailing singleton dimensions to 2D; try: ndims(ones(1,1,1))
[sphereX, sphereY, sphereZ] = ind2sub(nVoxelsSphere, find(sphere));    % (sub)indices pointing to sphere voxels
sphereSubs = [sphereX, sphereY, sphereZ];
centerCoord = nVoxelsSphere/2 + [0.5 0.5 0.5];            % center position (sphere necessarily has odd number of voxels in each dimension)
relativeCoords = sphereSubs - ones(size(sphereSubs,1),1)*centerCoord;  % center-relative sphere-voxel (sub)indices