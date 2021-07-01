function tauVals = searchlight_singleSubj(subjData, refRDM, pairsOfPairs, gmCoords, volSize, radius_mm, voxelSize_mm)
    
%%% Computes a Kendall-tau (type A) correlation between a reference RDM
%%% and dissimilarities from a subject's functional data

%%% INPUT:
%%% subjData: NxC cell of N voxels by C contrast values (or t-values)
%%% refRDM: Mx1 array of reference (dis)similarity values across all contrast pairs, 
%%%         in an order corresponding to that produced by pdist.m 
%%% pairsOfPairs = Kx2 array, indices of pairs of rows in refRDM that should be compared during RSA
%%%                (it is possible to compare only some of the pairs)
%%% gmCoords: Nx3 matrix of number, with (x,y,z) coordinates of voxels in subjData
%%% volSize: 1x3 array of numbers, size of brain volume (in voxels)
%%% radium_mm: 1x3 vector of numbers, radius of sphere in each direction (mm)
%%% voxelSize_mm = 1x3 vector of numbers, size of voxel dimensions (mm)

%%% OUTPUT:
%%% tauVals = Nx1 volume of correlations between the reference RDM and the
%%%           subject's RDM based on a sphere centered around each voxel
%%%           (NaN are inserted for voxels for which 50% of the sphere lies
%%%           outside gmCoords)

%%% Idan Blank, Feb 20, 2017
    
%% Generate spheres %%
relativeCoords = searchlight_makeSphere(radius_mm, voxelSize_mm);
disp('   Defining searchlight spheres');
sphereCoords = searchlight_sphereInds(relativeCoords, gmCoords, gmCoords, volSize);

%% Main %%
nSpheres = length(sphereCoords);
nVoxelsMin = 0.5*size(relativeCoords,1);    % at least half of the voxels in a sphere should be in subjData (i.e., in gray matter)
tauVals = zeros(nSpheres,1);
reverseStr = '';                            % for displaying progress
currPrcnt = 1;                              % for displaying progress
for s = 1:nSpheres  
    doNewLine = 0;
    if (100*s/nSpheres) >= currPrcnt
        fprintf(reverseStr);        
        msg = sprintf('   Sphere %d of %d (%d%%)', s, nSpheres, round(100*s/nSpheres));
        fprintf('   Sphere %d of %d (%d%%)', s, nSpheres, round(100*s/nSpheres));
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        currPrcnt = currPrcnt + 1;
        doNewLine = 1;
    end
    
    nVoxels = size(sphereCoords{s},1);
    if nVoxels >= nVoxelsMin
        %% Find indices of voxels of current sphere in subjData %%
        voxelInds = zeros(nVoxels,1);
        for v = 1:nVoxels
            xGood = gmCoords(:,1) == sphereCoords{s}(v,1);
            yGood = gmCoords(:,2) == sphereCoords{s}(v,2);
            zGood = gmCoords(:,3) == sphereCoords{s}(v,3);
            voxelInds(v) = find(xGood & yGood & zGood);            
        end
                
        %% Compute (dis)similarities %%
        currData = subjData(voxelInds,:)';          % Transpose, so that each row is a beta pattern across the sphere for one contrast
        spearmanDists = pdist(currData,'spearman');  % Spearman correlation between each pair of conditions (rows in currData)
        tauVals(s) = rankCorr_Kendall_taua(spearmanDists, refRDM, pairsOfPairs);
    else
        tauVals(s) = nan;
    end
    
    if doNewLine
        fprintf(1,'\n');
    end
end