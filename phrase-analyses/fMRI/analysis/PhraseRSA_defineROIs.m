function PhraseRSA_defineROIs(dataPath, ssNames, masks, fROIs, locConfig)

%%% INPUT:
%%% dataPath = string, full path to where subjects' data are stored
%%% ssNames = cell of strings (data for subject i have been saved in: ssNames{i}_data.mat
%%% masks = structure; for each system (field), two embedded fields:
%%%         voxelInds = v X 2 matrix;
%%%                      col 1: coordinate index for voxel v,
%%%                      col 2: integer denoting the mask to which voxel v belongs
%%%         voxelXYZ = v X 4 matrix; same as voxelInds, but with [x,y,z]
%%%                    coordinates in columns 1-3
%%% fROIs = r X 2 matrix. Col 1/2: bottom/top percentile of voxels to take from mask;
%%% locCOnfig = structure; for each system in masks (field), 4 embedded fields:
%%%             name = string, name of experiment that serves as localizer
%%%             conds = cell of strings, conditions (from colNames) for defining fROIs
%%%             weights = array of numbers, weights corresponding to each
%%%                       item in 'conds', for computing the localizer contrast
%%%             vals = 'tVals' or 'conVals', the data from which localizer contrast is computed

%%% OUTPUT:
%%% For each subject, a file in dataPath called <subjectName>_rois.mat,
%%% with a structure variable called roiInds. Fields: fROI definition criterion
%%% (e.g., "from90to100") -> system (m entires, one per mask) -> all/even/odd

systems = fieldnames(masks);
for ss = 1:length(ssNames)
    disp([num2str(ss), '. Defining ROIs for ', ssNames{ss}]);
    load(fullfile(dataPath, [ssNames{ss}, '_data']));
    vInds = data.voxelInds;                                 % indices of GM voxels included in the subject's datasets

    %% Define localizer contrasts %%
    locCon = struct;
    roiInds = struct;    
    for s = 1:length(systems)
        disp(['  ', systems{s}, ' ROIs']);
        locCurr = data.(locConfig.(systems{s}).name).(locConfig.(systems{s}).vals);   % localizer data
        colNames = data.(locConfig.(systems{s}).name).colNames;
        locCon.all = zeros(size(locCurr,1),1);
        locCon.even = zeros(size(locCurr,1),1);
        locCon.odd = zeros(size(locCurr,1),1);
        for c = 1:length(locConfig.(systems{s}).conds)
            colInd = strcmp(colNames, locConfig.(systems{s}).conds{c});
            locCon.all = locCon.all + locCurr(:,colInd)*locConfig.(systems{s}).weights(c);
                        
            colInd = strcmp(colNames, ['EVEN_', locConfig.(systems{s}).conds{c}]);
            locCon.even = locCon.even + locCurr(:,colInd)*locConfig.(systems{s}).weights(c);            
            
            colInd = strcmp(colNames, ['ODD_', locConfig.(systems{s}).conds{c}]);
            locCon.odd = locCon.odd + locCurr(:,colInd)*locConfig.(systems{s}).weights(c);                        
        end

        %% Define fROIs %%
        maskVoxels = masks.(systems{s}).voxelInds;
        nMasks = length(unique(maskVoxels(:,2)));
        for m = 1:nMasks
            mInds = maskVoxels(maskVoxels(:,2)==m,1);       % voxels belonging to mask m
            maskSize = length(mInds);                       % for computing how many voxels to use as an ROI
            [~,~,mIndsInV] = intersect(mInds,vInds);        % voxels in vInds that belong to mask m
            vIndsSubset = vInds(mIndsInV);                  % keeping track of which voxel coordinates are used
            allVals = locCon.all(mIndsInV);                 % localizer contrast values in mask m, all runs
            evenVals = locCon.even(mIndsInV);               % localizer contrast values in mask m, even runs
            oddVals = locCon.odd(mIndsInV);                 % localizer contrast values in mask m, odd runs
            
            %% Loop over threshold for fROI definition %%
            for f = 1:size(fROIs,1)                
                nBtm = (fROIs(f,1)/100)*maskSize;           % number of voxels corresponding to the bottom threshold
                nTop = (fROIs(f,2)/100)*maskSize;           % number of voxels corresponding to the top threshold
                nTotal = round(nTop-nBtm);
                currPrcnt = 100*(nTotal/length(allVals));   % percentage of voxels relative to the subject's number of voxels in this mask
                
                nAbove = maskSize-nTop;                     % number of voxels "above" the top threshold
                if length(allVals)-nAbove < 0
                    btmPrcnt = 0;
                    topPrcnt = min(100, currPrcnt);
                else
                    topPrcnt = 100*(1-(nAbove/length(allVals)));
                    btmPrcnt = max(topPrcnt-currPrcnt,0);
                end

                allInds = find((allVals >= prctile(allVals,btmPrcnt)) & (allVals <= prctile(allVals,topPrcnt)));
                evenInds = find((evenVals >= prctile(evenVals,btmPrcnt)) & (evenVals <= prctile(evenVals,topPrcnt)));
                oddInds = find((oddVals >= prctile(oddVals,btmPrcnt)) & (oddVals <= prctile(oddVals,topPrcnt)));
                
                eval(['roiInds.from', num2str(fROIs(f,1)), 'to', num2str(fROIs(f,2)), '.', ...
                    systems{s}, '(', num2str(m), ').all = vIndsSubset(allInds);']);
                                                            % e.g., for language IFG mask: roiInds.from90to100.language(3).all
                eval(['roiInds.from', num2str(fROIs(f,1)), 'to', num2str(fROIs(f,2)), '.', ...
                    systems{s}, '(', num2str(m), ').even = vIndsSubset(evenInds);']);
                eval(['roiInds.from', num2str(fROIs(f,1)), 'to', num2str(fROIs(f,2)), '.', ...
                    systems{s}, '(', num2str(m), ').odd = vIndsSubset(oddInds);']);                
            end
        end
    end
    
    save(fullfile(dataPath, [ssNames{ss}, '_rois']), 'roiInds');
end