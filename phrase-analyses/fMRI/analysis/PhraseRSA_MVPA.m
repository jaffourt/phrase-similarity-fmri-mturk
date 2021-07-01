function PhraseRSA_MVPA(dataPath, ssNames, mvpaConfig)

%%% MVPA analysis: correlation between multivariate response patterns to different sentences

%%% INPUT:
%%% dataPath = string, full path to where subjects' data are stored
%%% ssNames = cell of strings (data for subject i have been saved in: ssNames{i}_data.mat
%%% mvpaConfig = structure with 11 fields:
%%%             expt = string, experiment name
%%%             stim = cell of strings, stimulus classes
%%%             conds = cell of strings, condition classes
%%%             comps = c x 4 matrix of integers, each row contains two stim-cond indices
%%%                     (e.g., [1,2,1,4] = stim1-cond2 vs. stim1-cond4), for one MVPA comparison
%%%             dataFields = f x 3 cell, with each row f denoting:
%%%                          col1 = string, name of field in the data file where data are stored
%%%                          col2 = cell of strings, with distance measures to compute
%%%                                 (any pdist function, 'wbCorr' for Haxby within-between correlation difference,
%%%                                 LDt for Fisher linear discriminant, kendall for Kendall's tau type A)
%%%                          col3 = cell of cells of strings, where each cell is a different subset ('fold'/'split') of data
%%%                                 for cross validation (not used for pdist distances & Kendall), denoted by suffixes of 
%%%                                 the relevant column names (e.g., a cell for 'EVEN' and a cell for 'ODD';
%%%                                 a cell for 'Run1', 'Run3' and another for 'Run2', 'Run4');
%%%             doMeanSubtraction = boolean, whether to subtract from every beta of every voxel its average
%%%                                 across stimuli and conditions
%%%             plotOrder = 2D matrix, each row contains indices for rows in "comps"
%%%                         that will be grouped together in the barplot
%%%                         (in the particular order they are indexed here)
%%%             theLegend = legend for levels of within-group variable
%%%             xUnits = cell of strings, names for levels of grouping variable
%%%             colors = Lx3 array, with colors for each of the levels in theLegend
%%%             hemiSplit = array; for each cortical system, 1 to plot first half of masks separately
%%%                         from second half (i.e., split by hemisphere)

%%% OUTPUT:
%%% text files for mixed-effects linear modeling in R
%%% Plots of dissimilarities between different conditions and/or stimuli

%%% SANITY CHEKED, AUGUST 28, 2018 

%% Parameters for computing dissimilarities %%
nStim = numel(mvpaConfig.stim);
nConds = numel(mvpaConfig.conds);
n = nStim*nConds;
nSs = numel(ssNames);

nComps = size(mvpaConfig.comps,1);

dataType = 'indMed';
dispersionType = 'ci';


%% Compute relevant pairwise comparisons %%
disp('Multivariate pattern comparison:');
theDists = struct;
for d = 1:size(mvpaConfig.dataFields,1)
    disp(['  Analysis set ', num2str(d), ' (out of ', num2str(size(mvpaConfig.dataFields,1)), ')']);
    dists = mvpaConfig.dataFields{d,2};
    nDists = numel(dists);
    folds = mvpaConfig.dataFields{d,3};
    nFolds = numel(folds);
    trialNames = cell(n,1);  % for columns of the variable 'vals'
    
    %% Loop through subjects %%
    for ss = 1:nSs
        disp(['    ', num2str(ss), '. ', ssNames{ss}]);
        load(fullfile(dataPath, [ssNames{ss}, '_data']));
        x = data.(mvpaConfig.expt).(mvpaConfig.dataFields{d,1});        
        colNames = data.(mvpaConfig.expt).([mvpaConfig.dataFields{d,1}, 'ColNames']);
        
        load(fullfile(dataPath, [ssNames{ss}, '_rois']));
        if ss == 1
            fROI_crits = fieldnames(roiInds);
            systems = fieldnames(roiInds.(fROI_crits{1}));
        end
        v = data.voxelInds;                                         % all GM voxel coordinates for current subject   

       %% Get columns for different folds / splits %%
        foldCols = cell(nFolds,1);
        for f = 1:nFolds
            for c = 1:numel(folds{f})
                currCols = find(cellfun(@(x,y)(strcmp(folds{f}{c},x(y:(y+length(folds{f}{c})-1)))), ...
                    colNames, strfind(colNames,folds{f}{c})));
                           % (to avoid, e.g., indices for Run10 being found when searching for Run1)
                foldCols{f} = union(foldCols{f}, currCols);
                           % all columns whose name contains the string folds{f}{c}, regardless of condition / stimulus
            end
        end     
        nPerFold = cellfun(@numel,foldCols)/(nStim*nConds);
        foldTransitionInds = [0;cumsum(nPerFold)];
        
       %% Loop through fROI definition criteria + cortical systems %%
        for f = 1:length(fROI_crits)
            disp(['        fROI criterion: ', fROI_crits{f}]);
            for s = 1:length(systems)    
                nROI = numel(roiInds.(fROI_crits{f}).(systems{s}));   
                for r = 1:(nROI+2)      % last two entries: all ROIs in LH, all ROIs in RH
                    if ss == 1
                        for theDist = 1:nDists
                            fieldName = [mvpaConfig.dataFields{d,1}, '_', dists{theDist}];
                            theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).criticalComps = nan(nComps,nSs);
                            if (~strcmp(dists{theDist},'LDt')) && (~strcmp(dists{theDist}, 'wbCorr'))
                                theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).allPairwiseItems = nan(n,n,nSs); 
                            end
                        end                    
                    end

                    if r <= nROI
                        disp(['            ', systems{s}, ' ', num2str(r)]);                                            
                        inds = roiInds.(fROI_crits{f}).(systems{s})(r).all;
                        [~, ~, indsInV] = intersect(inds, v);   % all inds should be in v
                        if r == 1
                            allROIsIndsInV = cell(1,2);         % one entry per hemisphere
                            allROIsIndsInV{1} = indsInV;
                            hemi = 1;
                        elseif r == ((nROI/2)+1)
                            allROIsIndsInV{2} = indsInV;
                            hemi = 2;
                        else
                            allROIsIndsInV{hemi} = union(allROIsIndsInV{hemi}, indsInV);
                        end                        
                    else
                        disp('             all ROIs in each hemisphere');
                        indsInV = allROIsIndsInV{r-nROI};
                    end
                    vals = nan(length(indsInV),n,sum(nPerFold));
                    
                  %% Extract data for current subject & fROI (loop through ALL stimuli and conditions) %%
                    valsCol = 1;
                    for stim = 1:nStim
                        stimCols = find(~cellfun(@isempty, strfind(colNames, mvpaConfig.stim{stim})));                    
                        for c = 1:nConds           
                            if (ss==1) && (f==1) && (s==1) && (r==1)
                                ind = (stim-1)*nConds+c;
                                trialNames{ind} = [mvpaConfig.stim{stim}, '_', mvpaConfig.conds{c}];
                            end

                            condCols = find(~cellfun(@isempty, strfind(colNames, mvpaConfig.conds{c})));
                            stimCondcols = intersect(condCols, stimCols);                            
                            for fold = 1:nFolds
                                cols = intersect(stimCondcols,foldCols{fold});
                                firstLayer = foldTransitionInds(fold)+1;
                                lastLayer = foldTransitionInds(fold+1);
                                vals(:,valsCol,firstLayer:lastLayer) = permute(x(indsInV,cols),[1,3,2]);
                            end
                            valsCol = valsCol + 1;
                        end
                    end                    
                    nanRows = any(any(isnan(vals),2),3);
                    vals = vals(~nanRows,:,:); 
                    if mvpaConfig.doMeanSubtraction
                     %% Version 1: mass subtraction %%
                        vals = vals - repmat(mean(mean(vals,2),3),1,n,sum(nPerFold));   % subtracting the mean of each voxel's activation
                        
                     %% Version 2: subtraction per stimulus set %%
%                         cols = 1:nConds;
%                         for stim = 1:nStim
%                             currVals = vals(:,cols,:);
%                             currVals = currVals - repmat(mean(mean(currVals,2),3),1,nConds,sum(nPerFold));
%                             vals(:,cols,:) = currVals;
%                             cols = cols + nConds;
%                         end                            
                    end
                    
                  %% Compute all pairwise dissimilarities (using data folds) %%
                    distsForPDist = setdiff(dists, {'wbCorr', 'LDt', 'Kendall'});
                    for fold = 1:nFolds
                        foldFirstInd = foldTransitionInds(fold)+1;
                        foldLastInd = foldTransitionInds(fold+1);
                        foldInds = foldFirstInd:foldLastInd;                        
                                   % layers in vals for trials belonging to the current fold
                        foldVals = vals(:,:,foldInds);                       
                        restInds = setdiff(1:foldTransitionInds(end), foldInds);    
                                   % layers in vals for trials belonging to all other folds
                        restVals = vals(:,:,restInds);    
                        
                        for theDist = 1:numel(distsForPDist)
                            if fold == 1
                                allDists = zeros(n,n,1);
                                eval(['allDists_', distsForPDist{theDist}, ' = allDists;']);
                            end
                            eval(['allDists = allDists_', distsForPDist{theDist}, ';']);
                            allDists = allDists + ...
                                (1/nFolds)*pdist2((nanmean(foldVals,3))', (nanmean(restVals,3))', distsForPDist{theDist});
                            eval(['allDists_', distsForPDist{theDist}, ' = allDists;']);                                
                            if fold == nFolds
                                fieldName = [mvpaConfig.dataFields{d,1}, '_', dists{theDist}];                                
                                theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).allPairwiseItems(:,:,ss) = allDists;
                            end                                
                        end
                        
                        if any(strcmp(dists, 'Kendall'))
                            if fold == 1
                                kendallVals = zeros(n,n);
                            end
                            for item1 = 1:n
                                for item2 = item1:n
                                    kendallVal = 1 - rankCorr_Kendall_tauTypeA(nanmean(foldVals(:,item1,:),3), nanmean(restVals(:,item2,:),3));
                                    kendallVals(item1,item2) = kendallVals(item1,item2) + (1/nFolds)*kendallVal;
                                    if ~(item1 == item2)
                                        kendallVals(item2,item1) = kendallVals(item2,item1) + (1/nFolds)*kendallVal;
                                    end
                                end
                            end                            
                            if fold == nFolds
                                fieldName = [mvpaConfig.dataFields{d,1}, '_Kendall'];                                
                                theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).allPairwiseItems(:,:,ss) = kendallVals;
                            end
                        end                            
                    end
                    
                  %% Compute dissimilarities: critical comparisons %%
                    distsForPDist = setdiff(dists, {'wbCorr', 'LDt', 'Kendall'});
                    for c = 1:nComps
                        trial1 = strcmp(trialNames, ...
                            [mvpaConfig.stim{mvpaConfig.comps(c,1)}, '_', mvpaConfig.conds{mvpaConfig.comps(c,2)}]);
                        trial2 = strcmp(trialNames, ...
                            [mvpaConfig.stim{mvpaConfig.comps(c,3)}, '_', mvpaConfig.conds{mvpaConfig.comps(c,4)}]);
                        vals1 = nanmean(vals(:,trial1,:),3);
                        vals2 = nanmean(vals(:,trial2,:),3);                                               
                        
                     %% pdist dissimilarities: not using data folds %%
                        for theDist = 1:numel(distsForPDist)
                            currDist = pdist2(vals1', vals2', distsForPDist{theDist});
                            fieldName = [mvpaConfig.dataFields{d,1}, '_', distsForPDist{theDist}];
                            theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).criticalComps(c,ss) = currDist; 
                        end
                        
                     %% Kendall dissimilarities %%
                        if any(strcmp(dists, 'Kendall'))
                            rKendall = rankCorr_Kendall_tauTypeA(vals1,vals2);
                            fieldName = [mvpaConfig.dataFields{d,1}, '_Kendall'];                            
                            theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).criticalComps(c,ss) = 1-rKendall;
                        end
                        
                     %% Distances using data folds %%
                        otherDists = setdiff(dists, union(distsForPDist,'Kendall'));
                        if ~isempty(otherDists)
                            for fold = 1:nFolds
                                foldFirstInd = foldTransitionInds(fold)+1;
                                foldLastInd = foldTransitionInds(fold+1);
                                foldInds = foldFirstInd:foldLastInd;                        
                                           % layers in vals for trials belonging to the current fold
                                foldVals = vals(:,:,foldInds);                       
                                restInds = setdiff(1:foldTransitionInds(end), foldInds);    
                                           % layers in vals for trials belonging to all other folds
                                restVals = vals(:,:,restInds);

                                for theDist = 1:numel(otherDists)
                                    switch otherDists{theDist}
                                        case 'wbCorr'
                                            f1 = nanmean(foldVals(:,trial1,:),3);
                                            f2 = nanmean(foldVals(:,trial2,:),3);
                                            r1 = nanmean(restVals(:,trial1,:),3);
                                            r2 = nanmean(restVals(:,trial2,:),3);
                                            rBetween = corr([f1,r1],[f2,r2],'type','spearman');
                                            rBetween = tanh(mean(atanh(rBetween(:))));
                                            rWithin = [corr(f1,r1,'type','spearman'), corr(f2,r2,'type','spearman')];
                                            rWithin = tanh(mean(atanh(rWithin)));
                                            currDist = (rWithin-rBetween);                                            
                                            
                                        case 'LDt'
                                            f1 = permute(foldVals(:,trial1,:),[3,1,2]); % trials x voxels
                                            f2 = permute(foldVals(:,trial2,:),[3,1,2]);
                                            r1 = permute(restVals(:,trial1,:),[3,1,2]);
                                            r2 = permute(restVals(:,trial2,:),[3,1,2]);
                                            
                                       %% Train on current fold, test on rest %%
                                            mu_f1 = mean(f1,1);
                                            mu_f2 = mean(f2,1);
                                            sigma_f1 = covdiagFromRSAtoolbox(f1);
                                            sigma_f2 = covdiagFromRSAtoolbox(f2);
                                            sigmaPooled = size(f1,1)*sigma_f1 + size(f2,1)*sigma_f2;
                                            w = sigmaPooled\(mu_f1'-mu_f2');   % this is the linear discriminant: a vector
                                                                               % than when trials are projected onto it,
                                                                               % discriminates between trials of type 1 and 2 
                                            [~,~,~,stats] = ttest2(r1*w,r2*w); % projection onto w, and t-test of results
                                            currDist = stats.tstat;                                          
                                    end
                                    fieldName = [mvpaConfig.dataFields{d,1}, '_', otherDists{theDist}];
                                    if fold == 1
                                        theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).criticalComps(c,ss) = ...
                                            (1/nFolds)*currDist;
                                    else
                                        theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).criticalComps(c,ss) = ...
                                            theDists.(fROI_crits{f}).(systems{s})(r).(fieldName).criticalComps(c,ss) + ...
                                            (1/nFolds)*currDist;  
                                    end
                                end
                            end
                        end
                    end
                end   
                save RSA_ROI_results theDists                         
            end
        end
    end
end

dists = struct;
dists2 = struct;
for f = 1:length(fROI_crits)
    for s = 1:length(systems)
        dists.(fROI_crits{f}).(systems{s}) = theDists.(fROI_crits{f}).(systems{s})(1:(end-2));      % individual ROIs
        dists2.(fROI_crits{f}).(systems{s}) = theDists.(fROI_crits{f}).(systems{s})((end-1):end);   % all ROIs in LH, RH
    end
end
trialComps = mvpaConfig.comps;
save RSA_ROI_results dists trialNames fROI_crits systems trialComps

clear dists
dists = dists2;
clear dists2
save RSA_allROIsPerHemi_results dists trialNames fROI_crits systems trialComps


%% Create a table for mixed effects modeling in R %%
load RSA_ROI_results
theDists = dists;
load RSA_allROIsPerHemi_results
dists2 = dists;
for f = 1:length(fROI_crits)
    for s = 1:length(systems)
        theDists.(fROI_crits{f}).(systems{s})(end+1:end+2) = dists2.(fROI_crits{f}).(systems{s});     % individual ROIs
    end
end
distFields = fieldnames(theDists.(fROI_crits{1}).(systems{1})(1));
nROI_total = 0;
for s = 1:length(systems)    
    nROI_total = nROI_total + numel(theDists.(fROI_crits{1}).(systems{s}));
end
nRows = length(fROI_crits)*nROI_total*numel(distFields)*nSs*nComps;
data = cell(nRows,7);   % columns: distance measure, stimulus class, comparison, fROI criterion, system/fROI, subject, value
row = 1;

for f = 1:length(fROI_crits)
    for s = 1:length(systems)
        nROI = numel(theDists.(fROI_crits{f}).(systems{s})) - 2;
        for r = 1:(nROI+2)            
            %% Make table for R %%
            if r <= nROI
                roiNum = ['00', num2str(r)];
                roiNum = roiNum(end-2:end);
                roiName = [systems{s}, '_r', roiNum];     
            elseif r == (nROI + 1)
                roiName = [systems{s}, '_LH'];
            else
                roiName = [systems{s}, '_RH'];
            end
            
            for theDist = 1:length(distFields)
                vals = theDists.(fROI_crits{f}).(systems{s})(r).(distFields{theDist}).criticalComps;
                for ss = 1:nSs
                    ssName = ['00', num2str(ss)];
                    ssName = ssName(end-2:end);
                    ssName = ['ss', num2str(ssName)];
                    for c = 1:nComps
                        data{row,1} = distFields{theDist};
                        data{row,2} = mvpaConfig.stim{mvpaConfig.comps(c,1)};
                        data{row,3} = [mvpaConfig.conds{mvpaConfig.comps(c,2)}, '_VS_', mvpaConfig.conds{mvpaConfig.comps(c,4)}];
                        data{row,4} = fROI_crits{f};
                        data{row,5} = roiName;
                        data{row,6} = ssName;
                        data{row,7} = vals(c,ss);
                        row = row+1;
                    end
                end
            end
        end           
    end
end

t = cell2table(data, 'variablenames', {'DistanceType', 'StimClass', 'Comparison', 'fROI_Criterion', 'fROI', 'ID', 'Distance'});
save RSA_ROI_results t -append
writetable(t,'RSA_ROI_results_forR.txt', 'delimiter', '\t');           