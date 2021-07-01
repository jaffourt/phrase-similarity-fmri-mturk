function PhraseRSA_MVPA_relVox(dataPath, ssNames, mvpaConfig)

%%% MVPA analysis: correlation between multivariate response patterns to different sentences
%%% NOTE: this code chooses the top n voxels across the entire gray matter based on response reliability
%%% (for a code that works with pre-defined fROIs, use PhraseRSA_MVPA.m)

%%% INPUT:
%%% dataPath = string, full path to where subjects' data are stored
%%% ssNames = cell of strings (data for subject i have been saved in: ssNames{i}_data.mat
%%% mvpaConfig = structure with 4 fields:
%%%             expt = string, experiment name
%%%             stim = cell of strings, stimulus classes
%%%             conds = cell of strings, condition classes
%%%             chooseVoxels = cell of strings, reliability measures for choosing voxels
%%%                            (should be field names of the relVox variable in the <ssName>_relVox.mat file
%%%             nVoxels = integer array, number of top voxels to choose
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

%%% OUTPUT:
%%% text files for mixed-effects linear modeling in R
%%% Plots of dissimilarities between different conditions and/or stimuli

%% Parameters for computing dissimilarities %%
nStim = numel(mvpaConfig.stim);
nConds = numel(mvpaConfig.conds);
n = nStim*nConds;
nSs = numel(ssNames);

nComps = size(mvpaConfig.comps,1);
relMeas = mvpaConfig.chooseVoxels;  % reliability measures
nArr = mvpaConfig.nVoxels;

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
        colNames = data.(mvpaConfig.expt).colNames;        
        load(fullfile(dataPath, [ssNames{ss}, '_relVox']));

        %% Get columns for different folds / splits %%
        foldCols = cell(nFolds,1);
        for f = 1:nFolds
            for c = 1:numel(folds{f})
                currCols = find(cellfun(@(x,y)(length(folds{f}{c})==length(x(y:end))), ...
                    colNames, strfind(colNames,folds{f}{c})));
                           % to avoid, e.g., indices for Run10 being found when searching for Run1
                foldCols{f} = union(foldCols{f}, currCols);
                           % all columns whose name contain the string folds{f}{c}, regardless of condition / stimulus
            end
        end     
        nPerFold = cellfun(@numel,foldCols)/(nStim*nConds);
        foldTransitionInds = [0;cumsum(nPerFold)];
        
        %% Loop through relability measures + number of top voxels %%
        for m = 1:length(relMeas)
            disp(['        reliability measure: ', relMeas{m}]);
            currMeas = relVox.(relMeas{m});
            
            for nVox = 1:length(nArr)    
                if ss == 1
                    for theDist = 1:nDists
                        fieldName = [mvpaConfig.dataFields{d,1}, '_', dists{theDist}];
                        theDists.(['chooseBy', relMeas{m}]).(['top', num2str(nArr(nVox))]).(fieldName) = nan(nComps,nSs);
                    end
                end
                
                thePrctile = min(100, 100*(nArr(nVox)/numel(data.voxelInds)));
                voxelInds = currMeas >= prctile(currMeas, 100-thePrctile);                    
                vals = nan(sum(voxelInds),n,sum(nPerFold));
                    
                %% Extract data for current subject (loop through ALL stimuli and conditions) %%
                valsCol = 1;
                for stim = 1:nStim
                    stimCols = find(~cellfun(@isempty, strfind(colNames, mvpaConfig.stim{stim})));                    
                    for c = 1:nConds           
                        if (ss==1) && (m==1) && (nVox==1)
                            ind = (stim-1)*nConds+c;
                            trialNames{ind} = [mvpaConfig.stim{stim}, '_', mvpaConfig.conds{c}];
                        end

                        condCols = find(~cellfun(@isempty, strfind(colNames, mvpaConfig.conds{c})));
                        stimCondcols = intersect(condCols, stimCols);                            
                        for fold = 1:nFolds
                            cols = intersect(stimCondcols,foldCols{fold});
                            firstLayer = foldTransitionInds(fold)+1;
                            lastLayer = foldTransitionInds(fold+1);
                            vals(:,valsCol,firstLayer:lastLayer) = permute(x(voxelInds,cols),[1,3,2]);
                        end
                        valsCol = valsCol + 1;
                    end
                end
                nanRows = any(any(isnan(vals),2),3);
                vals = vals(~nanRows,:,:); 
                if mvpaConfig.doMeanSubtraction
                    vals = vals - repmat(mean(mean(vals,2),3),1,n,sum(nPerFold));
                end                    
                    
                %% Compute dissimilarities %%
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
                        theDists.(['chooseBy', relMeas{m}]).(['top', num2str(nArr(nVox))]).(fieldName)(c,ss) = currDist;  
                    end

                    %% Kendall dissimilarities %%
                    if any(strcmp(dists, 'Kendall'))
                        rKendall = rankCorr_Kendall_tauTypeA(vals1,vals2);
                        fieldName = [mvpaConfig.dataFields{d,1}, '_Kendall'];                            
                        theDists.(['chooseBy', relMeas{m}]).(['top', num2str(nArr(nVox))]).(fieldName)(c,ss) = 1-rKendall;
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
                                                                           % that when trials are projected onto it,
                                                                           % discriminates between trials of type 1 and 2 
                                        [~,~,~,stats] = ttest2(r1*w,r2*w); % projection onto w, and t-test of results
                                        currDist = stats.tstat;                                          
                                end
                                fieldName = [mvpaConfig.dataFields{d,1}, '_', otherDists{theDist}];
                                if fold == 1
                                    theDists.(['chooseBy', relMeas{m}]).(['top', num2str(nArr(nVox))]).(fieldName)(c,ss) = (1/nFolds)*currDist;
                                else
                                    theDists.(['chooseBy', relMeas{m}]).(['top', num2str(nArr(nVox))]).(fieldName)(c,ss) = ...
                                        ttheDists.(['chooseBy', relMeas{m}]).(['top', num2str(nArr(nVox))]).(fieldName)(c,ss) + ...
                                        (1/nFolds)*currDist;  
                                end
                            end
                        end
                    end
                end                
                save RSA_reliableVoxels_results theDists                
            end
        end
    end
end
dists = theDists;
trialComps = mvpaConfig.comps;
save RSA_reliableVoxels_results dists trialNames relMeas nVox trialComps


%% Create a table for mixed effects modeling in R + plot the results %%
load RSA_reliableVoxels_results
distFields = fieldnames(dists.(['chooseBy', relMeas{1}]).(['top', num2str(nArr(1))]));
nVoxelSets = numel(relMeas)*numel(nArr);
nRows = nVoxelSets*numel(distFields)*nSs*nComps;
data = cell(nRows,7);   % columns: distance measure, stimulus class, comparison, reliability measure, n voxels, subject, value

row = 1;
pdfNames = cell(nVoxelSets*numel(distFields),1);
figNum = 1;
for m = 1:length(relMeas)
    for nVox = 1:length(nArr)
        for theDist = 1:length(distFields)
            
            vals = dists.(['chooseBy', relMeas{m}]).(['top', num2str(nArr(nVox))]).(distFields{theDist});
            for ss = 1:nSs
                ssName = ['00', num2str(ss)];
                ssName = ssName(end-2:end);
                ssName = ['ss', num2str(ssName)];
                for c = 1:nComps
                    data{row,1} = distFields{theDist};
                    data{row,2} = mvpaConfig.stim{mvpaConfig.comps(c,1)};
                    data{row,3} = [mvpaConfig.conds{mvpaConfig.comps(c,2)}, '_VS_', mvpaConfig.conds{mvpaConfig.comps(c,4)}];
                    data{row,4} = relMeas{m};
                    data{row,5} = ['top', num2str(nArr(nVox))];
                    data{row,6} = ssName;
                    data{row,7} = vals(c,ss);
                    row = row+1;
                end
            end
                
            %% Plot %%
            nGroups = size(mvpaConfig.plotOrder,1);
            nWithinGroup = size(mvpaConfig.plotOrder,2);
            newVals = nan(nSs,nGroups,nWithinGroup);
            for g = 1:nGroups
                for w = 1:nWithinGroup
                    newVals(:,g,w) = vals(mvpaConfig.plotOrder(g,w),:)';
                end
            end
            theTitle = [relMeas{m}, ', top ', num2str(nArr(nVox)), ' voxels: ', distFields{theDist}];
            theFigureAxes = plotData(newVals, dataType, dispersionType, [], ...
                mvpaConfig.theLegend, theTitle, mvpaConfig.colors, [0.1 nGroups+0.9]);
            axes(theFigureAxes);
            theFigure = gcf;
            title(regexprep(theTitle,'_', ' '));
            ylabel(regexprep(distFields{theDist},'_',' '), 'fontname', 'calibri', 'fontsize', 12);
            set(theFigureAxes,'xtick',1:nGroups,'xticklabels',mvpaConfig.xUnits','fontname','calibri','fontsize',12);
            set(theFigure,'units','normalized','paperposition',[0 0 1 1]);
            orient('landscape');                
            print(theFigure, num2str(figNum), '-dpdf');
            pdfNames{figNum} = [num2str(figNum), '.pdf'];                   
            if theDist < length(distFields)
                figNum = figNum + nVoxelSets;         % so that PDFs will be ordered first by voxel set, 
                                                      % then by distance type (and not vice versa)
            else
                figNum = (m-1)*numel(nArr)+nVox;
            end
            close(theFigure)
        end
    end
end

t = cell2table(data, 'variablenames', {'DistanceType', 'StimClass', 'Comparison', 'reliabilityMeasure', 'nVoxels', 'ID', 'Distance'});
save RSA_reliableVoxels_results t -append
writetable(t,'RSA_reliableVoxels_results_forR.txt', 'delimiter', '\t');           

%% Merge PDFs together %%
append_pdfs('MVPA_reliableVoxels.pdf', pdfNames{:});      
for fig = 1:length(pdfNames)
    delete(pdfNames{fig});
end               