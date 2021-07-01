function PhraseRSA_findReliableVoxels(dataPath, ssNames, relConfig)

%%% Finding reliable voxels for each subject

%%% INPUT:
%%% dataPath = string, full path to where subjects' data are stored
%%% ssNames = cell of strings (data for subject i have been saved in: ssNames{i}_data.mat
%%% relConfig = structure with 4 fields:
%%%             expt = string, experiment name
%%%             dataField = string, name of field in the data file where data are stored
%%%             measures = cell of strings, with measures to compute ('correlation', 'reliability')
%%%             splitHalf = 2x1 cell of cells, where each cell is a different subset ('fold'/'split') of data,
%%%                         denoted by suffixes of the relevant column names (e.g., a cell for 'EVEN' and a 
%%%                         cell for 'ODD'; a cell for 'Run1', 'Run3' and another for 'Run2', 'Run4');
%%%                         NOTE: within each half, data are averaged across runs (if multiple runs exist)

%%% OUTPUT:
%%% For each subject, a file in dataPath called <subjectName>_relVox.mat,
%%% with a structure variable called relVox. Fields: one per measure
%%% specified in relConfig.measures, each with a column vector of
%%% reliability values per voxel

disp('Computing response reliability for each voxel');
nSs = numel(ssNames);
for ss = 1:nSs
    disp(['  ', num2str(ss), '. ', ssNames{ss}]);
    relVox = struct;

    load(fullfile(dataPath, [ssNames{ss}, '_data']));
    x = data.(relConfig.expt).(relConfig.dataField);        
    colNames = data.(relConfig.expt).colNames;

    %% Get columns for the two data halves %%
    shCols = cell(2,1); % split-half columns
    uniqueTrials = cell(2,1);
    for h = 1:2
        reducedColNames = {};
        for c = 1:numel(relConfig.splitHalf{h})
            currCols = find(cellfun(@(x,y)(length(relConfig.splitHalf{h}{c})==length(x(y:end))), ...
                colNames, strfind(colNames,relConfig.splitHalf{h}{c})));
                       % to avoid, e.g., indices for Run10 being found when searching for Run1
            shCols{h} = union(shCols{h}, currCols);
                       % all columns whose name contain the string folds{h}{c}, regardless of condition / stimulus
            reducedColNames = union(reducedColNames, regexprep(colNames(currCols), relConfig.splitHalf{h}{c}, ''));  
        end
        uniqueTrials{h} = unique(reducedColNames);
    end    
    ut = union(uniqueTrials{1}, uniqueTrials{2});
    uniqueTrials = ut;
    n = length(uniqueTrials);

    %% Extract data for current subject %%
    vals = nan(size(x,1),n,2);   % voxels X unique trials X two data halves
    for h = 1:2                 
        for t = 1:n
            trialCols = find(~cellfun(@isempty, strfind(colNames, uniqueTrials{t})));
            cols = intersect(trialCols,shCols{h}); 
            vals(:,t,h) = nanmean(x(:,cols),2);
        end
    end
    nanRows = any(any(isnan(vals),2),3);
    vals = vals(~nanRows,:,:); 

    %% Compute reliability measures %%
    for m = 1:numel(relConfig.measures)
        switch relConfig.measures{m}
            case 'correlation'
                zVals = zscore(vals,1,2);  % z-score each voxel across unique trials (separately for each data half)
                rVals = mean(zVals(:,:,1).*zVals(:,:,2),2);  % correlation = mean of z-score products
            case 'reliability'
                a = vals(:,:,1);
                b = vals(:,:,2);
                projAonB = (repmat(sum(b.*a,2),1,n).*b)./repmat(sum(b.*b,2),1,n);   % from Sam NH's 2015 Neuron paper
                aResid = a - projAonB;
                rVals = 1 - (sum(aResid.*aResid,2))./(sum(a.*a,2));
        end
        relVox.(relConfig.measures{m}) = rVals;
    end

    save(fullfile(dataPath, [ssNames{ss}, '_relVox']), 'relVox');
end