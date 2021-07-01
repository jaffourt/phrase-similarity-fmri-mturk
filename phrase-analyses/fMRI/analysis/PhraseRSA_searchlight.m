function PhraseRSA_searchlight(dataPath, ssNames, contrastNames, slConfig)
    
%%% This is the main code for running the searchlight analysis

%%% INPUT:
%%% dataPath = string, path to where subject data are stored
%%% ssNames = Sx1 cell of strings, with subject names (corresponding to file names of the form <name>_data.mat)
%%% contrastNames = Cx1 cell of strings, with names of contrasts in the order in which they appear in the subjects' data
%%% pairsOfPairs = Px2 cell of cells. Each row is a "pair of contrast pairs":
%%%                In a given row, each entry is a cell of two strings denoting two contrasts to compare.
%%%                Only pairs of pairs denoted here will be included in computing the RDMs
%%% slConfig = config file as defined in PhraseRSA_Main.m

%%% Idan Blank; Feb 22, 2017

%% Find the indices of each pairofPairs contrast in contrastNames %%
pairsOfPairs = slConfig.pairsOfPairs;
relevantContrasts = {};
for ii = 1:size(pairsOfPairs,1)
    for jj = 1:size(pairsOfPairs,2)
        relevantContrasts = union(relevantContrasts, pairsOfPairs{ii,jj});
    end
end

nContrasts = length(relevantContrasts);
conInds = zeros(nContrasts,1);
for c = 1:nContrasts
    conInds(c) = find(strcmp(relevantContrasts{c}, contrastNames));
end
[conInds, origPlaces] = sort(conInds,'ascend');
relevantContrasts = relevantContrasts(origPlaces);

%% Define indices correpsonding to pairsOfPairs in a way that can be used by the searchlight functions %%
pairsOfPairsInds = zeros(size(pairsOfPairs,1),2);
for p = 1:size(pairsOfPairs,1)
    pair1 = sort([find(strcmp(pairsOfPairs{p,1}{1},relevantContrasts)), ...
        find(strcmp(pairsOfPairs{p,1}{2},relevantContrasts))], 'ascend');
    n = pair1(1)-1;                 % number of items in arithmetic series below
    a1 = nContrasts-1;
    sum1 = (n/2)*(2*a1-1*(n-1));    % sum of arithmetic series: number of entires in pdist up to distances with pair1(1)
    pairsOfPairsInds(p,1) = sum1 + (pair1(2)-pair1(1));
    
    pair2 = sort([find(strcmp(pairsOfPairs{p,2}{1},relevantContrasts)), ...
        find(strcmp(pairsOfPairs{p,2}{2},relevantContrasts))], 'ascend');
    n = pair2(1)-1;
    sum2 = (n/2)*(2*a1-1*(n-1));
    pairsOfPairsInds(p,2) = sum2 + (pair2(2)-pair2(1));
end

%% Main %%
disp('Searclight analysis: ');
for ss = 1:length(ssNames)
    disp([num2str(ss), '. ', ssNames{ss}]);
    load(fullfile(dataPath, [ssNames{ss}, '_data']));
    gmCoords = data.voxelXYZ;
    
%     disp('  Analyzing t-values: ');
%     tData = data.(slConfig.expt).tVals(:,conInds);
%     tTauVals = searchlight_singleSubj(tData, slConfig.refRDM, pairsOfPairsInds, gmCoords, ...
%         slConfig.volSize, slConfig.radius_mm, slConfig.voxelSize_mm);

    disp('  Analyzing contrast values: ');    
    conData = data.(slConfig.expt).conVals(:,conInds);
    conTauVals = searchlight_singleSubj(conData, slConfig.refRDM, pairsOfPairsInds, gmCoords, ...
        slConfig.volSize, slConfig.radius_mm, slConfig.voxelSize_mm); 
    
    ssNum = ['0', num2str(ss)];
    ssNum = ssNum(end-1:end);
    save(['searchlight_ss', ssNum], 'conTauVals');
end