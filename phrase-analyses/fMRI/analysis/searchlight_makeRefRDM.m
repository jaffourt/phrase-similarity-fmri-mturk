function refRDM = searchlight_makeRefRDM(pairsOfPairs, colNames, dissimilarities)
    
%%% Creates a reference representational dissimilarity matrix (RDM)

%%% INPUT:
%%% pairsOfPairs = Px2 cell of cells. Each row is a "pair of contrast pairs":
%%%                In a given row, each entry is a cell of two strings denoting two contrasts to compare.
%%%                Only pairs of pairs denoted here will be included in computing the RDMs
%%% colNames = names of contrast columns in subject data
%%% dissimilarities = dissimilarities between pairs of contrasts
%%%                   (missing dissimilarities are treated as NaN)

%%% SANITY CHECKED (for current project at least) %%%

%% Find the indices of each pairofPairs contrast in colNames %%
relevantCols = {};
for ii = 1:size(pairsOfPairs,1)
    for jj = 1:size(pairsOfPairs,2)
        relevantCols = union(relevantCols, pairsOfPairs{ii,jj});
    end
end

nCols = length(relevantCols);
colInds = zeros(nCols,1);
for c = 1:nCols
    colInds(c) = find(strcmp(relevantCols{c}, colNames));
end
[colInds, origPlaces] = sort(colInds,'ascend');
relevantCols = relevantCols(origPlaces);

%% Fill in the dissimilarity matrix %%
refRDM = nan(nCols);
for ii = 1:nCols
    refRDM(ii,ii) = 0;
end
for d = 1:size(dissimilarities,1)
    goodCols = find(~cellfun(@isempty, strfind(relevantCols, dissimilarities{d,1})));
    for c = 1:length(goodCols)
        col1 = goodCols(c);
        underlineInd = strfind(relevantCols{col1},'_');
        condName = relevantCols{col1}(1:underlineInd);
        col2 = find((~cellfun(@isempty, strfind(relevantCols, condName))) & ...
            (~cellfun(@isempty, strfind(relevantCols, dissimilarities{d,2}))));
        refRDM(col1,col2) = dissimilarities{d,3};
        refRDM(col2,col1) = dissimilarities{d,3};
    end
end
refRDM = squareform(refRDM);
