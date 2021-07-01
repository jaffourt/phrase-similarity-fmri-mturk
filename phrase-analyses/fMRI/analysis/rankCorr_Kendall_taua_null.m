function nullDist = rankCorr_Kendall_taua_null(a,varargin)

    %%% Computes the null distribution of Kendall tau-a correlations

    %%% INPUT:
    %%% a: vector of pairwise (dis)similarities
    %%% pairs (optional): Nx2 array of numbers denoting pairs of pairs to compare across a and its permuted versions
    %%%                   (in case only a subset of the pairs should be compared)
    %%% subsets (optional) = Sx1 cell of vectors, each one with indices of
    %%%                      valuse that come from within the same subset (values will only be
    %%%                      shuffled within subsets)
    %%% NOTE: pair members could only be within the same subset

    %%% OUTPUT:
    %%% nullDist = px2 matrix; col1 = correlation, col2 = associated p value

    %% Initalize variables %%
    if isempty(varargin)
        nAll = nchoosek(length(a),2);
        pairs = zeros(nAll,2);
        ind1 = 1;
        ind2 = 2;
        for p = 1:nAll
            pairs(p,:) = [ind1, ind2];
            ind2 = ind2+1;
            if ind2 > length(a)
                ind1 = ind1+1;
                ind2 = ind1+1;
            end
        end
    else
        pairs = varargin{1};
    end

    if length(varargin) <= 1
        subsets = {1:length(a)};
    else
        subsets = varargin{2};
    end

    %% For each subset, compute its contribution to the nominator of tau_a for each permutation of items in that subset %%
    nPairs = size(pairs,1);
    nSubsets = length(subsets);
    nullDist = zeros(2*nPairs+1,2);                      % values between -1 and 1
    nullDist(:,1) = -1:(1/nPairs):1;
    nullDist(:,1) = round((10^7)*nullDist(:,1))/(10^7);    % because Matlab does not understand that 1/3 == 1/3

    nomVals = cell(nSubsets,1);        % nominator of tau-A for each subset = concordant pairs - discordant pairs
    for s = 1:nSubsets
        aSub1 = a(subsets{s});        
        pairInds = zeros(nPairs,1);

        %% Find pairs in current subset %%
        for p = 1:nPairs
            if any(pairs(p,1)==subsets{s}) && any(pairs(p,2)==subsets{s})
                pairInds(p)=1;
            end
        end
        pairInds = pairInds > 0;
        currPairs = pairs(pairInds,:);
        currPairs = currPairs - min(subsets{s}) + 1;    % shift indices to start from 1 (for rankCorr_Kendall_taua.m)
        
        %% Compute all possible permutations of indices in subsets{s} %%
        if factorial(length(subsets{s})) <= (5*10^(5/nSubsets))          % analytical solution: all data permutations 
                                                                         % (up to 5*10^6 permutations when multiplying across all subsets)
            allPerms = makeAllPerms(subsets{s});         
            nomVals{s} = zeros(size(allPerms,1),1);
            for p = 1:size(allPerms,1)
                aSub2 = a(allPerms(p,:));
                nomVals{s}(p) = rankCorr_Kendall_taua(aSub1,aSub2,currPairs)*size(currPairs,1);    % nominator of tau-A = concordant pairs - discordant pairs
            end
        else
            nPerms = floor(5*10^(6/nSubsets));
            nomVals{s} = zeros(nPerms,1);        
            for p = 1:nPerms
                aSub2 = a(permute(subsets{s}));
                nomVals{s}(p) = rankCorr_Kendall_taua(aSub1,aSub2,currPairs)*size(currPairs,1);    % nominator of tau-A = concordant pairs - discordant pairs
            end    
        end
    end
    
    %% Combinations of nominators across subsets, to produce tau_a values %%
    subsetInds = ones(nSubsets,1);  % counter for computing all possible combinations of values across nomVals  
    currSubset = nSubsets;
    stopWhile = 0;    
    while ~stopWhile %%% sum(subsetInds == cellfun(@length, nomVals)) <= nSubsets               
        tau_a = sum(cellfun(@(x,y) x(y), nomVals, num2cell(subsetInds))) / nPairs;
        tau_a = round((10^7)*tau_a)/(10^7);        
        nullDist(nullDist(:,1)==tau_a,2) = nullDist(nullDist(:,1)==tau_a,2) + 1;

        %% Update counters %%
        if subsetInds(currSubset)+1 <= length(nomVals{currSubset})
            subsetInds(currSubset) = subsetInds(currSubset) + 1;
        else
            isGood = 0;
            while ~isGood
                subsetInds(currSubset) = 1;
                currSubset = currSubset - 1;
                if currSubset > 0
                    subsetInds(currSubset) = subsetInds(currSubset)+1;
                    isGood = subsetInds(currSubset) <= length(nomVals{currSubset});
                else
                    stopWhile = 1;
                    break
                end
            end
            currSubset = nSubsets;
        end        
    end         
    nNullVals = prod(cellfun(@length,nomVals)); % from counts to probabilities
    nullDist(:,2) = nullDist(:,2)/nNullVals;
    
end

function allPerms = makeAllPerms(v)
    %%% Computes all possible permutations of vector v %%%
    
    v = sort(v,'ascend');
    if size(v,1) > 1
        v = v';
    end    
    n = length(v);

    nPerms = factorial(n)-1;
    allPerms = zeros(nPerms,n);
    goodInd = 1;
    row = 1;
    while ~isempty(goodInd)
        goodInd = find(v(1:n-1)-v(2:n)<0,1,'last');
        if ~isempty(goodInd)
            swapInd = find(v(goodInd)-v((goodInd+1):n)<0,1,'last')+goodInd;
            val1 = v(goodInd);
            val2 = v(swapInd);
            v(goodInd) = val2;
            v(swapInd) = val1;
            v((goodInd+1):n) = fliplr(v((goodInd+1):n));
            allPerms(row,:) = v;
            row = row + 1;
        end
    end
end