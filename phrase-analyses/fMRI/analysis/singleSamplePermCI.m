function ci = singleSamplePermCI(data, nPerms, varargin)

%%% INPUT:
%%% data = a vecor of datapoints
%%% nPerms = number of permutations
%%% theType = optional string, either 'mean' (defult) or 'med' (median)

%%% OUTPUT
%%% ci = 95% confidence interval based on a permutation test
%%%      (All means for which the null hypothesis would not be
%%%       rejected in a permutation test)

if isempty(varargin)
    theType = 'mean';
else
    theType = varargin{1};
end
    
ci = zeros(1,2);
data = data(~isnan(data));
n = length(data);

perms = zeros(nPerms,1);  
switch theType
    case 'mean'
        dataDemean = data-mean(data);        
        for i = 1:nPerms
            inds = ceil(n*rand(n,1));       % observation to sample (with replacement)
            sample = dataDemean(inds);      % permutation sample
            perms(i) = mean(sample)/(std(sample)/sqrt(n));
        end
        [~,inds] = sort(perms(:),'descend');
        perms = perms(inds,:);

        indL = floor(0.025*nPerms);
        ci(1) = mean(data) - perms(indL,1)*(std(data)/sqrt(n)); % Solve for X: (mean(data)-X)/(std(data)/sqrt(n)) = t value at 2.5 percentile

        indU = ceil(0.975*nPerms);
        ci(2) = mean(data) - perms(indU,1)*(std(data)/sqrt(n)); % Solve for X: (mean(data)-X)/(std(data)/sqrt(n)) = t value at 97.5 percentile        
        
    case 'med'
        for i = 1:nPerms
            inds = ceil(n*rand(n,1));
            sample = data(inds);
            perms(i) = prctile(sample,50);
        end
        ci = [prctile(perms,2.5), prctile(perms,97.5)];
end