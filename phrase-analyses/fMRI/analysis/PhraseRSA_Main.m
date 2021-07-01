%%% Main code for analyzing the Phrase RSA data %%%
%%% Idan Blank, Nov 13 2017; EvLab Rulz!        %%%

%% Parameters %%
dataPath = '/mindhive/evlab/u/iblank/Desktop/Projects/PhraseRSA/Data';
load(fullfile(dataPath, 'subjectInfo'));
ssNames = subjects.main;

stimSets = [1,2,14,30];     % the four sentence sets from the behavioral experiment used in the fMRI experiment
stimSetNames = {'MAN', 'CAT', 'PEN', 'GLASS'};
stimTypes = {2, 'MEANINGPRES';
    3, 'NCHANGE';
    4, 'PREPCHANGE';
    5, 'AJCHANGE'};

load(fullfile(dataPath, 'masks'));

fROIs = [0 100;
    90 100;
    80 100;
    70 100;
    80 90;
    70 80
    50 100];         % percentiles of voxels to take


%% Define ROIs %%
% locConfig.lang.name = 'langloc';
% locConfig.lang.conds = {'S', 'N'};
% locConfig.lang.weights = [1, -1];       % language: S > N
% locConfig.lang.vals = 'conVals';
% 
% locConfig.MD.name = 'langloc';
% locConfig.MD.conds = {'S', 'N'};        % MD: N > S
% locConfig.MD.weights = [-1, 1];
% locConfig.MD.vals = 'conVals';
% 
% locConfig.DMN.name = 'langloc';
% locConfig.DMN.conds = {'N'};
% locConfig.DMN.weights = -1;
% locConfig.DMN.vals = 'conVals';         % DMN: FIX > N
% 
% locConfig.vis.name = 'langloc';
% locConfig.vis.conds = {'S'};
% locConfig.vis.weights = 1;              % visual: S > FIX
% locConfig.vis.vals = 'conVals';
% 
% PhraseRSA_defineROIs(dataPath, ssNames, masks, fROIs, locConfig);


%% Univariate analysis %%
% uniConfig.langloc.cols = {'S', 'N'};
% uniConfig.langloc.avg = {};
% uniConfig.langloc.hemiSplit = [1, 1, 1, 0];                         % 1 = plot separately for LH and RH fROIs
% uniConfig.langloc.isloc = {'lang', 'MD', 'DMN', 'vis'};
% uniConfig.langloc.vals = 'conVals';
% uniConfig.complang.cols = {'MAN_BASE', 'MAN_NCHANGE', 'MAN_AJCHANGE', 'MAN_PREPCHANGE', 'MAN_MEANINGPRES', ... 
%     'CAT_BASE', 'CAT_NCHANGE', 'CAT_AJCHANGE', 'CAT_PREPCHANGE', 'CAT_MEANINGPRES', ...
%     'PEN_BASE', 'PEN_NCHANGE', 'PEN_AJCHANGE', 'PEN_PREPCHANGE', 'PEN_MEANINGPRES', ...
%     'GLASS_BASE', 'GLASS_NCHANGE', 'GLASS_AJCHANGE', 'GLASS_PREPCHANGE', 'GLASS_MEANINGPRES'};
% uniConfig.complang.avg = {'BASE', 'NCHANGE', 'AJCHANGE', 'PREPCHANGE', 'MEANINGPRES'};
% uniConfig.complang.hemiSplit = [1, 1, 1, 0];                        % 1 = plot separately for LH and RH fROIs
% uniConfig.complang.isloc = {};
% uniConfig.complang.vals = 'tVals';
% 
% PhraseRSA_univariate(dataPath, ssNames, uniConfig); 


%% Multivariate analysis of fROIs: extracting data for mixed-effects modeling %%
mvpaConfig.expt = 'complang';
mvpaConfig.stim = {'MAN', 'CAT', 'PEN', 'GLASS'};
mvpaConfig.conds = {'BASE', 'NCHANGE', 'AJCHANGE', 'PREPCHANGE', 'MEANINGPRES'};
mvpaConfig.dataFields = {'conVals', {'spearman', 'wbCorr'}, ...
    {{'Run1', 'Run3', 'Run5', 'Run7', 'Run9'};
    {'Run2', 'Run4', 'Run6', 'Run8', 'Run10'}}; ...
    'tVals', {'spearman', 'wbCorr'}, {{'EVEN'}; {'ODD'}}};
mvpaConfig.doMeanSubtraction = false;
mvpaConfig.comps = [1,1, 1,2;
    1,1, 1,3;
    1,1, 1,4;
    1,1, 1,5;
    2,1, 2,2;
    2,1, 2,3;
    2,1, 2,4;
    2,1, 2,5;
    3,1, 3,2;
    3,1, 3,3;
    3,1, 3,4;
    3,1, 3,5;
    4,1, 4,2;
    4,1, 4,3;
    4,1, 4,4;
    4,1, 4,5];      % cols 1,3 are stim, columns 2,4 conds    
mvpaConfig.plotOrder = [1,2,3,4;
    5,6,7,8;
    9,10,11,12;
    13,14,15,16];   % rows from comps
mvpaConfig.theLegend = {'N vs. BASE', 'ADJ vs. BASE', 'PREP vs. BASE', 'PRESERVED vs. BASE'};
mvpaConfig.xUnits = {'MAN', 'CAT', 'PEN', 'GLASS'};
mvpaConfig.colors = [255 178 102;
    51 153 255;
    255 51 51;
    0 153 0]/255;
mvpaConfig.hemiSplit = [1, 1, 1, 0];    % 1 = plot separately for LH and RH fROIs
PhraseRSA_MVPA(dataPath, ssNames, mvpaConfig);    % ONLY NEED TO DO THIS ONCE
return

%% Compute split-half reliability for each voxel %%
% relConfig.expt = 'complang';
% relConfig.dataField = 'conVals';
% relConfig.measures = {'correlation', 'reliability'};
% relConfig.splitHalf = {{'Run1', 'Run3', 'Run5', 'Run7', 'Run9'};
%     {'Run2', 'Run4', 'Run6', 'Run8', 'Run10'}};
% PhraseRSA_findReliableVoxels(dataPath, ssNames, relConfig);


%% Multivariate analysis of reliable voxels: extracting data for mixed-effects modeling %%
mvpaConfig.expt = 'complang';
mvpaConfig.stim = {'MAN', 'CAT', 'PEN', 'GLASS'};
mvpaConfig.conds = {'BASE', 'NCHANGE', 'AJCHANGE', 'PREPCHANGE', 'MEANINGPRES'};
mvpaConfig.dataFields = {'conVals', {'spearman'}, ...
    {{'Run1', 'Run3', 'Run5', 'Run7', 'Run9'};
    {'Run2', 'Run4', 'Run6', 'Run8', 'Run10'}}};
mvpaConfig.chooseVoxels = {'correlation', 'reliability'};
mvpaConfig.nVoxels = [100; 500; 1000; 2500; 5000; 10000];
mvpaConfig.comps = [1,1, 1,2;
    1,1, 1,3;
    1,1, 1,4;
    1,1, 1,5;
    2,1, 2,2;
    2,1, 2,3;
    2,1, 2,4;
    2,1, 2,5;
    3,1, 3,2;
    3,1, 3,3;
    3,1, 3,4;
    3,1, 3,5;
    4,1, 4,2;
    4,1, 4,3;
    4,1, 4,4;
    4,1, 4,5];      % cols 1,3 are stim, columns 2,4 conds    
mvpaConfig.doMeanSubtraction = false;
mvpaConfig.plotOrder = [1,2,3,4;
    5,6,7,8;
    9,10,11,12;
    13,14,15,16];   % rows from comps
mvpaConfig.theLegend = {'N vs. BASE', 'ADJ vs. BASE', 'PREP vs. BASE', 'PRESERVED vs. BASE'};
mvpaConfig.xUnits = {'MAN', 'CAT', 'PEN', 'GLASS'};
mvpaConfig.colors = [255 178 102;
    51 153 255;
    255 51 51;
    0 153 0]/255;
PhraseRSA_MVPA_relVox(dataPath, ssNames, mvpaConfig);    % ONLY NEED TO DO THIS ONCE

locConfig.expt = 'langloc';
locConfig.vals = 'conVals';
locConfig.conds = {'S', 'N'};
locConfig.colors = [150 0 0; 0 0 150; 200 200 200]/255;
locConfig.chooseVoxels = {'correlation', 'reliability'};
locConfig.nVoxels = [100; 500; 1000; 2500; 5000; 10000];
PhraseRSA_MVPA_relvox_locData(dataPath, ssNames, locConfig)





%% Searchlight analysis %%
%%% NOTE: at some point it is worth testing more pairs, such as N-vs.Prep = AJ change
%%% NOTE: also worth testing whether pairs from different stimulus sets are maximally different

% slConfig.expt = 'complang';
% slConfig.volSize = [79 95 69];
% slConfig.voxelSize_mm = [2 2 2];
% slConfig.radius_mm = 6;
% slConfig.pairsOfPairs = {{'MAN_BASE', 'MAN_MEANINGPRES'}, {'MAN_BASE', 'MAN_AJCHANGE'};
%     {'MAN_BASE', 'MAN_MEANINGPRES'}, {'MAN_BASE', 'MAN_NCHANGE'};
%     {'MAN_BASE', 'MAN_PREPCHANGE'}, {'MAN_BASE', 'MAN_AJCHANGE'};
%     {'MAN_BASE', 'MAN_PREPCHANGE'}, {'MAN_BASE', 'MAN_NCHANGE'};
%     {'MAN_BASE', 'MAN_AJCHANGE'}, {'MAN_BASE', 'MAN_NCHANGE'};   
% 
%     {'CAT_BASE', 'CAT_MEANINGPRES'}, {'CAT_BASE', 'CAT_AJCHANGE'};
%     {'CAT_BASE', 'CAT_MEANINGPRES'}, {'CAT_BASE', 'CAT_NCHANGE'};
%     {'CAT_BASE', 'CAT_PREPCHANGE'}, {'CAT_BASE', 'CAT_AJCHANGE'};
%     {'CAT_BASE', 'CAT_PREPCHANGE'}, {'CAT_BASE', 'CAT_NCHANGE'};
%     {'CAT_BASE', 'CAT_AJCHANGE'}, {'CAT_BASE', 'CAT_NCHANGE'};
%     
%     {'PEN_BASE', 'PEN_MEANINGPRES'}, {'PEN_BASE', 'PEN_AJCHANGE'};
%     {'PEN_BASE', 'PEN_MEANINGPRES'}, {'PEN_BASE', 'PEN_NCHANGE'};
%     {'PEN_BASE', 'PEN_PREPCHANGE'}, {'PEN_BASE', 'PEN_AJCHANGE'};
%     {'PEN_BASE', 'PEN_PREPCHANGE'}, {'PEN_BASE', 'PEN_NCHANGE'};
%     {'PEN_BASE', 'PEN_AJCHANGE'}, {'PEN_BASE', 'PEN_NCHANGE'}
% 
%     {'GLASS_BASE', 'GLASS_MEANINGPRES'}, {'GLASS_BASE', 'GLASS_AJCHANGE'};
%     {'GLASS_BASE', 'GLASS_MEANINGPRES'}, {'GLASS_BASE', 'GLASS_NCHANGE'};
%     {'GLASS_BASE', 'GLASS_PREPCHANGE'}, {'GLASS_BASE', 'GLASS_AJCHANGE'};
%     {'GLASS_BASE', 'GLASS_PREPCHANGE'}, {'GLASS_BASE', 'GLASS_NCHANGE'};
%     {'GLASS_BASE', 'GLASS_AJCHANGE'}, {'GLASS_BASE', 'GLASS_NCHANGE'}};
% dissimilarities = {'BASE', 'MEANINGPRES', 1;
%     'BASE', 'PREPCHANGE', 1;
%     'BASE', 'AJCHANGE', 2;
%     'BASE', 'NCHANGE', 3};
% refRDM = searchlight_makeRefRDM(slConfig.pairsOfPairs, colNames.(slConfig.expt), dissimilarities);
% slConfig.refRDM = refRDM;   % should be in the alphabetical order of contrast names in the subject data (see colNames.mat)
% PhraseRSA_searchlight(dataPath, ssNames, colNames.(slConfig.expt), slConfig);

%% Move this to its own function (I wrote this for the Sam G meeting)
% p = [98,99];
% clusterCenters = cell(length(ssNames),1);
% clusterTaus = cell(length(ssNames),1);
% clusterM = cell(length(ssNames),1);
% for ss = 1:length(ssNames)
%     load(fullfile(dataPath, [ssNames{ss}, '_data']));
%     ssNum = ['0', num2str(ss)];
%     load(['searchlight_ss', ssNum(end-1:end)]);
%     thres = prctile(conTauVals,98);
%     voxelInds = conTauVals > thres;
%     currData = data.langloc.conVals(voxelInds,[3 6]);
%     
%     voxelCoords = data.voxelXYZ(voxelInds,:);
%     voxelDists = squareform(pdist(voxelCoords));
%     voxelSims = max(voxelDists(:)) - voxelDists;
%     voxelSims = voxelSims /  max(voxelSims(:));
%     voxelNames = cellfun(@num2str,num2cell(find(voxelInds)),'uniformoutput',false);
%     [clusters,Q] = clusteringHierarchical(voxelSims, voxelNames, 0);
%     clusters = clusters(:, Q==max(Q));
%     m = zeros(length(unique(clusters)),2);
%     se = zeros(length(unique(clusters)),2);
%     for c = 1:length(unique(clusters))
%         m(c,1) = mean(currData(clusters==c,1));
%         se(c,1) = std(currData(clusters==c,1))/sqrt(sum(clusters==c));
%         m(c,2) = mean(currData(clusters==c,2));
%         se(c,2) = std(currData(clusters==c,2))/sqrt(sum(clusters==c));
%     end
%     
%     figure(ss)
%     hold on
%     b = bar(m);
%     set(b(1),'facecolor',[1 0.3 0.3]);
%     set(b(2),'facecolor', [0.3 0.3 1]);
%     errorbar((1:size(m,1))'-0.1429, m(:,1), se(:,1), '-k', 'linestyle', 'none');
%     errorbar((1:size(m,1))'+0.1429, m(:,2), se(:,2), '-k', 'linestyle', 'none'); 
% 
%     clusterM{ss} = m;
%     clusterCenters{ss} = zeros(size(m,1),3);
%     clusterTaus{ss} = zeros(size(m,1),1);
%     yLims = get(gca,'ylim');    
%     for ii = 1:size(m,1)
%         xyz = round(mean(voxelCoords(clusters==ii,:),1));
%         clusterCenters{ss}(ii,:) = xyz;        
%         currTau = conTauVals(voxelInds);
%         clusterTaus{ss}(ii) = mean(currTau(clusters==ii));
% 
%         text(ii, yLims(1)+1.05*(yLims(2)-yLims(1)), [num2str(xyz(1)),',',num2str(xyz(2)),',',num2str(xyz(3))], ...
%             'horizontalalignment','center','verticalalignment','top','fontsize',12);
%         text(ii, yLims(1)+1.1*(yLims(2)-yLims(1)), num2str(round(100*mean(currTau(clusters==ii)))/100), ...
%             'horizontalalignment', 'center', 'verticalalignment','top','fontsize',12);        
%     end
%     set(gca,'ylim',[yLims(1) yLims(1)+1.2*(yLims(2)-yLims(1))]);
%     set(gca,'xlim',[0.5 find(Q==max(Q))+0.5]);
%     
%     if ss > 4
%         close
%     end
% end
% save myJunk
% 
% load myJunk
% clusterCenters = cell2mat(clusterCenters);
% clusterM = cell2mat(clusterM);
% megaDists = squareform(pdist(clusterCenters));
% megaSims = max(megaDists(:)) - megaDists;
% megaSims = megaSims /  max(megaSims(:));
% megaNames = cellfun(@num2str,num2cell((1:size(megaSims,1))'),'uniformoutput',false);
% [clusters,Q] = clusteringHierarchical(megaSims, megaNames, 0);
% clusters = clusters(:,find(Q==max(Q))+2);
% megaM = zeros(length(unique(clusters)),2);
% megaCenters = zeros(length(unique(clusters)),3);
% figure
% hold on
% for c = 1:length(unique(clusters))
%     megaM(c,:) = mean(clusterM(clusters==c,2),1);
%     megaCenters(c,:) = round(mean(clusterCenters(clusters==c,:)));
%     mCurr = clusterM(clusters==c,:);
%     nCurr = size(mCurr,1);
%     plot(repmat([c-0.1429; c+0.1429],1,nCurr), mCurr', 'o-k');    
%     scatter(repmat(c-0.1429,1,nCurr), mCurr(:,1)', 'o', 'markeredgecolor', [1 0.3 0.3], 'markerfacecolor', [1 0.3 0.3]);    
%     scatter(repmat(c+0.1429,1,nCurr), mCurr(:,2)', 'o', 'markeredgecolor', [0.3 0.3 1], 'markerfacecolor', [0.3 0.3 1]); 
%     plot(c+[-0.1429,+0.1429], [mean(mCurr(:,1)), mean(mCurr(:,2))], 'o-k', 'markerfacecolor', [0 0 0], 'linewidth', 3);
%     centerCurr = round(mean(clusterCenters(clusters==c,:)));
%     text(c, 1.1*max(mCurr(:)), [num2str(centerCurr(1)),',',num2str(centerCurr(2)),',',num2str(centerCurr(3))], ...
%         'horizontalalignment', 'center');    
% end    
% yLims = get(gca,'ylim');
% set(gca,'ylim',[yLims(1), 1.2*yLims(2)]);


%% OLD VERSION (REMOVED MAY 30, 2018): Reference RDMs for MVPA %%
% rankings = readtable('/Users/iblank/Desktop/MIT/Experiments/PhraseRSA/BehavioralData/rankings.csv', ...
%     'ReadVariableNames', true);
% nStim = length(stimSets);
% nTypes = length(stimTypes);
% 
% rdm = nan(nStim*(nTypes+1), nStim*(nTypes+1));
% colNames = cell(nStim*(nTypes+1),1);
% row = 1;
% for s = 1:nStim
%     setInds = rankings.set==stimSets(s);    
%     colNames{row} = [stimSetNames{s}, '_BASE'];
%     for t = 1:nTypes
%         typeInds = rankings.type==stimTypes{t,1};
%         data = rankings.ranking(setInds & typeInds);
%         if ~isempty(data)
%             m = mean(data);
%             rdm(row,row+t) = m;
%             rdm(row+t,row) = m;
%         end
%         colNames{row+t} = [stimSetNames{s}, '_', stimTypes{t,2}];
%     end
%     row = row + nTypes + 1;
% end
% rdm(1:(size(rdm,1)+1):end) = 0;
% refRDM.specific4StimSets = rdm;
% refRDM.colNames = colNames;
% 
% rdm = nan((nTypes+1), (nTypes+1), nStim);
% for s = 1:nStim
%     setInds = rankings.set==stimSets(s);    
%     for t = 1:nTypes
%         typeInds = rankings.type==stimTypes{t,1};
%         data = rankings.ranking(setInds & typeInds);
%         if ~isempty(data)
%             m = mean(data);
%             rdm(1,t+1,s) = m;
%             rdm(t+1,1,s) = m;
%         end
%     end
% end
% rdm = mean(rdm,3);
% rdm(1:(size(rdm,1)+1):end) = 0;
% rdm = squareform(tiedrank(squareform(rdm)));
% rdmOrig = rdm;
% for s = 1:(nStim-1)
%     rdm = blkdiag(rdm,rdmOrig);
% end
% rdm(isnan(refRDM.specific4StimSets)) = nan;
% refRDM.all30StimSets = rdm;
% 
% 
% rdm = nan((nTypes+1), (nTypes+1));
% pairsOfPairs = {1, [2, 4];
%     2, [];
%     3, 2;
%     4, 2};  % the item on the left is significantly more similar to BASE than the items on the right
% for p = 1:size(pairsOfPairs,1)
%     rdm(pairsOfPairs{p,2}, pairsOfPairs{p,1}) = 1;  % columnwise organization: each column has 1s for items that are more dissimilar to BASE
% end
% rdm(1:(size(rdm,1)+1):end) = 0;
% rdmOrig = rdm;
% for s = 1:(nStim-1)
%     rdm = blkdiag(rdm,rdmOrig);
% end
% rdm(isnan(refRDM.specific4StimSets)) = nan;
% refRDM.all30StimSets_sigPairs = rdm;
% 
% save referenceRDM refRDM 

%% OLD VERSION (REMOVED MAY 30, 2018): Null distribution for following analyses %%
% pairsPerSet = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
% pairs = [];
% subsets = cell(4,1);
% for ii = 1:4
%     pairs = [pairs; pairsPerSet+(4*(ii-1))];
%     subsets{ii} = ((ii-1)*4+1):(ii*4);
% end
% nullDist = rankCorr_Kendall_taua_null(refRDM(~isnan(refRDM)),pairs,subsets);
% save KendallNullDist nullDist


%% OLDEST VERSION: generating reference RDMs %%
% pairsOfPairs = {{'MAN_BASE', 'MAN_MEANINGPRES'}, {'MAN_BASE', 'MAN_AJCHANGE'};
%     {'MAN_BASE', 'MAN_MEANINGPRES'}, {'MAN_BASE', 'MAN_NCHANGE'};
%     {'MAN_BASE', 'MAN_PREPCHANGE'}, {'MAN_BASE', 'MAN_AJCHANGE'};
%     {'MAN_BASE', 'MAN_PREPCHANGE'}, {'MAN_BASE', 'MAN_NCHANGE'};
%     {'MAN_BASE', 'MAN_AJCHANGE'}, {'MAN_BASE', 'MAN_NCHANGE'};   
% 
%     {'CAT_BASE', 'CAT_MEANINGPRES'}, {'CAT_BASE', 'CAT_AJCHANGE'};
%     {'CAT_BASE', 'CAT_MEANINGPRES'}, {'CAT_BASE', 'CAT_NCHANGE'};
%     {'CAT_BASE', 'CAT_PREPCHANGE'}, {'CAT_BASE', 'CAT_AJCHANGE'};
%     {'CAT_BASE', 'CAT_PREPCHANGE'}, {'CAT_BASE', 'CAT_NCHANGE'};
%     {'CAT_BASE', 'CAT_AJCHANGE'}, {'CAT_BASE', 'CAT_NCHANGE'};
%     
%     {'PEN_BASE', 'PEN_MEANINGPRES'}, {'PEN_BASE', 'PEN_AJCHANGE'};
%     {'PEN_BASE', 'PEN_MEANINGPRES'}, {'PEN_BASE', 'PEN_NCHANGE'};
%     {'PEN_BASE', 'PEN_PREPCHANGE'}, {'PEN_BASE', 'PEN_AJCHANGE'};
%     {'PEN_BASE', 'PEN_PREPCHANGE'}, {'PEN_BASE', 'PEN_NCHANGE'};
%     {'PEN_BASE', 'PEN_AJCHANGE'}, {'PEN_BASE', 'PEN_NCHANGE'}
% 
%     {'GLASS_BASE', 'GLASS_MEANINGPRES'}, {'GLASS_BASE', 'GLASS_AJCHANGE'};
%     {'GLASS_BASE', 'GLASS_MEANINGPRES'}, {'GLASS_BASE', 'GLASS_NCHANGE'};
%     {'GLASS_BASE', 'GLASS_PREPCHANGE'}, {'GLASS_BASE', 'GLASS_AJCHANGE'};
%     {'GLASS_BASE', 'GLASS_PREPCHANGE'}, {'GLASS_BASE', 'GLASS_NCHANGE'};
%     {'GLASS_BASE', 'GLASS_AJCHANGE'}, {'GLASS_BASE', 'GLASS_NCHANGE'}};
% dissimilarities = {'BASE', 'MEANINGPRES', 1;
%     'BASE', 'PREPCHANGE', 2;
%     'BASE', 'AJCHANGE', 3;
%     'BASE', 'NCHANGE', 4};
% [refRDM, refColNames] = searchlight_makeRefRDM(pairsOfPairs, colNames.complang, dissimilarities);
% save referenceRDM refRDM pairsOfPairs refColNames
