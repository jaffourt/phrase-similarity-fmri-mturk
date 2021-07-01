function PhraseRSA_MVPA_relvox_locData(dataPath, ssNames, locConfig)

%%% Plot the correlations between reliability measures and localizer effect sizes across voxels;
%%% as well as the effect sizes themselves for the top most reliable voxels

%%% INPUT:
%%% dataPath = string, full path to where subjects' data are stored
%%% ssNames = cell of strings (data for subject i have been saved in: ssNames{i}_data.mat
%%% locConfig = structure with 6 fields:
%%%             expt = string, experiment name
%%%             conds = cell of strings, conditions in expt for which correlation should be computed
%%%                     (note: correlation is also computed for the difference between the first two conditions)
%%%             colors = n x 3 array of integers, each row is an RGB for a condition in conds 
%%%                      (+ additional row for the difference between the first two conditions)
%%%             chooseVoxels = cell of strings, reliability measures for choosing voxels
%%%                            (should be field names of the relVox variable in the <ssName>_relVox.mat file
%%%             nVoxels = integer array, number of top voxels to choose
%%%             vals = string, name of field in the data file where data are stored

%%% OUTPUT:
%%% text files for mixed-effects linear modeling in R

%% Parameters %%
nConds = numel(locConfig.conds);
nSs = numel(ssNames);

relMeas = locConfig.chooseVoxels;  % reliability measures
nArrOrig = locConfig.nVoxels;
nArrNames = cell(numel(nArrOrig)+1,1);
for nVox = 1:numel(nArrOrig)
    nArrNames{nVox} = ['top', num2str(nArrOrig(nVox))];
end
nArrNames{end} = 'allVoxels';

condNames = [locConfig.conds, {[locConfig.conds{1}, 'vs', locConfig.conds{2}]}];

dataType = 'indMed';
dispersionType = 'ci';


%% Compute correlations %%
disp('Correlation between voxel reliability and localizer effect size:');
fisherCorrs = struct;
conVals = struct;
for ss = 1:nSs
    disp(['  ', num2str(ss), '. ', ssNames{ss}]);
    load(fullfile(dataPath, [ssNames{ss}, '_data']));
    x = data.(locConfig.expt).(locConfig.vals);    
    colNames = data.(locConfig.expt).colNames;        
    nArr = [nArrOrig; numel(data.voxelInds)];
    load(fullfile(dataPath, [ssNames{ss}, '_relVox'])); 

    %% Loop through relability measures + number of top voxels %%
    for m = 1:length(relMeas)
        disp(['    reliability measure: ', relMeas{m}]);
        currMeas = relVox.(relMeas{m});

        for nVox = 1:length(nArr)    
            if ss == 1
                fisherCorrs.(['chooseBy', relMeas{m}]).(nArrNames{nVox}) = nan(nConds+1,nSs);
                conVals.(['chooseBy', relMeas{m}]).(nArrNames{nVox}) = nan(nConds+1,nSs);
            end

            thePrctile = min(100, 100*(nArr(nVox)/numel(data.voxelInds)));
            voxelInds = currMeas >= prctile(currMeas, 100-thePrctile);                    
            vals = nan(sum(voxelInds),nConds+1);

            %% Extract data for current subject %%
            for c = 1:nConds
                if c < 3
                    eval(['col', num2str(c), ' = col;']);
                end
                vals(:,c) = x(voxelInds,col);
            end            
            vals(:,nConds+1) = x(voxelInds,col1) - x(voxelInds,col2);
            rels = currMeas(voxelInds);
            
            nanRows1 = any(isnan(vals),2);
            nanRows2 = any(isnan(rels),2);
            goodRows = (~nanRows1) & (~nanRows2);            
            vals = vals(goodRows,:,:);           
            rels = rels(goodRows,:,:);
            
            fisherCorrs.(['chooseBy', relMeas{m}]).(nArrNames{nVox})(:,ss) = atanh(corr(vals,rels));
            conVals.(['chooseBy', relMeas{m}]).(nArrNames{nVox})(:,ss) = nanmean(vals,1)';
        end
    end
end
save RSA_reliableVoxels_locEffectSize fisherCorrs conVals condNames relMeas nVox


%% Plot %%
load RSA_reliableVoxels_locEffectSize
pdfNames = cell(2*numel(relMeas),1);
for m = 1:numel(relMeas)
    fData = nan(nSs, numel(nArrNames), nConds+1);
    pVals_f = zeros(numel(nArrNames)*nConds, 3);

    cData = nan(nSs, numel(nArrNames), nConds+1);
    pVals_c = zeros(numel(nArrNames)*nConds, 3);
    row = 1;
    
    for nVox = 1:numel(nArrNames)        
        fData(:,nVox,:) = permute(fisherCorrs.(['chooseBy', relMeas{m}]).(nArrNames{nVox}), [2 3 1]);
        cData(:,nVox,:) = permute(conVals.(['chooseBy', relMeas{m}]).(nArrNames{nVox}), [2 3 1]);
        
        for c = 1:(nConds+1)
            [~,~,~,s] = ttest((fisherCorrs.(['chooseBy', relMeas{m}]).(nArrNames{nVox})(c,:))');
            p = tcdf(s.tstat,s.df,'upper');
            pVals_f(row,:) = [nVox, c, p];
            
            [~,~,~,s] = ttest((conVals.(['chooseBy', relMeas{m}]).(nArrNames{nVox})(c,:))');
            p = tcdf(s.tstat,s.df,'upper');
            pVals_c(row,:) = [nVox, c, p];                        
            row = row + 1;
        end
    end
        
    qVals_f = FDRcorrect(pVals_f(:,3));
    pVals_f = pVals_f(qVals_f<0.05, 1:2);
    if ~isempty(pVals_f)
        pVals_f = [pVals_f, zeros(size(pVals_f,1),1)];
    else
        pVals_f = [];
    end
    
    qVals_c = FDRcorrect(pVals_c(:,3));
    pVals_c = pVals_c(qVals_c<0.05, 1:2);
    if ~isempty(pVals_c)
        pVals_c = [pVals_c, zeros(size(pVals_c,1),1)];
    else
        pVals_c = [];
    end        
    
    axes1 = plotData(fData, dataType, dispersionType, pVals_f, ...
        condNames, [relMeas{m} ' vs. effect size'], locConfig.colors, [0 numel(nArrNames)+1]);
    axes(axes1);
    theFigure = gcf;
    ylabel('Fisher correlation', 'fontname', 'calibri', 'fontsize', 12);
    set(axes1,'xtick',1:numel(nArrNames),'xticklabels',nArrNames,'fontname','calibri','fontsize',12);
    set(theFigure,'units','normalized','paperposition',[0 0 1 1]);
    orient('landscape');                
    print(theFigure, num2str((m-1)*2+1), '-dpdf');
    pdfNames{(m-1)*2+1} = [num2str((m-1)*2+1), '.pdf'];  

    pVals = zeros(numel(nArrNames)*nConds, 3);    
    axes2 = plotData(cData, dataType, dispersionType, pVals_c, ...
        condNames, ['Effect size of voxels chosen by ', relMeas{m}], locConfig.colors, [0 numel(nArrNames)+1]);
    axes(axes2);
    theFigure = gcf;
    ylabel('Effect size', 'fontname', 'calibri', 'fontsize', 12);
    set(axes2,'xtick',1:numel(nArrNames),'xticklabels',nArrNames,'fontname','calibri','fontsize',12);
    set(theFigure,'units','normalized','paperposition',[0 0 1 1]);
    orient('landscape');                
    print(theFigure, num2str(m*2), '-dpdf');
    pdfNames{m*2} = [num2str(m*2), '.pdf'];  
end
append_pdfs('RSA_reliableVoxels_locEffectSize.pdf', pdfNames{:});      
for fig = 1:length(pdfNames)
    delete(pdfNames{fig});
end               



%% Create a table for mixed effects modeling in R %%
% load RSA_reliableVoxels_locEffectSize
% nVoxelSets = numel(relMeas)*numel(nArrNames);
% nRows = nVoxelSets*nSs*(nConds+1);
% data = cell(nRows,3);   % columns: reliability measure, nVoxels, subject, condition,
%                         % effect size for each condition + correlation value for each condition
% row = 1;
% for m = 1:length(relMeas)
%     for nVox = 1:length(nArrNames)
%         vals = fisherCorrs.(['chooseBy', relMeas{m}]).(nArrNames{nVox});
%         for ss = 1:nSs
%             ssName = ['00', num2str(ss)];
%             ssName = ssName(end-2:end);
%             ssName = ['ss', num2str(ssName)];         
%             for c = 1:(nConds+1)
%                 data{row,1} = relMeas{m};
%                 data{row,2} = nArrNames{nVox};
%                 data{row,3} = ssName;     
%                 data{row,4} = condNames{c};
%                 data{row,5} = dataForPlot(ss, nVox, c, m);
%                 data{row,6} = vals(c,ss);
%                 row = row+1;
%             end
%         end
%     end
% end
% t = cell2table(data, 'variablenames', {'reliabilityMeasure', 'nVoxels', 'ID', 'condition', 'beta', 'rWithReliability'});
% save RSA_reliableVoxels_locEffectSize t -append
% writetable(t,'RSA_reliableVoxels_locEffectSize_forR.txt', 'delimiter', '\t');                          