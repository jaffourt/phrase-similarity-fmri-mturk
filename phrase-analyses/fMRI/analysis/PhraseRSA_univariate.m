function PhraseRSA_univariate(dataPath, ssNames, uniConfig)

%%% Univariate analysis: response magnitude to different experimental conditions in different regions

%%% INPUT:
%%% dataPath = string, full path to where subjects' data are stored
%%% ssNames = cell of strings (data for subject i have been saved in: ssNames{i}_data.mat
%%% uniConfig = structure; for each experiment (field), 4 embedded fields:
%%%             cols = string cell, columns to analyze from colNames
%%%             avg = cell of strings; for each entry, data will be averaged across conditions
%%%                   containing that string
%%%             hemiSplit = array; for each cortical system, 1 to plot first half of masks separately
%%%                         from second half (i.e., split by hemisphere)
%%%             isloc = cell of strings, names of systems for which this
%%%                     experiment was used as a localizer
%%%             vals = tVals' or 'conVals', the data from which localizer contrast is computed 

%%% OUTPUT:
%%% Plots of univariate responses for each condition of each experiment, in each  mask,
%%% as well as data averaged across conditions (see above) and across masks

%%% NOTE: for the localizer experiment, data are automatically cross validated

%% Parameters %%
dataType = 'indMean';
dispersionType = 'ci';
expNames = fieldnames(uniConfig);

%% Extract data %%
disp('Univariate analysis:');
univarData = struct;
for ss = 1:length(ssNames)
    disp([num2str(ss), '. ', ssNames{ss}]);
    load(fullfile(dataPath, [ssNames{ss}, '_data']));
    load(fullfile(dataPath, [ssNames{ss}, '_rois']));
    if ss == 1
        fROIs = fieldnames(roiInds);
        systems = fieldnames(roiInds.(fROIs{1}));
    end
    v = data.voxelInds;                                         % all GM voxel coordinates for current subject
    
    for e = 1:length(expNames)
        x = data.(expNames{e}).(uniConfig.(expNames{e}).vals);
        colNames = data.(expNames{e}).([uniConfig.(expNames{e}).vals, 'ColNames']);
        colInds = cellfun(@(x)(find(strcmp(colNames,x))), uniConfig.(expNames{e}).cols);
        
        %% Loop through fROI definition criteria + cortical systems %%
        for f = 1:length(fROIs)
            for s = 1:length(systems)
                
                %% Should cross-validation be used? %%
                if sum(strcmp(systems{s}, uniConfig.(expNames{e}).isloc)) > 0
                    doCV = 1;
                else
                    doCV = 0;
                end
                
                %% Loop through masks %%
                nROI = length(roiInds.(fROIs{f}).(systems{s}));  
                if ss == 1
                    nCols = length(uniConfig.(expNames{e}).avg);
                    if nCols==0
                        nCols = length(colInds);
                    end
                    valsPlot = nan(length(ssNames), nROI, nCols);    % subjects x ROIs x conditions
                else
                    valsPlot = univarData.(expNames{e}).(fROIs{f}).(systems{s});
                end
                
                for r = 1:nROI
                    if ~doCV
                        inds = roiInds.(fROIs{f}).(systems{s})(r).all;
                        [~, ~, indsInV] = intersect(inds, v);
                        valsCurr = [];
                        for c = 1:length(uniConfig.(expNames{e}).avg)
                            currCols = ~cellfun(@isempty, strfind(uniConfig.(expNames{e}).cols, uniConfig.(expNames{e}).avg{c})); % all columns from condition c                            
                            valsCurr(1,c) = nanmean(nanmean(x(indsInV,colInds(currCols)),1),2);                        
                        end
                        if isempty(valsCurr)                 % in case no averaging is required
                            valsCurr = nanmean(x(indsInV, colInds),1);   
                        end                        

                    else   % cross validation
                        indsEven = roiInds.(fROIs{f}).(systems{s})(r).even;
                        [~, ~, evenInV] = intersect(indsEven, v);
                        indsOdd = roiInds.(fROIs{f}).(systems{s})(r).odd;
                        [~, ~, oddInV] = intersect(indsOdd, v);
                        valsCurr = [];

                        for c = 1:length(uniConfig.(expNames{e}).avg)
                            condCols = ~cellfun(@isempty, strfind(colNames, uniConfig.(expNames{e}).avg{c}));
                            oddCols = ~cellfun(@isempty, strfind(colNames, 'ODD'));
                            valsEven = x(evenInV, (condCols & oddCols));
                            
                            evenCols = ~cellfun(@isempty, strfind(colNames, 'EVEN'));
                            valsOdd = x(oddInV, (condCols & evenCols));
                            
                            valsCurr(1,c) = nanmean(nanmean(valsEven,1)/2 + nanmean(valsOdd,1)/2, 2);
                        end
                        if isempty(valsCurr)                 % in case no averaging is required
                            for c = 1:length(uniConfig.(expNames{e}).cols)
                                oddCol = strcmp(colNames, ['ODD_', uniConfig.(expNames{e}).cols{c}]);
                                valsEven = x(evenInV, oddCol);

                                evenCol = strcmp(colNames, ['EVEN_', uniConfig.(expNames{e}).cols{c}]);
                                valsOdd = x(oddInV, evenCol);
                                
                                valsCurr(1,c) = nanmean(nanmean(valsEven,1)/2 + nanmean(valsOdd,1)/2, 2);
                            end
                        end
                    end
                    valsPlot(ss,r,:) = permute(valsCurr, [1 3 2]);                    
                end
                univarData.(expNames{e}).(fROIs{f}).(systems{s}) = valsPlot;
            end
        end
    end
    save univariateResults univarData
end

save univariateResults univarData
        
%% Plot data %%
disp(['Warning: the plotting part of this code does not work on VNC; Idan Blank, August, 2018']);
load univariateResults
load(fullfile(dataPath, [ssNames{1}, '_data']));        % for colNames
for e = 1:length(expNames)
    fROIs = fieldnames(univarData.(expNames{e}));
    systems = fieldnames(univarData.(expNames{e}).(fROIs{1}));
    
    condLabels = uniConfig.(expNames{e}).avg;
    if isempty(condLabels)
        condLabels = uniConfig.(expNames{e}).cols;
    end
    nConds = length(condLabels);
    colors = repmat(linspace(0, 0.75, nConds)', 1, 3);
    
    for f = 1:length(fROIs)
        figNames = {};
        figNum = 1;        
        
        for s = 1:length(systems)               
            vals = univarData.(expNames{e}).(fROIs{f}).(systems{s});

            if ~uniConfig.(expNames{e}).hemiSplit(s)        
                theTitle = [expNames{e}, ': ', systems{s}, ' ROIs, voxel selection (%): ', fROIs{f}];            
                
                plotData(vals, dataType, dispersionType, [], condLabels, {theTitle; '(Contrast values)'}, colors, []);
                orient('landscape');
                print(gcf, num2str(figNum), '-dpdf');
                figNames{figNum} = [num2str(figNum), '.pdf'];              
                figNum = figNum+1;
                
            else
                n = size(vals,2)/2;
                theTitle = [expNames{e}, ': LH ', systems{s}, ' ROIs, voxel selection (%): ', fROIs{f}];                

                plotData(vals(:,1:n,:), dataType, dispersionType, [], condLabels, {theTitle; '(Contrast values)'}, colors, []);
                orient('landscape');                
                print(gcf, num2str(figNum), '-dpdf');
                figNames{figNum} = [num2str(figNum), '.pdf'];              
                figNum = figNum+1;
                
                theTitle = [expNames{e}, ': RH ', systems{s}, ' ROIs, voxel selection (%): ', fROIs{f}];
                plotData(vals(:,(n+1):end,:), dataType, dispersionType, [], condLabels, {theTitle; '(Contrast values)'}, colors, []);
                orient('landscape');                
                print(gcf, num2str(figNum), '-dpdf');
                figNames{figNum} = [num2str(figNum), '.pdf'];              
                figNum = figNum+1;               
            end
        end
        
        %% Save %%
        append_pdfs(['Univar_', expNames{e}, '_', fROIs{f}, '.pdf'], figNames{:});      
        for fig = 1:(figNum-1)
            delete(figNames{fig});
        end       
        close all
    end
end