function f = plotData(data, dataType, dispersionType, pVals, withinGroup, theTitle, colors, xLims)

%%% This function plots the contrast estimates for the the tasks of Experiments 1 & 2

%%% INPUT:
%%% data = 3D martix of numbers, subjects X grouping levels (usually, ROIs) X within group levels (usually, conditions)
%%% dataType = string, from among:
%%%            'bar': bars of sample means
%%%            'indMean': means + individual data points
%%%            'indMed': medians + individual data points 
%%% dispersionType = string, from among:
%%%                  'se': standard errors
%%%                  'sd': standard deviations
%%%                  'ci': 95% confidence intervals
%%% pVals = m x 3 array of significant results (can be empty)
%%%         col1: grouping level
%%%         col2,3: within group levels
%%%                 for single level, do: [level, 0]
%%%                 for pairwise comparisons: [level1, level2]
%%%                 for interactions: [(level1+level2)/2, (level3+level4)/2]
%%% withinGroup = 1 x n cell of strings, with names for each within-group level
%%% theTitle = string with title of experiment
%%% colors = n X 3 matrix of RGB values (in 0-1 range), one per withinGroup value
%%% xLims = 1X2 numerical array, for controlling the width of the x-axis
%%%         (if you want to make sure that several different figures have
%%%         bars of the same width; otherwise, can be an empty array)

%%% Idan Blank, 03.05.16; edited 03.22.2017; edited 09.29.2017

nWithin = size(data,3);
nGroups = size(data,2);
nSs = size(data,1);
figure

%% Extract means and dispersion measures %%
switch dataType
    case 'bar'
        m = permute(nanmean(data,1),[2 3 1]);
    case 'indMean'
        m = permute(nanmean(data,1),[2 3 1]);
    case 'indMed'
        m = permute(prctile(data,50,1),[2 3 1]);
end
m(isnan(m)) = 0;

switch dispersionType
    case 'se'
        dispData = permute(nanstd(data,0,1)./sqrt(sum(~isnan(data),1)), [2 3 1]);
        dispData = cat(3, m-dispData, m+dispData);
        widthWeight = 0.5;
    case 'sd'
        dispData = permute(nanstd(data,0,1), [2 3 1]);
        dispData = cat(3, m-dispData, m+dispData);  
        widthWeight = 0.5;
    case 'ci'
        nPerms = 1000;
        dispData = zeros(nGroups, nWithin, 2);
        for g = 1:nGroups
            for w = 1:nWithin
                switch dataType
                    case 'indMean'
                        ci = singleSamplePermCI(data(:,g,w), nPerms, 'mean');
                    case 'bar'
                        ci = singleSamplePermCI(data(:,g,w), nPerms, 'mean');
                    case 'indMed'
                        ci = singleSamplePermCI(data(:,g,w), nPerms, 'med');
                end
                dispData(g,w,1) = ci(1);
                dispData(g,w,2) = ci(2);
            end
        end
        widthWeight = 1;
end

if strcmp(dataType, 'bar')
    dispData(:,:,1) = m-dispData(:,:,1);    % turn back to dispersion data NOT relative to the mean
    dispData(:,:,2) = dispData(:,:,2)-m;
end      

dispData(isnan(dispData)) = 0;

%% Generate plot %%
m
dispData

if size(m,1) > 1
    b = bar(zeros(size(m)));        % just to get bar width + offsets
else
    b = bar(repmat(m,2,1));
end

halfWidth = 0.5*b(1).BarWidth*(b(2).XOffset-b(1).XOffset);
offsets = zeros(length(b),1);
for i = 1:length(b)
    offsets(i) = b(i).XOffset;
end

clf reset
% set(gcf, 'paperunits', 'centimeters', 'position', [0 0 3.5 3.5]);
set(gcf, 'paperunits', 'normalized', 'position', [0 0 1 1]);
set(gca, 'fontname', 'helvetica', 'fontsize', 8);
if size(m,1) > 1
    b = bar(zeros(size(m)));
else
    b = bar(zeros(2,size(m,2)));
end
for w = 1:nWithin
    set(b(w),'FaceColor',colors(w,:));
end

legend(cellfun(@regexprep, withinGroup, repmat({'_'},1,nWithin), repmat({' '},1,nWithin),'uniformoutput',false),...
    'location', 'northeast','box','off');
hold on

if strcmp(dataType, 'bar')
    if size(m,1) > 1
        b = bar(m);
    else
        b = bar([m; zeros(1,size(m,2))]);
    end
    for w = 1:nWithin
        b(w).FaceColor = colors(w,:);
    end

    allOffsets = repmat(offsets', nGroups, 1);        
    errorbar(repmat((1:nGroups)',1,nWithin)+allOffsets, m, dispData(:,:,1), dispData(:,:,2), '-k', 'linestyle', 'none');

    ySimple = m + dispData(:,:,2);      % for p-values of simple effects
    
else
    allOffsets = repmat(offsets', nGroups, 1);
    xLows = repmat(1:nGroups,1,nWithin)+(allOffsets(:))'-halfWidth*widthWeight;
    xHighs = repmat(1:nGroups,1,nWithin)+(allOffsets(:))'+halfWidth*widthWeight;
    xRand = repmat(1:nGroups, nSs, 1)+repmat((1-2*rand(nSs,1))*halfWidth*widthWeight*0.3,1,nGroups);

    d1 = dispData(:,:,1);
    d1 = (d1(:))';
    d2 = dispData(:,:,2);
    d2 = (d2(:))';

    f = fill([xLows; xLows; xHighs; xHighs; xLows], ...
        [d1; d2; d2; d1; d1], ...
        [0 0 0], 'facealpha', 0.3, 'edgecolor', 'none');     % dispersion        
    for w = 1:nWithin
        plot(xRand+offsets(w), data(:,:,w), 'o', 'markeredgecolor','none','markerfacecolor',[0 0 0],'markersize',2.5);   
    end

    p = plot([xLows; xHighs], repmat((m(:))',2,1), '-k', 'linewidth', 2);
    for w = 1:nWithin
        firstInd = (w-1)*nGroups + 1;
        lastInd = w*nGroups;
        set(f(firstInd:lastInd), 'FaceColor', colors(w,:));
        set(p(firstInd:lastInd), 'Color', colors(w,:));
    end        

    ySimple = permute(max(data,[],1), [2 3 1]);
end

%% Cosmetics %%
box off
if isempty(xLims)
    set(gca, 'xlim', [0, nGroups+1], 'xtick', []);
else
    set(gca, 'xlim', xLims, 'xtick', []);    
end

title(theTitle, 'fontname', 'calibri', 'fontsize', 12);
switch dataType
    case 'indMean'
        yMin = min(min(1.15*data(:)), 0);
        yMax = max(max(1.15*data(:)), 0.05);
    case 'indMed'
        yMin = min(min(1.15*data(:)), 0);
        yMax = max(max(1.15*data(:)), 0.05);        
    case 'bar'
        d1 = dispData(:,:,1);
        d2 = dispData(:,:,2);
        yMin = min(min(1.15*(m(:)-d1(:))), 0);
        yMax = max(max(1.15*(m(:)+d2(:))), 0.05);
end
set(gca, 'ylim', [yMin yMax]);


%% Significance %%
if ~isempty(pVals)
    stepUp = 0.015*(yMax-yMin);
    yMax = max(ySimple(:))+7*stepUp;
    set(gca, 'ylim', [yMin yMax]);
    
    for i = 1:size(pVals,1)
        g = pVals(i,1);
        w1 = pVals(i,2);
        w2 = pVals(i,3);
        if w2 == 0
            plot([g+offsets(w1)-halfWidth*widthWeight, g+offsets(w1)+halfWidth*widthWeight], repmat(max(ySimple(g,:))+stepUp,1,2), '-', ...
                'color', colors(w1,:), 'linewidth', 0.5);
            text(g+offsets(w1), max(ySimple(g,:))+1.6*stepUp, '*', ...
                'fontname', 'helvetica', 'fontsize', 8, 'color', colors(w1,:), 'horizontalalignment', 'center', 'verticalalignment', 'cap');
        else
            if w2 == round(w2)
                x = g+offsets([w1 w2]);
                y = repmat(max(ySimple(g,:))+3*stepUp,1,2);
                yText = max(ySimple(g,:))+3.6*stepUp;
                theLine = 1;
                theStar = 9;
            else
                x1 = g + offsets(w1-0.5) + (offsets(w1+0.5)-offsets(w1-0.5))*0.5;
                x2 = g + offsets(w2-0.5) + (offsets(w2+0.5)-offsets(w2-0.5))*0.5;   
                x = [x1 x2];
                y = repmat(max(ySimple(g,:))+5*stepUp,1,2);
                yText = max(ySimple(g,:))+5.6*stepUp;
                theLine = 1.5;
                theStar = 10;
            end
            
            plot(x, y, '-k', 'linewidth', theLine);
            text(x(1)+(x(2)-x(1))*0.5, yText, '*', ...
                'fontname', 'helvetica', 'fontsize', theStar, 'horizontalalignment', 'center', 'verticalalignment', 'cap');            
        end
    end
end


%% Add axis labels %%
ylabel('Beta weight', 'fontname', 'calibri', 'fontsize', 12);
goodInds = isstrprop(theTitle, 'alphanum');     % letters or digits
theTitle(~goodInds) = '_';
isStop = 0;
while ~isStop
    newTitle = regexprep(theTitle, '__', '_');
    if strcmp(newTitle, theTitle)
        isStop = 1;
    else
        theTitle = newTitle;
    end
end
saveas(gcf, theTitle, 'epsc');
f = gca;