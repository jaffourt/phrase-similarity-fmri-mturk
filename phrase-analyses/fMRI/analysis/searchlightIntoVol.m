
%% Parameters %%
dataPath = '/mindhive/evlab/u/iblank/Desktop/Projects/PhraseRSA';
load(fullfile(dataPath, 'subjectInfo'));
load('KendallNullDist');
subjects = subjects.main;
nSs = length(subjects);

pThres = 0.001;

theSize = [79 95 69];
nVoxels = prod(theSize);
vol = cell(nVoxels,4);        % columns: 1=voxel, 2=cell of Kendall tau values; 3=n total, 4=n significant

%% Create volumes %%
% for ii = 1:nSs
%     disp(num2str(ii));
%     ssNum = ['0', num2str(ii)];    
%     load(['searchlight_ss', ssNum(end-1:end)]);
%     corrs = conTauVals;
%     load(fullfile(dataPath, [subjects{ii}, '_data']));
%     voxels = data.voxelInds;
%     clear data
%     
%     goodInds = ~isnan(corrs);
%     corrs = corrs(goodInds);
%     voxels = voxels(goodInds);
%     
%     for v = 1:length(voxels)
%         if mod(v,1000) == 0
%             disp(['  ', num2str(v)]);
%         end
%                 
%         ind = voxels(v);
%         if ~isempty(vol{ind,1})
%             vol{ind,2}(ii) = corrs(v);
%         else
%             vol{ind,1} = ind;
%             vol{ind,2} = nan(nSs,1);
%             vol{ind,2}(ii) = corrs(v);
%             vol{ind,3} = 0;
%             vol{ind,4} = 0;
%         end
%         
%         vol{ind,3} = vol{ind,3} + 1;
%         pVal = sum(nullDist(nullDist(:,1)>=corrs(v),2));
%         if pVal <= pThres
%             vol{ind,4} = vol{ind,4} + 1;
%         end            
%     end
% end
% 
% goodInds = ~cellfun(@isempty, vol(:,1));
% vol = vol(goodInds,:);
% medVol = nan(theSize);
% nVol = nan(theSize);
% for v = 1:size(vol,1)
%     if vol{v,3} >= 0.5*nSs
%         medVol(vol{v,1}) = prctile(vol{v,2},50);
%     end
%     nVol(vol{v,1}) = vol{v,4};
% end
% 
% save searchlightResultsVol medVol nVol


%% Create parcels %%
%%% Taken from spm_ss_estimate_ROI.m
addpath('/users/evelina9/fMRI_PROJECTS/spm_ss_vE/');
load searchlightResultsVol
goodInds = find(medVol(:)>0.05);

imgFile = spm_vol(fullfile('/mindhive/evlab/u/Shared/SUBJECTS', subjects{1}, ...
    'firstlevel_complang_spat_matching/con_0001.img'));
medFile = struct('fname','medianTauVolume.img','descrip','median Kendall tau (A)',...
    'mat',imgFile.mat,'dim',imgFile.dim,'dt',[spm_type('float32'),spm_platform('bigend')]);
medVol = spm_write_vol(medFile,medVol);

spm_smooth(medVol,'medianTauSmoothed.img',6*[1,1,1]); % 6 = FWHM for smoothing kernel
medVolSmooth = spm_vol('medianTauSmoothed.img');
data = spm_read_vols(medVolSmooth);
figure(1)
hist(data(goodInds),20)
pause

% for ii = 1:3:size(data,3)
%     imagesc(data(:,:,ii))
%     caxis([-0.1 0.15])
%     colorbar
%     pause(1)
% end

parcelsRSA = spm_ss_watershed(-data,find(data > 0.11));
fprintf('Done. Defined %d regions\n',max(parcelsRSA(:)));
finalVol = struct('fname','allParcels_RSA.img','mat',medVolSmooth.mat,'dim',medVolSmooth.dim,...
    'dt',[spm_type('int16') spm_platform('bigend')],'pinfo',[1;0;0]);
spm_write_vol(finalVol,parcelsRSA);
return


%% Plot %%
t1= spm_vol('/mindhive/evlab/u/iblank/Desktop/MNI_Parcels/T1_resizedToEPI.nii');
t1 = spm_read_vols(t1);
t1 = permute(t1,[3,2,1]);
t1 = flipdim(t1,1);

load searchlightResultsVol
nVol = permute(nVol,[3,2,1]);
nVol = flipdim(nVol,1);
medVol = permute(medVol, [3,2,1]);
medVol = flipdim(medVol,1);
% caxisMed = [min(medVol(:)), max(medVol(:))];
% caxisN = [0, max(nVol(:))];

c = colormap('hot');
c = c(3:end,:);
medThres = linspace(min(medVol(:)), max(medVol(:)), size(c,1)+1);
nThres = linspace(min(nVol(:)), max(nVol(:)), size(c,1)+1);

for ii = 1:ceil(size(nVol,3)/2)
    figure(1)
    clf reset
    theLayer = (ii-1)*2+1;
    for jj = 1:4
        subplot(2,2,jj)
        anat = t1(:,:,theLayer);
        if mod(jj,2) == 0
            data = medVol(:,:,theLayer);
            theLayer = theLayer + 1;
            thres = medThres;
        else
            data = nVol(:,:,theLayer);
            thres = nThres;
        end
        
        if any(~(data(:)==0))
            im = repmat(anat,1,1,3);
            r = im(:,:,1);
            g = im(:,:,2);
            b = im(:,:,3);
            
            goodInds = ~(data==0);    
            colorInds = cellfun(@(x)(find(x<=thres,1)), num2cell(data(goodInds)),'uniformOutput',false);
            colorInds = cell2mat(colorInds)-1;
            colorInds(colorInds==0) = 1;
            
            r(goodInds) = c(colorInds,1);
            g(goodInds) = c(colorInds,2);
            b(goodInds) = c(colorInds,3);
            im = cat(3,r,g,b);
        else
            im = repmat(anat,1,1,3);
        end
        
        imagesc(im);
%         if mod(jj,2)==0
%             caxis(caxisMed)
%         else
%             caxis(caxisN);
%         end
       colormap(c);
       colorbar('YTick',linspace(0,1,6),'YTicklabel',thres(round(linspace(1,length(thres),6))));
    end
    pause
end