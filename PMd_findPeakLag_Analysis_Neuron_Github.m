%% The following code illustrates how the peak lag was determined for each individual escape attempt, here
% for jumps in the CO2 assay.

% see 'findpeaks' on line 73 for calcium peak definition.

clear all

onlyFirstJumpOfBout = 1; %remove jumps in preceded or followed by another jump by 10 seconds.
onlyCategorizedCells = 1; %use only cells in time lag analysis that were categorized as behavior cells.

preIdx = 7.5 .* 10;
postIdx = 7.5 .* 10;

%CO2 sessions
folders{1} = 'G:\PMD_Miniscope\668\CO2\2'; %more jump epochs ... 2, but better than 1
folders{2} = 'G:\PMD_Miniscope\665\CO2\2';
folders{3} = 'G:\PMD_Miniscope\663\CO2';
folders{4} = 'G:\PMD_Miniscope\20a\CO2\exposure_6_best\H17_M10_S10';
folders{5} = 'G:\PMD_Miniscope\611\CO2\2_best';
folders{6} = 'G:\PMD_Miniscope\612\CO2\2_best';
folders{7} = 'G:\PMD_Miniscope\613\CO2\2_best';
folders{8} = 'G:\PMD_Miniscope\616\CO2\2_best';
folders{9} = 'G:\PMD_Miniscope\617\CO2\4_best';
folders{10} = 'G:\PMD_Miniscope\618\CO2\1_best';
folders{11} = 'G:\PMD_Miniscope\620\CO2\1_best';
folders{12} = 'G:\PMD_Miniscope\822\CO2\exposure_2_best\H18_M58_S45';

for mouseNum = 1:length(folders)
   
    cd(folders{mouseNum})
    
    %find all cell indexes from all jump and climb hose categorization
    if exist('jumpSeg.mat','file')
        load('jumpSeg.mat'); load('jumpSegPreBehav.mat'); %load('jumpSegPreOnlyBehav.mat');
    end
    load('good_neurons.mat')
    load('output_CNMF-E.mat','neuron'); nrn = neuron.C_raw; clearvars neuron;
    load('BehaviorMS.mat','climbHoseFrameMS','jumpFrameMS')
    
    if onlyFirstJumpOfBout==1
        if exist('jumpFrameMS')
            jumpDiffRemoveIdx = diff(jumpFrameMS(:,1));
            jumpDiffRemoveIdx = find(jumpDiffRemoveIdx < 75); jumpDiffRemoveIdx = jumpDiffRemoveIdx + 1;
            jumpFrameMS(jumpDiffRemoveIdx,:) = [];
        end
    end
    
    cellIdx_jump = unique([find(jumpSeg==1) find(jumpSegPreBehav==1)]);
        
    cellIdx = find(good_neurons);
    
    clearvars coeffSeg coeffSeg_shuff jumpSeg jumpSegPreBehav jumpSegPreOnlyBehav climbHoseSeg climbHoseSegPreBehav climbHoseSegPreOnlyBehav
    
    %z-score calcium traces
    for cellNum = 1:size(nrn,1)
       nrn(cellNum,:) = zscore(nrn(cellNum,:)); 
    end
    
    %get mean activity for each good neuron
    for cellNum = 1:length(cellIdx)  
        for behavNum = 1:size(jumpFrameMS,1)
            if jumpFrameMS(behavNum,1) > preIdx
                jumpSig{mouseNum,cellIdx(cellNum)}(behavNum,:) = nrn(cellIdx(cellNum),jumpFrameMS(behavNum,1)-preIdx:jumpFrameMS(behavNum,1)+postIdx);
            else
                jumpSig{mouseNum,cellIdx(cellNum)}(behavNum,:) = nan(1,preIdx+postIdx+1);                
            end
        end
 
        % find the peak of each individual behavioral trial and save as variable
        key = linspace(-10,10,151);

            for behavNum = 1:size(jumpSig{mouseNum,cellIdx(cellNum)},1)
                tempSig = smoothdata(jumpSig{mouseNum,cellIdx(cellNum)}(behavNum,:));
                [pks,locs,w,p] = findpeaks(tempSig, 'MinPeakProminence',.4); %findpeaks(tempSig, 'MinPeakProminence',.4,'Annotate','extents')
                [temp, index] = max(p);
                lag = locs(index);
                lag = key(lag);
                if ~isempty(lag)
                    lag_jump{mouseNum,cellIdx(cellNum)}(behavNum) = lag;
                else
                    lag_jump{mouseNum,cellIdx(cellNum)}(behavNum) = nan;                
                end
            end
                lag_jump_se(mouseNum,cellIdx(cellNum)) = nanstd(lag_jump{mouseNum,cellIdx(cellNum)});
                lag_jump_mean(mouseNum,cellIdx(cellNum)) = nanmean(lag_jump{mouseNum,cellIdx(cellNum)});   
                
        jumpSigAll{mouseNum,cellIdx(cellNum)} = jumpSig{mouseNum,cellIdx(cellNum)};        
        jumpSig{mouseNum,cellIdx(cellNum)} = nanmean(jumpSig{mouseNum,cellIdx(cellNum)});
        
    end
    
    for cellNum = 1:length(cellIdx_jump)
        jumpSig_GLMAll_CO2{mouseNum,cellIdx_jump(cellNum)} = jumpSigAll{mouseNum,cellIdx_jump(cellNum)};
        jumpSig_GLM_CO2{mouseNum,cellIdx_jump(cellNum)} = jumpSig{mouseNum,cellIdx_jump(cellNum)};
    end
        
    if onlyCategorizedCells==1
        temp = zeros(1,length(lag_jump_mean(mouseNum,:))); temp(cellIdx_jump) = 1;
        lag_jump_mean(~temp) = 0;
        lag_jump_se(~temp) = 0;
        
        idxToKeep = find(temp);
        for cellNum = 1:length(idxToKeep)
            lag_jump_cat{mouseNum,idxToKeep(cellNum)} = lag_jump{mouseNum,idxToKeep(cellNum)};
        end
    end 
end

%compile all cell activity into matrix

cntr = 1;
for mouseNum = 1:size(jumpSig,1)
    for cellNum = 1:size(jumpSig,2)
        
        if isempty(jumpSig{mouseNum,cellNum})
            continue
        end
        
        jumpAll(cntr,:) = jumpSig{mouseNum,cellNum};
        cntr = cntr+1;
    end
end


[what, I_jump] = max(jumpAll,[],2);

[what,sortJump] = sort(I_jump);

%convert each matrix to min/max
jumpAllSort = jumpAll(sortJump,:);
for cellNum = 1:size(jumpAllSort,1)    
    minVal = min(jumpAllSort(cellNum,:));
    maxVal = max(jumpAllSort(cellNum,:));
    jumpAllSort(cellNum,:) = (jumpAllSort(cellNum,:) - minVal) / ( maxVal - minVal );
end

A = subplot(1,5,1)
imagesc(jumpAllSort)
title(['CO2: n=', num2str(size(jumpAllSort,1))])
%caxis(vals)
A.XTick = [1,preIdx+1,preIdx+postIdx+1]; 
A.XTickLabel = [-10,0,10];
xlabel('time (s)')
hold on;
plot([preIdx+1 preIdx+1],[0 size(jumpAll,1)],'r')


jumpSig_CO2 = jumpSig; 
climbHoseSig_CO2 = climbHoseSig;
jumpAll_CO2 = jumpAll;
climbHoseAll_CO2 = climbHoseAll;
