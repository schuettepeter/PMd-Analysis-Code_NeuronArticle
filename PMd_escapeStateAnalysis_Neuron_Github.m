%% THE FOLLOWING IS CODE TO ACCOMPANY THE NEURON ARTICLE:
%"Coordination of escape and spatial navigation circuits orchestrate
%versatile flight from threats," Figures 3G-H and S6 A-B 

%Can an escape 'state' predict escape in each assay?

clear all; close all;

%choose whether to form clusters by k means or hmm
hmmAnalysis=1;
kMeansAnalysis=0;

zscore_CA_activity = 1;
doPCA=1;
varThresh = 80;
maxClusterNum = 2; %for k-means only
hmmClusterNum = 2; %for HMM only
bottomPC = 1; %remove bottom PCs if desired.
iter = 1000; %bootstrap iterations
negEscapeLag = 0;

%Rat Climb sessions
folders{1,1} = 'G:\PMD_Miniscope\668\RatClimb';
folders{2,1} = 'G:\PMD_Miniscope\665\RatClimb\2';
folders{3,1} = 'G:\PMD_Miniscope\663\RatClimb';
folders{5,1} = 'G:\PMD_Miniscope\611\RatClimb\2_best';
folders{6,1} = 'G:\PMD_Miniscope\612\RatClimb\2_best';
folders{7,1} = 'G:\PMD_Miniscope\613\RatClimb\2_best';
folders{8,1} = 'G:\PMD_Miniscope\616\RatClimb\2_best';
folders{9,1} = 'G:\PMD_Miniscope\617\RatClimb\1';
folders{11,1} = 'G:\PMD_Miniscope\620\RatClimb\1_best';
folders{12,1} = 'G:\PMD_Miniscope\822\RatClimb\exposure_2_best\H18_M4_S58';

%Hotplate sessions
folders{1,2} = 'G:\PMD_Miniscope\668\Hotplate\2';
folders{2,2} = 'G:\PMD_Miniscope\665\Hotplate\1';
folders{3,2} = 'G:\PMD_Miniscope\663\Hotplate\1';
folders{4,2} = 'G:\PMD_Miniscope\20a\Hotplate\H17_M58_S22';
folders{5,2} = 'G:\PMD_Miniscope\611\Hotplate\1_best';
folders{6,2} = 'G:\PMD_Miniscope\612\Hotplate\1_best';
folders{7,2} = 'G:\PMD_Miniscope\613\Hotplate\1';
folders{8,2} = 'G:\PMD_Miniscope\616\Hotplate\4_best';
folders{11,2} = 'G:\PMD_Miniscope\620\Hotplate\1';

%CO2 sessions
folders{1,3} = 'G:\PMD_Miniscope\665\CO2\2';
folders{2,3} = 'G:\PMD_Miniscope\663\CO2';
folders{3,3} = 'G:\PMD_Miniscope\20a\CO2\exposure_6_best\H17_M10_S10';
folders{4,3} = 'G:\PMD_Miniscope\611\CO2\2_best';
folders{5,3} = 'G:\PMD_Miniscope\612\CO2\2_best';
folders{6,3} = 'G:\PMD_Miniscope\613\CO2\2_best';
folders{7,3} = 'G:\PMD_Miniscope\616\CO2\2_best';
folders{8,3} = 'G:\PMD_Miniscope\617\CO2\4_best';
folders{9,3} = 'G:\PMD_Miniscope\618\CO2\1_best';
folders{10,3} = 'G:\PMD_Miniscope\620\CO2\1_best';
folders{11,3} = 'G:\PMD_Miniscope\822\CO2\exposure_2_best\H18_M58_S45';

%Hotplate RT
folders{1,4} = 'G:\PMD_Miniscope\611\Hotplate_RT';
folders{2,4} = 'G:\PMD_Miniscope\612\Hotplate_RT';
folders{3,4} = 'G:\PMD_Miniscope\616\Hotplate_RT';
folders{4,4} = 'G:\PMD_Miniscope\617\Hotplate_RT';
folders{5,4} = 'G:\PMD_Miniscope\620\Hotplate_RT';

%Toy Rat Climb
folders{1,5} = 'G:\PMD_Miniscope\20a\RatClimb_Toy\H16_M33_S5';
folders{2,5} = 'G:\PMD_Miniscope\668\RatClimb_Toy';
folders{3,5} = 'G:\PMD_Miniscope\822\RatClimb_Toy\H17_M34_S47';
folders{4,5} = 'G:\PMD_Miniscope\611\RatClimb\2_best';
folders{5,5} = 'G:\PMD_Miniscope\612\RatClimb_Toy';
folders{6,5} = 'G:\PMD_Miniscope\616\RatClimb_Toy';
folders{7,5} = 'G:\PMD_Miniscope\620\RatClimb_Toy';

idx = 1;
titleAll = {'rat','hotplate','CO2','hotplate RT','toy rat'};
figure(1001)

for mouseNum = 1:size(folders,1)
    for assayNum = 1:size(folders,2)
        if isempty(folders{mouseNum,assayNum})
            A = subplot(size(folders,1),size(folders,2),idx)
            box off; A.XTickLabel = []; A.YTickLabel = [];
            if mouseNum == 1
                title(titleAll{assayNum})
            end
            idx = idx + 1;
            continue
        end
        
        cd(folders{mouseNum,assayNum})
        
        load('output_CNMF-E.mat','neuron'); load('good_neurons.mat'); 
        
        if assayNum==1 | assayNum==5
            load('BehaviorMS.mat','climbUFrameMS');
            escapeIndices = climbUFrameMS;
        end
        if assayNum==2 | assayNum==4
            load('BehaviorMS.mat','climbHoseFrameMS');
            escapeIndices = climbHoseFrameMS;
        end
        if assayNum==3
            load('BehaviorMS.mat','jumpFrameMS');
            escapeIndices = jumpFrameMS;
        end
                
        caActivity = neuron.C_raw(find(good_neurons),:);
        if zscore_CA_activity == 1
            for cellNum = 1:size(caActivity,1)
               caActivity(cellNum,:) = zscore(caActivity(cellNum,:)); 
            end
        end
        caActivity = caActivity';
        
        
        %run PCA such that caActivity is now the 'score' output --
        %projections.
        if doPCA==1
            % De-mean
            caActivity = bsxfun(@minus,caActivity,mean(caActivity));
            % Do the PCA
            [coeff,score,latent,tsquared,explained,mu] = pca(caActivity);
            %determine how many PCs to include
            temp = cumsum(explained);
            temp = find(temp > varThresh);
            numPCs = min(temp);
            caActivity = score(:,bottomPC:numPCs);
        end
        
        if hmmAnalysis==1
            numClusters(mouseNum,assayNum) = hmmClusterNum; %SET THE NUMBER OF CLUSTERS FOR HMM HERE!!!
            [Mu, Cov, P, Pi, LL] = hmm(caActivity, length(caActivity), numClusters(mouseNum,assayNum), 500)
            if isnan(Mu(1))
                A = subplot(size(folders,1),size(folders,2),idx)
                box off; A.XTickLabel = []; A.YTickLabel = [];
                if mouseNum == 1
                    title(titleAll{assayNum})
                end
                idx = idx + 1;
                continue
            end
            for numSamples = 1:size(caActivity,1)
                for numClust = 1:size(Mu,1)
                    tempDist(numClust) = pdist2(Mu(numClust,:), caActivity(numSamples,:),'euclidean');
                end
                clusterIdx(numSamples) = find(tempDist==min(tempDist)); clearvars tempDist
            end
        end
                
        if kMeansAnalysis==1
        for clusterNum = 1:maxClusterNum
            clusterIdx(clusterNum,:) = kmeans(caActivity,clusterNum);
        end
        va = evalclusters(caActivity,clusterIdx','CalinskiHarabasz');
        clusterIdx = clusterIdx(va.OptimalK,:);
        numClusters(mouseNum,assayNum) = va.OptimalK;
        end
                
        %save variables for statistical analysis later on:
        clusterIdx_All{mouseNum,assayNum} = clusterIdx;
        escapeInOut_All{mouseNum,assayNum} = escapeIndices;
        temp = zeros(1,length(clusterIdx));
        for escNum=1:size(escapeIndices,1)
           temp(escapeIndices(escNum,1):escapeIndices(escNum,2)) = 1; 
        end
        escapeIndices_All{mouseNum,assayNum} = temp;
        
        if 1==1 %plot the output of clustering and escape in/out
        %plot
        A = subplot(size(folders,1),size(folders,2),idx)
        plot(clusterIdx); hold on;
        
        for i = 1:size(escapeIndices,1)
            jbfill([escapeIndices(i,1) escapeIndices(i,2)], [3 3], [0 0], [1 0 0], [1 0 0], 1, .4);
        end
        xlim([1 length(caActivity)])
        ylim([0 numClusters(mouseNum,assayNum)])
        text(1,0,['#clus=',num2str(numClusters(mouseNum,assayNum))]); hold on;
        box off; A.XTickLabel = []; A.YTickLabel = [];
        idx = idx + 1;
        if mouseNum == 1
           title(titleAll{assayNum})
        end
        end
        
        clearvars clusterIdx
     end
end

%% PERFORM STATS ON THE CLUSTERING OUTPUT

%first, convert all cluster outputs to either '1' or '2'

%flip so '2' corresponds to 'more escape'
for mouseNum = 1:size(folders,1)
    for assayNum = 1:size(folders,2)
        
        if isempty(folders{mouseNum,assayNum})
            continue
        end
            
        for clusterNum = 1:max(max(numClusters)) 
            clustEsc(clusterNum) = sum(clusterIdx_All{mouseNum,assayNum} == clusterNum & escapeIndices_All{mouseNum,assayNum} == 1)
        end
        
        maxClust = find(clustEsc == max(clustEsc));
        maxClustIdx = clusterIdx_All{mouseNum,assayNum} == maxClust(1);
        otherClustIdx = ~maxClustIdx;
        clusterIdx_All{mouseNum,assayNum} = zeros(1,length(clusterIdx_All{mouseNum,assayNum}));
        clusterIdx_All{mouseNum,assayNum}(find(maxClustIdx)) = 2;
        clusterIdx_All{mouseNum,assayNum}(find(otherClustIdx)) = 1;        
    end
end
    
%concatenate sessions across each assay
clust_All = cell(1,size(folders,2));  
for assayNum = 1:size(folders,2)
    for mouseNum = 1:size(folders,1)   
       clust_All{assayNum} = [clust_All{assayNum}, clusterIdx_All{mouseNum,assayNum}];
    end
end
    
esc_All = cell(length(negEscapeLag),size(folders,2));
for lagNum = 1:length(negEscapeLag)
    for assayNum = 1:size(folders,2)
        for mouseNum = 1:size(folders,1)  
            if isempty(escapeIndices_All{mouseNum,assayNum})
                continue
            end
            if negEscapeLag(lagNum)==0
                esc_All{lagNum,assayNum} = [esc_All{lagNum,assayNum}, escapeIndices_All{mouseNum,assayNum}];   
            else
                esc_All{lagNum,assayNum} = [esc_All{lagNum,assayNum}, escapeIndices_All{mouseNum,assayNum}(1+negEscapeLag(lagNum):end), escapeIndices_All{mouseNum,assayNum}(1:negEscapeLag(lagNum))];                   
            end
        end
    end
end

%find actual fraction of escape indices that occurred during cluster '2'
%for each assay
for lagNum = 1:length(negEscapeLag)
    for assayNum = 1:size(folders,2)
        escCount = sum(esc_All{lagNum,assayNum});
        fracEscClust(lagNum,assayNum) = sum(esc_All{lagNum,assayNum}==1 & clust_All{assayNum}==2) ./ escCount;
    end
end

%build bootstrap to compare with 'fracEscClust' values for each assay
rng shuffle;
    for assayNum = 1:size(folders,2)
        for iterNum = 1:iter
            escCount = length(find(esc_All{1,assayNum}==1)); %find total number of escape indices
            idx = randi(length(esc_All{1,assayNum}),1,escCount); %choose random indices, same number as escape indices
            esc_rand = zeros(1,length(esc_All{1,assayNum})); esc_rand(idx) = 1; %build random escape indices for concat sessions

            fracEscClust_boot(iterNum,assayNum) = sum(esc_rand==1 & clust_All{assayNum}==2) ./ escCount;
        end
    end
    
%check it actual value is greater than 95% of bootstrapped values
figure(2002)
for assayNum = 1:size(folders,2)
    subplot(1,5,assayNum)
    hist(fracEscClust_boot(:,assayNum)); hold on;
    plot([fracEscClust(1,assayNum) fracEscClust(1,assayNum)], [0 300],'r')
    title(titleAll{assayNum}); box off;
    xlim([0 1]); ylim([0 300])
    
    pValAssay(assayNum) = 1 - (sum(fracEscClust(1,assayNum) > fracEscClust_boot(:,assayNum)) ./ 1000);
    text(.9, 280, ['p=', num2str(pValAssay(assayNum))],'Color','r')
    ylabel('bootstrap count (random esc times)')
    xlabel('frac. escape during escape cluster')
end


%% plot fraction correct relative to 95th percentile.

figure(3003)
subplot(1,2,1)
labels = {'rat','hotplate','CO2','RT','toy rat'};

%find 95% percentile value for each assay
for assayNum = 1:size(folders,2)
    thresh95(assayNum) = prctile(fracEscClust_boot(:,assayNum),95)
end

bar(fracEscClust(1,:).*100); hold on;
ylim([0 100]); xlim([.5 5.5])
ylabel('% accurately classified escape')
set(gca, 'XTickLabel', labels); box off;

for assayNum = 1:size(folders,2)

    plot([assayNum-.5 assayNum+.5], [thresh95(assayNum).*100 thresh95(assayNum).*100],'r'); hold on;
    
end

% plot the same as above, but just the amount above chance for each assay

subplot(1,2,2)
labels = {'rat','hotplate','CO2','RT','toy rat'};

fracEscClustAboveChance = fracEscClust(1,:) - thresh95;

bar(fracEscClustAboveChance.*100); hold on;
%ylim([0 .9]); 
xlim([.5 5.5]);
ylabel('% accurately classified escape above chance')
set(gca, 'XTickLabel', labels); box off;