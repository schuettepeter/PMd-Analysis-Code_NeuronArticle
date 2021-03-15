%% The following code was used to determine which putative neurons significantly encoded escape behavior.

clearvars -except folders

%behaviors to classify:
climbUBehav = 0; %Rat climb
climbHoseBehav = 0; %Heated plate
jumpBehav = 0; %CO2

includeBehav = 0;
includePreBehav = 1;
includeOnlyPreBehav = 0;
includeOnlyPostBehav = 0;
numSamplesPre = round(7.5 .* 5); %fps * # seconds

samples = 1000; %bootstrap distribution

load('output_CNMF-E.mat', 'neuron') %load output from CNMF-E
load('Tracking.mat'); 
load('good_neurons.mat')

load('BehaviorMS.mat','jumpFrameMS','jumpIndicesMS','climbHoseFrameMS','climbHoseIndicesMS','climbUFrameMS','climbUIndicesMS') 
    
good_neurons = find(good_neurons);

%%
if jumpBehav == 1 && exist('jumpFrameMS')

            jumpIndicesMS = zeros(length(jumpIndicesMS),1);  
            
            %include only behavior samples
            if includeBehav == 1
            for behavNum = 1:size(jumpFrameMS,1)
                jumpIndicesMS([jumpFrameMS(behavNum,1):jumpFrameMS(behavNum,2)]) = 1;
            end
            end
            
            %include prior x samples in model
            if includePreBehav == 1
                jumpFrameMS(:,1) = jumpFrameMS(:,1) - numSamplesPre;
                jumpFrameMS(find(jumpFrameMS<1)) = 1;
                jumpFrameMS(find(jumpFrameMS>length(jumpIndicesMS))) = length(jumpIndicesMS);
            for behavNum = 1:size(jumpFrameMS,1)
                jumpIndicesMS([jumpFrameMS(behavNum,1):jumpFrameMS(behavNum,2)]) = 1;
            end
            end
            
            %include only prior x samples in model
            if includeOnlyPreBehav == 1
                jumpFrameMS(:,2) = jumpFrameMS(:,1);
                jumpFrameMS(:,1) = jumpFrameMS(:,1) - numSamplesPre;
                jumpFrameMS(find(jumpFrameMS<1)) = 1;
                jumpFrameMS(find(jumpFrameMS>length(jumpIndicesMS))) = length(jumpIndicesMS);
            for behavNum = 1:size(jumpFrameMS,1)
                jumpIndicesMS([jumpFrameMS(behavNum,1):jumpFrameMS(behavNum,2)]) = 1;
            end
            end

            %include only post x samples in model
            if includeOnlyPostBehav == 1
                jumpFrameMS(:,2) = jumpFrameMS(:,1);
                jumpFrameMS(:,2) = jumpFrameMS(:,2) + numSamplesPre;
                jumpFrameMS(find(jumpFrameMS<1)) = 1;
                jumpFrameMS(find(jumpFrameMS>length(jumpIndicesMS))) = length(jumpIndicesMS);
            for behavNum = 1:size(jumpFrameMS,1)
                jumpIndicesMS([jumpFrameMS(behavNum,1):jumpFrameMS(behavNum,2)]) = 1;
            end
            end          
            
            Indices = jumpIndicesMS;
                
        for seg = 1:length(good_neurons)

            CalciumTrace = neuron.C_raw(good_neurons(seg),:)';

            while length(Indices) > length(CalciumTrace)
                Indices = Indices(1:end-1);
            end
            while length(Indices) < length(CalciumTrace)
                Indices = [Indices; 0];
            end
            
            tbl = table(Indices, CalciumTrace);
            lm = fitlm(tbl, 'CalciumTrace~Indices');

            coeffSeg(good_neurons(seg),:) = table2array(lm.Coefficients(2:end,[1,4]));

        end
    
%create bootstrap distributions for each neuron to compare with output from
%previous section -- 
% move freeze behaviors as epochs -- not individual samples

    disp(['building bootstrap distribution'])
    
    InOut = jumpFrameMS;
    BoutLength =  InOut(:,2) - InOut(:,1);
    
    for iter = 1:samples
        %n = length(FreezeIndices); %number of random integers to choose from

        l = size(BoutLength,1); %how many numbers to pick
        
        offsetFromEnd = (max(BoutLength)) + 7;
        
        startIndices = datasample([7:length(neuron.C)-offsetFromEnd],l);
        lengthOrder = BoutLength(randperm(size(InOut,1),size(InOut,1))); 
        
        Indices = zeros(length(neuron.C),1);
        
        for i = 1:length(startIndices)
           Indices(startIndices(i):startIndices(i)+lengthOrder(i)) = 1; 
        end
        
            for seg = 1:length(good_neurons)

                CalciumTrace = neuron.C_raw(good_neurons(seg),:)';

                tbl = table(Indices, CalciumTrace);
                lm = fitlm(tbl, 'CalciumTrace~Indices');

                coeffSeg_shuff(good_neurons(seg),iter) = table2array(lm.Coefficients(2:end,1));

            end
    end
    
%determine, for each seg, where it's coefficient falls in the distribution
%-- EITHER ABOVE OR BELOW (so either positively or negatively modulated by behavior)
    for seg = 1:length(good_neurons)
        %positively modulated by behavior
        temp = length(find(coeffSeg_shuff(good_neurons(seg),:) > coeffSeg(good_neurons(seg),1))) ./ samples;
        if temp <= .05 & coeffSeg(good_neurons(seg),2) <= .05 %p-value of bootstrap and coefficient.
           jumpSeg(good_neurons(seg)) = 1;
        else
           jumpSeg(good_neurons(seg)) = 0;
        end
        
        %negatively modulated by behavior
        temp = length(find(coeffSeg_shuff(good_neurons(seg),:) < coeffSeg(good_neurons(seg),1))) ./ samples; %fraction bootstrapped vals that are LESS THAN ACTUAL VALUE
        if temp <= .05 & coeffSeg(good_neurons(seg),2) <= .05 %p-value of bootstrap and coefficient.
           jumpSeg(good_neurons(seg)) = -1;
        %else
        %   jumpSeg(good_neurons(seg)) = 0;
        end
    end
    
    if includeBehav == 1
        save('jumpSeg.mat','jumpSeg','coeffSeg','coeffSeg_shuff')
    end
    if includePreBehav == 1
        jumpSegPreBehav = jumpSeg; clearvars jumpSeg
        save('jumpSegPreBehav.mat','jumpSegPreBehav','coeffSeg','coeffSeg_shuff')
    end
    if includeOnlyPreBehav == 1
        jumpSegPreOnlyBehav = jumpSeg; clearvars jumpSeg
        save('jumpSegPreOnlyBehav.mat','jumpSegPreOnlyBehav','coeffSeg','coeffSeg_shuff')
    end
    if includeOnlyPostBehav == 1
        jumpSegPostOnlyBehav = jumpSeg; clearvars jumpSeg
        save('jumpSegPostOnlyBehav.mat','jumpSegPostOnlyBehav','coeffSeg','coeffSeg_shuff')
    end
    
end

clearvars coeffSeg coeffSeg_shuff
%%
if climbHoseBehav == 1 && exist('climbHoseFrameMS')
    
            climbHoseIndicesMS = zeros(length(climbHoseIndicesMS),1);  
            
            %include only behavior samples
            if includeBehav == 1
            for behavNum = 1:size(climbHoseFrameMS,1)
                climbHoseIndicesMS([climbHoseFrameMS(behavNum,1):climbHoseFrameMS(behavNum,2)]) = 1;
            end
            end
            
            %include prior x samples in model
            if includePreBehav == 1
                climbHoseFrameMS(:,1) = climbHoseFrameMS(:,1) - numSamplesPre;
                climbHoseFrameMS(find(climbHoseFrameMS<1)) = 1;
                climbHoseFrameMS(find(climbHoseFrameMS>length(climbHoseIndicesMS))) = length(climbHoseIndicesMS);
            for behavNum = 1:size(climbHoseFrameMS,1)
                climbHoseIndicesMS([climbHoseFrameMS(behavNum,1):climbHoseFrameMS(behavNum,2)]) = 1;
            end
            end
            
            %include only prior x samples in model
            if includeOnlyPreBehav == 1
                climbHoseFrameMS(:,2) = climbHoseFrameMS(:,1);
                climbHoseFrameMS(:,1) = climbHoseFrameMS(:,1) - numSamplesPre;
                climbHoseFrameMS(find(climbHoseFrameMS<1)) = 1;
                climbHoseFrameMS(find(climbHoseFrameMS>length(climbHoseIndicesMS))) = length(climbHoseIndicesMS);
            for behavNum = 1:size(climbHoseFrameMS,1)
                climbHoseIndicesMS([climbHoseFrameMS(behavNum,1):climbHoseFrameMS(behavNum,2)]) = 1;
            end
            end

            %include only post x samples in model
            if includeOnlyPostBehav == 1
                climbHoseFrameMS(:,2) = climbHoseFrameMS(:,1);
                climbHoseFrameMS(:,2) = climbHoseFrameMS(:,2) + numSamplesPre;
                climbHoseFrameMS(find(climbHoseFrameMS<1)) = 1;
                climbHoseFrameMS(find(climbHoseFrameMS>length(climbHoseIndicesMS))) = length(climbHoseIndicesMS);
            for behavNum = 1:size(climbHoseFrameMS,1)
                climbHoseIndicesMS([climbHoseFrameMS(behavNum,1):climbHoseFrameMS(behavNum,2)]) = 1;
            end
            end
            
            
            Indices = climbHoseIndicesMS;
                
        for seg = 1:length(good_neurons)

            CalciumTrace = neuron.C_raw(good_neurons(seg),:)';

            while length(Indices) > length(CalciumTrace)
                Indices = Indices(1:end-1);
            end
            while length(Indices) < length(CalciumTrace)
                Indices = [Indices; 0];
            end
            
            tbl = table(Indices, CalciumTrace);
            lm = fitlm(tbl, 'CalciumTrace~Indices');

            coeffSeg(good_neurons(seg),:) = table2array(lm.Coefficients(2:end,[1,4]));

        end

    
%create bootstrap distributions for each neuron to compare with output from
%previous section -- 
% move freeze behaviors as epochs -- not individual samples

    disp(['building bootstrap distribution'])
    %Indices = climbHoseIndicesMS;
    InOut = climbHoseFrameMS;
    BoutLength =  InOut(:,2) - InOut(:,1);
    
    for iter = 1:samples
        %n = length(FreezeIndices); %number of random integers to choose from

        l = size(BoutLength,1); %how many numbers to pick
        
        offsetFromEnd = (max(BoutLength)) + 7;
        
        startIndices = datasample([7:length(neuron.C)-offsetFromEnd],l);
        lengthOrder = BoutLength(randperm(size(InOut,1),size(InOut,1))); 
        
        Indices = zeros(length(neuron.C),1);
        
        for i = 1:length(startIndices)
           Indices(startIndices(i):startIndices(i)+lengthOrder(i)) = 1; 
        end
        
            for seg = 1:length(good_neurons)

                CalciumTrace = neuron.C_raw(good_neurons(seg),:)';

                tbl = table(Indices, CalciumTrace);
                lm = fitlm(tbl, 'CalciumTrace~Indices');

                coeffSeg_shuff(good_neurons(seg),iter) = table2array(lm.Coefficients(2:end,1));

            end
    end
    
%determine, for each seg, where it's coefficient falls in the distribution
%-- EITHER ABOVE OR BELOW (so either positively or negatively modulated by behavior)
    for seg = 1:length(good_neurons)
        %positively modulated by behavior
        temp = length(find(coeffSeg_shuff(good_neurons(seg),:) > coeffSeg(good_neurons(seg),1))) ./ samples;
        if temp <= .05 & coeffSeg(good_neurons(seg),2) <= .05 %p-value of bootstrap and coefficient.
           climbHoseSeg(good_neurons(seg)) = 1;
        else
           climbHoseSeg(good_neurons(seg)) = 0;
        end
        
        %negatively modulated by behavior
        temp = length(find(coeffSeg_shuff(good_neurons(seg),:) < coeffSeg(good_neurons(seg),1))) ./ samples; %fraction bootstrapped vals that are LESS THAN ACTUAL VALUE
        if temp <= .05 & coeffSeg(good_neurons(seg),2) <= .05 %p-value of bootstrap and coefficient.
           climbHoseSeg(good_neurons(seg)) = -1;
        end
    end
    
    if includeBehav == 1
        save('climbHoseSeg.mat','climbHoseSeg','coeffSeg','coeffSeg_shuff')
    end
    if includePreBehav == 1
        climbHoseSegPreBehav = climbHoseSeg; clearvars climbHoseSeg
        save('climbHoseSegPreBehav.mat','climbHoseSegPreBehav','coeffSeg','coeffSeg_shuff')
    end
    if includeOnlyPreBehav == 1
        climbHoseSegPreOnlyBehav = climbHoseSeg; clearvars climbHoseSeg
        save('climbHoseSegPreOnlyBehav.mat','climbHoseSegPreOnlyBehav','coeffSeg','coeffSeg_shuff')
    end
    if includeOnlyPostBehav == 1
        climbHoseSegPostOnlyBehav = climbHoseSeg; clearvars climbHoseSeg
        save('climbHoseSegPostOnlyBehav.mat','climbHoseSegPostOnlyBehav','coeffSeg','coeffSeg_shuff')
    end

end
%%
clearvars coeffSeg coeffSeg_shuff

if climbUBehav == 1 && exist('climbUFrameMS')
    
            climbUIndicesMS = zeros(length(climbUIndicesMS),1);    
            
            %include only behavior samples
            if includeBehav == 1
            for behavNum = 1:size(climbUFrameMS,1)
                climbUIndicesMS([climbUFrameMS(behavNum,1):climbUFrameMS(behavNum,2)]) = 1;
            end
            end
            
            %include prior x samples in model
            if includePreBehav == 1
                climbUFrameMS(:,1) = climbUFrameMS(:,1) - numSamplesPre;
                climbUFrameMS(find(climbUFrameMS<1)) = 1;
                climbUFrameMS(find(climbUFrameMS>length(climbUIndicesMS))) = length(climbUIndicesMS);
            for behavNum = 1:size(climbUFrameMS,1)
                climbUIndicesMS([climbUFrameMS(behavNum,1):climbUFrameMS(behavNum,2)]) = 1;
            end
            end
            
            %include only prior x samples in model
            if includeOnlyPreBehav == 1
                climbUFrameMS(:,2) = climbUFrameMS(:,1);
                climbUFrameMS(:,1) = climbUFrameMS(:,1) - numSamplesPre;
                climbUFrameMS(find(climbUFrameMS<1)) = 1;
                climbUFrameMS(find(climbUFrameMS>length(climbUIndicesMS))) = length(climbUIndicesMS);
            for behavNum = 1:size(climbUFrameMS,1)
                climbUIndicesMS([climbUFrameMS(behavNum,1):climbUFrameMS(behavNum,2)]) = 1;
            end
            end

            %include only post x samples in model
            if includeOnlyPostBehav == 1
                climbUFrameMS(:,2) = climbUFrameMS(:,1);
                climbUFrameMS(:,2) = climbUFrameMS(:,2) + numSamplesPre;
                climbUFrameMS(find(climbUFrameMS<1)) = 1;
                climbUFrameMS(find(climbUFrameMS>length(climbUIndicesMS))) = length(climbUIndicesMS);
            for behavNum = 1:size(climbUFrameMS,1)
                climbUIndicesMS([climbUFrameMS(behavNum,1):climbUFrameMS(behavNum,2)]) = 1;
            end
            end            
            
            Indices = climbUIndicesMS;
                
        for seg = 1:length(good_neurons)

            CalciumTrace = neuron.C_raw(good_neurons(seg),:)';

            while length(Indices) > length(CalciumTrace)
                Indices = Indices(1:end-1);
            end
            while length(Indices) < length(CalciumTrace)
                Indices = [Indices; 0];
            end
            
            tbl = table(Indices, CalciumTrace);
            lm = fitlm(tbl, 'CalciumTrace~Indices');

            coeffSeg(good_neurons(seg),:) = table2array(lm.Coefficients(2:end,[1,4]));

        end


    %create bootstrap distributions for each neuron to compare with output from
    %previous section -- 
    % move freeze behaviors as epochs -- not individual samples

    disp(['building bootstrap distribution'])
    %Indices = climbUIndicesMS;
    InOut = climbUFrameMS;
    BoutLength =  InOut(:,2) - InOut(:,1);
    
    for iter = 1:samples
        %n = length(FreezeIndices); %number of random integers to choose from

        l = size(BoutLength,1); %how many numbers to pick
        
        offsetFromEnd = (max(BoutLength)) + 7;
        
        startIndices = datasample([7:length(neuron.C)-offsetFromEnd],l);
        lengthOrder = BoutLength(randperm(size(InOut,1),size(InOut,1))); 
        
        Indices = zeros(length(neuron.C),1);
        
        for i = 1:length(startIndices)
           Indices(startIndices(i):startIndices(i)+lengthOrder(i)) = 1; 
        end
        
            for seg = 1:length(good_neurons)

                CalciumTrace = neuron.C_raw(good_neurons(seg),:)';

                tbl = table(Indices, CalciumTrace);
                lm = fitlm(tbl, 'CalciumTrace~Indices');

                coeffSeg_shuff(good_neurons(seg),iter) = table2array(lm.Coefficients(2:end,1));

            end
    end
    
%determine, for each seg, where it's coefficient falls in the distribution
%-- EITHER ABOVE OR BELOW (so either positively or negatively modulated by behavior)
    for seg = 1:length(good_neurons)
        %positively modulated by behavior
        temp = length(find(coeffSeg_shuff(good_neurons(seg),:) > coeffSeg(good_neurons(seg),1))) ./ samples;
        if temp <= .05 & coeffSeg(good_neurons(seg),2) <= .05 %p-value of bootstrap and coefficient.
           climbUSeg(good_neurons(seg)) = 1;
        else
           climbUSeg(good_neurons(seg)) = 0;
        end
        
        %negatively modulated by behavior
        temp = length(find(coeffSeg_shuff(good_neurons(seg),:) < coeffSeg(good_neurons(seg),1))) ./ samples; %fraction bootstrapped vals that are LESS THAN ACTUAL VALUE
        if temp <= .05 & coeffSeg(good_neurons(seg),2) <= .05 %p-value of bootstrap and coefficient.
           climbUSeg(good_neurons(seg)) = -1;
        end
    end
    
    if includeBehav == 1
        save('climbUSeg.mat','climbUSeg','coeffSeg','coeffSeg_shuff')
    end
    if includePreBehav == 1
        climbUSegPreBehav = climbUSeg; clearvars climbUSeg
        save('climbUSegPreBehav.mat','climbUSegPreBehav','coeffSeg','coeffSeg_shuff')
    end
    if includeOnlyPreBehav == 1
        climbUSegPreOnlyBehav = climbUSeg; clearvars climbUSeg
        save('climbUSegPreOnlyBehav.mat','climbUSegPreOnlyBehav','coeffSeg','coeffSeg_shuff')
    end
    if includeOnlyPostBehav == 1
        climbUSegPostOnlyBehav = climbUSeg; clearvars climbUSeg
        save('climbUSegPostOnlyBehav.mat','climbUSegPostOnlyBehav','coeffSeg','coeffSeg_shuff')
    end

end

clearvars coeffSeg coeffSeg_shuff
