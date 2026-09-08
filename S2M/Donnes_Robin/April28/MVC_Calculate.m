
close all
clear;clc

%% Load

load("C:\Users\jerem\Documents\Projet de stage\DATA_Jeremy\April28\c3dFile_MVC.mat")

%% Options

graph_mvcs             = false;

%% Parameters

addpath('C:\Users\turne\Documents\MATLAB\ezc3d_matlab')
rootFolder = 'D:\turne\uleth\PhD\Thesis\Data_Collection\DATA';
extension    = '.c3d';
SubjectsToDo = {'P01' 'P02' 'P03' 'P04' 'P05' 'P06' 'P07' 'P08' 'P09' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21' 'P22' 'P23' 'P24'};

MVC_Trials = {'Bic_1';'Bic_2';
              'Brachio_1';'Brachio_2';...
              'Delt_1';'Delt_2';...
              'Trap_1';'Trap_2';...
              'Tri_1';'Tri_2';...
              'FlxRad_Finger_1';'FlxRad_Finger_2';...
              'FlxUln_Wrist_1';'FlxUln_Wrist_2'};

Keyboard   = {'Norm','CS60'};
Excerpts   = {'Rch'; 'Chp'; 'Deb'};
Cond       = {'L'; 'S'};
GoodTrials = {'A'; 'B'};

%%
clc
close all

for iS = 1:length(SubjectsToDo)
    if iS == 14 % No Brachialis MVC (participant complications)
        MVC_Trials = {'Bic_1';'Bic_2';
                      'Delt_1';'Delt_2';...
                      'Trap_1';'Trap_2';...
                      'Tri_1';'Tri_2';...
                      'FlxRad_Finger_1';'FlxRad_Finger_2';...
                      'FlxUln_Wrist_1';'FlxUln_Wrist_2'};
        Muscles = {'Flexor digitorum';'Wrist flexor';'Bicep';'Tricep';'Deltoid';'Trap'};
    else
        MVC_Trials = {'Bic_1';'Bic_2';
                      'Brachio_1';'Brachio_2';...
                      'Delt_1';'Delt_2';...
                      'Trap_1';'Trap_2';...
                      'Tri_1';'Tri_2';...
                      'FlxRad_Finger_1';'FlxRad_Finger_2';...
                      'FlxUln_Wrist_1';'FlxUln_Wrist_2'};
        Muscles = {'Brachioradialis';'Flexor digitorum';'Wrist flexor';'Bicep';'Tricep';'Deltoid';'Trap'};        
    end

    for iM = 1:length(MVC_Trials)
        dataset = c3dFile_MVC.(SubjectsToDo{iS}).(MVC_Trials{iM});        

        % Graphs
        if graph_mvcs == true
            figure(iM)
            sgtitle(MVC_Trials{iM})
            for i = 1:length(Muscles)
                subplot(3,3,i)
                plot(dataset(:,i))
                subtitle(Muscles{i})
            end
        end

        if iS == 14
            % Bicep MVCs
            if iM == 1 || iM == 2
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,3);
    
            % Deltoid MVCs
            elseif iM == 3 || iM == 4
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM})= dataset(:,5);
    
            % Trap MVCs    
            elseif iM == 5 || iM == 6
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,6);
    
            % Tricep MVCs    
            elseif iM == 7 || iM == 8
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,4);
    
            % Flexor digitorum MVCs    
            elseif iM == 9 || iM == 10
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,2);
    
            % Wrist flexor MVCs    
            elseif iM == 11 || iM == 12
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,1); 
            end

        else

            % Bicep MVCs
            if iM == 1 || iM == 2
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,4);
    
            % Brachio MVCs    
            elseif iM == 3 || iM == 4
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,1);
    
            % Deltoid MVCs
            elseif iM == 5 || iM == 6
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM})= dataset(:,6);
    
            % Trap MVCs    
            elseif iM == 7 || iM == 8
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,7);
    
            % Tricep MVCs    
            elseif iM == 9 || iM == 10
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,5);
    
            % Flexor digitorum MVCs    
            elseif iM == 11 || iM == 12
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,2);
    
            % Wrist flexor MVCs    
            elseif iM == 13 || iM == 14
                MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM}) = dataset(:,3); 
            end
        end
    end
end

%% Filtering Parameters

Freq = 2000;
Ts = 1/Freq;

[b,a] = butter(2,[10 400]/(Freq/2),"bandpass"); % Parametre du filtre BP (Band Pass) 10-400 Hz
wind_length = round(30*10^(-3)/Ts); % 30 ms == 60 frames

%% MVC - Filtering and computing of RMS envelope

for iS = 1:length(SubjectsToDo)
    if iS == 14 % No Brachialis MVC (participant complications)
        MVC_Trials = {'Bic_1';'Bic_2';
                      'Delt_1';'Delt_2';...
                      'Trap_1';'Trap_2';...
                      'Tri_1';'Tri_2';...
                      'FlxRad_Finger_1';'FlxRad_Finger_2';...
                      'FlxUln_Wrist_1';'FlxUln_Wrist_2'};
        Muscles = {'Flexor digitorum';'Wrist flexor';'Bicep';'Tricep';'Deltoid';'Trap'};
    else
        MVC_Trials = {'Bic_1';'Bic_2';
                      'Brachio_1';'Brachio_2';...
                      'Delt_1';'Delt_2';...
                      'Trap_1';'Trap_2';...
                      'Tri_1';'Tri_2';...
                      'FlxRad_Finger_1';'FlxRad_Finger_2';...
                      'FlxUln_Wrist_1';'FlxUln_Wrist_2'};
        Muscles = {'Brachioradialis';'Flexor digitorum';'Wrist flexor';'Bicep';'Tricep';'Deltoid';'Trap'};        
    end
    for iM = 1:length(MVC_Trials)
        MVCTrial = MVC_data.(SubjectsToDo{iS}).(MVC_Trials{iM});
        tp      = MVCTrial - mean(MVCTrial);
        tp_filt = filtfilt(b,a,tp);
        
        allMVCData.(SubjectsToDo{iS}).(MVC_Trials{iM}).RawData               = MVCTrial;
        allMVCData.(SubjectsToDo{iS}).(MVC_Trials{iM}).Centered              = tp;
        allMVCData.(SubjectsToDo{iS}).(MVC_Trials{iM}).DataFiltered          = tp_filt;
        [allMVCData.(SubjectsToDo{iS}).(MVC_Trials{iM}).DataEnvelope, lower] = envelope(tp_filt, wind_length, 'rms');
        
        %MVC_Length(iM,iS) = length(allMVCData.(SubjectsToDo{iS}).(Conditions{iC}).DataEnvelope);
    end
end

for iS = 1:length(SubjectsToDo)
    if iS == 14
        MVC_cat.(SubjectsToDo{iS}).Bicep = vertcat(allMVCData.(SubjectsToDo{iS}).Bic_1.DataEnvelope,             allMVCData.(SubjectsToDo{iS}).Bic_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).Delt = vertcat(allMVCData.(SubjectsToDo{iS}).Delt_1.DataEnvelope,             allMVCData.(SubjectsToDo{iS}).Delt_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).Trap = vertcat(allMVCData.(SubjectsToDo{iS}).Trap_1.DataEnvelope,             allMVCData.(SubjectsToDo{iS}).Trap_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).Tricep = vertcat(allMVCData.(SubjectsToDo{iS}).Tri_1.DataEnvelope,            allMVCData.(SubjectsToDo{iS}).Tri_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).FFlexor = vertcat(allMVCData.(SubjectsToDo{iS}).FlxRad_Finger_1.DataEnvelope, allMVCData.(SubjectsToDo{iS}).FlxRad_Finger_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).WFlexor = vertcat(allMVCData.(SubjectsToDo{iS}).FlxUln_Wrist_1.DataEnvelope,  allMVCData.(SubjectsToDo{iS}).FlxUln_Wrist_2.DataEnvelope);
    else
        MVC_cat.(SubjectsToDo{iS}).Bicep = vertcat(allMVCData.(SubjectsToDo{iS}).Bic_1.DataEnvelope,             allMVCData.(SubjectsToDo{iS}).Bic_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).Brachio = vertcat(allMVCData.(SubjectsToDo{iS}).Brachio_1.DataEnvelope,       allMVCData.(SubjectsToDo{iS}).Brachio_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).Delt = vertcat(allMVCData.(SubjectsToDo{iS}).Delt_1.DataEnvelope,             allMVCData.(SubjectsToDo{iS}).Delt_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).Trap = vertcat(allMVCData.(SubjectsToDo{iS}).Trap_1.DataEnvelope,             allMVCData.(SubjectsToDo{iS}).Trap_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).Tricep = vertcat(allMVCData.(SubjectsToDo{iS}).Tri_1.DataEnvelope,            allMVCData.(SubjectsToDo{iS}).Tri_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).FFlexor = vertcat(allMVCData.(SubjectsToDo{iS}).FlxRad_Finger_1.DataEnvelope, allMVCData.(SubjectsToDo{iS}).FlxRad_Finger_2.DataEnvelope);
        MVC_cat.(SubjectsToDo{iS}).WFlexor = vertcat(allMVCData.(SubjectsToDo{iS}).FlxUln_Wrist_1.DataEnvelope,  allMVCData.(SubjectsToDo{iS}).FlxUln_Wrist_2.DataEnvelope);
    end
end

% Sort from highest to smallest value
for iS = 1:length(SubjectsToDo)
    if iS == 14 % No Brachialis MVC (participant complications)
        MVC_Trials = {'Bic_1';'Bic_2';
                      'Delt_1';'Delt_2';...
                      'Trap_1';'Trap_2';...
                      'Tri_1';'Tri_2';...
                      'FlxRad_Finger_1';'FlxRad_Finger_2';...
                      'FlxUln_Wrist_1';'FlxUln_Wrist_2'};
        MVC_muscles = {'Bicep';'Delt';'Trap';'Tricep';'FFlexor';'WFlexor'};
    else
        MVC_Trials = {'Bic_1';'Bic_2';
                      'Brachio_1';'Brachio_2';...
                      'Delt_1';'Delt_2';...
                      'Trap_1';'Trap_2';...
                      'Tri_1';'Tri_2';...
                      'FlxRad_Finger_1';'FlxRad_Finger_2';...
                      'FlxUln_Wrist_1';'FlxUln_Wrist_2'};
        MVC_muscles = {'Bicep';'Brachio';'Delt';'Trap';'Tricep';'FFlexor';'WFlexor'};
    end
    for iM = 1:length(MVC_muscles)
        Sorted.(SubjectsToDo{iS}).(MVC_muscles{iM}) = sort(MVC_cat.(SubjectsToDo{iS}).(MVC_muscles{iM}), 'descend'); 
        MVCsorted = Sorted.(SubjectsToDo{iS}).(MVC_muscles{iM}); 
        MVC_value_Med.(SubjectsToDo{iS}).(MVC_muscles{iM}) = median(MVCsorted(1:2000,:)); % Caclulate median MVC value
    end
end
save('C:\Users\jerem\Documents\Projet de stage\MVC_value_Med.mat', 'MVC_value_Med');
disp('MVC_value_Med.mat sauvegardé !');