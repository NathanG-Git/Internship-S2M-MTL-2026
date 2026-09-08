%%%%%%%%% ------------------------------------- %%%%%%%%%        
%%%%%%%%% Compute joint data from X-sens data   %%%%%%%%%
%%%%%%%%% ------------------------------------- %%%%%%%%%

clear
clc

%% Options

% plot_axes = true; % Place a stop in the code if true!
% save_data = false;
% loadNewXsens = true;

cd('/Users/robin/Documents/MATLAB/Piano_Expression/XSens/Joint_Angles')

load('XSens_MIDI_cut_P01_P08.mat', 'XSens_final');
XSens_part1 = XSens_final;
load('XSens_MIDI_cut_P09_P14.mat', 'XSens_final')
XSens_part2 = XSens_final;

XSens = catstruct(XSens_part1, XSens_part2);
clearvars -except XSens

%% Define variables

Subject = {'P01' 'P02' 'P04' 'P05' 'P06' 'P07' 'P08' 'P09' 'P10' 'P11' 'P12' 'P13' 'P14'};

Condition  = {'Competition', 'Extrait', 'Chopin'};
Expression = {'IE', 'PS'};
Takes_default = {'A', 'B', 'C'};
Takes = {'A', 'B', 'C'};

index = [7 8 9 10 11 12 13 14];
names = {'R scap';'R shoulder';'R elbow';'R wrist';'L scap';'L shoulder';'L elbow';'L wrist'};
joints = {'R_scap';'R_shoulder';'R_elbow';'R_wrist';'L_scap';'L_shoulder';'L_elbow';'L_wrist'};

% Define segments and joints
Pelvis    =  1; L5        =  2; L3       =  3; T12   =  4; T8 = 5; Neck = 6; Head = 7; 
rShould =  8; rUpperArm =  9; rForeArm = 10; rHand = 11; 
lShould = 12; lUpperArm = 13; lForeArm = 14; lHand = 15; 
rUpperLeg = 16; rLowerLeg = 17; rFoot    = 18; rToe  = 19;  
lUpperLeg = 20; lLowerLeg = 21; lFoot    = 22; lToe  = 23; 

segmentsToUse = [Pelvis,           T8,     Head,... 
                 rShould, rUpperArm, rForeArm, rHand,...
                 lShould, lUpperArm, lForeArm, lHand,...
                 rUpperLeg, rLowerLeg, rFoot,...
                 lUpperLeg, lLowerLeg, lFoot];
             
jointsNames   = {'PL' ;  'TH';  'NK';       ...
                 'rST'; 'rGH'; 'rEL'; 'rWR';... 
                 'lST'; 'lGH'; 'lEL'; 'lWR';... 
                 'rCF'; 'rKN'; 'rAK';       ...
                 'lCF'; 'lKN'; 'lAK'};

Joints = [NaN           Pelvis; %PL (parent, child)
          Pelvis            T8; %TH 
          T8              Head; %NK %% See if instead of T8 we should use middle position of rShoulder and lShoulder
          
          T8           rShould; %rST 
          rShould    rUpperArm; %rGH 
          rUpperArm   rForeArm; %rEL 
          rForeArm       rHand; %rWR 
          
          T8           lShould; %lST 
          lShould    lUpperArm; %lGH 
          lUpperArm   lForeArm; %lEL 
          lForeArm       lHand; %lWR 
          
          Pelvis     rUpperLeg; %rCF 
          rUpperLeg  rLowerLeg; %rKN 
          rLowerLeg      rFoot; %rAK 
          
          Pelvis     lUpperLeg; %lCF 
          lUpperLeg  lLowerLeg; %lKN 
          lLowerLeg      lFoot; %lAK 
                             ];

% Define axis sequences 
Sequences = ['yxz'; %PL  
             'yxz'; %T8 
             'yxz'; %NK 
             'yxz'; %rST 
             'yxz'; %rGH 
             'yxz'; %rEL 
             'yxz'; %rWR 
             'yxz'; %lST 
             'yxz'; %lGH 
             'yxz'; %lEL
             'yxz'; %lWR       
             'yxz'; %rCF
             'yxz'; %rKN
             'yxz'; %rAK
             'yxz'; %lCF
             'yxz'; %lKN
             'yxz'; %lAK
                 ];


% Define XSens joints to compile
% L5S1 = 1:3; T9T8 = 10:12; C1Head = 16:18; rScapula = 19:21; rShoulder = 22:24; rElbow = 25:27; rWrist = 28:30; rHip = 43:45; lHip = 55:57; 
% XSensJointsIdx = [L5S1; T9T8; C1Head; rScapula; rShoulder; rElbow; rWrist; rHip; lHip;]; 
XSensJointNames = {'L5S1' 'T9T8' 'C1Head' 'rScapula' 'rShoulder' 'rElbow' 'rWrist' 'rHip' 'lScapula' 'lShoulder' 'lElbow' 'lWrist' 'lHip'};

%Structures to return

JointAngles = struct();
XSens_Angles = struct();
XSens_Segments = struct();

%% Compute and Vizualize the data

for iS = 1:length(Subject)
    fprintf([Subject{iS} '...\n'])

    for iC = 1:length(Condition)
        if ismember(Subject{iS}, {'P05', 'P11', 'P14'}) && strcmp((Condition{iC}), 'Extrait')
                Takes = {'A','B'};
    
            elseif ismember(Subject{iS}, {'P10'}) && strcmp((Condition{iC}), 'Chopin')
                Takes = {'B','C'};
            else
                Takes = Takes_default;
        end


        if strcmp(Condition{iC}, 'Competition')
            % Compile and Organize all Orientation, Position, and Joint Angle data into a table 

            % Extraire la matrice de rotation du pelvis à la frame 1 par rapport à la
            % tpose
            
            xsens_data_ref = XSens.(Subject{iS}).(Condition{iC}).segmentData;
            q_pelvis = quaternion(xsens_data_ref(1).tpose.orientation);
            R_pelvis = q_pelvis.RotationMatrix;
            
            theta_z = atan2(R_pelvis(2,1), R_pelvis(1,1));
            
            % Créer la matrice de rotation inverse pour corriger cette orientation
            R_correction = [cos(-theta_z), -sin(-theta_z), 0;
                            sin(-theta_z),  cos(-theta_z), 0;
                                      0   ,           0   , 1];

            xsens_data = XSens.(Subject{iS}).(Condition{iC}).segmentData;
            xsens_ja   = XSens.(Subject{iS}).(Condition{iC}).jointData;

            RT = compute_joint_angles( ...
                Subject{iS}, Condition{iC}, 'NA', 'NA', ...
                xsens_data, xsens_ja, ...
                Joints, Sequences, ...
                rUpperArm, rForeArm, rHand, ...
                lUpperArm, lForeArm, lHand, ...
                R_correction, JointAngles, XSens_Angles, XSens_Segments);
            % if plot_axes
            %     plot_joint_axes(RT, segmentsToUse, Joints, Subject{iS}, Condition{iC}, 'NA', 'NA');
            % end

        else
            for iE = 1:length(Expression)
                for iT = 1:length(Takes)


                    % Extraire la matrice de rotation du pelvis à la frame 1 par rapport à la
                    % tpose
                    
                    xsens_data_ref = XSens.(Subject{iS}).(Condition{iC}).(Expression{iE}).(Takes{iT}).segmentData;
                    q_pelvis = quaternion(xsens_data_ref(1).tpose.orientation);
                    R_pelvis = q_pelvis.RotationMatrix;
                    
                    theta_z = atan2(R_pelvis(2,1), R_pelvis(1,1));
                    
                    % Créer la matrice de rotation inverse pour corriger cette orientation
                    R_correction = [cos(-theta_z), -sin(-theta_z), 0;
                                    sin(-theta_z),  cos(-theta_z), 0;
                                              0   ,           0   , 1];
                                          
                    xsens_data = XSens.(Subject{iS}).(Condition{iC}).(Expression{iE}).(Takes{iT}).segmentData;
                    xsens_ja   = XSens.(Subject{iS}).(Condition{iC}).(Expression{iE}).(Takes{iT}).jointData;

                    RT = compute_joint_angles( ...
                        Subject{iS}, Condition{iC}, Expression{iE}, Takes{iT}, ...
                        xsens_data, xsens_ja, ...
                        Joints, Sequences, ...
                        rUpperArm, rForeArm, rHand, ...
                        lUpperArm, lForeArm, lHand, ...
                        R_correction, JointAngles, XSens_Angles, XSens_Segments);

                    % if plot_axes
                    %     plot_joint_axes(RT, segmentsToUse, Joints, Subject{iS}, Condition{iC}, Expression{iE}, Takes{iT});
                    % end
                end 
            end
        end
    end 
end



%% Saving

if save_data
    save('JointAngles', 'JointAngles', '-v7.3');
    save('XSens_Angles', 'XSens_Angles', '-v7.3');
    save('XSens_Segments', 'XSens_Segments', '-v7.3');
end

%% Joint Angles Computation

function [RT] = compute_joint_angles( ...
    subject, condition, expression, take, ...
    xsens_data, xsens_ja, ...
    Joints, Sequences, ...
    rUpperArm, rForeArm, rHand, ...
    lUpperArm, lForeArm, lHand, R_correction, ...
    JointAngles, XSens_Angles, XSens_Segments)

    % Position and orientation matrix initialization
    orientation_all = [];
    position_all    = [];
    jointangle_all  = [];
    for i = 1:23
        orientation_all  = [orientation_all, xsens_data(i).orientation];
        position_all     = [position_all,    xsens_data(i).position];
    end
    for i = 1:22
        jointangle_all = [jointangle_all, xsens_ja(i).jointAngle];
    end

    data.orientation_all = array2table(orientation_all);
    data.position_all    = array2table(position_all);
    data.jointangle_all  = array2table(jointangle_all);

    nbSeg    = size(data.orientation_all,2)/4;
    nbFRall  = size(data.orientation_all,1);
    RT = nan(4,4,nbFRall,nbSeg);
    RT(4,1,:)=0; RT(4,2,:)=0; RT(4,3,:)=0; RT(4,4,:)=1;

    for i = 1:nbSeg
        %if i == 1
         %   q_root = array2table(xsens_data(1).tpose.orientation);
         %   q = data.orientation_all(:, i*4-3:i*4);
         %   q(1,1:4) = q_root;
       % else
            q = data.orientation_all(:, i*4-3:i*4);
       % end
        q = quaternion(q{:,:});
        R = squeeze(q.RotationMatrix);
        R = multiprod(R_correction, R);

        if i == 0 % USED TO BE 1, 0 resets the pelvis angles to 0
             T_root = array2table(xsens_data(1).tpose.position);
             T = data.position_all(:, i*3-2:i*3);
             T(1,1:3) = T_root;
        
        else                    
             T = data.position_all(:, i*3-2:i*3);
        end 
        T = T{:,:}';
        T = reshape(T,3,1,[]);

        RT(1:3,1:4,:,i) = [R T];

    end

    % Joint Angles Computation
    nbJ = size(Joints,1);
    JointMat = nan(nbFRall,3,nbJ);
    for j = 1:nbJ
        if isnan(Joints(j,1))
            RJ = multiprod(RT(1:3,1:3,:,Joints(j,2)), multitransp(RT(1:3,1:3,1,Joints(j,2))));
        else
            RJ = multiprod(RT(1:3,1:3,:,Joints(j,2)), multitransp(RT(1:3,1:3,:,Joints(j,1))));
        end
        JointMat(:,:,j) = squeeze(angleRotation(RJ, Sequences(j,:)))';
        if ismember(Joints(j,2), [rUpperArm, lUpperArm, rForeArm, lForeArm])
            JointMat(:,:,j) = unwrap(JointMat(:,:,j));
        end
        JointMat(:,:,j) = rad2deg(JointMat(:,:,j));
    end

    % Stockage => PAS DE FRAMES (à demander)
    if strcmp(expression, 'NA') && strcmp(take, 'NA')
        JointAngles.(subject).(condition).JA = JointMat;
        XSens_Angles.(subject).(condition).JA = xsens_ja;
        XSens_Segments.(subject).(condition).SEG = xsens_data;
    else
        JointAngles.(subject).(condition).(expression).(take).JA = JointMat;
        XSens_Angles.(subject).(condition).(expression).(take).JA = xsens_ja;
        XSens_Segments.(subject).(condition).(expression).(take).SEG = xsens_data;
    end

    assignin('caller', 'JointAngles', JointAngles);
    assignin('caller', 'XSens_Angles', XSens_Angles);
    assignin('caller', 'XSens_Segments', XSens_Segments);
end


%% Plot function
function plot_joint_axes(RT, segmentsToUse, Joints, subject, condition, expression, take)

    frameTP = 1;

    Fig1 = figure('Name', [subject ' - ' condition ' - ' expression ' - ' take]);
    
    % Plot local axes
    for iSeg = 1:length(segmentsToUse)
        idx = segmentsToUse(iSeg);
        plotAxes(RT(:,:,frameTP,idx));
        hold on
    end

    % Plot skeleton links
    seg = Joints(2:end,:);
%                     seg = [1 5; 5 6; 6 7; 8 9; 9 10; 10 11; 12 13; 13 14; 14 15];
                %seg = [1 5; 5 6; 6 7;       9 10; 10 11;        13 14; 14 15];
    for i = 1:size(seg,1)
        hold on
        plot3(squeeze(RT(1,4,frameTP, seg(i,:))), ...
              squeeze(RT(2,4,frameTP, seg(i,:))), ...
              squeeze(RT(3,4,frameTP, seg(i,:))), ...
              'k','LineWidth', 3);
        plot3(squeeze(RT(1,4,frameTP, seg(i,1))), ...
              squeeze(RT(2,4,frameTP, seg(i,1))), ...
              squeeze(RT(3,4,frameTP, seg(i,1))), ...
              'ko','LineWidth', 3);
    end

    axis equal
    title([subject ' - ' condition ' - ' expression ' - ' take])
end


