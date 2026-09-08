%%%%%%%%% ------------------------------------- %%%%%%%%%        
%%%%%%%%% Compute joint data from X-sens data   %%%%%%%%%
%%%%%%%%% ------------------------------------- %%%%%%%%%

clear
clc
%clearvars -except XSens
close all

%% Options

plot_axes = false; % If true, place a stop in the code after plotting or else MATLAB will crash!!!!!
save_data = false;

%% Load Data

% Change the files names accordingly
cd('D:\turne\uleth\PhD\Thesis\Data_Collection\DATA\MATLAB_Data')
load('KeyTimes_sec.mat')
% load('XSens_Norm.mat')
% 
% XSens.P01 = XSens_Norm.P01;
% XSens.P06 = XSens_Norm.P06;
% XSens.P09 = XSens_Norm.P09;
% XSens.P11 = XSens_Norm.P11;
% XSens.P13 = XSens_Norm.P13;
% XSens.P14 = XSens_Norm.P14;
% XSens.P16 = XSens_Norm.P16;
% XSens.P19 = XSens_Norm.P19;
% XSens.P22 = XSens_Norm.P22;
% XSens.P23 = XSens_Norm.P23;
% clear XSens_Norm
% 
% load('XSens_Both_1.mat')
% XSens.P02 = XSens_Both_1.P02;
% XSens.P03 = XSens_Both_1.P03;
% XSens.P04 = XSens_Both_1.P04;
% XSens.P05 = XSens_Both_1.P05;
% XSens.P07 = XSens_Both_1.P07;
% XSens.P08 = XSens_Both_1.P08;
% XSens.P10 = XSens_Both_1.P10;
% clear XSens_Both_1
% 
% load('XSens_Both_2.mat')
% XSens.P12 = XSens_Both_2.P12;
% XSens.P15 = XSens_Both_2.P15;
% XSens.P17 = XSens_Both_2.P17;
% XSens.P18 = XSens_Both_2.P18;
% XSens.P20 = XSens_Both_2.P20;
% XSens.P21 = XSens_Both_2.P21;
% XSens.P24 = XSens_Both_2.P24;
% clear XSens_Both_2

%% Define variables

SavePath = 'D:\turne\uleth\PhD\Thesis\Data_Collection\DATA\MATLAB_Data';

SubjectsToDo = {'P01' 'P02' 'P03' 'P04' 'P05' 'P06' 'P07' 'P08' 'P09' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21' 'P22' 'P23' 'P24'};

Keyboard   = {'Norm','CS60'};
Excerpts   = {'Rch'; 'Chp'; 'Deb'};
Cond       = {'L'; 'S'};
GoodTrials = {'A'; 'B'};

index = [7 8 9 10 11 12 13 14];
names = {'R scap';'R shoulder';'R elbow';'R wrist';'L scap';'L shoulder';'L elbow';'L wrist'};
joints = {'R_scap';'R_shoulder';'R_elbow';'R_wrist';'L_scap';'L_shoulder';'L_elbow';'L_wrist'};

% add functions for processing of kinematic data 
addpath('D:\turne\uleth\PhD\Thesis\Expression du piano\MATLAB\Functions\XSens_Angles')

fs_xsens = 60;

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
          T8              Head; %NK
          
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

% Define and correction matrix for left and right GH
Rleft = [1  0  0 
         0  0 -1
         0  1  0];  
Rright = [1  0  0  
          0  0  1
          0 -1  0];

% Define and correction matrix for left and right EL and WR
Rleft2 = [0  1  0 
          0  0 -1
          -1  0  0]; 
Rright2 = [0  -1 0          
           0  0  1
           -1 0 0]; 

% Define XSens joints to compile
L5S1 = 1:3; T9T8 = 10:12; C1Head = 16:18; rScapula = 19:21; rShoulder = 22:24; rElbow = 25:27; rWrist = 28:30; rHip = 43:45; lHip = 55:57; 
XSensJointsIdx = [L5S1; T9T8; C1Head; rScapula; rShoulder; rElbow; rWrist; rHip; lHip;]; 
XSensJointNames = {'L5S1' 'T9T8' 'C1Head' 'rScapula' 'rShoulder' 'rElbow' 'rWrist' 'rHip' 'lHip'};

%% Compute and Vizualize the data

for iS = 1:length(SubjectsToDo)
    fprintf([SubjectsToDo{iS} '...\n']);
    subject = SubjectsToDo{iS};
    if ismember(subject, {'P01', 'P06', 'P09', 'P11', 'P13', 'P14', 'P16', 'P19' 'P22' 'P23'}) % No P14
        Keyboard   = {'Norm'};
    elseif ismember(subject, {'P02', 'P03', 'P04', 'P05', 'P07', 'P08', 'P10', 'P12', 'P15', 'P17', 'P18', 'P20' 'P21' 'P24'})
        Keyboard   = {'Norm','CS60'};
    end

    for iK = 1:length(Keyboard)
        for iE = 1:length(Excerpts)
            for iC = 1:length(Cond)
                if ismember(subject, {'P01'}) && iE == 1
                    GoodTrials = {'A'};
                elseif ismember(subject, {'P02'}) && iK == 1 && iE == 3 && iC == 1
                    GoodTrials = {'A'};                    
                elseif ismember(subject, {'P10'}) && iK == 1 && iE == 1 && iC == 1
                    GoodTrials = {'B', 'C'};       
                elseif ismember(subject, {'P12'}) && iK == 1 && iE == 3 && iC == 1
                    GoodTrials = {'A', 'C'};                        
                elseif ismember(subject, {'P01', 'P02', 'P03', 'P04', 'P05', 'P06'})
                    GoodTrials = {'A', 'B'};
                else
                    GoodTrials = {'A', 'B', 'C'};                   
                end                
                for iG = 1:length(GoodTrials)
                    
                    if iS == 2 && iK == 1 && iE == 3 && iC == 1 && iG == 2
                        continue
                    else
                        % Compile and Organize all Orientation, Position, and Joint Angle data into a table                        
                        xsens_data = XSens.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).segmentData;
                        xsens_ja   = XSens.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).jointData;

                        orientation_all = [];
                        position_all    = [];
                        jointangle_all  = [];
                        for i = 1:23
                            orientation_data = xsens_data(i).orientation;
                            orientation_all  = [orientation_all, orientation_data]; 
                            
                            position_data    = xsens_data(i).position;
                            position_all     = [position_all, position_data];
                        end
                        for i = 1:22
                            jointangle_data = xsens_ja(i).jointAngle;
                            jointangle_all  = [jointangle_all, jointangle_data]; 
                        end
        
                        data.orientation_all = array2table(orientation_all);
                        data.position_all    = array2table(position_all);
                        data.jointangle_all  = array2table(jointangle_all);                
                        
                        nbSeg    = size(data.orientation_all,2)/4;       
                        nbFRall  = size(data.orientation_all,1);  
        
                        RT = nan(4,4,nbFRall,nbSeg);
                        RT(4,1,:)=0; 
                        RT(4,2,:)=0; 
                        RT(4,3,:)=0; 
                        RT(4,4,:)=1; 
        
                        for i = 1:nbSeg
                            % Load and compute q data
                            if i == 1 % 0 resets the pelvis angles to 0
                                 q_root = array2table(xsens_data(1).tpose.orientation);                                 
                                 q = data.orientation_all(:, i*4-3:i*4);
                                 q(1,1:4) = q_root;
                                 q = quaternion(q{:,:});
                                 R = q.RotationMatrix;
                                 R = squeeze(R);   
                            else
                                q = data.orientation_all(:, i*4-3:i*4);
                                q = quaternion(q{:,:});
                                R = q.RotationMatrix;
                                R = squeeze(R);

                                % correct axis of Upper Limb joints of participants
                                if i==rUpperArm                % for right GH  
                                    R = multiprod(R,Rright);
                                elseif i==rForeArm || i==rHand % for right EL WR
                                    R = multiprod(R,Rright2); %
                                end
%             
                                if i==lUpperArm                % for left GH   
                                    R = multiprod(R,Rleft); 
                                elseif i==lForeArm || i==lHand % for left EL WR 
                                    R = multiprod(R,Rleft2); %
                                end
                            end
                            
                            if i == 1 % 0 resets the pelvis angles to 0
                                % Comment based on the visual options
                                % wanted - does not change how joint angles
                                % are calculated (only the translation part
                                % in global space)
                                % -- A. keep world-space coordinates (recommended) -------------
                                T = data.position_all(:,1:3);
                                T = T{:,:}.';                % 3×N
                                T = reshape(T,3,1,[]);
                            
                                % -- B. first-frame–centred  -----------------------------------
                                % T_raw   = data.position_all(:,1:3);
                                % T_shift = T_raw - T_raw(1,:);
                                % T       = T_shift.';         % 3×N
                                % T       = reshape(T,3,1,[]);
                            
                                % -- C. neutral-centred ----------------------------------------
                                % T_root  = xsens_data(1).tpose.position;
                                % T_raw   = data.position_all(:,1:3);
                                % T_shift = T_raw - T_root;
                                % T       = T_shift.';         % 3×N
                                % T       = reshape(T,3,1,[]);
                            else                    
                                T = data.position_all(:, i*3-2:i*3);
                                T = T{:,:}';
                                T = reshape(T,3,1,[]);
                            end
                            RT(1:3,1:4,:,i) = [R T];
                        end

                        % Plot the axes at each joint
                        if plot_axes == true
                            % Plot axis
                            frameTP    = 500;
                            Fig1 = figure;
                            for iSeg = 1:length(segmentsToUse)
                                % Define segment index
                                idx = segmentsToUse(iSeg);
                                % plot axis
                                plotAxes(RT(:,:,frameTP,idx))
                                hold on
                            end
                            
                            seg = Joints(2:end,:);
            
                            for i=1:size(seg,1)
                                hold on
                                plot3(squeeze(RT(1,4,frameTP, seg(i,:))), squeeze(RT(2,4,frameTP, seg(i,:))), squeeze(RT(3,4,frameTP, seg(i,:))), 'k','LineWidth', 3)
                                hold on
                                plot3(squeeze(RT(1,4,frameTP, seg(i,1))), squeeze(RT(2,4,frameTP, seg(i,1))), squeeze(RT(3,4,frameTP, seg(i,1))), 'ko','LineWidth', 3)
                                hold on
                            end
                            axis equal
                        end
        
%                         % 1) Check that R is non-identity
%                         R = RM.(SubjectsToDo{iS}).(Keyboard{iK});
%                         disp( ' R = ' ), disp(R)
%                         if norm(R - eye(3))<1e-6
%                             warning('Your R is effectively identity → no change will happen');
%                         end
%                         
%                         % 2) Check your RT size
%                         sz = size(RT);
%                         fprintf('RT is %dx%dx%dx%d\n', sz);

                        % Compute joint angles
                        nbJ = size(Joints,1);
                        JointAngles = nan(nbFRall,3,nbJ);     
                        for j = 1:nbJ
                            if isnan(Joints(j,1))
                                RJ = multiprod(RT(1:3,1:3,:,Joints(j,2)), multitransp(RT(1:3,1:3,1,Joints(j,2)))); %multitransp(RT(1:3,1:3,1,Joints(j,2))));
                            else
                                RJ = multiprod(RT(1:3,1:3,:,Joints(j,2)), multitransp(RT(1:3,1:3,:,Joints(j,1)))); %multitransp(RT(1:3,1:3,:,Joints(j,1))));
                            end
                            JointAngles(:,:,j) = squeeze(angleRotation(RJ, Sequences(j,:)))';
                            if Joints(j,2) == rUpperArm || Joints(j,2) == lUpperArm || Joints(j,2) == rForeArm || Joints(j,2) == lForeArm 
                                JointAngles(:,:,j) = unwrap(JointAngles(:,:,j)); % pour corriger les sauts de pi potentiels à l'épaule et au coude
                            end
                            JointAngles(:,:,j) = rad2deg(JointAngles(:,:,j));
                        end
        
                        % Compile computed angles
                        for iJ = 1:nbJ
                            JointAngles_Calculated.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).(jointsNames{iJ}).jointAngle = squeeze(JointAngles(:,:,iJ));
                        end
                        
                        % Compile XSens calculated Joint Angles
                        XSens_Angles.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).JA     = xsens_ja;

                        % Compile XSens Segment data
                        xsens_seg = XSens.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).segmentData;
                        XSens_Segments.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).SEG = xsens_seg;
                    
                        pl = squeeze(JointAngles(:,:,1));
                        th = squeeze(JointAngles(:,:,2));

                        pl_new = [];
                        th_new = [];

                        % Correct order to match XSens (Change order to Y, Z, X)
                        pl_new(:,1) = pl(:,2); % y
                        pl_new(:,2) = pl(:,3); % z
                        pl_new(:,3) = pl(:,1); % x

                        th_new(:,1) = th(:,2); % y
                        th_new(:,2) = th(:,3); % z
                        th_new(:,3) = th(:,1); % x                        

                        % Store the calculated thorax and pelvis angles with the XSens joint angles
                        JointAngles_Compiled.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).PL = pl_new;
                        JointAngles_Compiled.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).TH = th_new;
                        JointAngles_Compiled.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).NK = xsens_ja(6).jointAngle;
                        JointAngles_Compiled.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).rGH = xsens_ja(8).jointAngle;
                        JointAngles_Compiled.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).rEL = xsens_ja(9).jointAngle;
                        JointAngles_Compiled.(SubjectsToDo{iS}).(Keyboard{iK}).(Excerpts{iE}).(Cond{iC}).(GoodTrials{iG}).rWR = xsens_ja(10).jointAngle; 
                    end
                end
            end
        end
    end
end

%% Save

if save_data == true
    cd('D:\turne\uleth\PhD\Thesis\Data_Collection\DATA\Proximal_Kinematics')
    save('JointAngles_Compiled.mat','JointAngles_Compiled')
end
