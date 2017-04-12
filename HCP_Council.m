% HCP_Council.m
% Calcs the DSE var decomposition & different combinations of the DVARS
% inference. 
%
% This code was written for HCP data-sets. However, can be used as a
% general template for other datasets by changing section D of the code. 
%
% Dependencies: 
%               DVARSCalc.m
%               /Nifti_Util
%               DSEvars.m
%               FDCalc.m
%               MovPartextImport.m %pretty useless though, alt can be used
%
%               These files can be found in https://github.com/asoroosh/DVARS
%
%%%NOTE:
%       Results of this file for 115320 can be found here:
%       https://github.com/asoroosh/DVARS_Paper17/tree/master/Data/HCP_115320_fMRIDiag
%
%       And likewise for 118730:
%       https://github.com/asoroosh/DVARS_Paper17/tree/master/Data/HCP_118730_fMRIDiag
%
%%%REFERENCES
%
%   Afyouni S. & Nichols T.E., Insights and inference for DVARS, 2017
%   http://www.biorxiv.org/content/early/2017/04/06/125021.1               
%
% SA, 2017, UoW
% srafyouni@gmail.com


clear; 

close all;

%------- A.Where from?! ------------
Site={'HCP'};

ModeList={'Unproc','Pre_Fix','Post_Fix'};
GSRStatList={'noglobal'}; %%BE CAREFULLLLL****************
%------- B.Subjects ----------------
%21 Healthy 
%SubList={'100307','103414','105115','111312','113619','115320','117122','118730','123117','151526','187345','303624','132118','901139','171330','263436','191336','779370','195041','145127','172029'};

SubList={'115320'};%,'115320'};%,'118730'}; %of interests

%------- C.addpaths ----------------
% You need to clone the DVARS directory:
% https://github.com/asoroosh/DVARS.git
% and addpath the directory (+ subfolders)

%addpath ~/DVARS/
%addpath ~/DVARS/Nifti_Util

%Change this path to the path of the directory where you store your HCP
%data. 
rt='/Volumes/HCP_S900/HCP_10Unrel_Vols/';

for GSRSFlga=GSRStatList
    for s=SubList
        for m=ModeList
            %-----------------------------------D.Find and Read-----------------------
            disp(['#-------------' s{1} ' -- ' m{1} ' -- ' GSRSFlga{1} ' --------------------------#'])

            output_dir=['R/' Site{1} '_' s{1} '_fMRIDiag'];
            if exist(output_dir,'dir')~=7; mkdir(output_dir); end;

            if strcmp(m,'Pre_Fix')
                OneSub=[rt '/rfMRI_REST1_LR.nii.gz']; %Minimally Pre-processed
            elseif strcmp(m,'Post_Fix')
                OneSub=[rt '/rfMRI_REST1_LR_hp2000_clean.nii.gz']; %Fully pre-processed
            elseif strcmp(m,'Unproc')
                OneSub=[rt  s{1} '_3T_rfMRI_REST1_LR_Brain.nii.gz']; %RAW
            end

            if strcmpi(GSRSFlga{1},'global')
                gsrflag=1; 
                disp('**GSR flag is on!');
            elseif strcmpi(GSRSFlga{1},'noglobal')
                gsrflag=0;
            else
                error('Unknown gsr flag!')
            end;
            
            %read the dude
            V1 = load_untouch_nii(OneSub);
            V2 = V1.img; 
            X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
            I0 = prod([X0,Y0,Z0]);
            Y  = reshape(V2,[I0,T0]); clear V2 V1; 
            %
            %-----------------------------------E.FD-----------------------
            % Change this to path of where you keep your
            % Movement_Regressors.txt files
            Path_to_Movement_Reg = ['/Volumes/HCP_S900/HCP_10Unrel_Vols/' s{1} '/' s{1} '/MNINonLinear/Results/rfMRI_REST1_LR/Movement_Regressors.txt'];
            MovPar=MovPartextImport(Path_to_Movement_Reg);
            
            [FDts,FD_Stat]=FDCalc(MovPar);
            
            if strcmpi(m,'Unproc')
                InstyNrm={'norm',100};
            elseif strcmpi(m,'Pre_Fix') || strcmpi(m,'Post_Fix')
                InstyNrm={'scale',1/100};
            else
                error('DUDDDE!')
            end
            
            %------------------------------F.Var and Inf---------------------
            
            %DSE Variance decomposition
            [V,DSE_Stat]=DSEvars(Y, 'gsrflag',gsrflag, InstyNrm{1},InstyNrm{2} ,'verbose',1);

            %Different test combinations 
            [DVARS_X2_m1d3, DVARS_Stat_X2_m1d3] = DVARSCalc(Y,'TestMethod','X2' , 'MeanType',1 , 'TransPower',1/3, 'VarType','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1);
            [DVARS_X2_m1d1, DVARS_Stat_X2_m1d1] = DVARSCalc(Y,'TestMethod','X2' , 'MeanType',1 , 'TransPower',1  , 'VarType','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1);
            
            [DVARS_X2_m3d3, DVARS_Stat_X2_m3d3] = DVARSCalc(Y,'TestMethod','X2' , 'MeanType',3 , 'TransPower',1/3, 'VarType','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1);
            [DVARS_X2_m3d1, DVARS_Stat_X2_m3d1] = DVARSCalc(Y,'TestMethod','X2' , 'MeanType',3 , 'TransPower',1  , 'VarType','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1);            
            
            %This one does the Relative DVARS as well. 
            [DVARS_Z   , DVARS_Stat_Z   ] = DVARSCalc(Y,'TestMeth','Z'  , 'WhichExpVal',3 , 'TransPower',1  , 'VarRobIQR','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1 ,'RDVARS');
            
            
            %f_hndl=figure('position',[50,500,1600,500]); hold on;
            %BOLDImgPan(Y,'scale',1/100,'destdir',output_dir,'handle',f_hndl,'FD',FDts,'DVARS',DVARS_Stat_Z.NDVARS_Z,'ColRng',[-10 10],'gsrflag',gsrflag,'verbose',0);
            %hold off; clear f_hndl;
            
            %close all
            %------------------------------G.Save----------------------------
            save([output_dir '/' Site{1} '_' s{1} '_' GSRSFlga{1} '_' m{1} '_Council_Test.mat'],...
                'FDts' , 'MovPar'   , 'FD_Stat',...
                'V'    , 'DSE_Stat' ,...
                'DVARS_X2_m1d3' , 'DVARS_Stat_X2_m1d3',...
                'DVARS_X2_m1d1' , 'DVARS_Stat_X2_m1d1',...
                'DVARS_X2_m3d3' , 'DVARS_Stat_X2_m3d3',...
                'DVARS_X2_m3d1' , 'DVARS_Stat_X2_m3d1',...
                'DVARS_Z'       , 'DVARS_Stat_Z');

            clear Y *0 V* FDts MovPar DVARS_* *Stat* OneSub output_dir
        end
    end
end
