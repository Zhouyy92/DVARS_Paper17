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
%               
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
addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC
addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC/Nifti_Util

for GSRSFlga=GSRStatList
    for s=SubList
        for m=ModeList
            %-----------------------------------D.Find and Read-----------------------
            disp(['#-------------' s{1} ' -- ' m{1} ' -- ' GSRSFlga{1} ' --------------------------#'])

            output_dir=['R/' Site{1} '_' s{1} '_fMRIDiag'];
            if exist(output_dir,'dir')~=7; mkdir(output_dir); end;

            if strcmp(m,'Pre_Fix')
                OneSub=['/Volumes/HCP_S900/HCP_10Unrel_Vols/' s{1} '/rfMRI_REST1_LR.nii.gz'];
            elseif strcmp(m,'Post_Fix')
                OneSub=['/Volumes/HCP_S900/HCP_10Unrel_Vols/' s{1} '/rfMRI_REST1_LR_hp2000_clean.nii.gz'];
            elseif strcmp(m,'Unproc')
                OneSub=['/Volumes/HCP_S900/HCP_10Unrel_Vols/' s{1} '/' s{1} '_3T_rfMRI_REST1_LR_Brain.nii.gz'];
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
            MovPar=MovPartextImport(['/Volumes/HCP_S900/HCP_10Unrel_Vols/' s{1} '/' s{1} '/MNINonLinear/Results/rfMRI_REST1_LR/Movement_Regressors.txt']);
            
            [FDts,FD_Stat]=FDCalc(MovPar);
            
            if strcmpi(m,'Unproc')
                InstyNrm={'norm',100};
            elseif strcmpi(m,'Pre_Fix') || strcmpi(m,'Post_Fix')
                InstyNrm={'scale',1/100};
            else
                error('DUDDDE!')
            end
            
            %------------------------------F.Var and Inf---------------------
            [V,DSE_Stat]=DSEvars(Y, 'gsrflag',gsrflag, InstyNrm{1},InstyNrm{2} ,'verbose',1);

            [DVARS_X2_m1d3, DVARS_Stat_X2_m1d3] = DVARSCalc(Y,'TestMeth','X2' , 'WhichExpVal',1 , 'TransPower',1/3, 'VarRobIQR','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1 ,'nflag',0);
            [DVARS_X2_m1d1, DVARS_Stat_X2_m1d1] = DVARSCalc(Y,'TestMeth','X2' , 'WhichExpVal',1 , 'TransPower',1  , 'VarRobIQR','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1 ,'nflag',0);
            
            [DVARS_X2_m3d3, DVARS_Stat_X2_m3d3] = DVARSCalc(Y,'TestMeth','X2' , 'WhichExpVal',3 , 'TransPower',1/3, 'VarRobIQR','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1 ,'nflag',0);
            [DVARS_X2_m3d1, DVARS_Stat_X2_m3d1] = DVARSCalc(Y,'TestMeth','X2' , 'WhichExpVal',3 , 'TransPower',1  , 'VarRobIQR','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1 ,'nflag',0);            
            
            [DVARS_Z   , DVARS_Stat_Z   ] = DVARSCalc(Y,'TestMeth','Z'  , 'WhichExpVal',3 , 'TransPower',1  , 'VarRobIQR','hIQR' , 'gsrflag',gsrflag , InstyNrm{1},InstyNrm{2} , 'verbose',1 ,'nflag',1);
            
            
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
