clear; 

%close all;

Site={'HCP'};

%ModeList={'Unproc','Pre_Fix','Post_Fix'};
ModeList={'Pre_Fix','Post_Fix'};

%21 Healthy 
%SubList={'100307','103414','105115','111312','113619','115320','117122','118730','123117','151526','187345','303624','132118','901139','171330','263436','191336','779370','195041','145127','172029'};
SubList={'118730'};%,'115320'};%,'115320'};%,'118730'};

GSRStatList={'noglobal'}; %%BE CAREFULLLLL****************

TestMeth={'X2_m3d3'}; %'X2d3','Z'
%Choose:  X2d1 Z
addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC
%addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC/Nifti_Util
addpath /Users/sorooshafyouni/Home/matlab/Ext/cifti-matlab-master

for g=GSRStatList
    for s=SubList
        %close all
        for m=ModeList

            disp(['#-------------' s{1} ' -- ' m{1} ' -- ' g{1}])

            read_dir=['R/' Site{1} '_' s{1} '_fMRIDiag'];
            if exist(read_dir,'dir')~=7; error('Dir doesnt exists!'); end;

            CIFTIrt=['/Volumes/HCP_S900/HCP_10Unrel_Vols/' s{1} '/' s{1} '/MNINonLinear/Results/rfMRI_REST1_LR'];
            if strcmp(m,'Pre_Fix')
                OneSub=[CIFTIrt '/rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii'];
                TickScaler=1/2;
            elseif strcmp(m,'Post_Fix')
                OneSub=[CIFTIrt '/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'];
                TickScaler=1;
%            elseif strcmp(m,'Unproc')
%                OneSub=['/Volumes/HCP_S900/HCP_10Unrel_Vols/' s{1} '/' s{1} '_3T_rfMRI_REST1_LR_Brain.nii.gz'];
            end

            if strcmpi(g{1},'global')
                gsrflag=1; 
                disp('**GSR flag is on!');
                gsrflagstr='GSR';
            elseif strcmpi(g{1},'noglobal')
                gsrflag=0;
                gsrflagstr='NoGSR';
            else
                error('Unknown gsr flag!')
            end;
            
            %read the dude -- NIFTI
%             V1 = load_untouch_nii(OneSub);
%             V2 = V1.img; 
%             X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
%             I0 = prod([X0,Y0,Z0]);
%             Y  = reshape(V2,[I0,T0]); clear V2 V1; 
            %--CIFTI
            V1 = ft_read_cifti(OneSub);
            V2=V1.dtseries;
            I0=size(V2,1); T0=size(V2,2);
            Y=V2; clear V2 V1; 
            %
            load([read_dir '/' Site{1} '_' s{1} '_' g{1} '_' m{1} '_Council.mat']);
            %[FDts,FD_Stat]=FDCalc(MovPar);
            
            pvals  = eval(['DVARS_Stat_' TestMeth{1} '.pvals']);
            Idx=find(pvals<(0.05./(T0-1)));
            
            if      contains(TestMeth{1},'X2')
                NDVARS = eval(['DVARS_Stat_' TestMeth{1} '.NDVARS_X2']);
            elseif  contains(TestMeth{1},'Z')
                NDVARS = eval(['DVARS_Stat_' TestMeth{1} '.NDVARS_Z']);
            else
                error('Unknown test method: choose btwn: X2d3, X2d1, Z');
            end
            
            ColRng=[-3 3];
            %f_hndl=figure('position',[50,500,1600,1400]); hold on;
            f_hndl=figure('position',[50,500,600,500]); hold on;
            hsp0=suptitle([s{1} '-' m{1} '-' TestMeth{1}]);
            hsp0.Interpreter='none';
            hsp0.FontSize=9;
            %BOLDImgPan(Y,'scale',1/100,'destdir',read_dir,'handle',f_hndl,'FD',FDts,'DVARS',NDVARS,'ColRng',ColRng,'prefix',prefix,'Idx',Idx);
            
            %'AbsMov',[FD_Stat.AbsRot FD_Stat.AbsTrans]
            fMRIDiag_plot(V,'GrandMean',mean(V.MeanOrig),'TickScaler',TickScaler,'bold',Y,'ColRng',ColRng,'FD',FDts,'scale',1/100,'gsrflag',gsrflag,'Idx',Idx,'handle',f_hndl,'fontsize',14,'linewidth',1.2)
            
            hold off; 
            
            prefix=[s{1} '_' m{1} '_' TestMeth{1}];
            %export_fig(f_hndl,['fMRIDiag_' gsrflagstr '_ColRng' num2str(ColRng(2)) '_' prefix '.pdf'])
            
            %close all
            clear Y prefix I0 T0 f_hndl V FDts MovPar *FD* DVARS_* *Stat* pvals NDVARS
        end
    end
end