clear

%close all;

Site={'NYU'};

ModeList={'func_minimal','nofilt_noglobal'};%,'Unproc','nofilt_global'};
%ModeList={'nofilt_noglobal'};
%25 Healthy 
%SubList={'0051036','0051038','0051039','0051040','0051041','0051042','0051044','0051045','0051046'...
%    ,'0051047','0051048','0051049','0051050','0051051','0051052','0051053','0051054','0051055',...
%    '0051056','0051057','0051058','0051059','0051060','0051061','0051062'};

SubList={'0051050','0051055'};

WhoToBelieve='PCPWebsite';

%SubList={'0051055'};


GSRStatList={'noglobal'}; %%BE CAREFULLLLL****************
PPline0={'cpac'}; %
addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC
addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC/Nifti_Util
addpath /Users/sorooshafyouni/Home/matlab/fieldtrip-20160107/external/afni

TestMeth={'X2_m3d3'}; %'X2d3','Z'

rt='/Users/sorooshafyouni/Home/PCP/Vols/Outputs/';
%rt='/Volumes/HCP_S900/PCP/Vols/Outputs/';

addpath Nifti_Util/
addpath /Users/sorooshafyouni/Home/GitClone/DVARS
%addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC_NowGitClone

TickScaler=1;

for g=GSRStatList
    for s=SubList
        %close all
        for m=ModeList
            PPline1=PPline0;
            disp(['#-------------' s{1} ' -- ' m{1} ' -- ' g{1}])

            read_dir=['R/' Site{1} '_' s{1} '_fMRIDiag'];
            if exist(read_dir,'dir')~=7; error('Akh!'); end;

            if strcmp(m,'nofilt_noglobal')
                InstyNrm={'norm',100};
                OneSub=[rt PPline1{1} '/' m{1} '/func_preproc/' Site{1} '_' s{1} '_func_preproc.nii.gz'];
            elseif strcmp(m,'func_minimal')
                InstyNrm={'norm',100};
                OneSub=[rt PPline1{1} '/func_minimal/' Site{1} '_' s{1} '_func_minimal.nii.gz'];
            elseif strcmp(m,'Unproc')
                InstyNrm={'scale',1/10};
                PPline1={'raw'};
                OneSub=[rt 'raw/' Site{1} '_' s{1}(3:end) '/' Site{1} '_' s{1} '_rest_Brain.nii.gz'];
            else
                error('Unknown mode arg!')
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
            
            T0=176;
            %read the dude
            V1 = load_untouch_nii(OneSub);
            V2 = V1.img; 
            X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
            I0 = prod([X0,Y0,Z0]);
            Y  = reshape(V2,[I0,T0]); clear V2 V1; 
            
            %
            load([read_dir '/' Site{1} '_' s{1} '_' g{1} '_' m{1} '_Council_' WhoToBelieve '.mat']);
            %[FDts,FD_Stat]=FDCalc(MovPar);
            
            %mean(V.MeanOrig)
            
%             pvals  = eval(['DVARS_Stat_' TestMeth{1} '.pvals']);
%             Idx=find(pvals<(0.05./(T0-1)));
            
            pvals  = eval(['DVARS_Stat_' TestMeth{1} '.pvals']);
            Idxs=find(pvals<(0.05./(T0-1)));
            
            DpDvar  = eval(['DVARS_Stat_' TestMeth{1} '.DeltapDvar']);
            Idxp   = find(pvals<(0.05./(T0-1)) & DpDvar>5);
            
            Idxs=Idxs+1; Idxp=Idxp+1;
            
            if      contains(TestMeth{1},'X2')
                NDVARS = eval(['DVARS_Stat_' TestMeth{1} '.SDVARS_X2']);
            elseif  contains(TestMeth{1},'Z')
                NDVARS = eval(['DVARS_Stat_' TestMeth{1} '.SDVARS_Z']);
            else
                error('Unknown test method: choose btwn: X2d3, X2d1, Z');
            end
            
            ColRng=[-3 3];
            f_hndl=figure('position',[50,500,600,500]); hold on;
            hsp0=suptitle([s{1} '-' m{1} '-' TestMeth{1}]);
            hsp0.Interpreter='none';
            hsp0.FontSize=9;
            %BOLDImgPan(Y,'scale',1/100,'destdir',read_dir,'handle',f_hndl,'FD',FDts,'DVARS',NDVARS,'ColRng',ColRng,'prefix',prefix,'Idx',Idx);

            %TfMRIDiag_plot(V,'FD',FDts,InstyNrm{1},InstyNrm{2},'gsrflag',gsrflag,'Idx',Idx,'handle',f_hndl,'Thick',0.5,'fontsize',14,'linewidth',1.2)            
            
            %fMRIDiag_plot_PP(V,'TickScaler',1/2,'grandmean',100,'bold',Y,'ColRng',ColRng,'FD',FDts,InstyNrm{1},InstyNrm{2},'gsrflag',gsrflag,'Idx',Idx,'handle',f_hndl,'Thick',0.5,'fontsize',14,'linewidth',1.2)
            fMRIDiag_plot_4paper(V,'TickScaler',TickScaler,'bold',Y,'ColRng',ColRng,'FD',FDts,InstyNrm{1},InstyNrm{2},'gsrflag',gsrflag,'Stat_Idx',Idxs,'Pract_Idx',Idxp,'handle',f_hndl,'fontsize',14,'linewidth',1.2)

            %fMRIDiag_plot(V,'bold',Y,'ColRng',ColRng,'FD',FDts,'scale',1/100,'gsrflag',gsrflag,'Idx',Idx,'handle',f_hndl,'fontsize',14,'linewidth',1.2)

            hold off; 
            
            prefix=[s{1} '_' m{1} '_' TestMeth{1}];
            export_fig(f_hndl,['Figs/fMRIDiag_' gsrflagstr '_ColRng' num2str(ColRng(2)) '_' prefix '.pdf'])
            
            %close all
            clear Y prefix I0 T0 f_hndl V FDts MovPar *FD* DVARS_* *Stat* pvals NDVARS InstyNrm
        end
    end
end