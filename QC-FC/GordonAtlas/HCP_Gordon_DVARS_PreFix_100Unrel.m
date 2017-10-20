clear 
close all


Atlas={'Gordon'};
Mode={'PreFix'};

load(['HCP_100Unrel_SubList.mat'])
HCP_10Unrel_SubList(3)=[];

SubList=HCP_10Unrel_SubList;

addpath /Users/sorooshafyouni/Home/GitClone/DVARS/Nifti_Util
addpath /Users/sorooshafyouni/Home/GitClone/DVARS
addpath /Users/sorooshafyouni/Home/GitClone/DVARS/Aux

nRlz=1000;

thr_fd_con=0.2;
thr_fd_lib=0.5;
thr_dv=5;

scl=1/100;
%--------------------------------------------------------------
cnt_s=1;
for s=SubList
    disp([num2str(cnt_s) ' ' s{1} '------'])
    load(['/Volumes/HCP_S900/HCP_100Unrelated/' Mode{1} '_NetMats/' Atlas{1} '/ts/HCP_' Mode{1} '_' s{1} '_' Atlas{1} '_ts.mat'],'ats')
    %ats=ats(150,:);
    %----------------------------------------------------------    
    [DVARS,DVARS_Stat]  = DVARSCalc(ats','Scale',scl,'verbose',0);
    dv_rmv_st             = find(DVARS_Stat.pvals<(0.05./1199));
    dv_rmv_pr              = dv_rmv_st(ismember(dv_rmv_st,find(DVARS_Stat.DeltapDvar>=5)));
%     dv_rmv_dp           = unique([dv_rmv dv_rmv+1]);
    %----------------------------------------------------------
    [V,DSE_Stat]        = DSEvars(ats','Scale',scl,'verbose',0);
    %----------------------------------------------------------
%    load(['R/NetMats/HCP_MinimalF_' s{1} '_Power_NetMats.mat'],'FDts','FD_Stat'); %Because I left the Hard back in the office!!
     MovPar=MovPartextImport(['/Volumes/HCP_S900/HCP_100Unrelated/' Mode{1} '/' s{1} '/' s{1} '/MNINonLinear/Results/rfMRI_REST1_LR/Movement_Regressors.txt']);
     [FDts,FD_Stat]=FDCalc(MovPar,[0 0 0 1 1 1],50,0);
    %----------------------------------------------------------
    fd_rmv_con=find(FDts>thr_fd_con);
    fd_rmv_dp_con=unique([fd_rmv_con fd_rmv_con+1]);
    
    fd_rmv_lib=find(FDts>thr_fd_lib);
    fd_rmv_dp_lib=unique([fd_rmv_lib fd_rmv_lib+1]);
    %----------------------------------------------------------
    load(['/Volumes/HCP_S900/HCP_100Unrelated/' Mode{1} '_NetMats/' Atlas{1} '/mts/HCP_' Mode{1} '_' s{1} '_' Atlas{1} '_ROIs.mat'],'mts','ROI')
    %----------------------------------------------------------
    for r0=1:size(ROI.ts,2)
        c0=mean(ROI.XYZ{r0});
        for r1=1:size(ROI.ts,2)
            c1=mean(ROI.XYZ{r1});
            ED(r0,r1)=sqrt(sum((c0-c1).^2));
        end
    end
    %----------------------------------------------------------
    mat=corr(mts);
    %----------------------------------------------------------
    %rgs_key = unique([dv_rmv-3 dv_rmv-2 dv_rmv-1 dv_rmv dv_rmv+1 dv_rmv+2 dv_rmv+3]);
    rgs_key = unique([dv_rmv_pr dv_rmv_pr+1]);%unique([dv_rmv]);
    rgs_key(rgs_key>length(DVARS))=[];
    
    conNlib=intersect(fd_rmv_con,fd_rmv_lib);
    conDlib=setdiff(fd_rmv_con,fd_rmv_lib);
    
    dvthr_rmv=find(DVARS>thr_dv);
    
    dvQthr_rmv=DVARS_IQR(DVARS,'DVARS');
    
    mat_dvs_scrmbld=zeros(size(mat,1),size(mat,2),nRlz);
    mat_fd_con=[]; mat_fd_lib=[]; mat_dv_thr=[]; mat_dv_Qthr=[];
    %----------------------------------------------------------
    mts0                      = mts;
    if ~isempty(fd_rmv_dp_con); mts0(fd_rmv_dp_con,:)=[];   end;
    mat_fd_con = corr(mts0);
    clear mts0
    %----------------------------------------------------------
    mts0                      = mts;
    if ~isempty(fd_rmv_dp_lib); mts0(fd_rmv_dp_lib,:)=[];   end;
    mat_fd_lib = corr(mts0);
    clear mts0 
    %----------------------------------------------------------
    mts0                 = mts;
    if ~isempty(dvthr_rmv);     mts0(dvthr_rmv,:) = [];     end;
    mat_dv_thr = corr(mts0);
    clear mts0
%----------------------------------------------------------
    mts0 = mts;
    if ~isempty(dvQthr_rmv);    mts0(dvQthr_rmv,:) = [];    end;
    mat_dv_Qthr = corr(mts0);
    clear mts0
%----------------------------------------------------------    
    if ~isempty(rgs_key)
    
        rgsDM=[];
        for i=1:numel(rgs_key) 
            imhere=zeros(size(ats,1),1);
            imhere(rgs_key(i))=1;%DSE_Stat.DeltapDvar(rgs_key(i)); 
            rgsDM=[rgsDM,imhere];
            clear imhere
        end

        %rgs=[0 DSE_Stat.DeltapDvar()]';
        %zeros(size(ats,1),1)

        dvr_ats=RegDVARS(ats',rgsDM);%,'Scale',scl);
        [DVARS_DVC_reg,DVARS_Stat_DVC_reg]  = DVARSCalc(dvr_ats','Scale',scl,'verbose',0);
        [V_DVC_reg,DSE_Stat_DVC_reg]        = DSEvars(dvr_ats','Scale',scl,'verbose',0);
        %----------------------------------------------------------
        dvs_ats=ats;
        dvs_ats(rgs_key,:)=[];
        [DVARS_DVC,DVARS_Stat_DVC]  = DVARSCalc(dvs_ats','Scale',scl,'verbose',0);
        [V_DVC,DSE_Stat_DVC]        = DSEvars(dvs_ats','Scale',scl,'verbose',0);       
        %----------------------------------------------------------
        mts0                = mts;
        mts0(rgs_key,:)     = [];
        mat_dvs             = corr(mts0);
        clear mts0
        
        disp('Scrambling....')
        for si=1:nRlz
            mts0 = mts;
            rgs_key_scmbld=datasample(1:size(ats,1),numel(rgs_key));
            mts0(rgs_key_scmbld,:) = [];
            mat_dvs_scrmbld(:,:,si) = corr(mts0);
            clear mts0 rgs_key_scmbld
        end
        %----------------------------------------------------------
        mts0                = mts;
        mts0                = RegDVARS(mts0',rgsDM);
        mat_dvr             = corr(mts0);
        clear mts0
        %----------------------------------------------------------
    %     ats_DVC=ats;
    %     ats_DVC(dv_rmv_dp,:)=[];
        %----------------------------------------------------------
    else
        mat_dvr=[]; mat_dvs=[]; DSE_Stat_DVC=[];
        V_DVC=[]; DVARS_DVC=[]; DVARS_Stat_DVC=[];
        V_DVC_reg=[]; DSE_Stat_DVC_reg=[];
        DVARS_DVC_reg=[]; DVARS_Stat_DVC_reg=[]; rgsDM=[];
    end
    save(['/Volumes/HCP_S900/HCP_100Unrelated/' Mode{1} '_NetMats/' Atlas{1} '/NetMats/HCP_' Mode{1} '_' s{1} '_' Atlas{1} '_NetMats_DVC.mat'],'mat_dvs_scrmbld','rgs_key','mat_dv_Qthr','dvQthr_rmv','dvthr_rmv','rgsDM','V_DVC','DSE_Stat_DVC','DVARS_DVC','DVARS_Stat_DVC','thr_fd_con','thr_fd_lib','V_DVC_reg','DSE_Stat_DVC_reg','DVARS_DVC_reg','DVARS_Stat_DVC_reg','mat','mat_dvr','mat_dvs','mat_dv_thr','mat_fd_lib','fd_rmv_dp_lib','fd_rmv_lib','mat_fd_con','fd_rmv_dp_con','fd_rmv_con','dv_rmv_pr','dv_rmv_st','DVARS','DVARS_Stat','FDts','FD_Stat','ED','DSE_Stat','V','conNlib','conDlib')
    
%     hh0=figure('position',[50,500,1500,600]); 
%     subplot(2,1,1); hold on; box on;
%     title(['DVARS - ' s{1}])
%     plot(DVARS,'linewidth',1.2,'color',[.8 .8 .8])
%     if ~isempty(rgs_key)
%         plot(DVARS_DVC_reg,'linewidth',1.2,'color',[0.6 .1 .1])
%         %------
%         tt=zeros(1,numel(DVARS));
%         tt(rgs_key)=1;
%         dv_idx_b=find(~tt);
%         DVARS_DVC_b(dv_idx_b)=DVARS_DVC;
%         DVARS_DVC_b(rgs_key)=NaN;
%         plot(DVARS_DVC_b,'linewidth',1.2,'color',[.2 .2 .4])
%         legend({'DVARS','DVARS (Rgsd out)','DVARS (Scrubbed)'})
%         clear tt *_b
%     end
%     ylim0s1=ylim;
% %     if ~isempty(dvthr_rmv)
% %         for i=1:numel(dvthr_rmv)
% %             line([dvthr_rmv(i) dvthr_rmv(i)],[ylim0s1(1) ylim0s1(2)],'color','r')
% %         end
% %     end
%     
%     if ~isempty(dvQthr_rmv)
%         for i=1:numel(dvQthr_rmv)
%             line([dvQthr_rmv(i) dvQthr_rmv(i)],[ylim0s1(1) ylim0s1(2)],'color','r')
%         end
%     end
%  %---------------------------------
%     subplot(2,1,2); hold on; box on; 
%     title(['FD - ' s{1}]) 
%     plot(FDts,'linewidth',1.2,'color',[1 0 0])
%     ylim02=ylim;
%         
%     for i=1:numel(conNlib)
%         line([conNlib(i) conNlib(i)],[ylim02(1) ylim02(2)],'color','k')
%     end
%     
%     for i=1:numel(conDlib)
%         line([conDlib(i) conDlib(i)],[ylim02(1) ylim02(2)],'color','g')
%     end
%  %---------------------------------    
%     export_fig(hh0,['DVARSplots_100Unrel/DVARS_plots_' s{1} '.pdf'])
%     
%     pause(2)
%     
%     close all
%     
    cnt_s=cnt_s+1;
    clear ROI  ED FD* rmv* fd_* dv_*
    clear mat* fd_* dv_* mts* *_reg DVARS* DSE* rgs*
    clear mat_dvs_scrmbld
end