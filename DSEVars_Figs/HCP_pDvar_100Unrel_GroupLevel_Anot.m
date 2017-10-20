% This function was used to produce Figure 7.
%
% NB! You have to run HCP_Council.m for all 20/21 subjects beforehand and
% save the results. The results for all 21 HCP subjects are available to be
% shared upon request at srafyouni@gmail.com
%
% 
%
% SA, UoW, 2017
% srafyouni@gmail.com 
%
%%%REFERENCES
%
%   Afyouni S. & Nichols T.E., Insights and inference for DVARS, 2017
%   http://www.biorxiv.org/content/early/2017/04/06/125021.1  
%
%   SA, UoW, 2017
%   srafyouni@gmail.com
%

clear; close all;

Site={'HCP'};

PPLine={''};

addpath /Users/sorooshafyouni/Home/GitClone/DVARS/

%ModeList={'Unproc','Pre_Fix','Post_Fix'};
ModeList={'Unproc','PreFix','PostFix'};

%21 Healthy 
%SubList={'100307','103414','105115','111312','113619','115320','117122','118730','123117','151526','187345','303624','132118','901139','171330','263436','191336','779370','195041','145127','172029'};
%SubList={'118730'};%,'115320'};%,'118730'};
load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/HCP/QC-FC/PowerAtlas/HCP_100Unrel_SubList.mat'])
SubList=HCP_10Unrel_SubList;
%------------------------------------------------------------------
GSRStat={'noglobal'}; %%BE CAREFULLLLL****************

pl=15;
lw=2;
lfs=14;

addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC
addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/QC/Nifti_Util

m_cnt = 1;
for m = ModeList
    s_cnt = 1;
    for s=SubList
            %load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/' Site{1} '/R/' Site{1} '_' s{1} '_fMRIDiag/' Site{1} '_' s{1} '_noglobal_' m{1} '_Council.mat'])
            load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/HCP/DSE/R/' m{1} '/HCP_DSE_' m{1} '_' s{1} '_' GSRStat{1} '.mat'],'DSE_Stat','V')
%------------------------------------------------------------------
            
            %Nu_X2_m3d3(s_cnt,m_cnt) = DVARS_Stat_X2_m3d3.nu;
            %Nu_X2_m1d3(s_cnt,m_cnt) = DVARS_Stat_X2_m1d3.nu;
            
            %Nu_X2_m3d1(s_cnt,m_cnt) = DVARS_Stat_X2_m3d1.nu;
            %Nu_X2_m1d1(s_cnt,m_cnt) = DVARS_Stat_X2_m1d1.nu;
            
            Vxl(s_cnt,m_cnt) = DSE_Stat.dim(1);
%------------------------------------------------------------------
            
            Info_Stat{s_cnt,m_cnt}=DSE_Stat;
            Info_V{s_cnt,m_cnt}=V;
%------------------------------------------------------------------
            
            s_cnt = s_cnt + 1;
    end
        m_cnt = m_cnt+1;
end

%------------------------------------------------------------------

XLabTicks_Labels={'Raw','Minimal','Full'};
Col=get(groot,'defaultAxesColorOrder');
Acol=Col(5,:); % Green
Dcol=Col(1,:); % Blue
Scol=Col(3,:); % Yellow
FDcol=Col(2,:); % Red (/orange!)
%------------------------------------------------------------------
%-------------------------WHOLE VARS-------------------------------
%------------------------------------------------------------------

DSE=cellfun(@(x) x.Prntg(2:4),Info_Stat,'UniformOutput', false);
D=cellfun(@(x) double(x(1)),DSE);
S=cellfun(@(x) double(x(2)),DSE);
E=cellfun(@(x) double(x(3)),DSE);

whole_fh=figure('position',[50,500,400,550]); 
hold on;
%------------------------------------------------------------------
hdl0=subplot(12,1,[1 9]); 
hold on; box on;
%title(['Var Comp - ' Site{1} ' ' PPLine{1} ' ' GSRStat{1}])
[dsh,Dxxlim]=ScatterBoxPlots(D,'subplot',hdl0,'PointSize',pl,'Color',repmat(Dcol,3,1),'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);%,'line');
[ssh,Sxxlim]=ScatterBoxPlots(S,'subplot',hdl0,'PointSize',pl,'Color',repmat(Scol,3,1),'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);%,'line');
set(gca,'Xtick',1:size(D,2),'xticklabel',[],'FontSize',lfs-1)
%line([0 size(DSE,2)+1],[50 50],'color',[.5 .5 .5],'linestyle','-.','linewidth',2)
line(xlim,[50 50],'color',[.5 .5 .5],'linestyle','-.','linewidth',lw)
ylabel('\% of A-var','fontsize',lfs,'interprete','latex')
legend([dsh{1} ssh{1}],{'D-var','S-var'},'fontsize',lfs-2)

% text(Sxxlim(3,1),S(3,1),'\leftarrow 101107')
% text(Sxxlim(3,2),S(3,2),'\leftarrow 101107')
% text(Sxxlim(3,3),S(3,3),'\leftarrow 101107')
% 
% text(Sxxlim(45,1),S(45,1),'\leftarrow 135932')
% text(Sxxlim(45,2),S(45,2),'\leftarrow 135932')
% text(Sxxlim(45,3),S(45,3),'\leftarrow 135932')
%--------%
% text(Sxxlim(60,1),S(60,1),'\leftarrow 151627')
% text(Sxxlim(60,2),S(60,2),'\leftarrow 151627')
%text(Sxxlim(27,3),S(27,3),'\leftarrow 151627')
% 
% text(Sxxlim(27,1),  S(27,1),    '\leftarrow 122620')
% text(Sxxlim(27,2),  S(27,2),    '\leftarrow 122620')
%text(Sxxlim(60,3),S(60,3),'\leftarrow 122620')
%--------%--------%--------%--------%--------%--------
% text(Dxxlim(27,1),  D(27,1),  '\leftarrow 122620')
% text(Dxxlim(27,2),  D(27,2),  '\leftarrow 122620')
%text(Dxxlim(27,3),D(27,3),'\leftarrow 122620')

% text(Dxxlim(60,1),D(60,1),'\leftarrow 151627')
% text(Dxxlim(60,2),D(60,2),'\leftarrow 151627')
%text(Dxxlim(60,3),D(60,3),'\leftarrow 151627')
%--------%
% text(Dxxlim(3,1),D(3,1),'\leftarrow 101107')
% text(Dxxlim(3,2),D(3,2),'\leftarrow 101107')
% text(Dxxlim(3,3),D(3,3),'\leftarrow 101107')

% text(Dxxlim(45,1),D(45,1),'\leftarrow 135932')
% text(Dxxlim(45,2),D(45,2),'\leftarrow 135932')
% text(Dxxlim(45,3),D(45,3),'\leftarrow 135932')

%------------------------------------------------------------------
% text(0.2,D(45,1),'135932')
% text(0.2,D(3,1),'101107')
% text(0.2,D(60,1),'151627')
% text(0.2,D(27,1),'122620')


%NB!! 
ss_list=find(ismember(SubList,{'101107','135932','151627','122620'}));
MarkerList={'d','>','o','s'};
%MarkerList={'o','o','o','o'};
ss_cnt=1;
for ss=ss_list
    SubList(ss)
    plot(Dxxlim(ss,:),D(ss,:),'Color',[FDcol 0.5],'linestyle','-.')
    plot(Sxxlim(ss,:),S(ss,:),'Color',[FDcol 0.5],'linestyle','-.')
    for mm=1:3
    scatter(Dxxlim(ss,mm),D(ss,mm),'MarkerFaceColor',Dcol,'Marker',MarkerList{ss_cnt},'linewidth',1,'markeredgecolor','r')
    scatter(Sxxlim(ss,mm),S(ss,mm),'MarkerFaceColor',Scol,'Marker',MarkerList{ss_cnt},'linewidth',1,'markeredgecolor','r')
    end
    ss_cnt=ss_cnt+1;
end


%------------------------------------------------------------------
hdl1=subplot(12,1,[10 11]); 
hold on; box on;
%%%HERE
esh=ScatterBoxPlots(log10(E),'subplot',hdl1,'PointSize',pl,'Color',repmat([.5 .5 .5],3,1),'line');
ylabel('log$_{10}$(\% of A-var)','fontsize',lfs-2,'interprete','latex')
%xlabel('Level of Pre-processing','fontsize',12)
set(gca,'Xticklabel',XLabTicks_Labels,'XTickLabelRotation',45,'Xtick',1:size(D,2),'FontSize',lfs-1)%set(gca,'Xtick',1:size(D,2),'xticklabel',[])%'FontSize',lfs)
legend([esh{1}],{'E-var'},'fontsize',lfs-2,'location','southwest')
%------------------------------------------------------------------
% hdl2=subplot(12,1,[11 12]);
% hold on; box on; 
% ScatterBoxPlots(Nu_X2_m3d3./Vxl,'subplot',hdl2,'Col',FDcol);
% 
% set(gca,'Xticklabel',XLabTicks_Labels,'XTickLabelRotation',45,'Xtick',1:size(Nu_X2_m3d3,2),'FontSize',lfs-1);
% 
 xlabel('Level of Pre-processing','fontsize',lfs,'interprete','latex')
% legend({'Spatial EDF'},'fontsize',lfs-2)
% ylabel('\% of DF','fontsize',lfs,'interprete','latex')
% %set(gca,'YTickLabel',strcat(cellstr(num2str(round(sqrt(yticks))')),'^2'),'FontSize',lfs-1);
% 
 set(gcf,'color','w');

export_fig(whole_fh,['R/Group_WholeVarComps_' Site{1} '_' GSRStat{1} '_AnotSM.pdf'])
%------------------------------------------------------------------
%-------------------------GLOBAL VARS------------------------------
%------------------------------------------------------------------
gDSE=cellfun(@(x) x.Prntg(6:8),Info_Stat,'UniformOutput', false);
gD=cellfun(@(x) double(x(1)),gDSE);
gS=cellfun(@(x) double(x(2)),gDSE);
gE=cellfun(@(x) double(x(3)),gDSE);

figure('position',[50,500,400,550]); hold on;
%------------------------------------------------------------------
hdl0=subplot(12,1,[1 9]); hold on; box on;
%title(['Global Var Comp - ' Site{1} ' ' PPLine{1}  ' ' GSRStat{1}])
[gdsh,gDxxlim]=ScatterBoxPlots(gD,'subplot',hdl0,'PointSize',pl,'Color',repmat(Dcol,3,1),'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
[gssh,gSxxlim]=ScatterBoxPlots(gS,'subplot',hdl0,'PointSize',pl,'Color',repmat(Scol,3,1),'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
set(gca,'Xtick',1:size(gD,2))
set(gca,'xticklabel',[],'FontSize',lfs-1)
ylabel('\% of A-var','fontsize',lfs,'interprete','latex')
legend([gdsh{1} gssh{1}],{'D_G-var','S_G-var'},'fontsize',lfs-2)

ss_cnt=1;
for ss=ss_list
    plot(gDxxlim(ss,:),gD(ss,:),'Color',[FDcol 0.5],'linestyle','-.')
    plot(gSxxlim(ss,:),gS(ss,:),'Color',[FDcol 0.5],'linestyle','-.')
    for mm=1:3
    scatter(gDxxlim(ss,mm),gD(ss,mm),'MarkerFaceColor',Dcol,'Marker',MarkerList{ss_cnt},'linewidth',1,'markeredgecolor','r')
    scatter(gSxxlim(ss,mm),gS(ss,mm),'MarkerFaceColor',Scol,'Marker',MarkerList{ss_cnt},'linewidth',1,'markeredgecolor','r')
    end
    ss_cnt=ss_cnt+1;
end

%------------------------------------------------------------------
hdl1=subplot(12,1,[10 11]); hold on; box on;
gesh=ScatterBoxPlots(log10(gE),'subplot',hdl1,'PointSize',pl,'Color',repmat([.5 .5 .5],3,1),'line');
xlabel('Level of Pre-processing','fontsize',lfs,'interprete','latex')
ylabel('log$_{10}$(\% of A-var)','fontsize',lfs-2,'interprete','latex')
set(gca,'Xticklabel',XLabTicks_Labels,'XTickLabelRotation',45,'Xtick',1:size(gD,2),'FontSize',lfs-1);
legend([gesh{1}],{'E_G-var'},'fontsize',lfs-2,'location','southwest')
set(gcf,'color','w');

export_fig(['R/Group_GlobalVarComps_' Site{1} '_' GSRStat{1} '_AnotSM.pdf'])
%------------------------------------------------------------------