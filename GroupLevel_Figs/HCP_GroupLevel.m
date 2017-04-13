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

clear; 
%close all;

Site={'HCP'};

PPLine={''};

%ModeList={'Unproc','Pre_Fix','Post_Fix'};
ModeList={'Unproc','Pre_Fix','Post_Fix'};

%21 Healthy
SubList={'100307','103414','105115','111312','113619','115320','117122','118730','123117','151526','187345','303624','132118','901139','171330','263436','191336','779370','195041','145127','172029'};

GSRStat={'noglobal'}; %%BE CAREFULLLLL****************

% You need to clone the DVARS directory:
% https://github.com/asoroosh/DVARS.git
% and addpath the directory (+ subfolders)

%addpath ~/DVARS/
%addpath ~/DVARS/Nifti_Util

lw=2;
lfs=14;

m_cnt = 1;
for m = ModeList
    s_cnt = 1;
    for s=SubList
            load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/' Site{1} '/R/' Site{1} '_' s{1} '_fMRIDiag/' Site{1} '_' s{1} '_noglobal_' m{1} '_Council.mat'])
            
            Nu_X2_m3d3(s_cnt,m_cnt) = DVARS_Stat_X2_m3d3.nu;
            Nu_X2_m1d3(s_cnt,m_cnt) = DVARS_Stat_X2_m1d3.nu;
            
            Nu_X2_m3d1(s_cnt,m_cnt) = DVARS_Stat_X2_m3d1.nu;
            Nu_X2_m1d1(s_cnt,m_cnt) = DVARS_Stat_X2_m1d1.nu;
            
            Info_Stat{s_cnt,m_cnt}=DSE_Stat;
            Info_V{s_cnt,m_cnt}=V;
            
            Vxl(s_cnt,m_cnt) = DSE_Stat.dim(1);
            
            s_cnt = s_cnt + 1;
    end
        m_cnt = m_cnt+1;
end

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

whole_fh=figure('position',[50,500,250,550]); 
hold on;
%------------------------------------------------------------------
hdl0=subplot(12,1,[1 9]); 
hold on; box on;
%title(['Var Comp - ' Site{1} ' ' PPLine{1} ' ' GSRStat{1}])
dsh=ScatterBoxPlots(D,'subplot',hdl0,'PointSize',20,'Col',Dcol);
ssh=ScatterBoxPlots(S,'subplot',hdl0,'PointSize',20,'Col',Scol);
set(gca,'Xtick',1:size(D,2),'xticklabel',[],'FontSize',lfs-1)
%line([0 size(DSE,2)+1],[50 50],'color',[.5 .5 .5],'linestyle','-.','linewidth',2)
line(xlim,[50 50],'color',[.5 .5 .5],'linestyle','-.','linewidth',lw)
ylabel('\% of A-var','fontsize',lfs,'interprete','latex')
legend([dsh{1} ssh{1}],{'D-var','S-var'},'fontsize',lfs-2)
%------------------------------------------------------------------
hdl1=subplot(12,1,[10 11]); 
hold on; box on;
%%%HERE
esh=ScatterBoxPlots(log10(E),'subplot',hdl1,'PointSize',20,'Col',[.5 .5 .5]);
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

%export_fig(whole_fh,['R/Group_WholeVarComps_' Site{1} '_' GSRStat{1} '.pdf'])
%------------------------------------------------------------------
%-------------------------GLOBAL VARS------------------------------
%------------------------------------------------------------------
gDSE=cellfun(@(x) x.Prntg(6:8),Info_Stat,'UniformOutput', false);
gD=cellfun(@(x) double(x(1)),gDSE);
gS=cellfun(@(x) double(x(2)),gDSE);
gE=cellfun(@(x) double(x(3)),gDSE);

figure('position',[50,500,250,550]); hold on;
%------------------------------------------------------------------
hdl0=subplot(12,1,[1 9]); hold on; box on;
%title(['Global Var Comp - ' Site{1} ' ' PPLine{1}  ' ' GSRStat{1}])
gdsh=ScatterBoxPlots(gD,'subplot',hdl0,'PointSize',20,'Col',Dcol);
gssh=ScatterBoxPlots(gS,'subplot',hdl0,'PointSize',20,'Col',Scol);
set(gca,'Xtick',1:size(gD,2))
set(gca,'xticklabel',[],'FontSize',lfs-1)
ylabel('\% of A-var','fontsize',lfs,'interprete','latex')
legend([gdsh{1} gssh{1}],{'D_G-var','S_G-var'},'fontsize',lfs-2)
%------------------------------------------------------------------
hdl1=subplot(12,1,[10 11]); hold on; box on;
gesh=ScatterBoxPlots(log10(gE),'subplot',hdl1,'PointSize',20,'Col',[.5 .5 .5]);
xlabel('Level of Pre-processing','fontsize',lfs,'interprete','latex')
ylabel('log$_{10}$(\% of A-var)','fontsize',lfs-2,'interprete','latex')
set(gca,'Xticklabel',XLabTicks_Labels,'XTickLabelRotation',45,'Xtick',1:size(gD,2),'FontSize',lfs-1);
legend([gesh{1}],{'E_G-var'},'fontsize',lfs-2,'location','southwest')
set(gcf,'color','w');

%export_fig(['R/Group_GlobalVarComps_' Site{1} '_' GSRStat{1} '.pdf'])
%------------------------------------------------------------------