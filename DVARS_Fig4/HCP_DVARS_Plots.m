% This script produce Figure 4.
%
% NB! the two subjects which were used on the table are 8th and 6th element
% of the SubList below. 
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

clear

%Addpath /AuxDraw available on the same repository
addpath ~/DVARS_Paper17/AuxDraw
addpath ~/DVARS_Paper17/AuxDraw/concentricplots

Site={'HCP'};
SubList={'100307','103414','105115','111312','113619','115320','117122','118730','123117','151526','187345','303624','132118','901139','171330','263436','191336','779370','195041','145127','172029'};


%We used subject 115320 and 118730 (6th and 8th elements, repectively).
% Therefore, you have to choose between either 6 or 8. However, the plot is
% reproducible for all other subjects as long as you produce the
% HCP_Council files for them.
s=SubList{6}; %8 and 6

%Choose which stage of the pre-processing you are intrested in. Choose
%between 'Pre_Fix' and 'Post_Fix'. NB that by default the last sub-plot of
%the figure is the step further of the pre-processing under study! In other
%words, if you choose Pre_Fix, the last sub-plot is results of the
%Post_Fix.
m={'Pre_Fix'}; %Post_Fix Pre_Fix

% You can see the results of different testing techniques. 
% X2_m3d3: Chi-squared stat, empirical median, 1/3 power transfo on var
% X2_m3d1: Chi-squared stat, empirical median, no power transfo on var
% X2_m1d3: Chi-squared stat, mean as \mu^D_0, 1/3 power transfo on var
% Z      : Use the good old Z stats!
TestMeth={'X2_m3d3'}; %'X2_m3d1' 'X2_m1d1' 'Z'

%If you want to save the figures, set saveflag to 1. Only do this when you
%have export_fig already in your Matlab path 
saveflag=0;

%#####################################################

load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/' Site{1} '/R/' Site{1} '_' s '_fMRIDiag/' Site{1} '_' s '_noglobal_' m{1} '_Council_Test.mat'])
T=DSE_Stat.dim(2);
Time=1:T;
hTime=(1:(T-1))+0.5;
eTime=[1 T];

% pvals  = eval(['DVARS_Stat_' TestMeth{1} '.pvals']);
% Idx=find(pvals<(0.05./(T0-1)));
% Idx=Idx+1;

if    contains(TestMeth{1},'X2')
    NDVARS = eval(['DVARS_Stat_' TestMeth{1} '.NDVARS_X20']);
elseif  contains(TestMeth{1},'Z')
    NDVARS = eval(['DVARS_Stat_' TestMeth{1} '.NDVARS_Z']);
else
    error('Unknown test method: choose btwn: X2d3, X2d1, Z');
end

pvals  = eval(['DVARS_Stat_' TestMeth{1} '.pvals']);
Idx=find(pvals<(0.05./(T-1)));

% Idx=find(DVARS_Stat_Z.pvals<0.05/(T-1));
% NDVARS=DVARS_Stat_Z.NDVARS_Z;


% Idx=find(DVARS_Stat_X2_m1d3.pvals<0.05/(T-1));
% NDVARS=DVARS_Stat_X2_m1d3.NDVARS_X2;

%Idx=find(fdr_bh(DVARS_Stat_X2d3.pvals));

diffpNDVARS=(V.Dvar_ts-median(V.Dvar_ts))/mean(V.Avar_ts)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alevel=0.05/T;
z2p = @(p) -sqrt(2) * erfcinv(p*2);
zl=z2p(alevel);

lw=1;
Col=get(groot,'defaultAxesColorOrder');
Acol=Col(5,:); % Green
FDcol=Col(2,:); % Red (/orange!)
Dcol=Col(1,:); % Blue
Scol=Col(3,:); % Yellow
Ecol=Col(4,:); % Purple

lfs=12; %label fontsize
lfts=14;

nsp=20;

spOrd={[1 2],[3 4],[5 6],[7 8],[9 10],[11 12],[13 16],[17 nsp]};

ffh1=figure('position',[50,500,700,1400]);
hold on; 
aa=suptitle([Site{1} '-' num2str(s) '-' m{1} '_' TestMeth{1}]);
aa.Interpreter='none';
sph0=subplot(nsp,1,spOrd{5});%--------------------------------------------
hold on; box on; grid on; axis tight;
%    plot(hTime,FDts,'Color',FDcol,'linewidth',lw-0.5)
%    cvl0=line([0 T-1],[.2 .2],'Color','k','linewidth',lw+.5,'linestyle','-.'); 
%    cvl1=line([0 T-1],[.5 .5],'Color','k','linewidth',lw+.5,'linestyle','--');
%    legend([cvl0 cvl1],{'Strict FD threshold','Lenient FD threshold'})
%    ylabel('FD (mm)','fontsize',lfts,'interpreter','latex')

    plot(hTime,DVARS_Stat_Z.NDVARS,'Color',Dcol,'linewidth',lw)
    ylabel('RDVARS','fontsize',lfts-2,'interpreter','latex')
    set(sph0,'ycolor','k','XTick',[],'xlim',[1 T])
    
    PatchMeUp(Idx);
sph1=subplot(nsp,1,spOrd{1});%--------------------------------------------
hold on; box on; grid on; axis tight; 
    %yyaxis right
    plot(hTime,DVARS_Z,'Color',Dcol,'linewidth',lw)
    set(sph1,'ycolor','k','XTick',[],'xlim',[1 T])
    ylabel('DVARS','fontsize',lfts-2,'interpreter','latex')
    
    %yyaxis left
    %set(sph1,'ycolor','k','YTick',[])
    
    PatchMeUp(Idx);
sph2=subplot(nsp,1,spOrd{2});%--------------------------------------------
hold on; box on; grid on; axis tight; 
    plot(hTime,sqrt(V.Dvar_ts),'Color',Dcol,'linewidth',lw)
    %ylabel('$\sqrt{\mathrm{D}-\mathrm{var}}$','fontsize',lfts,'interpreter','latex')
    %ylabel('$\surd$D-var','fontsize',lfts,'interpreter','latex')
    ylabel('$\sqrt{\mathrm{D\textendash var}}$','fontsize',lfts,'interpreter','latex')
    
    
    
    set(sph2,'ycolor','k','XTick',[],'xlim',[1 T]) 
    
    PatchMeUp(Idx);
sph3=subplot(nsp,1,spOrd{3});%--------------------------------------------
hold on; box on; grid on; axis tight; 
    %yyaxis right
    plot(V.Dvar_ts./mean(V.Avar_ts)*100,'Color',Dcol,'linewidth',lw)
    ylabel('$\%$ D\textendash var','fontsize',lfts,'interpreter','latex');
    set(sph3,'ycolor','k','XTick',[],'xlim',[1 T])
    
    %yyaxis left
    %set(sph3,'ycolor','k','YTick',[]) 
    
    PatchMeUp(Idx);
sph4=subplot(nsp,1,spOrd{4});%--------------------------------------------
hold on; box on; grid on; axis tight;   
    plot(hTime,diffpNDVARS,'Color',Dcol,'linewidth',lw)
    set(sph4,'ycolor','k','XTick',[],'xlim',[1 T])
    ylabel('$\Delta \%$ D\textendash var','fontsize',lfts-1,'interpreter','latex')
    
    PatchMeUp(Idx);
sph5=subplot(nsp,1,spOrd{6});%--------------------------------------------
hold on; box on; grid on; axis tight;
    %yyaxis right
    plot(hTime,NDVARS,'Color',Dcol,'linewidth',lw)
    set(sph5,'ycolor','k','XTick',[],'xlim',[1 T])
    cvl=line([0 T-1],[zl zl]*-1,'Color','r','linewidth',lw,'linestyle','-.');
    legend([cvl],{'Adj critical value'},'location','southeast')
    ylabel('Z(D\textendash var)','fontsize',lfts,'interpreter','latex')
    
    %yyaxis left
    %set(sph5,'ycolor','k','YTick',[])
    PatchMeUp(Idx);
sph6=subplot(nsp,1,spOrd{7});%--------------------------------------------
hold on; box on;
yyaxis(sph6,'left')
    line(Time,sqrt(V.Avar_ts),'LineStyle','-','linewidth',lw,'color',Acol)
    line(Time,ones(1,T).*mean(sqrt(V.Avar_ts)),'LineStyle',':','linewidth',.5,'color',Acol)
    
    line(hTime,sqrt(V.Dvar_ts),'LineStyle','-','linewidth',lw,'color',Dcol)
    line(hTime,ones(1,T-1).*mean(sqrt(V.Dvar_ts)),'LineStyle',':','linewidth',.5,'color',Dcol)
    
    line(hTime,sqrt(V.Svar_ts),'LineStyle','-','linewidth',lw,'color',Scol)
    line(hTime,ones(1,T-1).*mean(sqrt(V.Svar_ts)),'LineStyle',':','linewidth',.5,'color',Scol)
    
    line(eTime,sqrt(V.Evar_ts),'LineStyle','none','Marker','o','markerfacecolor',Ecol,'linewidth',3,'color',Ecol)
    ylabel('$\sqrt{\mathrm{Variance}}$','fontsize',lfts,'interpreter','latex')  
    
    YLim2=ylim.^2/mean(V.Avar_ts)*100;
    
    %YLim2=sqrt((ylim/mean(V.Avar_ts)).^2*100);
    
    set(sph6,'ycolor','k','xlim',[1 T])
    
yyaxis(sph6,'right')
    YTick2=PrettyTicks_fnc(YLim2,1/4); 
    %YTick2 = PrettyTicksBetter(YLim2.^2);
    YTick=sqrt(YTick2);
    set(sph6,'Ylim',sqrt(YLim2),'YTick',sqrt(YTick2),'YtickLabel',num2str([YTick2']));
    ylabel('$\%$ of A \textendash var','fontsize',lfts,'interpreter','latex')
    h=abline('h',YTick);
    set(h,'linestyle','-','color',[.5 .5 .5]); %the grids!
    set(sph6,'ycolor','k','xlim',[1 T],'XTick',[])    
    
    PatchMeUp(Idx);
    
sph7=subplot(nsp,1,spOrd{8});%--------------------------------------------

PostFix=load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/' Site{1} '/R/' Site{1} '_' s '_fMRIDiag/' Site{1} '_' s '_noglobal_Post_Fix_Council.mat']);

PostFix_pvals  = eval(['PostFix.DVARS_Stat_' TestMeth{1} '.pvals']);
PostFix_Idx=find(PostFix_pvals<(0.05./(T-1)));

hold on; box on;
yyaxis(sph7,'left')
    line(Time,sqrt(PostFix.V.Avar_ts),'LineStyle','-','linewidth',lw,'color',Acol)
    line(Time,ones(1,T).*mean(sqrt(PostFix.V.Avar_ts)),'LineStyle',':','linewidth',.5,'color',Acol)
    
    line(hTime,sqrt(PostFix.V.Dvar_ts),'LineStyle','-','linewidth',lw,'color',Dcol)
    line(hTime,ones(1,T-1).*mean(sqrt(PostFix.V.Dvar_ts)),'LineStyle',':','linewidth',.5,'color',Dcol)
    
    line(hTime,sqrt(PostFix.V.Svar_ts),'LineStyle','-','linewidth',lw,'color',Scol)
    line(hTime,ones(1,T-1).*mean(sqrt(PostFix.V.Svar_ts)),'LineStyle',':','linewidth',.5,'color',Scol)
    
    line(eTime,sqrt(PostFix.V.Evar_ts),'LineStyle','none','Marker','o','markerfacecolor',Ecol,'linewidth',3,'color',Ecol)
    ylabel('$\sqrt{\mathrm{Variance}}$','fontsize',lfts,'interpreter','latex')  
    
    PostFix_YLim2=ylim.^2/mean(PostFix.V.Avar_ts)*100;
    
    %YLim2=sqrt((ylim/mean(V.Avar_ts)).^2*100);
    
    set(sph7,'ycolor','k','xlim',[1 T])
    
yyaxis(sph7,'right')
    PostFix_YTick2=PrettyTicks_fnc(PostFix_YLim2,1/2); 
    %YTick2 = PrettyTicksBetter(YLim2.^2);
    PostFix_YTick=sqrt(PostFix_YTick2);
    set(sph7,'Ylim',sqrt(PostFix_YLim2),'YTick',sqrt(PostFix_YTick2),'YtickLabel',num2str([PostFix_YTick2']));
    ylabel('$\%$ of A \textendash var','fontsize',lfts,'interpreter','latex')
    h=abline('h',PostFix_YTick);
    set(h,'linestyle','-','color',[.5 .5 .5]); %the grids!
    set(sph7,'ycolor','k','xlim',[1 T])    
    
    PatchMeUp(PostFix_Idx);
    
%--------------------------    
if saveflag; export_fig(ffh1,['Figs/VARS_' Site{1} '_' s '_' m{1} '_' TestMeth{1} '.pdf']); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t3_varn = {'Avar','Dvar','Svar','Evar'};
t3_rown = {'Whole','Global','non-Global'};
Row_labs = {'Avar','Dvar','Svar','Evar','A_{G}-var','D_{G}-var','S_{G}-var','E_{G}-var'};
Col_labs = {'RMS','Prcntg_of_whole','Rel_iid'};
    
Var_Tab = [V.w_Avar,V.w_Dvar,V.w_Svar,V.w_Evar;...
V.g_Avar,V.g_Dvar,V.g_Svar,V.g_Evar;...
V.ng_Avar,V.ng_Dvar,V.ng_Svar,V.ng_Evar];

disp([Site{1} '-' num2str(s) '-' m{1}])

disp('----------------------')
%disp('Sum-of-Mean-Squared (SMS) Table')
%disp(array2table(fix(Var_Tab),'VariableNames',t3_varn,'RowNames',t3_rown))
%disp('------------')
T2print=array2table([DSE_Stat.RMS',DSE_Stat.Prntg',DSE_Stat.RelVar'],'VariableNames',Col_labs,'RowNames',Row_labs);
disp(T2print)
writetable(T2print,['Figs/DVARS_AND_Table_' Site{1} '_' s '_' m{1} '.csv'],'WriteRowNames',1)
disp('----------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ffh2=figure('position',[50,500,700,400],'color','w'); 
hold on; 
%aa=suptitle([Site{1} '-' num2str(s) '-' m{1}]);
%aa.Interpreter='none';
%spa3=subplot(1,2,1);
hold on; box on; axis off;
title('Whole Var \%A \textendash var','interpreter','latex','fontsize',lfts)
concentricplots_draw(0,0,DSE_Stat,'whole','verbose',0,'figure',ffh2,'Col',[Dcol;Scol;Ecol]);

if saveflag; export_fig(ffh2,['Figs/ConCent_Whole_' Site{1} '_' s '_' m{1} '.pdf']); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ffh3=figure('position',[50,500,200,200],'color','w'); 
hold on; grid on; box on; 

%aa=title([Site{1} '-' num2str(s) '-' m{1}]);
%aa.Interpreter='none';
%spa0=subplot(1,2,2);
%hold on; grid on; box on; 

title('Whole','interpreter','latex','fontsize',lfts)
bp0=bar(1:4,diag(DSE_Stat.RelVar(1:4)),'stacked');
icnt=1;
for icol=[5 1 3 4]
    bp0(icnt).FaceColor=Col(icol,:);
    icnt=icnt+1;
end
line([0 5],[1 1],'color','r','linewidth',lw)
set(gca,'xlim',[0 5],'XTick',1:4,'XTickLabel',Row_labs(1:4),'fontsize',lfts,'XTickLabelRotation',45)
ylabel('Relative to IID','interpreter','latex','fontsize',lfts)

%export_fig(['DVARS_AND_Table_' Site{1} '_' s '_' m{1} '.pdf'])
if saveflag; export_fig(ffh3,['Figs/RMSRel_Whole_' Site{1} '_' s '_' m{1} '.pdf']); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ffh4=figure('position',[50,500,700,400],'color','w'); 
hold on; box on; axis off;

%aa=suptitle([Site{1} '-' num2str(s) '-' m{1}]);
%aa.Interpreter='none';
%spa4=subplot(1,2,1);
%hold on; box on; axis off;
title('Global Var \%A$_{Gt}$','interpreter','latex','fontsize',lfts)
concentricplots_draw(0,0,DSE_Stat,'global','verbose',0,'figure',ffh4,'Col',[Dcol;Scol;Ecol]);

if saveflag; export_fig(ffh4,['Figs/ConCent_Global_' Site{1} '_' s '_' m{1} '.pdf']); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ffh5=figure('position',[50,500,200,200],'color','w'); 
%spa1=subplot(1,2,2);
hold on; grid on; box on; 
title('Global','interpreter','latex','fontsize',lfts)
bp1=bar(1:4,diag(log10(DSE_Stat.RelVar(5:end))),'stacked');
icnt=1;
for icol=[5 1 3 4]
    bp1(icnt).FaceColor=Col(icol,:);
    icnt=icnt+1;
end
line([0 5],[1 1],'color','r','linewidth',lw)
set(gca,'xlim',[0 5],'XTick',1:4,'XTickLabel',Row_labs(5:end),'fontsize',lfts,'XTickLabelRotation',45)
ylabel('log$_{10}$(Relative to IID)','interpreter','latex','fontsize',lfts)

if saveflag; export_fig(ffh5,['Figs/RMSRel_Global__' Site{1} '_' s '_' m{1} '.pdf']); end;
%set(ffh2,'color','w')

