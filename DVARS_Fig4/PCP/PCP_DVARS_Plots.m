clear

close all

Site={'NYU'};
SubList={'0051036','0051038','0051039','0051040','0051041','0051042','0051044','0051045','0051046'...
    ,'0051047','0051048','0051049','0051050','0051051','0051052','0051053','0051054','0051055',...
    '0051056','0051057','0051058','0051059','0051060','0051061','0051062'};

s=SubList{18}; %13 and 18

m={'func_minimal'}; %nofilt_noglobal -- func_minimal

TestMeth={'X2_m3d3'}; %'X2_m3d1' 'X2_m1d1' 'Z'

Who2Believe='PCPWebsite'; %ppppffff

PracThr=5;

saveflag=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/PCP/R/' Site{1} '_' s '_fMRIDiag/' Site{1} '_' s '_noglobal_' m{1} '_Council_' Who2Believe '.mat'])
T=DSE_Stat.dim(2);
Time=1:T;
hTime=(1:(T-1))+0.5;
eTime=[1 T];

% pvals  = eval(['DVARS_Stat_' TestMeth{1} '.pvals']);
% Idx=find(pvals<(0.05./(T0-1)));
% Idx=Idx+1;

if    contains(TestMeth{1},'X2')
    NDVARS = eval(['DVARS_Stat_' TestMeth{1} '.SDVARS_X2']);
elseif  contains(TestMeth{1},'Z')
    NDVARS = eval(['DVARS_Stat_' TestMeth{1} '.SDVARS_Z']);
else
    error('Unknown test method: choose btwn: X2d3, X2d1, Z');
end

pvals  = eval(['DVARS_Stat_' TestMeth{1} '.pvals']);
Idxs=find(pvals<(0.05./(T-1)));

DpDvar  = eval(['DVARS_Stat_' TestMeth{1} '.DeltapDvar']);
Idxp=find(pvals<(0.05./(T-1)) & DpDvar>PracThr);

% Idx=find(DVARS_Stat_Z.pvals<0.05/(T-1));
% NDVARS=DVARS_Stat_Z.NDVARS_Z;


% Idx=find(DVARS_Stat_X2_m1d3.pvals<0.05/(T-1));
% NDVARS=DVARS_Stat_X2_m1d3.NDVARS_X2;

%Idx=find(fdr_bh(DVARS_Stat_X2d3.pvals));

diffpNDVARS=(V.Dvar_ts-median(V.Dvar_ts))/mean(V.Avar_ts)*100;


addpath('/Users/sorooshafyouni/Home/GitClone/DVARS_Paper17/AuxDraw/')
%addpath /Users/sorooshafyouni/Home/DVARS/fMRIDiag/HCP/AuxDraw
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

lfs  = 12; %label fontsize
lfts = 14;
nsp  = 20;

trans = 0.3;
IdxpCol     = [FDcol];
IdxsCol     = [0 0 0];
PtchWdth    = 2;

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

    plot(hTime,DVARS_Stat_Z.RDVARS,'Color',Dcol,'linewidth',lw)
    ylabel('RDVARS','fontsize',lfts-2,'interpreter','latex')
    set(sph0,'ycolor','k','XTick',[],'xlim',[1 T])
    
    %PatchMeUp(Idx);
    PatchMeUpLine(setdiff(Idxs,Idxp)+1,Idxp+1,PtchWdth,IdxsCol,IdxpCol,trans,sph0);

sph1=subplot(nsp,1,spOrd{1});%--------------------------------------------
hold on; box on; grid on; axis tight; 
    %yyaxis right
    plot(hTime,DVARS_Z,'Color',Dcol,'linewidth',lw)
    set(sph1,'ycolor','k','XTick',[],'xlim',[1 T])
    ylabel('DVARS','fontsize',lfts-2,'interpreter','latex')
    
    %yyaxis left
    %set(sph1,'ycolor','k','YTick',[])
    
    %PatchMeUp(Idx);
    PatchMeUpLine(setdiff(Idxs,Idxp)+1,Idxp+1,PtchWdth,IdxsCol,IdxpCol,trans,sph1);

sph2=subplot(nsp,1,spOrd{2});%--------------------------------------------
hold on; box on; grid on; axis tight; 
    plot(hTime,sqrt(V.Dvar_ts),'Color',Dcol,'linewidth',lw)
    %ylabel('$\sqrt{\mathrm{D}-\mathrm{var}}$','fontsize',lfts,'interpreter','latex')
    %ylabel('$\surd$D-var','fontsize',lfts,'interpreter','latex')
    ylabel('$\sqrt{\mathrm{D\textendash var}}$','fontsize',lfts,'interpreter','latex')
    
    
    set(sph2,'ycolor','k','XTick',[],'xlim',[1 T]) 
    
    %PatchMeUp(Idx);
    PatchMeUpLine(setdiff(Idxs,Idxp)+1,Idxp+1,PtchWdth,IdxsCol,IdxpCol,trans,sph2);
sph3=subplot(nsp,1,spOrd{3});%--------------------------------------------
hold on; box on; grid on; axis tight; 
    %yyaxis right
    plot(V.Dvar_ts./mean(V.Avar_ts)*100,'Color',Dcol,'linewidth',lw)
    ylabel('$\%$ D\textendash var','fontsize',lfts,'interpreter','latex');
    set(sph3,'ycolor','k','XTick',[],'xlim',[1 T])
    
    %yyaxis left
    %set(sph3,'ycolor','k','YTick',[]) 
    
    %PatchMeUp(Idx);
    PatchMeUpLine(setdiff(Idxs,Idxp)+1,Idxp+1,PtchWdth,IdxsCol,IdxpCol,trans,sph3);

sph4=subplot(nsp,1,spOrd{4});%--------------------------------------------
hold on; box on; grid on; axis tight;   
    plot(hTime,diffpNDVARS,'Color',Dcol,'linewidth',lw)
    set(sph4,'ycolor','k','XTick',[],'xlim',[1 T])
    ylabel('$\Delta \%$ D\textendash var','fontsize',lfts-1,'interpreter','latex')
    
    %PatchMeUp(Idx);
    PatchMeUpLine(setdiff(Idxs,Idxp)+1,Idxp+1,PtchWdth,IdxsCol,IdxpCol,trans,sph4);

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
    
    %PatchMeUp(Idx);
    PatchMeUpLine(setdiff(Idxs,Idxp)+1,Idxp+1,PtchWdth,IdxsCol,IdxpCol,trans,sph5);

sph6=subplot(nsp,1,spOrd{7});%--------------------------------------------
hold on; box on;  axis tight
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
    %YTick2=PrettyTicks_fnc(YLim2,1/4); 
    YTick2=PrettyTicks0(YLim2); 
    %YTick2 = PrettyTicksBetter(YLim2.^2);
    YTick=sqrt(YTick2);
    set(sph6,'Ylim',sqrt(YLim2),'YTick',sqrt(YTick2),'YtickLabel',num2str([YTick2']));
    ylabel('$\%$ of A \textendash var','fontsize',lfts,'interpreter','latex')
    h=abline('h',YTick);
    set(h,'linestyle','-','color',[.5 .5 .5]); %the grids!
    set(sph6,'ycolor','k','xlim',[1 T],'XTick',[])    
    
    %PatchMeUp(Idx);
    PatchMeUpLine(setdiff(Idxs,Idxp)+1,Idxp+1,PtchWdth,IdxsCol,IdxpCol,trans,sph6);

sph7=subplot(nsp,1,spOrd{8});%--------------------------------------------

mTemp={'nofilt_noglobal'};
PostFix=load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/PCP/R/' Site{1} '_' s '_fMRIDiag/' Site{1} '_' s '_noglobal_' mTemp{1} '_Council_' Who2Believe '.mat']);

PostFix_pvals  = eval(['PostFix.DVARS_Stat_' TestMeth{1} '.pvals']);
PostFix_Idx=find(PostFix_pvals<(0.05./(T-1)));

PostFix_DpDvar  = eval(['PostFix.DVARS_Stat_' TestMeth{1} '.DeltapDvar']);
PostFix_Idxp    = find(PostFix_pvals<(0.05./(T-1)) & PostFix_DpDvar>PracThr);


hold on; box on;  axis tight
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
    %PostFix_YTick2=PrettyTicks_fnc(PostFix_YLim2,1/4); 
    PostFix_YTick2=PrettyTicks0(PostFix_YLim2);
    %YTick2 = PrettyTicksBetter(YLim2.^2);
    PostFix_YTick=sqrt(PostFix_YTick2);
    set(sph7,'Ylim',sqrt(PostFix_YLim2),'YTick',sqrt(PostFix_YTick2),'YtickLabel',num2str([PostFix_YTick2']));
    ylabel('$\%$ of A \textendash var','fontsize',lfts,'interpreter','latex')
    h=abline('h',PostFix_YTick);
    set(h,'linestyle','-','color',[.5 .5 .5]); %the grids!
    set(sph7,'ycolor','k','xlim',[1 T])    
    
    %PatchMeUp(PostFix_Idx);
    PatchMeUpLine(setdiff(PostFix_Idx,PostFix_Idxp)+1,PostFix_Idxp+1,PtchWdth,IdxsCol,IdxpCol,trans,sph7);

%--------------------------    
if saveflag; export_fig(ffh1,['Figs/VARS_' Site{1} '_' s '_' m{1} '_' TestMeth{1} '.pdf']); end;

clear h sph*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
addpath /Users/sorooshafyouni/Home/DVARS/concentricplots
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