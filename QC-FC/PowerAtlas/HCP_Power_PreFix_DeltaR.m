clear 
close all


%SubList={'100307','103414','105115','111312','113619','115320','117122','118730','123117','151526','187345','303624','132118','901139','171330','263436','191336','779370','195041','145127','172029'};
%SubList={'115320'};
%SubList={'105115'};
%SubList={'117122'};

SubList={'118730','115320'};

addpath /Users/sorooshafyouni/Home/GitClone/DVARS/Nifti_Util
%--------------------------------------------------------------

Col=get(groot,'defaultAxesColorOrder');
Acol=Col(5,:); % Green
%FDcol=Col(2,:); % Red (/orange!)
FDcol=[.5 .5 .5];
Dcol=Col(1,:); % Blue
Scol=Col(3,:); % Yellow
Ecol=Col(4,:); % Purple


nRlz=1000;
alevel=0.025;

idx=find(triu(ones(264),1));
cnt_s=1;
for s=SubList
    
    disp([num2str(cnt_s) ' ' s{1} '  ------'])
    load(['R/NetMats/HCP_MinimalF_' s{1} '_Power_NetMats_DVC.mat'],'ED','mat','mat_dvr','FDts','DVARS','DSE_Stat','mat_dvs_scrmbld')

    sarr=[]; %BinCountsX_tmp=[]; BinEdgeX_tmp=[];
    for i=1:nRlz
        scmat_tmp=mat-mat_dvs_scrmbld(:,:,i);
        sarr=[sarr scmat_tmp(idx)]; 
        clear *_tmp
    end
    
    ED=2*ED(idx);
    [sED,ii] = sort(ED);
    %105115 is the best candidate!
    figure('position',[50,500,700,300]); hold on; box on; 
    
    subplot(1,2,1); hold on; box on; grid on;
    title([s{1}])
    %ScatterPlots---------------------------
    dmat_dv=mat-mat_dvr;
    dmat_dv=dmat_dv(idx);
    scatter(sED,dmat_dv(ii),'marker','.','markerfacecolor',Dcol,'markeredgecolor',Dcol,'MarkerFaceAlpha',0.1);

    scmat_tmp=sarr(:,randi(1000));
    scatter(sED,scmat_tmp(ii),'marker','.','markerfacecolor',FDcol,'markeredgecolor',FDcol,'MarkerFaceAlpha',0.1)
    
    %SmoothedLine---------------------------
    YY1      = smooth(sED,dmat_dv(ii),0.1,'lowess');
    h1       = line(sED,YY1,'linewidth',1.5,'color','r');
    
    smat_tmp = mean(sarr,2);
    YY0      = smooth(sED,smat_tmp(ii),0.1,'lowess');
    h0       = line(sED,YY0,'linewidth',1.5,'color','g');
    
    xlim([0 max(ED)+10])
    ylabel('$\Delta$ r','fontsize',12,'interpreter','Latex')
    xlabel('Inter-node Distance (mm)','fontsize',12,'interpreter','Latex')
    
    
    subplot(1,2,2); hold on; axis off
    title([s{1}])
    %find the affected edges
    H=prctile(abs(sarr),(1-alevel./numel(dmat_dv))*100,2)<abs(dmat_dv);
    sig_idx=idx(H);
    
    sig_mat=zeros(264);
    sig_mat(sig_idx)=1;
    
    spy(sig_mat)
    
    set(gcf,'Color','w')
    export_fig(hf0,[s{1} '_Power_PreFix_DeltaR.pdf']);
    
    cnt_s=cnt_s+1;
end