clear;

close all

Site={'HCP'};

AtlasList={'Power2011'};

MarkerAlp=0.1;

addpath /Users/sorooshafyouni/Home/GitClone/DVARS/Nifti_Util

%21 Healthy 
SubList=[100307,103414,105115,111312,113619,115320,117122,118730,123117,151526,187345,303624,132118,901139,171330,263436,191336,779370,195041,145127,172029];
GSRStatList={'NoGSR'}; %%BE CAREFULLLLL****************

StateList={'DVCleaned'};%'Minimal','FDCleaned'};% raw DVCleaned

Col=get(groot,'defaultAxesColorOrder');
Acol=Col(5,:); % Green
%Dcol=Col(1,:); % Blue
Dcol=[0.5 0.5 0.5];
Scol=Col(3,:); % Yellow
FDCol=[1 0 0];

smthspan=0.1;

lalp=0.5;

for st=StateList
    for a=AtlasList
        for g=GSRStatList
                DV=[]; FD_con=[];
                cnt_s=1;
                for s=SubList

                    load(['R/NetMats/HCP_MinimalF_' num2str(s) '_Power_NetMats_DVC.mat'])
                    nn  = size(mat,1);
                    fge=nn*(nn-1)./2;
                    idx=find(triu(ones(nn),1));
                    
                    if isempty(rgs_key)
                        disp([num2str(s)  ' exluded.']);
                        continue;
                    end
%---------------------------------------------              
                        liemat(:,cnt_s)=mat(idx);
                        %Porpottion of scrubbed-------
                        FD_Null(cnt_s)=numel(mat_fd_lib);                    
%---------------------------------------------              
                        liemat_fd_lib(:,cnt_s)=mat_fd_lib(idx);
                        %Porpottion of scrubbed-------
                        FD_lib(cnt_s)=numel(fd_rmv_dp_lib);
%---------------------------------------------              
                        liemat_fd_con(:,cnt_s)=mat_fd_con(idx);
                        %Porpottion of scrubbed-------
                        FD_con(cnt_s)=numel(fd_rmv_dp_con);
%---------------------------------------------                        
                        liemat_dv(:,cnt_s)=mat_dvs(idx);
                        %Porpottion of scrubbed-------
                        DV(cnt_s)=numel(dv_rmv);
%---------------------------------------------                        
                        liemat_dvQ(:,cnt_s)=mat_dv_Qthr(idx);
                        %Porpottion of scrubbed-------
                        DV_thrQ(cnt_s)=numel(dvQthr_rmv);
%---------------------------------------------                        
                        liemat_dv_thr(:,cnt_s)=mat_dv_thr(idx);
                        %Porpottion of scrubbed-------
                        DV_thr(cnt_s)=numel(dvthr_rmv);                          
%---------------------------------------------
                    cnt_s=cnt_s+1;
                    
                end
                figure; hold on; 
                bar(mean([FD_lib;FD_con;DV;DV_thrQ;DV_thr],2))
%--------------------------------------------------------------------------                
                ED=ED(idx);
%---------------------------------------------
                QCFC=diag(corr(repmat(FD_Null',1,fge),liemat'));


                DVvFC=diag(corr(repmat(DV',1,fge),liemat_dv'));
                DVqvFC=diag(corr(repmat(DV_thrQ',1,fge),liemat_dv'));
                DVthrvFC=diag(corr(repmat(DV_thr',1,fge),liemat_dv'));
                
                FDlibvFC=diag(corr(repmat(FD_lib',1,fge),liemat_fd_lib'));
                FDconvFC=diag(corr(repmat(FD_con',1,fge),liemat_fd_con'));
                
                clear *liemat*
%--------------------------------------------------------------------------
                %ScatterBoxPlots(FDvFC)
                fh0=figure; 
                title([st{1} ' QC-FC, FD $\&$ DVARS'],'interpreter','latex')
                hold on; box on; 
                
                [sED,idx_ED]=sort(ED);
                rndidx=datasample(1:numel(DVvFC),round(numel(DVvFC)./10));
%------------------------------------------
       
                sQCFC=QCFC(idx_ED);
                scatter(sED(rndidx),sQCFC(rndidx),'marker','.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
                
%--------------------                
                sDVvFC=DVvFC(idx_ED);
                scatter(sED(rndidx),sDVvFC(rndidx),'marker','.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
                
                sDVqvFC=DVqvFC(idx_ED);
                scatter(sED(rndidx),sDVqvFC(rndidx),'marker','.','MarkerFaceColor',Dcol,'MarkerEdgeColor',Dcol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
                
                sDVthrvFC=DVqvFC(idx_ED);
                scatter(sED(rndidx),sDVthrvFC(rndidx),'marker','.','MarkerFaceColor',Dcol,'MarkerEdgeColor',Dcol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);

%--------------------                                
                sFDvFC_lib=FDlibvFC(idx_ED);
                scatter(sED(rndidx),sFDvFC_lib(rndidx),'marker','.','MarkerFaceColor',FDCol,'MarkerEdgeColor',FDCol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
                
                sFDvFC_con=FDconvFC(idx_ED);
                scatter(sED(rndidx),sFDvFC_con(rndidx),'marker','.','MarkerFaceColor',FDCol,'MarkerEdgeColor',FDCol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
%------------------------------------------ 
                QCFC_line = smooth(sED,sQCFC,smthspan,'rlowess');
                lhn         = line(sED,QCFC_line,'linewidth',1.5,'color','b');
%--------------------  
                DVvFC_line = smooth(sED,sDVvFC,smthspan,'rlowess');
                lh0         = line(sED,DVvFC_line,'linewidth',1.5,'color','k');
                
                DVqvFC_line = smooth(sED,sDVqvFC,smthspan,'rlowess');
                lh1         = line(sED,DVqvFC_line,'linewidth',1.5,'color',[Dcol lalp]);
                
                DVthrvFC_line = smooth(sED,sDVthrvFC,smthspan,'rlowess');
                lh2         = line(sED,DVthrvFC_line,'linewidth',1.5,'color',[Dcol lalp]);
%--------------------                
                FDlibvFC_line = smooth(sED,sFDvFC_lib,smthspan,'rlowess');
                lh3         = line(sED,FDlibvFC_line,'linewidth',1.5,'color',[FDCol lalp]);                
                
                FDconvFC_line = smooth(sED,sFDvFC_con,smthspan,'rlowess');
                lh4         = line(sED,FDconvFC_line,'linewidth',1.5,'color',[FDCol lalp]);   
                
                ylabel(['QC-FC'],'fontsize',12)
                xlabel(['Inter-node Distance (voxels)'],'fontsize',12)
                %legend([hpQ0 hpQ1 hpL0 hpL1],{'FD Quad','%D-var Quad','FD Lienar','%D-var Linear'},'fontsize',12)
                
                set(gcf,'Color','w')
                export_fig(fh0,['FCfigs/' st{1} '/' st{1} '_QC-FC_FD_DVARS.pdf'])

%--------------------------------------------------------------------------                
                
%                 h1=figure; 
%                 title([st{1} ' QC-FC Global Var Comps'],'interpreter','latex')
%                 hold on; box on;  
%                 scatter(ED,gDVvFC,'marker','.','MarkerFaceColor',Dcol,'MarkerEdgeColor',Dcol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
%                 scatter(ED,gSVvFC,'marker','.','MarkerFaceColor',Scol,'MarkerEdgeColor',Scol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
%                 
%                 [xl_gDV,yl_gDV]=FitMeLine(ED,gDVvFC,[1 4]);
%                 [hpQ0]=plot(xl_gDV(2,:),yl_gDV(2,:),'color',Dcol,'linewidth',1.5);
%                 [hpL0]=plot(xl_gDV(1,:),yl_gDV(1,:),'color',Dcol,'linewidth',1.5,'linestyle','-.');
%                 
%                 [xl_gSV,yl_gSV]=FitMeLine(ED,gSVvFC,[1 4]);
%                 [hpQ1]=plot(xl_gSV(2,:),yl_gSV(2,:),'color',Scol,'linewidth',1.5);
%                 [hpL1]=plot(xl_gSV(1,:),yl_gSV(1,:),'color',Scol,'linewidth',1.5,'linestyle','-.');
%                 %line([min(ED) max(ED)],[0 0],'color','r','linewidth',1.3)
%                 
%                 ylabel(['QC-FC'],'fontsize',12)
%                 xlabel(['Inter-node Distance (voxels)'],'fontsize',12)
%                 legend([hpQ0 hpQ1 hpL0 hpL1],{'gD-var Quad','gS-var Quad','gD-var Linear','gS-var Linear'},'fontsize',12)
%                 
%                 set(gcf,'Color','w')
%                 export_fig(h1,['FCfigs/' st{1} '/' st{1} '_QC-FC_GlobalVarComps.pdf'])

%--------------------------------------------------------------------------                
                
%                 h2=figure; 
%                 title([st{1} ' QC-FC NonGlobal Var Comps'],'interpreter','latex')
%                 hold on; box on;  
%                 scatter(ED,ngDVvFC,'marker','.','MarkerFaceColor',Dcol,'MarkerEdgeColor',Dcol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
%                 scatter(ED,ngSVvFC,'marker','.','MarkerFaceColor',Scol,'MarkerEdgeColor',Scol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
%                 
%                 [xl_ngDV,yl_ngDV]=FitMeLine(ED,ngDVvFC,[1 4]);
%                 [hpQ0]=plot(xl_ngDV(2,:),yl_ngDV(2,:),'color',Dcol,'linewidth',1.5);
%                 [hpL0]=plot(xl_ngDV(1,:),yl_ngDV(1,:),'color',Dcol,'linewidth',1.5,'linestyle','-.');
%                 
%                 [xl_ngSV,yl_ngSV]=FitMeLine(ED,ngSVvFC,[1 4]);
%                 [hpQ1]=plot(xl_ngSV(2,:),yl_ngSV(2,:),'color',Scol,'linewidth',1.5);
%                 [hpL1]=plot(xl_ngSV(1,:),yl_ngSV(1,:),'color',Scol,'linewidth',1.5,'linestyle','-.');
%                 
%                 %line([min(ED) max(ED)],[0 0],'color','r','linewidth',1.3)
%                 ylabel(['QC-FC'],'fontsize',12)
%                 xlabel(['Inter-node Distance (voxels)'],'fontsize',12)
%                 legend([hpQ0 hpQ1 hpL0 hpL1],{'ngD-var Quad','ngS-var Quad','ngD-var Linear','ngS-var Linear'},'fontsize',12)
%                 
%                 set(gcf,'Color','w')
%                 export_fig(h2,['FCfigs/' st{1} '/' st{1} '_QC-FC_NonGlobalVarComps.pdf'])
                
%--------------------------------------------------------------------------                
                
                
%                 h3=figure; hold on; box on; 
%                 title([st{1} ' - Correlation to Functional Connectivity'],'interpreter','latex')
%                 ScatterBoxPlots([[ngDVvFC ngSVvFC gDVvFC gSVvFC DVvFC FDvFC]],'figure',h3,'PointSize',3);
%                 h3.Children.XTickLabel={'ngD-var v. FC';'ngS-var v. FC';'gD-var v. FC';'gS-var v. FC';'$\%$D-var v. FC';'FD v. FC'};
%                 yyy=h3.Children.XAxis; yyy.FontSize=12; yyy.TickLabelInterpreter='Latex';
%                 yyy.TickLabelRotation=45;
%                 
%                 set(gcf,'Color','w')
%                 export_fig(h3,['FCfigs/' st{1} '/' st{1} '_Corrs.pdf'])
                
                
%--------------------------------------------------------------------------                
                
%                 h4=figure('position',[50,500,300,900]); hold on; 
%                 subplot(3,1,1)
%                 hold on; box on;
%                 title([st{1} ' Corr(FC,{FD,%D-var})'],'fontsize',12)
%                 histogram(FDvFC,50,'normalization','probability')
%                 histogram(DVvFC,50,'normalization','probability')
%                 legend({'FD','%D-var'},'fontsize',12)
%                 line([0 0],[0 0.1],'color','r','linewidth',1.2)
%                 ylabel('Intensity','fontsize',12)
%                 xlabel('Corr Coeffs','fontsize',12)
%                 
%                 subplot(3,1,2)
%                 hold on; box on;
%                 title([st{1} ' Corr(FC,GlobalVarComps)'],'fontsize',12)
%                 histogram(gDVvFC,50,'normalization','probability')
%                 histogram(gSVvFC,50,'normalization','probability')
%                 legend({'gD-var','gS-var'},'fontsize',12)
%                 line([0 0],[0 0.1],'color','r','linewidth',1.2)
%                 ylabel('Intensity','fontsize',12)
%                 xlabel('Corr Coeffs','fontsize',12)
%                 
%                 subplot(3,1,3)
%                 hold on; box on;
%                 title([st{1} ' Corr(FC,NonGlobalVarComps)'],'fontsize',12)
%                 histogram(ngDVvFC,50,'normalization','probability')
%                 histogram(ngSVvFC,50,'normalization','probability')
%                 legend({'ngD-var','ngS-var'},'fontsize',12)
%                 line([0 0],[0 0.1],'color','r','linewidth',1.2)
%                 ylabel('Intensity','fontsize',12)
%                 xlabel('Corr Coeffs','fontsize',12)
%                 
%                 set(gcf,'Color','w')
%                 export_fig(h4,['FCfigs/' st{1} '/' st{1} '_CorrHists.pdf'])
%--------------------------------------------------------------------------                
                 %pause(2)
                 %close all;
                 clear *_fd *_dv dv_* fd_* DV FD
        end
    end
end