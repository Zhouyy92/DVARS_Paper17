clear;

close all

Mode={'PreFix'};

Site={'HCP'};
AtlasList={'Power'};
MarkerAlp=0.1;

addpath /Users/sorooshafyouni/Home/GitClone/DVARS/Nifti_Util
%----------------------------------------
%21 Healthy 
load(['HCP_100Unrel_SubList.mat'])
HCP_10Unrel_SubList(3)=[];
SubList=HCP_10Unrel_SubList;
%----------------------------------------
GSRStatList={'NoGSR'}; %%BE CAREFULLLLL****************
StateList={'DVCleaned'};%'Minimal','FDCleaned'};% raw DVCleaned
smthspan=0.05;
lalp=0.7;
%----------------------------------------
Col=get(groot,'defaultAxesColorOrder');
Acol=Col(5,:); % Green
%Dcol=Col(1,:); % Blue
Dcol=[0.5 0.5 0.5];
Scol=Col(3,:); % Yellow
FDCol=[1 0 0];
%----------------------------------------
for st=StateList
    for a=AtlasList
        for g=GSRStatList
                
                disp('Porportions & NetMats...')
            
                DV=[]; DV_thrQ=[]; DV_thr=[]; 
                FD_con=[]; FD_lib=[]; FD_Null=[];
                cnt_s=1;
                for s=SubList
                    load(['/Volumes/HCP_S900/HCP_100Unrelated/' Mode{1} '_NetMats/' a{1} '/NetMats/HCP_MinimalF_' s{1} '_' a{1} '_NetMats_DVC.mat'])
                                                                                                   
                    nn  = size(mat,1);
                    fge=nn*(nn-1)./2;
                    idx=find(triu(ones(nn),1));
                    
                    if isempty(rgs_key) || numel(rgs_key)>round(1200*0.30)
                        disp([s{1}  ' exluded.']);
                        continue;
                    end
%---------------------------------------------              
                        liemat(:,cnt_s)=mat(idx);
                        %Porpottion of scrubbed-------
                        %FD_Null(cnt_s)=numel(fd_rmv_lib); 
                        %FD_Null(cnt_s)=numel(dv_rmv_pr); 
                        FD_Null(cnt_s)=mean(FDts);
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
                        DV(cnt_s)=numel(dv_rmv); %check this later!
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
                
                disp('Correlations...')
                
%--------------------------------------------------------------------------                
                ED   = 2*ED(idx);
%---------------------------------------------
                smtp = 1;
                
                ED   = ED(1:smtp:fge);
                
                QCFC.QCFC       = diag(corr(repmat(FD_Null'   ,1,numel(1:smtp:fge)),liemat(1:smtp:end,:)'));

                QCFC.DVvFC      = diag(corr(repmat(DV'       ,1,numel(1:smtp:fge)),liemat_dv(1:smtp:end,:)'));
                QCFC.DVqvFC     = diag(corr(repmat(DV_thrQ' ,1,numel(1:smtp:fge)),liemat_dv(1:smtp:end,:)'));
                QCFC.DVthrvFC   = diag(corr(repmat(DV_thr',1,numel(1:smtp:fge)),liemat_dv(1:smtp:end,:)'));
                
                QCFC.FDlibvFC   = diag(corr(repmat(FD_lib',1,numel(1:smtp:fge)),liemat_fd_lib'));
                QCFC.FDconvFC   = diag(corr(repmat(FD_con',1,numel(1:smtp:fge)),liemat_fd_con'));
                
%--------------------------------------------------------------------------
                %ScatterBoxPlots(FDvFC)
                fh0=figure; 
                %title([st{1} ' QC-FC, FD $\&$ DVARS'],'interpreter','latex')
                hold on; box on; 
                
                for i=-0.15:0.05:0.15
                    line([0 max(ED)+10],[i i],'color',[0.5 0.5 0.5 0.5])
                end

                disp('Sorting & scatter plotting...')
                
                [sED,idx_ED]=sort(ED);
                
                rndidx=datasample(1:numel(QCFC.DVvFC),round(numel(QCFC.DVvFC)./10));
%------------------------------------------
                sQCFC=QCFC.QCFC(idx_ED);
                scatter(sED(rndidx),sQCFC(rndidx),'marker','.','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
                
%--------------------                
                sDVvFC=QCFC.DVvFC(idx_ED);
                scatter(sED(rndidx),sDVvFC(rndidx),'marker','.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
                
                sDVqvFC=QCFC.DVqvFC(idx_ED);
                scatter(sED(rndidx),sDVqvFC(rndidx),'marker','.','MarkerFaceColor',Dcol,'MarkerEdgeColor',Dcol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
                
                sDVthrvFC=QCFC.DVqvFC(idx_ED);
                scatter(sED(rndidx),sDVthrvFC(rndidx),'marker','.','MarkerFaceColor',Dcol,'MarkerEdgeColor',Dcol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);

%--------------------                                
                sFDvFC_lib=QCFC.FDlibvFC(idx_ED);
                scatter(sED(rndidx),sFDvFC_lib(rndidx),'marker','.','MarkerFaceColor',FDCol,'MarkerEdgeColor',FDCol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);
                
                sFDvFC_con=QCFC.FDconvFC(idx_ED);
                scatter(sED(rndidx),sFDvFC_con(rndidx),'marker','.','MarkerFaceColor',FDCol,'MarkerEdgeColor',FDCol,'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',MarkerAlp);

                                
%------------------------------------------ 
                
                disp('Smoothing & lines...')                
                
                Lines.QCFC_line     = smooth(sED,sQCFC,smthspan,'lowess');
                lhn_fd              = line(sED,Lines.QCFC_line,'linewidth',1.5,'color',[0 0 1 lalp]);
%--------------------  
                Lines.DVvFC_line    = smooth(sED,sDVvFC,smthspan,'lowess');
                lh0_dv              = line(sED,Lines.DVvFC_line,'linewidth',1.5,'color','k');
                
                Lines.DVqvFC_line   = smooth(sED,sDVqvFC,smthspan,'lowess');
                lh1_dv              = line(sED,Lines.DVqvFC_line,'linewidth',1.5,'color',[Dcol lalp],'linestyle','-.');
                
                Lines.DVthrvFC_line = smooth(sED,sDVthrvFC,smthspan,'lowess');
                lh2_dv              = line(sED,Lines.DVthrvFC_line,'linewidth',1.5,'color',[Dcol lalp],'linestyle',':');
%--------------------                
                Lines.FDlibvFC_line = smooth(sED,sFDvFC_lib,smthspan,'lowess');
                lh3_df              = line(sED,Lines.FDlibvFC_line,'linewidth',1.5,'color',[FDCol lalp],'linestyle','-.');                
                
                Lines.FDconvFC_line = smooth(sED,sFDvFC_con,smthspan,'lowess');
                lh4_df              = line(sED,Lines.FDconvFC_line,'linewidth',1.5,'color',[FDCol lalp],'linestyle',':');   
                
                ylabel(['Correlation Coefficients (Between \#Spikes and FC)'],'fontsize',12,'interpreter','Latex')
                xlabel(['Inter-node Distance (mm)'],'fontsize',12,'interpreter','Latex')
                legend([lhn_fd lh0_dv lh1_dv lh2_dv lh3_df lh4_df],{'Unscrubbed','DVARS Test','DVARS IQR','DVARS Thresholded','FD-Lenient Threshold','FD-Conservative Threshold'},'fontsize',12)
                
                xlim([0 max(ED)+10])
                
                set(gcf,'Color','w')
                title(['QC-FC: MPP data, 100 HCP Subjects, ' a{1} ' atlas'],'fontsize',12,'interpreter','latex')
                export_fig(fh0,['FCfigs/' st{1} '/' st{1} '_' a{1} '_' Mode{1} '_100Unrel_QC-FC_FD_DVARS.pdf'])

%--------------------------------------------------------------------------

                QCFC.DoF=[DV;DV_thrQ;DV_thr;FD_lib;FD_con;];

                SLineCol=[0 0 0; Dcol; Dcol; FDCol; FDCol];
                dof_fh=figure; hold on; box on; 
                ScatterBoxPlots(sqrt(QCFC.DoF'),'Color',SLineCol,'MedLineColor','b','figure',dof_fh);
                ylabel('$\sqrt{\mathrm{ \# of  Detected  Spikes}}$','fontsize',12,'interpreter','Latex')
                xlabel('Methods','fontsize',12,'interpreter','Latex')
                dof_fh.Children.XTickLabel={'DVARS Test';'DVARS IQR';'DVARS Thresholded';'FD-Lenient Threshold';'FD-Conservative Threshold'};
                dof_fh.Children.XTickLabelRotation=45;
                
                set(gcf,'Color','w')
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
%                 export_fig(h2,['FCfigs/' st{1} '/' st{1} '_' a{1} '_QC-FC_NonGlobalVarComps.pdf'])
                
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
                 
                 %save(['R/QCFC_100Unrelated_' a{1} '_MF.mat'],'sED','Lines','QCFC','ED')
                 
                 %clear *_fd *_dv dv_* fd_* DV FD *liemat*
        end
    end
end