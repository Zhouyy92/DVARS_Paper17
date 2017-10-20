%
% SA, UoW, 2017
% srafyouni@gmail.com 
%
%%%REFERENCES
%
%   Afyouni S. & Nichols T.E., Insights and inference for DVARS, 2017
%   http://www.biorxiv.org/content/early/2017/04/06/125021.1  
%

clear

I               = 9000; 
T_List          = [150,250,600,1200];
SpksRt_List     = [1/100,1/10,1/5,1/3]; %[20% 10% 2% 1%]
rhorng          = [0,0.2,0.4,0.6];
SD_List         = [1,10,100];
noiseflag_list  = [1];%,1];

nRlz=200;

%--TEST
%SD=1; T=150; nspike=2;
%----------------------------------------
for noiseflag=noiseflag_list
    for t_cnt=1:numel(T_List)
        clear T
        T = T_List(t_cnt);
        for ns_cnt=1:numel(SpksRt_List)
            clear SpksRt nspike
            SpksRt      = SpksRt_List(ns_cnt);
            nspike      = round(T*SpksRt);
            for s_cnt=1:numel(SD_List)
                SD = SD_List(s_cnt);

                for nrho=1:numel(rhorng)
                %----------------------------------------
                    rho     = rhorng(nrho);
                    ACMat   = full(spm_Q(rho,T)); %AC matrix
                    %disp(['NoiseFlag: ' num2str(noiseflag) ', T: ' num2str(T) ', rho: ' num2str(rho) ', SD: ' num2str(SD) ', I: ' num2str(I) ', #SpksRt: ' num2str(SpksRt)])
                    load(['/Users/sorooshafyouni/GoogleDrive/FromIDH/DSE_Sim/Code/R/DSE_Sim_' num2str(T) 'l_' num2str(rho*100) 'p_' num2str(SD) 'vh_' num2str(nspike) 'nspk_noisef' num2str(noiseflag) '.mat'])
                    for it=1:nRlz
                        clear dtctd_idx fp insrt_idx fp
                        insrt_idx   =   nonzeros(idxspike(it,:));
                        %-------------------------------------
                        dtctd_idx   =   find(spvals(it,:)<0.05./(T-1));% & sDeltapDvar(it,:)>5);
                        %-------------------------------------
                        fp=setdiff(dtctd_idx,insrt_idx);
                        dscnt=0;
                        if ~isempty(fp)
                            for i=1:numel(fp); if ismember(fp(i)+1,insrt_idx); dscnt=dscnt+1; end; end
                        end
                                                
                        bct_tp           =   numel(intersect(insrt_idx,dtctd_idx));
                        bct_fp           =   numel(setdiff(dtctd_idx,insrt_idx))-dscnt;
                        bct_fn           =   numel(setdiff(insrt_idx,dtctd_idx));
                        bct_tn           =   T-1-numel(insrt_idx);
                        
                        dvt_Sens_tmp(it) =   bct_tp/(bct_tp+bct_fn);
                        dvt_Spec_tmp(it) =   bct_tn/(bct_tn+bct_fp);
                        
%--------------------------------------------------------------------------                        
                        
                        %POWER----------------------------------------
                        clear dtctd_idx bct_* dscnt fp
                        dtctd_idx       =   sIdx_iqr{it};
                        
                        fp=setdiff(dtctd_idx,insrt_idx);
                        dscnt=0;
                        if ~isempty(fp)
                            for i=1:numel(fp); if ismember(fp(i)+1,insrt_idx); dscnt=dscnt+1; end; end
                        end
                        
                        bct_tp           =   numel(intersect(insrt_idx,dtctd_idx));
                        bct_fp           =   numel(setdiff(dtctd_idx,insrt_idx))-dscnt;
                        bct_fn           =   numel(setdiff(insrt_idx,dtctd_idx));
                        bct_tn           =   T-1-numel(insrt_idx);
                        
                        fsl_Sens_tmp(it) =   bct_tp/(bct_tp+bct_fn);
                        fsl_Spec_tmp(it) =   bct_tn/(bct_tn+bct_fp);
                    end
                    dvt_Sens{t_cnt,s_cnt}(ns_cnt,nrho)=mean(dvt_Sens_tmp);
                    dvt_Spec{t_cnt,s_cnt}(ns_cnt,nrho)=mean(dvt_Spec_tmp);
                    
                    fsl_Sens{t_cnt,s_cnt}(ns_cnt,nrho)=mean(fsl_Sens_tmp);
                    fsl_Spec{t_cnt,s_cnt}(ns_cnt,nrho)=mean(fsl_Spec_tmp);
                    
                    clear *_tmp
                 %----------------------------------------
                end
            end
        end

    end
end    
set(gcf,'Color','w')
%----------------------------------------
BCT={dvt_Sens,fsl_Sens;...
    dvt_Spec,fsl_Spec};
bct_List={'Sensitivity','Specificity'};
test_List={'DVARS Test','DVARS FSL'};
tl_cnt=1;
for tl=test_List
    bct_cnt=1;
    for bct=bct_List
        fh0=figure('position',[50,500,700,1500]); 
        hold on;
        
        suptitle(tl{1})
        
        sp_cnt=1;
        for t_cnt=1:numel(T_List)
            T = T_List(t_cnt);
            for s_cnt=1:numel(SD_List)
                SD = SD_List(s_cnt);

                sp0=subplot(4,3,sp_cnt);
                hold on; box on; grid on; 
                title(['T= ' num2str(T) '- VH=' num2str(SD)])

                for i=1:numel(rhorng)
                    plot(BCT{bct_cnt,tl_cnt}{t_cnt,s_cnt}(:,i)*100,'linewidth',1.2,'marker','d','color',[[0.15 0.15].*(5-i) 1])
                end
                if bct_cnt==1
                    ylim([-1 105])
                elseif bct_cnt==2
                    ylim([99 100])
                end
                    
                xlim([1 4])
                sp0.XTickLabel={'1%','10%','20%','30%'};

                if t_cnt==numel(T_List) 
                    xlabel('Spikes (% of T)')
                end

                if s_cnt==1 && bct_cnt==1
                    ylabel('Sensitivity %')
                elseif s_cnt==1 && bct_cnt==2
                    ylabel('Specificity %')
                end
                
                if t_cnt==numel(T_List) && s_cnt==numel(SD_List)
                   legend({'AC=0','AC=0.2','AC=0.4','AC=0.6'}) 
                end
                
                sp_cnt=sp_cnt+1;

            end
        end
        
        set(fh0,'color','w')
        export_fig(fh0,['SpkDtctn_' tl{1} '_' bct{1} '.pdf'])
        
        bct_cnt=bct_cnt+1;
    end
    tl_cnt=tl_cnt+1;
end
%----------------------------------------