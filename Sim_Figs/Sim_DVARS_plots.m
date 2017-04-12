%
% This file re-produce Fig 2 and 3 of the paper out of Sim_DVARS_Results.m
% files. 
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


TheseMeans=[1 2 3 5];

%load('R/RR_Zval_Pvals.mat')
nRlz=1000; I=90e3;
RRalpha=[0.1 0.05 0.01 0.005];
%TT=[100,200,600,1200];

TT=[100,200,600,1200];
TT_Clr=[0,0.4470,0.7410; 0.8500,0.3250,0.0980; 0.9290,0.6940,0.1250; 0.4940,0.1840,0.5560];
SDrng_List=[200 200; 200 250; 200 500];

MeanNms={'xbar','sig2bar','sig2med','sigbar2','med'};
%MeanNms_lab={'mean','$\bar{\sigma^2}$','$\tilde{\sigma^2}$','$\bar{\sigma}^2$','median'};

MeanNms_lab={'$\hat{\mu}_0^{\mbox{DVARS}}$','$\hat\mu_0^D$','$\tilde{\mu}_{0}^{D}$','$\tilde{\mu}_0^{\mbox{DVARS}}$'};
%

%MeanNms=MeanNms(TheseMeans);
%MeanNms_lab=MeanNms_lab(TheseMeans);

VarNms={'S2','IQR','hIQR','IRQ2','hIRQ2','IRQ3','hIRQ3','IRQ4','hIQR4'};

VarNms_lab={'$\hat\sigma_0^2$','IQR','hIQR','IRQ','hIRQ','IRQ','hIRQ','IRQ','hIQR'};

NormPNms={'I','rt2','rt3','rt4'};
RRNms={'p1','p05','p01','p005'};
RRNms_lab={'0.1','0.05','0.01','0.005'};

TestNms={'X2','Z'};
% Zstat =  @(x,m,s) (x-m)/s;
% Zp =  @(x,m,s) 1-normcdf(Zstat(x,m,s));
 X2stat =  @(x,m,s) 2*m/s^2*x;
 X2df =  @(m,s) 2*m^2/s^2;
% X2p =  @(x,m,s) 1-chi2cdf(X2stat(x,m,s),X2df(m,s));

Fig_Idx_pp=reshape(1:3*2,[2 3])';
Fig_Idx_bb=reshape(1:3*2,[2 3])';

lfts=14;

fig_cnt_pp=1;
fig_cnt_bb=50;

m2_cnt=1;
for m_cnt=[3 5] %[3 5] is only for the paper, [2 3 5] for SM
    v2_cnt=1;
    for v_cnt=[3 7] %[3 7] is only for the paper, 2:9 for SM
        %figure(fig_cnt)
%         figure
%         hold on; 
        %suptitle(['t ' num2str(T) ' -- ' MeanNms{m_cnt} '-' VarNms{v_cnt}])
        
        for vh_cnt=1:size(SDrng_List,1)
            SDrng=SDrng_List(vh_cnt,:);        
            
            for t_cnt=1:numel(TT)
                T=TT(t_cnt);

                load(['R/RR_Zval_Pvals_t' num2str(T) '_vh' num2str(diff(SDrng)) '_Paper.mat'],'Mean','MnTrue','NormP','Pval','RR','VaTrue','Var','Zval')
                
                Mean=Mean(:,TheseMeans); % This dude > sigbar2 is just for checking!
                
                m_bias=(mean(Mean,1)./MnTrue-1)*100;
                size(Mean)
                %v_bias=(sqrt(mean(Var,1)./VaTrue)-1)*100;
                v_bias=(mean((Var),1)./(VaTrue)-1)*100;
                %v_bias=(mean((Var),1)./(VaTrue)-1)*100;
                mnRR=squeeze(mean(RR{m2_cnt,v2_cnt})); % Dim: length(RRalpha) x length(TestNms)
                rr_bias=bsxfun(@times,bsxfun(@minus,mnRR,RRalpha'),1./RRalpha')*100;
                
                %MnTrue
                
                %VaTrue
                
                ffh=figure(fig_cnt_bb);
                
                %if m2_cnt==1 && v2_cnt==1


                    sp_1=subplot(3,2,Fig_Idx_bb(vh_cnt,1));
                    hold on; grid on; box on;
                    plot(m_bias,'color',TT_Clr(t_cnt,:),'marker','x','LineWidth',2)
                    set(sp_1,'XTick',1:numel(MeanNms_lab),'TickLabelInterpreter', 'latex','XTickLabel',MeanNms_lab,...
                        'FontSize',12,'XTickLabelRotation',45);
                    
                    ylim([-10 5])
                    
                    ylabel('Mean Bias \%','interpreter','latex','fontsize',lfts)
                    clear sp_1

                    sp_2=subplot(3,2,Fig_Idx_bb(vh_cnt,2));
                    hold on; grid on; box on;
                    plot(v_bias,'color',TT_Clr(t_cnt,:),'marker','x','LineWidth',2)
                    set(sp_2,'TickLabelInterpreter', 'latex');
                    sp_2.XTick=1:numel(VarNms_lab);
                    sp_2.XTickLabel=VarNms_lab;
                    sp_2.FontSize=12;
                    sp_2.XTickLabelRotation=90;
                    
                    ylim([-1 2])
                    
                    ylabel('Variance Bias \%','interpreter','latex','fontsize',lfts)


                    set(ffh,'Color','w')
                %end
                
%                 sp_3=subplot(3,3,Fig_Idx_bb(vh_cnt,3));
%                 hold on; grid on; box on;
%                 plot(rr_bias,'color',TT_Clr(t_cnt,:),'marker','x','LineWidth',2)
%                 set(sp_3,'TickLabelInterpreter', 'latex');
%                 sp_3.XTick=1:numel(RRNms);
%                 sp_3.XTickLabel=RRNms_lab;
%                 sp_3.FontSize=12;
%                 sp_3.XTickLabelRotation=90;
%                 ylabel('Rejection Rate Bias')
                
                
%                 fprintf('\n')
% 
%                 fprintf('\nSettings: nRlz=%d  I=%d  T=%d  SDrng=[%d,%d]\n',nRlz,I,T,SDrng);...
% 
%                 fprintf('\n  Mean Bias, Percent\n\n');...
%                 disp(array2table(m_bias,'VariableNames',MeanNms));...
% 
%                 fprintf('\n  sqrt(Variance Bias), Percent\n\n');...
%                 disp(array2table(v_bias,'VariableNames',VarNms));...
% 
%                 fprintf('\n  Reject Rate\n\n');...
%                 disp(array2table(mnRR,'VariableNames',TestNms,'RowNames',RRNms));...
% 
%                 fprintf('\n  Reject Rate Bias, Percent\n\n');...
%                 disp(array2table(rr_bias,'VariableNames',TestNms,'RowNames',RRNms));...
% 
%                 fprintf('\n  Normality P-value (mean)\n\n');...
%                 disp(array2table(mean(NormP,1),'VariableNames',NormPNms))

                figure(fig_cnt_pp)
%                 subplot(3,2,Fig_Idx(vh_cnt,1))
%                 plot((1:(T-1))/T,squeeze(median(sort(Pval{m2_cnt,v2_cnt},3),1)),'linewidth',2);
%                 hline=refline(1,0); hline.Color='k'; hline.LineStyle='-.';
%                 grid on
%                 axis tight
%                 legend(TestNms,'Location','SouthEast','Interpreter','none')
%                 %abline(0,1,'linestyle','-','color',[1 1 1]*.8,'linewidth',4)
%                 title('PP plot, median of all realisations')
                subplot(3,2,Fig_Idx_pp(vh_cnt,1))
                loglog((1:(T-1))/T,squeeze(median(sort(Pval{m2_cnt,v2_cnt},3),1)),'linewidth',2,'Color',TT_Clr(t_cnt,:));
                %loglog((1:(T-1))/T,squeeze(median(sort(Pval{m2_cnt,v2_cnt},3),1)),'linewidth',2,'Color',TT_Clr(t_cnt,:));
                hold on; grid on;axis tight
                hline=refline(1,0); hline.Color='k'; hline.LineStyle='-.';
                %legend(TestNms,'Location','SouthEast','Interpreter','none')
                %abline(0,1,'linestyle','-','color',[1 1 1]*.8,'linewidth',4)
                title('log P-P plot','interpreter','latex','fontsize',lfts)
% 
                %figure(2)
%                 subplot(2,2,3)
%                 h=histogram(DVARS2_all(:,1000),50,'Norm','pdf');
%                 M_DV2 = Mean(1000,3);             % Median of voxelwise sig2
%                 S_DV2 = sqrt(Var(1000,2));       % half-IQR of DVARS2
%                 c=X2stat(1,M_DV2,S_DV2);
%                 nu=X2df(M_DV2,S_DV2);
%                 xx=get(gca,'xlim');xx=linspace(xx(1),xx(2),100);
%                 hold on 
%                 h(2)=plot(xx,chi2pdf( xx   *c,nu)*c,'linewidth',4);
%                 hold off
%                 legend(h,{'DVARS','Chi2(s2m_Dh)'},'Interpreter','none')
%                 title('Emperical and fitted distribution, last realisation')
% 
                subplot(3,2,Fig_Idx_pp(vh_cnt,2))
                DefCol=get(groot,'defaultAxesColorOrder');
                h(1)=histogram(Zval{m2_cnt,v2_cnt},100,'Norm','pdf','DisplayStyle','Stairs','LineWidth',2,'edgecolor',TT_Clr(t_cnt,:));
                hold on; box on; grid on;
                %h(2)=histogram(Zval{t_cnt,vh_cnt}{m2_cnt,v2_cnt},100,'Norm','pdf','DisplayStyle','Stairs','LineWidth',2,'edgecolor',DefCol(2,:));
                %h(3)=histogram(Zval{t_cnt,vh_cnt}{m2_cnt,v2_cnt},100,'Norm','pdf','DisplayStyle','Stairs','LineWidth',2,'edgecolor',DefCol(3,:));
                xx=get(gca,'xlim');xx=linspace(xx(1),xx(2),100);
                h(2)=plot(xx,normpdf(xx),'linewidth',4,'color',[1 1 1]*.8,'linestyle','-.');

                %legend(h,['ZStat',{'N(0,1)'}],'Location','SouthEast','Interpreter','none')
                title('Emperical and null Z tests','interpreter','latex','fontsize',lfts)
                
                clear Mean NormP RR Var Zval *_bias mnRR
            end
        %export_fig(['fig_' MeanNms{m_cnt} '_' VarNms{v_cnt} '_' num2str(T) '.png'])
        end
        fig_cnt_pp=fig_cnt_pp+1;
        %fig_cnt_bb=fig_cnt_bb+1;
        set(gcf,'Color','w')
    
        v2_cnt=v2_cnt+1;
    end
    m2_cnt=m2_cnt+1;
end