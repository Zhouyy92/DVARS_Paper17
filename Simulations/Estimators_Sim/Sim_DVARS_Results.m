%
% This function produce test statistics out of the simulated moments
% already computed via Sim_DVARS_Inference.m function. 
%
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


nRlz=1000;

MeanNms={'xbar','sig2bar','sig2med','sigbar2','med'};
VarNms={'S2','IQR','hIQR','IRQ2','hIRQ2','IRQ3','hIRQ3','IRQ4','hIQR4'};


TT=[100,200,600,1200];
SDrng_List=[200 200; 200 250; 200 500];


Zstat =  @(x,m,s) (x-m)/s;
Zp =  @(x,m,s) 1-normcdf(Zstat(x,m,s));
X2stat =  @(x,m,s) 2*m/s^2*x;
X2df =  @(m,s) 2*m^2/s^2;
X2p =  @(x,m,s) 1-chi2cdf(X2stat(x,m,s),X2df(m,s));


for t_cnt=1:numel(TT)
    
    for vh_cnt=1:size(SDrng_List,1)
        
        T=TT(t_cnt); 
        SDrng=SDrng_List(vh_cnt,:);
        
        disp([num2str(T) ' - ' num2str(diff(SDrng))])
        
        load(['R/Sim_DVARS_Inf_t' num2str(T) '_vh' num2str(diff(SDrng)) '.mat']);
        
        %--Rejection Rate
        RRalpha=[0.1 0.05 0.01 0.005];
        %RRNms={'p1','p05','p01','p005'};
        %TestNms={'s2m_Dh','Z_s2m_Dh','s2b_Dh','Z_s2b_Dh','Dm_Dh','Z_Dm_Dh'}; % Every other is Z
        %RR=zeros(nRlz,length(RRalpha),length(TestNms));

        %--All pvalues, test stats
        %Pval=zeros(nRlz,length(TestNms),T-1);
        %ZStat=zeros(nRlz,length(TestNms)/2,T-1); 
        
        m2_cnt=1;
        % sig2bar sig2med med
        for m_cnt=[3 5] %[3 5] is only for the paper, [2 3 5] for SM
            v2_cnt=1;
            % 'IQR','hIQR','IRQ2','hIRQ2','IRQ3','hIRQ3','IRQ4','hIQR4'
            for v_cnt=[3 7] %[3 7] is only for the paper, 2:9 for SM
            
                disp([MeanNms{m_cnt} ' - ' VarNms{v_cnt}])
                
                %Mean_all{m2_cnt,v2_cnt}=Mean;
                %Var_all{m2_cnt,v2_cnt}=Var;
                
                for r=1:nRlz       
                    % Tests -- don't bother, you can do them later:
                    M_DV2=Mean(r,m_cnt);            % Median of voxelwise sig2
                    S_DV2=sqrt(Var(r,v_cnt));       % half-IQR of DVARS2
                    
                    DVARS2=DVARS2_all(:,r)';
                    
                    Zval{m2_cnt,v2_cnt}(r,:)=Zstat(DVARS2,M_DV2,S_DV2);
                    
                    Pval_tmp=[X2p(DVARS2,M_DV2,S_DV2);Zp(DVARS2,M_DV2,S_DV2)];
                    Pval{m2_cnt,v2_cnt}(r,:,:)=Pval_tmp;
                    
                    for i=1:length(RRalpha)
                        RR{m2_cnt,v2_cnt}(r,i,:)=mean(Pval_tmp<=RRalpha(i),2);
                    end
                    
                    
                    clear Pval_tmp DVARS2 p_Zp p_X2p
                end
                                
                v2_cnt=v2_cnt+1;
            end
            m2_cnt=m2_cnt+1;
        end
        save(['R/RR_Zval_Pvals_t' num2str(T) '_vh' num2str(diff(SDrng)) '_Paper.mat'],'-v7.3','NormP','RR','Pval','Zval','Mean','Var','DVARS2_all','VaTrue','MnTrue')
        clear RR Pval Zval Mean Var DVARS2_all VaTrue MnTrue
    end
end

