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

addpath DVARS


IntnstyScl  = @(Y,md,scl) (Y./md).*scl; 

splflag     = 1;

I           = 9000; 
T_List      = [150,250,600,1200];
SpksRt_List = [0.05 0.10 0.30];
rhorng      = 0:0.2:0.6;
SD_List     = [1 10 100]; %leave me alone for now!!
B0          = 0; %it is the grand mean of the time series

nRlz = 200;
%----------------------------------------
% I           = 9000;

T_List      = [150,250,600,1200];
SpksRt_List = [1/3 1/5 1/10 1/100]; %[20% 10% 2% 1%]

rhorng      = 0:0.2:0.6;
SD_List     = [1 10 100];
%nRlz        = 5;
%----------------------------------------

for noiseflag=[0 1]
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

                    disp(['NoiseFlag: ' num2str(noiseflag) ', T: ' num2str(T) ', rho: ' num2str(rho) ', SD: ' num2str(SD) ', I: ' num2str(I) ', #SpksRt: ' num2str(SpksRt)])

        %----------------------------------------
                     
                    clear *Dvar *Svar magspk *Avar Y sY sDeltapDvar spvals sDVARS rng0* idxspike_tmp idxspike idxspike_tmp0
                    idxspike = zeros(nRlz,nspike);
                    magspk   = zeros(nRlz,nspike);
                    
                    for it=1:nRlz
                        Y           = B0+sqrtm(ACMat)*randn(T,I)*SD+randn(T,I)*SD*noiseflag;
                        Y           = Y'; %because everything else, apart from MassAC, is in IxT form.
                        %---
                        %figure; imagesc(full(spm_Q(rho,T)))
                        %------------------
                        %ACv_T(:,it) = MassAC(Y',1);
                        %------------------
                        %Y=randn(I,T)+10; %white noise
                        %----------------------------------------
                        % Add outliers
                        %pBad = 0.2;  % Proportion of bad obs
                        %nBad = round(pBad*T);
                        %Y(randperm(T,nBad)',:) = randn(nBad,1)*SD*5;
                        %----------------------------------------
                        %scY=IntnstyScl(Y,median(mean(Y)),100); %

                        [V_tmp,DSE_Stat_tmp] = DSEvars(Y,'verbose',0);

                        %V_DSE{nrho,it}=V_tmp;
                        %DSE_Stat{nrho,it}=DSE_Stat_tmp;

                        gAvar(it,:) = V_tmp.g_Avar_ts;
                        gDvar(it,:) = V_tmp.g_Dvar_ts;
                        gSvar(it,:) = V_tmp.g_Svar_ts;

                        Avar(it,:)  = V_tmp.Avar_ts;
                        Dvar(it,:)  = V_tmp.Dvar_ts;
                        Svar(it,:)  = V_tmp.Svar_ts;
                        %----------------------------------------
                        sY=Y;

                        magspk_tmp      = SD*(.5+(1-.5)*rand(1,nspike)+B0)./(B0+1);
                        magspk(it,:)    = magspk_tmp;

                        rng0            = 0:round(1/SpksRt):T;
                        rng0s           = rng0(1:end-1);
                        rng0e           = rng0(2:end)-1;
                        if splflag 
                            for ns=1:numel(rng0e)
                                wgts=[]; idxspike_tmp0=[]; sScan=[];
                                cnstnt=1; 
                                %rng0s(ns):rng0e(ns)
                                %[0 ones(1,(rng0e(ns)-rng0s(ns))-2) 0]
                                if ns==1
                                    %disp('head')
                                    wgts=[0 0 ones(1,(rng0e(ns)-rng0s(ns))-2)   0];
                                    if ~sum(wgts); continue; end;
                                    idxspike_tmp0 = sort(datasample(rng0s(ns):rng0e(ns),1,'Weights',wgts));
                                elseif ns==nspike && rng0e(end)>T-2
                                    %disp('tail')
                                    wgts=[0   ones(1,(rng0e(ns)-rng0s(ns))-2) 0 0];
                                    if ~sum(wgts); continue; end;
                                    idxspike_tmp0 = sort(datasample(rng0s(ns):rng0e(ns),1,'Weights',wgts));
                                else
                                    %disp('middle')
                                    wgts=[0   ones(1,(rng0e(ns)-rng0s(ns))-1)   0];
                                    if ~sum(wgts); continue; end;
                                    idxspike_tmp0 = sort(datasample(rng0s(ns):rng0e(ns),1,'Weights',wgts));
                                end
                                sScan=Y(:,idxspike_tmp0);
                                sY(:,idxspike_tmp0)=sScan+magspk_tmp(ns)*cnstnt;
                                idxspike(it,ns)=idxspike_tmp0;
                            end
                        end
                        %-----------
                        [sDVARS_tmp,sDVARS_Stat_tmp] = DVARSCalc(sY,'verbose',0);
                        sDVARS(it,:)                 = sDVARS_tmp;
                        spvals(it,:)                 = sDVARS_Stat_tmp.pvals;
                        sDeltapDvar(it,:)            = sDVARS_Stat_tmp.DeltapDvar;
                        %-----------
                        [sV_tmp,sDSE_Stat_tmp]       = DSEvars(sY,'verbose',0);

                        %sV{nrho,it}=sV_tmp;
                        %sDSE_Stat{nrho,it}=sDSE_Stat_tmp;
                        
                        sgAvar(it,:) = sV_tmp.g_Avar_ts;
                        sgDvar(it,:) = sV_tmp.g_Dvar_ts;
                        sgSvar(it,:) = sV_tmp.g_Svar_ts;

                        sAvar(it,:) = sV_tmp.Avar_ts;
                        sDvar(it,:) = sV_tmp.Dvar_ts;
                        sSvar(it,:) = sV_tmp.Svar_ts;
                        %----------------------------------------
                        [idx_iqr_tmp,~,dv] = DVARS_IQR(sDVARS_tmp,'DVARS');
                        sIdx_iqr{it}                 = idx_iqr_tmp;
                    end
    %----------------------------------------
                    save(['R\DSE_Sim_' num2str(T) 'l_' num2str(rho*100) 'p_' num2str(SD) 'vh_' num2str(nspike) 'nspk_noisef' num2str(noiseflag) '.mat'],'sIdx_iqr','nspike','idxspike','magspk','sAvar','sSvar','sDvar','Avar','Svar','Dvar','sgAvar','sgSvar','sgDvar','gAvar','gSvar','gDvar','sDVARS','spvals','sDeltapDvar')
                    clear *Dvar *Svar magspk idxspike *Avar Y sY sDVARS spvals sDeltapDvar sDeltapDvar spvals sDVARS *_tmp
                end
            end
        end

    end
end    

load('R/DSE_Sim_100l_60p_100vh_5nspk.mat')
figure; hold on; 
wo=randi(100);
plot(sDVARS(wo,:))
idx=find(spvals(wo,:)<(0.05/numel(sDVARS(wo,:))) & sDeltapDvar(wo,:)>5);
scatter(idx,ones(1,numel(idx))*min(sDVARS(wo,:)),'ro')
scatter(idxspike(wo,:),ones(1,nspike)*min(sDVARS(wo,:)),'kx')