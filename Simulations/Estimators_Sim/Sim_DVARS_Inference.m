function Sim_DVARS_Inference(nRlz,T_Idx,SD_idx)
% This function is to compute the first two moments of DVARS2 and is called
% from the WorkBench on a multi-processor cluster for more than 1,000 for
% each combination of TT and SDrng_List (3x4=12).
%
%%%EXAMPLE
%
%  % This produce (and saves) 10e2 realisation of a DVARS2 of 90e3 timeseries of size 200
%  % with high level of variance hetrogeniety.
%
%  Sim_DVARS_Inference(10e2,2,3) 
%
%%%REFERENCES
%
%   Afyouni S. & Nichols T.E., Insights and inference for DVARS, 2017
%   http://www.biorxiv.org/content/early/2017/04/06/125021.1  
%
%   SA, UoW, 2017
%   srafyouni@gmail.com
%

% TT=[100,300,500,1200]; :: T_Idx
% SDrng_List=[200 200; 200 250; 200 500]; :: SD_idx

%nRlz=100;

I=90000; 
TT=[100,200,600,1200];
%TT=[100,200,800,1200];
T=TT(T_Idx); 

SDrng_List=[200 200; 200 250; 200 500];
SDrng=SDrng_List(SD_idx,:);

disp(['Settings: T:' num2str(T) ' SD:' num2str(diff(SDrng)) ' #Real:' num2str(nRlz) ' #voxels:' num2str(I)])

%% Body
MeanNms={'xbar','sig2bar','sig2med','sigbar2','med'};
Mean=zeros(nRlz,length(MeanNms));

VarNms={'S2','IQR','hIQR','IRQ2','hIRQ2','IRQ3','hIRQ3','IRQ4','hIQR4'};
Var=zeros(nRlz,length(VarNms));

NormPNms={'I','rt2','rt3','rt4'};
NormP=zeros(nRlz,length(NormPNms));

% --Rejection Rate
%RRalpha=[0.1 0.05 0.01 0.005];
%RRNms={'p1','p05','p01','p005'};
%TestNms={'s2m_Dh','Z_s2m_Dh','s2b_Dh','Z_s2b_Dh','Dm_Dh','Z_Dm_Dh'}; % Every other is Z
%RR=zeros(nRlz,length(RRalpha),length(TestNms));

% --All pvalues, test stats
% Pval=zeros(nRlz,length(TestNms),T-1);
% ZStat=zeros(nRlz,length(TestNms)/2,T-1); 

% --(Possibly) heterogeneous variance
SD=SDrng(1)+diff(SDrng)*randn(I,1);
KSD=spdiags(SD,0,I,I);

Y=randn(I,T); %white noise
Y=KSD*Y;

IQRsd = @(x) (quantile(x,0.75)-quantile(x,0.25))/1.349;
hIQRsd = @(x) (quantile(x,0.5)-quantile(x,0.25))/1.349*2;

% Zstat =  @(x,m,s) (x-m)/s;
% Zp =  @(x,m,s) 1-normcdf(Zstat(x,m,s));
% X2stat =  @(x,m,s) 2*m/s^2*x;
% X2df =  @(m,s) 2*m^2/s^2;
% X2p =  @(x,m,s) 1-chi2cdf(X2stat(x,m,s),X2df(m,s));

warning('off','stats:lillietest:OutOfRangePHigh')
warning('off','stats:lillietest:OutOfRangePLow')

for r=1:nRlz
  if ~mod(r,10); disp(num2str(r)); end
  
  
  DY=diff(Y,1,2);
  %DVARS=sqrt(sum(DY.^2)./I);
  DVARS2=mean(DY.^2);

  Rob_S_D = IQRsd(DY')';
  Mn = [ mean(DVARS2) mean(Rob_S_D.^2) median(Rob_S_D.^2) mean(Rob_S_D).^2 ];
  Va = var(DVARS2);
  NrmP  = [];
  for d=[1 1/2 1/3 1/4] 
    if d==1
      Mn = [ Mn median(DVARS2) ];
      Va = [ Va IQRsd(DVARS2)^2 hIQRsd(DVARS2)^2 ];
      [~,p]=lillietest(DVARS2);
      NrmP=[NrmP p];
    else				       
      Z = DVARS2.^d;
      M_Z = median(Z);
      S_Z = IQRsd(Z);
      hS_Z = hIQRsd(Z);

      Va = [Va (1/d*M_Z^(1/d-1)*S_Z)^2  (1/d*M_Z^(1/d-1)*hS_Z)^2 ];

      [~,p]=lillietest(DVARS2.^d);
      NrmP=[NrmP p];
    end
  end

  Mean(r,:) = Mn;
  Var(r,:)  = Va;
  NormP(r,:) =  [NrmP];


end

MnTrue=2*mean(SD.^2);
VaTrue=8*mean(SD.^4)/I;

save(['Sim_DVARS_Inf_t' num2str(T) '_vh' num2str(diff(SDrng)) '.mat'],'NormP','NormPNms','DVARS2','DVARS','nRlz','I','T','SDrng','MeanNms','DVARS','VarNms','MnTrue','VaTrue','Var','Mean')

%% Vis
% 
% fprintf('\n')
% 
% fprintf('\nSettings: nRlz=%d  I=%d  T=%d  SDrng=[%d,%d]\n',nRlz,I,T,SDrng);...
% fprintf('\n  Mean Bias, Percent\n\n');...
% disp(array2table((mean(Mean,1)./MnTrue-1)*100,'VariableNames',MeanNms));...
% fprintf('\n  sqrt(Variance Bias), Percent\n\n');...
% disp(array2table((sqrt(mean(Var,1)./VaTrue)-1)*100,'VariableNames',VarNms));...
% fprintf('\n  Reject Rate\n\n');...
% mnRR=squeeze(mean(RR)); % Dim: length(RRalpha) x length(TestNms)
% disp(array2table(mnRR,'VariableNames',TestNms,'RowNames',RRNms));...
% fprintf('\n  Reject Rate Bias, Percent\n\n');...
% disp(array2table(bsxfun(@times,bsxfun(@minus,mnRR,RRalpha'),1./RRalpha')*100,'VariableNames',TestNms,'RowNames',RRNms));...
% fprintf('\n  Normality P-value (mean)\n\n');...
% disp(array2table(mean(NormP,1),'VariableNames',NormPNms))
% 
% figure(1)
% subplot(1,2,1)
% plot((1:(T-1))/T,squeeze(median(sort(Pval(:,1:2:end,:),3),1))','linewidth',2);
% axis tight
% legend(TestNms(1:2:end),'Location','SouthEast','Interpreter','none')
% %abline(0,1,'linestyle','-','color',[1 1 1]*.8,'linewidth',4)
% title('PP plot, median of all realisations')
% subplot(1,2,2)
% loglog((1:(T-1))/T,squeeze(median(sort(Pval(:,1:2:end,:),3),1))','linewidth',2);
% axis tight
% legend(TestNms(1:2:end),'Location','SouthEast','Interpreter','none')
% %abline(0,1,'linestyle','-','color',[1 1 1]*.8,'linewidth',4)
% title('log PP plot, median of all realisations')
% 
% 
% figure(2)
% subplot(1,2,1)
% h=histogram(DVARS2,50,'Norm','pdf');
% M_DV2 = Mn(3);             % Median of voxelwise sig2
% S_DV2 = sqrt(Va(2));       % half-IQR of DVARS2
% c=X2stat(1,M_DV2,S_DV2);
% nu=X2df(M_DV2,S_DV2);
% xx=get(gca,'xlim');xx=linspace(xx(1),xx(2),100);
% hold on 
% h(2)=plot(xx,chi2pdf( xx   *c,nu)*c,'linewidth',4);
% hold off
% legend(h,{'DVARS','Chi2(s2m_Dh)'},'Interpreter','none')
% title('Emperical and fitted distribution, last realisation')
% 
% subplot(1,2,2)
% DefCol=get(groot,'defaultAxesColorOrder');
% h(1)=histogram(ZStat(:,1,:),100,'Norm','pdf','DisplayStyle','Stairs','LineWidth',2,'edgecolor',DefCol(1,:));
% hold on
% h(2)=histogram(ZStat(:,2,:),100,'Norm','pdf','DisplayStyle','Stairs','LineWidth',2,'edgecolor',DefCol(2,:));
% h(3)=histogram(ZStat(:,3,:),100,'Norm','pdf','DisplayStyle','Stairs','LineWidth',2,'edgecolor',DefCol(3,:));
% xx=get(gca,'xlim');xx=linspace(xx(1),xx(2),100);
% h(4)=plot(xx,normpdf(xx),'linewidth',4,'color',[1 1 1]*.8);
% hold off
% legend(h,[TestNms(2:2:end),{'N(0,1)'}],'Location','SouthEast','Interpreter','none')
% title('Emperical and null Z tests, all realisations')