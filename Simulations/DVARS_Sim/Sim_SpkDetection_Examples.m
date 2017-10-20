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

%----------------------------------------
rho = 0.6;
noiseflag= 1;

for ns_cnt=1:numel(SpksRt_List)
    sp_cnt=1;
    fh0=figure('position',[50,500,500,400]); 
    hold on;
    %suptitle(tl{1})
    for t_cnt=1:numel(T_List)
        T = T_List(t_cnt);
        nspike = round(T*SpksRt_List(ns_cnt));
        for s_cnt=1:numel(SD_List)
            SD = SD_List(s_cnt);  


            load(['/Users/sorooshafyouni/GoogleDrive/FromIDH/DSE_Sim/Code/R/DSE_Sim_' num2str(T) 'l_' num2str(rho*100) 'p_' num2str(SD) 'vh_' num2str(nspike) 'nspk_noisef' num2str(noiseflag) '.mat'])

            sp0=subplot(4,3,sp_cnt);
            hold on; box on; grid on; 
            title(['T= ' num2str(T) '- VH=' num2str(SD)])

            plot(sDeltapDvar(randi(nRlz),:),'Color','b')
            if s_cnt==1
                ylabel('$\Delta \%$ D-var','fontsize',12,'interpreter','latex')
            end

            if t_cnt==numel(T_List)
                xlabel('Scans','fontsize',12,'interpreter','latex')
            end

            ylim([-20 20])
            xlim([0 T-1])

            sp_cnt=sp_cnt+1;
        end
    end
    set(fh0,'color','w')
    export_fig(fh0,['eg/Sim_Examples_' num2str(SpksRt_List(ns_cnt)) 'spkrt_' num2str(rho*100) 'p_' num2str(noiseflag) 'nosflg.pdf'])
end
        

%----------------------------------------