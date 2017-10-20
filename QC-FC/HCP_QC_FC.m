clear; close all;

Site={'HCP'};

ModeList={'Unproc','Pre_Fix','Post_Fix'};

AtlasList={'Power2011'};

addpath /Users/sorooshafyouni/Home/GitClone/DVARS/Nifti_Util

matpath='/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/Netmats/All_nets/DirRL_Run1/';

%21 Healthy 
SubList_DVARS=[100307,103414,105115,111312,113619,115320,117122,118730,123117,151526,187345,303624,132118,901139,171330,263436,191336,779370,195041,145127,172029];
GSRStatList={'NoGSR'}; %%BE CAREFULLLLL****************

load(['PowerAtlas_EuclideanDistance.mat'],'ED')

m=ModeList(3);
for a=AtlasList
    for g=GSRStatList
        load([matpath 'Netmats_' Site{1} '_' a{1} '_' g{1} '.mat'],'netmat_raw','SubList')
        [SubList1,AI] = intersect(SubList,SubList_DVARS);  
        mat = netmat_raw(:,:,AI); clear netmat_raw;
        nn  = size(mat,1);
        idx=find(triu(ones(nn),1));
        ED=ED(idx);
%         cnt_m = 1;
%         for m = ModeList
            cnt_s=1;
            for s=1:numel(SubList1) 
                smat=mat(:,:,s);
                liemat(:,s)=smat(idx);
                
                load(['/Users/sorooshafyouni/Home/DVARS/fMRIDiag/HCP/R/HCP_' num2str(SubList1(s)) '_fMRIDiag/' Site{1} '_' num2str(SubList1(s)) '_noglobal_' m{1} '_Council.mat'],'FDts','DVARS_Stat_X2_m3d3')
                disp([num2str(SubList1(s)) ' -- ' m{1} ' -- ' g{1}])
                FD(cnt_s)=mean(FDts);
                
                %rdvars(cnt_s)=median(DVARS_Stat_X2_m3d3.Zval);
                
                %rdvars(cnt_s)=median(DVARS_Stat_X2_m3d3.Zval);
                
                cnt_s=cnt_s+1;
            end
            %FDvFC=diag(corr(repmat(FD',1,34716),liemat'));
            
            FDvFC=diag(corr(repmat(FD',1,34716),liemat'));
            
            scatter(ED,FDvFC,'marker','.')
            line([min(ED) max(ED)],[0 0],'color','r','linewidth',1.3)
            ylabel(['QC-FC'],'fontsize',12)
            xlabel(['Inter-node Distance (mm)'],'fontsize',12)
%             cnt_m=cnt_m+1;
%         end
    end
end