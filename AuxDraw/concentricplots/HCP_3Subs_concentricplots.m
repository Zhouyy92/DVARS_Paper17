clear

Site='HCP';
modee={'Pre_FIX','Post_FIX'}; %modee: FIX/Prep

sublist={'100307','105115','113619'};

figh_DES=figure; hold on; suptitle('DES variances (circles) - Global NonGlobal (Sectors)')
figh_Glob=figure; hold on; suptitle('Global DES variances')
figh_nGlob=figure; hold on; suptitle('NonGlobal DES variances')
figh_GlobAndNonGlob=figure; hold on; %suptitle('Global and NonGlobal DES variances')
sp_cnt=1;
for mm=modee
    g_Vall=[];
    for s=1:numel(sublist) 
        disp(['Mode:' mm{1} 'Sub: ' sublist{s} ', ' num2str(s)])
        
        load(['/Users/sorooshafyouni/Home/DVARS/R_HCP_3Vols/R_' mm{1} '/DVARSInf_DSEVars_' mm{1} '_' sublist{s} '.mat'],'VT','V')
        
        figure(figh_DES)
        subplot(2,3,sp_cnt);
        hold on; 
        title(['Mode:' mm{1} 'Sub: ' sublist{s}],'Interpreter', 'none')
        %concentricplots_draw(log10(V.w_Avar),log10(V.w_Avar),log10(VT),[1 2],1,figh);
        concentricplots_draw(V.w_Avar,V.w_Avar,VT,[1 2],1,figh_DES);
        hold off;
        
        figure(figh_Glob)
        subplot(2,3,sp_cnt);
        hold on; 
        title(['Mode:' mm{1} 'Sub: ' sublist{s}],'Interpreter', 'none')        
        concentricplots_draw(V.w_Avar,V.w_Avar,VT,1,1,figh_Glob);
        hold off;
        
        figure(figh_nGlob)
        subplot(2,3,sp_cnt);
        hold on; 
        title(['Mode:' mm{1} 'Sub: ' sublist{s}],'Interpreter', 'none')        
        concentricplots_draw(V.w_Avar,V.w_Avar,VT,2,1,figh_nGlob);
        hold off;
        
        figure(figh_GlobAndNonGlob)
        h_sp_1(sp_cnt)=subplot(2,1,1); hold on; axis tight; title('Global Variances')
        concentricplots_draw(2*sp_cnt*10e13,V.g_Avar,VT,1,1,figh_GlobAndNonGlob);
        h_sp_2(sp_cnt)=subplot(2,1,2); hold on; axis tight; title('Non-Global Variances')
        concentricplots_draw(2*sp_cnt*10e13,V.rg_Avar,VT,2,1,figh_GlobAndNonGlob);
        
        yl1(s,:)=ylim(h_sp_1(sp_cnt)); yl2(s,:)=ylim(h_sp_2(sp_cnt));
        
        sp_cnt=sp_cnt+1;
    end
    plot(h_sp_1(sp_cnt-1),[2*(s+0.5)*10e13,2*(s+0.5)*10e13],max(yl1),'k','linewidth',1.8)
    plot(h_sp_2(sp_cnt-1),[2*(s+0.5)*10e13,2*(s+0.5)*10e13],max(yl2),'k','linewidth',1.8)
end
h_sp_1(1).XTick=h_sp_1(1).XTick(1:2:end);
h_sp_1(1).XTickLabel=[strcat('PreFix-',sublist),strcat('PostFix-',sublist)]';
%h_sp_1(1).XTickLabelRotation=90;

h_sp_2(1).XTick=h_sp_2(1).XTick(1:2:end);
h_sp_2(1).XTickLabel=[strcat('PreFix-',sublist),strcat('PostFix-',sublist)]';