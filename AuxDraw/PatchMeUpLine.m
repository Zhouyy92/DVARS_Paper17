%function ph=PatchMeUp(Idx,varargin)
function [ph00,ph11]=PatchMeUpLine(Idxs,Idxp,stpjmp,sLcol,pLcol,faceAlpha,ax_hndl)
% setdiff(Idxs,Idxp),Idxp,PtchWdth,trans,sph5
% ph=PatchMeUp(Idx,varargin)
% patch around the significantly detected spikes.
%
% SA, NISOx, 2017

%Pos=get(ax_hndl,'position');
%set(ax_hndl,'Color','none')
%PatchAx=axes('position',Pos,'visible','on','Ytick',[],'Xtick',[]);
yLm=ylim;

subplot(ax_hndl)

ph00=[];%ED
for ii=1:numel(Idxs)
     %xtmp=[Idxs(ii)-stpjmp   Idxs(ii)-stpjmp   Idxs(ii)+stpjmp  Idxs(ii)+stpjmp];
     %ytmp=[yLm(1)               yLm(2)         yLm(2)        yLm(1)    ];
     %patch(xtmp,ytmp,sLcol,'edgecolor','none','facealpha',faceAlpha);
     ph=line([Idxs(ii) Idxs(ii)],[yLm(1) yLm(2)],'color',[sLcol faceAlpha],'linewidth',stpjmp);
     ph00=[ph00 ph];
     clear *tmp
end
%uistack(ph00,'down');

ph11=[];%ED
for ii=1:numel(Idxp)
     %xtmp=[Idxp(ii)-stpjmp   Idxp(ii)-stpjmp   Idxp(ii)+stpjmp  Idxp(ii)+stpjmp];
     %ytmp=[yLm(1)               yLm(2)         yLm(2)        yLm(1)    ];
     %patch(xtmp,ytmp,pLcol,'edgecolor','none','facealpha',faceAlpha);
     ph=line([Idxp(ii) Idxp(ii)],[yLm(1) yLm(2)],'color',[pLcol faceAlpha],'linewidth',stpjmp);
     ph11=[ph11 ph];
     clear *tmp ph
end
%uistack(ph11,'down');

%uistack(PatchAx,'bottom')

%FOR SIMPLE LINES, USE BELOW LINES:
% faceAlpha=1;
% for ii=1:numel(Idxs)
%    ph=line([Idxs(ii) Idxs(ii)],[yLm(1) yLm(2)],'color',[sLcol faceAlpha]);
%    %uistack(ph,'bottom')
%    clear *tmp
% end

return