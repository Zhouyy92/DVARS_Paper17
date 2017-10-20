%function ph=PatchMeUp(Idx,varargin)
function PatchMeUp0(Idxs,Idxp,stpjmp,sLcol,pLcol,faceAlpha,ax_hndl)
% setdiff(Idxs,Idxp),Idxp,PtchWdth,trans,sph5
% ph=PatchMeUp(Idx,varargin)
% patch around the significantly detected spikes.
%
% SA, NISOx, 2017

%Pos=get(ax_hndl,'position');
%set(ax_hndl,'Color','none')
%PatchAx=axes('position',Pos,'visible','on','Ytick',[],'Xtick',[]);
yLm=ylim;

%ph0=[];%ED
for ii=1:numel(Idxs)
     xtmp=[Idxs(ii)-stpjmp   Idxs(ii)-stpjmp   Idxs(ii)+stpjmp  Idxs(ii)+stpjmp];
     ytmp=[yLm(1)               yLm(2)         yLm(2)        yLm(1)    ];
     patch(xtmp,ytmp,sLcol,'edgecolor','none','facealpha',faceAlpha);
     clear *tmp
end

%ph1=[];%ED
for ii=1:numel(Idxp)
     xtmp=[Idxp(ii)-stpjmp   Idxp(ii)-stpjmp   Idxp(ii)+stpjmp  Idxp(ii)+stpjmp];
     ytmp=[yLm(1)               yLm(2)         yLm(2)        yLm(1)    ];
     patch(xtmp,ytmp,pLcol,'edgecolor','none','facealpha',faceAlpha);
     clear *tmp
end

%uistack(PatchAx,'bottom')

%FOR SIMPLE LINES, USE BELOW LINES:
%faceAlpha=1;
%for ii=1:numel(Idx)
%    ph=line([Idx(ii) Idx(ii)],[yyll(1) yyll(2)],'color',[Lcol faceAlpha]);
%    %uistack(ph,'bottom')
%    clear *tmp
%end

return