function ph=PatchMeUp(Idx)
stpjmp=1;
yyll=ylim;
for ii=1:numel(Idx)
    xtmp=[Idx(ii)-stpjmp   Idx(ii)-stpjmp   Idx(ii)+stpjmp  Idx(ii)+stpjmp];
    ytmp=[yyll(1)               yyll(2)         yyll(2)        yyll(1)    ];
    ph(ii)=patch(xtmp,ytmp,[.5 .5 .5],'FaceAlpha',0.3,'edgecolor','none');
    clear *tmp
end
return 