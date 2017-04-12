function [hndls]=concentricplots_sect_all(x2,y2,VT,figh)

rr=VT(1,2:end);
l_idx=VT(2:end,2:end)';

if ~exist('verbose','var'); verbose=1; end
if verbose
    t3_rown={'Dvar','Svar','Evar'};
    t3_varn={'Global','nonGlobal'};
    disp('----------------------')
    disp('Whole Var Components (Radius of Layers)')
    disp(array2table(rr,'VariableNames',{'Dvar','Svar','Evar'}))
    disp('------------')
    disp('Sum-of-squares (SS) Table (Sectors)')
    disp(array2table(fix(l_idx),'VariableNames',t3_varn,'RowNames',t3_rown))
    disp('----------------------')
end    

l_idx=round(l_idx./repmat(sum(l_idx,2),1,size(l_idx,2))*100); %row: compvar (D/E/S) col: sub compvar (g/ng)

%rr=[6 9 0.5];
%Clk=[.9 .9 .9;.6 .6 .6;.3 .3 .3];
%Clk=[0 0 1;0 1 0;1 0 0];
Clk=[0    0.4470    0.7410;0.9290    0.6940    0.1250;0.4660    0.6740    0.1880];

if exist('figh','var')
    figure(figh);
else
    figure;
end
hold on; axis equal; box on; grid on; 
[hndls]=nested_circ(x2,y2,rr,l_idx,Clk);
legend([hndls{1,1}.Patch hndls{1,2}.Patch],{'Global Var','non-Global Var'})

function [hndls]=nested_circ(xxx,yyy,rr,l_idx,Clk)
rr=fliplr(cumsum(rr)+0.5);
l_idx=flipud(l_idx);

%[rr,idx_rr]=sort(rr+0.5,'descend');
%Clk=Clk(idx_rr,:);
for r_cnt=1:numel(rr)
    %hold on;
    l_idx_l=l_idx(r_cnt,:);
    %xycum=cumsum(l_idx_l);
    %rad0=[0 cumsum(l_idx_l./100)];
    %rad1=cumsum(l_idx_l./100);
    rad0=[0 cumsum(l_idx_l./100)];
    rad1=cumsum(l_idx_l./100);
    
    for l_cnt=1:size(l_idx,2)
        %disp('---')
        %[rad0(l_cnt),rad1(l_cnt)]
        %[x,y]=draw_circ(xxx,yyy,rr(r_cnt),rad0(l_cnt),rad1(l_cnt),Clk(l_cnt,:));
        [~,~,hndls_tmp]=draw_circ(xxx,yyy,rr(r_cnt),rad0(l_cnt),rad1(l_cnt),Clk(l_cnt,:));
        hndls{r_cnt,l_cnt}=hndls_tmp;
        %plot([xxx x(xycum(l_cnt))],[yyy y(xycum(l_cnt))],'color',[1 1 0],'linewidth',1.5);
    end
end
draw_circ(xxx,yyy,0.5,0,1,'w');

function [x,y,hndls]=draw_circ(xx,yy,r,rad0,rad1,Clk)
a1=2*pi*rad0;
a2=2*pi*rad1;
%[a1,a2]
t=linspace(a1,a2);
x=r*cos(t)+xx; 
y=r*sin(t)+yy;

lin_h=plot(x,y,'color','k','linewidth',4);
pat_h=patch([xx,x,xx],[yy,y,yy],Clk,'edgecolor','none');
hndls.Line=lin_h;
hndls.Patch=pat_h;
%fill(x,y)
%patch(x,y,Clk)