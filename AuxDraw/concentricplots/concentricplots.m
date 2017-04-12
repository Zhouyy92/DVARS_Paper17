function concentricplots(x2,y2)

l_idx=randi(100,3);
l_idx=round(l_idx./repmat(sum(l_idx,2),1,3)*100) %row: compvar (D/E/S) col: sub compvar (g/ng)

find(l_idx)

vv=[6 9 0.5];

Clk=[.9 .9 .9;.6 .6 .6;.3 .3 .3];
%Clk=[0 0 1;1 0 0;0 1 0];
figure; hold on; axis equal; box on; grid on; 
nested_circ(x2,y2,vv,l_idx,Clk)

function nested_circ(xxx,yyy,rr,l_idx,Clk)
rr=fliplr(cumsum(rr)+0.5);
l_idx=flipud(l_idx);
rr
%[rr,idx_rr]=sort(rr+0.5,'descend');
%Clk=Clk(idx_rr,:);

for r_cnt=1:numel(rr)
    hold on;
    [x,y]=draw_circ(xxx,yyy,rr(r_cnt),1,Clk(r_cnt,:));
    xycum=cumsum(l_idx(r_cnt,:));
    for l_cnt=1:size(l_idx,2)
        plot([xxx x(xycum(l_cnt))],[yyy y(xycum(l_cnt))],'color',[1 1 0],'linewidth',1.5)
    end
end
draw_circ(xxx,yyy,0.5,1,'w');

function [x,y]=draw_circ(xx,yy,r,rad,Clk)
t=linspace(0, 2*pi./rad);
x=r*cos(t)+xx; 
y=r*sin(t)+yy;

plot(x,y,'color','k','linewidth',1.1)
patch(x,y,Clk)
