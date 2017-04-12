% [hndls]=concentricplots_draw(x2,y2,VT,globalvar_id,verbose,figh)
%
%
%INPUTS:
%   x2, y2: coordinates of the origin
%   VT: *raw* VT matrix from DSEvars
%   globalvar_id: 
%               1:      Global Variances
%               2:      Non-Global Variances
%               [1 2]:  Global and Non-Global Variances
%               
%               Considerations:
%               1. If the first two options (1 or 2) is selected then the
%               layout is as a single concentric circles, with maximume of
%               three sectors. Each sector represent one of the DSE 
%               variance components.
%
%               2. If the last option ([1 2]) is selected, then the layout
%               of the circles are changed. There will be max three
%               concentric circles. For each of the circles, a three-part
%               sector demonstrate the DSE variance components. 
%
%               3. For raw VT matrix, the global variance components and E
%               variance component is far small, therefore, they cannot be
%               visible in the visualisations. Therefore, it is suggested
%               to use log10 of the values. 
%
%
%OUTPUTS:
%   hndls: Handle of all illustrated patches/lines.
%
% Example:
% 
%
%%Sim specs
% I=10e3; T=1000; Y=randn(I,T);
% SDrng=[200 500];
% SD=SDrng(1)+diff(SDrng)*rand(I,1);
% KSD=spdiags(SD,0,I,I);
% Y=KSD*Y;
% 
% [~,V,~,VT]=DSEvars(Y,I,T)
% 
% concentricplots_draw(log10(V.w_Avar),log10(V.w_Avar),log10(VT),[1 2],1);
%
%
%   SA, 2017, UoW
%   srafyouni@gmail.com

function [hndls]=concentricplots_draw(x2,y2,Stat,varargin)

VT=Stat.VT;

%VarNameList={'Global','NonGlobal'};
%RowNameList={'Dvar','Svar','Evar'};

if sum(strcmpi(varargin,'whole'))
   globalvar_id      =  0;
   Clk=[0 0.447 0.741;0.929 0.694 0.125;0.466 0.674 0.188];
   rr=VT(1,1);
   l_idx=VT(globalvar_id+1,2:end);
elseif sum(strcmpi(varargin,'global'))
   globalvar_id      =  1;
   Clk=[0 0.447 0.741;0.929 0.694 0.125;0.466 0.674 0.188];
   rr=VT(1,1);
   l_idx=VT(globalvar_id+1,2:end);
elseif sum(strcmpi(varargin,'noglobal'))
    globalvar_id      =  2;
    Clk=[0.46 0.67 0.18;0.63 0.07 0.18;0.30 0.74 0.93];
    rr=VT(1,2:end);
    l_idx=VT(globalvar_id+1,2:end)';
else
    globalvar_id      =  [1 2];
end
if sum(strcmpi(varargin,'Col'))
   Clk      =   varargin{find(strcmpi(varargin,'Col'))+1};
end
if sum(strcmpi(varargin,'verbose'))
   verbose      =   varargin{find(strcmpi(varargin,'verbose'))+1};
end
if sum(strcmpi(varargin,'figure'))
   figh      =   varargin{find(strcmpi(varargin,'figure'))+1};
   figure(figh);
end
if sum(strcmpi(varargin,'subplot'))
   figh      =   varargin{find(strcmpi(varargin,'subplot'))+1};
   subplot(figh);
end
if sum(strcmpi(varargin,'subplot')) && sum(strcmpi(varargin,'figure'))
   figure;
end

% if numel(globalvar_id)==1
%   Clk=[0 0.447 0.741;0.929 0.694 0.125;0.466 0.674 0.188];
%   rr=VT(1,1);
%   l_idx=VT(globalvar_id+1,2:end);
% elseif numel(globalvar_id)==2
%   Clk=[0.46 0.67 0.18;0.63 0.07 0.18;0.30 0.74 0.93];
%   rr=VT(1,2:end);
%   l_idx=VT(globalvar_id+1,2:end)';
% else
%     error('Wrong dude!')
% end

% if verbose
%     t3_varn=VarNameList(globalvar_id);
%     t3_rown=RowNameList;
%     disp('----------------------')
%     disp('Whole Var Components (Radius of Layers)')
%     %disp(array2table(rr,'VariableNames',{'Dvar','Svar','Evar'}))
%     disp('------------')
%     %if numel(globalvar_id)==2
%         disp('Sum-of-squares (SS) Table (Sectors)')
%         fix(l_idx)
%         %disp(array2table(l_idx,'VariableNames',t3_varn,'RowNames',t3_rown))
%     %end
%     disp('----------------------')
% end    

l_idx = round(l_idx./repmat(sum(l_idx,2),1,size(l_idx,2)) * 100); %row: compvar (D/E/S) col: sub compvar (g/ng)

hold on; axis equal; box on; grid on; 
[hndls] = nested_circ(x2,y2,rr,l_idx,Clk);

function [hndls] = nested_circ(xxx,yyy,rr,l_idx,Clk,draw_rad)
% calls 'draw_circ' appropriately so it explains all the variance comps. 
% SA, 2017, UoW
%
if ~exist('draw_rad','var'); draw_rad=1; end

rr    = fliplr(cumsum(rr)+0.5);
l_idx = flipud(l_idx);

for r_cnt=1:numel(rr) %for each strip...
    
    l_idx_l = l_idx(r_cnt,:);
    rad0    = [0 cumsum(l_idx_l./100)];
    rad1    = cumsum(l_idx_l./100);
    
    for l_cnt=1:size(l_idx,2) %draw a sector
        [~,~,hndls_tmp]    = draw_circ(xxx,yyy,rr(r_cnt),rad0(l_cnt),rad1(l_cnt),Clk(l_cnt,:));
        hndls{r_cnt,l_cnt} = hndls_tmp; clear hndls_tmp;
        %plot([xxx x(xycum(l_cnt))],[yyy y(xycum(l_cnt))],'color',[1 1 0],'linewidth',1.5);
    end
end
if draw_rad
    xtmp = rr(1) * cos(linspace(0,2 * pi)) + xxx; 
    ytmp = rr(1) * sin(linspace(0,2 * pi)) + yyy;
    
    for i = 1:50:numel(xtmp)
        plot([xxx xtmp(i)],[yyy ytmp(i)],'color','k','linewidth',1.2,'linestyle','-.'); 
    end;
    clear xtmp ytmp
end
draw_circ(xxx,yyy,sum(rr)./20,0,1,'w');

function [x,y,hndls]=draw_circ(xx,yy,r,rad0,rad1,Clk)
%draws a sector of a give radiant, on a given radius, with a given fill
%colour and with/without radii grids.
%
% SA, 2017, UoW

%test circles--------
% t = linspace(0, 2*pi);
% r = 1;
% x = r*cos(t);
% y = r*sin(t);
% figure(1)
% patch(x, y, 'g')
% axis equal
%---------------------

a1      = 2 * pi * rad0;
a2      = 2 * pi * rad1;
t       = linspace(a1,a2);
x       = r * cos(t) + xx; 
y       = r * sin(t) + yy;
lin_h   = plot(x,y,'color','k','linewidth',4);
pat_h   = patch([xx,x,xx],[yy,y,yy],Clk,'edgecolor','none');

hndls.Line = lin_h;
hndls.Patch = pat_h;