function fMRIDiag_WholeVar(V,varargin)

FaceAlpha=0.2;

Col=get(groot,'defaultAxesColorOrder');
Acol=Col(5,:); % Green
%FDcol=Col(2,:); % Red (/orange!)
FDcol=[.5 .5 .5];
Dcol=Col(1,:); % Blue
Scol=Col(3,:); % Yellow
Ecol=Col(4,:); % Purple

Idxs = []; Idxp = []; lw = 2;  lfs = 12;
TickScaler = 1; Thickness = 1; 
%
IdxpCol     = Col(2,:);
IdxsCol     = [.5 .5 .5];
%
if sum(strcmpi(varargin,'subplot'))
   hndl  = varargin{find(strcmpi(varargin,'subplot'))+1};
   subplot(hndl)
elseif sum(strcmpi(varargin,'figure'))
   hndl  = varargin{find(strcmpi(varargin,'figure'))+1};
   figure(hndl)
else
    hndl=figure('position',[50,500,1600,1400]); 
    hold on; box on;
    figure(hndl)
end
%
if sum(strcmpi(varargin,'Stat_Idx'))
   Idxs     =   varargin{find(strcmpi(varargin,'Stat_Idx'))+1};
end
if sum(strcmpi(varargin,'Pract_Idx'))
   Idxp     =   varargin{find(strcmpi(varargin,'Pract_Idx'))+1};
end
if sum(strcmpi(varargin,'TickScaler'))
   TickScaler          =   varargin{find(strcmpi(varargin,'TickScaler'))+1};
end

if sum(strcmpi(varargin,'linewidth'))
   lw      =   varargin{find(strcmpi(varargin,'linewidth'))+1};   
end

if sum(strcmpi(varargin,'fontsize'))
   lfs     =   varargin{find(strcmpi(varargin,'fontsize'))+1};
end
if sum(strcmpi(varargin,'Thick'))
   Thickness      =   varargin{find(strcmpi(varargin,'Thick'))+1};
end

T=length(V.Avar_ts);
Time=1:T;

hTime=(1:(T-1))+0.5;
eTime=[1 T];

% hndl=subplot(1,1,1); 
% hold on; box on;

%title('DSE Variance Decomposition (RMS units)','fontsize',13)
yyaxis(hndl,'left')
    line(Time,sqrt(V.Avar_ts),'LineStyle','-','linewidth',lw,'color',Acol)
    line(Time,ones(1,T).*mean(sqrt(V.Avar_ts)),'LineStyle',':','linewidth',.5,'color',Acol)
    
    line(hTime,sqrt(V.Dvar_ts),'LineStyle','-','linewidth',lw,'color',Dcol)
    line(hTime,ones(1,T-1).*mean(sqrt(V.Dvar_ts)),'LineStyle',':','linewidth',.5,'color',Dcol)
    
    line(hTime,sqrt(V.Svar_ts),'LineStyle','-','linewidth',lw,'color',Scol)
    line(hTime,ones(1,T-1).*mean(sqrt(V.Svar_ts)),'LineStyle',':','linewidth',.5,'color',Scol)
    
    line(eTime,sqrt(V.Evar_ts),'LineStyle','none','Marker','o','markerfacecolor',Ecol,'linewidth',3,'color',Ecol)
    ylabel('$\sqrt{\mathrm{Variability}}$','fontsize',lfs,'interpreter','latex')  
    YLim2=ylim.^2/mean(V.Avar_ts)*100;
    set(hndl,'ycolor','k','xlim',[1 T])
yyaxis(hndl,'right')
    YTick2=PrettyTicks(YLim2,TickScaler); YTick=sqrt(YTick2);
    set(hndl,'Ylim',sqrt(YLim2),'YTick',sqrt(YTick2),'YtickLabel',num2str([YTick2']));
    ylabel('\% of A-var','fontsize',lfs,'interpreter','latex')
    %set(sph2,'ygrid','on')
    h=abline('h',YTick);
    set(h,'linestyle','-','color',[.5 .5 .5]); %the grids!
    set(hndl,'ycolor','k','xlim',[1 T])
    
    PatchMeUp(setdiff(Idxs,Idxp),Thickness,IdxsCol,FaceAlpha);
    PatchMeUp(Idxp,Thickness,IdxpCol,FaceAlpha);
    
    
% function T=PrettyTicks(Lim,varargin)
% % For a given axis limit, find pretty tick spacing; assumes 50 is always
% % in the plot (i.e. that rounded integers are always appropriate)
% % Ylim - Y axis limts
% %
% %   TEN & SA, 2017, UoW
% %   srafyouni@gmail.com
% %
% MinTick=3;  % Minimum number of tick locations
% if ~isempty(varargin)
%     TickSp=[15 5 2.5 1 0.5 0.2]./varargin{1};
% elseif isempty(varargin)
%     TickSp=[15 5 2.5 1 0.5 0.2];
% end
% ts=0;
% T=[];
% while length(T)<MinTick
%   ts=ts+1;
%   if ts>length(TickSp)
%     break
%   end
%   TS=TickSp(ts);
%   T = ceil(Lim(1)/TS)*TS : TS : floor(Lim(2)/TS)*TS;
% end
% 
% return
%###################################################################################
function h=abline(a,b,varargin)
% FORMAT h = abline(a,b,...)
% Plots y=a+b*x in dotted line
% FORMAT h = abline('h',y,...)
% Plots a horizontal line at y; y can be a vector, & then multiple lines plotted
% FORMAT h = abline('v',x,...)
% Plots a vertical line at x; x can be a vector, & then multiple lines plotted
%
% ...  Other graphics options, e.g. "'LineStyle','-'" or "'LineWidth',2" or
%      "'color',[1 0 0]",  etc
%
% Like Splus' abline.  Line is plotted and then moved behind all other
% points on the graph.
%
% $Id: abline.m,v 1.1 2013/06/04 10:38:11 nichols Exp $

if (nargin==2) && ischar(a)
  a = lower(a);
else

  if (nargin<1)
    a = 0;
  end
  if (nargin<2)
    b = 0;
  end
  
end

XX=get(gca,'Xlim');
YY=get(gca,'Ylim');
h_exist = get(gca,'children');

g = [];
if ischar(a) && (a=='h')

  for i=1:length(b)
    g=[g;line(XX,[b(i) b(i)],'LineStyle',':',varargin{:})];
  end

elseif ischar(a) && (a=='v')

  for i=1:length(b)
    g=[g;line([b(i) b(i)],YY,'LineStyle',':',varargin{:})];
  end

else

  g=line(XX,a+b*XX,'LineStyle',':',varargin{:});

end

uistack(h_exist,'top');

if (nargout>0)
  h=g;
end

set(gcf,'color','w');
return
%###################################################################################
% function ph=PatchMeUp(Idx,varargin)
% %   Draws a patch across the significantly identified scans on var plots
% %
% %   SA, 2017, UoW
% %   srafyouni@gmail.com
% if nargin == 1
%     stpjmp=1;
% elseif nargin == 2
%     stpjmp=varargin{1};
% end
% 
% yyll=ylim;
% for ii=1:numel(Idx)
%     xtmp=[Idx(ii)-stpjmp   Idx(ii)-stpjmp   Idx(ii)+stpjmp  Idx(ii)+stpjmp];
%     ytmp=[yyll(1)               yyll(2)         yyll(2)        yyll(1)    ];
%     ph(ii)=patch(xtmp,ytmp,[.5 .5 .5],'FaceAlpha',0.3,'edgecolor','none');
%     clear *tmp
% end
% return 
%###################################################################################
function gsrY=fcn_GSR(Y)
%Global Signal Regression
%Inspired from FSLnets
%For the fMRIDiag, it needs to be transposed. 
%
%   SA, 2017, UoW
%   srafyouni@gmail.com
%
Y=Y';
mgrot=mean(Y,2); 
gsrY=Y-(mgrot*(pinv(mgrot)*Y));
gsrY=gsrY';    