function T=PrettyTicks_fnc(Lim,varargin)
% For a given axis limit, find pretty tick spacing; assumes 50 is always
% in the plot (i.e. that rounded integers are always appropriate)
% Ylim - Y axis limts
MinTick=3;  % Minimum number of tick locations
if ~isempty(varargin)
    TickSp=[10 5 2.5 1 0.5 0.1]./varargin{1};
elseif isempty(varargin)
    TickSp=[10 5 2.5 1 0.5 0.1];
end
ts=0;
T=[];
while length(T)<MinTick
  ts=ts+1;
  if ts>length(TickSp)
    break
  end
  TS=TickSp(ts);
  T = ceil(Lim(1)/TS)*TS : TS : floor(Lim(2)/TS)*TS;
end


% if Lim(2)<=100
%   % Ideally, 0 and 100 will always be in the plot
%   T = [ceil(Lim(1)/10)*10]:10:[floor(Lim(2)/10)*10];
%   if length(T)>6
%     T = [0 25 50 75 100];
%   end
% elseif Lim(2)<=200
%   T = [0 25 50 75 100 150 200];
% elseif Lim(2)<=400
%   T = [0 50 100 200 300 400];
% else
%   T = [0 50 100*round(linspace(1,Lim(2)/100,4))];
% end
% T(T<Lim(1))=[];
% T(T>Lim(2))=[];

% if Lim(2)<=100
%   % Ideally, 0 and 100 will always be in the plot
%   T = [0 25 50 75 100];
% elseif Lim(2)<=200
%   T = [0 25 50 75 100 150 200];
% elseif Lim(2)<=400
%   T = [0 50 100 200 300 400];
% else
%   T = [0 50 100*round(linspace(1,Lim(2)/100,4))];
% end
% T(T<Lim(1))=[];
% T(T>Lim(2))=[];

return