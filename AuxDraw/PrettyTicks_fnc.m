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

return