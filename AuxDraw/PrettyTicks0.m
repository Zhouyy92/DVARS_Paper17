function T=PrettyTicks0(Lim)

if Lim(2)<=100
  % Ideally, 0 and 100 will always be in the plot
  T = [ceil(Lim(1)/10)*10]:10:[floor(Lim(2)/10)*10];
  if length(T)>6
    T = [0 25 50 75 100];
  end
elseif Lim(2)<=200
  T = [0 25 50 75 100 150 200];
elseif Lim(2)<=400
  T = [0 50 100 200 300 400];
else
  T = [0 50 100*round(linspace(1,Lim(2)/100,4))];
end
T(T<Lim(1))=[];
T(T>Lim(2))=[];