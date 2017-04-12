function Tics = PrettyTicksBetter(Ylim0)

MaxNtick=11;

TickA=[0:10:150 170:20:500];% 600:100:1000];  

TickB=[0:20:150 170:40:500];% 600:100:1000]; 

TicA=TickA(TickA>Ylim0(1) & TickA<Ylim0(2));

TicB=TickB(TickB>Ylim0(1) & TickB<Ylim0(2));

if length(TicA)<MaxNtick
  Tics=TicA;
else
  Tics=TicB;
end