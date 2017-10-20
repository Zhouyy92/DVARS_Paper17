# This function calls a mediator internal function Sub_Sim_DVARS.sh with time series and level of variance het. 
# 
# SA, UoW, 2017
# srafyouni@gmail.com
#
####REFERENCES
#
#   Afyouni S. & Nichols T.E., Insights and inference for DVARS, 2017
#   http://www.biorxiv.org/content/early/2017/04/06/125021.1

nRLz=1000;

wallclk=(1 2 5 12)
Mem=(3 3 5 9)

for t_cnt in {0..3}
do
  	for sd_cnt in {0..2}
        do
                sh Sub_Sim_DVARS.sh $((t_cnt+1)) $((sd_cnt+1)) ${wallclk[t_cnt]} ${nRLz} ${Mem[t_cnt]}
               echo $((t_cnt+1)) - $((sd_cnt+1)) - ${wallclk[t_cnt]}
        done
done
