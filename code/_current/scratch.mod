reset;
set Sectors ordered, = 1..4;


param temp_sum3 default 0;
param INV_SHR 'the importance of i in jss kapital' {Sectors, Sectors} >= 0,
  default 1;
data;
for {j in 1..card(Sectors)} {
  option randseed j;
  param: temp_sum3 := sum{Sectors} Uniform(.1, .9);
  option randseed j;
  param: INV_SHR[*, j] = Uniform(.1, .9) / temp_sum2;
}
display INV_SHR;
