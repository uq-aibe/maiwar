/*=============================================================================
This is NAIU (New Australian Inter* Model with Uncertainty, where * stands for
Sectoral, Time, Age and Regional). It is an extension of
Cai--Judd (2021) to include multiple sectors among among other things. If any
part of this code is re-used, please cite OCallaghan (2022).
=============================================================================*/

/*=============================================================================
THE MODEL starts here
=============================================================================*/
/*=============================================================================
Parameters for generating sets
=============================================================================*/ 
param TInf 'infimum of the set of times' default 0;
param LInf 'start period/infimum of the look forward set', default 0;
param LSup 'end period/supremum of the look forward set' > LInf, default 3e+1;
param PInf 'infimum of times on a path', default 0;
param PSup 'supremum of times on a path (T_star in CJ)'# eg 2050 - 2022 = 28 
  default 2 >= PInf; 
param TSup 'final time' >= PInf + LInf, default PSup + LSup;
/*=============================================================================
Sets
=============================================================================*/ 
set Regions; 
set Sectors ordered;
#-----------the planning horizon and corresponding set:
#-----------lower values of LFwd are more myopic; longer satisfy Euler Eqn.
set LookForward 'planning set' = {LInf .. LSup + LInf - 1} ordered;
set LookForwardClosure 'includes extra time for closing out the problem'
  = LookForward union {LSup + LInf} ordered;
#-----------Each path is a state of the world and has the following structure:
set PathTimes 'times along a path' = {PInf .. PSup} ordered;
#-----------The set of paths are indexed by
set PathIndices 'path space'
  # for default we adopt a two-state Markov chain with unique absorbing state
    default {PInf .. PSup + 1};  
set AllTimes 'all time periods associated with a path' = {TInf .. TSup};                      
/*=============================================================================
Parameters
=============================================================================*/
param StrtTm 'start time for each step on each path' {PathTimes, PathIndices}
  default time();
param EndTm 'start time for each step on each path' {PathTimes, PathIndices}
  default time();
param BETA 'discount factor', in (0, 1), default .98; 
param ALPHA 'importance of capital in production' {Sectors} default 33e-2 >= 0;
param DELTA 'rate of depreciation for kapital' {Sectors} default .025 >= 0;
param PHI_ADJ 'kapital adjustment cost' {Sectors} default 5e-1 >= 0;
param GAMMA 'intertemporal elasticity of subst.' {Regions} default 5e-1 >= 0;
param ETA 'Frisch elasticity of labour supply'{Regions} default 5e-1 >= 0;
param EoS_KAP 'elasticity of substitution for kapital' default 1e-1 > 0; 
param INV_MIN 'lower bound of investment' default 0 >= 0;
#param kmin 'smallest capital' default 1e-1 >= 0;
#param kmax 'largest capital' default 1e+1 >= 0;
param ZETA1 'TFP before shock' default 1 >= 0;
param ZETA2 'TFP after shock' default .95 >= 0;
param PROB2 'one period probability of jump in TFP' default 0.01 >= 0;
param TAIL_CON_SHR 'tail consumption share (of output)' default 0.75 >= 0;
#
param UInf 'infimum of interval for uniform dbn' default 0.49 in [0, .5);
param USup 'supremum of interval for uniform dbn' default 0.51 in (.5, 1];
param VInf 'infimum of interval for basic variables' default 1e-4;
param VSup 'supremum of interval for basic variables' default 1e+4;
param OInf 'infimum of interval for observed/actual values' default 1e-4;
param OSup 'supremum of interval for observed/actual values' default 1e+4;
#
param RAW_CON_FLW {Regions, Sectors} default Uniform(UInf, USup) >= 0;
param RAW_LAB_FLW {Regions, Sectors} default Uniform(UInf, USup) >= 0;
param RAW_INV_FLW {Regions, Sectors, Sectors} 
  default Uniform(UInf, USup) >= 0;
param CON 'observed consumption' {Regions, Sectors, PathTimes}
  default 1e+0 in (OInf, OSup);
param INV 'observed investment' {Regions, Sectors, Sectors, PathTimes}
  default 1e+0 in (OInf, OSup);
param INV_SEC 'observed investment' {Regions, Sectors, PathTimes}
  default 1e+0 in (OInf, OSup);
param INV_SUM 'observed total investment' {Regions, Sectors, PathTimes}
  default 1e+0 in (OInf, OSup); 
param LAB 'observed labour' {Regions, Sectors, PathTimes}
  default 1e+0 in (OInf, OSup);
param KAP 'observed kapital' {Regions, Sectors, PathTimes}
  default 1e+0 in (OInf, OSup);
param E_OUTPUT 'observed Exp. output' {Regions, Sectors, PathTimes}
  default 1e+0 in (OInf, OSup); 
param ADJ_COST_KAP 'observed adjustment costs of kapital'
  {Regions, Sectors, PathTimes} default 0 in [0, OSup);
param MKT_CLR 'observed output' {Sectors, PathTimes} default 0 in (-1e-4, 1e-4); 
param JAC_ID 'observed jacobi identity' {Regions, Sectors, Sectors, PathTimes}
  default 0 in (-1e-1, 1e-1);
/*=============================================================================
Computed parameters
=============================================================================*/
param GAMMA_HAT 'utility parameter' {r in Regions}, = 1 - GAMMA[r];
param ETA_HAT 'utility parameter'{r in Regions}, = 1 + ETA[r];
param RHO 'exponent of the ces function', = (EoS_KAP - 1) / EoS_KAP; 
param RHO_INV 'inverse of RHO', = 1 / RHO;
param A 'productivity trend' {i in Sectors}
  default (1 - (1 - DELTA[i]) * BETA) / (ALPHA[i] * BETA) >= 0;
param B 'relative weight of consumption and leisure in utility' 
  {r in Regions, i in Sectors}
    default 1 >= 0;
    # (1 - ALPHA[i]) * A[i] * (A[i] - DELTA[i]) ^ (-1 / GAMMA[r]) >= 0;
param  REG_WGHT 'regional (population) weights' {r in Regions}
  default 1 / card(Regions);

#-----------set the seed for the random number generator for weights
option randseed 12345;
#let USup := .6;
#-----------Share parameters in the aggregators for utility and production.
#-----------For default values, we draw from the uniform distribution.
#-----------consumption
param CON_FLW_SUM {r in Regions} = sum{i in Sectors} RAW_CON_FLW[r, i];
param CON_SHR 'consumption weights for each good in utility'
  {r in Regions, i in Sectors} = RAW_CON_FLW[r, i] / CON_FLW_SUM[r];
param CON_SHR_CES 'consumption weights for each good in utility'
  {r in Regions, i in Sectors} = RAW_CON_FLW[r, i] / CON_FLW_SUM[r] ^ EoS_KAP;
#-----------labour
param LAB_FLW_SUM {r in Regions} = sum{i in Sectors} RAW_LAB_FLW[r, i];
param LAB_SHR 'labour weights for each sector in utility'
  {r in Regions, j in Sectors} = RAW_LAB_FLW[r, j] / LAB_FLW_SUM[r];
#-----------investment (inputs of sector i to j's kapital).
param INV_FLW_RSUM {r in Regions, j in Sectors} 
  = sum{i in Sectors} RAW_INV_FLW[r, i, j];
param INV_SHR "the importance of i in j's kapital"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_INV_FLW[r, i, j] / INV_FLW_RSUM[r, j]);
param INV_SHR_CES "the importance of i in j's kapital"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_INV_FLW[r, i, j] / INV_FLW_RSUM[r, j]) ^ EoS_KAP;
param Pr_shk 'probability of SHK'
  {Regions, Sectors, t in LookForward} = (1 - PROB2) ^ t;
param E_shk 'expected shock (exogenous)'
  {r in Regions, i in Sectors, t in LookForward}
    = ZETA2 + Pr_shk[r, i, t] * (ZETA1 - ZETA2);

/*=============================================================================
Basic (economic) variables
=============================================================================*/
var inv 'investment flows'
  {r in Regions, Sectors, j in Sectors, LookForward, s in PathTimes}
  in [VInf, VSup] default DELTA[j] * KAP[r, j, s];
# inv 'investment flows'
#  {r in Regions, i in Sectors, j in Sectors, t in LookForward, s in PathTimes}
#    in [VInf, VSup]
#    if (i = j or j = 'PbSc') then {default DELTA[j] * KAP[r, j, s];}
#    else { = (inv[r, i, 'PbSc', t, s] / INV_SHR[r, i, 'PbSc'])
#      / (inv[r, i, 'PbSc', t, s]) / INV_SHR[r, i, 'PbSc'])
#      * (inv[r, 'PbSc', j, t, s] / INV_SHR[r, 'PbSc', j])
#      * INV_SHR[r, i, j, t, s];
#    };
var con 'consumption flows' {Regions, Sectors, LookForward, PathTimes}
  in [VInf, VSup] default 1e-0;
var lab 'labour flows' {Regions, Sectors, LookForward, PathTimes}
  in [VInf, VSup] default 1e-0;
#-----------kapital, the dynamic variable, is defined on LookForwardClosure
var kap 'kapital stocks (dynamic: defined on LookForwardClosure)'
  {r in Regions, j in Sectors, LookForwardClosure, s in  PathTimes}
  in [VInf, VSup] default KAP[r, j, s];
/*=============================================================================
Potential intermediate variables (substituted out during pre-solving)
=============================================================================*/
var con_sec_CD 'Cobb-Douglas consumption aggregate (across sectors)'
  {r in Regions, t in LookForward, s in PathTimes}
  = prod{i in Sectors} con[r, i, t, s] ^ CON_SHR[r, i];
var con_sec_CES 'Cobb-Douglas consumption aggregate (across sectors)'
  {r in Regions, t in LookForward, s in PathTimes}
  = (sum{i in Sectors} CON_SHR_CES[r, i] * con[r, i, t, s] ^ RHO) ^ RHO_INV;
var con_sec_SumPow 'Sum of power consumption aggregate (across sectors)'
  {r in Regions, t in LookForward, s in PathTimes}
  = sum{i in Sectors}
    CON_SHR[r, i] * con[r, i, t, s] ^ GAMMA_HAT[r] / GAMMA_HAT[r];
var con_sec_SumShr 'Sum of fractional powers from consumption shares'
  {r in Regions, t in LookForward, s in PathTimes}
  = sum{i in Sectors} con[r, i, t, s] ^ CON_SHR[r, i];
var inv_sec_CES 'Const. Elast. Subst. investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward, s in PathTimes}
  = (sum{i in Sectors} INV_SHR_CES[r, i, j] * inv[r, i, j, t, s] ^ RHO)
    ^ RHO_INV;
var inv_sec_CD 'Cobb-Douglas investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward, s in PathTimes}
  = prod{i in Sectors} (inv[r, i, j, t, s] ^ INV_SHR[r, i, j]);
var lab_sec_CD 'Cobb-Douglas labour aggregate (across sectors)'
  {r in Regions, t in LookForward, s in PathTimes}
  = prod{j in Sectors} B[r, j] * lab[r, j , t, s] ^ LAB_SHR[r, j];
var lab_sec_Q 'quadratic labour aggregate (across sectors)'
  {r in Regions, t in LookForward, s in PathTimes}
  = sum{j in Sectors} (lab[r, j , t, s] ^ 2);
var lab_sec_F 'Frisch labour aggregate (across sectors)'
  {r in Regions, t in LookForward, s in PathTimes}
  = sum{j in Sectors} B[r, j] * lab[r, j , t, s] ^ ETA_HAT[r] / ETA_HAT[r];
var adj_cost_kap_Q 'quadratic adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward, s in PathTimes}
  = PHI_ADJ[i] * kap[r, i, t, s]
    * (kap[r, i, t + 1, s] / kap[r, i, t, s] - 1) ^ 2;
var E_output_CD 'Cobb-Douglas output transformation'
  {r in Regions, i in Sectors, t in LookForward, s in PathTimes}
  = E_shk[r, i, t] 
    * kap[r, i, t, s] ^ ALPHA[i] * lab[r, i, t, s] ^ (1 - ALPHA[i]);
var utility_CD 'Cobb-Douglas instantaneous utility'
  {t in LookForward, s in PathTimes}
  = sum{r in Regions}
    REG_WGHT[r] * (con_sec_CD[r, t, s] - lab_sec_CD[r, t, s]);
var utility_logCD_Q 'log Cobb-Douglas and Quadratic instantaneous utility'
  {t in LookForward, s in PathTimes}
  = sum{r in Regions}
    REG_WGHT[r] * (log(con_sec_CD[r, t, s]) - lab_sec_Q[r, t, s]);
var utility_CES_Q 'Const. Elast. Subst. and Quadratic instantaneous utility'
  {t in LookForward, s in PathTimes}
  = sum{r in Regions}
    REG_WGHT[r] * (con_sec_CD[r, t, s] - lab_sec_Q[r, t, s]);
var utility_CD_F 'Cobb-Douglas and Frisch instantaneous utility'
  {t in LookForward, s in PathTimes}
  = sum{r in Regions}
    REG_WGHT[r] * (con_sec_CD[r, t, s] - lab_sec_F[r, t, s]);
var utility_SumShr_Q 'utility: SumShr for consumption, quadratic for labour'
  {t in LookForward, s in PathTimes}
  = sum{r in Regions}
    REG_WGHT[r] * (con_sec_SumShr[r, t, s] - lab_sec_Q[r, t, s]);
var utility_SumPow_Q'utility: SumPow for consumption and quadratic for labour'
  {t in LookForward, s in PathTimes}
  = sum{r in Regions}
    REG_WGHT[r] * (con_sec_SumPow[r, t, s] - lab_sec_Q[r, t, s]);
var tail_val_logCD 'SumShr continuation value from time LSup + LInf onwards'
  {s in PathTimes}
  = sum{r in Regions} REG_WGHT[r] * (sum{i in Sectors} CON_SHR[r, i] * (
      TAIL_CON_SHR * 5e+0 * ALPHA[i] * log( kap[r, i, LSup + LInf, s])
      ) / (1 - BETA) - .4);
/*=============================================================================
Current intermediate variables (substituted out during pre-solving)
=============================================================================*/
var E_output 'current intermediate variable for output'
  {r in Regions, i in Sectors, t in LookForward, s in PathTimes}
    = E_output_CD[r, i, t, s] * 5e+0;
var inv_sec 'current intermediate variable for aggregated investment'
  {r in Regions, j in Sectors, t in LookForward, s in PathTimes}
    = inv_sec_CES[r, j, t, s];
var adj_cost_kap 'current adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward, s in PathTimes}
    = adj_cost_kap_Q[r, i, t, s];
var utility 'current intermediate variable for utility'
  {t in LookForward, s in PathTimes} = utility_logCD_Q[t, s];
var tail_val 'current intermediate variable for tail value function'
  {s in PathTimes} = tail_val_logCD[s] * 6e-2;# for LSup=2e+2, tail_param = 3.789; 
/*=============================================================================
The objectives and constraints
=============================================================================*/
maximize pres_disc_val 'present discounted value of utilities'
  {s in PathTimes}:
    sum{t in LookForward} BETA ^ (t - LInf) * utility[t, s]
      + BETA ^ (LSup - LInf) * tail_val[s];
subject to accum_kap_eq 'accumulation of kapital'
  {r in Regions, j in Sectors, t in LookForward, s in PathTimes}:
    kap[r, j, t + 1, s]
      = (1 - DELTA[j]) * kap[r, j, t, s] + inv_sec[r, j, t, s];
subject to market_clearing_eq 'market clearing for each sector and time'
  {i in Sectors, t in LookForward, s in PathTimes}:
    sum{r in Regions}(
      con[r, i, t, s]
      + sum{j in Sectors}(inv[r, i, j, t, s])
      + adj_cost_kap[r, i, t, s]
      - E_output[r, i , t, s]
      ) <= 0;
/*=============================================================================
The_data
=============================================================================*/
update data Regions, Sectors;
data;
#-----------2x2 model
#set Regions := SEQ RoQ;
#set Sectors := Agrc PbSc;
#-----------2x3 model
#set Regions := SEQ RoQ;
#set Sectors := Agrc PbSc Mnfc;
#-----------2x4 model
#set Regions := SEQ RoQ;
#set Sectors := Agrc PbSc Mnfc Srvc;
#-----------3x4 model
#set Regions := SEQ RoQ RoA;
#set Sectors := Agrc PbSc Mnfc Srvc;
#-----------3x5 model
#set Regions := SEQ RoQ RoA;
#set Sectors := Agrc PbSc Mnfc Srvc Trns;
#-----------3x6 model
#set Regions := SEQ RoQ RoA;
#set Sectors := Agrc PbSc Mnfc Srvc Trns Utlt;
##-----------4x4
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := Agrc PbSc Mnfc Srvc;# Trns Utlt;
##-----------4x5
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := Agrc PbSc Mnfc Srvc Utlt;
##-----------4x6
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := Agrc PbSc Frst Mnfc Srvc Trns Utlt;
##-----------10x8
#set Regions := CWQ MIW DrD Ftz FNQ NWQ SEQ SWQ NrQ WBB;
#set Sectors := Agrc PbSc Elct Frst Mnfc Srvc Trns WtSw;
##-----------7x10
#set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA;# MIW DrD Ftz FNQ NWQ;# SWQ NrQ WBB;
#set Sectors := F G H I J K L M N PbSc;
##-----------3x10
set Regions := CnQ FNQ Mck;# MIW DrD Ftz FNQ NWQ;# SWQ NrQ WBB;
set Sectors := F G H I J K L M N PbSc;
##-----------2x11
#set Regions := CnQ FNQ;# MIW DrD Ftz FNQ NWQ;# SWQ NrQ WBB;
#set Sectors := E F G H I J K L M N PbSc;
##-----------3x11
#set Regions := CnQ FNQ Mck;# MIW DrD Ftz FNQ NWQ;# SWQ NrQ WBB;
#set Sectors := E F G H I J K L M N PbSc;
##-----------8x16
#set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA RoW;# MIW DrD Ftz FNQ NWQ;# SWQ NrQ WBB;
#set Sectors := A B C D E F G H I J K L M N PbSc P;
##-----------7x20
#set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA;# MIW DrD Ftz FNQ NWQ;# SWQ NrQ WBB;
#set Sectors := A B C D E F G H I J K L M N PbSc P Q R T U;
#-----------display some parameter values:
display  CON_SHR, LAB_SHR, INV_SHR;
/*=============================================================================
Solving the model
=============================================================================*/
param InstanceName symbolic;
option substout 1;
option solver conopt;
#option solver knitro;
#option solver baron;
#option baron_options trace;
option show_stats 1;
#-----------start the algorithm for iterating along a given path
for {s in PathTimes}{
  if s = PInf then {
    fix {r in Regions, j in Sectors}
      kap[r, j, LInf, s] := KAP[r, j, s];
#-----------set and solve the plan for start time s
    objective pres_disc_val[s];
    #let InstanceName := ("~/Documents/_uq-aibe/maiwar-solutions/"
    #  & time() & "maiwar" & card(Regions) & "x" & card(Sectors)
    #  & "x" & card(LookForward) & "x" & card(PathTimes) & "s" & s & "log");
    #write ("b" & InstanceName);
    solve;
    #solution (InstanceName & ".sol");
    display _ampl_elapsed_time, _total_solve_time;
  }
  else {
    display s;
    for {r in Regions, i in Sectors}{
#-----------save actual path values of variables to parameter
      let CON[r, i, s - 1] := con[r, i, LInf, s - 1];
      let INV_SEC[r, i, s - 1] := inv_sec[r, i, LInf, s - 1];
      let INV_SUM[r, i, s - 1] := sum{j in Sectors} inv[r, i, j, LInf, s - 1];
      let LAB[r, i, s - 1] := lab[r, i, LInf, s - 1];
      let E_OUTPUT[r, i, s - 1] := E_output[r, i, LInf, s - 1];
      let ADJ_COST_KAP[r, i, s - 1] := adj_cost_kap[r, i, LInf, s - 1];
      let KAP[r, i, s] := kap[r, i, 1, s - 1];
      for {j in Sectors}{
        let INV[r, i, j, s - 1] := inv[r, i, j, LInf, s - 1];
        let JAC_ID[r, i, j, s - 1] := max{ii in Sectors, jj in Sectors}(
            INV[r, i, j, s] / INV_SHR[r, i , j]
            - (INV[r, i, ii, s - 1] / INV_SHR[r, i, ii])
              / (INV[r, jj, ii, s - 1] / INV_SHR[r, jj, ii])
              * (INV[r, jj, j, s - 1] / INV_SHR[r, jj, j])
          );
      };
#-----------update kapital (CJ call this the simulation step)
      fix kap[r, i, LInf, s] := KAP[r, i, s];
#-----------and give the other variables a warm start
      for {t in LookForward}{
        let con[r, i, t, s] := con[r, i, t, s - 1];
        let lab[r, i, t, s] := lab[r, i, t, s - 1];
        let kap[r, i, t, s] := kap[r, i, t + 1, s - 1];
        for {j in Sectors}{let inv[r, i, j, t, s] := inv[r, i, j, t, s - 1];};
      };
    };
    for {i in Sectors}{
#-----------save actual path values of market clearing to parameter
      let MKT_CLR[i, s - 1] := sum{rr in Regions}(
        E_OUTPUT[rr, i, s - 1] 
        - CON[rr, i, s - 1]
        - INV_SUM[rr, i, s - 1 ] 
        - ADJ_COST_KAP[rr, i, s - 1]
      );
    };
      
#-----------set and solve the plan for start time s
    objective pres_disc_val[s];
    #let InstanceName := ("~/Documents/_uq-aibe/maiwar-solutions/"
    #  & time() & "maiwar" & card(Regions) & "x" & card(Sectors)
    #  & "x" & card(LookForward) & "x" & card(PathTimes) & "s" & s & "log");
    #write ("b" & InstanceName);
    solve;
    #solution (InstanceName & ".sol");
    display s, _ampl_elapsed_time, _total_solve_time,
      E_OUTPUT, CON, INV_SUM, ADJ_COST_KAP, MKT_CLR, LAB, KAP;
  };
};
/*=============================================================================
Post processing: Euler Equation and Jacobi identity check
=============================================================================*/
param jacobi_id 'Intertemporal constraints on investment'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward, s in PathTimes
    : j <> 'PbSc' and i <> j}
      = inv[r, i, j, t, s].val
           - (inv[r, i, 'PbSc', t, s].val / INV_SHR[r, i, 'PbSc'])
             / (inv[r, 'PbSc', 'PbSc', t, s].val / INV_SHR[r, 'PbSc', 'PbSc'])
             * (inv[r, 'PbSc', j, t, s].val / INV_SHR[r, 'PbSc', j])
             * INV_SHR[r, i, j];
         #   / (inv[r, j, j, t, s].val / INV_SHR[r, j, j])
         #   * (inv[r, j, j, t, s].val / INV_SHR[r, j, j])
param max_jacobi_id 'the worst of the jacobi identities'
  = max{r in Regions, i in Sectors, j in Sectors, t in LookForward, s in 
      PathTimes: j <> 'PbSc' and i <> j}  jacobi_id[r, i, j, t, s];
display max_jacobi_id, _ampl_time, _total_solve_time;
