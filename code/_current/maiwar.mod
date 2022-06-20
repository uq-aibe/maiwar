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
param LSup 'end period/supremum of the look forward set' > LInf, default 48e+0; 
param PInf 'infimum of times on a path', default 0;
param PSup 'supremum of times on a path (T_star in CJ)'# eg 2050 - 2022 = 28 
  default 3 >= PInf; 
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
set PathTimes 'times along a path' = {PInf .. PSup - 1} ordered;
set PathTimesClosure 'includes extra time for closing out the problem'
  = LookForward union {PSup} ordered;
#-----------The set of paths are indexed by
set PathSpace 'path space: the set of paths'
  # for default we adopt a two-state Markov chain with unique absorbing state
    default {PInf .. PSup + 1};  
set AllTimes 'all time periods associated with a path' = {TInf .. TSup};                      
/*=============================================================================
Parameters
=============================================================================*/
param StrtTm 'start time for each step on each path' {PathTimes, PathSpace}
  default time();
param EndTm 'start time for each step on each path' {PathTimes, PathSpace}
  default time();
param BETA 'discount factor', in (0, 1), default .98; 
param DELTA 'rate of depreciation for kapital' {Sectors} default .025 >= 0;
param PHI_ADJ 'kapital adjustment cost' {Sectors} default 5e-1 >= 0;
param GAMMA 'intertemporal elasticity of subst.' {Regions} default 5e-1 >= 0;
param ETA 'Frisch elasticity of labour supply'{Regions} default 5e-1 >= 0;
#-----------output shares of each component
param SHR_KAP_OUT 'importance of capital in production' {Sectors}
  default 33e-2 in (0, 50e-2);
param SHR_LAB_OUT 'share of labour in output' {Sectors}
  default 33e-2 in (0, 50e-2);
param EoS_CON 'elasticity of substitution for consumption flows'
  default 1e-1 > 0; 
param EoS_OUT 'elasticity of substitution for output components'
  default 5e-1 > 0; 
param EoS_INT 'elasticity of substitution for intermediate flows'
  default 1e-1 > 0; 
param EoS_INV 'elasticity of substitution for kapital flows' default 1e-1 > 0; 
param SCALE_INV 'investment scale parameter' default 90e-2; # in (0, 1);
param SCALE_CON 'consumption scale parameter' default 90e-2; # in (0, 1);
param SCALE_OUT 'output scale parameter' default 90e-2; # in (0, 1);
param SCALE_INT 'intermediate scale parameter' default 90e-2; # in (0, 1);

param INV_MIN 'lower bound of investment' default 0 >= 0;
#param kmin 'smallest capital' default 1e-1 >= 0;
#param kmax 'largest capital' default 1e+1 >= 0;
param ZETA1 'TFP before shock' default 1 >= 0;
#param ZETA2 'TFP after shock' default .95 >= 0;
param ZETA2 'TFP after shock' default 1 >= 0;
#param ZETA2 'TFP after shock' default __ZETA2__ >= 0;
param PROB2 'one period probability of jump in TFP' default 0.01 >= 0;
param TAIL_CON_SHR 'tail consumption share (of output)' default 0.45 >= 0;
#
param UInf 'infimum of interval for uniform dbn' default 0.4999 in [0, .5);
param USup 'supremum of interval for uniform dbn' default 0.5001 in (.5, 1];
param VInf 'infimum of interval for basic variables' default 1e-4;
param VSup 'supremum of interval for basic variables' default 1e+4;
param OInf 'infimum of interval for observed/actual values' default 1e-7;
param OSup 'supremum of interval for observed/actual values' default 1e+7;
#-----------raw flow parameters (these would be replaced with normalised data)
param RAW_CON_FLW {Regions, Sectors} default Uniform(UInf, USup) >= 0;
param RAW_LAB_FLW {Regions, Sectors} default Uniform(UInf, USup) >= 0;
param RAW_INV_FLW {Regions, Sectors, Sectors} 
  default Uniform(UInf, USup) >= 0;
param RAW_INT_FLW {Regions, Sectors, Sectors} 
  default Uniform(UInf, USup) >= 0;
#-----------parameters for storing (observable) path values
param CON 'observed consumption' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param INV_SEC 'observed investment' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param INV_SUM 'observed total investment' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param INT_SUM 'observed total intermediate flows' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param LAB 'observed labour' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param KAP 'observed kapital' {Regions, Sectors, PathTimesClosure}
  default 1e+0; # in (OInf, OSup);
param E_OUTPUT 'observed Exp. output' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param ADJ_COST_KAP 'observed adjustment costs of kapital'
  {Regions, Sectors, PathTimes} default 0; # in [0, OSup);
param MKT_CLR 'observed output' {Sectors, PathTimes}
  default 0; # in (-1e-4, 1e-4); 
param DUAL_KAP 'lagrange multiplier for kapital accumulation'
  {Regions, Sectors, PathTimes} default 1e+0; # in (-OSup, OSup);
param DUAL_MCL 'lagrange multiplier for market clearing constraint'
  {Sectors, PathTimes} default 1e+0;# in [0, OSup);
param GROWTH_KAP 'lagrange multiplier for market clearing constraint'
  {Regions, Sectors, PathTimes} default 5e-2; # in (-1, 1);
param EULER_INTEGRAND 'Euler error integrand'
  {Regions, Sectors, PathTimesClosure} default 1; # in (-OSup, OSup);
param EULER_RATIO 'Expected Euler ratio'
  {Regions, Sectors, PathTimes} default 1; # in (-1e+2, 1e+2);
/*=============================================================================
Computed parameters
=============================================================================*/
param GAMMA_HAT 'utility parameter' {r in Regions}, = 1 - 1 / GAMMA[r];
param ETA_HAT 'utility parameter'{r in Regions}, = 1 + 1 / ETA[r];
param RHO_INV 'exponent of the investment ces aggregator',
  = (EoS_INV - 1) / EoS_INV; 
param RHO_INV_HAT 'inverse of RHO_INV', = 1 / RHO_INV;
param RHO_INT 'exponent of the intermediate ces aggregator',
  = (EoS_INT - 1) / EoS_INT; 
param RHO_INT_HAT 'inverse of RHO_INT', = 1 / RHO_INT;
param RHO_OUT 'exponent of the output ces aggregator',
  = (EoS_OUT - 1) / EoS_OUT; 
param RHO_OUT_HAT 'inverse of RHO_OUT', = 1 / RHO_OUT;
param RHO_CON 'exponent of the CON ces aggregator', = (EoS_CON - 1) / EoS_CON; 
param RHO_CON_HAT 'inverse of RHO_CON', = 1 / RHO_OUT;

param SHR_KAPLAB_OUT 'combined importance of kapital and labour in output'
  {i in Sectors}, = SHR_KAP_OUT[i] + SHR_LAB_OUT[i];
param SHR_INT_OUT 'share of intermediates in output' {i in Sectors}
  = 1 - SHR_KAPLAB_OUT[i];
#-----------output component shares for CES functions
param SHR_KAP_OUT_CES 'importance of kapital in production' {i in Sectors}
  = SHR_KAP_OUT[i] ** (1 / EoS_OUT);
param SHR_LAB_OUT_CES 'importance of labour in production' {i in Sectors}
  = SHR_LAB_OUT[i] ** (1 / EoS_OUT);
param SHR_KAPLAB_OUT_CES 'combined importance of kapital and labour in prod'
  {i in Sectors} = SHR_KAPLAB_OUT[i] ** (1 / EoS_OUT);
param SHR_INT_OUT_CES 'importance of intermediates in production'
  {i in Sectors} = SHR_INT_OUT[i] ** (1 / EoS_OUT);
#-----------productivity and relative importance of labour in utility
param A 'productivity trend' {i in Sectors}
  default 1; #(1 - (1 - DELTA[i]) * BETA) / (SHR_KAP_OUT_CES[i] * BETA) >= 0;
param B 'importance of disutility of labour (weight in utility function)' 
  {r in Regions, i in Sectors}
    default 1;
  #(1 - SHR_KAP_OUT_CES[i]) * A[i] * (A[i] - DELTA[i]) ** (-1 / GAMMA[r]) >= 0;
param C 'importance of consumption in utility'
  default 1;
param CAL_FAC_TAIL 'Calibration factor for terminal value function'
  default 1;
param CAL_FAC_INV 'Calibration factor for investment'
  default 1;
param CAL_FAC_INT 'Calibration factor for intermediate bundle'
  default 1;
param  REG_WGHT 'regional (population) weights' {r in Regions}
  default 1 / card(Regions);
#-----------set the seed for the random number generator for weights
option randseed 12345;
#-----------Share parameters in the aggregators for utility and production.
#-----------For default values, we draw from the uniform distribution.
#-----------consumption
param CON_FLW_SUM {r in Regions} = sum{i in Sectors} RAW_CON_FLW[r, i];
param CON_SHR 'consumption weights for each good in utility (for Cobb-Doug)'
  {r in Regions, i in Sectors} = (RAW_CON_FLW[r, i] / CON_FLW_SUM[r]);
  #{r in Regions, i in Sectors} = RAW_CON_FLW[r, i];
param CON_SHR_CES 'CES consumption weights for each good in utility'
  {r in Regions, i in Sectors}
    = (RAW_CON_FLW[r, i] / CON_FLW_SUM[r]) ** (1 / EoS_CON);
#-----------labour
param LAB_FLW_SUM {r in Regions} = sum{i in Sectors} RAW_LAB_FLW[r, i];
param LAB_SHR 'labour weights for each sector in utility'
  {r in Regions, j in Sectors} = RAW_LAB_FLW[r, j] / LAB_FLW_SUM[r];
  #{r in Regions, j in Sectors} = RAW_LAB_FLW[r, j];
#-----------inputs of sector i into j's investment).
param INV_FLW_RSUM {r in Regions, j in Sectors} 
  = sum{i in Sectors} RAW_INV_FLW[r, i, j];
param SEC_SHR_INV_CES "sectoral share of i in j's CES investment aggregator"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_INV_FLW[r, i, j] / INV_FLW_RSUM[r, j]) ** (1 / EoS_INV);
    #= RAW_INV_FLW[r, i, j] ** (1 / EoS_INV);
#-----------inputs of sector i into j's intermediate bundle).
param INT_FLW_RSUM {r in Regions, j in Sectors} 
  = sum{i in Sectors} RAW_INT_FLW[r, i, j];
param SEC_SHR_INT_CES "sectoral share of i in j's CES intermediate aggregator"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_INT_FLW[r, i, j] / INT_FLW_RSUM[r, j]) ** (1 / EoS_INV);
    #= RAW_INT_FLW[r, i, j] ** (1 / EoS_INV);
/*-----------------------------------------------------------------------------
#-----------uncertainty parameters
-----------------------------------------------------------------------------*/
param Pr_shk 'probability of SHK'
  {Regions, Sectors, t in LookForward} = (1 - PROB2) ** t;
param E_shk 'expected shock (exogenous)'
  {r in Regions, i in Sectors, t in LookForward}
    = ZETA2 + Pr_shk[r, i, t] * (ZETA1 - ZETA2);
/*=============================================================================
Basic (economic) variables
=============================================================================*/
var con 'consumption flows' {Regions, Sectors, LookForward}
  in [VInf, VSup] default 1e-0;
var inv 'investment flows'
  {r in Regions, i in Sectors, j in Sectors, LookForward}
  in [VInf, VSup] default DELTA[j] * KAP[r, j, PInf];
var int 'intermediate flows'
  {Regions, Sectors, Sectors, LookForward}
    in [VInf, VSup] default 1e-0;
var lab 'labour flows' {Regions, Sectors, LookForward}
  in [VInf, VSup] default 1e-0;
#-----------kapital, the dynamic variable, is defined on LookForwardClosure
var kap 'kapital stocks (dynamic: defined on LookForwardClosure)'
  {r in Regions, j in Sectors, LookForwardClosure}
  in [VInf, VSup] default KAP[r, j, PInf]; 
/*=============================================================================
Computed variables for paths
=============================================================================*/
#var inv_sum 'sum over destinations for saving'
# {r in Regions, i in Sectors, s in PathTimes}
#   = sum{j in Sectors} inv[r, i, j, LFwd, s];
#var kap_growth {r in Regions, i in Sectors, s in PathTimes}
#  = (kap[r, i, LFwd + 1] - KAP[r, i, LFwd]) / KAP[r, i, s];
##-----------Euler integrand for CES production
#    let EULER_INTEGRAND[r, i, s] := DUAL_KAP[r, i, s] * (1 - DELTA[i]) 
#      + DUAL_MCL[i, s] * (
#        SHR_KAP_OUT_CES[i] * (KAP[r, i, s] / E_OUTPUT[r, i, s]) ** (RHO_INV - 1)
#        - PHI_ADJ[i] * (2 * GROWTH_KAP[r, i, s] + GROWTH_KAP[r, i, s] ** 2)
#      );
#    if s > PInf then 
#      let EULER_RATIO[r, i, s] 
#        := BETA * EULER_INTEGRAND[r, i, s] / DUAL_KAP[r, i, s - 1];
#  };
/*=============================================================================
Potential intermediate variables (substituted out during pre-solving)
=============================================================================*/
#-----------variety of consumption aggregator functions 
var con_sec_CD 'Cobb-Douglas consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = prod{i in Sectors} con[r, i, t] ** (CON_SHR[r, i] * SCALE_CON);
var con_sec_CES 'Const. Elast. Subst. consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = (sum{i in Sectors} CON_SHR_CES[r, i] * con[r, i, t] ** RHO_CON)
    ** (RHO_CON_HAT * SCALE_CON) * C;
var con_sec_SumPow 'Sum of power consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = sum{i in Sectors}
    CON_SHR[r, i] * con[r, i, t] ** GAMMA_HAT[r] / GAMMA_HAT[r];
var con_sec_SumShr 'Sum of fractional powers from consumption shares'
  {r in Regions, t in LookForward}
  = sum{i in Sectors} con[r, i, t] ** CON_SHR[r, i];
#-----------variety of investment aggregator functions 
var inv_sec_CES 'Const. Elast. Subst. investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SEC_SHR_INV_CES[r, i, j] * inv[r, i, j, t] ** RHO_INV)
    ** (RHO_INV_HAT * SCALE_INV);
var inv_sec_CD 'Cobb-Douglas investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = prod{i in Sectors}
    inv[r, i, j, t] ** (SEC_SHR_INV_CES[r, i, j] * SCALE_INV);
#-----------variety of intermediate aggregator functions 
var int_sec_CES 'Const. Elast. Subst. intermediate aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SEC_SHR_INT_CES[r, i, j] * int[r, i, j, t] ** RHO_INT)
    ** (RHO_INT_HAT * SCALE_INT);
#-----------variety of labour aggregator functions 
var lab_sec_CD 'Cobb-Douglas labour aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = prod{j in Sectors} B[r, j] * lab[r, j , t] ** LAB_SHR[r, j];
var lab_sec_Q 'quadratic labour aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = sum{j in Sectors} (lab[r, j , t] ** 2);
var lab_sec_F 'Frisch labour aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = (sum{j in Sectors} B[r, j] * lab[r, j , t] ** ETA_HAT[r]) / ETA_HAT[r];
#-----------variety of adjustment cost functions 
var adj_cost_kap_Q 'quadratic adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward}
  =PHI_ADJ[i] * kap[r, i, t]
      * (kap[r, i, t + 1] / kap[r, i, t] - 1) ** 2;
#-----------variety of production functions 
var E_output_CD 'Cobb-Douglas output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i]
      * (kap[r, i, t] ** SHR_KAP_OUT[i] * lab[r, i, t] ** (1 - SHR_KAP_OUT[i]))
        ** SCALE_OUT;
var E_output_CES 'Constant Elasticity of Substitution output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i] * (
    SHR_KAP_OUT_CES[i] * kap[r, i, t] ** RHO_OUT
    + SHR_LAB_OUT_CES[i] * lab[r, i, t] ** RHO_OUT 
    + SHR_INT_OUT_CES[i] * int_sec_CES[r, i, t] ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT);
var E_output_ATA 'Atalay output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i] * (
    SHR_KAPLAB_OUT_CES[i]
      * (kap[r, i, t] ** SHR_KAP_OUT[i] * lab[r, i, t] ** SHR_LAB_OUT[i])
        ** RHO_OUT
    + SHR_INT_OUT_CES[i] * int_sec_CES[r, i, t] ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT);
#-----------variety of utility functions
var utility_CD 'Cobb-Douglas instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_CD[r, t]);
var utility_CD_Q 'Cobb-Douglas instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_Q[r, t]);
var utility_pow_CD_Q 'Power of Cobb-Douglas and Quadratic instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      (REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_Q[r, t]))
        ** GAMMA_HAT[r] / GAMMA_HAT[r];
var utility_CD_F 'Cobb-Douglas and Frisch instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_F[r, t]);
var utility_CES_F 'Cobb-Douglas and Frisch instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CES[r, t] - lab_sec_F[r, t]);
var utility_CES_Q 'Const. Elast. Subst. and Quadratic instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CES[r, t] - lab_sec_Q[r, t]);
var utility_powCES_Q 'Const. Elast. Subst. and Quadratic instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      5e-2 * REG_WGHT[r] * ((con_sec_CES[r, t]) ** GAMMA_HAT[r] / GAMMA_HAT[r]
        - lab_sec_Q[r, t]);
var utility_SumShr_Q 'utility: SumShr for consumption, quadratic for labour'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_SumShr[r, t] - lab_sec_Q[r, t]);
var utility_SumPow_Q 'utility: SumPow for consumption and quadratic for labour'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_SumPow[r, t] - lab_sec_Q[r, t]);
#-----------variety of tail or terminal value functions
var tail_val_CD_F 'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions}
    REG_WGHT[r] * (
      prod{i in Sectors}(
        (TAIL_CON_SHR * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[i] * SCALE_OUT))
          ** (CON_SHR[r, i] * SCALE_CON)
      )
      - sum{i in Sectors} B[r, i] * 1 ** ETA_HAT[r] / ETA_HAT[r]
    )) / (1 - BETA);
var tail_val_CD_Q 'SumShr continuation value from time LSup + LInf onwards'
  = (sum{r in Regions}
    REG_WGHT[r] * (
      prod{i in Sectors}(
        (TAIL_CON_SHR * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[i] * SCALE_OUT))
          ** (CON_SHR[r, i] * SCALE_CON)
      )
      - 1 ** 2
    )) / (1 - BETA);
var tail_val_CESutl_Q_CDout 'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ((
    sum{i in Sectors} CON_SHR_CES[r, i] * (TAIL_CON_SHR 
      * A[i] * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[i] * SCALE_OUT)
    ) ** RHO_CON
    ) ** (RHO_CON_HAT * SCALE_CON) #(RHO_INV_HAT + GAMMA_HAT[r]) / GAMMA_HAT[r]
    - 1 ** 2
  )) / (1 - BETA);
var tail_val_CDutl_F_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * (
    prod{i in Sectors} (TAIL_CON_SHR * A[i] * (
        SHR_KAP_OUT_CES[i] * kap[r, i, LSup + LInf] ** RHO_OUT
        + SHR_INT_OUT_CES[i] * 1 ** RHO_OUT
        + SHR_LAB_OUT_CES[i] * 1 ** RHO_OUT
      ) ** (RHO_OUT_HAT * SCALE_OUT) 
    )  ** (CON_SHR[r, i] * SCALE_CON)
    - sum{i in Sectors} B[r, i] * 1 ** ETA_HAT[r] / ETA_HAT[r]
  )) / (1 - BETA);
var tail_val_CESutl_F_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( C * (
    sum{i in Sectors} CON_SHR_CES[r, i] * (TAIL_CON_SHR * A[i] * (
      SHR_KAP_OUT_CES[i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_INT_OUT_CES[i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[i] * 1 ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)) ** RHO_CON
  ) ** (RHO_CON_HAT * SCALE_CON) 
    - sum{i in Sectors} B[r, i] * 1 ** ETA_HAT[r] / ETA_HAT[r]
  )) / (1 - BETA);
/*=============================================================================
Current intermediate variables (substituted out during pre-solving)
=============================================================================*/
var E_output 'current intermediate variable for output'
  {r in Regions, i in Sectors, t in LookForward}
    = E_output_CES[r, i, t];
var inv_sec 'current intermediate variable for aggregated investment'
  {r in Regions, j in Sectors, t in LookForward}
    #= __INVSEC__[r, j, t];
    = inv_sec_CES[r, j, t] * CAL_FAC_INV;
var adj_cost_kap 'current adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward}
    = adj_cost_kap_Q[r, i, t];
var utility 'current intermediate variable for utility'
  {t in LookForward} = utility_CES_F[t]; 
var tail_val 'current intermediate variable for tail value function'
  = tail_val_CESutl_F_CESout * CAL_FAC_TAIL;
/*=============================================================================
The objectives and constraints
=============================================================================*/
maximize pres_disc_val 'present discounted value of utilities':
    sum{t in LookForward} BETA ** (t - LInf) * utility[t]
      + BETA ** (LSup - LInf) * tail_val;
subject to kap_transition 'equation for the accumulation of kapital'
  {r in Regions, j in Sectors, t in LookForward}:
    kap[r, j, t + 1] - (1 - DELTA[j]) * kap[r, j, t] - inv_sec[r, j, t] = 0;
subject to market_clearing 'market clearing for each sector and time'
  {i in Sectors, t in LookForward}:
    sum{r in Regions}(
      con[r, i, t]
      + sum{j in Sectors}(inv[r, i, j, t])
      + sum{j in Sectors}(int[r, i, j, t])
      + adj_cost_kap[r, i, t]
      - E_output[r, i , t] 
      ) = 0;
#subject to jacobi_id 'Intertemporal constraints on investment'
#  {r in Regions, i in Sectors, j in Sectors, t in LookForward
#    ii in Sectors: 1 < ord(j) and i <> j}:
#      inv[r, i, j, t]
#        = (inv[r, i, ii, t] / SEC_SHR_INV_CES[r, i, ii])
#          / (inv[r, ii, ii, t] / SEC_SHR_INV_CES[r, ii, ii])
#          * (inv[r, ii, j, t] / SEC_SHR_INV_CES[r, ii, j])
#          * SEC_SHR_INV_CES[r, i, j];
/*=============================================================================
The_data section
=============================================================================*/
data;
#-----------2x2 model
#set Regions := SEQ RoQ;
#set Sectors := Agrc Frst;
#-----------2x3 model
#set Regions := SEQ RoQ;
#set Sectors := A Mnfc PbSc;
#-----------2x4 model
#set Regions := SEQ RoQ;
#set Sectors := Agrc Frst Mnfc Srvc;
#-----------3x4 model
#set Regions := SEQ RoQ RoA;
#set Sectors := Agrc Frst Mnfc Srvc;
#-----------3x5 model
#set Regions := SEQ RoQ RoA;
#set Sectors := A Frst Mnfc Srvc PbSc;
#-----------3x6 model
#set Regions := SEQ RoQ RoA;
#set Sectors := A Frst Mnfc Srvc PbSc Utlt;
##-----------4x4
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := Agrc Frst Mnfc Srvc;# Trns Utlt;
##-----------4x5
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := A Frst Mnfc PbSc Utlt;
##-----------4x6
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := A Elct Frst Mnfc Srvc PbSc Utlt;
##-----------4x7
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := A B C D E F PbSc;
##-----------4x8
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := A B C D E F G PbSc;
##-----------4x10
#set Regions := SEQ RoQ RoA RoW;
#set Sectors := A B C D E F G H I PbSc;
#-----------7x18
#set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA;
#set Sectors := A B C D E F G H I J K L M N PbSc P Q R;
#-----------7x20
set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA;
set Sectors := A B C D E F G H I J K L M N PbSc P Q R T U;
#-----------7x21
#set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA;
#set Sectors := A B C D E F G H I J K L M N PbSc P Q R T U Al;
##-----------7x40
#set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA;
#set Sectors := A B C D E F G H I J K L M N PbSc P Q R T U
#  A1 B1 C1 D1 E1 F1 G1 H1 I1 J1 K1 L1 M1 N1 PbSc1 P1 Q1 R1 T1 U1;
##-----------7x60
#set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA;
#set Sectors := A B C D E F G H I J K L M N PbSc P Q R T U
#  A1 B1 C1 D1 E1 F1 G1 H1 I1 J1 K1 L1 M1 N1 PbSc1 P1 Q1 R1 T1 U1
#  A2 B2 C2 D2 E2 F2 G2 H2 I2 J2 K2 L2 M2 N2 PbSc2 P2 Q2 R2 T2 U2;
/*-----------------------------------------------------------------------------
#-----------set the horizon and length of paths
-----------------------------------------------------------------------------*/
let LSup := 18;
let PSup := 5;
/*-----------------------------------------------------------------------------
#-----------opportunity to tune the calibration factors (still part of data)
-----------------------------------------------------------------------------*/
display A;
for {i in Sectors}{
  let A[i] := 105e-2; #* (card(Sectors) / 20) ** (1 - 20e-2);
  let DELTA[i] := 08e-2;
  let PHI_ADJ[i] := 100e-2;
  let SHR_KAP_OUT[i] := 49e-2;
  };
for {r in Regions, i in Sectors}{
  let B[r, i] := 001e-2;
  let KAP[r, i, PInf] := 2;
};
let C := 2000e-2;
let CAL_FAC_TAIL := 100e-2;

let TAIL_CON_SHR := 045e-2;
let CAL_FAC_INV := 065e-2;
let CAL_FAC_INT := 100e-2;
let EoS_INV := 25e-2;
let EoS_INT := 10e-2;
let EoS_CON := 20e-2;
let EoS_OUT := 10e-2;
let SCALE_CON := 90e-2;
let SCALE_INV := 99e-2;
let SCALE_INT := 45e-2;
let SCALE_OUT := 110e-2;


update data;
/*=============================================================================
Solving the model
=============================================================================*/
param InstanceName symbolic;
option solver knitro;
option show_stats 1;
#-----------solve the model for a given point on a given path
for {s in PathTimes}{
  display s;
#  if s <= 6 then option solver knitro; else option solver conopt;
#-----------display some parameter values:
  display  A['PbSc'], CON_SHR['SEQ', 'PbSc'], LAB_SHR['SEQ', 'PbSc'],
  SEC_SHR_INV_CES['SEQ', 'A', 'PbSc'], DELTA['PbSc'];
  display C, TAIL_CON_SHR, CAL_FAC_INV, CAL_FAC_INT, CAL_FAC_TAIL, EoS_INV,
  RHO_INV, EoS_INT, RHO_INT, EoS_CON, RHO_CON, EoS_OUT, RHO_OUT,
  SCALE_CON, SCALE_INV, SCALE_INT, SCALE_OUT;
#-----------update kapital (CJ call this the simulation step)
  fix {r in Regions, j in Sectors} kap[r, j, LInf] := KAP[r, j, s];
#-----------in this algorithm other variables automatically get warm start
  display E_output['SEQ', 'PbSc', LInf], con['SEQ', 'PbSc', LInf];
#-----------set and solve the plan for start time s
  objective pres_disc_val;
  let InstanceName := ("maiwar"
    & card(Regions) & "x" & card(Sectors)
    & "x" & card(LookForward) & "x");
  #write ("b" & InstanceName);
#-----------call the solver
  solve;
  #solution (InstanceName & ".sol");
  display ctime(), _ampl_elapsed_time, _total_solve_time,
  _total_solve_system_time, _total_solve_user_time
  >> (InstanceName & "-results.txt");
#-----------display step values
  display E_output["SEQ", "PbSc", LInf], con["SEQ", "PbSc", LInf],
  inv_sec["SEQ", "PbSc", LInf], kap["SEQ", "PbSc", LInf],
  lab["SEQ", "PbSc", LInf], kap_transition["SEQ", "PbSc", LInf],
  market_clearing["PbSc", LInf];

  for {r in Regions, i in Sectors}{
#-----------save actual path values of variables to parameter
    let CON[r, i, s] := con[r, i, LInf];
    let INV_SEC[r, i, s] := inv_sec[r, i, LInf];
    let INV_SUM[r, i, s] := sum{j in Sectors} inv[r, i, j, LInf];
    let INT_SUM[r, i, s] := sum{j in Sectors} int[r, i, j, LInf];
    let LAB[r, i, s] := lab[r, i, LInf];
    let E_OUTPUT[r, i, s] := E_output[r, i, LInf];
    let ADJ_COST_KAP[r, i, s] := adj_cost_kap[r, i, LInf];
    let KAP[r, i, s + 1] := kap[r, i, LInf + 1];
    let DUAL_KAP[r, i, s] := kap_transition[r, i, LInf]
  };
#    let TAIL_CON_SHR := (sum{r in Regions, i in Sectors} CON[r, i, s])
#      / (sum{r in Regions, i in Sectors} E_OUTPUT[r, i, s]);
  #display E_OUTPUT, CON, INV_SUM, INT_SUM, ADJ_COST_KAP, LAB, KAP;
  #display max {r in Regions, i in Sectors} DUAL_KAP[r, i, s];
  #display min {r in Regions, i in Sectors} DUAL_KAP[r, i, s];
  for {i in Sectors}{
#-----------save actual path values of market clearing to parameter
    let MKT_CLR[i, s] := sum{rr in Regions}(
      E_OUTPUT[rr, i, s] 
      - CON[rr, i, s]
      - INV_SUM[rr, i, s] 
      - INT_SUM[rr, i, s]
      - ADJ_COST_KAP[rr, i, s]
      );
    let DUAL_MCL[i, s] := market_clearing[i, LInf];
  };
#-----------growth rate of capital as a parameter
  for {r in Regions, i in Sectors}{
    let GROWTH_KAP[r, i, s] := (KAP[r, i, s + 1] - KAP[r, i, s]) / KAP[r, i, s];
#-----------Euler integrand for Cobb-Douglas production
    #let EULER_INTEGRAND[r, i, s] :=  DUAL_KAP[r, i, s] * (1 - DELTA[i]) 
    #  + DUAL_MCL[i, s] * (
    #    SHR_KAP_OUT[i] * (KAP[r, i, s] / LAB[r, i, s]) ** (SHR_KAP_OUT[i] - 1)
    #    - PHI_ADJ[i] * (2 * GROWTH_KAP[r, i, s] + GROWTH_KAP[r, i, s] ** 2)
    #  );
#-----------Euler integrand for CES production
    let EULER_INTEGRAND[r, i, s] := DUAL_KAP[r, i, s] * (1 - DELTA[i]) 
      + DUAL_MCL[i, s] * (
        SCALE_OUT * SHR_KAP_OUT_CES[i] * A[i] ** (RHO_OUT/SCALE_OUT)
          * E_OUTPUT[r, i, s]  ** (1 - RHO_OUT/SCALE_OUT)
          * KAP[r, i, s] **(RHO_OUT - 1)
        - PHI_ADJ[i] * (2 * GROWTH_KAP[r, i, s] + GROWTH_KAP[r, i, s] ** 2)
      );
    if s > PInf then 
    let EULER_RATIO[r, i, s] 
        := BETA * EULER_INTEGRAND[r, i, s] / DUAL_KAP[r, i, s - 1];
  };
  display GROWTH_KAP, DUAL_KAP, EULER_INTEGRAND, EULER_RATIO;
  display 
    min{r in Regions, i in Sectors} E_OUTPUT[r, i, s],
    max{r in Regions, i in Sectors} E_OUTPUT[r, i, s],
    min{r in Regions, i in Sectors} CON[r, i, s],
    max{r in Regions, i in Sectors} CON[r, i, s],
    min{r in Regions, i in Sectors} INV_SUM[r, i, s],
    max{r in Regions, i in Sectors} INV_SUM[r, i, s],
    min{r in Regions, i in Sectors} INT_SUM[r, i, s],
    max{r in Regions, i in Sectors} INT_SUM[r, i, s],
    max{r in Regions, i in Sectors} ADJ_COST_KAP[r, i, s],
    max{i in Sectors} abs(MKT_CLR[i, s]),
    min{i in Sectors} DUAL_MCL[i, s],
    min{r in Regions, i in Sectors} LAB[r, i, s],
    max{r in Regions, i in Sectors} LAB[r, i, s],
    min{r in Regions, i in Sectors} KAP[r, i, s + 1],
    max{r in Regions, i in Sectors} KAP[r, i, s + 1],
    min{r in Regions, i in Sectors} GROWTH_KAP[r, i, s],
    max{r in Regions, i in Sectors} GROWTH_KAP[r, i, s],
    min{r in Regions, i in Sectors} DUAL_KAP[r, i, s];
  display utility, tail_val, pres_disc_val;
  display (sum{r in Regions, i in Sectors} con[r, i, LInf])
    / (sum{r in Regions, i in Sectors} E_output[r, i, LInf]);
  display s, _ampl_elapsed_time, _total_solve_time;
#  for {r in Regions, i in Sectors}{
#  if s > PInf then display EULER_RATIO[r, i, s] - EULER_RATIO[r, i, s - 1];
#  };
};


#display KAP;
