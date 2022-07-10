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
param LSup 'supremum of the look forward set' > LInf, default 48e+0; 
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
param BETA 'discount factor', in (0, 1), default .975; 
param DELTA 'rate of depreciation for kapital' {Sectors} default .025 >= 0;
param PHI_ADJ 'kapital adjustment cost' {Sectors} default 5e-1 >= 0;
param GAMMA 'intertemporal elasticity of subst.' {Regions} default 5e-1 >= 0;
param EPS 'Frisch elasticity of labour supply'{Regions} default 5e-1 >= 0;
#-----------output shares of each component
param SHR_KAP_OUT 'importance of capital in production' {Sectors}
  default 33e-2 in (0, 50e-2);
param SHR_LAB_OUT 'share of labour in output' {Sectors}
  default 33e-2 in (0, 50e-2);
param EPS_CON 'elasticity of substitution for consumption flows'
  default 1e-1 > 0; 
param EPS_OUT 'elasticity of substitution for output components'
  default 5e-1 > 0; 
param EPS_INT 'elasticity of substitution for intermediate flows'
  default 1e-1 > 0; 
param EPS_INV 'elasticity of substitution for kapital flows'
  default 1e-1 > 0; 
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
param DUAL_MKT_CLR 'lagrange multiplier for market clearing constraint'
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
param EPS_HAT 'utility parameter'{r in Regions}, = 1 + 1 / EPS[r];
param RHO_INV 'exponent of the investment ces aggregator',
  = (EPS_INV - 1) / EPS_INV; 
param RHO_INV_HAT 'inverse of RHO_INV', = 1 / RHO_INV;
param RHO_INT 'exponent of the intermediate ces aggregator',
  = (EPS_INT - 1) / EPS_INT; 
param RHO_INT_HAT 'inverse of RHO_INT', = 1 / RHO_INT;
param RHO_OUT 'exponent of the output ces aggregator',
  = (EPS_OUT - 1) / EPS_OUT; 
param RHO_OUT_HAT 'inverse of RHO_OUT', = 1 / RHO_OUT;
param RHO_CON 'exponent of the CON ces aggregator',
  = (EPS_CON - 1) / EPS_CON; 
param RHO_CON_HAT 'inverse of RHO_CON', = 1 / RHO_OUT;

param SHR_KAPLAB_OUT 'combined importance of kapital and labour in output'
  {i in Sectors}, = SHR_KAP_OUT[i] + SHR_LAB_OUT[i];
param SHR_INT_OUT 'share of intermediates in output' {i in Sectors}
  = 1 - SHR_KAPLAB_OUT[i];
#-----------output component shares for CES functions
param SHR_KAP_OUT_CES 'importance of kapital in production' {i in Sectors}
  = SHR_KAP_OUT[i] ** (1 / EPS_OUT);
param SHR_LAB_OUT_CES 'importance of labour in production' {i in Sectors}
  = SHR_LAB_OUT[i] ** (1 / EPS_OUT);
param SHR_KAPLAB_OUT_CES 'combined importance of kapital and labour in prod'
  {i in Sectors} = SHR_KAPLAB_OUT[i] ** (1 / EPS_OUT);
param SHR_INT_OUT_CES 'importance of intermediates in production'
  {i in Sectors} = SHR_INT_OUT[i] ** (1 / EPS_OUT);
#-----------productivity and relative importance of labour in utility
param A 'productivity trend' {i in Sectors}
  default 1; #(1 - (1 - DELTA[i]) * BETA) / (SHR_KAP_OUT_CES[i] * BETA) >= 0;
param B 'importance of disutility of labour (weight in utility function)' 
  {r in Regions}
    default 1;
  #(1 - SHR_KAP_OUT_CES[i]) * A[i] * (A[i] - DELTA[i]) ** (-1 / GAMMA[r]) >= 0;
param A_UTL 'importance of consumption in utility'
  default 1;
param A_VAL 'Calibration factor for terminal value function'
  default 1;
param A_INV 'Calibration factor for investment'
  default 1;
param A_INT 'Calibration factor for intermediate bundle'
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
    = (RAW_CON_FLW[r, i] / CON_FLW_SUM[r]) ** (1 / EPS_CON);
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
    = (RAW_INV_FLW[r, i, j] / INV_FLW_RSUM[r, j]) ** (1 / EPS_INV);
    #= RAW_INV_FLW[r, i, j] ** (1 / EPS_INV);
#-----------inputs of sector i into j's intermediate bundle).
param INT_FLW_RSUM {r in Regions, j in Sectors} 
  = sum{i in Sectors} RAW_INT_FLW[r, i, j];
param SEC_SHR_INT_CES "sectoral share of i in j's CES intermediate aggregator"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_INT_FLW[r, i, j] / INT_FLW_RSUM[r, j]) ** (1 / EPS_INV);
    #= RAW_INT_FLW[r, i, j] ** (1 / EPS_INV);
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
#      + DUAL_MKT_CLR[i, s] * (
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
    ** (RHO_CON_HAT * SCALE_CON) * A_VAL;
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
  = prod{j in Sectors} B[r] * lab[r, j , t] ** LAB_SHR[r, j];
var lab_sec_Q 'quadratic labour aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = sum{j in Sectors} (lab[r, j , t] ** 2);
var lab_sec_bangF 'Frisch labour aggregate (across sectors) flat level sets'
  {r in Regions, t in LookForward}
  #= (sum{j in Sectors} lab[r, j , t]) ** EPS_HAT[r] / EPS_HAT[r];
  = B[r] * (sum{j in Sectors} lab[r, j , t]) ** EPS_HAT[r];
var lab_sec_concF 'Frisch labour aggregate (across sectors) concave level sets'
  {r in Regions, t in LookForward}
  #= (sum{j in Sectors} lab[r, j , t]) ** EPS_HAT[r] / EPS_HAT[r];
  = B[r] * sum{j in Sectors} lab[r, j , t] ** EPS_HAT[r];
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
var utility_CES_concF 'Const. Elast. Subst. and conc. Frisch instant. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CES[r, t] - lab_sec_concF[r, t]);
var utility_CES_bangF 'Const. Elast. Subst. and conc. Frisch instant. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CES[r, t] - lab_sec_bangF[r, t]);
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
      - sum{i in Sectors} B[r] * 1 ** EPS_HAT[r] / EPS_HAT[r]
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
    - sum{i in Sectors} B[r] * 1 ** EPS_HAT[r] / EPS_HAT[r]
  )) / (1 - BETA);
var tail_val_CESutl_bangF_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( A_UTL * (
    sum{i in Sectors} CON_SHR_CES[r, i] * (TAIL_CON_SHR * A[i] * (
      SHR_KAP_OUT_CES[i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_INT_OUT_CES[i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[i] * 1 ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)) ** RHO_CON
  ) ** (RHO_CON_HAT * SCALE_CON) 
    #- (sum{i in Sectors} 1) ** EPS_HAT[r] / EPS_HAT[r]
    - B[r] * (sum{i in Sectors} 1) ** EPS_HAT[r]
  )) / (1 - BETA);
var tail_val_CESutl_concF_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( A_UTL * (
    sum{i in Sectors} CON_SHR_CES[r, i] * (TAIL_CON_SHR * A[i] * (
      SHR_KAP_OUT_CES[i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_INT_OUT_CES[i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[i] * 1 ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)) ** RHO_CON
  ) ** (RHO_CON_HAT * SCALE_CON) 
    - B[r] * sum{i in Sectors} 1 ** EPS_HAT[r]
  )) / (1 - BETA);
/*=============================================================================
Current intermediate variables (substituted out during pre-solving)
=============================================================================*/
var E_output 'current intermediate variable for output'
  {r in Regions, i in Sectors, t in LookForward}
    = E_output_CES[r, i, t];
var inv_sec 'current intermediate variable for aggregated investment'
  {r in Regions, j in Sectors, t in LookForward}
# for nimrod:
    = __INVSEC__[r, j, t] * A_INV;
    #= inv_sec_CES[r, j, t] * A_INV;
var adj_cost_kap 'current adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward}
    = adj_cost_kap_Q[r, i, t];
var utility 'current intermediate variable for utility'
  {t in LookForward} = utility_CES_F[t]; 
var tail_val 'current intermediate variable for tail value function'
  = tail_val_CESutl_concF_CESout * A_VAL;
