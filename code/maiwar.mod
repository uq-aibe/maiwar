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
  = {LInf .. LSup + LInf} ordered;
#-----------Each path is a state of the world and has the following structure:
set PathTimes 'times along a path' = {PInf .. PSup + PInf - 1} ordered;
set PathTimesClosure 'includes extra time for closing out the problem'
  = {PInf .. PSup + PInf} ordered;
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
param ALPHA 'trend growth factor', >= 1, default 100e-2; 
param ALPHA_0 'initial trend growth factor', >= 1, default 100e-2; 
param BETA 'discount factor', in (0, 1), default .975; 
param DELTA 'rate of depreciation for kapital' {Sectors} default .025 >= 0;
param PHI_ADJ 'kapital adjustment cost' {Sectors} default 5e-1 >= 0;
param GAMMA 'intertemporal elasticity of subst.' {Regions} default 5e-1 >= 0;
param EPS_LAB 'Frisch elasticity of labour supply' default 5e-1 >= 0;
#-----------output shares of each component
param SHR_KAP_OUT 'importance of capital in production' {Sectors}
  default 33e-2 in (0, 50e-2);
param SHR_LAB_OUT 'share of labour in output' {Sectors}
  default 33e-2 in (0, 50e-2);
param EPS_CON 'elasticity of substitution for consumption flows'
  default 1e-1 > 0; 
param EPS_OUT 'elasticity of substitution for output components'
  default 5e-1 > 0; 
param EPS_MED 'elasticity of substitution for intermediate flows'
  default 1e-1 > 0; 
param EPS_INV 'elasticity of substitution for kapital flows'
  default 1e-1 > 0; 
param SCALE_INV 'investment scale parameter' default 90e-2; # in (0, 1);
param SCALE_CON 'consumption scale parameter' default 90e-2; # in (0, 1);
param SCALE_OUT 'output scale parameter' default 90e-2; # in (0, 1);
param SCALE_MED 'intermediate scale parameter' default 90e-2; # in (0, 1);
param SCALE_LAB 'labour scale parameter' default 300e-2; # in (0, 1);

param INV_MIN 'lower bound of investment' default 0 >= 0;
#param kmin 'smallest capital' default 1e-1 >= 0;
#param kmax 'largest capital' default 1e+1 >= 0;
param ZETA1 'TFP before shock' default 1 >= 0;
#param ZETA2 'TFP after shock' default .95 >= 0;
param ZETA2 'TFP after shock' default 1 >= 0;
#param ZETA2 'TFP after shock' default __ZETA2__ >= 0;
param PROB2 'one period probability of jump in TFP' default 0.01 >= 0;
param TAIL_SHR_CON 'tail consumption share (of output)' default 0.45 >= 0;
#
param UInf 'infimum of interval for uniform dbn' default 0.4999 in [0, .5);
param USup 'supremum of interval for uniform dbn' default 0.5001 in (.5, 1];
param VInf 'infimum of interval for basic variables' default 1e-4;
param VSup 'supremum of interval for basic variables' default 1e+5;
param OInf 'infimum of interval for observed/actual values' default 1e-7;
param OSup 'supremum of interval for observed/actual values' default 1e+7;
param LabSup 'supremum of interval for labour values' default 66e-2;
#-----------raw flow parameters (these would be replaced with normalised data)
param RAW_CON_FLW {Regions, Sectors} default Uniform(UInf, USup) >= 0;
param RAW_LAB_FLW {Regions, Sectors} default Uniform(UInf, USup) >= 0;
param RAW_INV_FLW {Regions, Sectors, Sectors} 
  default Uniform(UInf, USup) >= 0;
param RAW_MED_FLW {Regions, Sectors, Sectors} 
  default Uniform(UInf, USup) >= 0;
#-----------parameters for storing (observable) path values
param CON 'observed consumption' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param INV_SEC 'observed investment' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param INV_SUM 'observed total investment' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param MED_SUM 'observed total intermediate flows' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param LAB 'observed labour' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param LAB_EXT 'observed laborforce participation' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param KAP 'observed kapital' {Regions, Sectors, PathTimesClosure}
  default 1e+0; # in (OInf, OSup);
param E_OUT 'observed Exp. output' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param ADJ_COST_KAP 'observed adjustment costs of kapital'
  {Regions, Sectors, PathTimes} default 0; # in [0, OSup);
param MKT_CLR 'observed output' {Sectors, PathTimes}
  default 0; # in (-1e-4, 1e-4); 
param DUAL_KAP 'lagrange multiplier for kapital accumulation'
  {Regions, Sectors, PathTimes} default 1e+0; # in (-OSup, OSup);
param DUAL_MKT_CLR 'lagrange multiplier for market clearing constraint'
  {Sectors, PathTimes} default 1e+0;# in [0, OSup);
param GROWTH_KAP 'observed growth rate for kapital'
  {Regions, Sectors, PathTimes} default 5e-2; # in (-1, 1);
param GROWTH_OUT 'observed growth rate for output'
  {Regions, Sectors, PathTimes} default 5e-2; # in (-1, 1);
param EULER_INTEGRAND 'Euler error integrand'
  {Regions, Sectors, PathTimesClosure} default 1; # in (-OSup, OSup);
param EULER_RATIO 'Expected Euler ratio'
  {Regions, Sectors, PathTimes} default 1; # in (-1e+2, 1e+2);
param NAIRE 'Non-accelerating rate of employment'
  {Regions, Sectors, PathTimes} default 95e-2;
param EXP_LAB_EXT 'Exponent of lab_ext_sec'
  {Regions, Sectors, PathTimes} default 200e-2;
param DOM 'actual path values for domestic output'
  {Regions, Sectors, PathTimes}
  default 100e-2;
/*=============================================================================
Computed parameters
=============================================================================*/
param GAMMA_HAT 'utility parameter' {r in Regions} = 1 - 1 / GAMMA[r];
param RHO_LAB 'labour exponent parameter' = 1 + 1 / EPS_LAB;
param RHO_LAB_HAT 'inverse of labour exponent parameter'
  = 1 / RHO_LAB;
param RHO_INV 'exponent of the investment ces aggregator'
  = 1 - 1 / EPS_INV; 
param RHO_INV_HAT 'inverse of RHO_INV', = 1 / RHO_INV;
param RHO_MED 'exponent of the intermediate ces aggregator'
  = 1 - 1 / EPS_MED; 
param RHO_MED_HAT 'inverse of RHO_MED', = 1 / RHO_MED;
param RHO_OUT 'exponent of the output ces aggregator'
  = 1 - 1 / EPS_OUT; 
param RHO_OUT_HAT 'inverse of RHO_OUT', = 1 / RHO_OUT;
param RHO_CON 'exponent of the CON ces aggregator'
  = 1 - 1 / EPS_CON; 
param RHO_CON_HAT 'inverse of RHO_CON', = 1 / RHO_OUT;
param SHR_KAPLAB_OUT 'combined importance of kapital and labour in output'
  {i in Sectors}, = SHR_KAP_OUT[i] + SHR_LAB_OUT[i];
param SHR_MED_OUT 'share of intermediates in output' {i in Sectors}
  = 1 - SHR_KAPLAB_OUT[i];
#-----------output component shares for CES functions
param SHR_KAP_OUT_CES 'importance of kapital in production' {i in Sectors}
  = SHR_KAP_OUT[i] ** (1 / EPS_OUT);
param SHR_LAB_OUT_CES 'importance of labour in production' {i in Sectors}
  = SHR_LAB_OUT[i] ** (1 / EPS_OUT);
param SHR_KAPLAB_OUT_CES 'combined importance of kapital and labour in prod'
  {i in Sectors} = SHR_KAPLAB_OUT[i] ** (1 / EPS_OUT);
param SHR_MED_OUT_CES 'importance of intermediates in production'
  {i in Sectors} = SHR_MED_OUT[i] ** (1 / EPS_OUT);
#-----------productivity and relative importance of labour in utility
param A 'productivity trend' {i in Sectors}
  default 1; #(1 - (1 - DELTA[i]) * BETA) / (SHR_KAP_OUT_CES[i] * BETA) >= 0;
param A_LAB 'importance of disutility of labour (weight in utility function)' 
  {Regions, LookForwardClosure}
    default -1;
param A_LAB_EXT 'disutility weight for labourforce deviations in utility'
  default -1;
  #(1 - SHR_KAP_OUT_CES[i]) * A[i] * (A[i] - DELTA[i]) ** (-1 / GAMMA[r]) >= 0;
param A_CON 'importance of consumption in utility'
  default 1;
param A_VAL 'Calibration factor for terminal value function'
  default 1;
param A_INV 'Calibration factor for investment'
  default 1;
param A_MED 'Calibration factor for intermediate bundle'
  default 1;
param  REG_WGHT 'regional (population) weights' {r in Regions}
  default 1 / card(Regions);
#-----------set the seed for the random number generator for weights
option randseed 12345;
#-----------Share parameters in the aggregators for utility and production.
#-----------For default values, we draw from the uniform distribution.
#-----------consumption
param CON_FLW_SUM {r in Regions}
  = sum{i in Sectors} RAW_CON_FLW[r, i];
#  = 1;
param SHR_CON 'consumption weights for each good in utility (for Cobb-Doug)'
  {r in Regions, i in Sectors} = (RAW_CON_FLW[r, i] / CON_FLW_SUM[r]);
  #{r in Regions, i in Sectors} = RAW_CON_FLW[r, i];
param SHR_CON_CES 'CES consumption weights for each good in utility'
  {r in Regions, i in Sectors}
    = (RAW_CON_FLW[r, i] / CON_FLW_SUM[r]) ** (1 / EPS_CON);
#-----------labour
param LAB_FLW_SUM {r in Regions}
  = sum{i in Sectors} RAW_LAB_FLW[r, i];
#  = 1;
param SHR_LAB 'labour weights for each sector in utility'
  {r in Regions, j in Sectors} = RAW_LAB_FLW[r, j] / LAB_FLW_SUM[r];
  #{r in Regions, j in Sectors} = RAW_LAB_FLW[r, j];
#-----------inputs of sector i into j's investment).
param INV_FLW_RSUM {r in Regions, j in Sectors} 
  = sum{i in Sectors} RAW_INV_FLW[r, i, j];
#  = 1;
param SHR_INV_CES "sectoral share of i in j's CES investment aggregator"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_INV_FLW[r, i, j] / INV_FLW_RSUM[r, j]) ** (1 / EPS_INV);
    #= RAW_INV_FLW[r, i, j] ** (1 / EPS_INV);
#-----------inputs of sector i into j's intermediate bundle).
param MED_FLW_RSUM {r in Regions, j in Sectors} 
  = sum{i in Sectors} RAW_MED_FLW[r, i, j];
#  = 1;
param SHR_MED_CES "sectoral share of i in j's CES intermediate aggregator"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_MED_FLW[r, i, j] / MED_FLW_RSUM[r, j]) ** (1 / EPS_INV);
    #= RAW_MED_FLW[r, i, j] ** (1 / EPS_INV);
param SHR_LAB_CES "sectoral share of i in the CES labour aggregator"
  {r in Regions, i in Sectors}
    = SHR_LAB[r, i] ** (- 1 / EPS_LAB);
/*-----------------------------------------------------------------------------
uncertainty parameters
-----------------------------------------------------------------------------*/
param Pr_shk 'probability of SHK'
  {Regions, Sectors, t in LookForward} = (1 - PROB2) ** t;
param E_shk 'expected shock (exogenous)'
  {r in Regions, i in Sectors, t in LookForward}
    = ZETA2 + Pr_shk[r, i, t] * (ZETA1 - ZETA2);
/*-----------------------------------------------------------------------------
Armington parameters
-----------------------------------------------------------------------------*/
#-----------foreign prices
param PRC_YSA 'import prices' {Sectors, LookForward} default 100e-2;
param PRC_EXA 'export prices' {Sectors, LookForward} default 100e-2;
#-----------calibration factors
param A_CCON 'calibration factor for composite consumption'
  default 100e-2;
param A_CMED 'scaling factor for composite intermediates'
  default 100e-2;
param A_CINV 'scaling factor for composite intermediates'
  default 100e-2;
#-----------economies of scale exponents
param SCALE_CINV 'economies of scale for composite intermediate'
  default 100e-2;
param SCALE_CMED 'economies of scale for composite intermediate'
  default 100e-2;
#-----------elasticities
param EPS_CINV 'elasticity of subst. for composite investment flows'
  {Regions, Sectors, Sectors}
  default 0200e-2;
param EPS_CMED 'elasticity of subst. for composite intermediates'
  {Regions, Sectors, Sectors}
  default 0200e-2;
#-----------shorthand elasticity parameters
param RHO_CINV 'CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 - 1 / EPS_CINV[r, i, j];
param RHO_CINV_HAT 'inverse of CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 / RHO_CINV[r, i, j];
param RHO_CMED 'CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 - 1 / EPS_CMED[r, i, j];
param RHO_CMED_HAT 'inverse of CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 / RHO_CMED[r, i, j];
#-----------share parameters
param SHR_CON_CCON 'share of domestic consumption in comp. consumption'
  {Regions, Sectors}
  default 50e-2;
param SHR_YSA_CMED 'import weight in composite intermediate flows'
  {r in Regions, i in Sectors, j in Sectors}
  default 50e-2;
param SHR_DOM_CMED 'import weight in composite intermediate flows'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 - SHR_YSA_CMED[r, i, j];
param SHR_YSA_CINV 'import weight in composite investment flows'
  {r in Regions, i in Sectors, j in Sectors}
  default 50e-2;
param SHR_DOM_CINV 'import weight in composite investment flows'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 - SHR_YSA_CINV[r, i, j];
/*=============================================================================
===============================================================================
The variables
===============================================================================
===============================================================================
Basic (economic) variables
=============================================================================*/
var dcon 'consumption flows' {Regions, Sectors, LookForward}
  in [VInf, VSup] default 1e-0;
var dinv 'investment flows'
  {r in Regions, i in Sectors, j in Sectors, LookForward}
  in [VInf, VSup] default DELTA[j] * KAP[r, j, PInf];
var dmed 'intermediate flows'
  {Regions, Sectors, Sectors, LookForward}
    in [VInf, VSup] default 1e-0;
/*-----------------------------------------------------------------------------
labour split into time spent working and laborforce both normalised
-----------------------------------------------------------------------------*/
var lab 'labour hours' {Regions, Sectors, LookForward}
  in [VInf, LabSup] default 33e-2;
#-----------to do: add occupations as a set in the following
var lab_ext 'active laborforce but can also be interpreted as effort'
  {r in Regions, i in Sectors, t in LookForward}
  in [0, 100e-2] default NAIRE[r, i, t];
/*-----------------------------------------------------------------------------
kapital, the dynamic variable, is defined on LookForwardClosure
-----------------------------------------------------------------------------*/
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
#        SHR_KAP_OUT_CES[i] * (KAP[r, i, s] / E_OUT[r, i, s]) ** (RHO_INV - 1)
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
var ccon 'Leontief composite consumption good (includes imports)'
  {r in Regions, i in Sectors, t in LookForward}
  = A_CCON * dcon[r, i, t] / SHR_CON_CCON[r, i];
var ccon_sec_CD 'Cobb--Douglas consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = prod{i in Sectors} ccon[r, i, t] ** (SHR_CON[r, i] * SCALE_CON);
var con_sec_CD 'Cobb--Douglas consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = prod{i in Sectors} dcon[r, i, t] ** (SHR_CON[r, i] * SCALE_CON);
var con_sec_CES 'Const. Elast. Subst. consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = (sum{i in Sectors} SHR_CON_CES[r, i] * dcon[r, i, t] ** RHO_CON)
    ** (RHO_CON_HAT * SCALE_CON); 
var con_sec_SumPow 'Sum of power consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = sum{i in Sectors}
    SHR_CON[r, i] * dcon[r, i, t] ** GAMMA_HAT[r] / GAMMA_HAT[r];
var con_sec_SumShr 'Sum of fractional powers from consumption shares'
  {r in Regions, t in LookForward}
  = sum{i in Sectors} dcon[r, i, t] ** SHR_CON[r, i];
/*-----------------------------------------------------------------------------
Armington variables
-----------------------------------------------------------------------------*/
#-----------domestic prices
var prc_dom_CDL 'domestic price: Cobb-Douglas in sectors Leontief in imports'
  {j in Sectors, t in LookForward}
  = SCALE_CON * SHR_CON['SEQ', j] * (A_CCON / SHR_CON_CCON['SEQ', j])
    * (A_CON * ccon_sec_CD['SEQ', t] / dcon['SEQ', j, t]);
var prc_dom  'domestic price / the Lagrange multiplier'
  {j in Sectors, t in LookForward}
  = prc_dom_CDL[j, t];
#-----------imports
var ymed 'intermediate imports'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = (SHR_YSA_CMED[r, i, j] / SHR_DOM_CMED[r, i, j])
    * (prc_dom[j, t] / PRC_YSA[j, t]) ** EPS_CMED[r, i, j]
    * dmed[r, i, j, t];
var cmed 'composite intermediate flows'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = A_CMED * (
    SHR_DOM_CMED[r, i, j] * dmed[r, i, j, t] ** RHO_CMED[r, i, j]
    + SHR_YSA_CMED[r, i, j] * ymed[r, i, j, t] ** RHO_CMED[r, i, j]
  ) ** (RHO_CMED_HAT[r, i, j] * SCALE_CMED);
var yinv 'investment imports'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = (SHR_YSA_CINV[r, i, j] / SHR_DOM_CINV[r, i, j])
    * (prc_dom[j, t] / PRC_YSA[j, t]) ** EPS_CINV[r, i, j]
    * dinv[r, i, j, t];
var cinv 'composite investment flows' 
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = A_CINV * (
    SHR_DOM_CINV[r, i, j] * dinv[r, i, j, t] ** RHO_CINV[r, i, j]
    + SHR_YSA_CINV[r, i, j] * yinv[r, i, j, t] ** RHO_CINV[r, i, j]
  ) ** (RHO_CINV_HAT[r, i, j] * SCALE_CINV);
#-----------exports (other variables appear after output)
param SHR_EXO_JOUT 'share of exports in output (joint production CET function)'
  {Regions, Sectors}
  default 010e-2;
param SHR_DOM_JOUT 'share of domestic uses in output (CET function)'
  {r in Regions, i in Sectors}
  = 1 - SHR_EXO_JOUT[r, i];
param PRC_EXO 'exogenous price of exports (in domestic currency units)'
  {i in Sectors, t in LookForward}
  default 100e-2;
param EPS_JOUT 'elasticity of subst. for export CET function'
  {Regions, Sectors}
  = 4;
param SHR_DOM 'observed share of output for domestic uses'
  {Regions, Sectors, PathTimes}
  default 100e-2;
var shr_dom 'share of output going to domestic uses (p_dom / (p_dom + P_EXO))'
  {r in Regions, i in Sectors, t in LookForward}
    = 1 / (
      1 + SHR_EXO_JOUT[r, i] / SHR_DOM_JOUT[r, i]
        * (PRC_EXO[i, t] / prc_dom[i, t]) ** EPS_JOUT[r, i]
    );
#var dom 'domestic production'
#  = (A_YSA[j] * SHR_DOM_OUT[j] * mkt_clr[j, t] / prc_dom[j, t])
#    ** (1 - 1 / RHO_YSA[j])
#    * out[r, j, t];
#var gdp 'gross domestic product'
#  {r in Regions, j in Sectors, t in LookForward}
#  = (
#    (A_EXA[j] ** RHO_EXA[j] * SHR_DOM_GDP[j] * (1 + TAX_GDP[j]) * mkt_clr[j, t])
#      / (A_YSA[j] ** RHO_YSA[j] * SHR_YSA_OUT[j] * mkt_clr[j, t])
#  ) ** (1  - 1 / RHO_YSA[j])
#    * out[r, j, t];
#-----------variety of investment aggregator functions 
var inv_sec_CD 'Cobb--Douglas investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = prod{i in Sectors}
    dinv[r, i, j, t] ** (SHR_INV_CES[r, i, j] * SCALE_INV);
var inv_sec_CES 'Const. Elast. Subst. investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SHR_INV_CES[r, i, j] * dinv[r, i, j, t] ** RHO_INV)
    ** (RHO_INV_HAT * SCALE_INV);
var cinv_sec_CES 'composite CES investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SHR_INV_CES[r, i, j] * cinv[r, i, j, t] ** RHO_INV)
    ** (RHO_INV_HAT * SCALE_INV);
#-----------variety of intermediate aggregator functions 
var med_sec_CES 'Const. Elast. Subst. intermediate aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SHR_MED_CES[r, i, j] * dmed[r, i, j, t] ** RHO_MED)
    ** (RHO_MED_HAT * SCALE_MED);
var cmed_sec_CES 'Const. Elast. Subst. intermediate aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SHR_MED_CES[r, i, j] * cmed[r, i, j, t] ** RHO_MED)
    ** (RHO_MED_HAT * SCALE_MED);
#-----------variety of labour aggregator functions 
var lab_sec_CD 'Cobb--Douglas labour aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = prod{j in Sectors} A_LAB[r, t] * lab[r, j , t] ** SHR_LAB[r, j];
var lab_sec_Q 'quadratic labour aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = sum{j in Sectors} (lab[r, j , t] ** 2);
var lab_sec_bangF 'Frisch labour aggregate (across sectors) flat level sets'
  {r in Regions, t in LookForward}
  #= (sum{j in Sectors} lab[r, j , t]) ** RHO_LAB / RHO_LAB;
  = (sum{j in Sectors} lab[r, j , t]) ** RHO_LAB;
var lab_sec_caveF 'Frisch labour aggregate (across sectors) concave level sets'
  {r in Regions, t in LookForward}
  = (sum{j in Sectors} SHR_LAB_CES[r, j] * lab[r, j , t] ** RHO_LAB)
    ** (RHO_LAB_HAT * SCALE_LAB);
var lab_ext_sec 'labourforce aggregator'
  {r in Regions, t in LookForward}
  = sum{j in Sectors}
    (NAIRE[r, j, t] - lab_ext[r, j, t]) ** EXP_LAB_EXT[r, j, t];
#-----------variety of adjustment cost functions 
var adj_cost_kap_Q 'quadratic adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward}
  = PHI_ADJ[i] * kap[r, i, t]
      * (kap[r, i, t + 1] / kap[r, i, t] - 1) ** 2;
#-----------variety of production functions 
var E_out_CD 'Cobb--Douglas output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i]
      * (kap[r, i, t] ** SHR_KAP_OUT[i] * lab[r, i, t] ** (1 - SHR_KAP_OUT[i]))
        ** SCALE_OUT;
var E_out_ATA 'Atalay output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i] * (
    SHR_KAPLAB_OUT_CES[i]
      * (kap[r, i, t] ** SHR_KAP_OUT[i] * lab[r, i, t] ** SHR_LAB_OUT[i])
        ** RHO_OUT
    + SHR_MED_OUT_CES[i] * med_sec_CES[r, i, t] ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT);
var E_out_CES 'Constant Elasticity of Substitution output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i] * (
    SHR_KAP_OUT_CES[i]
      * (kap[r, i, t] * lab_ext[r, i, t]) ** RHO_OUT
    + SHR_LAB_OUT_CES[i]
      * (lab[r, i, t] * ALPHA * ALPHA_0 ** (t + 1)) ** RHO_OUT 
    + SHR_MED_OUT_CES[i]
      * (med_sec_CES[r, i, t] * lab_ext[r, i, t]) ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT);
var E_cout_CES 'Constant Elasticity of Substitution output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i] * (
    SHR_KAP_OUT_CES[i]
      * (kap[r, i, t] * lab_ext[r, i, t]) ** RHO_OUT
    + SHR_LAB_OUT_CES[i]
      * (lab[r, i, t] * ALPHA * ALPHA_0 ** (t + 1)) ** RHO_OUT 
    + SHR_MED_OUT_CES[i]
      * (cmed_sec_CES[r, i, t] * lab_ext[r, i, t]) ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT);
#-----------variety of utility functions
var utility_CD 'Cobb--Douglas instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_CD[r, t]);
var utility_CD_Q 'Cobb--Douglas instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_Q[r, t]);
var utility_pow_CD_Q 'Power of Cobb-Douglas and Quadratic instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      (REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_Q[r, t]))
        ** GAMMA_HAT[r] / GAMMA_HAT[r];
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
var utility_CD_F 'Cobb--Douglas and Frisch instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_caveF[r, t]);
var utility_CES_bangF 'Const. Elast. Subst. and conc. Frisch instant. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CES[r, t] - lab_sec_bangF[r, t]);
var utility_CES_caveF 'Const. Elast. Subst. and conc. Frisch instant. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (
        A_CON * con_sec_CES[r, t]
        + A_LAB_EXT * lab_ext_sec[r, t]
        + A_LAB[r, t] * lab_sec_caveF[r, t]
      );
var utility_CD_caveF 'Cobb-Douglas and concave Frisch inst. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (
        A_CON * con_sec_CD[r, t]
        + A_LAB_EXT * lab_ext_sec[r, t]
        + A_LAB[r, t] * lab_sec_caveF[r, t]
      );
var cutility_CD_caveF 'Cobb-Douglas-Leontief and concave Frisch inst. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (
        A_CON * ccon_sec_CD[r, t]
        + A_LAB_EXT * lab_ext_sec[r, t]
        + A_LAB[r, t] * lab_sec_caveF[r, t]
      );
#-----------variety of tail or terminal value functions
var tail_val_CD_F 'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions}
    REG_WGHT[r] * (
      prod{i in Sectors}(
        (TAIL_SHR_CON * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[i] * SCALE_OUT))
          ** (SHR_CON[r, i] * SCALE_CON)
      )
      - sum{i in Sectors} A_LAB[r, LSup + LInf] * 1 ** RHO_LAB / RHO_LAB
    )) / (1 - BETA);
var tail_val_CD_Q 'SumShr continuation value from time LSup + LInf onwards'
  = (sum{r in Regions}
    REG_WGHT[r] * (
      prod{i in Sectors}(
        (TAIL_SHR_CON * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[i] * SCALE_OUT))
          ** (SHR_CON[r, i] * SCALE_CON)
      )
      - 1 ** 2
    )) / (1 - BETA);
var tail_val_CESutl_Q_CDout 'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ((
    sum{i in Sectors} SHR_CON_CES[r, i] * (TAIL_SHR_CON 
      * A[i] * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[i] * SCALE_OUT)
    ) ** RHO_CON
    ) ** (RHO_CON_HAT * SCALE_CON) #(RHO_INV_HAT + GAMMA_HAT[r]) / GAMMA_HAT[r]
    - 1 ** 2
  )) / (1 - BETA);
var tail_val_CDutl_F_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * (
    prod{i in Sectors} (TAIL_SHR_CON * A[i] * (
        SHR_KAP_OUT_CES[i] * kap[r, i, LSup + LInf] ** RHO_OUT
        + SHR_MED_OUT_CES[i] * 1 ** RHO_OUT
        + SHR_LAB_OUT_CES[i] * 1 ** RHO_OUT
      ) ** (RHO_OUT_HAT * SCALE_OUT) 
    )  ** (SHR_CON[r, i] * SCALE_CON)
    - sum{i in Sectors} A_LAB[r, LSup + LInf] * 1 ** RHO_LAB / RHO_LAB
  )) / (1 - BETA);
var tail_val_CESutl_bangF_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( A_CON * (
    sum{i in Sectors} SHR_CON_CES[r, i] * (TAIL_SHR_CON * A[i] * (
      SHR_KAP_OUT_CES[i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_MED_OUT_CES[i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[i] * 1 ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)) ** RHO_CON
  ) ** (RHO_CON_HAT * SCALE_CON) 
    #- (sum{i in Sectors} 1) ** RHO_LAB / RHO_LAB
    - A_LAB[r, LSup + LInf] * (sum{i in Sectors} 1) ** RHO_LAB
  )) / (1 - BETA);
var tail_val_CESutl_caveF_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( A_CON * (
    sum{i in Sectors} SHR_CON_CES[r, i] * (TAIL_SHR_CON * A[i] * (
      SHR_KAP_OUT_CES[i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_MED_OUT_CES[i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[i] * (1 * ALPHA * ALPHA_0 ** (LSup + LInf)) ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)) ** RHO_CON
  ) ** (RHO_CON_HAT * SCALE_CON) 
    + A_LAB_EXT * 0
    + A_LAB[r, LSup + LInf] * (
      sum{i in Sectors} SHR_LAB_CES[r, i] * 33e-2 ** RHO_LAB
      ) ** (RHO_LAB_HAT * SCALE_LAB)
  )) / (1 - BETA);
var tail_val_CDutl_caveF_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( A_CON * (
    prod{i in Sectors} (TAIL_SHR_CON * A[i] * (
      SHR_KAP_OUT_CES[i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_MED_OUT_CES[i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[i] * (1 * ALPHA * ALPHA_0 ** (LSup + LInf)) ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)
    ) ** (SHR_CON[r, i] * SCALE_CON)
    )
    + A_LAB_EXT * 0
    + A_LAB[r, LSup + LInf] * (
      sum{i in Sectors} SHR_LAB_CES[r, i] * 33e-2 ** RHO_LAB
      ) ** (RHO_LAB_HAT * SCALE_LAB)
  )) / (1 - BETA);
/*=============================================================================
Current intermediate variables (substituted out during pre-solving)
=============================================================================*/
var con 'consumption'
  {r in Regions, i in Sectors, t in LookForward}
  = ccon[r, i, t];
var inv 'domestic investment flows'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = dinv[r, i, j, t];
var med 'domestic intermedate flows'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = cmed[r, i, j, t];
var inv_sec 'current intermediate variable for aggregated investment'
  {r in Regions, j in Sectors, t in LookForward}
    = inv_sec_CES[r, j, t] * A_INV;
# for nimrod:
    #= __INVSEC__[r, j, t] * A_INV;
#var med_sec 'current intermediate variable for aggregated investment'
#  {r in Regions, j in Sectors, t in LookForward}
#    = cmed_sec_CES[r, j, t] * A_INV;
var E_out 'current intermediate variable for output'
  {r in Regions, i in Sectors, t in LookForward}
    = E_cout_CES[r, i, t];
var dom 'domestic uses'
  {r in Regions, i in Sectors, t in LookForward}
  = shr_dom[r, i, t] * E_out[r, i, t];
var exo 'exports'
  {r in Regions, i in Sectors, t in LookForward}
  = (1 - shr_dom[r, i, t]) * E_out[r, i, t];
var utility 'current intermediate variable for utility'
  {t in LookForward} = cutility_CD_caveF[t]; 
var adj_cost_kap 'current adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward}
    = adj_cost_kap_Q[r, i, t];
var tail_val 'current intermediate variable for tail value function'
  = tail_val_CDutl_caveF_CESout * A_VAL;
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
      + sum{j in Sectors}(med[r, i, j, t])
      + adj_cost_kap[r, i, t]
      - dom[r, i , t] 
      ) = 0;
#subject to jacobi_id 'Intertemporal constraints on investment'
#  {r in Regions, i in Sectors, j in Sectors, t in LookForward
#    ii in Sectors: 1 < ord(j) and i <> j}:
#      inv[r, i, j, t]
#        = (inv[r, i, ii, t] / SHR_INV_CES[r, i, ii])
#          / (inv[r, ii, ii, t] / SHR_INV_CES[r, ii, ii])
#          * (inv[r, ii, j, t] / SHR_INV_CES[r, ii, j])
#          * SHR_INV_CES[r, i, j];
/*=============================================================================
The_data section
=============================================================================*/
data;
#-----------1x20 model
set Regions := SEQ;
set Sectors := A B C D E F G H I J K L M N PbSc P Q R T U;
#set Sectors := Agrc Frst;
#-----------2x2 model
#set Regions := SEQ RoQ;
#set Sectors := Agrc Frst;
#-----------2x3 model
#set Regions := SEQ RoQ;
#set Sectors := A Mnfc PbSc;
#-----------2x4 model
#set Regions := SEQ RoQ;
#set Sectors := A Frst Mnfc PbSc;
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
#set Sectors := A Frst Mnfc PbSc;# Trns Utlt;
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
#set Regions := CnQ FNQ Mck NrQ SEQ WBB RoA;
#set Sectors := A B C D E F G H I J K L M N PbSc P Q R T U;
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
let LSup := 15;
let PSup := 61;
/*-----------------------------------------------------------------------------
#-----------opportunity to tune the calibration factors (still part of data)
-----------------------------------------------------------------------------*/
let ALPHA := 1;#271828182846e-11;
let ALPHA_0 := 1;#271828182846e-11;
let ALPHA := 271828182846e-11;
let ALPHA_0 := 271828182846e-11;
let BETA := 950e-3;
display A;
param aA := 094e-2;
param c20A := 090e-2;
param c10A := 100e-2;
for {i in Sectors}{
  let A[i] := c10A; #* (card(Sectors) / 20) ** (1 - 20e-2);
  let DELTA[i] := 05e-2;
  let PHI_ADJ[i] := 300e-2;
  let SHR_KAP_OUT[i] := 33e-2;
  };
for {r in Regions, i in Sectors}{
  let KAP[r, i, PInf] := 1;
};
let A_CON := 09100e-2; #increase this to increase labour
let A_INV := 0010e-2;
let A_MED := 0010e-2;
let A_VAL := 0001e-2;
let A_LAB_EXT := -0845e-2;
#let A_CMED := 1;
let TAIL_SHR_CON := 045e-2;

let EPS_INV := 0255e-3;
let EPS_MED := 0260e-3;
let EPS_CON := 0999e-3;
let EPS_OUT := 0100e-3;

let SCALE_CON := 200e-3;
let SCALE_INV := 999e-3;
let SCALE_MED := 999e-3;
let SCALE_OUT := 999e-3;
let SCALE_CMED := 990e-3;
let SCALE_CINV := 990e-3;
let EPS_LAB := 050e-2;
let SCALE_LAB := 1200e-2;
for {r in Regions, j in Sectors, t in LookForward}{
  fix lab[r, j, t] := 33e-2;
  let NAIRE[r, j, t] := 80e-2;
  let EXP_LAB_EXT[r, j, t] := 2;
};

update data;
/*=============================================================================
Solving the model
=============================================================================*/
param InstanceName symbolic;
option solver knitro;
#option solver conopt;
option show_stats 1;
param CH_GROWTH_OUT {Regions, Sectors, PathTimes} default 0;
#-----------solve the model for a given point on a given path
for {s in PathTimes}{
  display s, ctime();
#-----------update kapital (CJ call this the simulation step)
  fix {r in Regions, j in Sectors} kap[r, j, LInf] := KAP[r, j, s];
  if s > PInf then
  let ALPHA := ALPHA * ALPHA_0;
  for {r in Regions, t in LookForwardClosure}{
    #let A_LAB[r, t] := - 271828182846e-11 ** - ((s + t) * SCALE_LAB);
  }
  #if s = 5 then
  #  unfix {r in Regions, j in Sectors, t in LookForward}
  #    lab[r, j, t] := 31.4e-2;
  #if s = 20 then
  #  let A['A'] := (1 - 33e-2) * A['A'];
  ;
#  if s <= 6 then option solver knitro; else option solver conopt;
#-----------display some parameter values:
  display PHI_ADJ['A'], A_LAB_EXT,
  NAIRE['SEQ', 'PbSc', 0], EXP_LAB_EXT['SEQ', 'A', 0],
  LabSup, LSup, ALPHA_0, ALPHA, BETA, A['PbSc'],
  SHR_CON['SEQ', 'PbSc'], SHR_LAB['SEQ', 'PbSc'],
  SHR_INV_CES['SEQ', 'A', 'PbSc'], DELTA['PbSc'];
  display A_CON, A_INV, A_MED, A_VAL, A_LAB['SEQ', 0],
  EPS_INV, EPS_MED, EPS_CON, EPS_OUT, EPS_LAB,
  RHO_INV,  RHO_MED,  RHO_CON, RHO_OUT, RHO_LAB,
  SCALE_CON, SCALE_INV, SCALE_MED, SCALE_OUT, SCALE_LAB,
  SCALE_CMED, SCALE_CINV;
#-----------in this algorithm other variables automatically get warm start
  display E_out['SEQ', 'PbSc', LInf], con['SEQ', 'PbSc', LInf];
#-----------set and solve the plan for start time s
  objective pres_disc_val;
  let InstanceName := ("maiwar"
    & card(Regions) & "x" & card(Sectors)
    & "x" & card(LookForward) & "x");
  #write ("b" & InstanceName);
#-----------call the solver
  display s;
  solve;
  #solution (InstanceName & ".sol");
  display ctime(), _ampl_elapsed_time, _total_solve_time,
  _total_solve_system_time, _total_solve_user_time
  >> (InstanceName & "-results.txt");
#-----------display step values
  display E_out["SEQ", "PbSc", LInf], con["SEQ", "PbSc", LInf],
  inv_sec["SEQ", "PbSc", LInf], kap["SEQ", "PbSc", LInf],
  lab["SEQ", "PbSc", LInf], kap_transition["SEQ", "PbSc", LInf],
  market_clearing["PbSc", LInf];
  for {r in Regions, i in Sectors}{
#-----------save actual path values of variables to parameter
    let CON[r, i, s] := con[r, i, LInf];
    let INV_SEC[r, i, s] := inv_sec[r, i, LInf];
    let INV_SUM[r, i, s] := sum{j in Sectors} inv[r, i, j, LInf];
    let MED_SUM[r, i, s] := sum{j in Sectors} med[r, i, j, LInf];
    let LAB[r, i, s] := lab[r, i, LInf];
    let LAB_EXT[r, i, s] := lab_ext[r, i, LInf];
    let E_OUT[r, i, s] := E_out[r, i, LInf];
    let DOM[r, i, s] := dom[r, i, LInf];
    let SHR_DOM[r, i, s] := shr_dom[r, i, LInf];
    let ADJ_COST_KAP[r, i, s] := adj_cost_kap[r, i, LInf];
    let KAP[r, i, s + 1] := kap[r, i, LInf + 1];
    let DUAL_KAP[r, i, s] := kap_transition[r, i, LInf]
  };
#    let TAIL_SHR_CON := (sum{r in Regions, i in Sectors} CON[r, i, s])
#      / (sum{r in Regions, i in Sectors} E_OUT[r, i, s]);
  #display E_OUT, CON, INV_SUM, MED_SUM, ADJ_COST_KAP, LAB, KAP;
  #display max {r in Regions, i in Sectors} DUAL_KAP[r, i, s];
  #display min {r in Regions, i in Sectors} DUAL_KAP[r, i, s];
  for {i in Sectors}{
#-----------save actual path values of market clearing to parameter
    let MKT_CLR[i, s] := sum{rr in Regions}(
      DOM[rr, i, s] 
      - CON[rr, i, s]
      - INV_SUM[rr, i, s] 
      - MED_SUM[rr, i, s]
      - ADJ_COST_KAP[rr, i, s]
      );
    let DUAL_MKT_CLR[i, s] := market_clearing[i, LInf];
  };
#-----------growth rate of capital as a parameter
  for {r in Regions, i in Sectors}{
    let GROWTH_KAP[r, i, s] := (KAP[r, i, s + 1] - KAP[r, i, s]) / KAP[r, i, s];
    if s > PInf then 
    let GROWTH_OUT[r, i, s] :=
      (E_OUT[r, i, s] - E_OUT[r, i, s - 1]) / E_OUT[r, i, s - 1];
   
#-----------Euler integrand for Cobb--Douglas production
    #let EULER_INTEGRAND[r, i, s] :=  DUAL_KAP[r, i, s] * (1 - DELTA[i]) 
    #  + DUAL_MKT_CLR[i, s] * (
    #    SHR_KAP_OUT[i] * (KAP[r, i, s] / LAB[r, i, s]) ** (SHR_KAP_OUT[i] - 1)
    #    - PHI_ADJ[i] * (2 * GROWTH_KAP[r, i, s] + GROWTH_KAP[r, i, s] ** 2)
    #  );
#-----------Euler integrand for CES production
    let EULER_INTEGRAND[r, i, s] := DUAL_KAP[r, i, s] * (1 - DELTA[i]) 
      + DUAL_MKT_CLR[i, s] * (
        SHR_KAP_OUT_CES[i] * LAB_EXT[r, i, s] ** RHO_OUT
          * KAP[r, i, s] ** (RHO_OUT - 1)
        * SCALE_OUT * A[i] * SHR_DOM[r, i, s]
          * (DOM[r, i, s] / (A[i] * SHR_DOM[r, i, s]))
            ** (1 - RHO_OUT / SCALE_OUT)
        - PHI_ADJ[i] * (2 * GROWTH_KAP[r, i, s] + GROWTH_KAP[r, i, s] ** 2)
      );
    if s > PInf then 
    let EULER_RATIO[r, i, s] 
        := BETA * EULER_INTEGRAND[r, i, s] / DUAL_KAP[r, i, s - 1];
  };
  display GROWTH_KAP, GROWTH_OUT, EULER_INTEGRAND, EULER_RATIO,
    min{r in Regions, i in Sectors} E_OUT[r, i, s],
    max{r in Regions, i in Sectors} E_OUT[r, i, s],
    min{r in Regions, i in Sectors} CON[r, i, s],
    max{r in Regions, i in Sectors} CON[r, i, s],
    min{r in Regions, i in Sectors} INV_SUM[r, i, s],
    max{r in Regions, i in Sectors} INV_SUM[r, i, s],
    min{r in Regions, i in Sectors} MED_SUM[r, i, s],
    max{r in Regions, i in Sectors} MED_SUM[r, i, s],
    max{r in Regions, i in Sectors} ADJ_COST_KAP[r, i, s],
    max{i in Sectors} abs(MKT_CLR[i, s]),
    min{i in Sectors} DUAL_MKT_CLR[i, s],
    max{i in Sectors} DUAL_MKT_CLR[i, s],
    min{r in Regions, i in Sectors} LAB[r, i, s],
    max{r in Regions, i in Sectors} LAB[r, i, s],
    min{r in Regions, i in Sectors} LAB_EXT[r, i, s],
    max{r in Regions, i in Sectors} LAB_EXT[r, i, s],
    min{r in Regions, i in Sectors} KAP[r, i, s + 1],
    max{r in Regions, i in Sectors} KAP[r, i, s + 1],
    min{r in Regions, i in Sectors} GROWTH_KAP[r, i, s],
    max{r in Regions, i in Sectors} GROWTH_KAP[r, i, s],
    min{r in Regions, i in Sectors} DUAL_KAP[r, i, s],
    utility, tail_val, pres_disc_val,
    (sum{r in Regions, i in Sectors} con[r, i, LInf])
      / (sum{r in Regions, i in Sectors} E_out[r, i, LInf]),
    KAP['SEQ', 'PbSc', s] / E_OUT['SEQ', 'PbSc', s],
    s, _ampl_elapsed_time, _total_solve_time, ctime();
#  for {r in Regions, i in Sectors}{
#  if s > PInf then display EULER_RATIO[r, i, s] - EULER_RATIO[r, i, s - 1];
#  };
};


#display KAP;
