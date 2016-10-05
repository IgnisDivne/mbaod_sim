$PROBLEM    PK model
$INPUT      ID TIME DV DROP DROP AMT DROP DROP DROP AGE WT
$DATA      sim_data.csv IGNORE=@
$SUBROUTINE ADVAN1 TRANS2
$PK 

 TVCL  = EXP(THETA(1))*EXP(ETA(1))
 TVV   = EXP(THETA(2))*EXP(ETA(2))
 TM50  = EXP(THETA(3))
 HL    = 0.75
 HILL  = EXP(THETA(4))

 CL    = TVCL*(WT/70)**HL*((AGE**HILL)/(AGE**HILL + TM50**HILL))
 V     = TVV *(WT/70)
  SC    = V

$THETA  (0,1,100) ; TVCL
 (0,3,100) ; TVV
 (0,3.651,1000) ; TM50
 (0,1.5261,10) ; Maturation HILL
$OMEGA  0.05
 0.05
$SIGMA  0.015
 0.0001  FIX
$ERROR 
 IPRED = F
 Y= IPRED*(1+EPS(1))  + EPS(2)

$SIMULATION (6638008) ONLYSIM NSUBPROBLEM=1
$TABLE      ID TIME DV AMT AGE WT FILE=mc_sim_1.tab NOAPPEND NOPRINT
            NOHEADER
;$EST METHOD=1 INTER MAXEVAL=9999 SIG=3 PRINT=5 NOABORT POSTHOC
