$PROBLEM PKPD model
$INPUT ID TIME DV AMT
$DATA est.dat IGNORE=@
$SUBROUTINE ADVAN2 TRANS2



$PK

CL   = THETA(1) * EXP(ETA(1))
V    = THETA(2) * EXP(ETA(2))
KA   = THETA(3);* EXP(ETA(3))
SC    = V

$ERROR
BASE = THETA(4) ;* EXP(ETA(4))
EMAX = THETA(5) * EXP(ETA(3))
EC50 = THETA(6)* EXP(ETA(4))
GAM  = THETA(7) ;* EXP(ETA(7))

CONC = F
IPRED = BASE + ((EMAX*CONC**GAM)/(EC50**GAM + CONC**GAM))
Y= IPRED*(1+EPS(1))  + EPS(2)

$THETA
(0.15) FIX  ; 1. CL
(8)    FIX  ; 2. V
(1)    FIX  ; 3. Ka
(0,1.5)       ; 4. BASE
(0, 150)    ; 5. EMAX
(0, 10.5)      ; 6. EC50
(0, 3.5)      ; 7. Hill



$OMEGA
0.07  FIX  ; 1 CL
0.02  FIX  ; 2 V
;0.6   FIX  ; 3 KA
;0.0625    ; 4 Base
0.0625       ; 5 EMAX
0.0625        ; 6 EC50
;0.0625        ; 7 EC50



$SIGMA
 0.015
 0.0001 FIX

$TABLE ID TIME DV AMT FILE=mc_est_1.tab NOAPPEND NOPRINT NOHEADER
$EST METHOD=1 INTER MAXEVAL=9999 SIG=3 PRINT=5 NOABORT POSTHOC
$COV