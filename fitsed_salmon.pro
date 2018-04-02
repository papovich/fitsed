;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; This IDL script returns the calzetti k(l) for a given l, given in
;; angstroms 
;; Salmon et al 2016
;; 
;; returns A(lambda) which gives the extinction (in mag) at a
;; wavelength lambda, given an E(B-V): 
;; 
;; A(lambda) = k(lambda) * E(B-V) (lambda/lambda_V AA)^delta.
;; 
;; where lambda_V=5500 AA (?) and delta = 0.62 * alog10( E(B-V))  + 0.26

;;

FUNCTION FITSED_SALMON, L, EBV, delta=delta

;;
;; USAGE:
USAGE = 'X = SALMON(WAVELENGTH, EBV)'
;;

IF N_PARAMS(0) LT 2 THEN BEGIN
    PRINT,'% USAGE: '+USAGE
    RETURN,0.0
ENDIF

if EBV le 0.0 then return, 0.0

A_lambda = fitsed_calzetti( L , EBV ) 
DELTA = 0.62 * alog10(EBV) + 0.26
return, A_lambda * (l / 5500.)^delta

END

