;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; This IDL script returns the calzetti k(l) for a given l, given in
;; angstroms 
;; Calzetti et al 2000
;; 
;; returns A(lambda) which gives the extinction (in mag) given an
;; E(B-V), i.e., A(lambda) = k(lambda) * E(B-V).
;;
;; THIS IS THE EXTINCTION ON THE STELLAR CONTINUUM.
;; TO GET NEBULAR EXTINCTION, MULTIPLY BY 1.0/0.44

FUNCTION FITSED_CALZETTI, L, EBV,  GAS=GAS

;;
;; USAGE:
USAGE = 'X = CALZETTI(WAVELENGTH, EBV,  GAS=GAS)'
;;

if EBV eq 0 then return, 0.0

IF N_PARAMS(0) LT 1 THEN BEGIN
    PRINT,'% USAGE: '+USAGE
    RETURN,0.0
ENDIF

LL = L / 1.e4
;print,ll

if total(size(ll,/dim)) le 1 then begin
   IF LL GT 0.63 THEN CAL = ( 2.659*((-1.857) + 1.040d0/LL ) + 4.045) $
   ELSE CAL = 2.659*( -2.156 + 1.509/LL - 0.198/LL^2 + 0.011/LL^3 ) + 4.045
endif else begin
   CAL = LL*0.0
   tLL = where( LL gt 0.63, comp=cLL)
   if tLL[0] ne -1 then $
      CAL[tLL] = ( 2.659*((-1.857) + 1.040d0/LL[tLL] ) + 4.045)
   if cLL[0] ne -1 then $
       CAL[cLL] = 2.659*( -2.156 + 1.509/LL[cLL]- 0.198/LL[cLL]^2 + 0.011/LL[cLL]^3 ) + 4.045
endelse

if keyword_set(GAS) then CAL=CAL * 0.44

RETURN, CAL*EBV

END

