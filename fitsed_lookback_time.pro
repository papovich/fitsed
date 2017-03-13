
function efunc_time, z
; NAME:
;       [(1+z)*E(z)]^-1 function used for time integral.  
; PURPOSE:
;       Returns the inverse value of (1+z)E(z) for all the cosmology integrals.
; EXPLANATION:
;       Returns the inverse value of (1+z)E(z), E.g., H(z) = E(z)*H(z=0)
;
; CALLING SEQUENCE:
;       a = efunc_time(z)
;
; INPUT PARAMATERS:
;       z = redshift
;       Note that !omega, !h, and !lambda should already be defined
;
; OUTPUT PARAMETERS:
;
; EXAMPLE:

on_error,0

if ( N_params() LT 1 ) then begin
    print,'% Syntax - EFUNC_TIME(z)'
    return,0.0
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

defsysv, '!omega', exists=omega_exists
defsysv, '!lambda', exists=lambda_exists
defsysv, '!h', exists=h_exists
if (h_exists eq 0 and lambda_exists eq 0 and omega_exists eq 0) then begin
    print,'% ERROR: !omega, !h, and !lambda must be defined'
    return,0.0
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

omegak = 1.0 - !omega - !lambda
return, (!omega*(1.0+z)^3.0 + omegak*(1.0+z)^2.0 + !lambda)^(-0.5) / (1.0+z)

end



function fitsed_Lookback_Time, z1, z2, h=h, omega=omega, lambda=lambda
; NAME:
;       Lookback_Time in years
; PURPOSE:
;       Returns the Lookback_Time between z=z1 and z=z2
; EXPLANATION:
;
; CALLING SEQUENCE:
;       a = Lookback_Time(z1,z2,h=h,omega=omega,lambda=lambda)
;
; INPUT PARAMATERS:
;       z1 = redshift
;       z2 = redshift
;       h = hubble parameter (H0 / 100)
;       omega = dimensionless mass density parameter
;       lambda = dimensionless cosmological constant parameter
; OUTPUT PARAMETERS:
;
; EXAMPLE:

on_error,0

if ( N_params() LT 2 ) then begin
    print,'% Syntax - Lookback_Time(z1,z2,h=h,omega=omega,lambda=lambda)'
    return,0.0
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if keyword_set(h) then h=h $
else h=1.0
if keyword_set(omega) then omega=omega $
else omega=1.0
if keyword_set(lambda) then lambda=lambda $
else lambda=0.0

defsysv, '!omega', exists=omega_exists
if (omega_exists eq 0) then defsysv,'!omega', omega $
else !omega=omega

defsysv, '!lambda', exists=lambda_exists
if (lambda_exists eq 0) then defsysv,'!lambda', lambda $
else !lambda=lambda

defsysv, '!h', exists=h_exists
if (h_exists eq 0) then defsysv,'!h', h $
else !h=h

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

hub_distance = 2.9979e5 / 100.0 / !h

omegak = 1.0 - !omega - !lambda

hub_time = 1.0 / 100.0 / !h * 3.085e19 / 3600.0 / 24.0 / 365.25; in years

if size(/dim, z1) gt 1 then begin
   lbt = fltarr(n_elements(z1))
   for i=0,n_elements(lbt)-1 do begin
      use1 = z1[i]
      use2 = z2
      if z1[i] gt z2 then begin
         use2 = z1[i]
         use1 = z2
      endif
      lbt[i] = hub_time * qromb('efunc_time',use1,use2)
   endfor
endif else begin
   if (z1 gt z2) then begin
      temp = z2
      z2 = z1
      z1 = temp
   endif
   lbt = hub_time * qromb('efunc_time',z1,z2)
endelse

return, lbt

end

