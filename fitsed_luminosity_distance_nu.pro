


function efunc, z
; NAME:
;       [E(z)]^-1 function used for cosmology integrals.  
; PURPOSE:
;       Returns the inverse value of E(z) for all the cosmology integrals.
; EXPLANATION:
;       Returns the inverse value of E(z), E.g., H(z) = E(z)*H(z=0)
;
; CALLING SEQUENCE:
;       a = efunc(z)
;
; INPUT PARAMATERS:
;       z - redshift
;       Note that !omega, !h, and !lambda should already be defined
;
; OUTPUT PARAMETERS:
;
; EXAMPLE:

on_error,0

if ( N_params() LT 1 ) then begin
    print,'% Syntax - EFUNC(z)'
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
return, (!omega*(1.0+z)^3.0 + omegak*(1.0+z)^2.0 + !lambda)^(-0.5)

end



function proper_distance, z, z1, h=h, omega=omega, lambda=lambda
; NAME:
;       Proper_Distance
; PURPOSE:
;       Returns the Proper Distance between two nearby objects
; EXPLANATION:
;       Returns the distance between two objects in the Universe which
;       remains constant with epoch if the two objects are moving with
;       the Hubble flow.
;
; CALLING SEQUENCE:
;       answer = Proper_Distance(z,h=h,omega=omega,lambda=lambda)
;
; INPUT PARAMATERS:
;       z - redshift
;       h = hubble parameter (H0 / 100)
;       omega = dimensionless mass density parameter
;       lambda = dimensionless cosmological constant parameter
; OUTPUT PARAMETERS:
;
; EXAMPLE:

on_error,0

if ( N_params() LT 1 ) then begin
    print,'% Syntax - Proper_Distance(z,h=h,omega=omega,lambda=lambda)'
    return,0.0
endif

if n_params() lt 2 then z1=0.0

if z lt z1 then begin 
    tmp=z
    z = z1
    z1=tmp
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

return, hub_distance*QROMB('efunc',z1,z)

end


function comoving_distance, z, z1, h=h, omega=omega, lambda=lambda
; NAME:
;       Comoving_Distance
; PURPOSE:
;       Returns the Transverse Comoving Distance between two nearby objects
; EXPLANATION:
;       The distance between two objects at the same redshift, but
;       separated on the sky by some angle d(theta) is:
;              Comoving_Distance*d(theta)
;
; CALLING SEQUENCE:
;       a = Comoving_Distance(z,h=h,omega=omega,lambda=lambda)
;
; INPUT PARAMATERS:
;       z - redshift
;       h = hubble parameter (H0 / 100)
;       omega = dimensionless mass density parameter
;       lambda = dimensionless cosmological constant parameter
; OUTPUT PARAMETERS:
;
; EXAMPLE:

on_error,0

if ( N_params() LT 1 ) then begin
    print,'% Syntax - Comoving_Distance(z,h=h,omega=omega,lambda=lambda)'
    return,0.0
endif

if n_params() lt 2 then z1=0.0

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

if omegak gt 0.0 then $
  return, hub_distance / sqrt(omegak) * $
  sinh( sqrt(omegak) * proper_distance(z,z1, h=!h,omega=!omega,lambda=!lambda)$
        / hub_distance )
if omegak eq 0.0 then $
  return, proper_distance(z,z1, h=!h,omega=!omega,lambda=!lambda)
if omegak lt 0.0 then $
  return, hub_distance / sqrt(-omegak) * $
  sin( sqrt(-omegak) * proper_distance(z,z1, h=!h,omega=!omega,lambda=!lambda)$
       /hub_distance )

end




function fitsed_luminosity_distance_nu, z, h=h, omega=omega, lambda=lambda
; NAME:
;       Luminosity_Distance
; PURPOSE:
;       Returns the Luminosity Distance for flux in f_nu units
; EXPLANATION:
;
; CALLING SEQUENCE:
;       a = Luminosity_Distance_nu(z,h=h,omega=omega,lambda=lambda)
;
; INPUT PARAMATERS:
;       z = redshift
;       h = hubble parameter (H0 / 100)
;       omega = dimensionless mass density parameter
;       lambda = dimensionless cosmological constant parameter
; OUTPUT PARAMETERS:
;
; EXAMPLE:

on_error,0

if ( N_params() LT 1 ) then begin
    print,'% Syntax - Luminosity_Distance_nu(h=h,omega=omega,lambda=lambda)'
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

if size(/dim, z) gt 1 then begin
   ld = dblarr(n_elements(z))
   for i=0,n_elements(z)-1 do begin
      ld[i] =  comoving_distance(z[i],h=!h,omega=!omega,lambda=!lambda) * sqrt(1.0+z[i])
   endfor
   return, ld
endif else begin
   return, comoving_distance(z,h=!h,omega=!omega,lambda=!lambda) * sqrt(1.0+z)
endelse

end

