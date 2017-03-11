;+
; NAME: ZATTENUATE_MEIKSIN
; PURPOSE:
;     Attenuates a given spectrum by IGM at observed z=ZED. 
; EXPLANATION: 
;     Takes dz=0.1 interpolated tables of Meiksin (2006)
;     and applies IGM attenuation to a given spectrum at given redshit
; CALLING SEQUENCE:
;     ZATTENUATE_MEIKSIN(IGM_tables,wavelength_array,flux_array,redshift)
; INPUT PARAMETERS: 
;     IGM_tables - the structure "meiksin" from the file interpol_meiksin.sav
;     wavelength - BC array of wavelengths
;     spectrum - 1D spectrum array
;     zed - redshift of object
; OPTIONAL PARAMETERS:
;     verbose - tell me more
; OUTPUT PARAMETERS:
;     returns attenuated spectrum
; EXAMPLE:
;     resulting_spec = zattenuate_meiksin( meiksin, lambda, spectrum, 5.345)
;-

;; %% Read in the structure once outside the procedure to save significant time %%
;; %%                restore, '$RS/LUTFIT/interpol_Meiksin.sav'                 %%
FUNCTION FITSED_ZATTENUATE, ZSTRUCTURE, WAVELENGTH, SPECTRUM, ZED, VERBOSE=VERBOSE, silent=silent

USAGE = '% USAGE: ZATTENUATE_MEIKSIN(WAVELENGTHS, FLUXES, Z)'

if (N_params() LT 3) then begin
    print,usage
    return, 0.0
endif

if ZED lt 0.0 and not keyword_set(silent) then begin
    print,'% Error: Redshift must be >= 0.0'
    return,0.0
endif

if total(size(zstructure)) eq 0 then begin
   findpro,'interpol_Meiksin',dir=dir,/noprint
   if ~file_test(dir[0]+'/interpol_Meiksin.sav') then $
      message,'ERROR: cannot find interpol_Meiksin.sav file'
   restore,dir[0]+'/interpol_Meiksin.sav'
   zstructure = meiksin
endif

;; Read in table of interpolated IGM transmissions
;;    This contains a structure called "meiksin". I'll rename
;;    it to "zatt", ykno, for the hell of it.
zatt = zstructure
zdat=dblarr(2,8000)
zdat[0,*] = zatt.lam

;; Specify the redshift
RZED = ZED
;; Find the correct tau(z) for input z, and assign it to zdat.
if zed gt 7.05d then begin
   zdat[1,*] = reform(zatt.tau[n_elements(zatt.zed)-1])
   zdat[0,*] = (zdat[0,*]/(1.d + 7.0d) )*(1.d + RZED)
   if keyword_set(verbose) and not keyword_set(silent) then message, "redshift higher than IGM table values.. using z=7", /contin
endif else if zed lt 1.45 then begin
   zdat[1,*] = reform(zatt.tau[0])
   zdat[0,*] = (zdat[0,*]/(1.d + 1.5d) )*(1.d + RZED)
   if keyword_set(verbose) and not keyword_set(silent) then message, "redshift lower than IGM table values.. using z=1.5", /contin
endif else zdat[1,*] = reform(zatt.tau[findel(zatt.zed, rzed),*])

;; Change in to rest frame wavelengths
zdat[0,*] = zdat[0,*]/(1.0d +rzed)
lam =reform(zdat[0,*])

if n_elements(wavelength) ne n_elements(spectrum) then begin
    message, "% Error, wavelength and flux arrays must have the same number of elements",/continue
    return, 0.0
endif

;; Cycle through all wavelengths of the spectrum and attenuate each
tspectrum = spectrum
swave = where(wavelength le max(lam) and $
              wavelength ge min(lam))
swave_null = where(wavelength lt min(lam))

;; If the spectrum has wavelength < legal values in zatt.lam,
;;    then make the spectrum =0 there.
if (swave_null[0] ne -1) then $
  tspectrum[swave_null] = 0.0d

;; Value of transmitted attenuation at each BC wavelength. 
;;    Since lam is at 1A resolution, we need only to find the
;;    corresponding BC wavelength.. so no interpolation.
for i=0,n_elements(swave)-1 do begin
   ztrans = zdat[1,findel(lam, wavelength[swave[i]])]
  
;; ******* Steve F EDIT July 12 2012 ********
;; Only attenuate half of Lyman alpha line with IGM.
  if wavelength[swave[i]] eq 1215.0d0 then begin
     tau_1215=(-1.0d0)*alog(ztrans)
     tau_lya=(-1.0d0)*alog(0.5d0) - alog(1.0d0+exp(-1.0d0*tau_1215))
     ztrans=exp((-1.0d0)*tau_lya)
  endif
;;*******************************************

;; Actually apply attenuation ztrans at each wavelength step of ztrans
  tspectrum[swave[i]] = spectrum[swave[i]] * ztrans
endfor
return, tspectrum

end
