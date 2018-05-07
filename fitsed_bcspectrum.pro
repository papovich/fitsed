
; plot a bc model given an ised file, age, ebv, redshift, mass
; normalization.
; 
; also has a flag to include nebular emission lines (and continuum?
; not yet)

pro fitsed_bcspectrum, isedfile,  z, mass, tage, tebv, $
                       chabrier=chabrier, cb07=cb07, $
                       ;nebular_emission=nebular_emission, $
                       nebular_fesc=nebular_fesc, $
                       outspectrum=outspectrum, $
                       nebular_nolya=nebular_nolya, $
                       nebular_lya_offset=nebular_lya_offset,$
                       endian=endian, $
                       high_res=high_res, $
                       EXTINCTION_LAW=EXTINCTION_LAW

; outspectrum is returned and is 2D array with
; outspectrum[0,*]=wavelength in angstroms
; outspectrum[1,*]=flux f_nu in micro-Jy

if keyword_set(extinction_law) then $
  extinction_law = extinction_law $
else $
  extinction_law = 'fitsed_calzetti'

; z: redshift
; mass: mass in solar masses
; tage: age in years
; tebv: color excess

  if not keyword_set(nebular_fesc) then nebular_fesc=0.0 
  ;; default is that all LyC photons absorbed and converted to other
  ;; emission lines.
;  if not keyword_set(nebular_fesc) then nebular_fesc=0.0 else $
;     nebular_fesc=nebular_fesc
  
;  print, nebular_fesc
  ;; this default means that no emission lines are added.
  if not keyword_set(nebular_lya_offset) then nebular_lya_offset = 5000.
  ; the default here is to have lya offset by 5000 kms red.
; Other keyword: NEBULAR_NOLYA = flag.  If set to 1 then use Lya, else don't.

  fitsed_readbct,isedfile,age=bcage, spec=bcspec,color4=c4, color3=c3, $
                 chabrier=chabrier,cb07=cb07,high_res=high_res,$
                 metallicity=metalZ,endian=endian
  bciage = (where( abs(bcage[1:*] - tage) eq $
                   min(abs(bcage[1:*] - tage)) ) )[0]
  tspec = fltarr(2,n_elements(bcspec[0,*]))
  tspec[0,*] = bcspec[0,*]
  tspec[1,*] = bcspec[bciage+1,*]
  
; add nebular emission if requested:
 ; if keyword_set(nebular_emission) then begin
  if nebular_fesc lt 1.0 then begin
     print, 'adding emission lines'
     nly = double(c3.nly[bciage])
     if finite(nly) eq 0 then nly=0.d
     metalZ=float(metalZ)
     fitsed_calc_nebular, reform(tspec[0,*]), tmp, nolya=nebular_nolya, $
                          metallicity=metalZ,$
                          n_lyc=nly, $
                          lya_offset=nebular_lya_offset, $
                          continuum=nebular_continuum
     tspec[1,*] = tspec[1,*] + tmp
  endif
;  endif
  
  outspectrum = tspec
;; dust attenuation
  if 0 then   $
     for j=0,n_elements(tspec[1,*])-1 do begin
     outspectrum[1,j] = tspec[1,j] * 10^(-0.4 *  $
                                         call_function(extinction_law, tspec[0,j], tebv))
;calzetti(tspec[0,j]))
  endfor
  A_lambda = call_function(extinction_law, tspec[0,*], tebv)
  tspec[1,*] *= 10^(-0.4 * A_lambda)
  outspectrum[1,*] = tspec[1,*]
;;  stop
;; IGM attenuation:
  temp = outspectrum
  if keyword_set(USE_MADAU) then begin
     outspectrum[1,*] = fitsed_zattenuate(temp[0,*], temp[1,*], z)
  endif   else begin
     meiksin_zed = z
     if meiksin_zed gt 1.5 then $
        outspectrum[1,*] = fitsed_zattenuate(meiksin, temp[0,*], temp[1,*], meiksin_zed) $
     else outspec[1,*] = temp[1,*]
     ;; there is a known bug in meiksin
     ;; fitsed_zattenuate that if you feed a
     ;; redshift < 1.5, weird crap occurs.  DONT do this.
  endelse

  scale = (mass) * 3.839d33 / 4. / !dpi  / $
          double(fitsed_luminosity_distance_nu(z,h=!h, omega=!omega, $
                                               lambda=!lambda))^2 / $
          (3.085d24)^2  / 1d-29 / c4.mstar[bciage]


  outspectrum[1,*] = outspectrum[1,*] * outspectrum[0,*]^2 /$
                     2.998d18*scale
  outspectrum[0,*] = outspectrum[0,*]*(1+z)

return
end
  
