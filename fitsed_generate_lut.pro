function int_spectrum_noi, spectrum, bandpass, zed=zed, $
                           flambda=flambda, angstrom=angstrom
;
; Returns a spectrum integrated with a bandpass.  The spectrum is
; assumed to be given as an array of two columns correspond to the
; wavelength (in microns) and the flux density (in flux per unit Hz).
; If the flambda flag is set then the flux density is assumed to be
; flux per unit angstrom; if the angstrom flag is set, then the
; wavelengths are assumed to be in angstroms.
;
; If zed is given, then the spectrum is redshifted to z=zed and then
; integrated through the bandpass
;
; calls interp.pro
;

; this version (_noi) does *not* interpolate the bandpass.  Rather it
; assumes that the bandpass is already aligned in wavelength space.

on_error,2

if ( N_params() LT 2) then begin
    print,'% Syntax: ans = int_spectrum(spec,bandpass,<zed=z>,</flambda>,</angstrom>)'
    return,0.0
endif

ssz = size(spectrum)
bsz = size(bandpass)

if (ssz[0] lt 2 or bsz[0] lt 2) then begin
    print,'% Error: spectrum and bandpass arrays must have >=2 columns'
    return,0.0
endif

if KEYWORD_SET(zed) then zed=zed else zed=0.0

lambda = spectrum[0,*]*1.e4
flux = spectrum[1,*]
tlambda = bandpass[0,*]*1.e4
ttrans = bandpass[1,*]
xxx = where( ttrans lt 0.0)
if xxx[0] ne -1 then ttrans[xxx] = 0.0

;lambda = lambda * (1.0+zed)

if KEYWORD_SET(angstrom) then begin
    lambda=lambda/1.e4
    tlambda=tlambda/1.e4
endif

if KEYWORD_SET(flambda) then flux = lambda^2.0 / 2.998e18 * flux

norm=0.0
solution=0.0

tflux = flux
tnu = 2.998d18 / tlambda


if not keyword_set(transcut) then transcut=0.0
y=where( ttrans gt transcut*max(ttrans))
;solution = int_tabulated( (tnu[y]), (tflux[y]*ttrans[y]/tnu[y]),/double) / $
;           int_tabulated( (tnu[y]), (ttrans[y]/tnu[y]),/double)
solution = tsum( tnu[y], tflux[y]*ttrans[y]/tnu[y]) / $
           tsum( tnu[y], ttrans[y]/tnu[y])
; 
; checkm0=0 & checkm1=1
; for i=0,n_elements(ttrans)-2 do begin
;     lam0 = tlambda[i] & lam1 = tlambda[i+1]
;     flux0 = flux[i]
;     flux1 = flux[i+1]
; ;    flux0 = interp(lambda,flux,tlambda[i],m=checkm0)
; ;    flux1 = interp(lambda,flux,tlambda[i+1],m=checkm1)
; ;    print,i,checkm0,tlambda[i],lambda[checkm0],lambda[checkm0+1],$
; ;          flux0,flux[checkm0],flux[checkm0+1]
; ;    flux0 = flux0 * (lam0/ (1.0+zed))^2.0
; ;    flux1 = flux1 * (lam1 / (1.0+zed))^2.0
;     T0 = ttrans[i]/lam0
;     T1 = ttrans[i+1]/lam1
;     
;     norm = norm + (lam1-lam0)*T0 + 0.5*(T0-T1)*(lam0-lam1)
;     solution = solution + (lam1-lam0)*(flux0*T0) + $
;       0.5*(T0*flux0 - T1*flux1)*(lam0-lam1)
;     
; endfor

;solution=solution/norm

return,solution

end



function interp, xvec, yvec, x, m=m, quiet=quiet
;
; *Simple* linear interpolation of array yvec at point x
; Nothing fancy, just lineary interpolates between adjacent points.
; Function returns the interpolated value of yvec at point x 
; m is the index where interpolation occured
;

if ( N_params() LT 3) then begin
    print,'% Syntax: ans = interp(xvector, yvector, x)'
    return,0.0
endif

if (n_elements(xvec) ne n_elements(yvec)) then begin
    print,'% Error: xvector and yvector must have equal # of elements'
    return,0.0
endif

nmax=n_elements(xvec)

if (x lt xvec[0] or x gt xvec[nmax-1]) then begin
    if not keyword_set(quiet) then begin
        print,'% Error: x out-of-bounds {x < xvec[0] or x > xvec[N]}'
        print,'% ',x,' < ',xvec[0],' or ',x,' > ',xvec[nmax-1]
    endif
    return,0.0
endif

m=0
for i=0,nmax-1 do begin
;    print,i,m,x,xvec[i]
    if (x ge xvec[i]) then m=i $
    else break
endfor

if (xvec[m] eq x) then solution=yvec[m] $
else begin
    y0=yvec[m] & y1=yvec[m+1]
    x0=xvec[m] & x1=xvec[m+1]
    
    solution = (y0-y1) / (x0-x1) * (x-x1) + y1
endelse

return,solution

end




FUNCTION INTERPOLATE_FILTER, SPEC, FILTERS, Z

; this routine takes a structure filled with FILTER transmission
; curves, and interpolates them to match the wavelengths in SPEC for a
; given Z.  It then returns a new structure with the interpolated
; filters

; define the tag names for each filter:
  tags = strarr(n_tags(filters))
  for i=0,n_elements(tags)-1 do begin
     tags[i] = 'field'+string(FORMAT=('(I1)'),i)
     if i ge 10 then $
        tags[i] = 'field'+string(FORMAT=('(I2)'),i)
     if i ge 100 then $
        tags[i] = 'field'+string(FORMAT=('(I3)'),i)
  endfor

; define output filters

  for ifilt = 0, n_tags(filters)-1 do begin

     tfilt = spec
     for i=0,n_elements(spec[0,*])-1 do begin
        tfilt[1,i] = interp(FILTERS.(ifilt)[0,*]/(1+z),$
                            FILTERS.(ifilt)[1,*],$
                            spec[0,i],/quiet)
     endfor
     
     if ifilt eq 0 then $
        newfilters=create_struct(tags[ifilt], tfilt) $
     else $
        newfilters=create_struct(newfilters, tags[ifilt], tfilt)

  endfor

  return,newfilters

end
    
    


;------------------------------------------------------------
;
;

PRO FITSED_GENERATE_LUT, OUTPUTFILE, $
; stuff to read in models:
                         SSP=SSP, $ ; 'bc03','cb07','bc11',... BC03(default), CB07, BC11
                         IMF=IMF, $ ; 'chab','salp'(default), for now
                         MODELDIR=MODELDIR,$ ; where models are stored
                         _EXTRA=_EXTRA, $ ;
                         ;; _EXTRA passed to readbct and can include
                         ;;        any of the following:
                         ;; HIGH_RES=HIGH_RES, $
                         ;; ENDIAN=ENDIAN,$ ; assumed ENDIAN='native'
                             ;; by default (need to check!) 
                         ;; CB07=CB07, $
                         ;; in addition, this is SET by IMF: 
                         ;; CHABRIER=CHABRIER,$
; filter information:
                         BANDPASS_FILE=BANDPASS_FILE, $ ; .RES FILE
                         BANDPASS_INDEXES=BANDPASS_INDEXES, $ ; array of indexes in .RES file for filters (needs to be same order as in .cat catalog
                         CATALOG=CATALOG, $ ; header of [CATHDR].cat  and .translate file
;                        either CATALOG or bandpass_indexes must be
;                        specified, if both are specifed, code will
;                        use CATALOG. 
;
; informaiton for the grid set up:
                         ZMIN=ZMIN, $
                         ZMAX=ZMAX, $
                         ZSTEP=ZSTEP,$
                         ZLOG=ZLOG, $ ; if set, then code assumes ZSTEP is delta(log z), ZED will still be redshift (but spacing will be log)
                         ZED=ZED, $
;
                         EXTINCTION_LAW=EXTINCTION_LAW, $
                         EBVMIN=EBVMIN, $
                         EBVMAX=EBVMAX, $
                         EBVSTEP=EBVSTEP, $
                         EBVARR=EBVARR, $
; OR: 
                         NODUST=NODUST, $ ; sets ebvmin=ebvmax=0

                         ;; second (delta) parameter for extinction
                         ;; law, not always needed, but allowed for
                         DELTAMIN=DELTAMIN, $
                         DELTAMAX=DELTAMAX, $
                         DELTASTEP=DELTASTEP, $
                         DELTAARR=DELTAARR, $ 
;
                         LOG_AGEMIN=LOG_AGEMIN, $
                         LOG_AGEMAX=LOG_AGEMAX, $
                         LOG_AGESTEP=LOG_AGESTEP, $
                         LOG_AGEARR=LOG_AGEARR, $
; 
                         TAUARR=TAUARR, $ ; input the array that you want (must match SSP values)
; 
                         METALARR=METALARR, $ ; input the array that you want (must match SSP values)
;
; stuff for nebular emission
                         NEBULAR_FESC=NEBULAR_FESC, $
                         NEBULAR_CONTINUUM=NEBULAR_CONTINUUM, $
                         NEBULAR_NOLYA=NEBULAR_NOLYA,$
                         NEBULAR_LYA_OFFSET=NEBULAR_LYA_OFFSET,$
; 
; other control stuff: 
                         RESTFRAME=RESTFRAME, $
; if restframe set, then filters are assumed to be in the rest-frame
                         USE_MADAU=USE_MADAU, $ ; otherwise uses Meiksin
                         VERBOSE=VERBOSE


; see readbct notes for the above, it controls byte swapping

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

; CALLING SEQUENCE:
;       GENERATE_LUT(<isedfile>,<outputfile>,<array of bandpass filenames>)
;
; INPUT PARAMETERS:
;
; OUTPUT PARAMETERS:
;
; EXAMPLE:

  USAGE = '% USAGE: GENERATE_LUT, OUTFILE,...'

  if (N_params() LT 1) then begin
    print,usage
    return
 endif
  
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; ~~~ Math underflow errors are turned off ~~~  
  currentExcept= !Except
  !Except=0
  void = Check_Math()
  floating_point_underflow = 32
  status = Check_Math()
  IF(status AND NOT floating_point_underflow) NE 0 THEN $
     Message, 'IDL Check_Math() error: ' + StrTrim(status, 2)
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


; SET WHERE MODELS ARE LOCATED AND SET UP STRINGs:
  if not keyword_set(MODELDIR) then $
     modeldir='./' ;; probably superceded by the stuff below:
  ;
  if not keyword_set(SSP) then SSP='bc03'
  ;modeldir=modeldir ; getenv('HOME')+'/Source/bc/bc03/mymodels/'
  if not keyword_Set(IMF) then IMF='salp'
  modelHdr = modeldir+imf
  modelTail='.ised'


; SET EXTINCTION LAW
  if keyword_set(extinction_law) then $
     extinction_law = extinction_law $
  else $
     extinction_law = 'fitsed_calzetti' 
;; this should be the name of the IDL routine that will be called as extinction_law(wavelength)

; SET NEBULAR EMISSION: 
  ;if not keyword_set(nebular_fesc) then nebular_fesc=1.0
  if not keyword_set(nebular_continuum) then nebular_continuum = 0 

; place Lya at systematic redshift ?  Hell no.  Use 500 km/s as the
; default based on Shapley et al. 2003, ApJ, 588, 65
  if not keyword_set(nebular_lya_offset) then nebular_lya_offset = 500.
; positive nebular_lya_offset is redshift

;------------------------------------------------------------
;
; READ IN THE FILTERS: 
  filters = fitsed_read_filters(bandpass_file, $
                                lambda=lambda, $
                                catalog=catalog, $
                                bandpass_indexes=bandpass_indexes, $
                                VERBOSE=VERBOSE)

; get the number of bands:
  sz_bands = n_tags(filters)
; get the name of the filters:
  nm_bands = tag_names(filters)

;------------------------------------------------------------
;
; INITIALIZE AND BUILD DATA ARRAYS
  if not keyword_set(log_agemin) then log_agemin=5.0
  if not keyword_set(log_agemax) then log_agemax=10.0
  if not keyword_set(log_agestep) then log_agestep=0.2
  LOG_AGEARR= findgen( (log_agemax-log_agemin)/log_agestep + 1) * log_agestep + $
            log_agemin

;; create metallicity sting for models:
  metalStr = strarr(n_elements(metalArr))
  for i=0,n_elements(metalArr)-1 do begin
     case metalArr[i] of
        0.02 : metalStr[i] = 'zsol'
        0.008 : metalStr[i] = 'z0p4sol'
        0.004 : metalStr[i] = 'z0p2sol'
        0.0004 : metalStr[i] = 'z0p02sol'
        0.05 : metalStr[i]='z2p5sol'
        else  : begin
           message, 'ERROR, metalArr[i] value of '+string(metalArr[i])+' not recognized, stopping'
        end
     endcase
  endfor

;; create tau string for the tau models
  TAUSTR=strarr(n_elements(tauArr))
  for i=0,n_elements(tauArr)-1 do begin
     case tauArr[i] of
        100 : taustr[i] = '100g'
        10 : taustr[i] = '10g'
        1 : taustr[i] = '1g'
        0.1 : taustr[i] = '100m'
        0.01 : taustr[i] = '10m'
        0.001 : taustr[i] = '1m'
        2 : taustr[i] = '2g'
        0.2 : taustr[i] = '200m'
        0.02 : taustr[i] = '20m'
        3 : taustr[i] = '3g'
        0.3 : taustr[i] = '300m'
        0.03 : taustr[i] = '30m'
        5 : taustr[i] = '5g'
        0.5 : taustr[i] = '500m'
        0.05 : taustr[i] = '50m'
        7 : taustr[i] = '7g'
        0.7 : taustr[i] = '700m'
        0.07 : taustr[i] = '70m'
        else : begin
           message, 'ERROR, tauArr[i] value of '+string(tauArr[i])+' not recognized, stopping'
        end
     endcase
  endfor
  
;; Set up E(B-V) [EBV] array:
  if keyword_set(ebvmin) then ebvmin = ebvmin else ebvmin=0.0
  if keyword_set(ebvmax) then ebvmax = ebvmax else ebvmax=1.0
  if keyword_set(ebvstep) then ebvstep = ebvstep else ebvstep=0.05
  if keyword_set(nodust) then begin 
     ebvmax = 0.0
     ebvmin = 0.0
  endif
  EBVARR = findgen( (ebvmax-ebvmin)/ebvstep + 1) * ebvstep + ebvmin

;; DELTA ALLOWED FOR: 
;; -- EASY SOLUTION IS TO SET IT UP TO USE SALMON+16 DUST LAW WHERE DELTA
;; -- IS SET BY EBV
  if not keyword_set(deltamin) then deltamin=0.0
  if not keyword_set(deltamax) then deltamax=0.0
  if not keyword_set(deltastep) then deltastep=0.1
  deltaArr=[0.0]                   ; for now DELTA DELTA DELTA TO DO

; SET UP REDSHIFT ARRAY:
  if keyword_set(zmin) then zmin = zmin else zmin=0.0
  if keyword_set(zmax) then zmax = zmax else zmax=5.0
  if keyword_set(zstep) then zstep = zstep else zstep=0.05
  if keyword_set(zlog) then begin
     zed=zmin
     while max(zed) lt zmax do zed = [zed, max(zed)+zstep*(1+max(zed))]
;     zed = findgen( (alog10(zmax)-alog10(zmin))/zstep + 2) * zstep + alog10(zmin)
;     zed=10^zed
  endif else $
     zed = findgen( (zmax-zmin)/zstep + 1) * zstep + zmin

; set array for luminosity distance (squared):
  LDSQ = double( fitsed_luminosity_distance_nu(zed, h=!h, omega=!omega, $
                                               lambda=!lambda) *$
                 1d6*3.0857d18 )^2

; INITIALIZED (L)OOK-(U)P (T)ABLE
  LUT = fltarr( n_elements(zed),  $
                n_elements(log_ageArr), $
                n_elements(ebvArr), $
                n_elements(deltaArr), $
                n_elements(metalArr),$
                n_elements(tauArr),$
                sz_bands )
  ;; lutMstar and lutSFR hold values from the c4 file on mstar and sfr
  lutMstar = fltarr( n_elements(log_ageArr), $
                     n_elements(metalArr),$
                     n_elements(tauArr))
  lutSFR =  fltarr( n_elements(log_ageArr), $
                     n_elements(metalArr),$
                     n_elements(tauArr))

  if keyword_set(verbose) then timestart=systime(/seconds)

; load in IGM unless Madau is requested (TO DO?)
;  restore,'$HOME/idl/pro/interpol_Meiksin.sav', VERB=VERBOSE

  for imetal=0,n_elements(metalArr)-1 do begin
     for itau=0, n_elements(tauArr)-1 do begin
        
      ;; define filename to model for this tau and metallicity:
        isedfile = modelHdr+'_'+metalStr[imetal]+'_tau'+taustr[itau]+'yr.ised'

        if not file_test(isedfile, /read) then begin
           message,'% ERROR, File '+isedfile+$
                   ' does not exist or is unreadable, stopping'
        endif
        if imetal eq 0 and itau eq 0 then isedfiles=[isedfile] $
        else isedfiles=[isedfiles,isedfile]

; ================================================
; Track TIME
        if keyword_set(verbose) then FITSED_CHECKTIME, timestart,isedfile
; ================================================
        fitsed_readbct, isedfile, age=bcage, spec=spec, metallicity=metalZ, $
                        color3=c3, color4=c4, _EXTRA=_EXTRA
      
           
        ;; define temporary spectra
        tspec = fltarr(2,n_elements(spec[0,*]))
        tspec[0,*] = spec[0,*]

        zspec = fltarr(2,n_elements(spec[0,*]))
        zspec[0,*] = spec[0,*]

        ;; AGE LOOP:
        for iage=0, n_elements(log_ageArr)-1 do begin
         
           ;; find the index of the closest age: 
           age_index = findel(bcage,10d^log_ageArr[iage])

           lutMstar[iage, imetal, itau] = c4.mstar[age_index]
           lutSFR[iage, imetal, itau] = c4.sfr[age_index]
           
           ;; EBV LOOP: 
           for iebv=0,n_elements(ebvArr)-1 do begin
              for idelta=0,n_elements(deltaArr)-1 do begin
                 
                 tspec[1,*] = spec[age_index+1,*]
; ------------------------------------------------------------
; NEBULAR EMISSION: 
                 ;; Add nebular emission if requested:
                 if nebular_fesc lt 1.0 then begin
                    nly = double(c3.nly[age_index])
                    if nly eq !VALUES.F_INFINITY then nly=0.d
                    metalZ=float(metalZ) ; metallicity taken from ISED FILE
                    fitsed_calc_nebular, $
                       reform(tspec[0,*]), tmp, nolya=nebular_nolya, $
                       metallicity=metalZ,$
                       n_lyc=nly, $
                       lya_offset = nebular_lya_offset, $
                       continuum=nebular_continuum
                    tspec[1,*] = tspec[1,*] + tmp
                 endif
;------------------------------------------------------------
; APPLY EXTINCTION LAW: (need to edit for SALMON+16, etc.
                 A_lambda = call_function(extinction_law, tspec[0,*], ebvArr[iebv])
                 tspec[1,*]*=10d^(-0.4*A_lambda)
                  
;------------------------------------------------------------
                    
                 ;; ITERATE OVER ALL REDSHIFT VALUES: 
                 for izz=0,n_elements(zed)-1 do begin
                    z=zed[izz]
                    ;; ATTENUATE the spectrum with cosmic transmission
                    if keyword_set(USE_MADAU) then begin
                       zspec[1,*] = fitsed_zattenuate(tspec[0,*], tspec[1,*], z)
                    endif   else begin
                       meiksin_zed = z
                       if meiksin_zed gt 1.5 then $
                          zspec[1,*] = fitsed_zattenuate(meiksin, tspec[0,*], tspec[1,*], meiksin_zed) $
                       else zspec[1,*] = tspec[1,*]
                       ;; there is a known bug in meiksin
                       ;; fitsed_zattenuate that if you feed a
                       ;; redshift < 1.5, weird crap occurs.  DONT do this.
                    endelse

                    ;------------------------------------------------------------
                    ; LOOP OVER ALL BANDPASSES:
                    for ibands=0, sz_bands - 1 do begin
                       ; if this is the first time, then interpolate the filters so that they
                       ; match the ised files at each wavelength
                       ; 
                       ; int_filter does interpolate, but this is supposed to speed it up.  
                       ; the interpolateion should be the same for all ised files at the same
                       ; redshift:
                       if $  ;; first time then do this:
                          itau eq 0 and $
                          imetal eq 0 and $
                          idelta eq deltamin and $
                          iebv eq 0 and $
                          iage eq 0 then begin
                          if keyword_set(verbose) and ibands eq 0 then $
                             message,/cont, $
                                     "Interpolating filters at z="+$
                                     strtrim(string(z),2)
                          if keyword_set(restframe) then begin
                             tinterpfilters = $
                                interpolate_filter(zspec,filters,0.0)
                          endif else begin
                             tinterpfilters = $
                                interpolate_filter(zspec,filters,z)
                          endelse

                          if z eq zmin then begin
                              interpfilters = $
                                 replicate(tinterpfilters, n_elements(zed))
		          endif
                          interpfilters[izz] = tinterpfilters
                       endif


                   
                    
; integrate the spectrum with the bandpass
                             
                       if keyword_set(restFrame) then begin
                          aveflx = $
                             int_spectrum_noi( zspec, $
                                               interpfilters[izz].(ibands), $
                                               zed=0.0, /flambda, /angstrom)
                       end else begin
                          aveflx = $
                             int_spectrum_noi( zspec, $
                                               interpfilters[izz].(ibands), $
                                               zed=z, /flambda, /angstrom)
                       endelse
                       ;; at this point aveflx has the units: [ Lsol/Hz/Msol]
                       
                       
; the get a Mass from this, do:
; Mass = 934.83 * (luminosity_distnace_nu(z))^2 * FLUX / aveflux
; where FLUX Is the actual photometry.  The ratio of FLUX/aveflux is
; actually just the scale factor from the SED-fitting part
                       
; Store in the output array
; and Convert LUT into units of f_nu: nJy/m_solar
                       LUT[izz, $
                           iage, $
                           iebv, $
                           idelta, $
                           imetal,$
                           itau, $
                           ibands] = $
                             aveflx * 3.839d33/(4*!PI*ldsq[izz]) / 1d-32 / $
                              double(lutMstar[iage, imetal, itau]) 
                    endfor      ; bandpass loop
                 endfor         ; redshift loop
              endfor            ; delta loop
           endfor               ; EBV loop
           print, "                        Pew Pew!! I just did Age #"+string(iage+1)+" of "+string(n_elements(log_ageArr))
        endfor                  ; age loop
        print, "                        Well, shucks.. I'm only working on Tau #"+string(itau+1)+string(n_elements(tauArr))
     endfor                     ; TAU loop
     print, "                        Cmon fella.. I JUST finished Metal #"+string(imetal+1)+string(n_elements(metalArr))
  endfor                        ; METAL loop
  
  save,filename=outputfile,$
       LUT, zed, log_ageArr, ebvArr, tauArr, metalArr, deltaArr, $
       lutMstar, lutSFR, $
       catalog, bandpass_file, filters, lambda, $
       isedfiles

  if keyword_set(verbose) then fitsed_checktime,timestart,'all done!'

  RETURN
END
