
; This script fits photometry in phot with uncertainties dphot to the
; LUT generated from BC models.
;
; Revision history:  Written in Aug 2003, revised Mar 2005.
; revised mar 2005 to use monte carlo'ing of data to construct errors
; revised dec 2004 to calculate errors fro cumulative PDF
; revised jan 2005 to use min chisq for best value (instead of max
; pdf) and to calculate instantaneous SFR for each model. 
; revised apr 2009 to allow for "free age" fitting of 2nd component
; with constant sf
; revised Oct 2011 to take out scaling.  LUT is in units of nJy per
; solar mass.   For flux densities in nJy, "scale" is just mass now....
; revised Aug 2012 to correct PDF analysis using TSUM instead of TOTAL
; MAJOR REWRITE 2016 August.  Cleaning up, added delta, etc.
; UPDATE 2019 February.  Removed delta.  Added alpha/beta for SFH=dpl.  

;------------------------------------------------------------
;------------------------------------------------------------
; 
; functions called by fitchisq: 

function old_calc_mass, mass, pdf, mst=mst, $
                    mass_best = mass_best, bestIndex=bestIndex, $
                    mass_mean = mass_mean, meanIndex=meanIndex, $
                    mass_sigma=mass_sigma, $
                    mass_median=mass_median, medianIndex=medianIndex, $
                    mass_lo68=mass_lo68, mass_hi68=mass_hi68, $
                    mass_lo95=mass_lo95, mass_hi95=mass_hi95, $
                    lo68Index=lo68Index, hi68Index=hi68Index, $
                    lo95Index=lo95Index, hi95Index=hi95Index

  mst = sort(mass)
  ttt = where(~finite(mass),/null)
  mass[ttt]=0.0 & pdf[ttt]=0.0  ; just to be sure... zero out infinite/NaNs
  y_mass = pdf / tsum(mass[mst], pdf[mst])

  mass_mean = tsum(mass[mst], mass[mst]*y_mass[mst])
  meanIndex = (where(abs(mass_mean - mass) eq $
                      min(abs(mass_mean - mass))))[0]
  mass_sigma = sqrt( tsum(mass[mst],$
                          abs(mass[mst]-mass_mean)^2.d*y_mass[mst]) )
  
  sty_mass = y_mass[mst]
  cumulative_pdf, mass[mst], y_mass[mst], stmass_median, $
                  stmass_lo68, stmass_hi68, stmass_lo95, stmass_hi95, mcum, $
                  mass_median, mass_lo68, mass_hi68, mass_lo95, mass_hi95, $
                  badobj2=badobj2
  if total(size(badobj2)) ne 0 then if badobj2 eq 1 then badobj=1
  if stmass_lo95 lt 0 then stmass_lo95=0
  if stmass_hi95 lt 0 then stmass_hi95=0

  ;medianIndex = mst[stmass_median]
  ;lo68Index=mst[stmass_lo68]
  ;hi68Index=mst[stmass_hi68]
  ;lo95Index=mst[stmass_lo95]
  ;hi95Index=mst[stmass_hi95]

  bestIndex = (where(y_mass eq max(y_mass)))[0]
  mass_best = mass[bestIndex]

  medianIndex = (where( abs(mass_median-mass) eq min(abs(mass_median-mass))))[0]
  lo68Index = (where( abs(mass_lo68-mass) eq min(abs(mass_lo68-mass))))[0]
  hi68Index = (where( abs(mass_hi68-mass) eq min(abs(mass_hi68-mass))))[0]
  lo95Index = (where( abs(mass_lo95-mass) eq min(abs(mass_lo95-mass))))[0]
  hi95Index = (where( abs(mass_hi95-mass) eq min(abs(mass_hi95-mass))))[0]


  ;message,/cont,'lo68Index= '+strn(lo68Index)
  ;message, /cont,'mass[lo68Index]= '+strn(mass[lo68Index])
  ;message, /cont,'mass_lo68= '+ strn(mass_lo68)
  return, y_mass

END


;------------------------------------------------------------

pro getStats, rawx, y_x, log=log, $
              bestx=bestx,$
              bestIndex=bestIndex, $
              meanx=meanx, $
              meanIndex=meanIndex,$
              sigmax=sigmax, $
              medianx=medianx, $
              medianIndex=medianIndex,$
              lo68x=lo68x, hi68x=hi68x, $
              lo95x=lo95x, hi95x=hi95x, $
              lo68Index=lo68Index, hi68Index=hi68Index, $
              lo95Index=lo95Index, hi95Index=hi95Index

  if keyword_set(log) then x = 10d^rawx else x=rawx

  bestIndex = (where(y_x eq max(y_x)))[0]
  bestx = x[bestIndex]
  if n_elements(x) gt 1 then begin
     meanx = tsum(x,x*y_x)
     sigmax = sqrt(tsum(x, (x-meanx)^2*y_x))
     cumulative_pdf, x, y_x, medianIndex, lo68Index, hi68Index, $
                     lo95Index, hi95Index, acum, $
                     medianx, lo68x, hi68x, lo95x, hi95x, badobj2=badobj2
     if total(size(badobj2)) ne 0 then $
        if badobj2 eq 1 then badobj=1
  endif else begin
     meanx = x[0]
     sigmax = 0.0
     medianx=0.0
     lo68x=0.0
     hi68x=0.0
     lo95x=0.0
     hi95x=0.0
     medianIndex=0
     lo68Index=0
     hi68Index=0
     lo95Index=0
     hi95Index=0
     
  endelse
  meanIndex = findel(meanx,x)
  
  if medianIndex le 0 then medianIndex = 0
  if lo68Index le 0 then lo68Index=0
  if lo95Index le 0 then lo95Index=0
  if hi68Index le 0 then hi68Index=0
  if hi95Index le 0 then hi95Index=0

  if medianIndex ge n_elements(y_x) then medianIndex=n_elements(y_x)-1
  if hi68Index ge n_elements(y_x) then hi68Index=n_elements(y_x)-1
  if hi95Index ge n_elements(y_x) then hi95Index=n_elements(y_x)-1

  if keyword_set(log) then begin
     bestx = alog10(bestx)
     sigmax = mean( abs(alog10(bestx) - [alog10(bestx+sigmax),alog10(bestx-sigmax)]) )
     medianx = alog10(medianx)
     lo68x = lo68x > 0 ? alog10( lo68x) : min(rawx)
     hi68x = alog10( hi68x)
     lo95x = lo95x > 0 ? alog10( lo95x) : min(rawx)
     hi95x = alog10( hi95x)
  endif

END


;------------------------------------------------------------
;------------------------------------------------------------

; function fitchisqGetVal, arr, val, badobj1=badobj1
; ttt=where(arr gt 0 and arr le 1)
; if ttt[0] eq -1 then ttt=where(arr gt 0 and arr le 1.d +1d-15)
; if (ttt)[0] eq -1 or (ttt)[0] gt n_elements(arr)-1 or long(total(finite(arr))) ne n_elements(arr) then begin
;    print, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SALMON EDIT :: BAD OBJECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
;    badobj1 = 1
;    ttt = 1
;    return, ttt
; endif
; ans = (where( abs(arr[ttt] - val) eq min(abs(arr[ttt]-val))))[0]
; return, ttt[ans]
; END


;------------------------------------------------------------
;------------------------------------------------------------
; need to leave cumulative pdf as an integral to allow for non-uniform
; gridding of parameters (otherwise summing would be fine!) 

pro cumulative_pdf, x, pdf, tmedian, tlo68, thi68, tlo95, thi95, cum, badobj2=badobj2, $
                    median, lo68, hi68, lo95, hi95,reverse=reverse, stop=stop, sort=sort
; sort means to sort first (for ISED, for example)
; reverse means to do it in reverse, useful for EBV

; PDF must be sorted from lowest array value to highest array value for mass
; Renormalize in case: 
if keyword_set(sort) then begin
   stx=sort(x)
   usex=x[stx]
   usepdf = pdf[stx]
endif else begin
   stx=lindgen(n_elements(x))
   usex=x[stx]
   usepdf=pdf
endelse

usepdf /= tsum(usex, usepdf)                                                             
cum = pdf*0.0d                                                                           
; The cumulative integral for stellar mass and SFR
; is huge and time consuming. You can achieve the same
; calculation with better than ~0.05 dex accuracy by
; doing the integration over a ~100x more coarse x.
if n_elements(usex) gt 1e3 then begin
   coarsex = [usex[0:*:100]] ; every 100th model
   if coarsex[-1] ne usex[-1] then  begin    ; include the last model
      coarsex=[coarsex,usex[-1]] 
      ;;cum = dblarr(n_elements(pdf)/100 +1 +1)
   endif ;; else begin
      ;;cum = dblarr(n_elements(pdf)/100 +1)
   ;;endelse
   cum = dblarr(n_elements(coarsex)) ;; CJP added this to force cum to have the same number of elements as coarsex... a bug was showing up. 
   integral_x = coarsex
   ;; CHECK if integral_x and cum have the same size.  Needed for
   ;; intepol below.
endif else $
   integral_x = usex
                                                                                         
if not keyword_set(reverse) then begin                                                   
   for i=0L,n_elements(integral_x)-1 do begin                                                  
      if i eq 0 then t = indgen(2) else t=where(usex le integral_x[i])                         
      cum[i] = tsum(usex[t], usepdf[t])                                                  
      if i eq 0 then cum[i] /= 2.                                                        
   endfor                                                                                
endif else begin                                                                         
   ; I don't think we ever reverse integrate Mass 
   ; or SFR.. so I won't fuss with that here
   for i=n_elements(usex)-2, 0, -1 do begin                                              
      t=where(usex ge usex[i])                                                           
      cum[i] = tsum(usex[t], usepdf[t])                                                  
   endfor                                                                                
   cum = 1.0 - cum                                                                       
endelse                                                                                  
                                                                                         
;if cum[0] ne cum[0] then stop                                                           
median = interpol(integral_x, cum, 0.5d)                                                       
lo68 = interpol(integral_x, cum, 0.15865d)                                                     
hi68 = interpol(integral_x, cum, 0.84135d)                                                     
lo95 = interpol(integral_x, cum, 0.02275d)                                                     
hi95 = interpol(integral_x, cum, 0.97725d)                                                     

tmedian = (where( abs(median-usex) eq min(abs(median-usex))))[0]
tlo68 = (where( abs(lo68-usex) eq min(abs(lo68-usex))))[0]
thi68 = (where( abs(hi68-usex) eq min(abs(hi68-usex))))[0]
tlo95 = (where( abs(lo95-usex) eq min(abs(lo95-usex))))[0]
thi95 = (where( abs(hi95-usex) eq min(abs(hi95-usex))))[0]
if keyword_set(stop) then stop

; tmedian = fitchisqGetVal(cum, 0.5d, badobj1=badobj1)
; if total(size(badobj1)) ne 0 then if badobj1 eq 1 then badobj2=1
; tlo68 = fitchisqGetVal(cum, 0.15865d, badobj1=badobj1)
; if total(size(badobj1)) ne 0 then if badobj1 eq 1 then badobj2=1
; thi68 = fitchisqGetVal(cum, 0.84135d, badobj1=badobj1)
; if total(size(badobj1)) ne 0 then if badobj1 eq 1 then badobj2=1
; tlo95 = fitchisqGetVal(cum, 0.02275d, badobj1=badobj1)
; if total(size(badobj1)) ne 0 then if badobj1 eq 1 then badobj2=1
; thi95 = fitchisqGetVal(cum, 0.97725d, badobj1=badobj1)
; if total(size(badobj1)) ne 0 then if badobj1 eq 1 then badobj2=1
if finite(lo68) eq 0 or lo68 lt min(x) then begin
   lo68 = 0.0d
   tlo68 = 0
endif
if finite(lo95) eq 0 or lo95 lt min(x) then begin
   lo95 = 0.0d
   tlo95 = 0
endif
if finite(hi68) eq 0 or hi68 gt max(x) then begin
   hi68 = max(x)
   tlo68 = n_elements(x)-1
endif
if finite(hi95) eq 0 or hi95 gt max(x) then begin
   hi95 = max(x)
   tlo95 = n_elements(x)
endif
if keyword_set(stop) then stop

END



;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------

function mytsum, x, y

  if size(/dim, x) eq 1 then return, y $
  else return, tsum(x,y)

end

;; Calculate maxmium likelihood of each parameter by marginalizing in
;; all directions.  Then derive parameter likelihood ranges
;; (uncertainties).  Mass and SFR are handled separately. 
;
;; If calculating errors (almost always are) then derive cummulative
;; distributions.

function old_y_p, x, pdf, i_x, a=a, b=b, c=c, norm=norm
;; this function takes the PDF, and integrates over a, b, c, and d to
;; get y(x), the posterior for x (marginalized over the other
;; parameters).
;;
;; a = array for the first variable, 
;; b= array for the second variable
;; c= third
;; i_x is assumed to be dimension of the variable for x.
;; 
;; Example:
;; if x is the 0th dimension, then i_x = 0 and pdf would have the form:
;; pdf[*, n_elements(a), n_elements(b), n_elements(c), n_elements(d)]
;;
;; if x is the 1st dimension, then i_x = 1 and the pdf would have the form:
;; pdf[n_elements(a), *, n_elements(b), n_elements(c), n_elements(d)]
;; 
;; etc... 

  sz = size(pdf,/dim)
  szx = sz[i_x]
  indexes = bindgen(4)
  indexes=indexes[where(i_x ne indexes)]

  y_x = fltarr(szx)
  tmpA = fltarr(szx, $
                n_elements(b),$
                n_elements(c))
  tmpB = fltarr(szx, $
                n_elements(c))
;  tmpC = fltarr(szx, $
;                n_elements(d))

  for i=0,n_elements(y_x)-1 do begin
;     for l=0,n_elements(d)-1 do begin
     for k=0,n_elements(c)-1 do begin
        for j=0,n_elements(b)-1 do begin
           case i_x of  
              0 : tmpA[i,j,k] = mytsum(a,pdf[i,*,j,k])
              1 : tmpA[i,j,k] = mytsum(a,pdf[*,i,j,k])
              2 : tmpA[i,j,k] = mytsum(a,pdf[*,j,i,k])
              3 : tmpA[i,j,k] = mytsum(a,pdf[*,j,k,i])
;              4 : tmpA[i,j,k] = mytsum(a,pdf[*,j,k,i])
              else : begin
                 message, 'i_x must be 0...3, cannot be '+strn(i_x)
                 stop
              end
           endcase
        endfor
        tmpB[i,k] = mytsum(b,tmpA[i,*,k])
     endfor
     y_x[i] = mytsum(c,tmpB[i,*])
  endfor

  norm = mytsum(x,y_x)
  y_x /= norm

  return, y_x

END

function old_y_p6, x, pdf, i_x, a=a, b=b, c=c, d=d, e=e, norm=norm
;; Same as y_p, but allows 5 variables (age, ebv, metal, tau, alpha,
;; beta) this function takes the PDF, and integrates over a, b, c, and d to
;; get y(x), the posterior for x (marginalized over the other
;; parameters).
;;
;; a = array for the first variable, 
;; b= array for the second variable
;; c, d,e = third, fourth fifth variable
;; i_x is assumed to be dimension of the variable for x.
;; 
;; Example:
;; if x is the 0th dimension, then i_x = 0 and pdf would have the form:
;; pdf[*, n_elements(a), n_elements(b), n_elements(c), n_elements(d)]
;;
;; if x is the 1st dimension, then i_x = 1 and the pdf would have the form:
;; pdf[n_elements(a), *, n_elements(b), n_elements(c), n_elements(d)]
;; 
;; etc... 

  sz = size(pdf,/dim)
  szx = sz[i_x]
  indexes = bindgen(6)
  indexes=indexes[where(i_x ne indexes)]

  y_x = fltarr(szx)
  tmpA = fltarr(szx, $
                n_elements(b),$
                n_elements(c),$
                n_elements(d),$
                n_elements(e))
  tmpB = fltarr(szx, $
                n_elements(c),$
                n_elements(d),$
                n_elements(e))
  tmpC = fltarr(szx, $
                n_elements(d),$
                n_elements(e))
  tmpD = fltarr(szx, $
                n_elements(e))

  for i=0,n_elements(y_x)-1 do begin
     for m=0,n_elements(e)-1 do begin
        for l=0,n_elements(d)-1 do begin
           for k=0,n_elements(c)-1 do begin
              for j=0,n_elements(b)-1 do begin
                 
                 case i_x of  
                    0 : tmpA[i,j,k,l,m] = mytsum(a,pdf[i,*,j,k,l,m])
                    1 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,i,j,k,l,m])
                    2 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,j,i,k,l,m])
                    3 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,j,k,i,l,m])
                    4 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,j,k,l,i,m])
                    5 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,j,k,l,m,i])
                    else : begin
                       message, 'i_x must be 0...5, cannot be '+strn(i_x)
                       stop
                    end
                 endcase
              endfor
              tmpB[i,k,l,m] = mytsum(b,tmpA[i,*,k,l,m])
           endfor
           tmpC[i,l,m] = mytsum(c,tmpB[i,*,l,m])
        endfor
        tmpD[i,m] = mytsum(d,tmpC[i,*,m])
     endfor
     y_x[i] = mytsum(e,tmpD[i,*])
  endfor
  
  norm = mytsum(x,y_x)
  y_x /= norm

  return, y_x

END




;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------

; REMOVE ALL 2-MASS STUFF.  IT"S JUST NUTS.
; CHANGE AGE to use AGEARR: 

PRO FITSED_FITCHISQ,  $
   ;; this is the photometry, errors, and redshift for the galaxy to
   ;; be fit:
   phot, dphot, z, $
   ;;
   ;; this is the central wavelengths of the fitted filters
   lambda_filters,$
   ;;
   ;; the following are from generate_lut in the same order as in lut[...]
   lut, zed, log_ageArr, ebvArr, metalArr, tauArr, alphaArr, betaArr, lutMstar, lutSFR, $
   nomass=nomass, $ ;; if set it skips the mass derivation; BS says this makes the code go much faster
   ;;
   ;; these control the code: 
   mfactor=mfactor, $
   calc_errors=calc_errors, $
;   num_iter=num_iter, $
;   monte_carlo=monte_carlo, $
;   run_monte_carlo=run_monte_carlo,$
   verbose=verbose, $
   sfh=sfh, $
   tltuniverse=tltuniverse, $
   badobj=badobj, $
   agelimit=agelimit, $ ; can be used to limit ages
   taulimit=taulimit, $ ; can be used to limit taus
   rest_lowerlimit=rest_lowerlimit, $ ; rest-frame wavelength lower limit to the bands you fit
   chisq=chisq, $
   scale=scale, $
   pdf=pdf, $
   mass=mass, $
   sfr=sfr, result=result 

; ------------------------------------------------------------
  
  if tauArr[0] lt 0.0 then tauArr = abs(tauArr)
  
  if not keyword_set(agelimit) then agelimit=0.0
  if not keyword_set(taulimit) then taulimit=0.0
  delvarx,  badobj
  if keyword_set(num_iter) then num_iter=num_iter else num_iter=100
  if not keyword_set(rest_lowerlimit) then rest_lowerlimit=0.0

;; MFACTOR allows the user to scale the masses as needed.  The
;; default is to assume the photometry is in uJy. For nJy, set
;; mfactor=1.0. For uJy, set mfactor=1000.
  if keyword_set(mfactor) then mfactor=mfactor $
  else mfactor=1000.d

  if strcmp(sfh,'dpl') then begin
     sz = size(lut,/dim)
     n_zed = sz[0]
     n_age = sz[1]
     n_ebv = sz[2]
;  n_delta = sz[3]
     n_metal = sz[3]
     n_tau = sz[4]
     n_alpha = sz[5]
     n_beta = sz[6]
     n_filters = sz[7]
  endif else begin
     sz = size(lut,/dim)
     n_zed = sz[0]
     n_age = sz[1]
     n_ebv = sz[2]
;  n_delta = sz[3]
     n_metal = sz[3]
     n_tau = sz[4]
     n_filters = sz[5]
     n_alpha=1
     n_beta=1
  endelse
  

  if n_tau+n_alpha+n_beta eq 1 then begin
      chisq = dblarr(n_age, n_ebv, n_metal)*0.0 + 1d62
      scale = dblarr(n_age, n_ebv, n_metal)
      pdf = dblarr(n_age, n_ebv, n_metal)
      mass = dblarr(n_age, n_ebv, n_metal)
      sfr =  dblarr(n_age, n_ebv, n_metal)
   endif else begin
      if strcmp(sfh,'dpl') then begin
         chisq = dblarr(n_age, n_ebv, n_metal, n_tau,n_alpha,n_beta)*0.0 + 1d62
         scale = dblarr(n_age, n_ebv, n_metal, n_tau,n_alpha,n_beta)
         pdf = dblarr(n_age, n_ebv, n_metal, n_tau,n_alpha,n_beta)
         mass = dblarr(n_age, n_ebv, n_metal, n_tau,n_alpha,n_beta)
         sfr =  dblarr(n_age, n_ebv, n_metal, n_tau,n_alpha,n_beta)
      endif else begin
         chisq = dblarr(n_age, n_ebv, n_metal, n_tau)*0.0 + 1d62
         scale = dblarr(n_age, n_ebv, n_metal, n_tau)
         pdf = dblarr(n_age, n_ebv,  n_metal, n_tau)
         mass = dblarr(n_age, n_ebv, n_metal, n_tau)
         sfr =  dblarr(n_age, n_ebv, n_metal, n_tau)
      endelse         
  endelse

  ;; find the index in the zed array for the galaxies' redshift
  z_from_catalog=z
  if z lt min(zed) then begin
     ;result=0.
     ;mass=0.
     ;sfr=0.
     if z lt min(zed) then z = min(zed)
     message,$
        ' redshift of galaxy ('+strn(z_from_catalog)+')'+$
        ' is outside range of models, forcing z=z_min, ('+strn(z)+')...',/continue
        ;' is outside range of models, skipping...',/continue
     ;return
  endif else if z gt max(zed) then begin
     if z gt max(zed) then z = max(zed)
     message,$
        ' redshift of galaxy '+strn(z)+$
        ' is outside range of models, forcing z=z_max...',/continue
        ;' is outside range of models, skipping...',/continue
  endif
  ; BWS EDIT here
  if 0 then begin
  endif else begin

     zind = (where( abs(zed-z) eq min(abs(zed-z))))[0]

     ;; calculate lookback time for object 
     lbt = fitsed_lookback_time(z,1100.d,h=!h,omega=!omega,lambda=!lambda)
     if not keyword_set(tltuniverse) then lbt = 1d100 ; really big number

     myvec = dblarr(n_filters)

;     if keyword_set(calc_errors) and keyword_set(run_monte_carlo) then begin
;        Message,'% WARNING. fitchisq_metal is not set up for Monte Carlo, due to the change in variables', continue
;        stop
;     endif

     ;print, "%%%% start chisq"
     ;timestart00=systime(/seconds)
     for b=0,n_beta-1 do begin
        for a=0,n_alpha-1 do begin
           for k=0,n_tau-1 do begin
              for l=0,n_metal-1 do begin
;                 for d=0l, n_delta-1 do begin
                 for j=0,n_ebv-1 do begin
                    for i=0,n_age-1 do begin
                       if strcmp(sfh,'dpl') then begin
                          myvec[*] = lut[zind,i,j,l,k,a,b,*]
                       endif else begin
                          myvec[*] = lut[zind,i,j,l,k,*]
                       endelse

                       sel_fin=where( finite(phot) eq 1 and finite(dphot) eq 1 and $
                                      lambda_filters gt rest_lowerlimit and $ ; BWS addition - ignore bands at low lam_rest
                                      finite(myvec) eq 1 and $ ; ignore bad LUT filters, just in case
                                      phot gt -98 and dphot gt -98)
                    
                       if n_elements(sel_fin) ge 2 then begin
                          tscale = $
                             svdfit( myvec[sel_fin], phot[sel_fin], 1, $
                                     a=ta, chisq=xtchisq, $
                                     measure_errors=dphot[sel_fin], $
                                     sigma=tsigma, $
                                     yfit=tyfit,/double, funct='fitsed_fitfunc')
                          tscale=tscale[0]
                          ;;----------------------------------------
                       ;; if the age is greater than the age of the universe, 
                          ;; then max out the chisq (zero out the pdf)
                          if 10d^log_ageArr[i] le lbt and $
                             10d^log_ageArr[i] ge agelimit then $
                                tchisq =total( (myvec[sel_fin]*tscale - $
                                                phot[sel_fin])^2 / $
                                               dphot[sel_fin]^2 ) 
                          
;                              chisq[i,j,d,l,k] = $
;                              total( (myvec[sel_fin]*scale[i,j,d,l,k] - $
;                                      phot[sel_fin])^2 / $
;                                     dphot[sel_fin]^2 ) 
                       endif  else begin ;; there are <= 2 good photometry points:
                          tscale=0.0d    ; scale[i,j,d,l,k] = 0.0000d
                          ;; leave chisq, etc, alone
                       endelse
                    
                    ;;----------------------------------------
                    ;; calculate the stellar mass and star formation rate

                       if strcmp(sfh,'dpl') then begin
                          chisq[i,j,l,k,a,b] = tchisq
                          scale[i,j,l,k,a,b] = tscale
                          mass[i,j,l,k,a,b] = scale[i,j,l,k,a,b]*mfactor
                          sfr[i,j,l,k,a,b] = lutSFR[i,l,k,a,b] / $
                                               lutMstar[i,l,k,a,b] * mass[i,j,l,k,a,b]
                       endif else begin
                          chisq[i,j,l,k] = tchisq
                          scale[i,j,l,k] = tscale
                          mass[i,j,l,k] = scale[i,j,l,k]*mfactor
                          sfr[i,j,l,k] = lutSFR[i,l,k] / $
                                           lutMstar[i,l,k] * mass[i,j,l,k]
                       endelse
                    
                    endfor  ;; age
                 endfor;; ebv
              endfor   ;; metallicity
           endfor      ;; tau
        endfor         ;; alpha
     endfor            ;; beta
   ;print, "%% All chisq "
   ;fitsed_checktime,timestart00

;; ------------------------------------------------------------
;; Define likelihood as PDF = exp(-chisq/2) 
;;
;; TODO:  allow for additional priors here.

     minchisq  = min(chisq)
     pdf = exp( - (chisq - min(chisq)) / 2.0 )
     ;minchisq  = min(chisq)
     ;pdf = exp( - (chisq) / 2.0 )
;; ===============================
;; Integrate over everything:

     if strcmp(sfh,'dpl') then begin
        y_age = fitsed_yp6(10d^log_ageArr, pdf, 0, a=ebvArr, b=metalArr, c=tauArr, d=alphaArr, e=betaArr,$
                     norm=norm)
;apply the normalization:
        pdf /= norm
        
        y_ebv = fitsed_yp6(ebvArr, pdf, 1, a=10d^log_ageArr, b=metalArr, c=tauArr, d=alphaArr, e=betaArr)
        y_metal = fitsed_yp6(metalArr, pdf, 2, a=10d^log_ageArr, b=ebvArr, c=tauArr, d=alphaArr, e=betaArr)
        if n_tau gt 1 then $
           y_tau = fitsed_yp6(tauArr, pdf, 3, a=10d^log_ageArr, b=ebvArr, c=metalArr, d=AlphaArr, e=betaArr) $
        else   y_tau = [1.0]
        if n_alpha gt 1 then $
           y_alpha = fitsed_yp6(alphaArr, pdf, 4, a=10d^log_ageArr, b=ebvArr, c=metalArr, d=tauArr, e=betaArr) $
        else y_alpha = [1.0]
        if n_beta gt 1 then $
           y_beta = fitsed_yp6(betaArr, pdf, 5, a=10d^log_ageArr, b=ebvArr, c=metalArr, d=tauArr, e=alphaArr) $
        else y_beta = [1.0]
        
     endif else begin
        y_age = fitsed_yp(10d^log_ageArr, pdf, 0, a=ebvArr, b=metalArr, c=tauArr, $
                    norm=norm)
;apply the normalization:
        pdf /= norm

        y_ebv = fitsed_yp(ebvArr, pdf, 1, a=10d^log_ageArr, b=metalArr, c=tauArr)
        y_metal = fitsed_yp(metalArr, pdf, 2, a=10d^log_ageArr, b=ebvArr, c=tauArr)
        if n_tau gt 1 then begin
           y_tau = fitsed_yp(tauArr, pdf, 3, a=10d^log_ageArr, b=ebvArr, c=metalArr)
        endif else begin
           y_tau = [1.0]
        endelse
     endelse

     ;; log Age
     getStats, log_ageArr, y_age, /log, $
               bestx = log_age_best, $
               bestIndex=tlog_age_best, $
               meanx = log_age_mean, $
               meanIndex=tlog_age_mean, $
               sigmax = log_age_sigma, $
               medianx = log_age_median, $
               medianIndex= tlog_age_median, $
               lo68x = log_age_lo68, $
               hi68x = log_age_hi68, $
               lo95x = log_age_lo95, $
               hi95x = log_age_hi95, $
               lo68Index = tlog_age_lo68, $
               hi68Index = tlog_age_hi68, $
               lo95Index = tlog_age_lo95, $
               hi95Index  = tlog_age_hi95

     ;; E(B-V)
     getStats, ebvArr, y_ebv,  $
               bestx = ebv_best, $
               bestIndex=tebv_best, $
               meanx = ebv_mean, $
               meanIndex=tebv_mean, $
               sigmax = ebv_sigma, $
               medianx = ebv_median, $
               medianIndex= tebv_median, $
               lo68x = ebv_lo68, $
               hi68x = ebv_hi68, $
               lo95x = ebv_lo95, $
               hi95x = ebv_hi95, $
               lo68Index = tebv_lo68, $
               hi68Index = tebv_hi68, $
               lo95Index = tebv_lo95, $
               hi95Index  = tebv_hi95

     ;; Metallicity
     getStats, metalArr, y_metal, $
               bestx = metal_best, $
               bestIndex=tmetal_best, $
               meanx = metal_mean, $
               meanIndex=tmetal_mean, $
               sigmax = metal_sigma, $
               medianx = metal_median, $
               medianIndex= tmetal_median, $
               lo68x = metal_lo68, $
               hi68x = metal_hi68, $
               lo95x = metal_lo95, $
               hi95x = metal_hi95, $
               lo68Index = tmetal_lo68, $
               hi68Index = tmetal_hi68, $
               lo95Index = tmetal_lo95, $
               hi95Index  = tmetal_hi95

     ;; Tau
     if n_tau gt 1 then begin
     getStats, tauArr, y_tau, $
               bestx = tau_best, $
               bestIndex=ttau_best, $
               meanx = tau_mean, $
               meanIndex=ttau_mean, $
               sigmax = tau_sigma, $
               medianx = tau_median, $
               medianIndex= ttau_median, $
               lo68x = tau_lo68, $
               hi68x = tau_hi68, $
               lo95x = tau_lo95, $
               hi95x = tau_hi95, $
               lo68Index = ttau_lo68, $
               hi68Index = ttau_hi68, $
               lo95Index = ttau_lo95, $
               hi95Index  = ttau_hi95
     endif else begin
               tau_best     = 0.0   
               ttau_best    = 0.0   
               tau_mean     = 0.0   
               ttau_mean    = 0.0   
               tau_sigma    = 0.0   
               tau_median   = 0.0   
               ttau_median  = 0.0   
               tau_lo68     = 0.0   
               tau_hi68     = 0.0   
               tau_lo95     = 0.0   
               tau_hi95     = 0.0   
               ttau_lo68    = 0.0   
               ttau_hi68    = 0.0   
               ttau_lo95    = 0.0   
               ttau_hi95    = 0.0   
     endelse                        

          ;; Alpha
     if n_alpha gt 1 then begin
     getStats, alphaArr, y_alpha, $
               bestx = alpha_best, $
               bestIndex=talpha_best, $
               meanx = alpha_mean, $
               meanIndex=talpha_mean, $
               sigmax = alpha_sigma, $
               medianx = alpha_median, $
               medianIndex= talpha_median, $
               lo68x = alpha_lo68, $
               hi68x = alpha_hi68, $
               lo95x = alpha_lo95, $
               hi95x = alpha_hi95, $
               lo68Index = talpha_lo68, $
               hi68Index = talpha_hi68, $
               lo95Index = talpha_lo95, $
               hi95Index  = talpha_hi95
     endif else begin
               alpha_best     = 0.0   
               talpha_best    = 0.0   
               alpha_mean     = 0.0   
               talpha_mean    = 0.0   
               alpha_sigma    = 0.0   
               alpha_median   = 0.0   
               talpha_median  = 0.0   
               alpha_lo68     = 0.0   
               alpha_hi68     = 0.0   
               alpha_lo95     = 0.0   
               alpha_hi95     = 0.0   
               talpha_lo68    = 0.0   
               talpha_hi68    = 0.0   
               talpha_lo95    = 0.0   
               talpha_hi95    = 0.0   
     endelse                        

               ;; Beta
     if n_beta gt 1 then begin
     getStats, betaArr, y_beta, $
               bestx = beta_best, $
               bestIndex=tbeta_best, $
               meanx = beta_mean, $
               meanIndex=tbeta_mean, $
               sigmax = beta_sigma, $
               medianx = beta_median, $
               medianIndex= tbeta_median, $
               lo68x = beta_lo68, $
               hi68x = beta_hi68, $
               lo95x = beta_lo95, $
               hi95x = beta_hi95, $
               lo68Index = tbeta_lo68, $
               hi68Index = tbeta_hi68, $
               lo95Index = tbeta_lo95, $
               hi95Index  = tbeta_hi95
     endif else begin
               beta_best     = 0.0   
               tbeta_best    = 0.0   
               beta_mean     = 0.0   
               tbeta_mean    = 0.0   
               beta_sigma    = 0.0   
               beta_median   = 0.0   
               tbeta_median  = 0.0   
               beta_lo68     = 0.0   
               beta_hi68     = 0.0   
               beta_lo95     = 0.0   
               beta_hi95     = 0.0   
               tbeta_lo68    = 0.0   
               tbeta_hi68    = 0.0   
               tbeta_lo95    = 0.0   
               tbeta_hi95    = 0.0   
     endelse                        

     
; Calculate stellar mass and SFR; Treated differently.  
; Need to sort-list masses over all of PDF, and
; start to add-up total PDF until confidence regions found.
     if not keyword_set(nomass) then begin
        y_mass = fitsed_calc_pmass(mass, pdf, mst=mst, $
                           mass_best = mass_best, bestIndex=tmass_best, $
                           mass_mean = mass_mean, meanIndex=tmass_mean, $
                           mass_sigma=mass_sigma, $
                           mass_median=mass_median, medianIndex=tmass_median, $
                           mass_lo68=mass_lo68, mass_hi68=mass_hi68, $
                           mass_lo95=mass_lo95, mass_hi95=mass_hi95, $
                           lo68Index=tmass_lo68, hi68Index=tmass_hi68, $
                           lo95Index=tmass_lo95, hi95Index=tmass_hi95)
        
        y_sfr = fitsed_calc_pmass(sfr, pdf, mst=sft, $
                          mass_best = sfr_best, bestIndex=tsfr_best, $
                          mass_mean = sfr_mean, meanIndex=tsfr_mean, $
                          mass_sigma=sfr_sigma, $
                          mass_median=sfr_median, medianIndex=tsfr_median, $
                          mass_lo68=sfr_lo68, mass_hi68=sfr_hi68, $
                          mass_lo95=sfr_lo95, mass_hi95=sfr_hi95, $
                          lo68Index=tsfr_lo68, hi68Index=tsfr_hi68, $
                          lo95Index=tsfr_lo95, hi95Index=tsfr_hi95)

        ;; TO-DO:  SFR SFR SFR -- interupt this step to average over
        ;;                        sfh with some timescale (probably
        ;;                        100 Myr) and report that!!
        
     endif ;; nomass

;;------------------------------------------------------------
;; calculate minimum value of chisq directly and get indexes:
     x = (where(chisq eq min(chisq)))[0]
;; Define indices of the best-fit model for each index:

     if strcmp(sfh,'dpl') then begin
        ibeta= x / n_metal / n_tau / n_alpha / n_ebv / n_age
        ialpha= (x-ibeta*n_age*n_ebv*n_metal*n_tau*n_alpha) / n_tau / n_metal / n_ebv / n_age
        itau = (x - ibeta*n_age*n_ebv*n_metal*n_tau*n_alpha $
                - ialpha*n_age*n_ebv*n_tau*n_metal) $
               / n_ebv / n_age / n_metal
        imetal = (x -ibeta*n_age*n_ebv*n_metal*n_tau*n_alpha $
                  - ialpha*n_age*n_ebv*n_tau*n_metal $
                  -itau*n_age*n_ebv*n_metal) / n_ebv / n_age
        iebv = (x - ibeta*n_age*n_ebv*n_metal*n_tau*n_alpha $
                - ialpha*n_age*n_ebv*n_tau*n_metal $
                -itau*n_age*n_ebv*n_metal $
                -imetal*n_age*n_ebv) /  n_age
        iage = x - ibeta*n_age*n_ebv*n_metal*n_tau*n_alpha $
               - ialpha*n_age*n_ebv*n_tau*n_metal $
               - itau*n_age*n_ebv*n_metal $
               - imetal*n_age*n_ebv $
               - iebv*n_age
     endif else begin
        itau = x / n_metal / n_ebv / n_age
        imetal = (x - itau*n_age*n_ebv*n_metal) / n_ebv / n_age
        iebv = (x - itau*n_age*n_ebv*n_metal $
                - imetal*n_age*n_ebv) / n_age
        iage = x - itau*n_age*n_ebv*n_metal $
               - imetal*n_age*n_ebv $
               - iebv*n_age
     endelse
     
     
;; make result structure:
     delvarx, result

     if total(size(badobj)) gt 0 then begin
        if badobj eq 1 then result=0 
     endif else begin
        badobj=0

        if keyword_set(nomass) then begin
           mass_struct = 0B
           sfr_struct = 0B
        endif else begin
           mass_struct = { mean: mass_mean, sigma: mass_sigma, $
                           mean_index: tmass_mean, $
                           minchisq_ind: x, minchisq: mass[x], $
                           maxpdf_ind: tmass_best, $
                           maxpdf: mass_best, $
                           ;mass_best: mass_best, $
                           median: mass_median,$ 
                           lo68: mass_lo68, hi68: mass_hi68, $
                           lo95: mass_lo95, $
                           hi95: mass_hi95, $
                           median_ind: tmass_median, $
                           lo68_ind: tmass_lo68, $
                           hi68_ind: tmass_hi68, $
                           lo95_ind: tmass_lo95, $
                           hi95_ind: tmass_hi95 }
           
           sfr_struct = { mean: sfr_mean, sigma: sfr_sigma, $
                          mean_ind: tsfr_mean, $
                          minchisq_ind: x, minchisq: sfr[x], $
                          maxpdf_ind: tsfr_best, $
                          maxpdf: sfr[tsfr_best], $
                          median: sfr_median,$ 
                          lo68: sfr_lo68, hi68: sfr_hi68, $
                          lo95: sfr_lo95, $
                          hi95: sfr_hi95, $
                          median_ind: tsfr_median, $
                          lo68_ind: tsfr_lo68, $
                          hi68_ind: tsfr_hi68, $
                          lo95_ind: tsfr_lo95, $
                          hi95_ind: tsfr_hi95 }
        endelse

        if strcmp(sfh,'dpl') then $
           result = {y_age: y_age, $
                     y_ebv:y_ebv, $
                     y_metal: y_metal, $
                     y_tau: y_tau, $
                     y_alpha:y_alpha, $
                     y_beta:y_beta, $
                     y_mass: y_mass, $
                     y_sfr: y_sfr, $
                     minchisq: chisq[x], $
                     norm: norm,$
                     log_age: {mean: log_age_mean, sigma: log_age_sigma, $
                               mean_ind: tlog_age_mean,$
                               minchisq_ind: iage, minchisq: log_ageArr[iage], $
                               maxpdf_ind: tlog_age_median, maxpdf: log_ageArr[tlog_age_best], $
                               median: log_age_median, median_ind: tlog_age_median, $
                               lo68_ind: tlog_age_lo68, hi68_ind: tlog_age_hi68, $
                               lo95_ind: tlog_age_lo95, hi95_ind: tlog_age_hi95, $
                               lo68: log_age_lo68, hi68: log_age_hi68, $
                               lo95: log_age_lo95, hi95: log_age_hi95 }, $
                     ebv: {mean: ebv_mean, sigma: ebv_sigma, $
                           mean_ind: tebv_mean,$
                           minchisq_ind: iebv, minchisq: ebvArr[iebv], $
                           maxpdf_ind: tebv_median, maxpdf: ebvArr[tebv_best], $
                           median: ebv_median, median_ind: tebv_median, $
                           lo68_ind: tebv_lo68, hi68_ind: tebv_hi68, $
                           lo95_ind: tebv_lo95, hi95_ind: tebv_hi95, $
                           lo68: ebv_lo68, hi68: ebv_hi68, $
                           lo95: ebv_lo95, hi95: ebv_hi95 }, $
                     metal: {mean: metal_mean, sigma: metal_sigma, $
                             mean_ind: tmetal_mean,$
                             minchisq_ind: imetal, minchisq: metalArr[imetal], $
                             maxpdf_ind: tmetal_median, maxpdf: metalArr[tmetal_best], $
                             median: metal_median, median_ind: tmetal_median, $
                             lo68_ind: tmetal_lo68, hi68_ind: tmetal_hi68, $
                             lo95_ind: tmetal_lo95, hi95_ind: tmetal_hi95, $
                             lo68: metal_lo68, hi68: metal_hi68, $
                             lo95: metal_lo95, hi95: metal_hi95 }, $
                     tau: {mean: tau_mean, sigma: tau_sigma, $
                           mean_ind: ttau_mean,$
                           minchisq_ind: itau, minchisq: tauArr[itau], $
                           maxpdf_ind: ttau_median, maxpdf: tauArr[ttau_best], $
                           median: tau_median, median_ind: ttau_median, $
                           lo68_ind: ttau_lo68, hi68_ind: ttau_hi68, $
                           lo95_ind: ttau_lo95, hi95_ind: ttau_hi95, $
                           lo68: tau_lo68, hi68: tau_hi68, $
                           lo95: tau_lo95, hi95: tau_hi95 }, $
                     alpha: {mean: alpha_mean, sigma: alpha_sigma, $
                             mean_ind: talpha_mean,$
                             minchisq_ind: ialpha, minchisq: alphaArr[ialpha], $
                             maxpdf_ind: talpha_median, maxpdf: alphaArr[talpha_best], $
                             median: alpha_median, median_ind: talpha_median, $
                             lo68_ind: talpha_lo68, hi68_ind: talpha_hi68, $
                             lo95_ind: talpha_lo95, hi95_ind: talpha_hi95, $
                             lo68: alpha_lo68, hi68: alpha_hi68, $
                             lo95: alpha_lo95, hi95: alpha_hi95 }, $
                     beta: {mean: beta_mean, sigma: beta_sigma, $
                             mean_ind: tbeta_mean,$
                             minchisq_ind: ibeta, minchisq: betaArr[ibeta], $
                             maxpdf_ind: tbeta_median, maxpdf: betaArr[tbeta_best], $
                             median: beta_median, median_ind: tbeta_median, $
                             lo68_ind: tbeta_lo68, hi68_ind: tbeta_hi68, $
                             lo95_ind: tbeta_lo95, hi95_ind: tbeta_hi95, $
                             lo68: beta_lo68, hi68: beta_hi68, $
                             lo95: beta_lo95, hi95: beta_hi95 }, $
                     mass: mass_struct, sfr: sfr_struct $
                    } $
        else $
           result = {y_age: y_age, $
                     y_ebv:y_ebv, $
                     y_metal: y_metal, $
                     y_tau: y_tau, $
                     y_mass: y_mass, $
                     y_sfr: y_sfr, $
                     minchisq: chisq[x], $
                     norm: norm,$
                     log_age: {mean: log_age_mean, sigma: log_age_sigma, $
                               mean_ind: tlog_age_mean,$
                               minchisq_ind: iage, minchisq: log_ageArr[iage], $
                               maxpdf_ind: tlog_age_median, maxpdf: log_ageArr[tlog_age_best], $
                               median: log_age_median, median_ind: tlog_age_median, $
                               lo68_ind: tlog_age_lo68, hi68_ind: tlog_age_hi68, $
                               lo95_ind: tlog_age_lo95, hi95_ind: tlog_age_hi95, $
                               lo68: log_age_lo68, hi68: log_age_hi68, $
                               lo95: log_age_lo95, hi95: log_age_hi95 }, $
                     ebv: {mean: ebv_mean, sigma: ebv_sigma, $
                           mean_ind: tebv_mean,$
                           minchisq_ind: iebv, minchisq: ebvArr[iebv], $
                           maxpdf_ind: tebv_median, maxpdf: ebvArr[tebv_best], $
                           median: ebv_median, median_ind: tebv_median, $
                           lo68_ind: tebv_lo68, hi68_ind: tebv_hi68, $
                           lo95_ind: tebv_lo95, hi95_ind: tebv_hi95, $
                           lo68: ebv_lo68, hi68: ebv_hi68, $
                           lo95: ebv_lo95, hi95: ebv_hi95 }, $
                     metal: {mean: metal_mean, sigma: metal_sigma, $
                             mean_ind: tmetal_mean,$
                             minchisq_ind: imetal, minchisq: metalArr[imetal], $
                             maxpdf_ind: tmetal_median, maxpdf: metalArr[tmetal_best], $
                             median: metal_median, median_ind: tmetal_median, $
                             lo68_ind: tmetal_lo68, hi68_ind: tmetal_hi68, $
                             lo95_ind: tmetal_lo95, hi95_ind: tmetal_hi95, $
                             lo68: metal_lo68, hi68: metal_hi68, $
                             lo95: metal_lo95, hi95: metal_hi95 }, $
                     tau: {mean: tau_mean, sigma: tau_sigma, $
                           mean_ind: ttau_mean,$
                           minchisq_ind: itau, minchisq: tauArr[itau], $
                           maxpdf_ind: ttau_median, maxpdf: tauArr[ttau_best], $
                           median: tau_median, median_ind: ttau_median, $
                           lo68_ind: ttau_lo68, hi68_ind: ttau_hi68, $
                           lo95_ind: ttau_lo95, hi95_ind: ttau_hi95, $
                           lo68: tau_lo68, hi68: tau_hi68, $
                           lo95: tau_lo95, hi95: tau_hi95 }, $
                     mass: mass_struct, sfr: sfr_struct $
                    } 
           
     endelse ;; badobj=0
     
     
  endelse ;; back to BWS edit... don't know what that is

  ;stop

end



