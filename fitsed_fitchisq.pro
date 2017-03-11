
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

;------------------------------------------------------------
;------------------------------------------------------------
; 
; functions called by fitchisq: 

function calc_mass, mass, pdf, mst=mst, $
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

if not keyword_set(reverse) then begin
   for i=0L,n_elements(usex)-1 do begin
      if i eq 0 then t = indgen(2) else t=where(usex le usex[i])
      cum[i] = tsum(usex[t], usepdf[t])
      if i eq 0 then cum[i] /= 2.
   endfor
endif else begin
   for i=n_elements(usex)-2, 0, -1 do begin
      t=where(usex ge usex[i])
      cum[i] = tsum(usex[t], usepdf[t])
   endfor
   cum = 1.0 - cum
endelse

;if cum[0] ne cum[0] then stop 
median = interpol(usex, cum, 0.5d)
lo68 = interpol(usex, cum, 0.15865d)
hi68 = interpol(usex, cum, 0.84135d)
lo95 = interpol(usex, cum, 0.02275d)
hi95 = interpol(usex, cum, 0.97725d)

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

function y_p, x, pdf, i_x, a=a, b=b, c=c, d=d, norm=norm
;; this function takes the PDF, and integrates over a, b, c, and d to
;; get y(x), the posterior for x (marginalized over the other
;; parameters).
;;
;; a = array for the first variable, 
;; b= array for the second variable
;; c, d= third, fourth variable
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
  indexes = bindgen(5)
  indexes=indexes[where(i_x ne indexes)]

  y_x = fltarr(szx)
  tmpA = fltarr(szx, $
                n_elements(b),$
                n_elements(c),$
                n_elements(d))
  tmpB = fltarr(szx, $
                n_elements(c),$
                n_elements(d))
  tmpC = fltarr(szx, $
                n_elements(d))

  for i=0,n_elements(y_x)-1 do begin
     for l=0,n_elements(d)-1 do begin
        for k=0,n_elements(c)-1 do begin
           for j=0,n_elements(b)-1 do begin
              case i_x of  
                 0 : tmpA[i,j,k,l] = mytsum(a,pdf[i,*,j,k,l])
                 1 : tmpA[i,j,k,l] = mytsum(a,pdf[*,i,j,k,l])
                 2 : tmpA[i,j,k,l] = mytsum(a,pdf[*,j,i,k,l])
                 3 : tmpA[i,j,k,l] = mytsum(a,pdf[*,j,k,i,l])
                 4 : tmpA[i,j,k,l] = mytsum(a,pdf[*,j,k,l,i])
                 else : begin
                    message, 'i_x must be 0...4, cannot be '+strn(i_x)
                    stop
                 end
              endcase
           endfor
           tmpB[i,k,l] = mytsum(b,tmpA[i,*,k,l])
        endfor
        tmpC[i,l] = mytsum(c,tmpB[i,*,l])
     endfor
     y_x[i] = mytsum(d,tmpC[i,*])
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
   ;; the following are from generate_lut in the same order as in lut[...]
   lut, zed, log_ageArr, ebvArr, deltaArr, metalArr, tauArr, lutMstar, lutSFR, $
   nomass=nomass, $ ;; if set it skips the mass derivation; BS says this makes the code go must faster
   ;;
   ;; these control the code: 
   mfactor=mfactor, $
   calc_errors=calc_errors, $
;   num_iter=num_iter, $
;   monte_carlo=monte_carlo, $
;   run_monte_carlo=run_monte_carlo,$
   verbose=verbose, $
   tltuniverse=tltuniverse, $
   badobj=badobj, $
   agelimit=agelimit, $ ; can be used to limit ages
   taulimit=taulimit, $ ; can be used to limit taus
   chisq=chisq, $
   scale=scale, $
   pdf=pdf, $
   mass=mass, $
   sfr=sfr, result=result

  
  if not keyword_set(agelimit) then agelimit=0.0
  if not keyword_set(taulimit) then taulimit=0.0
  delvarx,  badobj
  if keyword_set(num_iter) then num_iter=num_iter else num_iter=100

;; MFACTOR allows the user to scale the masses as needed.  The
;; default is to assume the photometry is in uJy. For nJy, set
;; mfactor=1.0. For uJy, set mfactor=1000.
  if keyword_set(mfactor) then mfactor=mfactor $
  else mfactor=1000.d

  sz = size(lut,/dim)
  n_zed = sz[0]
  n_age = sz[1]
  n_ebv = sz[2]
  n_delta = sz[3]
  n_metal = sz[4]
  n_tau = sz[5]
  n_filters = sz[6]

  chisq = dblarr(n_age, n_ebv, n_delta, n_metal, n_tau)*0.0 + 1d62
  scale = dblarr(n_age, n_ebv, n_delta, n_metal, n_tau)
  pdf = dblarr(n_age, n_ebv, n_delta, n_metal, n_tau)
  mass = dblarr(n_age, n_ebv, n_delta, n_metal, n_tau)
  sfr =  dblarr(n_age, n_ebv, n_delta, n_metal, n_tau)

  ;; find the index in the zed array for the galaxies' redshift
  if z lt min(zed) or z gt max(zed) then begin
     result=0.
     mass=0.
     sfr=0.
     message,$
        ' redshift of galaxy '+strn(z)+$
        ' is outside range of models, skipping...',/continue
     return
  endif else begin

     zind = (where( abs(zed-z) eq min(abs(zed-z))))[0]
     ;; calculate lookback time for object 
     lbt = lookback_time(z,1100.d,h=!h,omega=!omega,lambda=!lambda)
     if not keyword_set(tltuniverse) then lbt = 1d100 ; really big number

     myvec = dblarr(n_filters)

;     if keyword_set(calc_errors) and keyword_set(run_monte_carlo) then begin
;        Message,'% WARNING. fitchisq_metal is not set up for Monte Carlo, due to the change in variables', continue
;        stop
;     endif

     for k=0,n_tau-1 do begin
        for l=0,n_metal-1 do begin
           for d=0l, n_delta-1 do begin
              for j=0,n_ebv-1 do begin
                 for i=0,n_age-1 do begin
                    myvec[*] = lut[zind,i,j,d,l,k,*]
                    sel_fin=where( finite(phot) eq 1 and finite(dphot) eq 1 and $
                                   phot gt -98 and dphot gt -98)

                    if n_elements(sel_fin) ge 2 then begin
                       scale[i,j,d,l,k] = $
                          svdfit( myvec[sel_fin], phot[sel_fin], 1, $
                                  a=ta, chisq=xtchisq, $
                                  measure_errors=dphot[sel_fin], $
                                  sigma=tsigma, $
                                  yfit=tyfit,/double, funct='fitsed_fitfunc')

                       ;;----------------------------------------
                       ;; if the age is greater than the age of the universe, 
                       ;; then max out the chisq (zero out the pdf)
                       if 10d^log_ageArr[i] le lbt and $
                          10d^log_ageArr[i] ge agelimit then $
                          chisq[i,j,d,l,k] = $
                          total( (myvec[sel_fin]*scale[i,j,d,l,k] - $
                                  phot[sel_fin])^2 / $
                                 dphot[sel_fin]^2 )
                    endif  else begin ;; there are <= 2 good photometry points:
                       scale[i,j,d,l,k] = 0.0000d
                       ;; leave chisq, etc, alone
                    endelse
               
                    ;;----------------------------------------
                    ;; calculate the stellar mass and star formation rate

                    mass[i,j,d,l,k] = scale[i,j,d,l,k]*mfactor
                    sfr[i,j,d,l,k] = lutSFR[i,l,k] / $
                                     lutMstar[i,l,k] * mass[i,j,d,l,k]
               
                 endfor  ;; age
              endfor ;; ebv
           endfor ;; delta
        endfor ;; metallicity
     endfor ;; tau

;; ------------------------------------------------------------
;; Define likelihood as PDF = exp(-chisq/2) 
;;
;; TODO:  allow for additional priors here.

     minchisq  = min(chisq)
     pdf = exp( - (chisq ) / 2.0 )
;; ===============================
;; Integrate over everything: 
     y_age = y_p(10d^log_ageArr, pdf, 0, a=ebvArr, b=deltaArr, c=metalArr, d=tauArr, $
                 norm=norm)
;apply the normalization:
     pdf /= norm

     y_ebv = y_p(ebvArr, pdf, 1, a=10d^log_ageArr, b=deltaArr, c=metalArr, d=tauArr)
     y_delta = y_p(deltaArr, pdf, 2, a=10d^log_ageArr, b=ebvArr, c=metalArr, d=tauArr)
     y_metal = y_p(metalArr, pdf, 3, a=10d^log_ageArr, b=ebvArr, c=deltaArr, d=tauArr)
     y_tau = y_p(tauArr, pdf, 4, a=10d^log_ageArr, b=ebvArr, c=deltaArr, d=metalArr)
     
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

     ;; Delta
     getStats, deltaArr, y_delta, $
               bestx = delta_best, $
               bestIndex=tdelta_best, $
               meanx = delta_mean, $
               meanIndex=tdelta_mean, $
               sigmax = delta_sigma, $
               medianx = delta_median, $
               medianIndex= tdelta_median, $
               lo68x = delta_lo68, $
               hi68x = delta_hi68, $
               lo95x = delta_lo95, $
               hi95x = delta_hi95, $
               lo68Index = tdelta_lo68, $
               hi68Index = tdelta_hi68, $
               lo95Index = tdelta_lo95, $
               hi95Index  = tdelta_hi95

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

     


; Calculate stellar mass and SFR; Treated differently.  
; Need to sort-list masses over all of PDF, and
; start to add-up total PDF until confidence regions found.

     if not keyword_set(nomass) then begin
   
        y_mass = calc_mass(mass, pdf, mst=mst, $
                           mass_best = mass_best, bestIndex=tmass_best, $
                           mass_mean = mass_mean, meanIndex=tmass_mean, $
                           mass_sigma=mass_sigma, $
                           mass_median=mass_median, medianIndex=tmass_median, $
                           mass_lo68=mass_lo68, mass_hi68=mass_hi68, $
                           mass_lo95=mass_lo95, mass_hi95=mass_hi95, $
                           lo68Index=tmass_lo68, hi68Index=tmass_hi68, $
                           lo95Index=tmass_lo95, hi95Index=tmass_hi95)
        
        y_sfr = calc_mass(sfr, pdf, mst=sft, $
                          mass_best = sfr_best, bestIndex=tsfr_best, $
                          mass_mean = sfr_mean, meanIndex=tsfr_mean, $
                          mass_sigma=sfr_sigma, $
                          mass_median=sfr_median, medianIndex=tsfr_median, $
                          mass_lo68=sfr_lo68, mass_hi68=sfr_hi68, $
                          mass_lo95=sfr_lo95, mass_hi95=sfr_hi95, $
                          lo68Index=tsfr_lo68, hi68Index=tsfr_hi68, $
                          lo95Index=tsfr_lo95, hi95Index=tsfr_hi95)

     endif ;; nomass

;;------------------------------------------------------------
;; calculate minimum value of chisq directly and get indexes:
     x = (where(chisq eq min(chisq)))[0]
;; Define indices of the best-fit model for each index:

     itau = x / n_metal / n_delta / n_ebv / n_age
     imetal = (x - itau*n_age*n_ebv*n_delta*n_metal) / n_delta / n_ebv / n_age
     idelta = (x - itau*n_age*n_ebv*n_delta*n_metal $
               - imetal*n_age*n_ebv*n_delta) $
              / n_ebv / n_age
     iebv = (x - itau*n_age*n_ebv*n_delta*n_metal $
             - imetal*n_age*n_ebv*n_delta $
             - idelta*n_age*n_ebv) / n_age
     iage = x - itau*n_age*n_ebv*n_delta*n_metal $
            - imetal*n_age*n_ebv*n_delta $
            - idelta*n_age*n_ebv $
            - iebv*n_age
     
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
         
        result = {y_age: y_age, $
                  y_ebv:y_ebv, $
                  y_metal: y_metal, $
                  y_delta:y_delta, $
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
                  delta: {mean: delta_mean, sigma: delta_sigma, $
                          mean_ind: tdelta_mean,$
                          minchisq_ind: idelta, minchisq: deltaArr[idelta], $
                          maxpdf_ind: tdelta_median, maxpdf: deltaArr[tdelta_best], $
                          median: delta_median, median_ind: tdelta_median, $
                          lo68_ind: tdelta_lo68, hi68_ind: tdelta_hi68, $
                          lo95_ind: tdelta_lo95, hi95_ind: tdelta_hi95, $
                          lo68: delta_lo68, hi68: delta_hi68, $
                          lo95: delta_lo95, hi95: delta_hi95 }, $
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

     endelse
  endelse

end



