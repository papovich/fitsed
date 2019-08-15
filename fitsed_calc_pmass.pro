; take a lut file and save file from FITSED and return the P(Mass)

; this is cut out of: 
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

; EXAMPLE:

; restore,'fitsed_z1_delayed_calzetti_Zsol/3647.sav',/verb
; y = fitsed_calc_pmass(mass, pdf, mst=mst)  
; plot, mass[mst], gauss_smooth(y[mst],20), /xlog, xr=[1e10, 1e12],/xst


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


; returns Y(Mass) = P(Mass)
function fitsed_calc_pmass, mass, pdf, mst=mst, $
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
