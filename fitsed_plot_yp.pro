; take: 
; param file (models used)
; filters res file (while filters used)
; fitsed_output/*.sav (best fit for object) - get from id
; 
; file and plot the best fitting model spectrum over the photometry


function findel, x, arr
  
  return, (where( abs( x - arr) eq min(abs(x - arr))))[0]
  
end


pro fitsed_plot_yp, id, $
                    par=par, $
                    result=result, $
                    noplot=noplot, $ ; if true then do not plot!
                    ps=ps, label=label, $
                    _EXTRA=_EXTRA, $
                    stop=stop

  ;; modelspec is returned and contains the best-fit model spectrum
  ;; lambda contains the wavelengths of the filters
  ;; phot and dphot contain the photometry and errors 

  fitsed_read_param, par, params=p
  data = fitsed_read_catalog(p.catalog,$
                             name_zphot=p.name_zphot,$
                             name_zspec=p.name_zspec,$
                             ab_zeropoint=p.ab_zeropoint)

  filters = fitsed_read_filters(p.bandpass_file,$
                                lambda=lambda, $
                                catalog=p.catalog, $
                                bandpass_indexes=bandpass_indexes)
  
  x = where( data.id eq id)
  if x[0] eq -1 then $
     message,'ERROR: id= '+strn(id)+' not found in catalog: '+p.catalog

  ;; restore .sav file:
  restore,p.outdir+strn(long(id))+'.sav'
  ;; restore lut file:
  restore,p.lutfile, /verb
  
   print, '% Plotting P(x) for object= '+strn(data[x].id)+', z= '+strn(data[x].z)
   print,'% Best-fit results:'
  print, '%    Log M/Msol= '+strn(alog10(result.mass.minchisq))
  print, '%    Log age/yr= '+strn(result.log_age.minchisq)
  print, '%    tau/Gyr= '+strn(result.tau.minchisq)
  if strcmp(p.sfh,'dpl') then begin
     print, '%    alpha= '+strn(result.alpha.minchisq)
     print, '%    beta= '+strn(result.beta.minchisq)
  endif
  print, '%    metal= '+strn(result.metal.minchisq)
  print, '%    E(B-V)= '+strn(result.ebv.minchisq)
  print, '%    min Chi-squared= '+strn(result.minchisq)
;  print, '%    delta= '+strn(result.delta.minchisq)

  ;; figure out model photometry from LUT file:
  restore,p.lutfile; ,/skip_existing
  ;; should give you lut, zed, etc... 
                                ; need to figure out best model here!

;------------------------------------------------------------
;------------------------------------------------------------
; CALCULATE Y VALUES: 
;------------------------------------------------------------
;------------------------------------------------------------
  ymass = fitsed_calc_pmass(mass, pdf, mst=mst)
  ysfr = fitsed_calc_pmass(sfr, pdf, mst=sft)

  y_age = result.y_age
                                ; fitsed_yp6(10d^log_ageArr, pdf, 0, a=ebvArr, b=metalArr, c=tauArr, d=alphaArr, e=betaArr,$
                                ;                 norm=norm)
; ;apply the normalization:
                                ;    pdf /= norm
;
  y_ebv=result.y_ebv
  y_metal=result.y_metal
  y_tau=result.y_tau
  if strcmp(p.sfh,'dpl') then begin
     y_alpha = result.y_alpha
     y_beta = result.y_beta
  endif
;------------------------------------------------------------

  
  if not keyword_set(noplot) then begin
     !p.charsize=1.5
     if strcmp(p.sfh,'dpl') then begin
        !p.multi=[0,3,3]
     endif else begin
        !p.multi=[0,3,2]
     endelse
     if keyword_set(ps) then begin
        plot, log_ageArr, y_age, tit='!5', xtit='log age/yr', ytit='P(age)'
     endif else begin
        plot, log_ageArr, y_age, xtit='log age/yr', ytit='P(age)'
     endelse
     plots, result.log_age.median, !y.crange, color=djs_icolor('red'), line=5
     plots, result.log_age.lo68, !y.crange, color=djs_icolor('red'), line=3
     plots, result.log_age.hi68, !y.crange,color=djs_icolor('red'), line=3

     ;; 
     plot, ebvArr, y_ebv, xtit='E(B-V)', ytit='P(E[B-V])'
     plots, result.ebv.median, !y.crange, color=djs_icolor('red'), line=5
     plots, result.ebv.lo68, !y.crange, color=djs_icolor('red'), line=3
     plots, result.ebv.hi68, !y.crange,color=djs_icolor('red'), line=3

     ;; 
     plot, metalArr, y_metal, xtit='Z', ytit='P(Z)'
     plots, result.metal.median, !y.crange, color=djs_icolor('red'), line=5
     plots, result.metal.lo68, !y.crange, color=djs_icolor('red'), line=3
     plots, result.metal.hi68, !y.crange,color=djs_icolor('red'), line=3

     ;; 
     plot, tauArr, y_tau, xtit='tau / Gyr', ytit='P(tau)', /xlog
     plots, result.tau.median, !y.crange, color=djs_icolor('red'), line=5
     plots, result.tau.lo68, !y.crange, color=djs_icolor('red'), line=3
     plots, result.tau.hi68, !y.crange,color=djs_icolor('red'), line=3
     if strcmp(p.sfh,'dpl') then begin
        plot, alphaArr, y_alpha, xtit='alpha', ytit='P(beta)', /xlog
        plots, result.alpha.median, !y.crange, color=djs_icolor('red'), line=5
        plots, result.alpha.lo68, !y.crange, color=djs_icolor('red'), line=3
        plots, result.alpha.hi68, !y.crange,color=djs_icolor('red'), line=3
        
        plot, betaArr, y_beta, xtit='beta', ytit='P(beta)', /xlog
        plots, result.beta.median, !y.crange, color=djs_icolor('red'), line=5
        plots, result.beta.lo68, !y.crange, color=djs_icolor('red'), line=3
        plots, result.beta.hi68, !y.crange,color=djs_icolor('red'), line=3
endif

     mr= [result.mass.lo68 / 5, result.mass.hi68 * 5]
     plot, mass[mst], gauss_smooth(ymass[mst],20), /xlog, xr=mr,/xst, $
           xtit='mass / Msol', ytit='P(Mass)'
     plots, result.mass.median, !y.crange, color=djs_icolor('red'), line=5
     plots, result.mass.lo68, !y.crange, color=djs_icolor('red'), line=3
     plots, result.mass.hi68, !y.crange,color=djs_icolor('red'), line=3

     sfrr= [result.sfr.lo68 / 5, result.sfr.hi68 * 5]
     plot, sfr[sft], gauss_smooth(ysfr[sft],20), /xlog, xr=sfrr,/xst, $
           xtit='sfr / Msol yr^-1', ytit='P(sfr)'
     plots, result.sfr.median, !y.crange, color=djs_icolor('red'), line=5
     plots, result.sfr.lo68, !y.crange, color=djs_icolor('red'), line=3
     plots, result.sfr.hi68, !y.crange,color=djs_icolor('red'), line=3

     if keyword_set(stop) then stop
     
  endif
  
  
end
  
 
