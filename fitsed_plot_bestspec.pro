; take: 
; param file (models used)
; filters res file (while filters used)
; fitsed_output/*.sav (best fit for object) - get from id
; 
; file and plot the best fitting model spectrum over the photometry


function findel, x, arr
  
  return, (where( abs( x - arr) eq min(abs(x - arr))))[0]
  
end


pro fitsed_plot_bestspec, id, $
                          par=par, modelspec=modelspec, $
                          modelphot=modelphot, $
                          lambda=lambda, $
                          phot=phot, dphot=dphot, $
                          data=data, $
                          result=result, $
                          nophot=nophot, $ ; dont plot real phot
                          nomodelphot=nomodelphot, $ ; dont plot model phot
                          noplot=noplot, $           ; if true then do not plot!
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
  print, '% Plotting best fit for object= '+strn(data[x].id)+', z= '+strn(data[x].z)
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

  ;; generate best-fit spectrum from result file:
  modelHdr = p.modeldir+p.imf
  modelTail='.ised'

  if strcmp(p.ssp,'bc03') then begin
     if strcmp(p.imf,'chab') then begin
        chabrier=1
        endian="big"
     endif else begin
        chabrier=0
        delvarx, endian
     endelse

     case result.metal.minchisq of
        0.02 : metalStr = 'zsol'
        0.008 : metalStr = 'z0p4sol'
        0.004 : metalStr = 'z0p2sol'
        0.0004 : metalStr = 'z0p02sol'
        0.05 : metalStr='z2p5sol'
        else  : begin
           message, 'ERROR, result.metal.maxpdf value of '+$
                    string(result.metal.maxpdf)+' not recognized, stopping'
        end
     endcase

                                ; ALPHA and BETA needed for DPL SFH
     if strcmp(p.SFH,'dpl') then begin
        tstr = strtrim(string(format='(f5.2)', result.tau.minchisq),2)
        taustr = repstr(tstr, '.','p')
        if result.alpha.minchisq lt 2 then begin
              alphastr = strtrim(string(format='(f3.1)', result.alpha.minchisq),2)
              alphastr = repstr(alphastr,'.','p')
           endif else if result.alpha.minchisq le 1000 then begin
              alphastr[i] = strtrim(string(format='(I4)', result.alpha.minchisq),2)
           endif
           if result.beta.minchisq lt 2 then begin
              betastr = strtrim(string(format='(f3.1)', result.beta.minchisq),2)
              betastr = repstr(betastr,'.','p')
           endif else if result.beta.minchisq le 1000 then begin
              betastr = strtrim(string(format='(I4)', result.beta.minchisq),2)
           endif
        endif else begin
        case result.tau.minchisq of
           -100 : taustr = 'neg100g'
           -10 : taustr = 'neg10g'
           -1 : taustr = 'neg1g'
           -70 : taustr = 'neg70g'
           -7 : taustr = 'neg7g'
           -0.7 : taustr = 'neg700m'
           -50 : taustr = 'neg50g'
           -5 : taustr = 'neg5g'
           -0.5 : taustr = 'neg500m'
           -30 : taustr = 'neg30g'
           -3 : taustr = 'neg3g'
           -0.3 : taustr[i] = 'neg300m'
           100 : taustr = '100g'
           10 : taustr = '10g'
           1 : taustr = '1g'
           0.1 : taustr = '100m'
           0.01 : taustr = '10m'
           0.001 : taustr = '1m'
           2 : taustr = '2g'
           0.2 : taustr = '200m'
           0.02 : taustr = '20m'
           3 : taustr = '3g'
           0.3 : taustr = '300m'
           0.03 : taustr = '30m'
           5 : taustr = '5g'
           0.5 : taustr = '500m'
           0.05 : taustr = '50m'
           7 : taustr = '7g'
           0.7 : taustr = '700m'
           0.07 : taustr = '70m'
           else : begin
              message, 'ERROR, result.tau.maxpdf value of '+string(result.tau.maxpdf)+' not recognized, stopping'
           end
        endcase
     endelse

     if strcmp(p.sfh,'exp') then $
        isedfile = modelHdr+'_'+metalStr+'_tau'+taustr+'yr.ised' $
     else if strcmp(p.sfh,'delayed') then $
        isedfile = modelHdr+'_'+metalStr+'_delayedtau'+taustr+'yr.ised' $
     else begin
        if strcmp(p.sfh,'dpl') then begin
           isedfile = modelHdr+'_'+metalStr+'_dpl_tau'+taustr+$
                      '_alpha'+alphastr+'_beta'+betastr+'.ised'
        endif else begin
           print,'dying... no SFH specified or recognized.  sfh=',sfh
           stop
        endelse
     endelse

     ;isedfile = modelHdr+'_'+metalStr+'_tau'+taustr+'yr.ised'
        
     fitsed_bcspectrum, isedfile, data[x].z, $
                        ;result.mass.maxpdf, $
                        result.mass.minchisq, $
                        10d^result.log_age.minchisq, $
                        result.ebv.minchisq, $
                        chabrier=chabrier, endian=endian, $
                        nebular_fesc=p.nebular_fesc, $
                        /nebular_nolya, $
                        extinction_law=p.extinction_law, $
                        outspectrum=modelspec
     
  endif else begin
     message,'Unrecognized SSP in param file: SSP= ',p.ssp
     return
  endelse

  ;; print photometry, filter wavelengths and die
  ;print, data[x].phot, data[x].dphot
  ;print, lambda

  ;; figure out model photometry from LUT file:
  restore,p.lutfile; ,/skip_existing
  ;; should give you lut, zed, etc... 
                                ; need to figure out best model here!

  if strcmp(p.sfh,'dpl') then begin
     zind = (where( abs(zed-data[x].z) eq min(abs(zed-data[x].z))))[0]
     i = result.log_age.minchisq_ind
     j = result.ebv.minchisq_ind
;  d =result.delta.minchisq_ind
     k = result.tau.minchisq_ind
     l = result.metal.minchisq_ind
     a=result.alpha.minchisq_ind
     b=result.beta.minchisq_ind
     myvec = (lut[zind,i,j,l,k,a,b,*])[*]
  endif else begin
     zind = (where( abs(zed-data[x].z) eq min(abs(zed-data[x].z))))[0]
     i = result.log_age.minchisq_ind
     j = result.ebv.minchisq_ind
;  d =result.delta.minchisq_ind
     k = result.tau.minchisq_ind
     l = result.metal.minchisq_ind
     myvec = (lut[zind,i,j,l,k,*])[*]
  endelse

  mfactor = 1000d               ;* 10d^(0.4*(23.9-p.AB_ZEROPOINT))               ;                             ;for uJy
  myscale = result.mass.minchisq / mfactor
  model_phot = myvec*myscale

  if not keyword_set(noplot) then begin
  
     ang='!3'+string("305b) & ang='!3'+string("305b) ;; repeated just to get " level correct

     yrange = minmax( modelspec[1,where(modelspec[0,*] ge 1300*(1+data[x].z) and $
                                        modelspec[0,*] le 80000)]) * [0.5,5.0]
     if keyword_set(ps) then begin
        plot, modelspec[0,*], modelspec[1,*], tit='!5', $
              xr=[1000, 100000],/xlog,/ylog, yr=yrange, $
              xtit='observed wavelength [!3'+ang+'!5]', ytit='flux density [!7l!5Jy]', $
              _EXTRA=_EXTRA
     endif else begin
        plot, modelspec[0,*], modelspec[1,*], $
              xr=[1000, 100000],/xlog,/ylog, yr=yrange, $
              xtit='observed wavelength ['+ang+']', ytit='flux density [!4l!3Jy]', $
              _EXTRA=_EXTRA
     endelse
     if not keyword_set(ps) then modelcolor='yellow' else modelcolor='black'
     if not keyword_set(nomodelphot) then $
        oplot, lambda, model_phot, $
               color=djs_icolor(modelcolor), psym=6, symsize=3
     modelphot=model_phot
     if not keyword_set(nophot) then $
        oploterror, lambda, data[x].phot, data[x].dphot, $
                    color=djs_icolor('red'), psym=6, symsize=2

     if keyword_set(label) then begin
        if strcmp(p.sfh,'dpl') then begin
           al_legend, [    'ID '+strn(data[x].id)+', z= '+strn(data[x].z), $
                           'Log M/Msol= '+strn(alog10(result.mass.minchisq)), $
                           'Log age/yr= '+strn(result.log_age.minchisq), $
                           'tau/Gyr= '+strn(result.tau.minchisq), $
                           '!7a!5= '+strn(result.alpha.minchisq), $
                           '!7b!5= '+strn(result.beta.minchisq), $
                           'metal= '+strn(result.metal.minchisq), $
                           'E(B-V)= '+strn(result.ebv.minchisq), $
                           '!7v!5!s!u2!n!r!d0!n='+strn(result.minchisq)], charsize=1
        endif else begin
           al_legend, [    'ID '+strn(data[x].id)+', z= '+strn(data[x].z), $
                           'Log M/Msol= '+strn(alog10(result.mass.minchisq)), $
                           'Log age/yr= '+strn(result.log_age.minchisq), $
                           'tau/Gyr= '+strn(result.tau.minchisq), $
                           'metal= '+strn(result.metal.minchisq), $
                           'E(B-V)= '+strn(result.ebv.minchisq), $
                           '!7v!5!s!u2!n!r!d0!n='+strn(result.minchisq)], charsize=1
                                ;'delta=
                                ;'+strn(result.delta.minchisq)],
                                ;charsize=1.0
        endelse
        
     endif
     
     data=data[x]
     phot=data[x].phot
     dphot=data[x].dphot
     
     if keyword_set(stop) then stop
     
  endif
  
  
end
  
 
