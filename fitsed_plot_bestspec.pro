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
  print, '%    metal= '+strn(result.metal.minchisq)
  print, '%    E(B-V)= '+strn(result.ebv.minchisq)
  print, '%    delta= '+strn(result.delta.minchisq)

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

     case result.tau.minchisq of
        -100 : taustr[i] = 'neg100g'
        -10 : taustr[i] = 'neg10g'
        -1 : taustr[i] = 'neg1g'
        -70 : taustr[i] = 'neg70g'
        -7 : taustr[i] = 'neg7g'
        -0.7 : taustr[i] = 'neg700m'
        -50 : taustr[i] = 'neg50g'
        -5 : taustr[i] = 'neg5g'
        -0.5 : taustr[i] = 'neg500m'
        -30 : taustr[i] = 'neg30g'
        -3 : taustr[i] = 'neg3g'
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


     isedfile = modelHdr+'_'+metalStr+'_tau'+taustr+'yr.ised'
        
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
  restore,lutfile; ,/skip_existing
  ;; should give you lut, zed, etc... 
                                ; need to figure out best model here!
  zind = (where( abs(zed-data[x].z) eq min(abs(zed-data[x].z))))[0]
  i = result.log_age.minchisq_ind
  j = result.ebv.minchisq_ind
  d =result.delta.minchisq_ind
  k = result.tau.minchisq_ind
  l = result.metal.minchisq_ind
  myvec = (lut[zind,i,j,d,l,k,*])[*]

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

     if keyword_set(label) then $
        al_legend, [    'ID '+strn(data[x].id)+', z= '+strn(data[x].z), $
                        'Log M/Msol= '+strn(alog10(result.mass.minchisq)), $
                        'Log age/yr= '+strn(result.log_age.minchisq), $
                        'tau/Gyr= '+strn(result.tau.minchisq), $
                        'metal= '+strn(result.metal.minchisq), $
                        'E(B-V)= '+strn(result.ebv.minchisq), $
                        'delta= '+strn(result.delta.minchisq)], charsize=1.0

     
     data=data[x]
     phot=data[x].phot
     dphot=data[x].dphot

     if keyword_set(stop) then stop
     
  endif
  
end
