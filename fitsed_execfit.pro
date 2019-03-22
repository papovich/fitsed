
pro fitsed_execfit, p, data=data
  ;; p is the param structure, loaded from 
  ;; fitsed_read_param, paramfile, params=p
  
  if ~file_test(p.outdir,/directory) then file_mkdir,p.outdir
;
; read in the catalog: 
  data = fitsed_read_catalog(p.catalog,$
                             name_zphot=p.name_zphot,$
                             name_zspec=p.name_zspec,$
                             ab_zeropoint=p.ab_zeropoint)

  ;; uncomment this eventually:    delvarx(version)
  restore,p.lutfile
  ;; comment   this out later:
  version=p.version
  ;; 
  if p.version ne version then begin
     print,'FITSED::   WHOA: Running version=',p.version, ' but LUT created using version=',version
     print,'FITSED::   Dying...'
     stop
  endif
  lambda_filters= lambda

  rest_lowerlimit = p.rest_lowerlimit
;;====
;; turn off floating precison error reporting:
  currentExcept= !Except
  !Except=0
  void = Check_Math()
  floating_point_underflow = 32
  status = Check_Math()
  IF(status AND NOT floating_point_underflow) NE 0 THEN $
     Message, 'IDL Check_Math() error: ' + StrTrim(status, 2)
;;===

  for i=0,n_elements(data)-1 do begin

     id = data[i].id
     phot = data[i].phot
     dphot = data[i].dphot
     dphot = sqrt( dphot^2 + (phot*p.err_flux_factor)^2)

     usez = data[i].z
   
     print,'working on object '+strn(id)+' with z='+strn(usez)

     fitsed_fitchisq, phot,dphot, usez, lambda_filters, $
                      lut, zed, log_ageArr, ebvArr, metalArr, tauArr, alphaArr, betaArr, $
                      lutMstar, lutSFR, $
                      tltuniverse=p.ageltuniverse,$
                      pdf=pdf, mass=mass, sfr=sfr, sfh=p.sfh, result=result, rest_lowerlimit=rest_lowerlimit
     filename = p.outsaveHeader+strn(id)+'.sav'
;   print, outdir+filename, phot, phot/dphot
     lutfile=p.lutfile
     if p.write_sav then save,file=p.outdir+filename, result, mass, sfr, pdf, lutfile
endfor

end
