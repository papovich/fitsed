
pro fitsed_write_catalog, data=data,$
                          catalog=catalog, fitsed_cat, outdir, $
                          bestfit=bestfit, $ ; write only best fit values
                          paramfile=paramfile, $ ; copied to header of fitsed_cat
                          fileHeader=fileHeader, $ ; header used to make .sav files
                          name_zphot=name_zphot, $ ; passed to read_catalog
                          name_zspec=name_zspec, $ ; passed to read_catalog
                          ab_zeropoint=ab_zeropoint ; passed to read_catalog
  ;; last three keywords (all passed to _read_catalog) not really
  ;; needed... 
  ;;
  ;; outdir must have '/' trailing or things will go nuts... 
  if not keyword_set(fileHeader) then fileHeader=''

  if keyword_set(catalog) then $
     data = fitsed_read_catalog(catalog,$
                                name_zphot=name_zphot,$
                                name_zspec=name_zspec,$
                                ab_zeropoint=ab_zeropoint)
  if total(size(data)) eq 0 then $
     message,'ERROR: either data, or catalog keyword need to be specified.'
        
  id = data.id & usez=data.z

  ; open output file and populate header:

  openw,lun, /get_lun,fitsed_cat
  printf,lun, '# FITSED v2.0 output file'
  printf,lun, '#'
  if keyword_set(paramfile) then begin
     n=file_lines(paramfile)
     lines=strarr(n)
     openr,lun2, paramfile, /get_lun
     readf,lun2,lines
     close, lun2
     free_lun, lun2
     for i=0,n_elements(lines)-1 do begin
        printf,lun,'# '+lines[i]
     endfor
  endif
  ; printf,lun, '# all parameters will be written here'
  if not keyword_set(bestfit) then begin
     columns = ['id', 'z', 'tau', 'tau_l68', 'tau_h68', $
                'metal', 'met_l68', 'met_h68', $
                'lg_age', 'lg_a_l68', 'lg_a_h68', $
                'ebv', 'ebv_l68', 'ebv_h68', $
                'delta', 'del_l68', 'del_h68', $
                'lg_mass', 'lg_mass_l68', 'lg_mass_h68', $
                'lg_sfr', 'lg_sfr_l68', 'lg_sfr_h68']
     myformat=['a6', 'a5', 'a7','a9','a9', $ ; id, z, tau
               'a7', 'a9', 'a9', $           ; metal
               'a7', 'a9', 'a9',$            ; lg_age
               'a6','a8','a8', $             ; ebv
               'a6','a8','a8', $             ; delta
               'a8','a12','a12', $           ; lg_mass
               'a7','a11','a11']             ; SFR
  endif else begin
     columns = ['id', 'z', 'tau', $
                'metal', $
                'lg_age', $
                'ebv', $
                'delta', $
                'lg_mass', $
                'lg_sfr']
     myformat=['a6', 'a5', 'a8', $ ; id, z, tau
               'a9',  $           ; metal
               'a8', $            ; lg_age
               'a8', $             ; ebv
               'a8', $             ; delta
               'a9', $           ; lg_mass
               'a9']             ; SFR
  endelse     
  colstr = ''
  for i=0,n_elements(columns)-1 do colstr = colstr+string(format='('+myformat[i]+',x)',columns[i])
; FIX THE ABOVE
          
  printf,lun, '# '+colstr

  if not keyword_set(bestfit) then begin
     format = '(i8,x,f6.3,x, 3(f7.3,x), 3(f7.4,x), 3(f7.2,x), 3(f6.2,x), 3(f7.2,x), 3(f8.2,x), 3(f9.3,x))'
  endif else begin
     format = '(i8,x,f6.3,x, 1(f7.3,x), 1(f7.4,x), 1(f7.2,x), 1(f6.2,x), 1(f7.2,x), 1(f8.2,x), 1(f9.3,x))'
  endelse

  for i=0,n_elements(id)-1 do begin
     delvarx, result
     restore,outdir+fileHeader+strn(id[i])+'.sav'

     if not keyword_set(bestfit) then begin
        if total(size(result)) gt 5 then begin
           printf,lun,format=format, $
                  id[i], usez[i],  $
                  result.tau.median, result.tau.lo68, result.tau.hi68, $
                  result.metal.median, result.metal.lo68, result.metal.hi68, $
                  result.log_age.median,result.log_age.lo68, result.log_age.hi68, $
                  result.ebv.median, result.ebv.lo68, result.ebv.hi68, $
                  result.delta.median, result.delta.lo68, result.delta.hi68, $
                  alog10(result.mass.median), alog10( result.mass.lo68), alog10(result.mass.hi68), $
                  alog10(result.sfr.median), alog10( result.sfr.lo68), alog10(result.sfr.hi68)
        endif else begin
           printf,lun,format=format, $
                  id[i], -1,  $
                  replicate(0.0,21)
        endelse     
     endif else begin
        if total(size(result)) gt 5 then begin
           printf,lun,format=format, $
                  id[i], usez[i],  $
                  result.tau.maxpdf, $
                  result.metal.maxpdf, $
                  result.log_age.maxpdf,$
                  result.ebv.maxpdf, $
                  result.delta.maxpdf, $
                  alog10(result.mass.maxpdf),$
                  alog10(result.sfr.maxpdf)
        endif else begin
           printf,lun,format=format, $
                  id[i], -1,  $
                  replicate(0.0,7)
        endelse     
     endelse


  endfor

  close,lun
  free_lun,lun

end

