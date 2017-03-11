
; read in the filters from the .RES file
; 
; butchered from FAST

function fitsed_read_filters, filters_res, lambda=lambda, catalog=catalog, $
                              bandpass_indexes=bandpass_indexes, $
                              stop=stop,$
                              VERBOSE=VERBOSE
  ;; read in the filters_res file, extract the filters you need, then
  ;; dump return the filters in the format needed by fitsed.pro

  filtersDone = 0
  if keyword_set(CATALOG) then begin
     if file_test(CATALOG) then begin
        num_filt = fitsed_read_translate(CATALOG,fnames=tempFnames)
        filtersDone=1
     endif 
  endif 
  if ~filtersDone and keyword_set(bandpass_indexes) then begin
     ;; check if bandpass_indexes exists and is integer array or
     ;; strings
     szbands = size(bandpass_indexes)
     if szbands[2] eq 2 then $
        num_filt = string(bandpass_indexes) $
     else if szbands[2] eq 7 then num_filt = bandpass_indexes $
     else begin
        print,'% FITSED_READ_FILTERS: ERROR -'
        print,'% FITSED_READ_FILTERS: BANDPASS_INDEXES MUST BE ARRAY OF INTEGERS OR STRINGS'
        print,'% FITSED_READ_FILTERS: stopping... '
        stop
     endelse
     filtersDone=1
  endif 

  if ~filtersDone then begin
     message,'KEYWORD CATALOG OR BANDPASS_INDEXES MUST BE SPECIFIED...AND EXIST/BE LEGAL. STOPPING....'
     return,0
  endif

  n_filt     = n_elements(num_filt)
  tags = strarr(n_filt) ;; and array with the name for each filter.

;   if not file_test(FILTERS_RES) then begin
;      print,"ERROR: defined filter file is not available" & exit
;   endif

  ;; read in the whole fugly file
  READCOL,FILTERS_RES,num,wl,tr,FORMAT='I,F,F',/silent
  infofile=FILTERS_RES+'.info'
  if strmatch(FILTERS_RES,'*.latest') then begin
     ;; filters.res ends in .latest, so change it to .info:
     infofile=repstr(FILTERS_RES,'.latest','.info')
  endif
  ;message,'% ERROR: cannot locate filter info file '+FILTERS_RES+'.info'
  if file_test(infofile) then begin
     READCOL,infofile,flen,fnames, $
             delimit=' =', FORMAT='i,a',/silent
  endif else fnames=tempFnames
  ;; 
  st_filt    = where(num eq 1,n_filt_res)

  en_filt    = [st_filt[1:n_filt_res-1]-1,n_elements(no)-1]
  nel_filt   = en_filt-st_filt+1
  lambda     = fltarr(n_elements(num_filt))
  c_num_filt  = num_filt-1
;  tags = fnames[fix(num_filt)]
  translate=intarr(n_elements(num_filt))
  for i=0,n_filt-1 do begin
     if file_test(infofile) then begin
        translate[i] = (where( fix(num_filt[i]) eq flen))[0]
        tstr = fnames[translate[i]]
     endif else $
        tstr=fnames[i]
     tstr = repstr(tstr,'/','')
     tstr = repstr(tstr,'+','')
     tstr = repstr(tstr,'_','')
     tstr = repstr(tstr,'-','')
     tstr = repstr(tstr,'.dat','')
     tstr = repstr(tstr,'.txt','')
     tstr = repstr(tstr,'ccd','')
     tstr = repstr(tstr,'2cols','')
     tstr = repstr(tstr,'.res','')
     tstr = repstr(tstr,'.raw','')
     tstr = repstr(tstr,'.','')
     tags[i] = tstr
  endfor
  ;; if any tag starts with a number, put a '_' in front of it to make
  ;; it a legimate structure tag:
  t = where(strmatch(tags,'*[1234567890]'))
  if t[0] ne -1 then begin
     for i=0,n_elements(t)-1 do begin
        tags[t[i]] = '_'+tags[t[i]]
     endfor
  endif

  for i=0,n_filt-1 do begin
     tmp_filt = [[wl[st_filt[c_num_filt(i)]:en_filt[c_num_filt[i]]]],$
                 [tr[st_filt[c_num_filt(i)]:en_filt[c_num_filt[i]]]]]
     ;; create structure to hold all the filters:
     if keyword_set(VERBOSE) then print, '% FITSED_READ_FILTERS: READING FILTER #, name= '+num_filt[i]+' into tag= '+tags[i]
     if i eq 0 then $
        filters=create_struct(tags[i], transpose(tmp_filt,[1,0])) $
     else $
        filters=create_struct(filters, tags[i], transpose(tmp_filt,[1,0]))


     ;; create array with effective wavelengths for each filter
     ;; TODO - SHOULD MAKE THIS USE TSUM, NOT TOTAL
     lambda[i] = TOTAL(tmp_filt[*,0]*tmp_filt[*,0]*tmp_filt[*,1]) / $
                 TOTAL(tmp_filt[*,0]*tmp_filt[*,1])
  endfor

  if n_elements(n_filt) eq 0 then n_filt  = 0
  if keyword_set(stop) then stop
  return, filters
end

