
; read in the filters from the .RES file
; 
; butchered from FAST

function num_filters, CATALOG

  ;; this routine strips out the filter #s you need
  ;; from the .translate file

  comment = '#'
  l_com   = strlen(comment)
  line    = ''
  w_head  = 0
  w_data  = 0

  TRANSLATE = repstr(CATALOG,'.cat','.translate')

  ;;...determine the number of data columns
  OPENR, lun, catalog, /GET_LUN
  while w_data eq 0 do begin
     READF, lun, line
     col   = strsplit(line, /extract)
     if strmid(col(0),0,l_com) ne comment then begin
        n_col  = n_elements(col)
        w_data = 1
    endif
  endwhile
  CLOSE, lun
  FREE_LUN, lun

  ;;...read header with same number of columns as the data
  OPENR, unit, CATALOG, /GET_LUN
  while w_head eq 0 do begin
     READF, unit, line
     col = strsplit(line, /extract)
     if strmid(col(0),0,l_com) eq comment and n_elements(col) gt 1 $
     then begin
                                ;...first entry may be attached to '#'
        if col(0) eq comment then tmp_head = col(1:(n_elements(col)-1)) $
        else tmp_head = [strmid(col(0),1),col(1:(n_elements(col)-1))]
        if n_elements(tmp_head) eq n_col then begin
            if n_elements(header) eq 0 then header = tmp_head else $
              header = [[header],[tmp_head]]
        endif
        
     endif else begin
        w_head = 1
     endelse
  endwhile
  CLOSE, unit
  FREE_LUN, unit

  ;;...change header file in case CATALOG.TRANSLATE is provided
  if FILE_TEST(translate) then begin
     readcol,TRANSLATE,tr1,tr2,format='(A,A)',/silent
     for i=0,n_elements(tr1)-1 do begin
        rep_h = where(header eq tr1(i),n_rep_h)
        if n_rep_h eq 1 then begin
           if (size(header))[0] eq 1 then begin
              header(rep_h) = tr2(i)
           endif else begin
              rep_h2 = array_indices(header,rep_h)
              header(rep_h2(0),rep_h2(1)) = tr2(i)
           endelse
        endif
     endfor
  endif

  ;;...determine filters from the header info
  if (size(header))[0] eq 0 then begin
     print,"ERROR: no header found in "+CATALOG
     print,"       Check that number of data columns is the same as header columns"
     exit
  endif
  if (size(header))[0] eq 1 then n_line_h = 1
  if (size(header))[0] eq 2 then n_line_h = (size(header))[2]
  header = REFORM(header,(size(header))[1],n_line_h)
  i_line = (array_indices(header,where(strmatch(header,'F*[1234567890]') eq 1)))[1,0]
  fl_ind = where(strmatch(header,'F*[1234567890]') eq 1,n_filt)
  filt   = strmid(header(fl_ind,i_line),1)
  
  return,filt

end

;
;
;

function fitsed_read_filters, filters_res, lambda=lambda, catalog=catalog, $
                              bandpass_indexes=bandpass_indexes, $
                              VERBOSE=VERBOSE
  ;; read in the filters_res file, extract the filters you need, then
  ;; dump return the filters in the format needed by fitsed.pro

  filtersDone = 0
  if keyword_set(CATALOG) then begin
     if file_test(CATALOG) then begin
        num_filt = num_filters(CATALOG)
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
  READCOL,FILTERS_RES,flen,fnames,fcom,flams,fabvegas,FORMAT='i,a,a,x,f,a',/silent

  ;; 
  st_filt    = where(num eq 1,n_filt_res)
  en_filt    = [st_filt[1:n_filt_res-1]-1,n_elements(no)-1]
  nel_filt   = en_filt-st_filt+1
  lambda     = fltarr(n_elements(num_filt))
  c_num_filt  = num_filt-1
;  tags = fnames[fix(num_filt)]
  for i=0,n_filt-1 do begin
     tstr = fnames[fix(num_filt[i])]
     tstr = repstr(tstr,'/','')
     tstr = repstr(tstr,'+','')
     tstr = repstr(tstr,'_','')
     tstr = repstr(tstr,'.dat','')
     tags[i] = tstr
  endfor

  for i=0,n_filt-1 do begin
     tmp_filt = [[wl[st_filt[c_num_filt(i)]:en_filt[c_num_filt[i]]]],$
                 [tr[st_filt[c_num_filt(i)]:en_filt[c_num_filt[i]]]]]
     ;; create structure to hold all the filters:
     if keyword_set(VERBOSE) then print, '% FITSED_READ_FILTERS: READING FILTER #, name= '+num_filt[i]+', '+fnames[i]+', into tag= '+tags[i]
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
  return, filters
end

