
function fitsed_read_translate, catalog, fl_ind=fl_ind, err_ind=err_ind, $
                                fnames=fnames

; based on old : ;function num_filters, CATALOG
; 
; this also will return fl_ind, which contains the indexes of which
; columns in the catalog file belong to the filters in the translate
; file.   err_ind contains the indexes of the error columns

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
  endif else begin
     message,'Cannot find translate file: ',translate
  endelse
 

  ;; capitalize everything:
  ;;header=strupcase(header)

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
  err_ind = where(strmatch(header,'E*[1234567890]') eq 1,n_filt)
  filt   = strmid(header[fl_ind,i_line],1)
  
  ; make an array with the filter names:
  fnames= strmid(tr1[where(strmatch(strupcase(tr1),'F_*'))],2)

  return,filt
end

;
;
;
