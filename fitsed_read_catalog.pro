
function fitsed_read, catalog, header=header

; this routine reads in the catalog file, strips out header and data
; and returns a big data matrix (all type float).

; borrowed heavily from fast... 


; % TO HERE - REPLACE WITH NUM_FILT ROUTINE FROM READ_FILTERS_ TO USE
;             TRANSLATE FILE TO FIGURE OUT FILTERS AND ORDER TO USE
;             FOR PHOTOMETRY... 
print, "Reading catalog: "+catalog

header=''
n=file_lines(catalog)
lines=strarr(n)
openr,lun, catalog, /get_lun
readf,lun,lines
close, lun
free_lun, lun

; read header into header lines.  Assume that the last header line will contain
; column information.

first_char = strmid(lines,0,1)
ihdr = where(first_char eq '#', nhdr, comp=idata, ncomp=ndata)

if ndata eq 0 then message, 'ERROR: '+catalog+' contains no data (or something else is messed up)'
ncol = n_elements(strsplit(lines[idata[0]]))
; strip last line of header to get column information:

; parse header
if nhdr gt 0 then begin
  h = strmid(lines[ihdr],1)   ; strip comment symbol
  for i=0,n_elements(h)-1 do begin
      head_line = strsplit(h[i],/extract)
      nhcol     = n_elements(head_line)
      if nhcol eq ncol then begin
          if n_elements(header) eq 0 or n_elements(header) eq 1 then $
            header = head_line else header = [[header],[head_line]]
      endif
  endfor
endif

data = make_array(ncol, ndata, type=type)
for i=0L, ndata-1 do begin
    if n_elements(strsplit(lines[idata[i]],/extract)) ne ncol then begin
       print,'ERROR: when reading '+catalog
       print,'       number of columns is not constant throughout file'
       exit
    endif
    data[*,i] = fix( strsplit(lines[idata[i]],/extract), type=4) ; type 4 
endfor
 
return, data

END

;------------------------------------------------------------

function fitsed_read_catalog, catalog, header=header, $
                              ab_zeropoint=ab_zeropoint, $
                              name_zphot=name_zphot, $
                              name_id=name_id, $
                              name_zspec=name_zspec

if KEYWORD_SET(NAME_ZPHOT) then name_zphot=strupcase(name_zphot) $
   else name_zphot = strupcase('z_phot')
message,/cont,'will use column name '+name_zphot+' for photo-z'

if keyword_set(NAME_ZSPEC) then begin
   NAME_ZSPEC=strupcase(name_zspec)
   message,/cont,'but, will use column name '+name_zspec+' for redshift if '+name_zspec+'>0.0'
endif
; if name_zspec is set, then code will use this column for redshift,
; if name_zspec column is > 0.0


if not KEYWORD_SET(NAME_ID) then name_id='ID'
if not KEYWORD_SET(AB_ZEROPOINT) then AB_ZEROPOINT = 25.  ;; for consistency with EAZY and FAST-like catalogs

transcat = repstr(catalog,'.cat','.translate')
zphotcat = repstr(catalog,'.cat','.zout')

; find the ID column and set that type to long.  All others are float
;type = replicate('flt',ncol)
;type[where( strcmp(strupcase(header),name_id),/null)]='lon'

cat = fitsed_read(catalog, header=header)
zcat = fitsed_read(zphotcat, header=hzcat)

; check that ID columns match up:
cid = (where( strcmp(strupcase(header),name_id) eq 1,/null))[0]
zid = (where( strcmp(strupcase(hzcat),name_id) eq 1, /null))[0]
if total(abs(cat[cid]-zcat[zid])) ne 0 then $
   message,'ERROR: flux catalog and zphot catalog have different (or unmatched) rows'

; create structure to hold photometry.  Place F* values in phot array
; and E* values in dphot array

; find zphot column:
zid = (where( strcmp(strupcase(hzcat),name_zphot) eq 1))[0]
if keyword_set(name_zspec) then $
   zspid = (where( strcmp(strupcase(header),name_zspec) eq 1))[0]

if zid[0] eq -1 then $
   message,'ERROR: cannot find zphot column with label '+name_zphot+' in photo-z catalog '+zphotcat

; find all the photometry columns:
n_filt = fitsed_read_translate(catalog, fl_ind=f1, err_ind=e1)
;f1 = where( strcmp(strupcase(strmid(header,0,2)),'F_') eq 1)
;e1 = where( strcmp(strupcase(strmid(header,0,2)),'E_') eq 1)

if n_elements(f1) ne n_elements(e1) or $
   (n_elements(f1) eq 0) then $
      message, 'ERROR: number of flux and error columns are not the same (or zero) in catalog file '+catalog

; apply ab_zeropoint to get all fluxes to uJy (assume uJy): 
mag2ab = 10d^(-0.4*(ab_zeropoint-23.9))
MESSAGE,/cont,'applying factor of '+strn(mag2ab)+' to convert catalog fluxes to uJy'

; create data structure for fitting:
szcat = size(/dim,cat)
if size(/dim, szcat) eq 1 then useSize=1 else useSize=szcat[1]
data = replicate( { id: 0L, z: 0.0, $
                    phot: fltarr(n_elements(f1)), $
                    dphot: fltarr(n_elements(e1))}, useSize)
data.id = reform(cat[cid,*])
data.z = reform(zcat[zid,*])
if keyword_set(name_zspec) then begin
   if zspid ne -1 then begin
      x = where( cat[zspid,*] gt 0.0,/null)
      if x ne !NULL then $
         data[x].z = reform(cat[zspid,x])
   endif
endif
data.phot = cat[f1,*] * mag2ab
data.dphot = cat[e1,*] * mag2ab

return, data

end


