

function parseflt, str, type=type
  ;; take a str, strip out '[' and ',' and return an array of
  ;; type=type

  if not keyword_set(type) then type=4 ; flt, 
  ;; see: http://www.harrisgeospatial.com/docs/size.html#S_820040301_1082572
  
  spl = strsplit(str,'[,]',/ext)
  sz=size(/dim,spl)
  
  return, fix(spl,type=type)
end

;------------------------------------------------------------

pro fitsed_read_param, paramfile, params=params

  
  if file_test(paramfile) eq 0 then $
     message,'ERROR: parameter file, '+paramfile+', does not exist or unreadable'

  print, "Reading Parameter File: "+paramfile

  header=''
  n=file_lines(paramfile)
  lines=strarr(n)
  openr,lun, paramfile, /get_lun
  readf,lun,lines
  close, lun
  free_lun, lun

  bestoutput=''

  for i=0,n_elements(lines)-1 do begin
     first_char = strmid(lines[i],0,1)
     if first_char eq '#' then continue ; skip commmented lines

     ; split line, 1st var is keyword, 2nd is '=' and 3rd is value
     vars = strsplit(lines[i], /extr)
     if strcmp(vars[0],'') then continue ; blank line
     key = strupcase(vars[0])
     case 1 of 
        strcmp(key,'CATALOG') : catalog=repstr(vars[2],"'","")
        strcmp(key,'CATALOG_PATH') : catdir=repstr(vars[2],"'","")
        strcmp(key,'AB_ZEROPOINT') : ab_zeropoint=fix(vars[2],type=4)
        strcmp(key,'FILTERS_RES') : bandpass_file=repstr(vars[2],"'","")
        strcmp(key,'ERR_FLX_FACTOR') : err_flux_factor=fix(vars[2],type=4)
        strcmp(key,'NAME_ZPHOT') : name_zphot=repstr(vars[2],"'","")
        strcmp(key,'NAME_ZSPEC') : name_zspec=repstr(vars[2],"'","")

        strcmp(key,'LUT_FILE') : lutfile=repstr(vars[2],"'","")
        strcmp(key,'SAVEFILE_DIR') : outdir=repstr(vars[2],"'","")
        strcmp(key,'SAVEFILE_HEAD') : outsaveHeader=repstr(vars[2],"'","")
        strcmp(key,'WRITE_SAVE') : write_sav=fix(vars[2],type=2)
        strcmp(key,'OUTPUT_FILE') : output=repstr(vars[2],"'","")
        strcmp(key,'BESTFIT_OUTPUT_FILE') : bestoutput=repstr(vars[2],"'","")
        ;; strcmp(key,'BEST_FIT') : bestfit=fix(vars[2],type=2)
        strcmp(key,'SKIP_FIT') : skipfit=fix(vars[2],type=2)

        strcmp(key,'LIBRARY') : ssp=repstr(vars[2],"'","")
        strcmp(key,'LIBRARY_DIR') : modeldir=repstr(vars[2],"'","")
        strcmp(key,'IMF') : imf=repstr(vars[2],"'","")
        strcmp(key,'EXTINCTION_LAW') : extinction_law=repstr(vars[2],"'","")
        strcmp(key,'NEBULAR_FESC') : nebular_fesc=fix(vars[2],type=4)
        strcmp(key,'RESTLAM_LOWERLIMIT') : rest_lowerlimit=fix(vars[2],type=4)

        strcmp(key,'TAU') : tau = parseflt(vars[2],type=4)
        strcmp(key,'METAL') : metal=parseflt(vars[2],type=4)
        strcmp(key,'LOG_AGE_MIN') : log_agemin=fix(vars[2],type=4)
        strcmp(key,'LOG_AGE_MAX') : log_agemax=fix(vars[2],type=4)
        strcmp(key,'LOG_AGE_STEP') : log_agestep=fix(vars[2],type=4)
        strcmp(key,'AGELTUNIVERSE') : ageltuniverse=fix(vars[2],type=2)
        strcmp(key,'Z_MIN') : zmin=fix(vars[2],type=4)
        strcmp(key,'Z_MAX') : zmax=fix(vars[2],type=4)
        strcmp(key,'Z_STEP') : zstep=fix(vars[2],type=4)
        strcmp(key,'Z_LOG') : zlog=fix(vars[2],type=2)
        strcmp(key,'Z_CUSTOM') : zcustom=parseflt(vars[2],type=4)
        strcmp(key,'EBV_MIN') : ebvmin=fix(vars[2],type=4)
        strcmp(key,'EBV_MAX') : ebvmax=fix(vars[2],type=4)
        strcmp(key,'EBV_STEP') : ebvstep=fix(vars[2],type=4)
        strcmp(key,'DELTA_MIN') : deltamin=fix(vars[2],type=4)
        strcmp(key,'DELTA_MAX') : deltamax=fix(vars[2],type=4)
        strcmp(key,'DELTA_STEP') : deltastep=fix(vars[2],type=4)

        strcmp(key,'H0') : h=fix(vars[2],type=4)/100.
        strcmp(key,'OMEGA_M') : omega=fix(vars[2],type=4)
        strcmp(key,'OMEGA_L') : lambda=fix(vars[2],type=4)

        else : message, /cont, 'ignoring unknown keyword '+vars[0]
     endcase
     
  endfor

  ;; set extinction law: 
  if strcmp(extinction_law,'salmon') then begin
     print,"%% Salmon extinction law assumed"
     extinction_law='fitsed_salmon' 
  endif else $
     extinction_law='fitsed_calzetti' ;; default is calzetti

  if strcmp(output,'') then begin
     output=repstr(catalog,'.cat','.fitsed.cat')
  endif
  if strcmp(output,catalog) then begin
     print,'% FITSED_READ_PARAM: error output cat and input cat have same name, shame, shame, crashing.... '
     stop
  endif
  ;; set output - try to guess what might have to happen:
  

  if strcmp(bestoutput,'') then begin
     bestoutput = repstr(output,'.cat','.bestfit.cat')
     if strcmp(output,bestoutput) then begin
        print,'% FITSED_READ_PARAM: output and best-fit output file names have same name.  Make sure your ouput catalog ends in .cat (or it will crash)'
        stop
     endif
  endif

  if strcmp(lutfile,'') then $
     lutfile= 'fitsed_lut_'+repstr(catalog,'.cat','')+'_'+$
              ssp+'_'+imf+'_'+extinction_law+'.sav'

  defsysv, '!omega', exists=omega_exists
  if (omega_exists eq 0) then defsysv,'!omega', omega $
  else !omega=omega

  defsysv, '!lambda', exists=lambda_exists
  if (lambda_exists eq 0) then defsysv,'!lambda', lambda $
  else !lambda=lambda

  defsysv, '!h', exists=h_exists
  if (h_exists eq 0) then defsysv,'!h', h $
  else !h=h

;outdir='fitsed_output/'

if ~file_test(outdir,/directory) then file_mkdir,outdir

params = { ssp: ssp, $
           imf: imf, $
           modeldir: modeldir, $
           bandpass_file: bandpass_file, $
           zmin: zmin,  $
           zmax: zmax, $
           zstep: zstep, $
           zlog: zlog, $
           zcustom: zcustom, $
           ;; 
           extinction_law: extinction_law, $
           ebvmin: ebvmin, $
           ebvmax: ebvmax, $
           ebvstep: ebvstep, $
           ;;
           log_agemin: log_agemin, $
           log_agemax: log_agemax, $
           log_agestep: log_agestep, $
           ageltUniverse: ageltUniverse, $
           ;;
           tau: tau, $
           metal: metal, $
           ;; 
           nebular_fesc: nebular_fesc, $
           rest_lowerlimit: rest_lowerlimit,$
           ;; 
           err_flux_factor: err_flux_factor, $
           name_zphot: name_zphot, $
           name_zspec: name_zspec, $
           ab_zeropoint: ab_zeropoint, $
           ;;
           write_sav: write_sav, $
           outsaveHeader: outsaveHeader, $
           catdir: catdir, $
           catalog: catdir+catalog, $
           lutfile: lutfile, $
           fitsed_cat: catdir+output, $
           bestfit_cat: catdir+bestoutput, $
;           bestfit: bestfit, $
           skipfit: skipfit, $
           outdir: outdir $
         }
           

end
