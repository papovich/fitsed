
;; compile some useful functions: 


function findel, x, arr
  
  return, (where( abs( x - arr) eq min(abs(x - arr))))[0]
  
end


;; main exec file.  Calls everything else, given a paramfile

pro fitsed, paramfile

version = 3.0

print,'Running FITSED version = v'+strn(version)

fitsed_read_param, paramfile, params=p, version=version

if ~file_test(p.lutfile) then begin
   ;; generate Look Up Table (LUT): 

   delvarx, chabrier, endian
;   if strcmp(p.ssp,'bc03') and strcmp(p.imf,'chab') then begin
;      chabrier=1
;      endian='big'
;   endif
   endian='native'

   fitsed_generate_lut,p.lutfile, /verbose, $
                       version=p.version, $
                       ssp=p.ssp, $
                       imf=p.imf,$
                       modeldir=p.modeldir,$
                       chabrier=chabrier, endian=endian, $
                       bandpass_file=p.bandpass_file, $
                       catalog=p.catalog, $
                       zmin=p.zmin,$
                       zmax=p.zmax,$
                       zstep=p.zstep,$
                       zlog=p.zlog,$
                       zcustom=p.zcustom,$
                       extinction_law=p.extinction_law, $
                       ebvmin=p.ebvmin, $
                       ebvmax=p.ebvmax, $
                       ebvstep=p.ebvstep, $
                       log_agemin=p.log_agemin, $
                       log_agemax=p.log_agemax, $
                       log_agestep=p.log_agestep, $
                       sfh=p.sfh, $
                       alpha=p.alpha, $
                       beta=p.beta, $
                       nebular_fesc=p.nebular_fesc, $
                       nebular_continuum=0, $
                       nebular_nolya=1, $
                       log_ageArr=log_age, ebvArr=ebv, deltaArr=delta, $
                       tauArr=p.tau, metalArr=p.metal
endif 


;; fit photometry against LUT: 
if ~p.skipfit then $
   fitsed_execfit, p, data=data
;;
;; write out results to output file
if ~p.skipfit then $
   ;; always write out catalog, first full catalog, then bestfit catalog
   fitsed_write_catalog,p,  data=data, paramfile=paramfile $
                         ; p.fitsed_cat, p.outdir $
else $
   fitsed_write_catalog, p, data=data, $
                         paramfile=paramfile
                         ; catalog=p.catalog, paramfile=paramfile, $
                         ;ab_zeropoint=p.ab_zeropoint, $
                         ;name_zphot=p.name_zphot, $
                         ;name_zspec=p.name_zspec, $
                         ;p.fitsed_cat, p.outdir

message,/cont, 'All done!'

end
