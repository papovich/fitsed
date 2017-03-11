
;; compile some useful functions: 


function findel, x, arr
  
  return, (where( abs( x - arr) eq min(abs(x - arr))))[0]
  
end


;; main exec file.  Calls everything else, given a paramfile

pro fitsed, paramfile

fitsed_read_param, paramfile, params=p



if ~file_test(p.lutfile) then begin
   ;; generate Look Up Table (LUT): 

   delvarx, chabrier, endian
   if strcmp(p.ssp,'bc03') and strcmp(p.imf,'chab') then begin
      chabrier=1
      endian='big'
   endif

   fitsed_generate_lut,p.lutfile, /verbose, $
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
                       extinction_law=p.extinction_law, $
                       ebvmin=p.ebvmin, $
                       ebvmax=p.ebvmax, $
                       ebvstep=p.ebvstep, $
                       log_agemin=p.log_agemin, $
                       log_agemax=p.log_agemax, $
                       log_agestep=p.log_agestep, $
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
   fitsed_write_catalog, data=data, paramfile=paramfile, $
                         bestfit=p.bestfit, $
                         p.fitsed_cat, p.outdir $
else $
   fitsed_write_catalog, data=data, catalog=p.catalog, paramfile=paramfile, $
                         bestfit=p.bestfit, $
                         ab_zeropoint=p.ab_zeropoint, $
                         name_zphot=p.name_zphot, $
                         name_zspec=p.name_zspec, $
                         p.fitsed_cat, p.outdir

message,/cont, 'All done!'

end
