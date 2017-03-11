
paramfile='fitsed.param'

fitsed_read_param, paramfile, params=p

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

end
