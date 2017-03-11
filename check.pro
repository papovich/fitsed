
; read in fitsed results and fast results and plot mass v mass

pro read_fitsed_cat, fitsedcat, best=best, cat=cat


if keyword_set(best) then begin
   readcol,fitsedcat, $
           format='l,f, f, f, f, f, f, f', $
           id, z, tau, $
           metal, $
           lg_age, $
           ebv, $
           delta, $
           lg_mass
   cat = makearrstruct({ id: id, z: z, $
                         tau: tau, $
                         metal: metal, $
                         lg_age: lg_age, $
                         ebv: ebv, $
                         delta: delta, $
                         lg_mass: lg_mass, $
                         lg_sfr: lg_mass*0.0})

endif else begin
   readcol,fitsedcat, $
           format='l,f, f,f,f, f,f,f, f,f,f, f,f,f, f,f,f, f,f,f, f,f,f', $
           id, z, tau, taul68, tauh68, $
           metal, metall68, metalh68, $
           lg_age, lg_agel68, lg_ageh68, $
           ebv, ebvl68, ebvh68, $
           delta, deltal68, deltah68, $ 
           lg_mass, lg_massl68, lg_massh68, $
           lg_sfr, lg_sfrl68, lg_sfrh68
   cat = makearrstruct({ id: id, z: z, $
                         tau: tau, taul68: taul68, tauh68: tauh68, $
                         metal: metal, metall68: metall68, metalh68: metalh68, $
                         lg_age: lg_age, lg_agel68: lg_agel68, lg_ageh68: lg_ageh68, $
                         ebv: ebv, ebvl68: ebvl68, ebvh68: ebvh68, $
                         delta: delta, deltal68: deltal68, deltah68: deltah68, $
                         lg_mass: lg_mass, lg_massl68: lg_massl68, lg_massh68: lg_massh68, $
                         lg_sfr: lg_sfr, lg_sfrl68: lg_sfrl68, lg_sfrh68: lg_sfrh68})
endelse

end


pro check, fitsedcat, fastcat,stop=stop, best=best

read_fitsed_cat, fitsedcat, best=best, cat=fsed

readcol, fastcat, $
         format='l,f, f,f,f, f,f,f, f,f,f', $
         id, z, ltau, metal, lage, av, lmass, lsfr, lssfr, la2t, chi2
fast = makearrstruct( { id: id, z: z, ltau: ltau, metal: metal, lage: lage, av: av, $
                        lmass: lmass, lsfr: lsfr, lssfr: lssfr, $
                        la2t: la2t, chi2: chi2})

plot, fsed.lg_mass, fast.lmass, xtit='Log Mass (FITSED)', ytit='Log Mass (FAST)', $
      psym=3, xr=[7,12], yr=[7,12],/yst
plots, !x.crange, !x.crange, line=1, color=cjp_icolor('red')

if keyword_set(stop) then stop

end
   
;fsfile=findfile('example_phot/hdfn*fitsed.cat')
fsfile=findfile('example_phot/hdfn*fitsed.bestfit.cat')
fastfile=findfile('example_phot/hdfn*.fout')
    
check, fsfile[0], fastfile[0], /stop, /best

;;     compare median to best fits: 
fsfile=findfile('example_phot/hdfn*fitsed.cat')
read_fitsed_cat, fsfile[0],cat=fsed 
fsfile=findfile('example_phot/hdfn*fitsed.bestfit.cat')
read_fitsed_cat, fsfile[0],cat=fbest,/best

plot, fsed.lg_mass, fbest.lg_mass, xtit='Log Mass (median)', ytit='Log Mass (best)',  psym=3, xr=[7,12], yr=[7,12],/yst
plots, !x.crange, !x.crange, line=1, color=cjp_icolor('red')


end   
