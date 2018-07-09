
; updated for BC03 2016 update... 

PRO READCOLORFILE, color, isedroot, metallicity=metallicity, $
                   colorext=colorext, verbose=verbose

if keyword_set(colorext) then colorext=colorext else colorext='.1color'

;; read in the associated bc color file and return a structure.

openr,lun,/get_lun, isedroot+colorext
done=0
count=0
line=''
while ~done do begin
  readf,lun,line
  count++
  ;;print,line
  if strcmp(line,'#',1) ne 1 then done=1
  if ~done then begin
     if count gt 1 then begin
        if strcmp(lastline,'#') eq 0 and $
        strcmp(lastline,'#                                                                                                                      M*_tot=') eq 0 and $
        strcmp(lastline,'#             <---------------- SDSS AB mag --------------->       <---------------------- CFHT MegaCam AB mag --------------------->  '+$
        '    <---------------- GALEX AB mag and flux --------------->') eq 0 then penlastline=lastline 
     endif
     lastline=line
  endif
endwhile
close,lun
free_lun,lun

tcolor = read_ascii(isedroot+colorext,data_start=count-1,header=header)
; guess at the tag name for the structure:
; (THIS PRODUCES AN ANNOYING TYPE-CONVERSION ERROR....)
;tag = tag_names(tcolor)
tcolor = tcolor.(0)

; establish array indices for various quantities:
; the crappy thing is the 2003 and 2016 lines are swapped.  WTF BC!!
;
; test them to figure out which one has "(1)" as the second element in
; the array...hopefully that'll do it:
check = (strsplit(lastline,' #',/extract))[0]
if strcmp(check,'(1)') eq 1 then bc2003=1 else bc2003=0

if ~bc2003 then begin
   tags = strsplit(lastline,' #',/extract)
   indices = strsplit(penlastline,' #()',/extract)
endif else begin
   tags = strsplit(penlastline,' #',/extract)
   indices = strsplit(lastline,' #()',/extract)
endelse

color=create_struct('file', isedroot+colorext, 'header', header)

metalline=strsplit(header[5],' =',/ext)
metallicity = float( metalline[n_elements(metalline)-2])

for aa=0,n_elements(tags)-1 do begin
;;    print,aa,n_elements(tags)-1,tags[aa]
    case tags[aa] of
        "log-age-yr" : color=create_struct(color, 'logage', tcolor[long(indices[aa])-1, *])
        "Mbol" : color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "Umag" : color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "Bmag" : color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "Rmag" : color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "J2mag": color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "Vmag" : color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "Kmag" : color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        '14-V' : color=create_struct(color, 'm14mv', tcolor[long(indices[aa])-1,*])
        '17-V' : color=create_struct(color, 'm17mv', tcolor[long(indices[aa])-1,*])
        '22-V' : color=create_struct(color, 'm22mv', tcolor[long(indices[aa])-1,*])
        '27-V' : color=create_struct(color, 'm27mv', tcolor[long(indices[aa])-1,*])
        "U-J"  : color=create_struct(color, 'umj', tcolor[long(indices[aa])-1,*])
        "J-F"  : color=create_struct(color, 'jmf', tcolor[long(indices[aa])-1,*])
        "F-N"  : color=create_struct(color, 'fmn', tcolor[long(indices[aa])-1,*])
        "U-B"  : color=create_struct(color, 'umb', tcolor[long(indices[aa])-1,*])
        "B-V"  : color=create_struct(color, 'bmv', tcolor[long(indices[aa])-1,*])
        "V-R"  : color=create_struct(color, 'vmr', tcolor[long(indices[aa])-1,*])
        "V-I"  : color=create_struct(color, 'vmi', tcolor[long(indices[aa])-1,*])
        "V-J"  : color=create_struct(color, 'vmj', tcolor[long(indices[aa])-1,*])
        "V-K"  : color=create_struct(color, 'vmk', tcolor[long(indices[aa])-1,*])
        "R-I"  : color=create_struct(color, 'rmi', tcolor[long(indices[aa])-1,*])
        "J-H"  : color=create_struct(color, 'jmh', tcolor[long(indices[aa])-1,*])
        "H-K"  : color=create_struct(color, 'hmk', tcolor[long(indices[aa])-1,*])
        "V-K'" : color=create_struct(color, 'vmkp', tcolor[long(indices[aa])-1,*])
        "(J-H)2M" : color=create_struct(color, 'jmh2m', tcolor[long(indices[aa])-1,*])
        "(J-Ks)2M" : color=create_struct(color, 'jmks2m', tcolor[long(indices[aa])-1,*])
        "B(4000)":  color=create_struct(color, 'b4000', tcolor[long(indices[aa])-1,*])
        "B4_VN":  color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "B4_SDSS":  color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "B(912)":  color=create_struct(color, 'b912', tcolor[long(indices[aa])-1,*])
        "NLy":  color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "SNR/yr/Lo":  color=create_struct(color, 'snr_yr_lo', tcolor[long(indices[aa])-1,*])
        "PNBR/yr/Lo":  color=create_struct(color, 'pnbr_yr_lo', tcolor[long(indices[aa])-1,*])
        "N(BH)":  color=create_struct(color, 'n_bh', tcolor[long(indices[aa])-1,*])
        "N(NS)":  color=create_struct(color, 'n_ns', tcolor[long(indices[aa])-1,*])
        "N(WD)":  color=create_struct(color, 'n_wd', tcolor[long(indices[aa])-1,*])
        "M(Remnants)":  color=create_struct(color, 'm_remnants', tcolor[long(indices[aa])-1,*])
        "M_remnants":  color=create_struct(color, 'm_remnants', tcolor[long(indices[aa])-1,*])
        "M*/Lb":  color=create_struct(color, 'mstar_Lb', tcolor[long(indices[aa])-1,*])
        "M*/Lv":  color=create_struct(color, 'mstar_Lv', tcolor[long(indices[aa])-1,*])
        "M*_liv/Lb":  color=create_struct(color, 'mstar_Lb', tcolor[long(indices[aa])-1,*])
        "M*_liv/Lv":  color=create_struct(color, 'mstar_Lv', tcolor[long(indices[aa])-1,*])
        "M*_liv/Lk":  color=create_struct(color, 'mstar_Lk', tcolor[long(indices[aa])-1,*])
        "M*":  color=create_struct(color, 'mstar', tcolor[long(indices[aa])-1,*])
        "M*_liv":  color=create_struct(color, 'mstar', tcolor[long(indices[aa])-1,*])
        "Mgas":  color=create_struct(color, 'mgas', tcolor[long(indices[aa])-1,*])
        "M_ret_gas":  color=create_struct(color, 'mgas', tcolor[long(indices[aa])-1,*])
        "Mgalaxy":  color=create_struct(color, 'mgalaxy', tcolor[long(indices[aa])-1,*])
        "M_galaxy":  color=create_struct(color, 'mgalaxy', tcolor[long(indices[aa])-1,*])
        "SFR/yr":  color=create_struct(color, 'SFR', tcolor[long(indices[aa])-1,*])
        "b(t)*'s/yr":  color=create_struct(color, 'bstar_yr', tcolor[long(indices[aa])-1,*])
        "B(t)/yr/Lo":  color=create_struct(color, 'b_yr_lo', tcolor[long(indices[aa])-1,*])
        "Turnoff_mass":  color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "BPMS/BMS":  color=create_struct(color, 'bpms_bms', tcolor[long(indices[aa])-1,*])
        "g_AB": color=create_struct(color, tags[aa], tcolor[long(indices[aa])-1,*])
        "(u-g)AB": color=create_struct(color, 'umg_ab', tcolor[long(indices[aa])-1,*])
        "(g-r)AB": color=create_struct(color, 'gmr_ab', tcolor[long(indices[aa])-1,*])
        "(g-i)AB": color=create_struct(color, 'gmi_ab', tcolor[long(indices[aa])-1,*])
        "(g-z)AB": color=create_struct(color, 'gmz_ab', tcolor[long(indices[aa])-1,*])
        else :
    endcase
 endfor

end


PRO FITSED_READBCT, ISEDFILE, age=age, spec=spec, na=na,$
             color4=color4, color5=color5, colorab=colorab,$
             color1=color1, color2=color2, color3=color3, $
             metallicity=metallicity, $
             endian=endian,$
             ignorecolorfiles=ignorecolorfiles, $
             _EXTRA=_EXTRA

;!!!! IMPORTANT !!!!
; The endian keyword controls the byte ordering, see help for
; read_binary.  Here are the important notes: ENDIAN Set this keyword
; to one of three string values: "big", "little" or "native" which
; specifies the byte ordering of the file to be read. If the computer
; running READ_BINARY uses byte ordering that is different than that
; of the file, READ_BINARY will swap the order of bytes in multi-byte
; data types read from the file. (Default: "native" = perform no byte
; swapping.) See Files and I/O for more information on byte order and
; cross-platform programming. ;
; 
; Some of this might be mitigated by setting "/swap_if_big_endian".
; 
; _EXTRA is fed below

if not keyword_set(endian) then endian="native"

; This procedure reads in a Bruzual & Charlot file from a binary .ised
; file and returns two arrays, one for age, and one for the spectra
;
; Input:
; ISEDFILE = Bruzual & Charlot .ised file (binary raw file for BC03
; release) *** AND ONLY FOR THE BC03 RELEASE?! ***
; na = number of time steps, (for most models it is 221)
;
; Returns:
; AGE array = fltarr(na) lists the time steps in the BC models (in yrs)
; SPEC array = fltarr(na+1,1221)
;            index[0,*] stores the wavelength information (angstroms)
;            index[1:na,*] stores the spectra that correspond to the
;                           na time steps.
; COLOR* arrays are the *.*color files for the particular ISED
;
; EXAMPLE:
; Read in the age, spectra, and color files for a ised file called,
; 'krap.ised':
; 
; IDL> readbct,'krap.ised',age=age,spec=spec,na=na,color1=c1,color2=c2,color3=c3, $
; IDL> color4=c4,color5=c5,colorab=cab
;
; IDL> help,c1,c2,c3,c4,c5,cab,age,spec
; C1              STRUCT    = -> <Anonymous> Array[1]
; C2              STRUCT    = -> <Anonymous> Array[1]
; C3              STRUCT    = -> <Anonymous> Array[1]
; C4              STRUCT    = -> <Anonymous> Array[1]
; C5              STRUCT    = -> <Anonymous> Array[1]
; CAB             STRUCT    = -> <Anonymous> Array[1]
; AGE             FLOAT     = Array[221]
; SPEC            FLOAT     = Array[222, 1221]
;
; See the Bruzual & Charlot documentation for a description of the
; tags in the color* structures.
;
; Calls:
; READCOLORFILE (see below)
;
; Written by Casey Papovich (some time ago...)
; Modified by A Stutz, May 2004.
; Modified by C Papovich, June 2004 to read the *color files
; Modified by C Papovich, June 2005 to read the Charbrier IMFs
; Modified more by CJP, Oct 2009 to do cb07 models Chabrier and
; Salpeter...
; Modified more by CJP, Jun 2017 for 2016 files: 
; Heavily modified/chaged by CJP to try and guess format of .ised and
; color files.  This should work... send bugs to CJP (papovich@tamu.edu).

if (not file_test(isedfile,/read)) then begin
    print,'% File ',isedfile,' does not exist or is unreadable'
endif

if keyword_set(verbose) and n_elements(na) eq 0 then  print, 'Using default value for number of time steps: 221'
if (n_elements(na) eq 0) then na = 221

;; make update to figure out number of time steps, wavelength steps,
;; etc:
; scan to find first real line
addonestep=0

REPEAT1 :  

RightEndian=0
count=0

while ~RightEndian do begin
   openr,lun,/get_lun, isedfile

   done=0
   while ~done do begin
      npre = read_binary(lun,data_type=2,data_dims=1,endian=endian)
      done=(npre ne 0)

   endwhile
   ;print,endian, npre
   if npre gt 1000 and npre lt 2000 then break
   ;; endian needs to be swapped - no way around it:

   close,lun
   free_lun, lun
   
   ;; cycle through them, but this really depends on what machine made
   ;; and is using the models. If you're using the same machine
   ;; to make the .ised files and to read them in, then it should be "native".
   if strcmp(endian,'native') then endian='big' $
   else $
      if strcmp(endian,'big') then endian='little' $
      else if strcmp(endian,'little') then endian='native'
   count++
   if count gt 3 then begin
      print,'% Crap, cannot figure out endian.  Died with endian='+endian
      stop
   endif
endwhile
nmark=npre
npre = read_binary(lun,data_type=2,data_dims=2,endian=endian)
nage = npre[1]

;; THIS IS IMPORTANT __ THERE IS SOMETHING DIFFERENT PHYSICALLY ABOUT
;; THE STRUCTURE FROM THE 2003 TO THE 2016 MODELS.  THIS CHECKS FOR
;; THAT AND TRIES TO CORRECT IT. 
if addonestep gt 0 then begin
   for i=1,addonestep do begin
      chuck=read_binary(lun,data_type=2,data_dims=1,endian=endian)
      ;print,addonestep, chuck
   endfor
endif

age = read_binary(lun,data_type=4,data_dims=nage,endian=endian)
;print, age[0], age[1], addonestep
if age[0] ne 0 or finite(alog10(age[1])) eq 0 or  abs(alog10(age[1])-5) gt 1.0 then begin
   addonestep+=1
   if addonestep gt 5 then begin
      print,'% Crap.  Cannot find age array, crashing'
      stop
   endif
   close,lun
   free_lun,lun
   
   GOTO,REPEAT1
endif

; scan foward until you find "1184" nmark=1184 for the 2016 models:
done=0
while ~done do begin
   res=read_binary(lun,data_type=3,data_dims=1,endian=endian)
   ;;done = (read_binary(lun,data_type=3,data_dims=1)==npre[0])
   ;;sprint,res,nmark
   done = (res eq nmark)
   ;print, res, npre[0], done
endwhile

; now get wavelength steps:
nstep = read_binary(lun,data_type= 3,data_dims=2,endian=endian)
nwave=nstep[1]

; now set up arrays to read fluxes
spec = fltarr(nage+1, nwave) ;; 0th is wavelength
;; read in wavelengths:
tflux = read_binary(lun,data_type=4,data_dims=nwave,endian=endian)
spec[0,*] = tflux

;print,'Reading in ised file with '+strn(nage)+' age steps'
;print,'                     and '+strn(nwave)+' wavelength steps'

; now iterate to get fluxes: 
for i=1,nage do begin
   done=0
   while ~done do begin
      done = (read_binary(lun,data_type=3,data_dims=1,endian=endian) eq nwave)
   endwhile
   tflux = read_binary(lun,data_type=4,data_dims=nwave,endian=endian)
   spec[i,*] = tflux
   ;; code has nind, but ignore for now:
   nind = read_binary(lun,data_type=3,data_dims=1,endian=endian)
   indexes = read_binary(lun,data_type=4,data_dims=nind,endian=endian)
endfor

close,lun
free_lun,lun

;; read in the color files if requested
isedroot = (strsplit(isedfile,".ised",/extract,/regex))[0]

if keyword_set(verbose) then $
  print,"%%%% IGNORE 'Type conversion error' errors (there will be six of them)"
if not keyword_set(ignorecolorfiles) then begin
   readcolorfile,color1,isedroot,colorext='.1color', $
                 metallicity=metallicity
   readcolorfile,color2,isedroot,colorext='.2color'
   readcolorfile,color3,isedroot,colorext='.3color'
   readcolorfile,color4,isedroot,colorext='.4color'
   readcolorfile,color5,isedroot,colorext='.5color'
   readcolorfile,colorab,isedroot,colorext='.1ABmag'
endif
if keyword_set(verbose) then $
  print,"%%%% IGNORE 'Type conversion error' errors (there will be six of them)"

RETURN
END

