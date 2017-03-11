
PRO READCOLORFILE, color, isedroot, bcheaderlines, metallicity=metallicity, $
                   colorext=colorext, verbose=verbose

if keyword_set(colorext) then colorext=colorext else colorext='.1color'

;; read in the associated bc color file and return a structure.

tcolor = read_ascii(isedroot+colorext,data_start=bcheaderlines,header=header)

; guess at the tag name for the structure:
; (THIS PRODUCES AN ANNOYING TYPE-CONVERSION ERROR....)
;tag = tag_names(tcolor)
tcolor = tcolor.(0)

; establish array indices for various quantities:
tags = strsplit(header[bcheaderlines-2],' #',/extract)
indices = strsplit(header[bcheaderlines-1],' #()',/extract)
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
        "M*/Lb":  color=create_struct(color, 'mstar_Lb', tcolor[long(indices[aa])-1,*])
        "M*/Lv":  color=create_struct(color, 'mstar_Lv', tcolor[long(indices[aa])-1,*])
        "M*":  color=create_struct(color, 'mstar', tcolor[long(indices[aa])-1,*])
        "Mgas":  color=create_struct(color, 'mgas', tcolor[long(indices[aa])-1,*])
        "Mgalaxy":  color=create_struct(color, 'mgalaxy', tcolor[long(indices[aa])-1,*])
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


PRO FITSED_READBCT, ISEDFILE, age=age, spec=spec, na=na, high_res=high_res, $
                    color4=color4, color5=color5, colorab=colorab,$
                    color1=color1, color2=color2, color3=color3, $
                    metallicity=metallicity, $
                    bcheaderlines=bcheaderlines,chabrier=chabrier,$
                    endian=endian,cb07=cb07,ignorecolorfiles=ignorecolorfiles, $
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

if keyword_set(bcheaderlines) then bcheaderlines=bcheaderlines $
else bcheaderlines=29

if (not file_test(isedfile,/read)) then begin
    print,'% File ',isedfile,' does not exist or is unreadable'
endif

if keyword_set(verbose) and n_elements(na) eq 0 then  print, 'Using default value for number of time steps: 221'
if (n_elements(na) eq 0) then na = 221

tempbin = read_binary(isedfile,data_type=4,endian=endian)

if keyword_set(high_res) then begin
   if keyword_set(cb07) then begin
      print,'% ERROR: I have not hardcoded the number of wavelength elements for the HR CB07 models yet'
      message,'% Sorry'
   endif
; in high-resolution mode:
    ns=long(6900)
    offa = 2        ; offset from beginning of file to first age index
    offb = 84            ; spaces between age info and wavelength info
    if keyword_set(chabrier) then offb=84-6
    offc = 4 ; spaces between wavlength info and the first spectrum
    offs = long(56)             ; space between spectra

endif else begin
    ns=long(1221)
    offa = 2        ; offset from beginning of file to first age index
    offb = 84            ; spaces between age info and wavelength info
    offc = 1281                 ; spaces between wavlength info and the first spectrum
    offs = long(56)             ; space between spectra
    if keyword_set(chabrier) then offb=84-6

    if keyword_set(cb07) then begin 
       ns=long(1238)
       offc=1281+17
       if not keyword_set(chabrier) then begin
          offb=84
       endif
    endif
    
endelse
 
index_first_spec = na + offa + offb - 1 
ifs = index_first_spec
if (tempbin(ifs) ne 91.0) then begin
   print,'% tempbin(ifs) != 91.0, it is ='+strtrim(string(tempbin(ifs),2))+' instead!'
   print,'% MAYBE ITS A PROBLEM WITH THE ENDIAN ?'
   message, 'There is a problem with extracting the wavelength array'
endif

index_second_spec = ifs + ns - 1 + offc
ifs2 = index_second_spec
age = fltarr(na)
age[*] = tempbin[offa:offa+na-1]

spec = fltarr(na+1,ns)
spec[0,*] = tempbin[ifs:ifs+ns-1]

for i=1L,na do begin
    i1 = long(ifs2)+long((i-1L)*(ns+offs))
    i2 = long(ifs2)+long((i-1L)*(ns+offs)+ns-1L)
    spec[i,*] = tempbin[i1:i2]
    ;; this is a debugging 'stop' set i to the spectrum number
    ;; (normally 1-221) you want to view).
    if i eq -200 then stop
endfor

;; read in the color files if requested
isedroot = (strsplit(isedfile,".ised",/extract,/regex))[0]

if keyword_set(verbose) then $
  print,"%%%% IGNORE 'Type conversion error' errors (there will be six of them)"
if not keyword_set(ignorecolorfiles) then begin
   readcolorfile,color1,isedroot,bcheaderlines,colorext='.1color', $
                 metallicity=metallicity
   readcolorfile,color2,isedroot,bcheaderlines,colorext='.2color'
   readcolorfile,color3,isedroot,bcheaderlines,colorext='.3color'
   readcolorfile,color4,isedroot,bcheaderlines,colorext='.4color'
   if not keyword_set(chabrier) then $
      readcolorfile,color5,isedroot,bcheaderlines,colorext='.5color'
   readcolorfile,colorab,isedroot,bcheaderlines,colorext='.1ABmag'
endif
if keyword_set(verbose) then $
  print,"%%%% IGNORE 'Type conversion error' errors (there will be six of them)"

    

RETURN
END

