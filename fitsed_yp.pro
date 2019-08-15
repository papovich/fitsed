; EXAMPLE usage: 

;         y_age = fitsed_yp(10d^log_ageArr, pdf, 0, a=ebvArr, b=metalArr, c=tauArr, $
;                     norm=norm)
; ;apply the normalization:
;         pdf /= norm
; 
;         y_ebv = fitsed_yp(ebvArr, pdf, 1, a=10d^log_ageArr, b=metalArr, c=tauArr)
;         y_metal = fitsed_yp(metalArr, pdf, 2, a=10d^log_ageArr, b=ebvArr, c=tauArr)
;         if n_tau gt 1 then begin
;            y_tau = fitsed_yp(tauArr, pdf, 3, a=10d^log_ageArr, b=ebvArr, c=metalArr)
;         endif else begin
;            y_tau = [1.0]
;         endelse
; 

function fitsed_yp, x, pdf, i_x, a=a, b=b, c=c, norm=norm
;; this function takes the PDF, and integrates over a, b, c, and d to
;; get y(x), the posterior for x (marginalized over the other
;; parameters).
;;
;; a = array for the first variable, 
;; b= array for the second variable
;; c= third
;; i_x is assumed to be dimension of the variable for x.
;; 
;; Example:
;; if x is the 0th dimension, then i_x = 0 and pdf would have the form:
;; pdf[*, n_elements(a), n_elements(b), n_elements(c), n_elements(d)]
;;
;; if x is the 1st dimension, then i_x = 1 and the pdf would have the form:
;; pdf[n_elements(a), *, n_elements(b), n_elements(c), n_elements(d)]
;; 
;; etc... 

  sz = size(pdf,/dim)
  szx = sz[i_x]
  indexes = bindgen(4)
  indexes=indexes[where(i_x ne indexes)]

  y_x = fltarr(szx)
  tmpA = fltarr(szx, $
                n_elements(b),$
                n_elements(c))
  tmpB = fltarr(szx, $
                n_elements(c))
;  tmpC = fltarr(szx, $
;                n_elements(d))

  for i=0,n_elements(y_x)-1 do begin
;     for l=0,n_elements(d)-1 do begin
     for k=0,n_elements(c)-1 do begin
        for j=0,n_elements(b)-1 do begin
           case i_x of  
              0 : tmpA[i,j,k] = mytsum(a,pdf[i,*,j,k])
              1 : tmpA[i,j,k] = mytsum(a,pdf[*,i,j,k])
              2 : tmpA[i,j,k] = mytsum(a,pdf[*,j,i,k])
              3 : tmpA[i,j,k] = mytsum(a,pdf[*,j,k,i])
;              4 : tmpA[i,j,k] = mytsum(a,pdf[*,j,k,i])
              else : begin
                 message, 'i_x must be 0...3, cannot be '+strn(i_x)
                 stop
              end
           endcase
        endfor
        tmpB[i,k] = mytsum(b,tmpA[i,*,k])
     endfor
     y_x[i] = mytsum(c,tmpB[i,*])
  endfor

  norm = mytsum(x,y_x)
  y_x /= norm

  return, y_x

END
