
; example usage:
; 
;        y_age = fitsed_yp6(10d^log_ageArr, pdf, 0, a=ebvArr, b=metalArr, c=tauArr, d=alphaArr, e=betaArr,$
;                      norm=norm)
; ;apply the normalization:
;         pdf /= norm
;         
;         y_ebv = fitsed_yp6(ebvArr, pdf, 1, a=10d^log_ageArr, b=metalArr, c=tauArr, d=alphaArr, e=betaArr)
;         y_metal = fitsed_yp6(metalArr, pdf, 2, a=10d^log_ageArr, b=ebvArr, c=tauArr, d=alphaArr, e=betaArr)
;         if n_tau gt 1 then $
;            y_tau = fitsed_yp6(tauArr, pdf, 3, a=10d^log_ageArr, b=ebvArr, c=metalArr, d=AlphaArr, e=betaArr) $
;         else   y_tau = [1.0]
;         if n_alpha gt 1 then $
;            y_alpha = fitsed_yp6(alphaArr, pdf, 4, a=10d^log_ageArr, b=ebvArr, c=metalArr, d=tauArr, e=betaArr) $
;         else y_alpha = [1.0]
;         if n_beta gt 1 then $
;            y_beta = fitsed_yp6(betaArr, pdf, 5, a=10d^log_ageArr, b=ebvArr, c=metalArr, d=tauArr, e=alphaArr) $
;         else y_beta = [1.0]
; 

function fitsed_yp6, x, pdf, i_x, a=a, b=b, c=c, d=d, e=e, norm=norm
;; Same as y_p, but allows 5 variables (age, ebv, metal, tau, alpha,
;; beta) this function takes the PDF, and integrates over a, b, c, and d to
;; get y(x), the posterior for x (marginalized over the other
;; parameters).
;;
;; a = array for the first variable, 
;; b= array for the second variable
;; c, d,e = third, fourth fifth variable
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
  indexes = bindgen(6)
  indexes=indexes[where(i_x ne indexes)]

  y_x = fltarr(szx)
  tmpA = fltarr(szx, $
                n_elements(b),$
                n_elements(c),$
                n_elements(d),$
                n_elements(e))
  tmpB = fltarr(szx, $
                n_elements(c),$
                n_elements(d),$
                n_elements(e))
  tmpC = fltarr(szx, $
                n_elements(d),$
                n_elements(e))
  tmpD = fltarr(szx, $
                n_elements(e))

  for i=0,n_elements(y_x)-1 do begin
     for m=0,n_elements(e)-1 do begin
        for l=0,n_elements(d)-1 do begin
           for k=0,n_elements(c)-1 do begin
              for j=0,n_elements(b)-1 do begin
                 
                 case i_x of  
                    0 : tmpA[i,j,k,l,m] = mytsum(a,pdf[i,*,j,k,l,m])
                    1 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,i,j,k,l,m])
                    2 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,j,i,k,l,m])
                    3 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,j,k,i,l,m])
                    4 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,j,k,l,i,m])
                    5 : tmpA[i,j,k,l,m] = mytsum(a,pdf[*,j,k,l,m,i])
                    else : begin
                       message, 'i_x must be 0...5, cannot be '+strn(i_x)
                       stop
                    end
                 endcase
              endfor
              tmpB[i,k,l,m] = mytsum(b,tmpA[i,*,k,l,m])
           endfor
           tmpC[i,l,m] = mytsum(c,tmpB[i,*,l,m])
        endfor
        tmpD[i,m] = mytsum(d,tmpC[i,*,m])
     endfor
     y_x[i] = mytsum(e,tmpD[i,*])
  endfor
  
  norm = mytsum(x,y_x)
  y_x /= norm

  return, y_x

END


