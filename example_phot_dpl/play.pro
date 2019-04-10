
a=''
for i=1,1067 do begin
    fitsed_plot_bestspec,i,par='fitsed.param'
    read,a,prompt='Enter any key (and return): '
 endfor

end
   
