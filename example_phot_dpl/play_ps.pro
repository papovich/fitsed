

xsize=7.5
ysize=6.5
charsize=1.8


set_plot,'ps'
device,filename='play.ps',encap=0
device,xsize=xsize,ysize=ysize,yoff=0.5, xoff=0.5,/inches
device,/isolatin
device,/color
   
!p.thick=4 & !p.charthick=4 & !x.thick=4 & !y.thick=4
!p.font=-1 & !p.charsize=charsize


for i=1,1067 do begin
    fitsed_plot_bestspec,i,par='fitsed.param', /label
;    read,a,prompt='Enter any key (and return): '
 endfor
device,/close
!p.thick=1 & !p.charthick=1 & !x.thick=1 & !y.thick=1
!p.font=-1
set_plot,'x'
end
   
