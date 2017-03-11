pro FITSED_CHECKTIME, timestart, comment

  if total(size(comment)) eq 0 then comment=0

  temptime= systime(/seconds) - timestart ; get current run time
  days= LONG(temptime/60./60./24.)
  hours= LONG(((temptime/60./60./24.)-days)*60.)
  minutes= LONG( ((((temptime/60./60./24.)-days)*60.)-hours)*60 )
  seconds= (((((temptime/60./60./24.)-days)*60.)-hours)*60 - minutes)*60
  days= STRCOMPRESS(days)
  hours= STRCOMPRESS(hours)
  minutes=STRCOMPRESS(minutes)
  seconds=STRCOMPRESS(seconds)
  
  if comment then print, '%  Working on ',comment
  print, '% Runtime= '+days+'d'+hours+'h'+minutes+'m'+seconds+'s'

end
