function fire_darkcorrect,im,dark

  ;; Correct for Dark current
  newim = im
  
  ;; Scale by exptime
  exptime = sxpar(newim.head,'exptime')
  darkflux = dark.flux * exptime
  ;; Only correct "good" pixels
  gdpix = where(newim.mask eq 1,ngdpix)
  newim.flux[gdpix] -= darkflux[gdpix]

  
  ; Add processing information to header
  ;--------------------------------------
  ;  Current timestamp information
  ;  Sun Oct  7 15:38:23 2012
  date = systime(0)
  datearr = strtrim(strsplit(date,' ',/extract),2)
  timarr = strsplit(datearr[3],':',/extract)
  datestr = datearr[1]+' '+datearr[2]+' '+strjoin(timarr[0:1],':')
  head = im.head
  fxaddpar,head,'DARK',datestr+' DARK is '+dark.file
  newim = fire_updateheader(newim,head)
  
  return,newim

  end
