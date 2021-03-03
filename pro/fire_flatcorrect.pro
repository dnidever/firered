function fire_flatcorrect,im,flat

  ;; Correct for variations in QA
  newim = im
  
  ;; Only correct "good" pixels
  medflat = median(flat.flux)
  gdpix = where(flat.flux gt 0.001 and im.mask eq 1,ngdpix)
  newim.flux[gdpix] /= flat.flux[gdpix]/medflat       ; divide by flat field image

  ; Add processing information to header
  ;--------------------------------------
  ;  Current timestamp information
  ;  Sun Oct  7 15:38:23 2012
  date = systime(0)
  datearr = strtrim(strsplit(date,' ',/extract),2)
  timarr = strsplit(datearr[3],':',/extract)
  datestr = datearr[1]+' '+datearr[2]+' '+strjoin(timarr[0:1],':')
  head = im.head
  fxaddpar,head,'FLATCOR',datestr+' Flat is '+flat.file+', scale '+strtrim(string(medflat,format='(F20.2)'),2)
  newim = fire_updateheader(newim,head)
  
  return,newim

  end
