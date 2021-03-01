function fire_bpmcorrect,im,bpm

  ;; Correct for BPM

  bdpix = where(bpm.flux eq 1,nbdpix)
  if nbdpix gt 0 then begin
    temp = im.flux
    temp[bdpix] = !values.f_nan
    mn = smooth(temp,21,/edge_truncate,/nan)
    ;med = median(im.flux,5)
    im.flux[bdpix] = mn[bdpix]
    im.err[bdpix] = 1e30
    im.mask[bdpix] = 0
  endif

  ; Add processing information to header
  ;--------------------------------------
  ;  Current timestamp information
  ;  Sun Oct  7 15:38:23 2012
  date = systime(0)
  datearr = strtrim(strsplit(date,' ',/extract),2)
  timarr = strsplit(datearr[3],':',/extract)
  datestr = datearr[1]+' '+datearr[2]+' '+strjoin(timarr[0:1],':')
  fxaddpar,head,'BPM',datestr+' BPM is '+bpm.file+', '+strtrim(nbdpix,2)+' bad pixels'
  im = fire_updateheader(im,head)

  return,im

  end
