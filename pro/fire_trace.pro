pro fire_trace,bstr,im,err,tstr

  ;; Trace orders of a stellar spectrum

  ;; Loop over orders
  norders = n_elements(bstr)
  tstr = replicate({order:0,bndx0:0,bndx1:0,bndy0coef:fltarr(4),bndy1coef:fltarr(4),$
                    tycoef:fltarr(4),tsigcoef:fltarr(4)},norders)
  ;; PUT IN BOUNDARY INFORMATION AS WELL
  for i=0,norders-1 do begin
    FIRE_TRACE_ORDER,bstr[i],im,err,tstr1,tcoef,sigcoef
    ;; save it in a large structure
    tstr[i].order = i+1
    tstr[i].bndx0 = bstr[i].lo
    tstr[i].bndx1 = bstr[i].hi
    tstr[i].bndy0coef = bstr[i].y0coef
    tstr[i].bndy1coef = bstr[i].y1coef    
    ;tstr[i].data = ptr_new(tstr1)
    tstr[i].tycoef = tcoef
    tstr[i].tsigcoef = sigcoef    
  endfor

  ;stop

  end
