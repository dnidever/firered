pro fire_trace,bstr,im,err

  ;; Trace orders of a stellar spectrum

  ;; Loop over orders
  norders = n_elements(bstr)
  tstr = relicate({order:0,data:ptr_new(),tcoef:fltarr(4)},norders)
  for i=0,norders-1 do begin
    FIRE_TRACE_ORDER,bstr[i],im,err,tstr1,tcoef
    ;; save it in a large structure
    tstr[i].order = i+1
    tstr[i].data = ptr_new(tstr1)
    tstr[i].tcoef = tcoef
  endfor

  stop

  end
