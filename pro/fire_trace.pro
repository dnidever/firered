pro fire_trace,objfile,bndfile,bpmfile

  ;; Trace orders of a stellar spectrum

  objfile = 'ut131222/fire_0084.fits'
  bndfile = 'fire_boundary_0011.fits'
  bpmfile = 'bpm2.fits'

  im = fire_makeimage(objfile)
  bpm = fire_makeimage(bpmfile)
  bstr = mrdfits(bndfile,1)
  norders = n_elements(bstr)
  npix = 2048

  ;; Correct for bpm
  im = fire_bpmcorrect(im,bpm)

  ;; Loop over orders
  norders = n_elements(bstr)
  tstr = replicate({order:0,bndx0:0,bndx1:0,bndy0coef:fltarr(4),bndy1coef:fltarr(4),$
                    tycoef:fltarr(4),tsigcoef:fltarr(4),thwhmcoef:fltarr(4),$
                    tmoffcoef:fltarr(4)},norders)
  ;; PUT IN BOUNDARY INFORMATION AS WELL
  for i=0,norders-1 do begin
    print,'Tracing order ',strtrim(i+1,2)
    FIRE_TRACE_ORDER,bstr[i],im,tstr1,coeffs
    ;; save it in a large structure
    tstr[i].order = i+1
    tstr[i].bndx0 = bstr[i].lo
    tstr[i].bndx1 = bstr[i].hi
    tstr[i].bndy0coef = bstr[i].y0coef
    tstr[i].bndy1coef = bstr[i].y1coef    
    ;tstr[i].data = ptr_new(tstr1)
    tstr[i].tycoef = coeffs.tcoef
    tstr[i].tsigcoef = coeffs.sigcoef
    tstr[i].thwhmcoef = coeffs.hwhmcoef
    tstr[i].tmoffcoef = coeffs.moffcoef
  endfor

  ;;mwrfits,tstr,'fire_trace_0084.fits',/create
  
  stop

  end
