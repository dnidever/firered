pro fire_extract,objfile,arcfile,tracefile,bpmfile

  ;objfile = 'ut131222/fire_0046.fits'
  objfile = 'ut131222/fire_0084.fits'
  arcfile = 'ut131222/fire_0047.fits'
  boundaryfile = 'fire_boundary_0011.fits'
  tracefile = 'fire_trace_0084.fits'  
  bpmfile = 'bpm2.fits'

  fits_read,objfile,imobj,objhead
  fits_read,arcfile,imarc,archead
  fits_read,bpmfile,imbpm,bpmhead
  bstr = mrdfits(boundaryfile,1)
  tstr = mrdfits(tracefile,1)  
  norders = n_elements(bstr)
  npix = 2048

  ;; Correct for bpm
  bd = where(imbpm eq 1,nbd)
  medobj = median(imobj,5)
  imobj[bd] = medobj[bd]
  medarc = median(imarc,5)  
  imarc[bd] = imarc[bd]

  ;; Fix "flipped" pixels
  bd = where(imobj lt -1000,nbd)
  if nbd gt 0 then imobj[bd] += 65536
  bd = where(imarc lt -1000,nbd)
  if nbd gt 0 then imarc[bd] += 65536  

  ;; Correct for Dark current
  ;; scale by exptime
  
  ;; Do flat field correction

  ;; Construct error arrays
  gain = sxpar(objhead,'gain')  ; electrons/DU
  errobj = sqrt((imobj>1)/gain)
  errarc = sqrt((imarc>1)/gain)  

  stop
  ;fire_trace,bstr,imobj,errobj,tstr 
  
  ;; Order loop
  x = findgen(npix)
  outobj = replicate({order:0,data:ptr_new()},norders)
  outobj.order = lindgen(norders)+1
  outarc = replicate({order:0,data:ptr_new()},norders)  
  outarc.order = lindgen(norders)+1
  ;For i=1,norders-1 do begin
  For i=20,20 do begin     
    print,'order = ',i
    outobj1 = FIRE_EXTRACT_ORDER(bstr[i],imobj,errobj)
    outobj[i].data = ptr_new(outobj1)
    outarc1 = FIRE_EXTRACT_ORDER(bstr[i],imarc,errarc,/arc)
    outarc[i].data = ptr_new(outarc1)    
    ;stop
  Endfor

  stop

  ;; Don't use the RECTIFIED images, use the masked aperture image instead
  ;; use trace to fit, shift slightly in Y
  ;; make a boxcar version too
  
  ;; Write a program that will fit the trace using a bright star.

  data = *outobj[20].data
  FIRE_LINEFIT2D,data.recim,data.recerr,linestr1,model1
  ;; sometimes can pick up some more sky lines with a second round
  FIRE_LINEFIT2D,data.recim-model1,data.recerr,linestr2,model2

  
;; For an object spectrum, rectify the order, sum over all columns and
;; find the peak, use this to "recenter" the trace for extraction.
  
  
  data = *outarc[20].data
  FIRE_LINEFIT2D,data.recim,data.recerr,linestr1,model1
  
  ;; put in 2D spectrum array?
  
  stop

  end
