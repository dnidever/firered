pro fire_extract,objfile,arcfile,tracefile,bpmfile

  objfile = 'ut131222/fire_0046.fits'
  arcfile = 'ut131222/fire_0047.fits'
  tracefile = 'fire_trace_0011.fits'
  bpmfile = 'bpm2.fits'

  fits_read,objfile,imobj,objhead
  fits_read,arcfile,imarc,archead
  fits_read,bpmfile,imbpm,bpmhead
  tstr = mrdfits(tracefile,1)
  norders = n_elements(tstr)
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
  
  ;; Order loop
  x = findgen(npix)
  outobj = replicate({order:0,data:ptr_new()},norders)
  outobj.order = lindgen(norders)+1
  outarc = replicate({order:0,data:ptr_new()},norders)  
  outarc.order = lindgen(norders)+1
  ;For i=1,norders-1 do begin
  For i=20,20 do begin     
    print,'order = ',i
    outobj1 = FIRE_EXTRACT_ORDER(tstr[i],imobj,errobj)
    outobj[i].data = ptr_new(outobj1)
    outarc1 = FIRE_EXTRACT_ORDER(tstr[i],imarc,errarc,/arc)
    outarc[i].data = ptr_new(outarc1)    
    ;stop
  Endfor

  stop
  
  ;; Write a program that will fit the trace using a bright star.

  data = *outobj[20].data
  FIRE_LINEFIT2D,data.recim,data.recerr,linestr1,model1
  ;; sometimes can pick up some more sky lines with a second round
  FIRE_LINEFIT2D,data.recim-model1,data.recerr,linestr2,model2

  
  data = *outarc[20].data
  FIRE_LINEFIT2D,data.recim,data.recerr,linestr1,model1
  
  ;; put in 2D spectrum array?
  
  stop

  end
