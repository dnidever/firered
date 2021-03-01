pro fire_extract,objfile,arcfile,tracefile,bpmfile

  ;objfile = 'ut131222/fire_0046.fits'
  objfile = 'ut131222/fire_0084.fits'
  arcfile = 'ut131222/fire_0047.fits'
  boundaryfile = 'fire_boundary_0011.fits'
  tracefile = 'fire_trace_0084.fits'  
  bpmfile = 'bpm3.fits'

  obj = fire_makeimage(objfile)
  arc = fire_makeimage(arcfile)
  bpm = fire_makeimage(bpmfile)
  bstr = mrdfits(boundaryfile,1)
  tstr = mrdfits(tracefile,1)  
  norders = n_elements(bstr)
  npix = 2048
  
  ;; Correct for bpm
  obj = fire_bpmcorrect(obj,bpm)
  arc = fire_bpmcorrect(arc,bpm)  

  ;; Correct for Dark current
  ;; scale by exptime
  
  ;; Do flat field correction
  
  ;; Order loop
  x = findgen(npix)
  outobj = replicate({order:0,data:ptr_new()},norders)
  outobj.order = lindgen(norders)+1
  outarc = replicate({order:0,data:ptr_new()},norders)  
  outarc.order = lindgen(norders)+1
  flux = fltarr(obj.nx,norders)
  err = fltarr(obj.nx,norders)
  mask = intarr(obj.nx,norders)
  undefine,arclines
  undefine,skylines
  For i=1,norders-1 do begin
  ;For i=8,20 do begin     
    print,'order = ',i
    ;; Recenter aperture 
    tstr1 = tstr[i]
    yrecenter = FIRE_RECENTER(tstr[i],obj)
    tstr1.tycoef[0] += yrecenter
    print,'Recenter = ',strtrim(yrecenter,2)
    ;; Get arc lines
    FIRE_LINEFIT2D,tstr1,arc,subarc,alinestr,amodel,aresidim,almodel,count=nalines,/arc
    if nalines gt 0 then push,arclines,alinestr
    print,strtrim(nalines,2),' arc lines found'
    ;; Remove sky lines
    FIRE_LINEFIT2D,tstr1,obj,subobj,slinestr,smodel,sresidim,slmodel,count=nslines
    if nslines gt 0 then push,skylines,slinestr
    print,strtrim(nslines,2),' sky lines found'
    ;; Extract the object spectrum
    spec = FIRE_EXTRACT_ORDER(tstr1,subobj)
    outobj[i].data = ptr_new(spec)
    
    ;; Fill in the spectrum information
    flux[tstr1.bndx0:tstr1.bndx1,i] = spec.flux
    err[tstr1.bndx0:tstr1.bndx1,i] = spec.err
    mask[tstr1.bndx0:tstr1.bndx1,i] = 1

    ;; The Gaussian PSF fits aren't that great
    ;; maybe scale PSF when doing recenter?
    
    ;stop
  Endfor  ;; order loop

  ;; 123 sky lines, 244 arc lines
  
  ;; Get LSF for each order

  ;; Get Wavelength for each order

stop
  
  ;; Write out the information
  outfile = 'xxxx.fits'
  MWRFITS,flux,outfile,head,/create
  MWRFITS,err,outfile
  MWRFITS,mask,outfile  

  stop

  end
