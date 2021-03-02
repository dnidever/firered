pro fire_extract,objfile,arcfile,bndfile=bndfile,tracefile=tracefile,bpmfile=bpmfile,$
                 moffat=moffat,outdir=outdir

  if n_elements(outdir) eq 0 then outdir='finalspec/'
  
  ;objfile = 'ut131222/fire_0046.fits'
  objfile = 'ut131222/fire_0084.fits'
  arcfile = 'ut131222/fire_0047.fits'
  if n_elements(bndfile) eq 0 then bndfile = 'fire_boundary_0011.fits'
  if n_elements(tracefile) eq 0 then tracefile = 'fire_trace_0084.fits'  
  if n_elements(bpmfile) eq 0 then bpmfile = 'bpm3.fits'
  if n_elements(moffat) eq 0 then moffat=1

  ;; Load the data
  obj = fire_makeimage(objfile)
  arc = fire_makeimage(arcfile)
  bpm = fire_makeimage(bpmfile)
  bstr = mrdfits(bndfile,1)
  tstr = mrdfits(tracefile,1)  
  norders = n_elements(bstr)
  
  ;; Correct for bpm
  obj = fire_bpmcorrect(obj,bpm)
  arc = fire_bpmcorrect(arc,bpm)  

  ;; Correct for Dark current
  ;; scale by exptime
  
  ;; Do flat field correction
  
  ;; Order loop
  outobj = replicate({order:0,data:ptr_new()},norders)
  outobj.order = lindgen(norders)+1
  outarc = replicate({order:0,data:ptr_new()},norders)  
  outarc.order = lindgen(norders)+1
  flux = fltarr(obj.nx,norders)
  err = fltarr(obj.nx,norders)
  mask = intarr(obj.nx,norders)
  undefine,arclines
  undefine,skylines
  objmodel = obj.flux*0
  For i=1,norders-1 do begin
  ;For i=20,20 do begin     
    print,'order = ',i
    ;; Recenter/scale aperture 
    tstr1 = FIRE_SCALE_TRACE(tstr[i],obj) ; recenter/scale
    ;; Get arc lines
    FIRE_LINEFIT2D,tstr1,arc,subarc,alinestr,amodel,aresidim,almodel,count=nalines,/arc
    if nalines gt 0 then push,arclines,alinestr
    print,strtrim(nalines,2),' arc lines found'
    ;; Remove sky lines
    FIRE_LINEFIT2D,tstr1,obj,subobj,slinestr,smodel,sresidim,slmodel,count=nslines
    if nslines gt 0 then push,skylines,slinestr
    print,strtrim(nslines,2),' sky lines found'
    ;; Extract the object spectrum
    apim = FIRE_GETORDERIMAGE(tstr1,obj)
    FIRE_EXTRACT_ORDER,tstr1,subobj,spec,objmodel1,moffat=moffat
    objmodel[min(apim.x):max(apim.x),min(apim.y):max(apim.y)] += objmodel1
    outobj[i].data = ptr_new(spec)
    
    ;; Fill in the spectrum information
    flux[tstr1.bndx0:tstr1.bndx1,i] = spec.flux
    err[tstr1.bndx0:tstr1.bndx1,i] = spec.err
    mask[tstr1.bndx0:tstr1.bndx1,i] = 1

    ;; The Gaussian PSF fits aren't that great
    ;; maybe scale PSF when doing recenter?
    ;; Maye try Moffat function instead
    
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
