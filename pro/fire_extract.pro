pro fire_extract,objfile,arcfile,bndfile=bndfile,tracefile=tracefile,$
                 moffat=moffat,outdir=outdir,clobber=clobber

  ;; Not enough inputs
  if n_elements(objfile) eq 0 or n_elements(arcfile) eq 0 then begin
    print,'Syntax - fire_extract,objfile,arcfile,bndfile=bndfile,tracefile=tracefile,'
    print,'                      moffat=moffat,outdir=outdir,clobber=clobber'
    return
  endif

  ;; Defaults
  if n_elements(outdir) eq 0 then outdir='finalspec/'
  ;;;objfile = 'ut131222/fire_0046.fits'
  ;objfile = 'ut131222/fire_0084.fits'
  ;arcfile = 'ut131222/fire_0047.fits'
  ;objfile = 'red/obj0084.fits'
  ;arcfile = 'red/arc0047.fits'
  if n_elements(bndfile) eq 0 then bndfile = 'fire_boundary_0011.fits'
  if n_elements(tracefile) eq 0 then tracefile = 'fire_trace_0084.fits'  
  if n_elements(moffat) eq 0 then moffat=1

  ;; Output filename
  base = file_basename(objfile,'.fits')
  outfile = outdir+base+'_spec.fits'
  if file_test(outfile) eq 1 and not keyword_set(clobber) then begin
    print,outfile,' exists and clobber not set'
    return
  endif

  print,'------------------------------'
  print,'Extracting ',objfile
  print,'------------------------------'  
  
  ;; Load the arc linelist
  repodir = '/Users/nidever/projects/firered/'
  linelist = importascii(repodir+'linelists/thar_linelist.txt',/header)
  
  ;; Load the data
  head = headfits(objfile)
  obj = fire_readimage(objfile)
  arc = fire_readimage(arcfile)
  bstr = mrdfits(bndfile,1)
  tstr = mrdfits(tracefile,1)  
  norders = n_elements(bstr)
  
  ;; Order loop
  ;;-----------
  outobj = replicate({order:0,data:ptr_new()},norders)
  outobj.order = lindgen(norders)+1
  outarc = replicate({order:0,data:ptr_new()},norders)  
  outarc.order = lindgen(norders)+1
  flux = fltarr(obj.nx,norders)
  err = fltarr(obj.nx,norders)
  mask = intarr(obj.nx,norders)
  wave = dblarr(obj.nx,norders)
  wcoef = fltarr(4,norders)
  wrms = dblarr(norders)
  lsfcoef = dblarr(3,norders)
  undefine,arclines
  undefine,skylines
  objmodel = obj.flux*0
  ;; skip first order, it has issues
  For i=1,norders-1 do begin
  ;For i=20,20 do begin     
    print,'order = ',strtrim(i+1,2)
    ;; Recenter/scale aperture 
    tstr1 = FIRE_SCALE_TRACE(tstr[i],obj) ; recenter/scale

    ;; Get arc lines
    FIRE_LINEFIT2D,tstr1,arc,subarc,alinestr,amodel,aresidim,almodel,count=nalines,/arc
    if nalines gt 0 then begin
      add_tag,alinestr,'type','arc',alinestr
      add_tag,alinestr,'order',i,alinestr        
      push,arclines,alinestr
    endif
    print,'  ',strtrim(nalines,2),' arc lines found'
    sxaddhist,'order '+strtrim(i+1,2)+':  '+strtrim(nalines,2)+' arc lines found',head
    
    arcrecim = FIRE_RECTIFY_ORDER(tstr1,arc,/exact)
    arcspec = arcrecim.flux[*,arcrecim.ny/2]
    
    ;; Remove sky lines
    FIRE_LINEFIT2D,tstr1,obj,subobj,slinestr,smodel,sresidim,slmodel,count=nslines
    if nslines gt 0 then begin
      add_tag,slinestr,'type','sky',slinestr
      add_tag,slinestr,'order',i,slinestr 
      push,skylines,slinestr
    endif
    print,'  ',strtrim(nslines,2),' sky lines found'
    sxaddhist,'order '+strtrim(i+1,2)+':  '+strtrim(nslines,2)+' sky lines found',head
    
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
    gline = where(linelist.order eq i+1,ngline)
    if ngline gt 0 then begin
      linelist1 = linelist[gline]
      ;; match the arc lines
      srcmatch,alinestr.xtrace,alinestr.xtrace*0,linelist1.xpix,linelist1.xpix*0,3.0,ind1,ind2,count=nmatch
      print,'  Matched ',strtrim(nmatch,2),' arc lines'
      sxaddhist,'order '+strtrim(i+1,2)+':  Matched '+strtrim(nmatch,2)+' arc lines',head
      
      if nmatch gt 0 then begin
        alinestr2 = alinestr[ind1]
        linelist2 = linelist1[ind2] 
        coef2 = poly_fit(alinestr2.xtrace,linelist2.wave,2)
        coef3 = poly_fit(alinestr2.xtrace,linelist2.wave,3)
        wrms1 = sqrt(mean((linelist2.wave-poly(alinestr2.xtrace,coef3))^2))
        wrms[i] = wrms1
        print,'  Wave solution RMS = ',stringize(wrms1,ndec=5),' pixels'
        sxaddhist,'order '+strtrim(i+1,2)+': Wave solution RMS = '+stringize(wrms1,ndec=5)+' pixels',head
        
        xorder = lindgen(tstr1.bndx1-tstr1.bndx0+1)+tstr1.bndx0
        wave[tstr1.bndx0:tstr1.bndx1,i] = poly(xorder,coef3)
        wcoef[*,i] = coef3
        ;if wrms1 gt 0.1 then stop,'bad solution!'
      endif
    endif else print,'  Nothing for this order in the linelist'

    ;; Get LSF coef each order
    lsfcoef1 = robust_poly_fit(alinestr.xtrace,alinestr.pars[1],1)
    ;glsf = where(alinestr.perror[1] gt 0.0,nglsf)
    ;if nglsf gt 2 then lsfcoef1 = robust_poly_fit(alinestr[glsf].xtrace,alinestr[glsf].pars[1],1)
    ;plot,alinestr.xtrace,alinestr.pars[1],ps=8
    ;x = findgen(2048)
    ;oplot,x,poly(x,lsfcoef1),co=250
    lsfcoef[0,i] = lsfcoef1
    sxaddhist,'order '+strtrim(i+1,2)+': mean LSF Gaussian sigma '+stringize(median(alinestr.pars[1]),ndec=3)+' pixels',head
    
    ;stop
  Endfor  ;; order loop

  ;lines = [arclines,skylines]
  ;gd = where(lines.pars[0] gt 0 and lines.status gt 0,ngd)
  ;lines = lines[gd]
  ;;;mwrfits,lines,'fire_lines.fits',/create
  ;base = file_basename(objfile,'.fits')
  ;MWRFITS,arclines,outdir+'/'+base+'_lines.fits',/create
  ;MWRFITS,skylines,outdir+'/'+base+'_lines.fits',/silent

  sxaddhist,'HDU0: Flux',head
  sxaddhist,'HDU1: Error',head  
  sxaddhist,'HDU2: Mask',head
  sxaddhist,'HDU3: Wavelength',head
  sxaddhist,'HDU4: Wavelength coefficients',head
  sxaddhist,'HDU5: LSF coefficients',head
  
  ;; Write out the information
  ;FIRE_WRITESPEC,spec,outfile
  print,'Writing spectrum to ',outfile
  MWRFITS,flux,outfile,head,/create
  MWRFITS,err,outfile,/silent
  MWRFITS,mask,outfile,/silent
  MWRFITS,wave,outfile,/silent
  MWRFITS,wcoef,outfile,/silent
  MWRFITS,lsfcoef,outfile,/silent
  ;stop

end
