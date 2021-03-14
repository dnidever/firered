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
  nlsforder = 1
  nwaveorder = 3
  schema = {order:0,xlo:0L,xhi:0L,tycoef:fltarr(4),tsigcoef:fltarr(4),thwhmcoef:fltarr(4),tmoffcoef:fltarr(4),$
            trecenter:0.0,trescale:0.0,flux:fltarr(obj.nx),err:fltarr(obj.nx),mask:intarr(obj.nx),$
            wave:dblarr(obj.nx),wcoef:dblarr(nwaveorder+1),wrms:0.0,lsfcoef:dblarr(nlsforder+1),$
            nskylines:0,narclines:0}
  outstr = replicate(schema,norders)
  undefine,arclines
  undefine,skylines
  objmodel = obj.flux*0
  ;; skip first order, it has issues
  For i=1,norders-1 do begin
    print,'order = ',strtrim(i+1,2)
    outstr[i].order = i+1 
    ;; Recenter/scale aperture 
    tstr1 = FIRE_SCALE_TRACE(tstr[i],obj,/pl) ; recenter/scale
    outstr[i].xlo = tstr1.bndx0
    outstr[i].xhi = tstr1.bndx1
    outstr[i].trecenter = tstr1.recenter
    outstr[i].trescale = tstr1.rescale
    outstr[i].tycoef = tstr1.tycoef
    outstr[i].tsigcoef = tstr1.tsigcoef
    outstr[i].thwhmcoef = tstr1.thwhmcoef
    outstr[i].tmoffcoef = tstr1.tmoffcoef    

    ;; Get arc lines
    FIRE_LINEFIT2D,tstr1,arc,subarc,alinestr,amodel,aresidim,almodel,count=nalines,/arc
    if nalines gt 0 then begin
      add_tag,alinestr,'type','arc',alinestr
      add_tag,alinestr,'order',i,alinestr        
      push,arclines,alinestr
    endif
    print,'  ',strtrim(nalines,2),' arc lines found'
    sxaddhist,'order '+strtrim(i+1,2)+':  '+strtrim(nalines,2)+' arc lines found',head
    outstr[i].nskylines = nalines
    
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
    outstr[i].nskylines = nslines
    
    ;; Extract the object spectrum
    apim = FIRE_GETORDERIMAGE(tstr1,obj)
    FIRE_EXTRACT_ORDER,tstr1,subobj,spec,objmodel1,moffat=moffat
    objmodel[min(apim.x):max(apim.x),min(apim.y):max(apim.y)] += objmodel1
    outobj[i].data = ptr_new(spec)
    
    ;; Fill in the spectrum information
    outstr[i].flux[tstr1.bndx0:tstr1.bndx1] = spec.flux
    outstr[i].err[tstr1.bndx0:tstr1.bndx1] = spec.err
    outstr[i].mask[tstr1.bndx0:tstr1.bndx1] = 1

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
        wcoef = reform( poly_fit(alinestr2.xtrace,linelist2.wave,nwaveorder) )
        wrms = sqrt(mean((linelist2.wave-poly(alinestr2.xtrace,wcoef))^2))
        outstr[i].wrms = wrms
        print,'  Wave solution RMS = ',stringize(wrms,ndec=5),' pixels'
        sxaddhist,'order '+strtrim(i+1,2)+': Wave solution RMS = '+stringize(wrms,ndec=5)+' pixels',head
        
        xorder = lindgen(tstr1.bndx1-tstr1.bndx0+1)+tstr1.bndx0
        outstr[i].wave[tstr1.bndx0:tstr1.bndx1] = poly(xorder,wcoef)
        outstr[i].wcoef = wcoef
        ;if wrms1 gt 0.1 then stop,'bad solution!'
      endif
    endif else print,'  Nothing for this order in the linelist'

    ;; Get LSF coef each order
    lsfcoef = reform( robust_poly_fit(alinestr.xtrace,alinestr.pars[1],nlsforder) )
    ;glsf = where(alinestr.perror[1] gt 0.0,nglsf)
    ;if nglsf gt 2 then lsfcoef = robust_poly_fit(alinestr[glsf].xtrace,alinestr[glsf].pars[1],nlsforder)
    ;plot,alinestr.xtrace,alinestr.pars[1],ps=8
    ;x = findgen(2048)
    ;oplot,x,poly(x,lsfcoef),co=250
    sxaddhist,'order '+strtrim(i+1,2)+': mean LSF Gaussian sigma '+stringize(median(alinestr.pars[1]),ndec=3)+' pixels',head
    outstr[i].lsfcoef = lsfcoef
    
    ;stop
  Endfor  ;; order loop
  
  ;; Write out the information
  ;FIRE_WRITESPEC,spec,outfile
  print,'Writing spectrum to ',outfile
  MWRFITS,0,outfile,head,/create
  MWRFITS,outstr,outfile,/silent

end
