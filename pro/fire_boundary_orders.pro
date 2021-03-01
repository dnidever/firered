pro fire_boundary_orders

  ;; Get boundaries for the orders

  ;;file = 'ut131222/fire_0001.fits'
  file = 'ut131222/fire_0011.fits'  
  fits_read,file,im,head

  norders = 22
  npix = 2048

  step = 25
  nbin = 100

  undefine,loall,hiall
  for i=0,npix/step-1 do begin
    xmed = i*step
    lo = i*step-nbin/2 > 0
    hi = i*step+nbin/2 < (npix-1) 
    ;hi =  lo+nbin-1
    ;xmed = lo+nbin/2
    med = median(im[lo:hi,*],dim=1)
    ;med = median(im[950:1050,*],dim=1)
    sm = smooth(med,3,/edge_truncate)
  
    ;; find edges
    dmed = slope(sm,/acc)
    sigdmed = mad(dmed)
    ;; mad maxima and minima
    maxima,dmed,allminarr,allmaxarr
    ;; sigma for each pixel
    bindata,findgen(npix),dmed,xbin,ybin,binsize=100,/mad
    interp,xbin,ybin,findgen(npix),sig
    nsig = 3 ;10
    
    ;; minima
    gmin = where(dmed[allminarr] lt -nsig*sig)
    minarr = allminarr[gmin]

    ;; neighboring pixels must also be negative
    gdmin = where(dmed[minarr-1] lt 0 and dmed[minarr+1] lt 0,ngdmin)
    minarr = minarr[gdmin]
    
    ;; deal with minima that are too close to each other
    bd = where(slope(minarr) lt 53,nbd)
    while (nbd gt 0) do begin
      bd = where(slope(minarr) lt 53,nbd)
      if nbd gt 0 then begin
        ind = [bd[0],bd[0]+1]
        si = sort(abs(dmed[minarr[ind]]))
        remove,ind[si[0]],minarr
        ;print,'removing ',si[0]
      endif
    endwhile
  
    ;; maxima
    gmax = where(dmed[allmaxarr] gt nsig*sig)
    maxarr = allmaxarr[gmax]  

    ;; neighboring pixels must also be positive
    gdmax = where(dmed[maxarr-1] gt 0 and dmed[maxarr+1] gt 0,ngdmax)
    maxarr = maxarr[gdmax]
    
    ;; deal with maxima that are too close to each other
    bd = where(slope(maxarr) lt 53,nbd)
    while (nbd gt 0) do begin
      bd = where(slope(maxarr) lt 53,nbd)
      if nbd gt 0 then begin
        ind = [bd[0],bd[0]+1]
        si = sort(abs(dmed[maxarr[ind]]))
        remove,ind[si[0]],maxarr
        ;print,'removing ',si[0]
      endif
    endwhile

    plot,dmed
    oplot,minarr,dmed[minarr],ps=1,co=250
    oplot,maxarr,dmed[maxarr],ps=1,co=150

    temp = replicate({x:0L,y:0L},n_elements(minarr))
    temp.x = xmed
    temp.y = minarr
    push,hiall,temp
    temp = replicate({x:0L,y:0L},n_elements(maxarr))
    temp.x = xmed
    temp.y = maxarr
    push,loall,temp
    
    ;if xmed eq 1000 then stop
    
    ;wait,0.5
    ;stop
    
  endfor

  displayc,im,min=0,max=2000
  oplot,loall.x,loall.y,ps=1,co=150
  oplot,hiall.x,hiall.y,ps=1,co=250

  ;; Create trace structure
  x = findgen(npix/step)*step
  trace = replicate({order:0,x:x,y0:fltarr(n_elements(x))-1,y1:fltarr(n_elements(x))-1,y0coef:fltarr(4),y1coef:fltarr(4),lo:0L,hi:0L},norders)
  trace.order = lindgen(22)+1
  
  ;; Initialize at x=1000
  ind = where(x eq 1000)
  indlo = where(loall.x eq 1000,nindlo)
  indhi = where(hiall.x eq 1000,nindhi)
  for i=0,norders-1 do trace[i].y0[ind] = loall[indlo[i]].y
  for i=0,norders-1 do trace[i].y1[ind] = hiall[indhi[i]].y
  
  ;; Now match them up
  ;; middle to right first
  ghi = where(x gt 1000.0,nghi)
  for i=0,nghi-1 do begin
    xind = ghi[i]
    xprev = xind-1 
    x1 = x[xind]
    ;; order loop
    for j=0,norders-1 do begin
      ;; Low
      if trace[j].y0[xprev] ge 0.0 then begin
        indlo = where(loall.x eq x[xind] and abs(loall.y-trace[j].y0[xprev]) lt 20,nindlo)
        if nindlo gt 1 then begin
           bestind = first_el(minloc(abs(loall[indlo].y-trace[j].y0[xprev])))
           indlo = indlo[bestind]
        endif
        if nindlo gt 0 then trace[j].y0[xind] = loall[indlo].y
      endif
      ;; High
      if trace[j].y1[xprev] ge 0.0 then begin
        indhi = where(hiall.x eq x[xind] and abs(hiall.y-trace[j].y1[xprev]) lt 20,nindhi)         
        if nindhi gt 1 then begin
           bestind = first_el(minloc(abs(hiall[indhi].y-trace[j].y1[xprev])))
           indhi = indhi[bestind]
        endif
        if nindhi gt 0 then trace[j].y1[xind] = hiall[indhi].y
      endif
    endfor
  endfor

  ;; middle to left
  glo = where(x lt 1000.0,nglo)
  for i=nglo-1,0,-1 do begin
    xind = glo[i]
    xprev = xind+1 
    x1 = x[xind]
    ;; order loop
    for j=0,norders-1 do begin
      ;; Low
      if trace[j].y0[xprev] ge 0.0 then begin
        indlo = where(loall.x eq x[xind] and abs(loall.y-trace[j].y0[xprev]) lt 20,nindlo)
        if nindlo gt 1 then begin
           bestind = first_el(minloc(abs(loall[indlo].y-trace[j].y0[xprev])))
           indlo = indlo[bestind]
        endif
        if nindlo gt 0 then trace[j].y0[xind] = loall[indlo].y
      endif
      ;; High
      if trace[j].y1[xprev] ge 0.0 then begin
        indhi = where(hiall.x eq x[xind] and abs(hiall.y-trace[j].y1[xprev]) lt 20,nindhi)         
        if nindhi gt 1 then begin
           bestind = first_el(minloc(abs(hiall[indhi].y-trace[j].y1[xprev])))
           indhi = indhi[bestind]
        endif
        if nindhi gt 0 then trace[j].y1[xind] = hiall[indhi].y
      endif
    endfor
  endfor

  ;; Now fit polynomials to the y positions
  for i=0,norders-1 do begin
    ;; lower boundary
    indlo = where(trace[i].y0 ge 0.0,nindlo)
    if i gt 0 then indlo = where(trace[i].y0 ge 10.0,nindlo)
    coef0 = robust_poly_fit(trace[i].x[indlo],trace[i].y0[indlo],3)
    trace[i].y0coef = coef0
    
    ;; upper boundary
    indhi = where(trace[i].y1 ge 0.0,nindhi)
    coef1 = robust_poly_fit(trace[i].x[indhi],trace[i].y1[indhi],3)
    trace[i].y1coef = coef1

    yr = minmax([trace[i].y0[indlo],trace[i].y1[indhi]])
    plot,trace[i].x[indlo],trace[i].y0[indlo],ps=8,yr=yr
    oplot,trace[i].x[indlo],poly(trace[i].x[indlo],coef0),co=250
    oplot,trace[i].x[indhi],trace[i].y1[indhi],ps=8
    oplot,trace[i].x[indhi],poly(trace[i].x[indhi],coef1),co=150

    trace[i].lo = max( [min(trace[i].x[indlo]),min(trace[i].x[indhi])] )
    trace[i].hi = min( [max(trace[i].x[indlo]),max(trace[i].x[indhi])] )
    ;stop
  endfor

  displayc,im,min=0,max=1000
  x = findgen(npix)
  for i=0,norders-1 do begin
    xx = x[trace[i].lo:trace[i].hi]
    oplot,xx,poly(xx,trace[i].y0coef),co=150
    oplot,xx,poly(xx,trace[i].y1coef),co=250  
  endfor

  ;mwrfits,trace,'fire_boundary_0011.fits',/create
  
  stop

  end
