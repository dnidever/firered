function fire_makeimage,input,err=err,mask=mask,x=x,y=y,head=head

  ;; Make an image object

  ;;   flux, err, mask (good mask), nx, ny, x (global), y (global),
  ;;   type (arc, object, quartz)

  ;; Input filename
  if size(input,/type) eq 7 then begin
    file = input
    if file_test(input) eq 0 then begin
      print,input,' NOT FOUND'
      return,-1
    endif
    FITS_READ,input,flux,head
  endif else begin
     flux = input
     file = ''
  endelse
    
  sz = size(flux)
  nx = sz[1]
  ny = sz[2]

  ;; Fix "flipped" pixels
  bd = where(flux lt -1000,nbd)
  if nbd gt 0 then flux[bd] += 65536
  
  ;; Flux
  im = {file:file,flux:flux}
  ;; Gain
  gain = 1.0
  if n_elements(head) gt 0 then begin
    hdgain = sxpar(head,'gain',count=ngain)  ; electrons/DU     
    if ngain gt 0 then gain=hdgain
  endif
  ;; Read noise
  rdnoise = 0.0
  if n_elements(head) gt 0 then begin
    hdrdnoise = sxpar(head,'enoise',count=nrdnoise)  ; electrons/read
    if nrdnoise gt 0 then rdnoise=hdrdnoise
  endif  
  ;; Error array
  if n_elements(err) eq 0 then begin
     ;; poisson error in electrons = sqrt((flux>1)*gain)
     ;; rdnoise (in electrons)
     ;; add in quadrature
     ;; err[e] = sqrt( (flux>1)*gain) + rdnoise^2)
     ;; convert back to ADU
     ;; err[ADU] = sqrt( (flux>1)*gain) + rdnoise^2)/gain
     ;;          = sqrt( (flux>1)/gain + (rdnoise/gain)^2 )
     err = sqrt((flux>1)/gain + (rdnoise/gain)^2)
  endif
  im = create_struct(im,'err',err)
  ;; Mask
  if n_elements(mask) eq 0 then mask=bytarr(nx,ny)+1
  im = create_struct(im,'mask',mask)
  ;; X and Y arrays
  if n_elements(x) eq 0 then x=lindgen(nx)
  if n_elements(y) eq 0 then y=lindgen(ny)
  im = create_struct(im,'x',x,'y',y,'nx',nx,'ny',ny)
  ;; Header
  if n_elements(head) eq 0 then mkhdr,head,flux
  im = create_struct(im,'head',head)
  ;; Exposure type
  exptype = 'object'
  if n_elements(head) gt 0 then begin
    hdexptype = sxpar(head,'exptype',count=nexptype)
    if nexptype gt 0 then begin
      exptype = strlowcase(hdexptype)
    endif else begin
       if stregex(input,'bpm',/boolean,/fold_case) eq 1 then exptype='bpm'
    endelse
  endif
  if strlowcase(exptype) eq 'science' then exptype='object'
  im = create_struct(im,'exptype',exptype)
  ;; Subimage/rectified
  im = create_struct(im,'subimage',0,'rectified',0)
  
  return,im
  
  end
