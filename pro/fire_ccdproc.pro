pro fire_ccdproc,input,bpm=bpm,dark=dark,flat=flat,outdir=outdir

  ;; Load the input
  LOADINPUT,input,files,count=nfiles
  
  ;; Main CCD reduction steps
  nfiles = n_elements(files)
  print,strtrim(nfiles,2),' files to process'

  ;; Defaults
  if n_elements(bpm) eq 0 then bpm = 'bpm3.fits'  
  if n_elements(dark) eq 0 then dark = 'dark.fits'
  if n_elements(flat) eq 0 then flat = 'flat.fits'  
  if n_elements(outdir) eq 0 then outdir = 'red/'
  if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
  print,'BPM file = ',bpm
  print,'Dark file = ',dark
  print,'Flat file = ',flat
  print,'Output directory = ',outdir
  
  ;; Load calibration products
  imbpm = FIRE_READIMAGE(bpm)
  imdark = FIRE_READIMAGE(dark)
  imflat = FIRE_READIMAGE(flat)  
  
  ;; File loop
  For i=0,nfiles-1 do begin
    file = files[i]
    if file_test(file) eq 0 then begin
      print,file,' NOT FOUND'
      goto,BOMB
    endif

    ;; Load the image
    FITS_READ,file,flux,head
    sz = size(flux)
    nx = sz[1]
    ny = sz[2]

    object = sxpar(head,'object')
    exptype = sxpar(head,'exptype')
    exptime = sxpar(head,'exptime')    
    print,strtrim(i+1,2),' ',file,' ',object,' ',exptype,' ',exptime

    ;; Fix "flipped" pixels
    bd = where(flux lt -1000,nbd)
    if nbd gt 0 then flux[bd] += 65536

    ;; Reference pixel correct
    flux = FIRE_REFCORRECT(flux)
    flux[*,0:3] = 0.0
    flux[*,2044:*] = 0.0
    flux[0:3,*] = 0.0
    flux[2044:*,*] = 0.0    

    ;; Make the error array
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

    ;; Make the image object
    im = FIRE_MAKEIMAGE(flux,err=err,head=head)
    
    ;; Bad pixel mask
    im = FIRE_BPMCORRECT(im,imbpm)
    
    ;; Dark correction
    im = FIRE_DARKCORRECT(im,imdark)
    
    ;; Flat correction
    im = FIRE_FLATCORRECT(im,imflat)

    ;; Write output file
    base = file_basename(file,'.fits')
    num = strmid(base,5)
    ;; use useful prefix, obj, arc, tell
    suff = 'im'
    if strlowcase(exptype) eq 'science' then suff='obj'
    if strlowcase(exptype) eq 'arc' then suff='arc'
    if strlowcase(exptype) eq 'telluric' then suff='tell'    
    ;if strlowcase(object) eq 'thne' or strlowcase(object) eq 'thar' then suff='arc'
    ;if stregex(object,'telluric',/boolean,/fold_case) eq 1 then suff='tell'
    outfile = outdir+'/'+suff+num+'.fits'
    
    print,'Writing reduced file to ',outfile
    FIRE_WRITEIMAGE,im,outfile
    
    BOMB:
  Endfor  ; file loop
     
  ;stop

  end
