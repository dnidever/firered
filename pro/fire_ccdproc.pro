pro fire_ccdproc,files,bpm=bpm,dark=dark,flat=flat

  ;; Main CCD reduction steps
  nfiles = n_elements(files)
  print,strtrim(nfiles,2),' files to process'

  ;; Defaults
  if n_elements(bpm) eq 0 then bpm = 'bpm3.fits'  
  if n_elements(dark) eq 0 then dark = 'dark.fits'
  if n_elements(flat) eq 0 then flat = 'flat.fits'  

  ;; Load calibration products
  imbpm = FIRE_MAKEIMAGE(bpm)
  imfldark = FIRE_MAKEIMAGE(dark)
  imflat = FIRE_MAKEIMAGE(flat)  
  
  ;; File loop
  For i=0,nfiles-1 do begin
    file = files[i]
    if file_test(file) eq 0 then begin
      print,file,' NOT FOUND'
      goto,BOMB
    endif
     
    ;; Read in the file
    im = FIRE_MAKEIMAGE(file)

    ;; Fix wraparound pixels
    ;; this is done in makeimage
    ;bd = where(im.flux lt -1000,nbd)
    ;if nbd gt 0 then im.flux[bd] += 65536

    ;; Bad pixel mask
    im = FIRE_BPMCORRECT(im,imbpm)

    ;; Dark correction
    im = FIRE_DARKCORRECT(im,imdark)

    ;; Flat correction
    im = FIRE_FLATCORRECT(im,imflat)

    ;; Write output file

    BOMB:
  Endfor  ; file loop
     
  stop

  end
