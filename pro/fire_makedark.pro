pro fire_makedark,files,outfile

  ;;  --- making dark ---

  files1 = 'ut131222/fire_00'+['18','19','20','21','22']+'.fits'
  files2 = 'ut131223/fire_010'+['5','6','7','8','9']+'.fits'
  files3 = 'ut131224/fire_020'+['0','1','2']+'.fits'    
  files = [files1,files2,files3]
  outfile = 'dark.fits'
  
  nfiles = n_elements(files)
  npix = 2048
  cube = fltarr(npix,npix,nfiles)
  exptime = fltarr(nfiles)
  for i=0,nfiles-1 do begin
    fits_read,files[i],im,head
    exptime[i] = sxpar(head,'exptime')
    print,files[i],exptime[i],median(im),mad(im)
    newim = FIRE_REFCORRECT(im)
    ;; scale by exptime
    cube[*,*,i] = float(newim)/exptime[i]
    displayc,float(newim)/exptime[i],/z,tit=files[i]
    wait,0.3
  endfor
  med = median(cube,dim=3) > 0
  ;; must be positive


  
  ;med = median(im,11)
  ;npix = 2048
  ;;; fix the edges
  ;med[0:5,*] = replicate(1,6)#med[6,*]
  ;med[npix-6:npix-1,*] = replicate(1,6)#med[npix-7,*]
  ;med[*,0:5] = med[*,6]#replicate(1,6)
  ;med[*,npix-6:npix-1] = med[*,npix-7]#replicate(1,6)
  ;sig = mad(im-med)
  ;nsig = 3
  ;mask = abs(im-med) gt nsig*sig
  ;;; There are also bad pixels in of a "cross" that have low values
  ;cross = fltarr(3,3)
  ;cross[1,0] = 1
  ;cross[1,2] = 1  
  ;cross[0,1] = 1
  ;cross[2,1] = 1
  ;crossmask = convol(float(mask),cross)
  ;mask = mask or (crossmask eq 4)

  print,'Writing to ',outfile
  mkhdr,outhead,med
  sxaddpar,outhead,'exptime',1.0
  sxaddpar,outhead,'exptype','Dark'
  mwrfits,med,outfile,outhead,/create

  ;stop

end
