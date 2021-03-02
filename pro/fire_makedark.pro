pro fire_makedark,files,outfile

  ;;  --- making dark ---

  files = 'ut131222/fire_00'+['18','19','20','21','22']+'.fits'
  outfile = 'dark.fits'
  
  nfiles = n_elements(files)
  npix = 2048
  cube = fltarr(npix,npix,nfiles)
  for i=0,nfiles-1 do begin
    fits_read,files[i],im,head
    cube[*,*,i] = float(im)
  endfor
  med = median(cube,dim=2)

  ;; scale by exptime
  
  stop
  
  med = median(im,11)
  npix = 2048
  ;; fix the edges
  med[0:5,*] = replicate(1,6)#med[6,*]
  med[npix-6:npix-1,*] = replicate(1,6)#med[npix-7,*]
  med[*,0:5] = med[*,6]#replicate(1,6)
  med[*,npix-6:npix-1] = med[*,npix-7]#replicate(1,6)
  sig = mad(im-med)
  nsig = 3
  mask = abs(im-med) gt nsig*sig
  ;; There are also bad pixels in of a "cross" that have low values
  cross = fltarr(3,3)
  cross[1,0] = 1
  cross[1,2] = 1  
  cross[0,1] = 1
  cross[2,1] = 1
  crossmask = convol(float(mask),cross)
  mask = mask or (crossmask eq 4)
  
  ;mwrfits,mask,outfile,/create

stop

end
