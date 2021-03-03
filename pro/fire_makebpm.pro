pro fire_makebpm,file,outfile

  file = 'ut131222/fire_0019.fits'
  outfile = 'bpm3.fits'
  
  ;;  --- making bpm ---
  fits_read,file,im,head
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
  
  print,'Writing to ',outfile
  mkhdr,outhead,mask
  sxaddpar,outhead,'exptype','BPM'
  mwrfits,mask,outfile,outhead,/create

stop

end
