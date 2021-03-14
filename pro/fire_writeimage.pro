pro fire_writeimage,im,outfile

  ;; Write an image object to a file

  ;;   flux, err, mask (good mask), nx, ny, x (global), y (global),
  ;;   type (arc, object, quartz)

  MWRFITS,im.flux,outfile,im.head,/create
  MWRFITS,im.err,outfile,/silent
  MWRFITS,im.mask,outfile,/silent  

  ;; add other things like, exptype, subimage, rectified to header
  
  ;stop
  
  end
