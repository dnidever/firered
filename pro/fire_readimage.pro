function fire_readimage,file

  ;; Read an image

  ;;   flux, err, mask (good mask), nx, ny, x (global), y (global),
  ;;   type (arc, object, quartz)
  
  ;; File filename
  if file_test(file) eq 0 then begin
    print,file,' NOT FOUND'
    return,-1
  endif

  ;; Get number of extensions
  nextensions = 0
  fits_open,file,fcb
  nrextensions=fcb.nextend
  fits_close,fcb
  
  ;; Multiple extensions, load error, mask, etc.
  if nextensions gt 0 then begin
    FITS_READ,file,flux,head,exten=0
    FITS_READ,file,err,ehead,exten=1
    FITS_READ,file,mask,mhead,exten=2
    im = FIRE_MAKEIMAGE(flux,err=err,mask=mask,head=head)
    
  ;; Single image
  endif else begin
    FITS_READ,file,flux,head
    im = FIRE_MAKEIMAGE(flux,head=head)
  endelse
  im.file = file

  return,im
  
  end
