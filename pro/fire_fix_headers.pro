pro fire_fix_headers

  ;; On the first two nights the object name was NOT put in the
  ;; header.  It was on the last night
  ;; But the coordinates are there. Match to the target catalog


  ;; There are offsets in the coordinates, use the offset patterns to
  ;; figure out what's what

  setdisp
  
  ;; Load the target catalog
  readcol,'Nidever_firetargets.cat',num,name,sra,sdec,format='i,a,a,a',skip=2
  ntargets = n_elements(name)
  targets = replicate({num:0,name:'',sra:'',sdec:'',ra:0.0d0,dec:0.0d0},ntargets)
  targets.num = long(num)
  targets.name = name
  targets.sra = sra
  targets.sdec = sdec
  targets.ra = sexig2ten(sra)*15
  targets.dec =sexig2ten(sdec)

  info = mrdfits('fire_info.fits',1)
  info.exptype = strtrim(info.exptype,2)
  gdsci = where(info.exptype eq 'Science',ngdsci)
  infosci = info[gdsci]
  
  ;; Night loop
  undefine,all
  schema = {file:'',night:'',expnum:0L,exptype:'',exptime:0.0,object:'',ra:0.0d0,dec:0.0d0,ra2:0.0d0,dec2:0.0d0}
  nights = ['ut131222','ut131223','ut131224']
  For i=0,2 do begin
    files = file_search(nights[i]+'/fire_????.fits',count=nfiles)
    ;; exposure loop
    For j=0,nfiles-1 do begin
      head = headfits(files[j])
      base = file_basename(files[j],'.fits')
      expnum = long((strsplit(base,'_',/extract))[1])
      ra = sxpar(head,'ra')
      dec = sxpar(head,'dec')
      exptype = sxpar(head,'exptype')
      exptime = sxpar(head,'exptime')
      object = sxpar(head,'object')
      temp = schema
      temp.file = files[j]
      temp.night = nights[i]
      temp.expnum = expnum
      temp.exptype = exptype
      temp.exptime = exptime
      temp.object = object
      temp.ra = ra
      temp.dec = dec
      temp.ra2 = ra - 0.06
      temp.dec2 = dec - 0.02
      push,all,temp
      if exptype eq 'Science' then begin
        print,files[j],' ',ra,' ',dec,' ',exptype,' ',exptime,' ',object
        ;; coordinates are shifted a bit
        ;;; ut131222
        ;ra -= 0.06
        ;dec -= 0.02
        ;; some of ut131223
        ;ra += 0.2
        ;dec -= 0.02        
        dist = sphdist(ra,dec,targets.ra,targets.dec,/deg)*3600
        ind = first_el(minloc(dist))
        plot,targets.ra,targets.dec,ps=1,xr=[-1,1]+ra,yr=[-1,1]+dec,xs=1,ys=1
        oplot,info.ra2,info.dec2,ps=4,co=150
        ;oplot,info.ra+0.2,info.dec,ps=4,co=150        
        oplot,[ra],[dec],ps=5,co=250,sym=2
        oplot,[targets[ind].ra],[targets[ind].dec],ps=6,co=200,sym=3
        if object eq 'Dark' then begin 
          txt = ''
          read,'Okay? ',txt
          if txt eq 'y' then begin
            sxaddpar,head,'object',targets[ind].name
            modfits,files[j],0,head
          endif
          ;;stop
        endif
      endif
      ;; Fix Arcs
      if exptype eq 'Arc' then begin
        if object eq 'FireStarter' or object eq 'Dark' then begin
          dist = sphdist(infosci.ra,infosci.dec,ra,dec,/deg)
          ind = first_el(minloc(dist))
          if min(dist) lt 0.1 and infosci[ind].night eq nights[i] and abs(expnum-infosci[ind].expnum) lt 5 then begin
             object = 'ThAr '+infosci[ind].object
             sxaddpar,head,'object',object
             modfits,files[j],0,head
             ;print,'match for arc ',files[j]
          endif else begin
             object = 'NO MATCH'
             sxaddpar,head,'object','ThAr'
             modfits,files[j],0,head             
             ;print,'NO science match for arc ',files[j]
          endelse
          print,files[j],' ',ra,' ',dec,' ',exptype,' ',exptime,' ',object
          ;stop
        endif
      endif
      ;; Fix quartz lamps
      if sxpar(head,'lamp') eq 'QuartzH'  or sxpar(head,'lamp') eq 'QuartzL' then begin
         print,files[j],' ',object,' QUARTZ!'
         if object eq 'FireStarter' or object eq 'Dark' then begin
            sxaddpar,head,'object',sxpar(head,'lamp')
            modfits,files[j],0,head
         endif
      endif
      ;; Fix tellurics
      if exptype eq 'Telluric' then begin
        if object eq 'FireStarter' or object eq 'Dark' then begin
            sxaddpar,head,'object','Telluric'
            modfits,files[j],0,head
        endif
      endif
    Endfor
    ;stop
  Endfor

; stopped on ut131223/fire_0170.fits
  
  ;mwrfits,all,'fire_info.fits',/create


  ;; Fix extra exposures based on the logsheets
  fix = importascii('fix_headers.txt',delim=string(9B),fieldnames=['num','object'])
  for i=0,n_elements(fix)-1 do begin
    file = file_search('ut13122?/fire_'+string(fix[i].num,format='(i04)')+'.fits',count=nfile)
    if nfile eq 0 then stop,'no file'
    head = headfits(file)
    object = sxpar(head,'object')
    print,num[i],' ',object,' -> ',fix[i].object
    sxaddpar,head,'object',fix[i].object
    modfits,file,0,head
  endfor

  
  
;; a handful of exposures from night 2 don't have IDs
;;ut131223/fire_0110.fits        42.677975       -19.506407 Science       10.6000 Dark  fixed by hand
;;ut131223/fire_0111.fits        42.677975       -19.506407 Science       21.2000 Dark  fixed by hand
;;ut131223/fire_0112.fits        42.677975       -19.506407 Science       42.4000 Dark  fixed by hand
;;ut131223/fire_0117.fits        57.769347       -42.523680 Science       10.6000 Dark  fixed by hand
;;ut131223/fire_0118.fits        57.769264       -42.523736 Science       10.6000 Dark  fixed by hand

;; got matches here
;; ut131224/fire_0204.fits        42.678038       -19.506412 Science       10.6000 twilight spectrum
;; ut131224/fire_0214.fits        57.770143       -42.523220 Science       10.6000 HD24331


;; Telluric stars from night 1
;;fire_0034.fits[2048,2048][real][Telluric][0][]: FireStarter        1.00000  2013-12-22  75.122283  -57.587020
;;fire_0042.fits[2048,2048][real][Telluric][0][]: FireStarter        10.6000  2013-12-22  77.328430  -54.980226
;;fire_0043.fits[2048,2048][real][Telluric][0][]: FireStarter        21.2000  2013-12-22  77.328472  -54.980226
;;fire_0049.fits[2048,2048][real][Telluric][0][]: FireStarter        10.6000  2013-12-22  42.339576  -43.836260
;;fire_0058.fits[2048,2048][real][Telluric][0][]: FireStarter        10.6000  2013-12-22  42.340157  -43.836900
;;fire_0064.fits[2048,2048][real][Telluric][0][]: FireStarter        10.6000  2013-12-22  92.993116  -60.835478
;;fire_0071.fits[2048,2048][real][Telluric][0][]: FireStarter        10.6000  2013-12-22  76.321636  -60.120272
;;fire_0077.fits[2048,2048][real][Telluric][0][]: FireStarter        10.6000  2013-12-22  85.798500  -60.862470

;; Telluric stars from night 2
;;fire_0125.fits[2048,2048][real][Telluric][0][]: Dark        10.6000  2013-12-23  70.979599  -49.877734   HD30252?
;;fire_0126.fits[2048,2048][real][Telluric][0][]: Dark        21.2000  2013-12-23  70.979682  -49.877734   HD30252?
;;fire_0138.fits[2048,2048][real][Telluric][0][]: Dark        10.6000  2013-12-23  120.68327  -65.631501
;;fire_0139.fits[2048,2048][real][Telluric][0][]: Dark        21.2000  2013-12-23  120.68339  -65.631529
;;fire_0150.fits[2048,2048][real][Telluric][0][]: Dark        10.6000  2013-12-23  145.73165  -51.229057
;;fire_0151.fits[2048,2048][real][Telluric][0][]: Dark        21.2000  2013-12-23  145.73169  -51.229002
;;fire_0162.fits[2048,2048][real][Telluric][0][]: Dark        10.6000  2013-12-23  142.10292  -71.718162
;;fire_0163.fits[2048,2048][real][Telluric][0][]: Dark        21.2000  2013-12-23  142.10313  -71.718190
;;fire_0174.fits[2048,2048][real][Telluric][0][]: Dark        10.6000  2013-12-23  95.734256  -70.954884
;;fire_0180.fits[2048,2048][real][Telluric][0][]: Dark        10.6000  2013-12-23  125.14269  -80.754328
  
;; Telluric from night 3
;;fire_0221.fits[2048,2048][real][Science][0][]: Telluric, HD42646        10.6000  2013-12-24  92.280630  -53.864380
;;fire_0226.fits[2048,2048][real][Telluric][0][]: Telluric, HD16636        10.6000  2013-12-24  39.624925  -52.929386
;;fire_0231.fits[2048,2048][real][Telluric][0][]: Telluric, HD23722        10.6000  2013-12-24  56.313567  -52.911322
;;fire_0241.fits[2048,2048][real][Telluric][0][]: Telluric, HD81484        10.6000  2013-12-24  140.65924  -67.572035
;;fire_0250.fits[2048,2048][real][Telluric][0][]: Telluric, HD30252        10.6000  2013-12-24  70.988125  -48.250704
;;fire_0259.fits[2048,2048][real][Telluric][0][]: Teeluric, HD62139        10.6000  2013-12-24  113.31385  -77.042429
;;fire_0265.fits[2048,2048][real][Telluric][0][]: Telluric, HD11213        10.6000  2013-12-24  204.85043  -73.507042
  
  stop

  end
