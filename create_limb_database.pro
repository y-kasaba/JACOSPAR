  ;------------------------------------------------------------------------------------------------------------------------
function WHERE_XYZ, Array_expression, Count, XIND=xind, YIND=yind, ZIND=zind
  ;------------------------------------------------------------------------------------------------------------------------
  ; works for 1, 2 or 3 dimensional arrays
  ;
  ; Returns the 1D indices (same as WHERE)
  ;
  ; ARGUMENTS
  ;  - same as WHERE (see WHERE)
  ;
  ; KEYWORDS
  ; - Optionally returns X, Y, Z locations through:
  ;
  ; XIND: Output keyword, array of locations
  ; along the first dimension
  ; YIND: Output keyword, array of locations
  ; along the second dimension (if present)
  ; ZIND: Output keyword, array of locations
  ; along the third dimension (if present)
  ;
  ; If no matches where found, then XIND returns -1
  ;
  index_array=where(Array_expression, Count, /L64)
  dims=size(Array_expression,/dim)
  xind=index_array mod dims[0]
  case n_elements(dims) of
    2: yind=index_array / dims[0]
    3: begin
      yind=index_array / dims[0] mod dims[1]
      zind=index_array / dims[0] / dims[1]
    end
    else:
  endcase
  return, index_array
end

Pro create_limb_database

;++++++++++++++++++++++++
;Written by SHOHEI AOKI
;++++++++++++++++++++++++

;================================================
;path
;================================================
path_data = '/home/pat1/WORK02/omega/data/'
path_save = '/home/pat1/WORK02/omega/limb_data/'

;================================================
;read_list
;================================================
files = strarr(1160)

name = ''
a=0
b=0
openr, lun, '/home/pat1/WORK02/omega/list_Limb.txt' ,/get_lun
for i = 0, 1160-1 do begin
  readf,lun,name,a,b,format='(A9,I,I)'
  files(i) = name
endfor
free_lun,lun

;================================================
;create complete list
;================================================
openw, lun, path_save + 'List_Limb_observations_OMEGA.txt' ,/get_lun
printf, lun, 'Filemame, Orbit, MY, Long, Lat, Ls, Slant-D, FOV, Total#-of-spectra, #-of-spectra at 0-5, 5-10, 10-15, 15-20, 20-25, 25-30, 30-35, 35-40, 40-45, 45-50 km'


;================================================
;search file
;================================================
;files = file_search(path_data+'*.sav',count=count)
;total_filen = count

;================================================
;Loop start
;================================================
for Loop = 0l, 1160-1 do begin

  file = files(Loop)
  print, file
  restore, path_data+file+'.sav'
  ip = n_elements(LATI(*,0))
  io = n_elements(LATI(0,*))

  ;reform(geocube(*,12,*)); in m
  tgh = reform(geocube(*,12,*))/1000. ;km

  a = WHERE_XYZ(tgh ge 65.535, count, XIND=XIND, YIND=YIND)
  total_number_of_spectra = count

  b = where(tgh ge 65.535 and tgh le 65.535 + 5., count)
  nsp_0_5 = count
  b = where(tgh ge 65.535 + 5. and tgh le 65.535 + 10., count)
  nsp_5_10 = count
  b = where(tgh ge 65.535 + 10. and tgh le 65.535 + 15., count)
  nsp_10_15 = count
  b = where(tgh ge 65.535 + 15. and tgh le 65.535 + 20., count)
  nsp_15_20 = count
  b = where(tgh ge 65.535 + 20. and tgh le 65.535 + 25., count)
  nsp_20_25 = count
  b = where(tgh ge 65.535 + 25. and tgh le 65.535 + 30., count)
  nsp_25_30 = count
  b = where(tgh ge 65.535 + 30. and tgh le 65.535 + 35., count)
  nsp_30_35 = count
  b = where(tgh ge 65.535 + 35. and tgh le 65.535 + 40., count)
  nsp_35_40 = count
  b = where(tgh ge 65.535 + 40. and tgh le 65.535 + 45., count)
  nsp_40_45 = count
  b = where(tgh ge 65.535 + 45. and tgh le 65.535 + 50., count)
  nsp_45_50 = count
  
  a = WHERE_XYZ(tgh ge 65.535 and tgh le 65.535 + 50., count, XIND=XIND, YIND=YIND)

  lati = reform(lati(XIND,YIND))
  longi = reform(longi(XIND,YIND))
  dmars = reform(dmars(XIND,YIND))
  Solar_longitude = reform(Solar_longitude(XIND,YIND))
  SUB_SOLAR_LONGITUDE = reform(SUB_SOLAR_LONGITUDE(XIND,YIND))

  jdat_save = fltarr(total_number_of_spectra,352)
  geocube_save = fltarr(total_number_of_spectra,51)

;  for i = 0l, total_number_of_spectra-1 do jdat_save(i,*) = reform(jdat(XIND(i),*,YIND(i)))
;  for i = 0l, total_number_of_spectra-1 do geocube_save(i,*) = reform(geocube(XIND(i),*,YIND(i)))
  for i = 0l, count-1 do jdat_save(i,*) = reform(jdat(XIND(i),*,YIND(i)))
  for i = 0l, count-1 do geocube_save(i,*) = reform(geocube(XIND(i),*,YIND(i)))
  jdat = jdat_save
  geocube = geocube_save

  slant_distance = reform(geocube(*,11)*1.e-3)
  FOV = slant_distance*tan(1.2/(1000.*2.))*2.

  save, wvl, jdat, lati, longi, geocube, dmars, Solar_longitude, SUB_SOLAR_LONGITUDE, filename=path_save+file_basename(file,'.sav')+'.sav'

  ;'Filemame, Orbit, MY, Long, Lat, Ls, Slant-D, FOV, Total#-of-spectra, #-of-spectra at 0-5, 5-10, 10-15, 15-20, 20-25, 25-30, 30-35, 35-40, 40-45, 45-50 km'
  printf, lun, STRCOMPRESS(file) + '  ' + $
    STRCOMPRESS(median(longi)) + '  ' + $
    STRCOMPRESS(median(lati)) + '  ' + $
    STRCOMPRESS(median(Solar_longitude)) + '  ' + $
    STRCOMPRESS(median(slant_distance)) + '  ' + $
    STRCOMPRESS(median(FOV)) + '  ' + $
    STRCOMPRESS(total_number_of_spectra) + '  ' + $
    STRCOMPRESS(nsp_0_5) + '  ' + $
    STRCOMPRESS(nsp_5_10) + '  ' + $
    STRCOMPRESS(nsp_10_15) + '  ' + $
    STRCOMPRESS(nsp_15_20) + '  ' + $
    STRCOMPRESS(nsp_20_25) + '  ' + $
    STRCOMPRESS(nsp_25_30) + '  ' + $
    STRCOMPRESS(nsp_30_35) + '  ' + $
    STRCOMPRESS(nsp_35_40) + '  ' + $
    STRCOMPRESS(nsp_40_45) + '  ' + $
    STRCOMPRESS(nsp_45_50)

endfor
free_lun,lun
close,/all
stop
end
