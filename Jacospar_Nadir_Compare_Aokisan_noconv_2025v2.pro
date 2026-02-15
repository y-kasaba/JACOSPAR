pro Jacospar_Nadir_Compare_Aokisan_noconv_2025v2
  ;JACOSPARで計算した0.02cm-1刻みのスペクトルと比較する
  ;sp2:青木さん
  ;arnaud_sp;あるのーさん
  xsize=1000
  ysize=600
  window,xsize=xsize,ysize=ysize
  device, retain=2, decomposed=0,SET_FONT='DejaVuSans', /TT_FONT
  loadct,39
  !p.multi=[0,0,2,0,0]
  !p.background = 255
  !p.color = 0
  !p.charsize=2

  bin=1; Binning parameter to aokisan results : 0:no binning, 1:binning

  path_save='OMEGA/Results/448320251118_dwn1_0-80km_limbCase3_gas_test1_dz=1_ice2';'/Users/juice/Documents/Results_picture/Jacospar_Nadir_Compare_Aokisan/';

  ;dataName='/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4483_1_0216_newdustdatabase_newicedatabase_3500-3800_nogas_dust1.5ice2.0_H2O70000ppm'
  ;dataName=['/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4483_1_0217_newdustdatabase_newicedatabase_3500-3750_nogas_dust1.5ice2.0_dwn=1','/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4483_1_0217_newdustdatabase_newicedatabase_3750-30_nogas_dust1.5ice2.0_dwn=1']
  dataName='OMEGA/Results/448320251118_dwn1_0-80km_limbCase3_gas_test1_dz=1_ice2/Orbit_step00_4483_00_01.res'

  ;;;;;;;===== shohei results===========
  sp2 = read_binary('OMEGA/Results/0588_1_limb_ex1_dwn=001_inputfile=shoehi/Orbit_step00_0588_1_00_00.res',Data_type=4)
  wn2=3500d + dindgen(50000)*0.02d;data(0,*);[um]
  nwn2=n_elements(wn)

  ; change to radiance
  wn_pfsolspec_hr = fltarr(875482)
  rad_pfsolspec_hr = fltarr(875482)
  file_pfsolspec_hr = 'Data/Sun/pfsolspec_hr.dat'

  openr, 1, file_pfsolspec_hr
  for i = 0l, 875481 do begin
    readf, 1, a, b
    wn_pfsolspec_hr(i) = a
    rad_pfsolspec_hr(i) = b*!dpi/(1.5d^2.d) ;covert @ Mars [1.5 AU]
  endfor
  close,1

  solar_flux = interpol(rad_pfsolspec_hr,wn_pfsolspec_hr,wn2)
  sp2 = sp2 * solar_flux
  ;;;;==============

  ;;;;;;;===== RISEI results===========
  sp = read_binary(dataname,Data_type=4)
  wn=3500d + dindgen(1000);data(0,*);[um]
  nwn=n_elements(wn)

  ; change to radiance
  wn_pfsolspec_hr = fltarr(875482)
  rad_pfsolspec_hr = fltarr(875482)
  file_pfsolspec_hr = 'Data/Sun/pfsolspec_hr.dat'

  openr, 1, file_pfsolspec_hr
  for i = 0l, 875481 do begin
    readf, 1, a, b
    wn_pfsolspec_hr(i) = a
    rad_pfsolspec_hr(i) = b*!dpi/(1.5d^2.d) ;covert @ Mars [1.5 AU]
  endfor
  close,1

  solar_flux = interpol(rad_pfsolspec_hr,wn_pfsolspec_hr,wn)
  sp = sp * solar_flux
  ;;;;==============


  ;  ;---------------- my calculation results ---------------
  ;  ;data=read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4484_0_0807_nadirtest_4100-4200_nogas_nodust_80km/rad')
  ;  data=read_ascii(dataName[0]+'/Radnoconv.txt')
  ;  data=data.field1
  ;  ndata=n_elements(data)
  ;  arnaud_sp=data(1,*)
  ;  wn=data(0,*);[um]
  ;
  ;  ;  wn=1/wl*1E4
  ;  nwn=n_elements(wn)
  ;
  ;  ;  if n_elements(dataName) eq 2 then begin
  ;  ;    data2=read_ascii(dataName[1]+'/rad')
  ;  ;    data2=data2.field1
  ;  ;    ndata2=n_elements(data2)
  ;  ;    sp2=data2(0:ndata2/2-1)
  ;  ;    wl2=data2(ndata2/2:-1);[um]
  ;  ;    wn2=1/wl2*1E4
  ;  ;    nwn2=n_elements(wn2)
  ;  ;  endif
  ;
  ;
  ;  ; change to radiance
  ;  wn_pfsolspec_hr = fltarr(875482)
  ;  rad_pfsolspec_hr = fltarr(875482)
  ;  file_pfsolspec_hr = '/Users/juice/jacospar/Inversion/Data/Sun/pfsolspec_hr.dat'
  ;
  ;  openr, 1, file_pfsolspec_hr
  ;  for i = 0l, 875481 do begin
  ;    readf, 1, a, b
  ;    wn_pfsolspec_hr(i) = a
  ;    rad_pfsolspec_hr(i) = b*!dpi/(1.5d^2.d) ;covert @ Mars [1.5 AU]
  ;  endfor
  ;  close,1
  ;
  ;  solar_flux = interpol(rad_pfsolspec_hr,wn_pfsolspec_hr,wn)
  ;  arnaud_sp = arnaud_sp * solar_flux
  ;  arnaud_sp= arnaud_sp/!pi

  ;  wnnew=wn
  ;  spnew=sp
  ;
  ;  ;=======================
  ;  if n_elements(dataName) eq 2 then begin
  ;    solar_flux2 = interpol(rad_pfsolspec_hr,wn_pfsolspec_hr,wn2)
  ;    sp2 = sp2 * solar_flux2
  ;    sp2= sp2/!pi
  ;    wnnew=[wn,wn2]
  ;    spnew=[sp,sp2]
  ;  endif
  ;
  ;------------------- aoki-san results --------
  ;  A_data=read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Gas/Nadir_test_Pha30.dat')
  ;  A_data=A_data.field1
  ;  A_wn=A_data(0,25:-27)
  ;  A_nwn=n_elements(A_wn)
  ;  A_rad=A_data(1,0:-2)
  ;A_radnew=interpol(A_rad,A_wn,wnnew)
  if bin eq 1 then begin
    sp2_new=dblarr(1000)
    wn2_new=dblarr(1000)
    for i=0,999 do begin
      sp2_new(i)=0.02*(total(sp2(i*50l:i*50l+49l)))
      wn2_new(i)=mean(wn2(i*50l:i*50l+49l))
    endfor
  endif
  ;  ;===== あるのーさん、青木さんの計算結果を書き出す
  ;  ;あるのーさん
  ;  data1=dblarr(2,nwn)
  ;  data1(0,*)=wn
  ;  data1(1,*)=arnaud_sp
  ;  openw,lun,path_save+'Arnaud_database',/get_lun
  ;  printf,lun,data1
  ;  close,lun
  ;  free_lun,lun
  ;  ;青木さん積分後
  ;  data2=dblarr(2,nwn-1)
  ;  data2(0,*)=wn2_new
  ;  data2(1,*)=sp2_new
  ;  openw,lun,path_save+'Shohei_database_1cm-1',/get_lun
  ;  printf,lun,data2
  ;  close,lun
  ;  free_lun,lun
  ;  ;青木さん積分前
  ;  data3=dblarr(2,50000)
  ;  data3(0,*)=wn2
  ;  data3(1,*)=sp2
  ;  openw,lun,path_save+'Shohei_database_002cm-1',/get_lun
  ;  printf,lun,data3
  ;  close,lun
  ;  free_lun,lun

  ;========- プロット　========--
  xrange=[3500,4500];[min(wnnew),max(wnnew)]
  ;------------------- plot figure -------------
  plot,wn,sp2_new,xrange=xrange,yrange=[0,0.1],/nodata,color=0,xtitle='wavenumber',ytitle='W/sr/cm-1/m2'
  ;xyouts,0.5,0.5,'test',charsize=2.5,color=0,/normal
  ;plots,wn2_new,sp2_new,color=200
  plots,wn,sp,color=100


  differ=(sp-sp2_new)/sp*1e2
  plot,wn,differ,xrange=xrange,yrange=[-5,5],/nodata,color=0,xtitle='wavenumber',ytitle='%',title='Difference(1-DISORT/JACOSPAR)'
  ;xyouts,0.5,0.5,'test',charsize=2.5,color=0,/normal
  oplot,wn,differ,color=150
  oplot,[min(sp),max(sp)],[0,0]

  ;------------------- integrate the radiance -----
  ; integrate my radiance by wavenumber
  wn_reso=(max(wn)-min(wn))/(nwn-1)
  energy=total(wn_reso*sp)
  if n_elements(dataName) eq 2 then begin
    wn_reso2=(max(wn2)-min(wn2))/(nwn2-1)
    energy=total(wn_reso*sp)+total(wn_reso2*sp2)
  endif


  ;integrate aoki-san radiance by wavenumber
  A_wnreso=(max(wn2_new)-min(wn2_new))/(nwn-1)

  ;  lower=where(A_wn eq min(wn))
  ;  upper=where(A_wn eq max(wn))
  min1=min((wn - min(wn2_new))^2,lower)
  min2=min((wn - max(wn2_new))^2,upper)
  if n_elements(dataName) eq 2 then begin
    min3=min((wn2 - min(wn2))^2,lower2)
    min4=min((wn2 - max(wn2))^2,upper2)
  endif

  A_energy=total(A_wnreso*sp2_new(lower:upper))
  if n_elements(dataName) eq 2 then A_energy=total(A_wnreso*A_rad(lower:upper))+total(A_wnreso*A_rad(lower2:upper2))

  ;xyouts,0.2,0.87,'Aoki-san '+'energy(W/sr/m2)='+strmid(A_energy,6,5),color=200,/normal
  xyouts,0.2,0.9,'Risei '+'energy(W/sr/m2)='+strmid(energy,6,5),color=100,/normal




  ; ================================================ to compare jacospar different input ======================================================
  ;  dataName=['/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4484_0_0107_nadir_gas_validation_3500-3800_Pa','/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4484_0_0107_nadir_gas_validation_3950-4200_Pa']
  ;
  ;
  ;;data=read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4484_0_0807_nadirtest_4100-4200_nogas_nodust_80km/rad')
  ;data=read_ascii(dataName[0]+'/rad')
  ;data=data.field1
  ;ndata=n_elements(data)
  ;sp=data(0:ndata/2-1)
  ;wl=data(ndata/2:-1);[um]
  ;wn=1/wl*1E4
  ;
  ;data2=read_ascii(dataName[1]+'/rad')
  ;data2=data2.field1
  ;ndata2=n_elements(data2)
  ;sp2=data2(0:ndata2/2-1)
  ;wl2=data2(ndata2/2:-1);[um]
  ;wn2=1/wl2*1E4
  ;
  ;
  ;wn_pfsolspec_hr = fltarr(875482)
  ;rad_pfsolspec_hr = fltarr(875482)
  ;file_pfsolspec_hr = '/Users/juice/jacospar/Inversion/Data/Sun/pfsolspec_hr.dat'
  ;
  ;openr, 1, file_pfsolspec_hr
  ;for i = 0l, 875481 do begin
  ;  readf, 1, a, b
  ;  wn_pfsolspec_hr(i) = a
  ;  rad_pfsolspec_hr(i) = b*!dpi/(1.5d^2.d) ;covert @ Mars [1.5 AU]
  ;endfor
  ;close,1
  ;
  ;solar_flux = interpol(rad_pfsolspec_hr,wn_pfsolspec_hr,wn)
  ;sp = sp * solar_flux
  ;sp= sp/!pi
  ;
  ;solar_flux2 = interpol(rad_pfsolspec_hr,wn_pfsolspec_hr,wn2)
  ;sp2 = sp2 * solar_flux2
  ;sp2= sp2/!pi
  ;
  ;wnnew=[wn,wn2]
  ;spnew=[sp,sp2]
  ;
  ;plots,wnnew,spnew,color=200
  snapshot = TVRD(True=1)
  Write_JPEG, path_save+'/jacospar_nadir_compare_v2.jpg', snapshot, True=1, Quality=100
  if n_elements(dataName) eq 2 then Write_JPEG, path_save+strmid(dataname[0],58,60)+'_noconv.jpg', snapshot, True=1, Quality=100
end