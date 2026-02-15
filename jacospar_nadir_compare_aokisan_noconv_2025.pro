pro Jacospar_Nadir_Compare_Aokisan_noconv_2025

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

  path_save='OMEGA/Results/003720251001_dwn10_0-100km_limbCase1/';

  ;dataName='/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4483_1_0216_newdustdatabase_newicedatabase_3500-3800_nogas_dust1.5ice2.0_H2O70000ppm'
  ;dataName=['/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4483_1_0217_newdustdatabase_newicedatabase_3500-3750_nogas_dust1.5ice2.0_dwn=1','/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4483_1_0217_newdustdatabase_newicedatabase_3750-30_nogas_dust1.5ice2.0_dwn=1']
  dataName='OMEGA/Results/003720251001_dwn10_0-100km_limbCase1'

  ;---------------- my calculation results ---------------
  ;data=read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/4484_0_0807_nadirtest_4100-4200_nogas_nodust_80km/rad')
  data=read_ascii(dataName[0]+'/rad')
  data=data.field1
  ndata=n_elements(data)
  sp=data(1,*)
  wn=data(0,*);[um]
  wn=1/wn*10000
  ;  sp = read_binary('/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/0685_1NADIRTEST_GASABSORPTIONTOSHOHEI_0519_3500-4500_dwn=002/Orbit_step00_0685_1_00_00.res',Data_type=4)
  ;  wn=3500d + dindgen(50000)*0.02d;data(0,*);[um]
  ;  wn=1/wl*1E4
  nwn=n_elements(wn)

  ;  if n_elements(dataName) eq 2 then begin
  ;    data2=read_ascii(dataName[1]+'/rad')
  ;    data2=data2.field1
  ;    ndata2=n_elements(data2)
  ;    sp2=data2(0:ndata2/2-1)
  ;    wl2=data2(ndata2/2:-1);[um]
  ;    wn2=1/wl2*1E4
  ;    nwn2=n_elements(wn2)
  ;  endif


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
  sp= sp/!pi

  wnnew=wn
  spnew=sp

  if n_elements(dataName) eq 2 then begin
    solar_flux2 = interpol(rad_pfsolspec_hr,wn_pfsolspec_hr,wn2)
    sp2 = sp2 * solar_flux2
    sp2= sp2/!pi
    wnnew=[wn,wn2]
    spnew=[sp,sp2]
  endif


  ;------------------- aoki-san results --------
  A_data=read_ascii('Data/Gas/Nadir_test_Pha30.dat')
  A_data=A_data.field1
  A_wn=A_data(0,25:-27)
  A_nwn=n_elements(A_wn)
  A_rad=A_data(1,0:-2)
  ;A_radnew=interpol(A_rad,A_wn,wnnew)
  if bin eq 1 then begin
    A_wn=A_data(0,25:-27)
    A_nwn=n_elements(A_wn)
    A_rad=A_data(1,25:-27)
    A_wn=rebin(A_wn,1,A_nwn/50)
    A_rad=rebin(A_rad,1,A_nwn/50)
    A_nwn=n_elements(A_wn)
    A_rad_compare=interpol(A_rad,A_wn,wnnew)
  endif

  xrange=[3500,4500];[min(wnnew),max(wnnew)]
  ;------------------- plot figure -------------
  plot,wnnew,spnew,xrange=xrange,yrange=[0,0.6],/nodata,color=0,xtitle='wavenumber',ytitle='W/sr/cm-1/m2'
  ;xyouts,0.5,0.5,'test',charsize=2.5,color=0,/normal
  plots,A_wn,A_rad,color=200
  plots,wnnew,spnew,color=100


  differ=(spnew-A_rad_compare)/spnew*1e2
  plot,wnnew,differ,xrange=xrange,yrange=[-3,3],/nodata,color=0,xtitle='wavenumber',ytitle='%',title='Difference(1-DISORT/JACOSPAR)'
  ;xyouts,0.5,0.5,'test',charsize=2.5,color=0,/normal
  oplot,wnnew,differ,color=150
  oplot,[min(A_wn),max(A_wn)],[0,0]

  ;------------------- integrate the radiance -----
  ; integrate my radiance by wavenumber
  wn_reso=(max(wn)-min(wn))/(nwn-1)
  energy=total(wn_reso*sp)
  if n_elements(dataName) eq 2 then begin
    wn_reso2=(max(wn2)-min(wn2))/(nwn2-1)
    energy=total(wn_reso*sp)+total(wn_reso2*sp2)
  endif


  ;integrate aoki-san radiance by wavenumber
  A_wnreso=(max(A_wn)-min(A_wn))/(A_nwn-1)

  ;  lower=where(A_wn eq min(wn))
  ;  upper=where(A_wn eq max(wn))
  min1=min((A_wn - min(wn))^2,lower)
  min2=min((A_wn - max(wn))^2,upper)
  if n_elements(dataName) eq 2 then begin
    min3=min((A_wn - min(wn2))^2,lower2)
    min4=min((A_wn - max(wn2))^2,upper2)
  endif

  A_energy=total(A_wnreso*A_rad(lower:upper))
  if n_elements(dataName) eq 2 then A_energy=total(A_wnreso*A_rad(lower:upper))+total(A_wnreso*A_rad(lower2:upper2))

  xyouts,0.2,0.9,'Aoki-san '+'energy(W/sr/m2)='+strmid(A_energy,6,5),color=200,/normal
  xyouts,0.2,0.87,'Risei '+'energy(W/sr/m2)='+strmid(energy,6,5),color=100,/normal




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
  Write_JPEG, path_save+strmid(dataname,58,60)+'.jpg', snapshot, True=1, Quality=100
  if n_elements(dataName) eq 2 then Write_JPEG, path_save+strmid(dataname[0],58,60)+'_noconv.jpg', snapshot, True=1, Quality=100
end