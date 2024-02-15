FUNCTION CREATE_AEROSOLS_OPTICAL_kogure, specie, RMean = RMean, INDICES = indices, FORCE = force, _EXTRA = extra

  IF ~KEYWORD_SET(force) THEN force = 0
  pathnames = PATH_MANAGEMENT()
  
  CASE STRLOWCASE(specie) OF
    'dust': dustdata = READ_ASCII(pathnames.aerosols + 'Dust_rfi_Wolff_2009_in_wn.dat')
    'ice': data = READ_ASCII(pathnames.aerosols + 'waterice_warren.txt')
    ;'so2': data = TRANSPOSE([[11000D, 28600D],[1.34D, 1.34D],[0D, 0D]])
    ;'s2': data = TRANSPOSE([[11000D, 28600D],[1.96D, 1.96D],[0D, 0D]])
  ENDCASE
  
  
  
  IF STRCMP(specie, 'dust') OR STRCMP(specie, 'ice') THEN dustdata = DOUBLE(dustdata.FIELD1)

  IF ~KEYWORD_SET(Rmean) THEN Rmean = [DINDGEN(9) * 0.001D + 0.001D, DINDGEN(20) * 0.01D + 0.01D, DINDGEN(18) * 0.1D + 0.3D, DINDGEN(6) * 0.5D + 2.5D, DINDGEN(5) * 5D + 10D] * 1D-3
  IF KEYWORD_SET(indices) THEN Rmean = Rmean(indices)
  FOR i = 0, N_ELEMENTS(Rmean) - 1 DO BEGIN
    IF NOT(FILE_TEST(pathnames.aerosols + specie + '_' + STRING(Rmean(i), "(E8.2)") + '.aero')) OR force THEN BEGIN
      data = GET_AEROSOLS_OPTICAL_kogure(specie, Rmean(i), 0.3d, OPTICALFILE = dustdata, _EXTRA = extra)
     ; SAVE, data, FILENAME = pathnames.aerosols + specie + '_' + STRING(Rmean(i), "(E8.2)") + '.aero'
     labe=['   lamda [um] ','  AbsorpCrossSection ',' ScatCrossSection ']
     openw,lun,'/Users/juice/Documents/Results_picture/'+strmid(RMean,3,16)+'Rs=0.3d.txt',/get_lun
     printf,lun,labe
     printf,lun,data
     close,lun
     free_lun,lun

     ;     ------------ plot ------------
     erase
     loadct,39
     ;window;,xsize=xsize,ysize=ysize
     !p.background = 255
     !p.color = 0
     !p.charsize=2.5
     !p.multi=[0,6,0,0,0]
     
     wavelengs=1/data(0,*)*1e4
     wavelength2=1/dustdata(0,*)*1e4
     
     plot,wavelength2,dustdata(1,*),position=[0.1,0.7,0.3,0.95],yrange=[1.3,1.6],xrange=[0,3.0],color=0,charsize=2.0,title='refractive index Re(um)'
     ; plots,wavelengs,dustdata(1,*),color=100
     
     plot,wavelength2,dustdata(2,*),position=[0.4,0.7,0.6,0.95],yrange=[0,0.015],xrange=[0,2.5],color=0,charsize=2.0,title='refractive index Im(um)'
     ;plots,wavelengs,dustdata(2,*),color=100
     
     plot,wavelengs,data(2,*),position=[0.7,0.7,0.9,0.95],/nodata,color=0,charsize=2.0,title='sca/abs cross section'
     plots,wavelengs,data(2,*),color=254
     plots,wavelengs,data(1,*),color=0
     
     plot,data(0,*),data(2,*),position=[0.5,0.35,0.9,0.6],/nodata,color=0,charsize=2.0,title='sca/abs cross section'
     plots,data(0,*),data(2,*),color=254
     plots,data(0,*),data(1,*),color=0
     
     ratio=data(2,*)/data(1,*)
     plot,data(0,*),ratio,position=[0.1,0.35,0.4,0.6],yrange=[0.1,1E5],/nodata,color=0,charsize=2.0,/ylog,title='sca/abs ratio'
     plots,data(0,*),ratio,color=254
 
     albedo=data(2,*)/(data(2,*)+data(1,*))
     
     plot,wavelengs,albedo,position=[0.1,0.05,0.3,0.3],xrange=[0.2,3.1],yrange=[0.7,1.0],color=0,charsize=2.0,title='Albedo (um)'
;     plots,wavelengs,dustdata(1,*),color=200
     
     plot,data(0,*),ratio,position=[0.5,0.05,0.9,0.3],yrange=[0.5,1.1],/nodata,color=0,charsize=2.0,title='Albedo(cm-1)'
     plots,data(0,*),albedo,color=0
     
     snapshot = TVRD(True=1)
     Write_JPEG, '/Users/juice/Documents/Results_picture/'+strmid(RMean,3,16)+'Rs=0.3d_um.jpg', snapshot, True=1, Quality=100
     
     
;     --------------- plot end ----------

     
    ENDIF
  ENDFOR
END