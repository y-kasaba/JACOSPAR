pro JACOSPAR_read_output_dusticedensradius_newplot_2024

  ;  device, retain=2, decomposed=0,SET_FONT='DejaVuSans', /TT_FONT
  ;  window,xsize=900,ysize=600
  ;  loadct,39
  ;  !P.Color = 0;'000000'xL
  ;  !P.Background = 'FFFFFF'xL
  ;  !p.charsize=2.5

;;-----------------------------------------resultファイルとアプリオリ、正解データの設定------------------------------------------------------------
  path='OMEGA/Results/003620240708_invsa_7000_7000_490_490_r1d-4_dodepend'
;  times = 1d;;;ダスト密度何倍か
  
  ;dens_correct = read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Apriori/apriori_4_2_2_GEMZ_wz_LOGND/mcd_dust_ave.dat') 
  ;dens_correct = read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Apriori/apriori_4_2_2_GEMZ_wz_LOGND/mcd_dust_cold.dat')
  dens_correct = read_ascii('Data/Apriori/apriori_4_2_2_GEMZ_wz_LOGND/LimbTest_Dust_one-tenth.dat')
  ;dens_correct = read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Apriori/apriori_4_2_2_GEMZ_wz_LOGND/mcd_dust_storm100times.dat')
  
  ice_dens_cor = read_ascii('Data/Apriori/apriori_4_2_2_GEMZ_wz_LOGND/NadirTest_H2O_ice.dat')
  
  radi = 0.8d-4;粒径正解データ[cm]
;;----------------------------------------------------------------------------------------------------------------------------------------
  
   
  path_result=path + '/Results.sav'
  path_resultZTPN=path + '/ResultsZTPN.sav'
  paths_TauAbs=file_search(path+'/Orbit_step*tau_abs.txt',count=count)
  paths_TauSca=file_search(path+'/Orbit_step*tau_aero_sca.txt',count=count)
  paths_abscross=file_search(path+'/Orbit_step*Aero_k_abs.txt',count=count)
  paths_scacross=file_search(path+'/Orbit_step*Aero_k_sca.txt',count=count)

  orbit = strmid(path,14,4)
  altitudepath = 'OMEGA/GEOMETRY/Geom_' + orbit + '.dat'

  restore,path_result
  ;restore,path_resultZTPN

  nlayer = nlayer
  iter = n_elements(convx)
  elements = nlayer*iter
  NSPECTRA = NSPECTRA
  NWN = NWN
  ngroups =n_elements(jac(*,0,0,0,0))
  nmembers = 1
  pathoutput = path + '/'
  naerosols = 2
  nojacobian = 0
  ;;=========== ジオメトリ ====================
  altitude=read_ascii(altitudepath,data_start=1)
  altitude=altitude.field1
  altitude=altitude(0,0:-1)

  ;;============ スペクトル　=========================
  file = file_search(path + '/Results.sav')

  wvl = file_search(path + '/rad')
  ;for loop = 0, count-1 do begin
  ;file = paths(loop)
  restore,file

  rad = rad(*,*,*)
  rad_measured = rad_measured
  wvl = read_ascii(wvl)
  wvl = wvl.field1(0,*)
  nalt = n_elements(rad(0,*,0))


;  count = iter

  ;!p.multi=[0,5,2,0,0]
  ;for loop=0, count-1 do begin
  ;  plot,wvl,rad_measured(*,0),color=0,yrange=[0,0.5],/nodata
  ;  for alt = 0, nalt-1 do begin
  ;    color=float(alt)/float(nalt)*254.0
  ;    oplot,wvl,rad_measured(*,alt),color=color,linestyle=2
  ;    oplot,wvl,rad(*,alt,loop), color = color
  ;    xyouts,0.2,0.8-alt*0.03,'alt step:'  + strmid(alt,6,2) ,color=color,/normal
  ;
  ;  endfor
  ;endfor
  ;snapshot = TVRD(True=1)
  ;Write_JPEG, path + '/Figures/Spectral_step' +strmid(loop,6,2)   + '.jpg', snapshot, True=1, Quality=100
  ;;========================= データ　==========================
  iter =n_elements(paths_TauAbs) - 1
  for loop=iter,iter do begin
    savename = 'Orbit_step' + STRING(iter, "(I02)") + '_' + orbit + '_00'
    R_std = dblarr(101,7)

    path_TauAbs=paths_TauAbs(loop)
    path_TauSca=paths_TauSca(loop)
    path_abscross=paths_abscross(loop)
    path_scacross=paths_scacross(loop)

    Tau_abs = read_ascii(path_TauAbs,data_start=1)
    tau_abs=tau_abs.field1
    TauSca = read_ascii(path_TauSca,data_start=1)
    TauSca=TauSca.field1
    abscross = read_ascii(path_abscross,data_start=1)
    abscross=abscross.field1
    scacross = read_ascii(path_scacross,data_start=1)
    scacross=scacross.field1

    ;;============== ヤコビアン　===========================
    ;  window,xsize=1400,ysize=700
    ;  device, retain=2, decomposed=0,SET_FONT='DejaVuSans', /TT_FONT
    ;  !p.multi=[0,1,0,0,0]
    ;  !p.charsize=2.5
    ;;  plot,Kmat(11,*,0),ztpn(0,*,0)
    ;  plot,Kmat(2*0,*,loop),ztpn(0,*,0),title='kmat(' + strmid(0,6,2) + ',*,loop)' + '_loop = ' + strmid(loop,7,1),/nodata,xrange=[-0.4,1.0]
    ;;  for i=0,n_elements(Kmat(*,0,0))/2-1 do begin
    ;;    oplot,Kmat(2*i,*,loop),ztpn(0,*,0),color=double(i)/n_elements(jac(0,0,*,0,0))*254d
    ;;    xyouts,0.06,0.85-0.03*i,string(round(altitude(i)))+'km',/normal,color=double(i)/n_elements(jac(0,0,*,0,0))*254d
    ;;  endfor
    ;  for i=0,n_elements(Kmat(*,0,0))/4-1 do begin
    ;    oplot,Kmat(4*i,*,loop),ztpn(0,*,0),color=double(i)/n_elements(jac(0,0,*,0,0))*254d
    ;    xyouts,0.06,0.85-0.03*i,string(round(altitude(i)))+'km',/normal,color=double(i)/n_elements(jac(0,0,*,0,0))*254d
    ;  endfor
    ;  snapshot = TVRD(True=1)
    ;  Write_JPEG, path+'/Figures/Kmat' + '_loop=' + strmid(loop,7,1) + '.jpg', snapshot, True=1, Quality=100
    ;  erase
    ;  !p.multi=[0,1,0,0,0]
    ;  !p.charsize=2.5
    ;;  plot,jac(0,*,0,1,loop),ztpn(0,*,0),title='jac(0,*,' + strmid(0,6,2) +',1,loop)' + '_loop = ' + strmid(loop,7,1),/nodata,xrange=[-0.02,0.05]
    ;;  for i=0,n_elements(jac(0,0,*,0,0))-1 do begin
    ;;    oplot,jac(0,*,i,1,loop),ztpn(0,*,0),color=double(i)/n_elements(jac(0,0,*,0,0))*254d
    ;;    xyouts,0.1,0.85-0.03*i,string(round(altitude(i)))+'km',/normal,color=double(i)/n_elements(jac(0,0,*,0,0))*254d
    ;;  endfor
    ;  plot,jac(0,*,0,1,loop),ztpn(0,*,0),title='jac(0,*,' + strmid(0,6,2) +',1,loop)' + '_loop = ' + strmid(loop,7,1),/nodata,xrange=[-0.2,0.2]
    ;  for i=0,n_elements(jac(0,0,*,0,0))-1 do begin
    ;    oplot,jac(0,*,i,1,loop),ztpn(0,*,0),color=double(i)/n_elements(jac(0,0,*,0,0))*254d
    ;    xyouts,0.1,0.85-0.03*i,string(round(altitude(i)))+'km',/normal,color=double(i)/n_elements(jac(0,0,*,0,0))*254d
    ;  endfor
    ;  snapshot = TVRD(True=1)
    ;  Write_JPEG, path+'/Figures/jac' + '_loop=' + strmid(loop,7,1) + '.jpg', snapshot, True=1, Quality=100
    ;  erase

    ;;============== Result data取得　===========================
    openw,1,path + '/checkData_step' + STRING(iter, "(I02)") + '.txt'
    ;  writef,1,'iter:' + strmid(loop)
    printf,1,'altitude'
    printf,1,format='(F20.10)',ztpn(0,*,0);,exp(xk(loop,2)),exp(xk(loop,4)),exp(xk(loop,8));リトリーバル変数 (イタレーション回数、ターゲット数＊大気数）
    printf,1,'retrieval variance'
    printf,1,format='(F20.10)',exp(xk(loop,*));,exp(xk(loop,2)),exp(xk(loop,4)),exp(xk(loop,8));リトリーバル変数 (イタレーション回数、ターゲット数＊大気数）
    printf,1,'output radiance'
    printf,1,format='(F20.10)',Fxk(loop,0);出力輝度 （イタレーション回数、大気数＊波長数）
    printf,1,'rad'
    printf,1,format='(F20.10)',rad(0,0,loop);出力輝度 （イタレーション回数、大気数＊波長数）
    printf,1,'Kmat'
    printf,1,format='(F20.17)',kmat(0,0,loop);,kmat(0,15,loop),kmat(0,30,loop),kmat(0,45,loop);Kmat : ヤコビアンの値（大気数＊波長数、ターゲット数＊大気数、イタレーション回数）
    printf,1,'convX'
    printf,1,format='(F20.10)',convX(loop);Xの事後確率密度関数（イタレーション回数）
    printf,1,'invsa'
    printf,1,format='(F20.10)',INVSA(0,0);,INVSA(15,15),INVSA(30,30),INVSA(45,45);
    printf,1,'invse'
    printf,1,format='(F20.10)',INVSE(*);,INVSE(15,15),INVSE(30,30),INVSE(45,45);
    printf,1,'Jac'
    printf,1,format='(F20.10)',Jac(0,0,0,*,loop);
    ;  printf,1,'reterror'
    ;  if loop gt 0 then printf,1,format='(F20.10)',reterror(loop-1,0);
    printf,1,'RMS'
    printf,1,format='(F20.10)',RMS(*,loop);
    printf,1,'S'
    printf,1,format='(F20.10)',S(0,0);,S(15,15),S(30,30),S(45,45);
    ;  printf,1,'SM'
    ;  printf,1,format='(F20.10)',SM(0,0);,SM(15,15),SM(30,30),SM(45,45);
    ;  printf,1,'SS'
    ;  printf,1,format='(F20.10)',SS(0,0);,SS(15,15),SS(30,30),SS(45,45);
    printf,1,'S_albedo'
    printf,1,format='(F20.10)',S_albedo(*,loop);
    printf,1,'X'
    printf,1,format='(F20.10)',X(*,0,loop);
    printf,1,'XA'
    printf,1,format='(F20.10)',XA(0);,XA(15),XA(30),XA(45);
    printf,1,'XK'
    printf,1,format='(F20.10)',XK(loop,0);,XK(loop,15),XK(loop,30),XK(loop,45);
    printf,1,'ZTPN'
    printf,1,format='(F20.10)',ZTPN(*,0,loop);
    printf,1,'Tau_abs'
    printf,1,format='(F20.10)',Tau_abs(0,0);
    ;  printf,1,'TauSca'
    ;  printf,1,format='(F20.10)',TauSca(1:2,0);
    printf,1,'abscross'
    printf,1,format='(F20.10)',abscross(1:2,0);
    printf,1,'scacross'
    printf,1,format='(F20.10)',scacross(1:2,0);


    R = DBLARR(ngroups * nmembers + 1, nspectra)
    ;  IF ~nojacobian THEN BEGIN
    dRdalb = DBLARR(ngroups * nmembers + 1, nspectra)
    dRdalb_std = DBLARR(ngroups * nmembers + 1, nspectra)
    dRdabs = DBLARR(ngroups * nmembers + 1, nlayer, nspectra)
    dRdabs_std = DBLARR(ngroups * nmembers + 1, nlayer, nspectra)
    dRdsca = DBLARR(ngroups * nmembers + 1, naerosols, nlayer, nspectra)
    dRdsca_std = DBLARR(ngroups * nmembers + 1, naerosols, nlayer, nspectra)
    ;  ENDIF
    FOR ifile = 0, nspectra - 1 DO BEGIN
      filename = savename + "_" + STRING(ifile, FORMAT = "(I02)") + '.res'
      data = READ_BINARY(pathoutput + filename, DATA_TYPE = 4)
      ind = 0
      indvec = INDGEN(ngroups)
      FOR i = 0, nmembers - 1 DO BEGIN
        R(i + indvec * nmembers, ifile) = data(ind * ngroups + indvec)
        ind += 1
        ;      R_std(i + indvec * nmembers, ifile) = data(ind * ngroups + indvec)
        ind += 1
        IF ~nojacobian THEN BEGIN
          dRdalb(i + indvec * nmembers, ifile) = data(ind * ngroups + indvec)
          ind += 1
          dRdalb_std(i + indvec * nmembers, ifile) = data(ind * ngroups + indvec)
          ind += 1
          FOR j = 0, nlayer - 1L DO BEGIN
            dRdabs(i + indvec * nmembers, j, ifile) = data(ind * ngroups + indvec)
            ind += 1
          ENDFOR
          FOR j = 0, nlayer - 1L DO BEGIN
            dRdabs_std(i + indvec * nmembers, j, ifile) = data(ind * ngroups + indvec)
            ind += 1
          ENDFOR
          FOR k = 0, naerosols - 1 DO BEGIN
            FOR j = 0, nlayer - 1L DO BEGIN
              dRdsca(i + indvec * nmembers, k, j, ifile) = data(ind * ngroups + indvec)
              ind += 1
            ENDFOR
            FOR j = 0, nlayer - 1L DO BEGIN
              dRdsca_std(i + indvec * nmembers, k, j, ifile) = data(ind * ngroups + indvec)
              ind += 1
            ENDFOR
          ENDFOR
        ENDIF
      ENDFOR
      R(ngroups * nmembers, ifile) = R(ngroups * nmembers - 1, ifile)
      ;    R_std(ngroups * nmembers, ifile) = R_std(ngroups * nmembers - 1, ifile)
    ENDFOR

    IF ~nojacobian THEN BEGIN
      dRdalb *= R
      ;    dRdalb_std *= R_std
      FOR i = 0, nlayer - 1L DO BEGIN
        dRdabs(*, i, *) *= R
        ;      dRdabs_std(*, i, *) *= R_std
        FOR k = 0, naerosols - 1 DO BEGIN
          dRdsca(*, k, i, *) *= R
          ;        dRdsca_std(*, k, i, *) *= R_std
        ENDFOR
      ENDFOR
      dRdalb(ngroups * nmembers, *) = dRdalb(ngroups * nmembers - 1, *)
      dRdalb_std(ngroups * nmembers, *) = dRdalb_std(ngroups * nmembers - 1, *)
      dRdabs(ngroups * nmembers, *, *) = dRdabs(ngroups * nmembers - 1, *, *)
      dRdabs_std(ngroups * nmembers, *, *) = dRdabs_std(ngroups * nmembers - 1, *, *)
      dRdsca(ngroups * nmembers, *, *, *) = dRdsca(ngroups * nmembers - 1, *, *, *)
      dRdsca_std(ngroups * nmembers, *, *, *) = dRdsca_std(ngroups * nmembers - 1, *, *, *)

      dRdalb(WHERE(~FINITE(dRdalb), /NULL)) = 0D
      dRdalb_std(WHERE(~FINITE(dRdalb_std), /NULL)) = 0D
      dRdabs(WHERE(~FINITE(dRdabs), /NULL)) = 0D
      dRdabs_std(WHERE(~FINITE(dRdabs_std), /NULL)) = 0D
      dRdsca(WHERE(~FINITE(dRdsca), /NULL)) = 0D
      dRdsca_std(WHERE(~FINITE(dRdsca_std), /NULL)) = 0D
    ENDIF
    printf,1,'R'
    printf,1,format='(F20.10)',R(0,0);
    printf,1,'R_std'
    printf,1,format='(E10.2)',R_std(0,0);
    printf,1,'drdalb'
    printf,1,format='(F20.10)',drdalb(0,0);
    printf,1,'dRdalb_std'
    printf,1,format='(E10.2)',dRdalb_std(0,0);
    printf,1,'dRdabs'
    printf,1,format='(F20.10)',dRdabs(0,0);
    printf,1,'dRdabs_std'
    printf,1,format='(E10.2)',dRdabs_std(0,0,0);
    printf,1,'dRdsca'
    printf,1,format='(F20.10)',dRdsca(0,0);
    printf,1,'dRdsca_std'
    printf,1,format='(E10.2)',dRdsca_std(0,0,0);
    ; writef,1,format='(3i5)',
    close,1
  endfor
  
  DEVICE, GET_SCREEN_SIZE = screen_size

;  if converged eq 1 then begin
  !p.multi=[0,2,4,0,0]
  dens_correct = dens_correct.Field1  &  dens_correct = interpol(dens_correct(1,*), dens_correct(0,*), ztpn(0,*))
  ice_dens_cor = ice_dens_cor.Field1  &  ice_dens_cor = interpol(ice_dens_cor(1,*), ice_dens_cor(0,*), ztpn(0,*))
  radi_correct = dblarr(n_elements(dens_correct)) & radi_correct(*) = radi

    
  p1 = PLOT(dens_correct, ztpn(0,*,0), 'g', XTITLE = 'Dust density(/cm3)', YTITLE = 'Altitude(km)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
    LAYOUT = [2, 4, 1], THICK = 2, NAME = 'Correct', /xlog, FONT_SIZE = 10)
  ;    p11 = PLOT(lambda_measured(*), rad_measured(*, 0) + drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT)
  ;    p12 = PLOT(lambda_measured(*), rad_measured(*, 0) - drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT)
  p1_1 = PLOT(exp(xk(iter,0:14)), ztpn(0,*,0), 'b', /OVERPLOT, THICK = 2, NAME = 'Retrieved', LINESTYLE = '-', FONT_SIZE = 10)
  p1_2 = PLOT(exp(xk(0,0:14)), ztpn(0,*,0), 'g', /OVERPLOT, THICK = 2, NAME = 'A priori', LINESTYLE = ':', FONT_SIZE = 10)
  l1 = LEGEND(TARGET = [p1, p1_1, p1_2], POSITION = [0.18, 0.93], FONT_SIZE = 7)
  p2 = PLOT(((exp(xk(iter,0:14)) - dens_correct)/dens_correct * 100),ztpn(0,*,0), 'b', /CURRENT, xrange = [-100, 200], XTITLE = 'Residual density(%)', YTITLE = 'Altitude(km)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
    LAYOUT = [2, 4, 2], THICK = 2, NAME = 'Measured', FONT_SIZE = 10)



  p3 = PLOT(ice_dens_cor, ztpn(0,*,0), 'g', /CURRENT, XTITLE = 'Ice density(/cm3)', YTITLE = 'Altitude(km)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
    LAYOUT = [2, 4, 3], THICK = 2, NAME = 'Correct', /xlog, FONT_SIZE = 10)
  ;    p11 = PLOT(lambda_measured(*), rad_measured(*, 0) + drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT)
  ;    p12 = PLOT(lambda_measured(*), rad_measured(*, 0) - drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT)
  p3_1 = PLOT(exp(xk(iter,15:29)), ztpn(0,*,0), 'b', /OVERPLOT, THICK = 2, NAME = 'Retrieved', LINESTYLE = '-', FONT_SIZE = 10)
  p3_2 = PLOT(exp(xk(0,15:29)), ztpn(0,*,0), 'g', /OVERPLOT, THICK = 2, NAME = 'A priori', LINESTYLE = ':', FONT_SIZE = 10)
  ;l3 = LEGEND(TARGET = [p3, p3_1, p3_2], POSITION = [0.38, 0.5])
  p2 = PLOT(((exp(xk(iter,15:29)) - ice_dens_cor)/ice_dens_cor * 100), ztpn(0,*,0), 'b', /CURRENT, xrange = [-100, 200], XTITLE = 'Residual density(%)', YTITLE = 'Altitude(km)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
    LAYOUT = [2, 4, 4], THICK = 2, NAME = 'Measured',  FONT_SIZE = 10)


  if n_elements(xk(iter,*)) gt 30 then begin
    
      p4 = PLOT(radi_correct, ztpn(0,*,0), 'g', /CURRENT, XTITLE = 'Dust Radius(cm)', YTITLE = 'Altitude(km)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
        LAYOUT = [2, 4, 5], THICK = 2, NAME = 'Correct', /xlog, FONT_SIZE = 10)
      ;    p11 = PLOT(lambda_measured(*), rad_measured(*, 0) + drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT)
      ;    p12 = PLOT(lambda_measured(*), rad_measured(*, 0) - drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT)
      p4_1 = PLOT(exp(xk(iter,30:44)), ztpn(0,*,0), 'b', /OVERPLOT, THICK = 2, NAME = 'Retrieved', LINESTYLE = '-', FONT_SIZE = 10)
      p4_2 = PLOT(exp(xk(0,30:44)), ztpn(0,*,0), 'g', /OVERPLOT, THICK = 2, NAME = 'A priori', LINESTYLE = ':', FONT_SIZE = 10)
     ; l4 = LEGEND(TARGET = [p1, p1_1, p1_2], POSITION = [0.38, 0.55], FONT_SIZE = 10)
      p5 = PLOT(((exp(xk(iter,30:44)) - radi_correct)/radi_correct * 100),ztpn(0,*,0), 'b', /CURRENT, xrange = [-100, 200], XTITLE = 'Residual density(%)', YTITLE = 'Altitude(km)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
        LAYOUT = [2, 4, 6], THICK = 2, NAME = 'Measured', FONT_SIZE = 10)
        
        
        
      p6 = PLOT(radi_correct, ztpn(0,*,0), 'g', /CURRENT, XTITLE = 'Ice Radius(cm)', YTITLE = 'Altitude(km)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
        LAYOUT = [2, 4, 7], THICK = 2, NAME = 'Correct', /xlog, FONT_SIZE = 10)
      ;    p11 = PLOT(lambda_measured(*), rad_measured(*, 0) + drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT)
      ;    p12 = PLOT(lambda_measured(*), rad_measured(*, 0) - drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT)
      p6_1 = PLOT(exp(xk(iter,45:59)), ztpn(0,*,0), 'b', /OVERPLOT, THICK = 2, NAME = 'Retrieved', LINESTYLE = '-', FONT_SIZE = 10)
      p6_2 = PLOT(exp(xk(0,45:59)), ztpn(0,*,0), 'g', /OVERPLOT, THICK = 2, NAME = 'A priori', LINESTYLE = ':', FONT_SIZE = 10)
     ; l6 = LEGEND(TARGET = [p1, p1_1, p1_2], POSITION = [0.38, 0.55], FONT_SIZE = 10)
      p7 = PLOT(((exp(xk(iter,45:59)) - radi_correct)/radi_correct * 100),ztpn(0,*,0), 'b', /CURRENT, xrange = [-100, 200], XTITLE = 'Residual density(%)', YTITLE = 'Altitude(km)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
        LAYOUT = [2, 4, 8], THICK = 2, NAME = 'Measured', FONT_SIZE = 10)
  endif


  p1.SAVE, path + "/Figures/iterated_density_3.jpg" , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
  p1.CLOSE
;  endif

  ;------------------ SPECTRA -------------------

  IF ~KEYWORD_SET(addname) THEN addname = ''

  RMS = DBLARR(nspectra)
  nwn = N_ELEMENTS(wn_measured)

  lambda_measured = 1D4 / wn_measured

  ;OPENW, lun, path + "Figures/Spectra_allAltitude.dat", /get_lun, WIDTH = 1000L

  p1 = PLOT(lambda_measured(*), rad_measured(*, 0), 'g', XTITLE = 'Wavelength ($\lambda [\mu m]$)', YTITLE = 'Radiance / Sun ', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
    TITLE = STRING(ztpn(0,1,0), "(F5.1)") + " km " , LAYOUT = [5, 3, 1], XRANGE = [MIN(lambda_measured), MAX(lambda_measured)], THICK = 2, NAME = 'Measured', FONT_SIZE = 10)
  p11 = PLOT(lambda_measured(*), rad_measured(*, 0) + drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT, FONT_SIZE = 10)
  p12 = PLOT(lambda_measured(*), rad_measured(*, 0) - drad_measured(*, 0), 'g', THICK = .5, /OVERPLOT, FONT_SIZE = 10)
  p3 = PLOT(lambda_measured(*), rad(*, 0, iter), 'b', /OVERPLOT, THICK = 2, NAME = 'Synthetic', LINESTYLE = '-', FONT_SIZE = 10)

  FOR i = 1, 13 DO BEGIN
    RMS(i) = SQRT(TOTAL((rad_measured(*, i) - rad(*, i, iter))^2))
    p1 = PLOT(lambda_measured(*), rad_measured(*, i), 'g', /CURRENT, XTITLE = 'Wavelength ($\lambda [\mu m]$)', YTITLE = 'Radiance / Sun ', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
      TITLE = STRING(ztpn(0,i+1,0), "(F5.1)") + " km " , LAYOUT = [5, 3, i+1], XRANGE = [MIN(lambda_measured), MAX(lambda_measured)], THICK = 2, NAME = 'Measured', FONT_SIZE = 10)
    p11 = PLOT(lambda_measured(*), rad_measured(*, i) + drad_measured(*, i), 'g', THICK = .5, /OVERPLOT, FONT_SIZE = 10)
    p12 = PLOT(lambda_measured(*), rad_measured(*, i) - drad_measured(*, i), 'g', THICK = .5, /OVERPLOT, FONT_SIZE = 10)
    p3 = PLOT(lambda_measured(*), rad(*, i, iter), 'b', /OVERPLOT, THICK = 2, NAME = 'Synthetic', LINESTYLE = '-', FONT_SIZE = 10)
  ENDFOR
  ;p5 = PLOT(lambda_measured(*), rad_measured(*, 1) - rad(*, 1, iter), 'g', /CURRENT, XTITLE = 'Wavelength ($\lambda [\mu m]$)', YTITLE = 'Residual radiance [CGS]', $
  ;    LAYOUT = [7, 2, 2], XRANGE = [MIN(lambda_measured), MAX(lambda_measured)], THICK = 2, NAME = 'RMS = ' + STRING(RMS(i), "(E10.2)"), LINESTYLE = '-', SYMBOL = 'x')

  l1 = LEGEND(TARGET = [p1, p3], POSITION = [0.98, 0.2])
  p1.SAVE, path + "/Figures/Spectra_allAltitude.jpg" , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
  p1.CLOSE

  ;------ ヤコビアン
  figProfiles = OBJARR(n_elements(kmat(*,0,0))/ngroups, 4)
  cc = COLORTABLE(39, NCOLORS = n_elements(kmat(*,0,0))/ngroups + 1, /TRANSPOSE)

  ;    figProfiles(0, 0) = PLOT(Kmat(ngroups,0:14,iter)*1.5, ztpn(0,*,0), 'g', XTITLE = 'Dust Jacobian', YTITLE = 'd(Radiance)/d(tau)', DIMENSIONS = ROUND(screen_size * .9), LINESTYLE = '-', $
  ;      LAYOUT = [2, 1, 1], THICK = 2, NAME = 'Dust Jacobian', /BUFFER, /nodata)

  for i=0,n_elements(kmat(*,0,0))/ngroups - 1 do begin
    figProfiles(i, 0) = PLOT(Kmat(ngroups*i,0:14,0)*1.5, ztpn(0,*,0), COLOR = cc(*,i), $
      LAYOUT = [2, 2, 1], TITLE = 'Dust Density Jacobian', DIMENSIONS = ROUND(screen_size * .9), FONT_SIZE = 10, $
      XTITLE = 'd(Radiance)/d(tau)', YTITLE = 'Altitude(km)', THICK = 2, CURRENT = (i GT 0) ? 1 : 0, OVERPLOT = (i GT 0) ? 1 : 0, NAME = [string(round(altitude(i)))+'km'], SYMBOL = 'x')
  endfor
  l = LEGEND(TARGET = figProfiles(*, 0), POSITION = [.55, .9])
  
  for i=0,n_elements(kmat(*,0,0))/ngroups - 1 do begin
    figProfiles(i, 1) = PLOT(Kmat(ngroups*i,15:29,0)*1.5, ztpn(0,*,0), COLOR = cc(*,i), $
      LAYOUT = [2, 2, 2], TITLE = 'Ice Density Jacobian', DIMENSIONS = ROUND(screen_size * .9), FONT_SIZE = 10, $
      XTITLE = 'd(Radiance)/d(tau)', YTITLE = 'Altitude(km)', THICK = 2, /CURRENT, OVERPLOT = (i GT 0) ? 1 : 0, NAME = [string(round(altitude(i)))+'km'], SYMBOL = 'x')
  endfor
  
  if n_elements(xk(iter,*)) gt 30 then begin
  
      for i=0,n_elements(kmat(*,0,0))/ngroups - 1 do begin
        figProfiles(i, 2) = PLOT(Kmat(ngroups*i,30:44,0),ztpn(0,*,0), COLOR = cc(*,i), $
          LAYOUT = [2, 2, 3], TITLE = 'Dust Radius Jacobian', FONT_SIZE = 10,  $
          xTITLE = 'd(Radiance)/d(Radius)', YTITLE = 'Altitude(km)', THICK = 2, /CURRENT, OVERPLOT = (i GT 0) ? 1 : 0, NAME = [string(round(altitude(i)))+'km'], SYMBOL = 'x')
      endfor
      ;l = LEGEND(TARGET = figProfiles(*, 1), POSITION = [.535, .9])
      
      for i=0,n_elements(kmat(*,0,0))/ngroups - 1 do begin
        figProfiles(i, 3) = PLOT(Kmat(ngroups*i,45:59,0),ztpn(0,*,0), COLOR = cc(*,i), $
          LAYOUT = [2, 2, 4], TITLE = 'Ice Radius Jacobian', FONT_SIZE = 10,  $
          xTITLE = 'd(Radiance)/d(Radius)', YTITLE = 'Altitude(km)', THICK = 2, /CURRENT, OVERPLOT = (i GT 0) ? 1 : 0, NAME = [string(round(altitude(i)))+'km'], SYMBOL = 'x')
      endfor
  endif
  
  figProfiles(0, 0).SAVE, path + '/Figures/Jacobian_iter0.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
  figProfiles(0, 0).CLOSE



end