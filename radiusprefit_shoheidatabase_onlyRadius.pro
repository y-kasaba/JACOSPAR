FUNCTION RADIUSPREFIT_shoheidatabase_onlyRadius, orbit, alt_range, ztpn, geometry, FOV, wnmeas, radmeas, StartVariance, $
  DOPLOTPREFIT = doplotprefit, $
  NAMEPLUSDIR = nameplusdir, $
  DOPRINTPREFIT = doprintprefit, $
  PRINTSTAGE = printstage, $
  _EXTRA = extra
;
; ver_onlyRadius : only fit radius, not variance.20241012 kogure.
;

  IF ~KEYWORD_SET(doplotprefit) THEN doplotprefit = 0L
  IF ~KEYWORD_SET(nameplusdir) THEN nameplusdir = ''
  IF ~KEYWORD_SET(doprintprefit) THEN doprintprefit = 0L
  IF ~KEYWORD_SET(printstage) THEN printstage = 0L

  IF doplotprefit THEN DEVICE, GET_SCREEN_SIZE = screen_size

  pathnames = PATH_MANAGEMENT()
  pathoutput = pathnames.results + orbit + nameplusdir + '/Figures/'
  IF ~FILE_TEST(pathoutput, /DIRECTORY) THEN FILE_MKDIR, pathoutput

  radmeas = radmeas(WHERE(wnmeas GT 4000D AND wnmeas LT 20000D), *)
  wnmeas = wnmeas(WHERE(wnmeas GT 4000D AND wnmeas LT 20000D))

  dumb = MIN(ABS(wnmeas - 9000D), indref)

  meanrad = DINDGEN(30L) * 0.1D-4 + 0.1D-4
  sigrad_fixed = 0.3D;StartVariance(0)  ; 分散を1つの固定値に設定
  nmean = N_ELEMENTS(meanrad)

  meanradm = meanrad * EXP(-2.5D * (ALOG(sigrad_fixed + 1D))^2)  ; 分散は固定して計算

  nwn = N_ELEMENTS(wnmeas)
  nangles = 181 ; Number of angles for the phase functions
  nspec = N_ELEMENTS(radmeas(0, *))

  costheta_s = cos(geometry(3,*)/180d*!pi)
  cosphi_s = cos(geometry(2,*)/180d*!pi)
  sintheta_s = sin(geometry(3,*)/180d*!pi)
  sinphi_s = sin(geometry(2,*)/180d*!pi)
  theta = costheta_s * cosphi_s / ((costheta_s^2 * cosphi_s^2 + costheta_s^2 * sinphi_s^2 + sintheta_s^2)^0.5d);内角を取って、位相角を導く
  theta = acos(theta)/!pi * 180d

  abscoeff = DBLARR(nwn, nmean)
  scacoeff = DBLARR(nwn, nmean)
  PP = DBLARR(nwn, nmean, nspec)

  ; 粒径のみを探索するループ
  FOR i = 0, nmean - 1 DO BEGIN
    IF printstage THEN PRINT, 'Radius: ' + STRING(i, "(I02)") + '/' + STRING(nmean, "(I02)")
    data = INTERPOL_AEROSOLS_OPTICAL_shoheidatabase_20241019('dust', meanradm(i), sigrad_fixed, /SAVECALC)
    abscoeff(*, i) = INTERPOL(data(*, 2), data(*, 0), wnmeas)
    scacoeff(*, i) = INTERPOL(data(*, 1), data(*, 0), wnmeas)
    A = INTERPOL_AEROSOLS_Phasefunction_ShoheiDatabase('dust', meanradm(i), sigrad_fixed, theta, wnvec = wnmeas, /SAVECALC)
    PP(*, i, *) = A
  ENDFOR
  
  radcalc = DBLARR(n_elements(radmeas(*,0)),nmean,nspec)
  
  ; abscoeffの2次元目（40要素）を逆にする
  FOR j = 0, nmean - 1 DO BEGIN
    FOR i = 0, nspec - 1 DO radcalc(*, j, i) =  scacoeff(*, j) * PP(*, j, i)
  ENDFOR
  
  radcalc = REVERSE(radcalc,1)
  radcalc = REVERSE(radcalc,2)

  indpos = INTARR(nspec)
  meanradval = DBLARR(nspec)
  sigradval = DBLARR(nspec)
  meanradmval = DBLARR(nspec)
  sigradmval = DBLARR(nspec)

  ; フィッティング結果の計算
  FOR i = 0, nspec - 1 DO BEGIN
    RMS = DBLARR(nmean)
    ;FOR j = 0, nmean - 1 DO RMS(j) = SQRT(TOTAL((radmeas(*, i) / radmeas(indref, i) - (scacoeff(*, j) * PP(*, j, i)) / (scacoeff(indref, j) * PP(indref, j, i)))^2))
    FOR j = 0, nmean - 1 DO RMS(j) = SQRT(TOTAL((radmeas(*, i) / radmeas(indref, i) - radcalc(*, j,i) / radcalc(indref,j, i))^2))
    val = MIN(RMS, ind)
    indpos(i) = ind
    meanradval(i) = meanrad(indpos(i))
    sigradval(i) = sigrad_fixed  ; 分散は固定値を使用
    meanradmval(i) = meanradm(indpos(i))
    sigradmval(i) = EXP(SQRT(ALOG(sigrad_fixed + 1D)))  ; 固定した分散で計算
    IF doprintprefit THEN PRINT, 'Spectrum ' + STRING(i, "(I02)") + '(' + STRING(geometry(0, i), "(F5.2)") + 'km): Mean = ' + STRING(meanradval(i) * 1D4, "(E8.2)") + ', Sigma = ' + STRING(sigradval(i), "(E8.2)")
  ENDFOR

  IF doplotprefit THEN BEGIN
    DEVICE, GET_SCREEN_SIZE = screen_size
    cc = COLORTABLE(39, NCOLORS = nspec + 5L, /TRANSPOSE)
    p = OBJARR(nspec, 2)
    i = 0
    p(0, 0) = PLOT(wnmeas, scacoeff(*, indpos(i)) * PP(*, indpos(i), 0) / scacoeff(indref, indpos(i)) / PP(indref, indpos(i), 0), THICK = 2, /BUFFER, $
      NAME = 'Spectrum ' + STRING(i, "(I02)") + ': Mean = ' + STRING(meanradval(i) * 1D4, "(E8.2)"), COLOR = cc(*, i))
    p(0, 1) = PLOT(wnmeas, radmeas(*, i) / radmeas(indref, i), COLOR = cc(*, i), LINESTYLE = '--', /OVERPLOT, THICK = 2)
    FOR i = 1, nspec - 1 DO BEGIN
      p(i, 0) = PLOT(wnmeas, scacoeff(*, indpos(i)) * PP(*, indpos(i), i) / scacoeff(indref, indpos(i)) / PP(indref, indpos(i), i), THICK = 2, /BUFFER, /OVERPLOT, $
        NAME = 'Spectrum ' + STRING(i, "(I02)") + ': Mean = ' + STRING(meanradval(i) * 1D4, "(E8.2)"), COLOR = cc(*, i))
      p(i, 1) = PLOT(wnmeas, radmeas(*, i) / radmeas(indref, i), COLOR = cc(*, i), /OVERPLOT, LINESTYLE = '--', THICK = 2)
    ENDFOR
    p(0, 0).SAVE, pathoutput + 'Prefit_Radius_Results.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
    p(0, 0).CLOSE
  ENDIF

  ; 高度ごとのフィッティング結果をまとめる
  nz = N_ELEMENTS(ztpn(0, *))
  meanr = DBLARR(nz)
  sigr = DBLARR(nz)
  meanrval_z = DBLARR(nz)
  sigrval_z = DBLARR(nz)
  dz = ztpn(0, 1) - ztpn(0, 0)
  FOR i = 0, nz - 1L DO BEGIN
    inds = WHERE(ABS(geometry(0, *) - ztpn(0, i)) LE dz / 2D, /NULL)
    IF N_ELEMENTS(inds) GT 0 THEN BEGIN
      meanr(i) = MEDIAN(meanradmval(inds))
      sigr(i) = MEDIAN(sigradmval(inds))
      meanrval_z(i) = MEDIAN(meanradval(inds))
      sigrval_z(i) = MEDIAN(sigradval(inds))
    ENDIF
  ENDFOR
  inds = WHERE(meanr EQ 0, /NULL)
  IF N_ELEMENTS(inds) GT 0 THEN BEGIN
    i = 0L
    WHILE meanr(i) EQ 0 DO i += 1L
    meanr(inds) = meanr(i)
    sigr(inds) = sigr(i)
    meanrval_z(inds) = meanrval_z(i)
    sigrval_z(inds) = sigrval_z(i)
  ENDIF

  RETURN, [[meanr], [sigr], [meanrval_z], [sigrval_z]]
END
