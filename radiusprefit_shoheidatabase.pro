FUNCTION RADIUSPREFIT_shoheidatabase, orbit, alt_range, ztpn, geometry, FOV, wnmeas, radmeas, $
  DOPLOTPREFIT = doplotprefit, $
  NAMEPLUSDIR = nameplusdir, $
  DOPRINTPREFIT = doprintprefit, $
  PRINTSTAGE = printstage, $
  _EXTRA = extra

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

  ;  dumb = MIN(ABS(wnmeas - 5800D), indref)
  dumb = MIN(ABS(wnmeas - 9000D), indref)

  meanrad = DINDGEN(30L) * 0.1D-4 + 0.1D-4
  sigrad = DINDGEN(10L) * 0.1D + 0.1D
  nmean = N_ELEMENTS(meanrad)
  nsig = N_ELEMENTS(sigrad)

  sigradm = DBLARR(nmean, nsig)
  meanradm = DBLARR(nmean, nsig)
  FOR i = 0, nmean - 1L DO sigradm(i, *) = EXP(SQRT(ALOG(sigrad + 1D)))
  FOR i = 0, nmean - 1L DO meanradm(i, *) = meanrad(i) * EXP(-2.5D * (ALOG(sigradm(i, *)))^2)

  nwn = N_ELEMENTS(wnmeas)
  nangles = 181 ; Number of angles for the phase functions
  nspec = N_ELEMENTS(radmeas(0, *))
  
  costheta_s = cos(geometry(3,*)/180d*!pi)
  cosphi_s = cos(geometry(2,*)/180d*!pi)
  sintheta_s = sin(geometry(3,*)/180d*!pi)
  sinphi_s = sin(geometry(2,*)/180d*!pi)
  theta = costheta_s * cosphi_s / ((costheta_s^2 * cosphi_s^2 + costheta_s^2 * sinphi_s^2 + sintheta_s^2)^0.5d);内角を取って、位相角を導く
  theta = acos(theta)/!pi * 180d

  abscoeff = DBLARR(nwn, nmean, nsig)
  scacoeff = DBLARR(nwn, nmean, nsig)
  gcoeff = DBLARR(nwn, nmean, nsig)
  PPpre = DBLARR(nwn, nmean, nsig)
  PP = DBLARR(nwn, nmean, nsig, nspec)
  FOR i = 0, nmean - 1 DO BEGIN
    FOR j = 0, nsig - 1 DO BEGIN
      IF printstage THEN PRINT, 'Radius: ' + STRING(i, "(I02)") + '/' + STRING(nmean, "(I02)") + ', Sigma: ' + STRING(j, "(I02)") + '/' + STRING(nsig, "(I02)")
      ;      data = INTERPOL_AEROSOLS_OPTICAL('dust', meanradm(i, j), sigradm(i, j), /SAVECALC)
      data = INTERPOL_AEROSOLS_OPTICAL_shoheidatabase('dust', meanradm(i, j), sigradm(i, j), /SAVECALC);shoheidatabaseに変更した。2021/12/01 Kogure Risei
      abscoeff(*, i, j) = INTERPOL(data(*, 2), data(*, 0), wnmeas)
      scacoeff(*, i, j) = INTERPOL(data(*, 1), data(*, 0), wnmeas)
      A = INTERPOL_AEROSOLS_Phasefunction_ShoheiDatabase('dust', meanradm(i, j), sigradm(i, j), theta, wnvec = wnmeas, /SAVECALC);入力：エアロゾル種、粒径、分散　　出力：その位相関数を波長、角度(0-180)ごと
      PP(*, i, j, *) = A
      ;      gcoeff(*, i, j) = INTERPOL(data(*, nangles + 7), data(*, 0), wnmeas);shoheidatabase用にコメントアウトした。2021/12/01 Kogure Risei
    ENDFOR
  ENDFOR
  

;  PP = DBLARR(nwn, nmean, nsig, nspec);上に移動した。shoheidatabase用に。
  ;FOR i = 0, nwn - 1 DO FOR j = 0, nmean - 1 DO FOR k = 0, nsig - 1 DO FOR ii = 0, nspec - 1 DO PP(i, j, k, ii) = hgfunc(COS(!DTOR * geometry(2, ii)),[gcoeff(i, j, k), 1D])

  ;  scaTg = scacoeff * gcoeff

  ;  IF doplotprefit THEN BEGIN
  ;    DEVICE, GET_SCREEN_SIZE = screen_size
  ;    cc = COLORTABLE(39, NCOLORS = MIN([256, nmean * nsig]), /TRANSPOSE)
  ;    p = OBJARR(nmean, nsig)
  ;    p(0, 0) = PLOT(wnmeas, scaTg(*, 0, 0) / scaTg(indref, 0, 0), THICK = 2, /BUFFER, NAME = 'M=' + STRING(meanrad(0) * 1D4, "(E8.2)") + 'µm, S=' + STRING(sigrad(0), "(E8.2)"), /YLOG, COLOR = cc(*, 0))
  ;    FOR i = 0, nmean - 1 DO FOR j = 0, nsig - 1 DO IF ~(i EQ 0 AND j EQ 0) THEN p(i, j) = PLOT(wnmeas, scaTg(*, i, j) / scaTg(0, i, j), /OVERPLOT, $
  ;      NAME = 'M=' + STRING(meanrad(i) * 1D4, "(E8.2)") + 'µm, S=' + STRING(sigrad(j), "(E8.2)"), COLOR = cc(*, LONG(j + nsig * i) MOD 255L), THICK = 2)
  ;    p(0, 0).SAVE, pathoutput + 'Prefit_Radius.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
  ;    p(0, 0).CLOSE
  ;  ENDIF

  indpos = INTARR(nspec, 2)
  meanradval = DBLARR(nspec)
  sigradval = DBLARR(nspec)
  meanradmval = DBLARR(nspec)
  sigradmval = DBLARR(nspec)
  FOR i = 0, nspec - 1 DO BEGIN
    RMS = DBLARR(nmean, nsig)
    FOR j = 0, nmean - 1 DO FOR k = 0, nsig - 1 DO RMS(j, k) = SQRT(TOTAL((radmeas(*, i) / radmeas(indref, i) - scacoeff(*, j, k) * PP(*, j, k, i) / scacoeff(indref, j, k) / PP(indref, j, k, i))^2))
    val = MIN(RMS, ind)
    indpos(i, *) = ARRAY_INDICES(RMS, ind)
    meanradval(i) = meanrad(indpos(i, 0))
    sigradval(i) = sigrad(indpos(i, 1))
    meanradmval(i) = meanradm(indpos(i, 0), indpos(i, 1))
    sigradmval(i) = sigradm(indpos(i, 0), indpos(i, 1))
    IF doprintprefit THEN PRINT, 'Spectrum ' + STRING(i, "(I02)") + '(' + STRING(geometry(0, i), "(F5.2)") + 'km): Mean = ' + STRING(meanradval(i) * 1D4, "(E8.2)") + ', Sigma = ' + STRING(sigradval(i), "(E8.2)")
  ENDFOR

  IF doplotprefit THEN BEGIN
    DEVICE, GET_SCREEN_SIZE = screen_size
    cc = COLORTABLE(39, NCOLORS = nspec + 5L, /TRANSPOSE)
    p = OBJARR(nspec, 2)
    i = 0
    p(0, 0) = PLOT(wnmeas, scacoeff(*, indpos(i, 0), indpos(i, 1)) * PP(*, indpos(i, 0), indpos(i, 1), 0) / scacoeff(indref, indpos(i, 0), indpos(i, 1)) / PP(indref, indpos(i, 0), indpos(i, 1), 0), THICK = 2, /BUFFER, $
      NAME = 'Spectrum ' + STRING(i, "(I02)") + ': Mean = ' + STRING(meanradval(i) * 1D4, "(E8.2)") + ', Sigma = ' + STRING(sigradval(i), "(E8.2)"), COLOR = cc(*, i))
    p(0, 1) = PLOT(wnmeas, radmeas(*, i) / radmeas(indref, i), COLOR = cc(*, i), LINESTYLE = '--', /OVERPLOT, THICK = 2)
    FOR i = 1, nspec - 1 DO BEGIN
      p(i, 0) = PLOT(wnmeas, scacoeff(*, indpos(i, 0), indpos(i, 1)) * PP(*, indpos(i, 0), indpos(i, 1), i) / scacoeff(indref, indpos(i, 0), indpos(i, 1)) / PP(indref, indpos(i, 0), indpos(i, 1), i), THICK = 2, /BUFFER, /OVERPLOT, $
        NAME = 'Spectrum ' + STRING(i, "(I02)") + ': Mean = ' + STRING(meanradval(i) * 1D4, "(E8.2)") + ', Sigma = ' + STRING(sigradval(i), "(E8.2)"), COLOR = cc(*, i))
      p(i, 1) = PLOT(wnmeas, radmeas(*, i) / radmeas(indref, i), COLOR = cc(*, i), /OVERPLOT, LINESTYLE = '--', THICK = 2)
    ENDFOR
    p(0, 0).SAVE, pathoutput + 'Prefit_Radius_Results.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
    p(0, 0).CLOSE
  ENDIF

  nz = N_ELEMENTS(ztpn(0, *))
  meanr = DBLARR(nz)
  sigr = DBLARR(nz)
  meanrval = DBLARR(nz)
  sigrval = DBLARR(nz)
  dz = ztpn(0, 1) - ztpn(0, 0)
  FOR i = 0, nz - 1L DO BEGIN
    inds = WHERE(ABS(geometry(0, *) - ztpn(0, i)) LE dz / 2D, /NULL)
    IF N_ELEMENTS(inds) GT 0 THEN BEGIN
      meanr(i) = MEDIAN(meanradmval(inds))
      sigr(i) = MEDIAN(sigradmval(inds))
      meanrval(i) = MEDIAN(meanradval(inds))
      sigrval(i) = MEDIAN(sigradval(inds))
    ENDIF
  ENDFOR
  inds = WHERE(meanr EQ 0, /NULL)
  IF N_ELEMENTS(inds) GT 0 THEN BEGIN
    i = 0L
    WHILE meanr(i) EQ 0 DO i += 1L
    meanr(inds) = meanr(i)
    sigr(inds) = sigr(i)
    meanrval(inds) = meanrval(i)
    sigrval(inds) = sigrval(i)
  ENDIF

  RETURN, [[meanr], [sigr], [meanrval], [sigrval]]
END