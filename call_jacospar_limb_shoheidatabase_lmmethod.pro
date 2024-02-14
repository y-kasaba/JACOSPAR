FUNCTION CALL_JACOSPAR_LIMB_shoheidatabase_LMmethod, orbit, $
  alt_range, $

  CHECKJACOBIANJACOSPAR = checkjacobianjacospar, $
  CHECKJACOBIANK= checkjacobiank, $
  CHECKJACOBIANRODGER = checkjacobianrodger, $
  CLEAN = clean, $
  DORADIUSPREFIT = doradiusprefit, $
  DWN = dwn, $
  FACTDUST = factdust, $
  FACTORINI = factorini, $
  FACTORIRADIUSM = factoriradiusm, $
  FACTORIRADIUSS = factoriradiuss, $
  FORCE = force, $
  GEOMETRY = geometry, $
  GETMEASUREDSPECTRA = getmeasuredspectra, $
  LOADZTPNFROMINAF = loadztpnfrominaf, $
  NODELETE = nodelete, $
  NOFIT = nofit, $
  NOMEASINTERP = nomeasinterp, $
  NOMOD = nomod, $
  NOREGION = noregion, $
  ONLYCOMPUTERAD = onlycomputerad, $
  SETALBEDO = setalbedo, $
  SETNOISE = setnoise, $
  SETPROFILEDUST = setprofiledust, $
  SETPROFILERADIUS = setprofileradius, $
  STARTFROM = startfrom, $
  STARTRADIUS = StartRadius, $
  STARTVARIANCE = StartVariance, $
  SUBSTRACTLASTSPECTRUM = substractlastspectrum, $
  THRESHOLD = threshold, $
  USESYNTHETIC = usesynthetic, $
  WN_MEASURED = wn_measured, $
  Z = z, $

  _EXTRA = extra

  !EXCEPT = 0

  T0=systime(1)
  ;;;;;;  kogure 2021/4/14
  ;;;;;;  with DOSORT input file(dust&ice's phase function,cross section.)
  ; -----------------------
  ; Setting up the keywords
  ; -----------------------
  IF ~KEYWORD_SET(checkjacobianjacospar) THEN checkjacobianjacospar = 0
  IF ~KEYWORD_SET(checkjacobiank) THEN checkjacobiank = 0
  IF ~KEYWORD_SET(checkjacobianrodger) THEN checkjacobianrodger = 0
  IF ~KEYWORD_SET(clean) THEN clean = 0
  IF ~KEYWORD_SET(factorini) THEN BEGIN
    factorini = DBLARR(5)
    factorini(*) = 1D
  ENDIF
  IF ~KEYWORD_SET(factoriradiusm) THEN BEGIN
    factoriradiusm = DBLARR(2)
    factoriradiusm(*) = 1D
  ENDIF
  IF ~KEYWORD_SET(factoriradiuss) THEN BEGIN
    factoriradiuss = DBLARR(5)
    factoriradiuss(*) = 1D
  ENDIF
  IF ~KEYWORD_SET(force) THEN force = 0
  IF ~KEYWORD_SET(getmeasuredspectra) THEN getmeasuredspectra = 0
  IF ~KEYWORD_SET(nodelete) THEN nodelete = 0
  IF ~KEYWORD_SET(nofit) THEN nofit = 0
  IF ~KEYWORD_SET(nomeasinterp) THEN nomeasinterp = 0L
  IF ~KEYWORD_SET(nomod) THEN nomod = 0
  IF ~KEYWORD_SET(doradiusprefit) THEN doradiusprefit = 0
  IF ~KEYWORD_SET(onlycomputerad) THEN onlycomputerad = 0
  IF ~KEYWORD_SET(StartRadius) THEN StartRadius = [1.5D-4, 2D-4]
  IF ~KEYWORD_SET(StartVariance) THEN StartVariance = [.3D, .3D]
  IF ~KEYWORD_SET(substractlastspectrum) THEN substractlastspectrum = 0
  IF ~KEYWORD_SET(usesynthetic) THEN usesynthetic = 0

  fitalbedo = 0L
  fitradius = 0L
  fitsigma = 0L
  makeplot = 0L
  makeprint = 0L
  onlyfitspecies = INDGEN(5)
  IF KEYWORD_SET(extra) THEN BEGIN
    names = TAG_NAMES(extra)
    FOR i = 0, N_ELEMENTS(names) - 1 DO BEGIN
      CASE STRLOWCASE(names(i)) OF
        'dz': dz = extra.dz
        'fitalbedo': fitalbedo = 1L
        'fitsigma': fitsigma = 1L
        'makeplot': makeplot = extra.makeplot
        'makeprint': makeprint = extra.makeprint
        'nameplusdir': nameplusdir = extra.nameplusdir
        'fitradius': fitradius = 1L
        'onlyfitspecies': onlyfitspecies = extra.onlyfitspecies
        'rm': Rm = extra.Rm
        'solar': solar = extra.solar
        'wnrange': wnrange = extra.wnrange
        ELSE:
      ENDCASE
    ENDFOR
  ENDIF
  IF N_ELEMENTS(WHERE(onlyfitspecies GT 4, /NULL)) GT 0 OR N_ELEMENTS(onlyfitspecies) EQ 0 THEN BEGIN
    PRINT, "The indexes of the fitted species exceeds the total number of species or no fitted species. Aborting."
    RETURN, -1
  ENDIF

  IF makeplot THEN DEVICE, GET_SCREEN_SIZE = screen_size;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;スクリーンさいずはどこで定義している？

  ; ------------------------
  ; Setting up the save path
  ; ------------------------
  itermax = 60L; plotする際は2以上に設定
  pathnames = PATH_MANAGEMENT()
  path_solar = pathnames.solar ; path of the solar radiance
  path_output = pathnames.results + orbit + '/'
  IF KEYWORD_SET(nameplusdir) THEN path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/" ELSE nameplusdir = ""
  IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
  IF clean THEN BEGIN
    RESULT = DIALOG_MESSAGE("Are you sure you want to clean the directory?", /CENTER, /QUESTION, /CANCEL)
    IF RESULT EQ "No" THEN clean = 0L
    IF RESULT EQ "Cancel" THEN RETURN, -1
    IF clean THEN BEGIN
      IF FILE_TEST(path_output + "/Figures") THEN BEGIN
        liste = FILE_SEARCH(path_output + '/Figures/*.*')
        IF N_ELEMENTS(liste) GT 1 THEN FOR i = 0, N_ELEMENTS(liste) - 1 DO FILE_DELETE, liste(i), /QUIET
      ENDIF
      liste = FILE_SEARCH(path_output + '*.*')
      IF N_ELEMENTS(liste) GT 1 THEN FOR i = 0, N_ELEMENTS(liste) - 1 DO FILE_DELETE, liste(i), /QUIET
    ENDIF
  ENDIF
  IF makeplot THEN IF ~FILE_TEST(path_output + '/Figures/', /DIRECTORY) THEN FILE_MKDIR, path_output + '/Figures/'
  IF ~clean AND ~nodelete THEN BEGIN
    listtodelete = FILE_SEARCH(path_output + '*.res', /TEST_ZERO_LENGTH)
    anydeleted = 0
    FOR i = 0, N_ELEMENTS(listtodelete) - 1 DO BEGIN
      IF STRLEN(listtodelete(i)) THEN BEGIN
        FILE_DELETE, listtodelete(i), /ALLOW_NONEXISTENT, /NOEXPAND_PATH, /QUIET
        FILE_DELETE, STRMID(listtodelete(i), 0, STRLEN(listtodelete(i)) - 4) + '.log', /ALLOW_NONEXISTENT, /NOEXPAND_PATH, /QUIET
        anydeleted = 1
      ENDIF
    ENDFOR
    IF anydeleted THEN PRINT, "Some .res and .log files have been deleted because they were empty."
  ENDIF

  ; ----------------------------------
  ; Check if run has been already done
  ; ----------------------------------
  IF FILE_TEST(path_output + 'Results.sav') AND ~force THEN BEGIN
    RESTORE, path_output + 'Results.sav'
    converged = iter LT itermax - 1
    PRINT, 'This orbit has already been computed. Convergence = ' + STRING(converged, "(I1)")
    RETURN, converged
  ENDIF

  ; --------------------------
  ; Setting up some parameters
  ; --------------------------
  gas_species = ['co2', 'co', 'h2o']
  aerosols_species = ['dust', 'ice']
  gas_species_all = gas_species
  aerosols_species_all = aerosols_species
  naerosols_max = N_ELEMENTS(aerosols_species)
  nspecies = N_ELEMENTS(onlyfitspecies)
  IF N_ELEMENTS(onlyfitspecies(WHERE(onlyfitspecies LE 2, /NULL))) GT 0 THEN $
    gas_species = gas_species(onlyfitspecies(WHERE(onlyfitspecies LE 2, /NULL))) $
  ELSE gas_species = !NULL
  IF N_ELEMENTS(onlyfitspecies(WHERE(onlyfitspecies - 3 GE 0, /NULL))) GT 0 THEN $
    aerosols_species = aerosols_species(onlyfitspecies(WHERE(onlyfitspecies - 3 GE 0, /NULL)) - 3) $
  ELSE aerosols_species = !NULL
  ngas = N_ELEMENTS(gas_species)
  naerosols = N_ELEMENTS(aerosols_species)

  kB = 1.38D-16 ; [cgs] Boltzmann's constant[erg/K]([m^2*kg/s^2/K])
  gM = 371.1D ; [cgs]表面重力[cm/s^2]
  Nav = 6.02D23;アボガドロ定数[/mol]
  MM_CO2 = 44.01 ; [g/mol]
                                                                                                                                                       ;                                         H2O -         CO -            CO2 -           H2O -           DUST -          H2O -          CO2 -           H2O -           DUST
  IF ~KEYWORD_SET(wnrange) THEN wnrange = [[2660D,2740D], [2820D,2880D],[2880D,3040D],[3185D,3245D],[3390D,3450D],[3885D,3945D],[4480D,4540D],[5570D,5630D],[6570D,6630D],[7750D,7810D]];[; [[3780D, 4020D],[4020D, 4500D], [4500D, 5180D], [5180D, 5800D], [5800D, 6600D], [6600D, 7700D],[7700D, 8500D], [8500D, 9000D], [9000D, 10000D]]
  IF ~KEYWORD_SET(threshold) THEN threshold = 1D-4;  20210930 kogure change to 1d-3
  resol = 31.2461D
  IF ~KEYWORD_SET(dwn) THEN dwn = 1D-2
  nu_ref = 8500D

  ; -------------
  ; Read geometry
  ; -------------
  geometry = READ_GEOMETRY(orbit)
  IF N_ELEMENTS(geometry) EQ 0 THEN BEGIN
    PRINT, "Orbit is not found. Aborting."
    RETURN, 0
  ENDIF
  sequence = WHERE((geometry(0, *) GE alt_range(0)) AND (geometry(0, *) LE alt_range(1)) AND (geometry(0, *) GT 0D), /NULL)
  IF N_ELEMENTS(sequence) EQ 0 THEN BEGIN
    PRINT, "There are no measurement in the selected altitude region. Aborting."
    RETURN, 0
  ENDIF
  geometry = geometry(*, sequence)
  nspectra = N_ELEMENTS(geometry(0, *))
  IF nspectra LE 1 THEN BEGIN
    print, "Not enough spectra in this orbit."
    RETURN, 0
  ENDIF

  FOV = 1.2D-3 / !DPI * 180D

  ; --------------------------------------------------------
  ; Reshape altitude range depending on the scanned geometry
  ; --------------------------------------------------------
  IF KEYWORD_SET(dz) THEN z = [DINDGEN(CEIL(geometry(0, nspectra - 1) / dz) + 1) * dz] ELSE z = [[0D], [DOUBLE(geometry(0, *))]]
  ;z = [DINDGEN(CEIL(geometry(0, nspectra - 1) / dz) + 1) * dz]; add kogure ,because sometimes upper IF sentense cause error.
  nlayer = N_ELEMENTS(z)
  ; -----------------------------------
  ; Defining the radii of the particles
  ; -----------------------------------
  Radiusm = DBLARR(naerosols_max, nlayer, itermax) ; Radius mean modal value
  Radiuss = DBLARR(naerosols_max, nlayer, itermax) ; Radius variance modal value
  IF KEYWORD_SET(setprofileradius) THEN BEGIN
    Radiusm(*, *, 0) = radiusM_mod;;;;;;;;;;;;;;;;;どこでradiusM_modを定義？？
    Radiuss(*, *, 0) = radiusS_mod
    Radiusfinalm(*, *, 0) = StartRadius
    Radiusfinals(*, *, 0) = StartVariance
  ENDIF ELSE BEGIN
    IF N_ELEMENTS(StartRadius(0, *)) EQ 1 THEN FOR i = 0, nlayer - 1 DO Radiusm(*, i, 0) = StartRadius ELSE Radiusm(*, *, 0) = radiusM_mod ; [cm]
    IF N_ELEMENTS(StartRadius(0, *)) EQ 1 THEN FOR i = 0, nlayer - 1 DO Radiuss(*, i, 0) = StartVariance ELSE Radiuss(*, *, 0) = radiusS_mod
  ENDELSE
;  Radius_ap = read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Apriori/apriori_4_2_2_GEMZ_wz_LOGND/Dust_Radius_Risei_1d-4.dat', data_start = 1)  & Radius_ap = Radius_ap.field1;;;;;;;2022/1/1 kogure risei 線形に粒径が減少するプロファイル
;  Radiusm(0,*,0) = interpol(Radius_ap(1,*), Radius_ap(0,*), z(0,*));;;;;;;2022/1/1 kogure risei 線形に粒径が減少するプロファイル

  ; -------------------
  ; Calculating the ILS
  ; -------------------
  wn_ILS = -50D + DINDGEN(ROUND(100 / dwn) + 1) * dwn
  ILS = GAUSSIAN(wn_ILS, 1D, resol, 0D)
  ILS /= INTEGRATE(wn_ILS, ILS)

  ; ---------------------
  ; Read measured spectra
  ; ---------------------
  nwindows = N_ELEMENTS(wnrange(0, *))
  IF usesynthetic THEN BEGIN
    rad_measured = READ_BINARY(pathnames.results + orbit + '_Synth/Spectra.dat', DATA_TYPE = 5)
    rad_measured = REFORM(rad_measured(2 : N_ELEMENTS(rad_measured) - 1), rad_measured(0), rad_measured(1) + 1L)
    sequence = INDGEN(nspectra)
    wn_measured = rad_measured(*, 0)
    rad_measured = rad_measured(*, 1 : nspectra)
    drad_measured = rad_measured
    drad_measured(*) = .01D
  ENDIF ELSE BEGIN
    rad_measured = READ_SPECTRA(orbit, WN = wn_measured, ERROR = drad_measured, _EXTRA = extra)
    IF substractlastspectrum THEN FOR i = 0, N_ELEMENTS(rad_measured(0, *)) - 1L DO rad_measured(*, i) -= rad_measured(*, N_ELEMENTS(rad_measured(0, *)) - 1L)
    IF N_ELEMENTS(rad_measured) EQ 0 THEN RETURN, !NULL
  ENDELSE
  inds = INTARR(N_ELEMENTS(wn_measured))
  FOR i = 0, nwindows - 1 DO inds = inds OR (wn_measured GE wnrange(0, i) AND wn_measured LE wnrange(1, i))
  inds = WHERE(inds);, /NULL)
  rad_measured = rad_measured(*, sequence)
  drad_measured = drad_measured(*, sequence)
  wnall = wn_measured
  radall = rad_measured
  wn_measured = wn_measured(inds)
  rad_measured = rad_measured(inds, *)
  drad_measured = drad_measured(inds, *)
  nspec = N_ELEMENTS(rad_measured(0, *))
  nwn = N_ELEMENTS(wn_measured)
  lambda_measured = 1D4 / wn_measured

  IF getmeasuredspectra THEN RETURN, [[wn_measured], [rad_measured]]

  IF makeplot THEN BEGIN
    cc = COLORTABLE(39, NCOLORS = nspec + 1, /TRANSPOSE)
    p = OBJARR(nspec)
    FOR i = 0, nspec - 1 DO  BEGIN
      p(i) = PLOT(lambda_measured, rad_measured(*, i), /OVERPLOT, XTITLE = 'Wavelength ($\lambda [\mu m]$)', YTITLE = 'Radiance [CGS]', COLOR = cc(*, i), $
        NAME = " - Alt: " + STRING(geometry(0, i), "(F6.3)") + "km", DIMENSIONS = ROUND(screen_size * .9), TITLE = orbit, /BUFFER, THICK = 3)
      k = PLOT(lambda_measured, rad_measured(*, i) + drad_measured(*, i), COLOR = cc(*, i), /OVERPLOT)
      k = PLOT(lambda_measured, rad_measured(*, i) - drad_measured(*, i), COLOR = cc(*, i), /OVERPLOT)
    ENDFOR
    L = LEGEND(TARGET = p)
    p(0).SAVE, path_output + 'Figures/Measured_spectra.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
    p(0).CLOSE
  ENDIF

  ; ---------------------------
  ; Get the a-priori atmosphere
  ; ---------------------------
  X = DBLARR(nspecies, nlayer, itermax)

  ztpn_ini = GETZTPN_LimbTest_new(MEAN(geometry(4, *)), MEAN(geometry(5, *)), MEAN(geometry(6, *)), MEAN(geometry(7, *)), z, cov, StartRadius);T[K],P[pa],CO2ndens,,
  openw,lun,path_output+'ztpn',/get_lun;;;;kogure 20200505
  printf,lun,ztpn_ini
  free_lun,lun;;;;;;;;;;;;;;;;;kogure 20200505
  IF KEYWORD_SET(startfrom) THEN BEGIN
    IF ~FILE_TEST(pathnames.results + startfrom + '/ResultsZTPN.sav') THEN BEGIN
      PRINT, 'Run ' + startfrom + ' not found. Aborting.'
      RETURN, -1
    ENDIF
    onlyfitspeciessave = onlyfitspecies
    Radiusmsave = Radiusm
    Radiusssave = Radiuss
    fitradiussave = fitradius
    fitsigmasave = fitsigma
    RESTORE, pathnames.results + startfrom + '/ResultsZTPN.sav'
    ztpn_ini(5 + onlyfitspecies, *) = ztpn(5 + onlyfitspecies, *, iter);;;;;;;;;;;;;;;;CO2,CO,H2O,H2Oice数密度
    Radiusmsave(*, *, 0) = Radiusm(*, *, iter)
    Radiusm = Radiusmsave
    Radiusssave(*, *, 0) = Radiuss(*, *, iter)
    Radiuss = Radiusssave
    onlyfitspecies = onlyfitspeciessave
    fitradius = fitradiussave
    fitsigma = fitsigmasave
  ENDIF
  Ht = kb * Nav * ztpn_ini(1, *) / gM / MM_CO2 * 1D-5 ; [km]スケールハイトHt=RT/Mg
  corrlength = MEAN(Ht)
  FOR i = 0, 4 DO ztpn_ini(5 + i, *) *= factorini(i)
  ztpn = DBLARR(N_ELEMENTS(ztpn_ini(*, 0)), N_ELEMENTS(ztpn_ini(0, *)), itermax)
  ztpn(*, *, 0) = ztpn_ini
  FOR i = 0, 1 DO Radiusm(i, *, 0) *= factoriradiusm(i)
  FOR i = 0, 1 DO Radiuss(i, *, 0) *= factoriradiuss(i)
  IF KEYWORD_SET(factdust) THEN ztpn(8, *, 0) *= [factdust(0), factdust]
  IF KEYWORD_SET(setprofiledust) THEN ztpn(8, *, 0) = setprofiledust

  ; ------------------------------------
  ; Calculate the apriori surface albedo
  ; ------------------------------------
  ;  albedoFile = READ_ASCII(pathnames.albedo + 'true_albedo_orb_0937_4_p_63_scan_482_ARS.dat'); ;2020/2/13 to set albedo = 0.15f
  ;  albedoFile = albedoFile.FIELD1
  S_albedo = DBLARR(nwn, itermax)
  wn_albedo = wn_measured
  ;  IF KEYWORD_SET(setalbedo) THEN BEGIN   ;2020/2/13 to set albedo = 0.15
  ;    CASE N_ELEMENTS(setalbedo) OF
  ;      1: S_albedo(*, 0) = setalbedo
  ;      2: S_albedo(*, 0) = INTERPOL(setalbedo, [wn_albedo(0), wn_albedo(N_ELEMENTS(wn_albedo) - 1L)], wn_albedo)
  ;      ELSE: S_albedo(*, 0) = INTERPOL(setalbedo(*, 1), setalbedo(*, 0), wn_albedo)
  ;    ENDCASE
  ;  ENDIF ELSE S_albedo(*, 0) = INTERPOL(albedoFile(1, *), albedoFile(0, *), wn_albedo)
  S_albedo(*,0) = 0.10d;makearr(n_elements(S_albedo(*,0)),Value=0.1500000)
  ; --------------------
  ; Prefit of the radius
  ; --------------------
  IF doradiusprefit THEN BEGIN
    resultsprefit = RADIUSPREFIT_shoheidatabase(orbit, alt_range, ztpn(*, *, 0), geometry, FOV, wnall, radall, _EXTRA = extra)
    Radiusm(0, *, 0) = resultsprefit(*, 0)
;    Radiuss(0, *, 0) = resultsprefit(*, 1)
  ENDIF

  ; -----------------------
  ; Set up the VMR profiles
  ; -----------------------
  X(*, *, 0) = GET_X(ztpn(*, *, 0), nlayer, onlyfitspecies)
  IF makeplot THEN p = PLOT_PROFILES(ztpn(*, *, 0), gas_species, aerosols_species, onlyfitspecies, 0, itermax, path_output, Radiusm(*, *, 0), Radiuss(*, *, 0), fitradius, fitsigma, wn_measured, S_albedo, _EXTRA = extra)
  openw,lun,path_output+'VMR',/get_lun
  printf,lun,X(*,*,0)
  free_lun,lun
  ; -----------------------
  ; Calculate the a-prioris
  ; -----------------------
  radtemp = CREATE_SYNTHETIC_SPECTRA_Limb_shoheiDatabase(orbit, alt_range, wn_apriori, z, ztpn(*, *, 0), Ht, geometry, dwn, FOV, Radiusm(*, *, 0), $
    Radiuss(*, *, 0), wn_albedo, S_albedo(*, 0), jactemp, ILS = ILS, NAMEPLUSFILE = 'step00_', $
    TAU_ABS = tau_abs, TAU_AERO_SCA = tau_aero_sca, RADNOCONV = radtempnoconv, JACNOCONV = jactempnoconv, DRDABS = dRdabs, DRDSCA = dRdsca, _EXTRA = extra)

  lambda_apriori = 1D4 / wn_apriori;;;;;;;;;;;;;波長に直している

  IF nomeasinterp OR checkjacobiank THEN BEGIN
    rad = radtemp
    jac = jactemp
    nwn = N_ELEMENTS(wn_apriori)
  ENDIF ELSE BEGIN
    rad = DBLARR(nwn, nspectra, itermax)
    jac = DBLARR(nwn, nlayer, nspectra, 1L + nspecies + naerosols * 2L, itermax)
    FOR i = 0, nspectra - 1 DO BEGIN
;      rad(*, i, 0) = MEANINTERP(wn_apriori, radtemp(*, i), wn_measured)
;      FOR j = 0, nlayer - 1L DO FOR k = 0,  1L + nspecies + naerosols * 2L - 1 DO jac(*, j, i, k, 0) = MEANINTERP(wn_apriori, jactemp(*, j, i, k), wn_measured)
      rad(*, i, 0) = INTERPOL(radtemp(*, i), wn_apriori, wn_measured);;;;;;/////kogure risei 2021/12/14 MEANINTERPという関数を使うことによって波長が飛んだ時に（例えばガス吸収波長を抜かした時に）、飛ばした波長の両端のきど値がおかしくなるので、INTERPOLに修正した。
      FOR j = 0, nlayer - 1L DO FOR k = 0,  1L + nspecies + naerosols * 2L - 1 DO jac(*, j, i, k, 0) = INTERPOL(jactemp(*, j, i, k), wn_apriori, wn_measured);;;;;/////kogure risei 2021/12/14
    ENDFOR
  ENDELSE

  IF makeplot THEN BEGIN
    PathSave = path_output + 'Figures/'
    IF ~FILE_TEST(PathSave, /DIRECTORY) THEN FILE_MKDIR, PathSave
    cc = COLORTABLE(39, NCOLORS = nspec + 1, /TRANSPOSE)
    p = OBJARR(nspec)
    FOR i = 0, nspec - 1 DO p(i) = PLOT(lambda_apriori, radtemp(*, i), /OVERPLOT, XTITLE = 'Wavelength ($\lambda [\mu m]$)', YTITLE = 'Radiance [CGS]', COLOR = cc(*, i), $
      NAME = " - Alt: " + STRING(geometry(0, i), "(F6.3)") + "km", DIMENSIONS = ROUND(screen_size * .9), TITLE = orbit, /BUFFER)
    L = LEGEND(TARGET = p)
    p(0).SAVE, path_output + 'Figures/Apriori_spectra_inf.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
    p(0).CLOSE
    p = OBJARR(nspec)
    FOR i = 0, nspec - 1 DO p(i) = PLOT(lambda_measured, rad(*, i, 0), /OVERPLOT, XTITLE = 'Wavelength ($\lambda [\mu m]$)', YTITLE = 'Radiance [CGS]', COLOR = cc(*, i), $
      NAME = " - Alt: " + STRING(geometry(0, i), "(F6.3)") + "km", DIMENSIONS = ROUND(screen_size * .9), TITLE = orbit, /BUFFER)
    L = LEGEND(TARGET = p)
    p(0).SAVE, path_output + 'Figures/Apriori_spectra.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
    p(0).CLOSE
  ENDIF

  ; 2020/2/13 followings is same as line 369
  ;  IF makeplot THEN BEGIN
  ;    p = OBJARR(nspec)
  ;    FOR i = 0, nspec - 1 DO p(i) = PLOT(lambda_apriori, radtemp(*, i), /OVERPLOT, XTITLE = 'Wavelength ($\lambda [\mu m]$)', YTITLE = 'Radiance [CGS]', COLOR = cc(*, i), $
  ;      NAME = " - Alt: " + STRING(geometry(0, i), "(F6.3)") + "km", DIMENSIONS = ROUND(screen_size * .9), TITLE = orbit, /BUFFER)
  ;    L = LEGEND(TARGET = p)
  ;    p(0).SAVE, path_output + 'Figures/Apriori_spectra_inf.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
  ;    p(0).CLOSE
  ;  ENDIF

  IF onlycomputerad THEN RETURN, rad(*, *, 0)
  filename='rad'
  data=dblarr(2,nwn)
  data(0,*)=lambda_measured
  data(1,*)=rad(*,0,0)
  openw,lun,path_output+filename,/get_lun
  printf,lun,data
  close,lun
  free_lun,lun

  ; ------------------------------------------
  ; Select the wavenumber regions to be fitted
  ; ------------------------------------------
  regions = INTARR(nwn)
  regions(*)= 1L
  IF KEYWORD_SET(noregion) THEN FOR i = 0, N_ELEMENTS(noregion(0, *)) - 1 DO regions(WHERE(wn_measured GT noregion(0, i) AND wn_measured LT noregion(1, i), /NULL)) = 0
  regions = WHERE(regions, /NULL)
  nregions = N_ELEMENTS(regions)
  IF regions EQ !NULL THEN RETURN, -1

  T1=systime(1)
  print,'Time :',T1-T0
  ; ----------------------------------------------
  ; PLOT THE MEASURED, A-PRIORI AND FITTED SPECTRA
  ; ----------------------------------------------
  IF makeplot THEN BEGIN
    RMS = DBLARR(nspectra, itermax)
    RMS(*, 0) = PLOT_SPECTRA_JACOSPAR(wn_measured, rad_measured, rad, drad_measured, regions, nspectra, geometry(0, *), sequence, path_output, 0)
  ENDIF

  IF checkjacobianjacospar THEN BEGIN
    CHECKJACOBIANJACOSPARFUNCTION_shoheidatabase, wn_apriori, nlayer, dwn, orbit, alt_range, z, ztpn(*, *, 0), geometry, dwn, FOV, Radiusm, Radiuss, wn_albedo, S_albedo(*, 0), solar, ILS, $
      radtemp, radtempnoconv, dRdabs, dRdsca, tau_abs, tau_aero_sca, path_output + 'Figures/', extra
    RETURN, -1
  ENDIF

  IF checkjacobianrodger THEN BEGIN
    CHECKJACOBIANRODGERFUNCTION, wn_apriori, nlayer, nspectra, dwn, orbit, alt_range, z, ztpn(*, *, 0), geometry, dwn, FOV, Radiusm, Radiuss, wn_albedo, S_albedo(*, 0), solar, ILS, $
      radtempnoconv, jactempnoconv, dRadnoconv, djactempnoconv, path_output + 'Figures/', extra
    RETURN, -1
  ENDIF

  ;====== kogure 2021/5/24 ====
  radiance=rad(*,*,0)
  wvl=lambda_measured
  alt=geometry(0, *)
  save,radiance,wvl,alt,filename='/Volumes/HD-PNFU3/JACOSPAR/Inversion/OMEGA/Results/' + orbit + NAMEPLUSDIR + '/Rad_orb'+orbit+'.sav'
  ;=======================

;  ;;----------------------- 20200818 kogure to stop to just calculate forward ------------
;  RESULT = DIALOG_MESSAGE("Do you want to do iteration?", /CENTER, /QUESTION)
;  IF RESULT EQ "No" THEN stop
;  ;;--------------------------------------------------------------------------------------

  ; ----------------------------------------
  ; Calculate the Rodger's a-priori matrices
  ; ----------------------------------------
;  drad_measured(where(drad_measured le 1d-3)) = 1d-3;;;;;;;;kogure 202107921
  Y = GET_Y(rad_measured, nspectra, regions);２次元配列の観測輝度[波長数、スペクトル数]を１次元配列Yに直す
  invSe = GET_INVSE(drad_measured, nspectra, regions) ;観測輝度の誤差の逆数についても同様に1時現に直す
  xa = GET_XVECTOR(X(*, *, 0), nlayer, onlyfitspecies, fitradius, fitsigma, Radiusm(*, *, 0), Radiuss(*, *, 0), S_albedo(*, 0), fitalbedo)  ;X:大気情報[リトリーバルターゲット、1kmごと大気層]　xa:リトリーバル大気ターゲットアプリオリを大気層で１列に並べ替えたもの。dustとiceなら２種類*粒径2種類*大気層nの要素数で、0:n-1までダスト、n:2n-1までアイス、と入っている。値はlog(molec/cm3)とlog(粒径)
  invSa = GET_INVSA(X(*, *, 0), z, gas_species_all, aerosols_species_all, onlyfitspecies, corrlength, fitradius, fitsigma, Radiusm(*, *, 0), Radiuss(*, *, 0), cov, S_albedo(*, 0), wn_albedo, MEAN(wn_measured(1 : nwn - 1L) - wn_measured(0 : nwn - 2L)), _EXTRA = extra);アプリオリの分散・この中でSaを変えている。
;  invSa = invSa*0.01d;;;;;;;;   一時的にテスト。パラメータ誤差を大きくする必要がある（北海道情報大、佐藤さん　20210907）
  ;invSa(where(invSa ne 0.0)) = 0.01;;;;;;;;;   一時的にテスト。パラメータ誤差を大きくする必要がある（北海道情報大、佐藤さん　20210907）
  Kmat = DBLARR(nregions * nspectra, nwn * fitalbedo + nlayer * (ngas + (1L + fitradius + fitsigma) * naerosols), itermax);ヤコビアンの値（大気数＊波長数、ターゲット数＊大気数、イタレーション回数）。入り方はxk,Fxkと同じ
  Kmat(*, *, 0) = GET_K(jac(*, *, *, *, 0), ztpn(*, *, 0), Radiusm(*, *, 0), Radiuss(*, *, 0), S_albedo(*, 0), onlyfitspecies, regions, wn_measured, fitradius, fitsigma, nlayer, nspectra, fitalbedo);ヤコビアンを入れる。ちなみにヤコビアンはcreate_syntheticの250行目あたりで計算されている
  Fxk = DBLARR(itermax, nregions * nspectra);出力輝度記録用配列
  Fxk(0, *) = GET_Y(rad(*, *, 0), nspectra, regions);最初のforward計算輝度の1次元配列。0:n-1まで大気層１の波長region個の輝度、n:2n-1まで大気層2の波長region個の輝度

  IF checkjacobiank THEN BEGIN
    CHECKJACOBIANKFUNCTION, wn_apriori, nlayer, nspectra, orbit, alt_range, z, ztpn(*, *, 0), geometry, dwn, FOV, Radiusm, Radiuss, solar, ILS, $
      radtemp, jactemp(*, *, *, *, 0), regions, Kmat(*, *, 0), xa, path_output + 'Figures/', ngas, naerosols, onlyfitspecies, extra
    RETURN, -1
  ENDIF
  print,string(7B)
  ; -----------------------------
  ; COMPUTE THE AVERAGING KERNELS
  ; -----------------------------
  IF makeprint THEN PRINT, "DOF for a_priori:"
  A = GET_A(z, Kmat(*, *, 0), invSa, invSe, nlayer, nwn, gas_species, aerosols_species, onlyfitspecies, 0, path_output, S, invSaKinvSeK, fitradius, fitsigma, _EXTRA = extra);Averaging kernels
;  IF ~FINITE(TOTAL(TOTAL(A,2))) THEN BEGIN
;    PRINT, "The A matrix contains NaNs."
;    RETURN, -1
;  ENDIF
  IF makeprint THEN PRINT, ""

  ; -----------------------------
  ; START THE RODGERS CALCULATION
  ; -----------------------------
  convX = DBLARR(itermax, 1)
  convX(0) = 1D6
  epsilonX = threshold * (nspecies + (fitradius + fitsigma) * naerosols) * nlayer
  xk = DBLARR(itermax, N_ELEMENTS(xa));リトリーバルされる変数x
  xk(0, *) = xa
  iter = 1L
  v = 1d1;20211101 KogureRisei. Parameter of Levenberg - Marquardt Method
  ratio = 1d;20211101 KogureRisei. Parameter of Levenberg - Marquardt Method
  kai2 = DBLARR(itermax, 1);20211119 kogure risei 評価関数のため加えた
  kai2(0) = 1d
  
  SAVE, Kmat, S, invSa, invSe, A, X, Y, Fxk, xk, rad, jac, xa, convX, wn_measured, rad_measured, drad_measured, nlayer, nwn, nspectra, alt_range, wnrange, ztpn, onlyfitspecies, S_albedo, RMS, iter, fitradius, fitsigma, FILENAME = path_output + 'Results.sav';add this line. 2021/10/9 risei kogure
  SAVE, ztpn, Radiusm, Radiuss, S_albedo, iter, onlyfitspecies, fitradius, fitsigma, FILENAME = path_output + 'ResultsZTPN.sav';add this line. 2021/10/9 risei kogure
  
;  WHILE (iter LT itermax) AND (convX(iter - 1) GE epsilonX) AND ~nofit DO BEGIN
  WHILE (iter LT itermax) AND (convX(iter - 1) GE epsilonX) AND ~nofit OR (Ratio lt 0) DO BEGIN;;;;20211109 KogureRisei.for L-M method.
    IF makeprint THEN PRINT, "Iter " + STRING(iter, "(I02)") + ":"
    TIC
    ;xk(iter, *) = xa + S # TRANSPOSE(Kmat(*, *, iter - 1L)) # DIAG_MATRIX(invSe) # (Y - Kmat(*, *, iter - 1L) # xa);2.30式　xの更新, Sは更新した共分散行列　この式を、豊岡さん修論のA4.5式にしなくてはいけないのでは？
    ;xk(iter, *) = xa + S # TRANSPOSE(Kmat(*, *, iter - 1L)) # DIAG_MATRIX(invSe) # (Y - Fxk(iter - 1l,*) + Kmat(*, *, iter - 1L) # (reform(xk(iter-1,*)) - xa));rodgerさん本の5.09式に従った式に書き換えてみた。ただし、y-F(x)の項は小さいとして無視
    ;xk(iter,*) = xk(iter-1,*) + LA_INVERT(((1d + v)*invSa + TRANSPOSE(Kmat(*, *, iter - 1L)) # DIAG_MATRIX(invSe) # Kmat(*, *, iter - 1L)) , /DOUBLE, STATUS = status) # (TRANSPOSE(Kmat(*, *, iter - 1L)) # DIAG_MATRIX(invSe) # (Y - Fxk(iter - 1l,*)) - invSa # (reform(xk(iter-1,*)) - xa));20211101 KogureRisei. Levenberg - Marquardt Method[Inverse methods for atmospheric sounding. page92]
    II = invSa;20211101 KogureRisei.for L-M method.
    II(where(invSa ne 0)) = 1;20211101 KogureRisei.for L-M method.
    xk(iter,*) = xk(iter-1,*) + LA_INVERT(((1d + v)*invSa + TRANSPOSE(Kmat(*, *, iter - 1L)) # DIAG_MATRIX(invSe) # Kmat(*, *, iter - 1L)) , /DOUBLE, STATUS = status) # (TRANSPOSE(Kmat(*, *, iter - 1L)) # DIAG_MATRIX(invSe) # (Y - Fxk(iter - 1l,*)) - invSa # (reform(xk(iter-1,*)) - xa));20211101 KogureRisei. Levenberg - Marquardt Method[Inverse methods for atmospheric sounding. page92]
    ;xk(iter, *) = xk(iter - 1, *) + S # (MATRIX_MULTIPLY(Kmat(*, *, iter - 1L), DIAG_MATRIX(invSe), /ATRANSPOSE) # (Y - Fxk(iter - 1, *)) - MATRIX_MULTIPLY(invSa, xk(iter - 1, *) - xa, /BTRANSPOSE))
    IF fitalbedo THEN S_albedo(*, iter) = EXP(xk(iter, 0 : nwn - 1L)) ELSE S_albedo(*, iter) = S_albedo(*, iter - 1L)
    FOR i = 0, nspecies - 1 DO X(i, *, iter) = xk(iter, nwn * fitalbedo + i * nlayer : nwn * fitalbedo + (i + 1) * nlayer - 1)
    Radiusm(*, *, iter) = Radiusm(*, *, iter - 1)
    Radiuss(*, *, iter) = Radiuss(*, *, iter - 1)
    IF fitradius THEN FOR i = 0, naerosols - 1 DO Radiusm(i, *, iter) = EXP(xk(iter, nwn * fitalbedo + (i + nspecies) * nlayer : nwn * fitalbedo + (i + nspecies + 1) * nlayer - 1))
    IF fitsigma THEN FOR i = 0, naerosols - 1 DO Radiuss(i, *, iter) = EXP(xk(iter, nwn * fitalbedo + (i + nspecies + naerosols * fitradius) * nlayer : nwn * fitalbedo + (i + nspecies + naerosols * fitradius + 1) * nlayer - 1))
    ztpn(*, *, iter) = UPDATE_ZTPN(ztpn(*, *, iter - 1), X(*, *, iter), onlyfitspecies)
    IF makeplot THEN p = PLOT_PROFILES(ztpn, gas_species, aerosols_species, onlyfitspecies, iter, itermax, path_output, Radiusm, Radiuss, fitradius, fitsigma, wn_measured, S_albedo, _EXTRA = extra)
    radtemp = CREATE_SYNTHETIC_SPECTRA_Limb_shoheiDatabase(orbit, alt_range, wn_apriori, z, ztpn(*, *, iter), Ht, geometry, dwn, FOV, Radiusm(*, *, iter), Radiuss(*, *, iter), wn_albedo, S_albedo(*, iter), jactemp, $
      ILS = ILS, NAMEPLUSFILE = 'step' + STRING(iter, "(I02)") + '_', _EXTRA = extra)
    FOR i = 0, nspectra - 1 DO BEGIN
;      rad(*, i, iter) = MEANINTERP(wn_apriori, radtemp(*, i), wn_measured)
;      FOR j = 0, nlayer - 1L DO FOR k = 0, 1L + nspecies + naerosols * 2L - 1L DO jac(*, j, i, k, iter) = MEANINTERP(wn_apriori, jactemp(*, j, i, k), wn_measured)
      rad(*, i, iter) = INTERPOL(radtemp(*, i), wn_apriori, wn_measured);;;;;;/////kogure risei 2021/12/14
      FOR j = 0, nlayer - 1L DO FOR k = 0,  1L + nspecies + naerosols * 2L - 1 DO jac(*, j, i, k, iter) = INTERPOL(jactemp(*, j, i, k), wn_apriori, wn_measured);;;;;/////kogure risei 2021/12/14
    ENDFOR

    Kmat(*, *, iter) = GET_K(jac(*, *, *, *, iter), ztpn(*, *, iter), Radiusm(*, *, iter), Radiuss(*, *, iter), S_albedo(*, iter), onlyfitspecies, regions, wn_measured, fitradius, fitsigma, nlayer, nspectra, fitalbedo);ヤコビアンの更新
    Fxk(iter, *) = GET_Y(rad(*, *, iter), nspectra, regions);出力輝度の記録

    IF makeprint THEN PRINT, "DOF for a_priori:"
    A = GET_A(z, Kmat(*, *, iter), invSa, invSe, nlayer, nwn, gas_species, aerosols_species, onlyfitspecies, iter, path_output, S, invSaKinvSeK, fitradius, fitsigma, _EXTRA = extra)
;    convX(iter) = ABS((xk(iter - 1, *) - xk(iter, *)) # MATRIX_MULTIPLY(S, xk(iter - 1, *) - xk(iter, *), /BTRANSPOSE));xの事後確率密度関数 2.25式  または　conversion test（豊岡さん修論のA5.1式）
    convX(iter) = ABS((xk(iter - 1, *) - xk(iter, *)) # MATRIX_MULTIPLY(invSaKinvSeK, xk(iter - 1, *) - xk(iter, *), /BTRANSPOSE));上の式からSをS-1に変更。20210929 小暮。xの事後確率密度関数 2.25式  または　conversion test（豊岡さん修論のA5.1式）
    
  ;  Sy = dblarr(1,nspectra*nregions)
    IF N_ELEMENTS(invSa) GT 1 THEN Sy = DIAG_MATRIX(1d/invse) # (1d/(Kmat(*, *, iter) # LA_INVERT(invSa, /DOUBLE, STATUS = status) # transpose(Kmat(*, *, iter)) + DIAG_MATRIX(1d/invse))) # DIAG_MATRIX(1d/invse) else Sy = DIAG_MATRIX(1d/invse) # (1d/(Kmat(*, *, iter) # (1d/invSa) # transpose(Kmat(*, *, iter)) + DIAG_MATRIX(1d/invse))) # DIAG_MATRIX(1d/invse) ;20211119 kogure risei. rodger 5.27式
    IF N_ELEMENTS(Sy) GT 1 THEN kai2(iter) = ABS(transpose(y - Fxk(iter, *)) # La_INVERT(Sy) # (y - Fxk(iter, *))) else kai2(iter) = ABS(transpose(y - Fxk(iter, *)) # (1d/Sy) # (y - Fxk(iter, *)));;20211119 kogure risei 評価関数のため加えた. rodger 5.32式    
    
    t = TOC()

    IF makeplot THEN RMS(*, iter) = PLOT_SPECTRA_JACOSPAR(wn_measured, rad_measured, rad, drad_measured, regions, nspectra, geometry(0, *), sequence, path_output, iter)
    IF makeplot THEN BEGIN
      DEVICE, GET_SCREEN_SIZE = screen_size
      layout = GET_LAYOUT(nspectra)
      FOR i = 0, nspectra - 1 DO p = PLOT(INDGEN(iter + 1), RMS(i, 0 : iter), $
        /CURRENT, XTITLE = 'Iterations', YTITLE = 'RMS', LAYOUT = [layout, i + 1], $
        COLOR = 'k', /YLOG, THICK = 2, TITLE = "Spectra " + STRING(i, "(I02)"), DIMENSIONS = ROUND(screen_size * .9), /BUFFER)

      p.SAVE, path_output + '/Figures/RMS.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
      p.CLOSE
    ENDIF
        
    Ratio = 1d - kai2(iter)/kai2(iter-1);20211101 KogureRisei.for L-M method.
    if ratio ge 0 then v = v/2d;20211101 KogureRisei.for L-M method.
    if ratio lt 0 then v = v*2d ;20211101 KogureRisei.for L-M method.
    if v gt 1d7 then v = 1d7
    ;if v lt 800d then v = 800
    print,'kai=', kai2(iter) , 'ratio=',ratio,', v=', v;20211101 KogureRisei.for L-M method.
    
;    IF makeprint THEN BEGIN
;      PRINT, "ConvX: " + STRING(convX(iter), "(E10.3)") + "(" + STRING(epsilonX, "(E10.3)") + ")"
;      PRINT, STRING(t, "(F9.1)") + "s"
;      PRINT, ""
;    ENDIF
    SAVE, xk, FILENAME = path_output + 'Step_' + STRING(iter, "(I02)") + '.sav'
    iter = iter + 1L
  ENDWHILE
  iter = iter - 1L
  converged = iter LT itermax AND ~nofit
  IF makeplot AND ~nofit THEN RMS(*, iter) = PLOT_SPECTRA_JACOSPAR(wn_measured, rad_measured, rad, drad_measured, regions, nspectra, geometry(0, *), sequence, path_output, iter)

  ; -------------------------
  ; CALCULATION OF THE ERRORS
  ; -------------------------
  ; Retrieval error
  G = S # MATRIX_MULTIPLY(Kmat(*, *, iter), DIAG_MATRIX(invSe), /ATRANSPOSE);3.27式
  RetErrorMat = ABS((xk(iter - 1, *) - xa) # (A - IDENTITY(N_ELEMENTS(A(*, 0))))) + ABS(G # invSe^(-1))
  IF fitalbedo EQ 1 THEN RetErrorAlbedo = RetErrorMat(0 : nwn - 1L) ELSE RetErrorAlbedo = DBLARR(N_ELEMENTS(regions))
  RetError = DBLARR(nspecies + (fitradius + fitsigma) * naerosols, nlayer)
  FOR i = 0, nspecies + (fitradius + fitsigma) * naerosols - 1 DO RetError(i, *) = RetErrorMat(nwn * fitalbedo + i * nlayer : nwn * fitalbedo + (i + 1) * nlayer - 1)

  ; A-priori and retrieval noise covariance
  Ss = S # invSa # S
  Sm = MATRIX_MULTIPLY(S, Kmat(*, *, iter), /BTRANSPOSE) # DIAG_MATRIX(invSe) # Kmat(*, *, iter) # S

  invSdiag = DIAG_MATRIX(S)
  IF fitalbedo EQ 1 THEN dXalbedo = invSdiag(0 : nwn - 1L)
  dX = DBLARR(nspecies + (fitradius + fitsigma) * naerosols, nlayer)
  FOR i = 0, nspecies + (fitradius + fitsigma) * naerosols - 1 DO dX(i, *) = invSdiag(nwn * fitalbedo + i * nlayer : nwn * fitalbedo + (i + 1) * nlayer - 1)
  IF makeplot THEN p = PLOT_PROFILES(ztpn, gas_species, aerosols_species, onlyfitspecies, iter, itermax, path_output, Radiusm, Radiuss, fitradius, fitsigma, wn_measured, S_albedo, $
    ERRORDEN = RetError, ERRORALBEDO = RetErrorAlbedo, ADDNAME = '_conv', _EXTRA = extra)

  ; ------------
  ; PLOT THE RMS
  ; ------------
  IF makeplot THEN BEGIN
    DEVICE, GET_SCREEN_SIZE = screen_size
    layout = GET_LAYOUT(nspectra)
    FOR i = 0, nspectra - 1 DO p = PLOT(INDGEN(iter + 1), RMS(i, 0 : iter), $
      /CURRENT, XTITLE = 'Iterations', YTITLE = 'RMS', LAYOUT = [layout, i + 1], $
      COLOR = 'k', /YLOG, THICK = 2, TITLE = "Spectra " + STRING(i, "(I02)"), DIMENSIONS = ROUND(screen_size * .9), /BUFFER)

    p.SAVE, path_output + '/Figures/RMS.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
    p.CLOSE
  ENDIF
  IF makeprint THEN BEGIN
    IF iter EQ itermax - 1L THEN PRINT, "Inversion has NOT converged." ELSE PRINT, "Inversion has converged."
  ENDIF

  ; --------------------
  ; Save results to file
  ; --------------------
  SAVE, Kmat, S, invSa, invSe, A, X, dX, Y, Fxk, xk, rad, jac, xa, convX, wn_measured, rad_measured, drad_measured, RetError, RetErrorAlbedo, $
    nlayer, nwn, nspectra, alt_range, wnrange, ztpn, onlyfitspecies, S_albedo, RMS, Sm, Ss, iter, fitradius, fitsigma, converged, $ ;Radiusfinalm, Radiusfinals,
    FILENAME = path_output + 'Results.sav'
  SAVE, ztpn, RetError, Radiusm, Radiuss, S_albedo, RetErrorAlbedo, iter, onlyfitspecies, fitradius, fitsigma, FILENAME = path_output + 'ResultsZTPN.sav'
  RETURN, converged
END
