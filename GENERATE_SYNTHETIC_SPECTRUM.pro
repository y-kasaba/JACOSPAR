FUNCTION GENERATE_SYNTHETIC_SPECTRUM, orbit, alt_range, $
                                      STARTRADIUS = StartRadius, $
                                      STARTVARIANCE = StartVariance, $
                                      MAKEPRINT = makeprint, $
                                      DWN = dwn, $
                                      NORADIUSPREFIT = noradiusprefit, $
                                      CLEAN = clean, NAMEPLUS = nameplus, $
                                      FIGVEC = figvec, NONOISE = nonoise, FACTORINI = factorini, _EXTRA = extra
                                      
  dz = 5D
  
  IF ~KEYWORD_SET(makeprint) THEN makeprint = 0
  IF ~KEYWORD_SET(clean) THEN clean = 0
  IF ~KEYWORD_SET(nonoise) THEN nonoise = 0
  IF ~KEYWORD_SET(StartRadius) THEN StartRadius = [1.5D-4, 3D-4]
  IF ~KEYWORD_SET(StartVariance) THEN StartVariance = [.3D, .1D]
  IF ~KEYWORD_SET(noradiusprefit) THEN noradiusprefit = 0
  IF ~KEYWORD_SET(factorini) THEN BEGIN
    factorini = DBLARR(5) 
    factorini(*) = 1D
  ENDIF
  IF KEYWORD_SET(extra) THEN BEGIN
    names = TAG_NAMES(extra)
    FOR i = 0, N_ELEMENTS(names) - 1 DO BEGIN
      CASE STRLOWCASE(names(i)) OF
        'dz': dz = extra.dz
        ELSE:
      ENDCASE
    ENDFOR
  ENDIF

  ; ------------------------
  ; Setting up the save path
  ; ------------------------
  pathnames = PATH_MANAGEMENT()
  path_output = pathnames.results + orbit + '_Synth/'
  orbitnumber = orbit
  orbitstring = STRING(orbit, "(I05)")
  IF KEYWORD_SET(nameplus) THEN path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplus + "/" ELSE nameplus = ""
  IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
  IF clean THEN BEGIN
    RESULT = DIALOG_MESSAGE("Are you sure you want to clean the directory?", /CENTER, /QUESTION)
    IF RESULT EQ "No" THEN RETURN, 0
    IF FILE_TEST(path_output + "/Figures") THEN BEGIN
      liste = FILE_SEARCH(path_output + '/Figures/*.*')
      IF N_ELEMENTS(liste) GT 1 THEN FOR i = 0, N_ELEMENTS(liste) - 1 DO FILE_DELETE, liste(i), /QUIET
    ENDIF
    liste = FILE_SEARCH(path_output + '*.*')
    IF N_ELEMENTS(liste) GT 1 THEN FOR i = 0, N_ELEMENTS(liste) - 1 DO FILE_DELETE, liste(i), /QUIET
  ENDIF
  IF ~FILE_TEST(path_output + '/Figures/', /DIRECTORY) THEN FILE_MKDIR, path_output + '/Figures/'

  ; ----------------
  ; Wavenumber range
  ; ----------------
  wnrange = [[3780D, 4020D],[4020D, 4500D], [4500D, 5180D], [5180D, 5800D], [5800D, 6600D], [6600D, 7700D],[7700D, 8500D], [8500D, 9000D], [9000D, 10000D]] ; CO - CO2 - H2O - NULL - H2O - CO2 - H2O
  gas_species = ['co2', 'co', 'h2o']
  aerosols_species = ['dust', 'ice']
  naerosols_max = N_ELEMENTS(aerosols_species)
  onlyfitspecies = INDGEN(5)
  nspecies = N_ELEMENTS(onlyfitspecies)

  ; --------------------------
  ; Setting up some parameters
  ; --------------------------
  ngas = N_ELEMENTS(gas_species)
  naerosols = N_ELEMENTS(aerosols_species)
  FOV = 1.2D-3 * 2D / !DPI * 180D
  IF ~KEYWORD_SET(dwn) THEN dwn = 1D-2
  resol = 31.2461D
  kB = 1.38D-16 ; [cgs]
  gM = 371.1D ; [cgs]
  Nav = 6.02D23
  MM_CO2 = 44.01D-3
      
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
  
  ; --------------------------------------------------------
  ; Reshape altitude range depending on the scanned geometry
  ; --------------------------------------------------------
  IF KEYWORD_SET(dz) THEN z = DINDGEN(CEIL(alt_range(1) / dz) + 1) * dz ELSE z = DOUBLE(geometry(0, *))
  nlayer = N_ELEMENTS(z)  
  
  ; -----------------------------------
  ; Defining the radii of the particles
  ; -----------------------------------
  radiusS_mod = EXP(SQRT(ALOG(StartVariance + 1D)))
  radiusM_mod = StartRadius * EXP(-2.5D * (ALOG(radiusS_mod))^2)
  Radiusm = DBLARR(naerosols_max, nlayer) ; Radius mean modal value
  Radiuss = DBLARR(naerosols_max, nlayer) ; Radius variance modal value
  IF N_ELEMENTS(StartRadius(0, *)) EQ 1 THEN FOR i = 0, nlayer - 1 DO Radiusm(*, i) = radiusM_mod ELSE Radiusm = radiusM_mod ; [cm]
  IF N_ELEMENTS(StartRadius(0, *)) EQ 1 THEN FOR i = 0, nlayer - 1 DO Radiuss(*, i) = radiusS_mod ELSE Radiuss = radiusS_mod
  Radiusfinalm = DBLARR(naerosols_max, nlayer)
  Radiusfinals = DBLARR(naerosols_max, nlayer)
  IF N_ELEMENTS(StartRadius(0, *)) EQ 1 THEN FOR i = 0, nlayer - 1 DO Radiusfinalm(*, i) = StartRadius ELSE Radiusfinalm = radiusM_mod / EXP(-2.5D * (ALOG(radiusS_mod))^2) ; [cm]
  IF N_ELEMENTS(StartRadius(0, *)) EQ 1 THEN FOR i = 0, nlayer - 1 DO Radiusfinals(*, i) = StartVariance ELSE Radiusfinals = EXP(ALOG(radiusS_mod)^2) - 1D

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
  rad_measured = READ_SPECTRA(orbit, WN = wn_measured, ERROR = drad_measured)
  drad_measured(*) = .01D
  IF N_ELEMENTS(rad_measured) EQ 0 THEN RETURN, !NULL
  inds = INTARR(N_ELEMENTS(wn_measured))
  FOR i = 0, nwindows - 1 DO inds = inds OR (wn_measured GE wnrange(0, i) AND wn_measured LE wnrange(1, i))
  inds = WHERE(inds, /NULL)
  rad_measured = rad_measured(*, sequence)
  drad_measured = drad_measured(*, sequence)
  wnall = wn_measured
  radall = rad_measured
  wn_measured = wn_measured(inds)
  rad_measured = rad_measured(inds, *)
  drad_measured = drad_measured(inds, *)
  snr = rad_measured / drad_measured
  nspec = N_ELEMENTS(rad_measured(0, *))
  nwn = N_ELEMENTS(wn_measured)

  ; ---------------------------
  ; Get the a-priori atmosphere
  ; ---------------------------
  X = DBLARR(nspecies, nlayer)
  
  ztpn = GETZTPN(MEAN(geometry(4, *)), MEAN(geometry(5, *)), MEAN(geometry(6, *)), MEAN(geometry(7, *)), z, cov, StartRadius(1))
  IF KEYWORD_SET(startfrom) THEN BEGIN
    IF ~FILE_TEST(pathnames.results + startfrom + '/ResultsZTPN.sav') THEN BEGIN
      PRINT, 'Run ' + startfrom + ' not found. Aborting.'
      RETURN, -1
    ENDIF
    onlyfitspeciessave = onlyfitspecies
    Radiusmsave = Radiusm
    Radiusssave = Radiuss
    RESTORE, pathnames.results + startfrom + '/ResultsZTPN.sav'
    ztpn(5 + onlyfitspecies, *) = ztpn(5 + onlyfitspecies, *, iter)
    Radiusmsave(*, *, 0) = Radiusm(*, *, iter)
    Radiusm = Radiusmsave
    Radiusssave(*, *, 0) = Radiuss(*, *, iter)
    Radiuss = Radiusssave
    onlyfitspecies = onlyfitspeciessave
  ENDIF
  Ht = kb * Nav * ztpn(1, *) / gM * MM_CO2 * 1D-5 ; [km]
  corrlength = MEAN(Ht) / 2D
  FOR i = 0, 4 DO ztpn(5 + i, *) *= factorini(i)
  
  IF ~noradiusprefit THEN BEGIN
    resultsprefit = RADIUSPREFIT(orbit, alt_range, ztpn, geometry, FOV, wnall, radall, _EXTRA = extra)
    Radiusm(0, *) = resultsprefit(*, 0)
    Radiuss(0, *) = resultsprefit(*, 1)
    Radiusfinalm(0, *) = resultsprefit(*, 2)
    Radiusfinals(0, *) = resultsprefit(*, 3)
  ENDIF
  
  ; ------------------------------------
  ; Calculate the apriori surface albedo
  ; ------------------------------------
  albedoFile = READ_ASCII(pathnames.albedo + 'true_albedo_orb_0937_4_p_63_scan_482_ARS.dat')
  albedoFile = albedoFile.FIELD1
  S_albedo = DBLARR(nwn)
  wn_albedo = wn_measured
  C = 0.4D
  IF KEYWORD_SET(setalbedo) THEN BEGIN
    CASE N_ELEMENTS(setalbedo) OF
      1: S_albedo(*) = setalbedo
      2: S_albedo = INTERPOL(setalbedo, [4000D, 8000D], wn_albedo)
      ELSE:
    ENDCASE
  ENDIF ELSE S_albedo = INTERPOL(albedoFile(1, *), albedoFile(0, *), wn_albedo)

  p = PLOT_PROFILES(ztpn, gas_species, aerosols_species, onlyfitspecies, 0, 1, path_output, Radiusm, Radiuss, 1, wn_measured, S_albedo, _EXTRA = extra)
  
  ; -----------------------
  ; Calculate the radiances
  ; -----------------------
  rad = DBLARR(nwn, nspectra)
  R = CREATE_SYNTHETIC_SPECTRA(orbit, alt_range, wn_apriori, z, ztpn, geometry, dwn, FOV, Radiusm, Radiuss, wn_albedo, S_albedo(*, 0), ILS = ILS, DWNAERO = 1D, /NOJACOBIAN, NAMEPLUSDIR = '_Synth', _EXTRA = extra)
  FOR i = 0, nspectra - 1 DO rad(*, i, 0) = MEANINTERP(wn_apriori, R(*, i), wn_measured)
  drad_measured = rad(*, *, 0) / SNR
  
  ; -------------
  ; Add the noise
  ; -------------
  IF NOT nonoise THEN FOR i = 0, nspectra - 1 DO rad(*, i) = rad(*, i) + RANDOMU(1, N_ELEMENTS(drad_measured), /DOUBLE, /NORMAL) * drad_measured
  
  ; ----------------
  ; Plot the spectra
  ; ----------------
  DEVICE, GET_SCREEN_SIZE = screen_size
  layout = GET_LAYOUT(nspectra)
  IF ~KEYWORD_SET(figvec) THEN BEGIN
    figvec = OBJARR(nspectra) 
    addfig = 0 
  ENDIF ELSE addfig = 1
  FOR i = 0, nspectra - 1 DO $
    IF addfig THEN figvec(i) = PLOT(wn_measured, rad(*, i), LAYOUT = figvec(i)) $
    ELSE figvec(i) = PLOT(wn_measured, rad(*, i), XTITLE = 'Wavenumber [$cm^{-1}$]', YTITLE = 'Radiance [CGS]', DIMENSIONS = ROUND(screen_size * .9), $
        TITLE = STRING(z(i), "(F5.1)") + " km (#" + STRING(sequence(i), "(I3)") + ")", $
        CURRENT = (i GT 0) ? 1 : 0, LAYOUT = [layout, i + 1], XRANGE = [wnrange(0, 0), wnrange(1, nwindows - 1)])
  
  figvec(0).SAVE, path_output + '/Figures/Spectra.jpg' , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
  figvec(0).CLOSE
  
  ; ----------------
  ; Save the spectra
  ; ----------------
  OPENW, lun, path_output + 'Spectra.dat', /GET_LUN
  WRITEU, lun, DOUBLE(nwn)
  WRITEU, lun, DOUBLE(nspectra)
  WRITEU, lun, DOUBLE(wn_measured)
  WRITEU, lun, DOUBLE(rad)
  CLOSE, lun  
  
  OPENW, lun, path_output + "Spectra.txt", /get_lun, WIDTH = 1000L
  PRINTF, lun, STRING(wn_measured, "(F8.2)")
  FOR i=0, nspectra - 1 DO PRINTF, lun, STRING(rad(*, i), "(F8.6)")
  FREE_LUN, lun
  
  RETURN, figvec
END