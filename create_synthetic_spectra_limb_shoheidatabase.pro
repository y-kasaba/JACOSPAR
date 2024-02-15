FUNCTION CREATE_SYNTHETIC_SPECTRA_Limb_shoheiDatabase, orbit, $
  alt_range, $
  wn, $
  z, $
  ztpn, $
  Ht, $
  geometry, $
  dwn, $
  FOV, $
  Radiusm, $
  Radiuss, $
  wn_albedo, $
  S_albedo, $
  Jacconv, $
  dRadconv, $

  NAMEPLUSDIR = nameplusdir, $
  NAMEPLUSFILE = nameplusfile, $
  TOFIT = tofit, $
  RADNOCONV = Radnoconv, $
  JACNOCONV = Jacnoconv, $
  DRADNOCONV = dRadnoconv, $
  DJACNOCONV = dJacnoconv, $
  WNRANGE = wnrange, $

  CONVOLVEALLINONE = convolveallinone, $
  NUMERICALJACOBIAN = numericaljacobian, $

  TAU_ABS = tau_abs, $
  TAU_AERO_SCA = tau_aero_sca, $
  BETA_AERO = Beta_aero, $
  DRDABS = dRdabs, $
  DRDSCA = dRdsca, $
  DRDALB = dRdalb, $
  DELTADRDABS = deltadRdabs, $
  DELTADRDSCA = deltadRdsca, $
  GETFILENAMES = getfilenames, $

  _EXTRA = extra

  !EXCEPT = 0

  ; ------------------------------------------
  ; Species to be fitted and wavenumber ranges
  ; ------------------------------------------
  IF ~KEYWORD_SET(wnrange) THEN wnrange = [[3780D, 4020D],[4020D, 4500D], [4500D, 5180D], [5180D, 5800D], [5800D, 6600D], [6600D, 7700D],[7700D, 8500D], [8500D, 9000D], [9000D, 10000D]] ; CO - CO2 - H2O - NULL - H2O - CO2 - H2O
  gas_species = ['co2', 'co', 'h2o']
  aerosols_species = ['dust', 'waterice'];;kogure 20210812 'ice' to 'waterice'
  onlyfitspecies = [INDGEN(N_ELEMENTS(gas_species)), N_ELEMENTS(gas_species) + INDGEN(N_ELEMENTS(aerosols_species))]

  ; -----------------------
  ; Setting up the keywords
  ; -----------------------
  nojacobian = 0
  nodelete = 0
  dwnaero = dwn
  fitradius = 0
  fitsigma = 0
  IF ~KEYWORD_SET(makeprint) THEN makeprint = 0L
  IF ~KEYWORD_SET(nameplusfile) THEN nameplusfile = ""
  IF ~KEYWORD_SET(getfilenames) THEN getfilenames = 0L
  IF ~KEYWORD_SET(convolveallinone) THEN convolveallinone = 0L
  IF ~KEYWORD_SET(numericaljacobian) THEN numericaljacobian = 0L
  IF KEYWORD_SET(extra) THEN BEGIN
    names = TAG_NAMES(extra)
    FOR i = 0, N_ELEMENTS(names) - 1L DO BEGIN
      CASE STRLOWCASE(names(i)) OF
        'nojacobian': nojacobian = extra.nojacobian
        'ils': ILS = extra.ILS
        'onlyfitspecies': onlyfitspecies = extra.onlyfitspecies
        'dwnaero': dwnaero = extra.dwnaero
        'fitradius': fitradius = extra.fitradius
        'fitsigma': fitsigma = extra.fitsigma
        ELSE:
      ENDCASE
    ENDFOR
  ENDIF
  nspecies = N_ELEMENTS(onlyfitspecies)
  IF N_ELEMENTS(onlyfitspecies(WHERE(onlyfitspecies LE 2, /NULL))) GT 0 THEN $
    gas_species_fitted = gas_species(onlyfitspecies(WHERE(onlyfitspecies LE 2, /NULL))) $
  ELSE gas_species_fitted = !NULL
  IF N_ELEMENTS(onlyfitspecies(WHERE(onlyfitspecies - 3 GE 0, /NULL))) GT 0 THEN $
    aerosols_species_fitted = aerosols_species(onlyfitspecies(WHERE(onlyfitspecies - 3 GE 0, /NULL)) - 3) $
  ELSE aerosols_species_fitted = !NULL
  ngas = N_ELEMENTS(gas_species)
  naerosols = N_ELEMENTS(aerosols_species)
  ngas_fitted = N_ELEMENTS(gas_species_fitted)
  naerosols_fitted = N_ELEMENTS(aerosols_species_fitted)
  nspectra = N_ELEMENTS(geometry(0, *))
  nlayer = N_ELEMENTS(z)
  epsilon = 1D-1

  ; ------------------------
  ; Setting up the save path
  ; ------------------------
  pathnames = PATH_MANAGEMENT()
  path_output = pathnames.results + orbit
  IF KEYWORD_SET(nameplusdir) THEN path_output = path_output + nameplusdir
  IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output

  ; --------------------------
  ; Setting up some parameters
  ; --------------------------
  kB = 1.38D-16 ; [cgs]
  gM = 371.1D ; [cgs]
  FBord = 0D
  resol = 31.2461D

  ; --------------------------------------------------------
  ; Reshape altitude range depending on the scanned geometry
  ; --------------------------------------------------------
  hmin = 0D
  hmax = MAX(geometry(0, *))
  alt_range = [hmin, hmax]

  ; ---------------------------
  ; Spectral windows definition
  ; ---------------------------
  nwindows = N_ELEMENTS(wnrange(0, *))
  nwn = ROUND((wnrange(1L, *) - wnrange(0, *)) / dwn) + 1L
  wn = DBLARR(TOTAL(nwn))
  ind = 0L
  FOR i = 0, nwindows - 1L DO BEGIN
    wn(ind : ind + nwn(i) - 1L) = wnrange(0, i) + DINDGEN(nwn(i)) * dwn
    ind += nwn(i)
  ENDFOR

  ; ------------
  ; Load the ILS
  ; ------------
  IF ~KEYWORD_SET(ILS) THEN BEGIN
    wn_ILS = -60D + DINDGEN(ROUND(120D / dwn) + 1L) * dwn
    ILS = GAUSSIAN(wn_ILS, 1D, resol, 0D)
    ILS = ILS / INTEGRATE(wn_ILS,ILS)
  ENDIF

  ; -----------------------------------------
  ; Calculate the gas absorption coefficients
  ; -----------------------------------------
  nwn_k = LONARR(nwindows)
  FOR i = 0, nwindows - 1L DO nwn_k(i) = ROUND((wnrange(1L, i) - wnrange(0, i) + 2 * FBord) / dwn) + 1L
  IF ngas GT 0 THEN BEGIN
    ACS = {k : DBLARR(MAX(nwn_k), nlayer, ngas), nwn : 0L, wn : DBLARR(MAX(nwn_k))}
    ACS = REPLICATE(ACS, nwindows)
    FOR i = 0, nwindows - 1L DO BEGIN
      FOR j = 0, ngas - 1L DO BEGIN
        ACS[i].nwn = nwn_k(i)
        temp = INTERPOL_ACS(gas_species(j), [wnrange(0, i) - FBord, wnrange(1L, i) + FBord], wnAbs, ztpn(2, *)*1d-2, ztpn(1L, *), dwn);Unit??
        IF N_ELEMENTS(temp) EQ 0 THEN CONTINUE
        ACS[i].k(0 : ACS[i].nwn - 1L, *, j) = temp
        ACS[i].wn[0 : ACS[i].nwn - 1L] = wnAbs(0) + DINDGEN(nwn_k(i)) * dwn
      ENDFOR
    ENDFOR
  ENDIF ELSE BEGIN
    ACS = {k : DBLARR(MAX(nwn_k), nlayer), nwn : 0L, wn : DBLARR(MAX(nwn_k))}
    ACS = REPLICATE(ACS, nwindows)
    FOR i = 0, nwindows - 1L DO FOR j = 0, ngas - 1L DO ACS[i].nwn = nwn_k(i)
  ENDELSE

  ; ----------------------------------------------------------------------------------------
  ; Prepare the matrix that will contain the aerosols absorption and scattering coefficients
  ; ----------------------------------------------------------------------------------------
  nwn_aero = LONARR(nwindows)
  FOR i = 0, nwindows - 1L DO nwn_aero(i) = ROUND((wnrange(1L, i) - wnrange(0, i) + 2 * FBord) / dwnaero) + 1L
  Beta_aero = {wn : DBLARR(MAX(nwn_k)), $
    nwn : 0L, $
    kabs : DBLARR(MAX(nwn_aero), naerosols, nlayer), $
    ksca : DBLARR(MAX(nwn_aero), naerosols, nlayer), $
    dkabsdrm : DBLARR(MAX(nwn_aero), naerosols, nlayer), $
    dkscadrm : DBLARR(MAX(nwn_aero), naerosols, nlayer), $
    dkabsds : DBLARR(MAX(nwn_aero), naerosols, nlayer), $
    dkscads : DBLARR(MAX(nwn_aero), naerosols, nlayer), $
    Qext : DBLARR(MAX(nwn_aero), naerosols, nlayer)}
  Beta_aero = REPLICATE(Beta_aero, nwindows, 2)
  FOR i = 0, nwindows - 1L DO Beta_aero[i].nwn = nwn_aero(i)

  ; -----------------------
  ; Calculate the a-prioris
  ; -----------------------
  IF numericaljacobian THEN FileNames = STRARR(nwindows * nspectra * (1L + nlayer * (TOTAL(onlyfitspecies LT 3) + TOTAL(onlyfitspecies GE 3) * (1L + fitradius + fitsigma)))) ELSE $
    FileNames = STRARR(nwindows * nlayer * (1L + naerosols * nlayer))
  ind = 0
  FOR i = 0, nwindows - 1L DO BEGIN
    IF numericaljacobian THEN BEGIN
      temp = CREATE_JACOSPAR_INPUT_FILES_NUMERICAL(orbit, wnrange(*, i) + [-FBord, FBord], gas_species, aerosols_species, alt_range, dwn, FOV, z, ztpn, Ht, geometry, Radiusm, Radiuss, wn_albedo, S_albedo, $
        'Orbit_' + nameplusfile + orbit + "_" + STRING(i, "(I02)"), PATH_OUTPUT = path_output + '/', epsilon, fitradius, fitsigma, $
        GAS_ACS = ACS[i].k(0 : ACS[i].nwn - 1L, *, *), TAU_AERO_SCA = tau_aero_sca, TAU_ABS = tau_abs, _EXTRA = extra)
    ENDIF ELSE BEGIN
      temp = CREATE_JACOSPAR_INPUT_FILES_Limb_shoheidatabase(orbit, wnrange(*, i) + [-FBord, FBord], gas_species, aerosols_species, alt_range, dwn, FOV, z, ztpn, Ht, geometry, Radiusm, Radiuss, wn_albedo, S_albedo, $
        PATH_OUTPUT = path_output + '/', SAVENAME = 'Orbit_' + nameplusfile + orbit + "_" + STRING(i, "(I02)"), $
        WN = temp0, AERO_K_ABS = temp1, AERO_K_SCA = temp2, DAERODRM_K_ABS = temp3, DAERODRM_K_SCA = temp4,  DAERODS_K_ABS = temp5, DAERODS_K_SCA = temp6, $
        GAS_ACS = ACS[i].k(0 : ACS[i].nwn - 1L, *, *), TAU_AERO_SCA = tau_aero_sca, TAU_ABS = tau_abs, QEXT = temp7, _EXTRA = extra)
      Beta_aero[i].wn(0 : Beta_aero[i, 0].nwn - 1L) = temp0
      Beta_aero[i].kabs(0 : Beta_aero[i].nwn - 1L, *, *) = temp1
      Beta_aero[i].ksca(0 : Beta_aero[i].nwn - 1L, *, *) = temp2
      Beta_aero[i].dkabsdrm(0 : Beta_aero[i].nwn - 1L, *, *) = temp3
      Beta_aero[i].dkscadrm(0 : Beta_aero[i].nwn - 1L, *, *) = temp4
      Beta_aero[i].dkabsds(0 : Beta_aero[i].nwn - 1L, *, *) = temp5
      Beta_aero[i].dkscads(0 : Beta_aero[i].nwn - 1L, *, *) = temp6
      Beta_aero[i].Qext(0 : Beta_aero[i].nwn - 1L, *, *) = temp7
    ENDELSE
    FileNames(ind : ind + N_ELEMENTS(temp) - 1L) = temp
    ind = ind + N_ELEMENTS(temp)
  ENDFOR
  FileNames = FileNames(0 : ind - 1L)
  IF getfilenames EQ 1L THEN RETURN, FileNames
  torun = LONARR(N_ELEMENTS(FileNames))
  FOR i = 0, N_ELEMENTS(FileNames) - 1L DO torun(i) = ~FILE_TEST(FileNames(i) + '.res')
  IF TOTAL(torun) GT 0L THEN status = PARALLEL_JACOSPAR(FileNames(WHERE(torun, /NULL)), _EXTRA = extra)

  Radnoconv = DBLARR(TOTAL(nwn) + ROUND(2D * FBord / dwn), nspectra)
  dRadnoconv = DBLARR(TOTAL(nwn) + ROUND(2D * FBord / dwn), nspectra)
  Jacnoconv = DBLARR(TOTAL(nwn) + ROUND(2D * FBord / dwn), nlayer, nspectra, 1L + ngas_fitted + naerosols_fitted * 3L)
  dJacnoconv = DBLARR(TOTAL(nwn) + ROUND(2D * FBord / dwn), nlayer, nspectra, 1L + ngas_fitted + naerosols_fitted * 3L)
  Radconv = DBLARR(TOTAL(nwn), nspectra)
  dRadconv = DBLARR(TOTAL(nwn), nspectra)
  Jacconv = DBLARR(TOTAL(nwn), nlayer, nspectra, 1L + ngas_fitted + naerosols_fitted * 3L)

  IF numericaljacobian THEN BEGIN
    indwn = 0L
    FOR i = 0, nwindows - 1L DO BEGIN
      ; Get the radiances and Jacobians calculated by JACOSPAR
      Rtemp = READ_JACOSPAR_OUTPUT_FILE(orbit, wnrange(*, i) + [-FBord, FBord], gas_species, aerosols_species, alt_range, dwn, dwnaero, z, geometry, wnloc, $
        SAVENAME = 'Orbit_' + nameplusfile + orbit + "_" + STRING(i, "(I02)"), PATH_OUTPUT = path_output + '/', /NOJACOBIAN, _EXTRA = extra)

      IF i EQ 0 THEN inds = WHERE(wnloc GE wnrange(0, i) - FBord AND wnloc LE wnrange(1L, i)) ELSE $
        IF i EQ nwindows - 1L THEN inds = WHERE(wnloc GE wnrange(0, i) AND wnloc LE wnrange(1L, i) + FBord) ELSE $
        inds = WHERE(wnloc GE wnrange(0, i) AND wnloc LE wnrange(1L, i))
      nwntemp = N_ELEMENTS(inds)

      IF ~convolveallinone THEN BEGIN
        FOR j = 0, nspectra - 1L DO BEGIN
          IF nwn_k(i) GE N_ELEMENTS(ILS) THEN BEGIN
            temp = CONVOL(Rtemp(*, j), ILS, /CENTER, /EDGE_WRAP) * (wn(1L) - wn(0))
            Radconv(indwn : indwn + nwntemp - 1L, j) = temp(inds)
          ENDIF ELSE Radconv(indwn : indwn + nwntemp - 1L, j) = Rtemp(inds, j)
        ENDFOR
      ENDIF

      ; Calculate the absorption Jacobians, convolve them and set into the final vectors
      IF ~nojacobian THEN BEGIN
        Jactemp = DBLARR(N_ELEMENTS(inds), nlayer, nspectra, 1L + ngas_fitted + naerosols_fitted * 3L);;;;inds：計算対象波数
        dJactemp = DBLARR(N_ELEMENTS(inds), nlayer, nspectra, 1L + ngas_fitted + naerosols_fitted * 3L)

        ; Jacobians to the gas VMR
        FOR k = 0, ngas_fitted - 1L DO BEGIN
          FOR ii = 0, nspectra - 1L DO BEGIN
            Rvmr = READ_JACOSPAR_OUTPUT_FILE(orbit, wnrange(*, i) + [-FBord, FBord], gas_species, aerosols_species, alt_range, dwn, dwnaero, z, geometry, wnloc, $
              SAVENAME = 'Orbit_' + nameplusfile + orbit + "_" + STRING(i, "(I02)") + '_' + gas_species(k) + '_L' + STRING(j, "(I03)"), PATH_OUTPUT = path_output + '/', /NOJACOBIAN, _EXTRA = extra)
            FOR j = 0, nlayer - 1L DO Jactemp(*, j, ii, k + 1L) = (Rvmr(*, j) - Rtemp(*, ii)) / (epsilon * ztpn(k + 5, j))
          ENDFOR
        ENDFOR

        ; Jacobians to the aerosols VMR and particle radius
        FOR k = 0, naerosols_fitted - 1L DO BEGIN
          FOR j = 0, nlayer - 1L DO BEGIN
            ; VMR
            Rvmr = READ_JACOSPAR_OUTPUT_FILE(orbit, wnrange(*, i) + [-FBord, FBord], gas_species, aerosols_species, alt_range, dwn, dwnaero, z, geometry, wnloc, $
              SAVENAME = 'Orbit_' + nameplusfile + orbit + "_" + STRING(i, "(I02)") + '_' + aerosols_species(k) + '_L' + STRING(j, "(I03)"), PATH_OUTPUT = path_output + '/', /NOJACOBIAN, _EXTRA = extra)
            FOR ii = 0, nspectra - 1L DO Jactemp(*, j, ii, ngas_fitted + k + 1L) = (Rvmr(*, ii) - Rtemp(*, ii)) / (epsilon * ztpn(k + 5 + ngas, j))

            IF fitradius THEN BEGIN
              ; Aerosols mean radius
              Rradius = READ_JACOSPAR_OUTPUT_FILE(orbit, wnrange(*, i) + [-FBord, FBord], gas_species, aerosols_species, alt_range, dwn, dwnaero, z, geometry, wnloc, $
                SAVENAME = 'Orbit_' + nameplusfile + orbit + "_" + STRING(i, "(I02)") + '_' + aerosols_species(k) + '_M_L' + STRING(j, "(I03)"), PATH_OUTPUT = path_output + '/', /NOJACOBIAN, _EXTRA = extra)
              FOR ii = 0, nspectra - 1L DO Jactemp(*, j, ii, ngas_fitted + naerosols_fitted + k + 1L) = (Rradius(*, ii) - Rtemp(*, ii)) * 1D4
            ENDIF

            IF fitsigma THEN BEGIN
              ; Aerosols variance radius
              Rsigma = READ_JACOSPAR_OUTPUT_FILE(orbit, wnrange(*, i) + [-FBord, FBord], gas_species, aerosols_species, alt_range, dwn, dwnaero, z, geometry, wnloc, $
                SAVENAME = 'Orbit_' + nameplusfile + orbit + "_" + STRING(i, "(I02)") + '_' + aerosols_species(k) + '_S_L' + STRING(j, "(I03)"), PATH_OUTPUT = path_output + '/', /NOJACOBIAN, _EXTRA = extra)
              FOR ii = 0, nspectra - 1L DO Jactemp(*, j, ii, ngas_fitted + 2 * naerosols_fitted + k + 1L) = (Rsigma(*, ii) - Rtemp(*, ii)) / (epsilon * Radiuss(k, j))
            ENDIF
          ENDFOR
        ENDFOR

        Jacnoconv(indwn : indwn + nwntemp - 1L, *, *, *) = Jactemp(inds, *, *, *)

        IF nwn_k(i) GE N_ELEMENTS(ILS) THEN BEGIN
          FOR j = 0, nspectra - 1L DO BEGIN
            FOR k = 0, nlayer - 1L DO BEGIN
              FOR ii = 0, 1L + (ngas_fitted + 3 * naerosols_fitted) - 1L DO BEGIN
                temp = CONVOL(Jactemp(*, k, j, ii), ILS, /CENTER, /EDGE_WRAP) * (wn(1L) - wn(0))
                Jacconv(indwn : indwn + nwntemp - 1L, k, j, ii) = temp(inds)
              ENDFOR
            ENDFOR
          ENDFOR
        ENDIF ELSE BEGIN
          Jacconv(indwn : indwn + nwntemp - 1L, *, *, *) = Jactemp(inds, *, *, *)
        ENDELSE
      ENDIF

      ; Update the wavenumber indices
      indwn = indwn + nwntemp - 1L
    ENDFOR
  ENDIF ELSE BEGIN
    dz = [(z(1L : nlayer - 1L) - z(0 : nlayer - 2)) * 1D5, 1D7 - z(nlayer - 1L) * 1D5]
    dRdalb = DBLARR(TOTAL(nwn) + ROUND(2D * FBord / dwn), nspectra)
    dRdabs = DBLARR(TOTAL(nwn) + ROUND(2D * FBord / dwn), nlayer, nspectra)
    dRdsca = DBLARR(TOTAL(nwn) + ROUND(2D * FBord / dwn), naerosols, nlayer, nspectra)
    indwn = 0L
    FOR i = 0, nwindows - 1L DO BEGIN
      ; Get the radiances and Jacobians calculated by JACOSPAR
      Rtemp = READ_JACOSPAR_OUTPUT_FILE(orbit, wnrange(*, i) + [-FBord, FBord], gas_species, aerosols_species, alt_range, dwn, dwnaero, z, geometry, wnloc, $
        dRdalbtemp, dRdtauabstemp, dRdtauscatemp, deltaRtemp, deltadRdalbtemp, deltadRdtauabstemp, deltadRdtauscatemp, $
        SAVENAME = 'Orbit_' + nameplusfile + orbit + "_" + STRING(i, "(I02)"), PATH_OUTPUT = path_output + '/', _EXTRA = extra)

      ; Get the index range in the wavenumber vector
      IF i EQ 0 THEN inds = WHERE(wnloc GE wnrange(0, i) - FBord AND wnloc LE wnrange(1L, i)) ELSE $;;;;inds：計算対象波数
        IF i EQ nwindows - 1L THEN inds = WHERE(wnloc GE wnrange(0, i) AND wnloc LE wnrange(1L, i) + FBord) ELSE $
        inds = WHERE(wnloc GE wnrange(0, i) AND wnloc LE wnrange(1L, i))
      nwntemp = N_ELEMENTS(inds)
      
;      dRdtauabstemp(*,-1,*)=0d;;;;;;;;;;;;;;;;;;;あとで消す！！！！！！！！！！！！！！！！！！！！！！！！小暮
;      dRdtauscatemp(*,0,-1,*)=0d;;;;;;;;;;;;;;;;;;;あとで消す！！！！！！！！！！！！！！！！！！！！！！！！小暮

      Radnoconv(indwn : indwn + nwntemp - 1L, *) = Rtemp(inds, *)
      dRadnoconv(indwn : indwn + nwntemp - 1L, *) = deltaRtemp(inds, *)
      dRdalb(indwn : indwn + nwntemp - 1L, *) = dRdalbtemp(inds, *)
      dRdabs(indwn : indwn + nwntemp - 1L, *, *) = dRdtauabstemp(inds, *, *)
      dRdsca(indwn : indwn + nwntemp - 1L, *, *, *) = dRdtauscatemp(inds, *, *, *)

;      ;------- 2021/4/7 kogure Risei -----------------
;      path_output = pathnames.results
;      openw,lun,path_output + orbit + nameplusdir + '/' +'Radnoconv.txt',/get_lun
;      data=fltarr(2,nwn)
;      data(0,*)=wn & data(1,*)=transpose(Radnoconv(*,0))
;      printf,lun,data
;      close,lun
;      free_lun,lun
;      ;----------------------kogure end --------------

      IF ~convolveallinone THEN BEGIN
        FOR j = 0, nspectra - 1L DO BEGIN
          IF nwn_k(i) GE N_ELEMENTS(ILS) THEN BEGIN
            temp = CONVOL(Rtemp(*, j), ILS, /CENTER, /EDGE_WRAP) * (wn(1L) - wn(0))
            Radconv(indwn : indwn + nwntemp - 1L, j) = temp(inds)
            temp = CONVOL(deltaRtemp(*, j), ILS, /CENTER, /EDGE_WRAP) * (wn(1L) - wn(0))
            dRadconv(indwn : indwn + nwntemp - 1L, j) = temp(inds)
          ENDIF ELSE BEGIN
            Radconv(indwn : indwn + nwntemp - 1L, j) = Rtemp(inds, j)
            dRadconv(indwn : indwn + nwntemp - 1L, j) = deltaRtemp(inds, j)
          ENDELSE
        ENDFOR
      ENDIF

      ; Calculate the absorption Jacobians, convolve them and set into the final vectors
      IF ~nojacobian THEN BEGIN
        ind = DBLARR(2, nlayer)
        weights = DBLARR(2, nlayer)
        FOR j = 0, nlayer - 1L DO ind(*, j) = [j - 1L, j]
        ind(0, 0) = 0L
        weights(*) = .5D
        weights(0, 0) = 0D
        Jactemp = DBLARR(N_ELEMENTS(dRdtauabstemp(*, 0, 0)), nlayer, nspectra, 1L + ngas_fitted + naerosols_fitted * 3L);Jactemp(*, 大気層n, *, 0)=アルベド、（*,*,*,1）＝ガス、(*,*,*,2:3)=エアロゾル
        dJactemp = DBLARR(N_ELEMENTS(dRdtauabstemp(*, 0, 0)), nlayer, nspectra, 1L + ngas_fitted + naerosols_fitted * 3L)

        ; Jacobian to the albedo
        FOR ii = 0, nspectra - 1L DO Jactemp(*, 0, ii, 0) = dRdalbtemp(*, ii)
        FOR ii = 0, nspectra - 1L DO dJactemp(*, 0, ii, 0) = deltadRdalbtemp(*, ii)

        ; Jacobians to the gas VMR
        FOR k = 0, ngas_fitted - 1L DO BEGIN
          ACSloc = ACS[i].k(0 : ACS[i].nwn - 1L, *, onlyfitspecies(k))
          FOR ii = 0, nspectra - 1L DO BEGIN
            FOR j = 0, nlayer - 1L DO Jactemp(*, j, ii, k + 1L) = (weights(0, j) * dRdtauabstemp(*, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * dRdtauabstemp(*, ind(1L, j), ii) * dz(ind(1L, j))) * ACSloc(*, j)
            FOR j = 0, nlayer - 1L DO dJactemp(*, j, ii, k + 1L) = (weights(0, j) * deltadRdtauabstemp(*, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * deltadRdtauabstemp(*, ind(1L, j), ii) * dz(ind(1L, j))) * ACSloc(*, j)
          ENDFOR
        ENDFOR

        ; Jacobians to the aerosols VMR and particle radius
        FOR k = 0, naerosols_fitted - 1L DO BEGIN
          indloc = onlyfitspecies(ngas_fitted + k) - ngas
          abscoeff = DBLARR(N_ELEMENTS(wnloc), nlayer)
          FOR j = 0, nlayer - 1L DO abscoeff(*, j) = INTERPOL(Beta_aero[i].kabs(0 : Beta_aero[i].nwn - 1L, indloc, j), Beta_aero[i].wn(0 : Beta_aero[i].nwn - 1L), wnloc)
          scacoeff = DBLARR(N_ELEMENTS(wnloc), nlayer)
          FOR j = 0, nlayer - 1L DO scacoeff(*, j) = INTERPOL(Beta_aero[i].ksca(0 : Beta_aero[i].nwn - 1L, indloc, j), Beta_aero[i].wn(0 : Beta_aero[i].nwn - 1L), wnloc)
          dabscoeffdrm = DBLARR(N_ELEMENTS(wnloc), nlayer)
          FOR j = 0, nlayer - 1L DO dabscoeffdrm(*, j) = INTERPOL(Beta_aero[i].dkabsdrm(0 : Beta_aero[i].nwn - 1L, indloc, j), Beta_aero[i].wn(0 : Beta_aero[i].nwn - 1L), wnloc)
          dscacoeffdrm = DBLARR(N_ELEMENTS(wnloc), nlayer)
          FOR j = 0, nlayer - 1L DO dscacoeffdrm(*, j) = INTERPOL(Beta_aero[i].dkscadrm(0 : Beta_aero[i].nwn - 1L, indloc, j), Beta_aero[i].wn(0 : Beta_aero[i].nwn - 1L), wnloc)
          dabscoeffds = DBLARR(N_ELEMENTS(wnloc), nlayer)
          FOR j = 0, nlayer - 1L DO dabscoeffds(*, j) = INTERPOL(Beta_aero[i].dkabsds(0 : Beta_aero[i].nwn - 1L, indloc, j), Beta_aero[i].wn(0 : Beta_aero[i].nwn - 1L), wnloc)
          dscacoeffds = DBLARR(N_ELEMENTS(wnloc), nlayer)
          FOR j = 0, nlayer - 1L DO dscacoeffds(*, j) = INTERPOL(Beta_aero[i].dkscads(0 : Beta_aero[i].nwn - 1L, indloc, j), Beta_aero[i].wn(0 : Beta_aero[i].nwn - 1L), wnloc)

          FOR ii = 0, nspectra - 1L DO BEGIN
            ; VMR
            FOR j = 0, nlayer - 1L DO Jactemp(*, j, ii, ngas_fitted + k + 1L) = ((weights(0, j) * dRdtauabstemp(*, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * dRdtauabstemp(*, ind(1L, j), ii) * dz(ind(1L, j))) * abscoeff(*, j) + $
              (weights(0, j) * dRdtauscatemp(*, indloc, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * dRdtauscatemp(*, indloc, ind(1L, j), ii) * dz(ind(1L, j))) * scacoeff(*, j))
            FOR j = 0, nlayer - 1L DO dJactemp(*, j, ii, ngas_fitted + k + 1L) = ((weights(0, j) * deltadRdtauabstemp(*, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * deltadRdtauabstemp(*, ind(1L, j), ii) * dz(ind(1L, j))) * abscoeff(*, j) + $
              (weights(0, j) * deltadRdtauscatemp(*, indloc, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * deltadRdtauscatemp(*, indloc, ind(1L, j), ii) * dz(ind(1L, j))) * scacoeff(*, j))

             ;Aerosols mean radius
            FOR j = 0, nlayer - 1L DO Jactemp(*, j, ii, ngas_fitted + naerosols_fitted + k + 1L) = ((weights(0, j) * dRdtauabstemp(*, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * dRdtauabstemp(*, ind(1L, j), ii) * dz(ind(1L, j))) * dabscoeffdrm(*, j) + $
              (weights(0, j) * dRdtauscatemp(*, indloc, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * dRdtauscatemp(*, indloc, ind(1L, j), ii) * dz(ind(1L, j))) * dscacoeffdrm(*, j))
            FOR j = 0, nlayer - 1L DO dJactemp(*, j, ii, ngas_fitted + naerosols_fitted + k + 1L) = ((weights(0, j) * deltadRdtauabstemp(*, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * deltadRdtauabstemp(*, ind(1L, j), ii) * dz(ind(1L, j))) * dabscoeffdrm(*, j) + $
              (weights(0, j) * deltadRdtauscatemp(*, indloc, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * deltadRdtauscatemp(*, indloc, ind(1L, j), ii) * dz(ind(1L, j))) * dscacoeffdrm(*, j))
;            ; Aerosols mean radius  小暮粒径フィットテスト2021/11/26 こことget_K
;            FOR j = 0, nlayer - 1L DO Jactemp(*, j, ii, ngas_fitted + naerosols_fitted + k + 1L) = ((weights(0, j) * dRdtauabstemp(*, ind(0, j), ii) * (dz(ind(0, j))^2) + weights(1L, j) * dRdtauabstemp(*, ind(1L, j), ii) * (dz(ind(1L, j))^2))* abscoeff(*, j) * dabscoeffdrm(*, j) + $
;              (weights(0, j) * dRdtauscatemp(*, indloc, ind(0, j), ii) * (dz(ind(0, j))^2) + weights(1L, j) * dRdtauscatemp(*, indloc, ind(1L, j), ii) * (dz(ind(1L, j))^2)) * scacoeff(*, j) * dscacoeffdrm(*, j))
;            FOR j = 0, nlayer - 1L DO dJactemp(*, j, ii, ngas_fitted + naerosols_fitted + k + 1L) = ((weights(0, j) * deltadRdtauabstemp(*, ind(0, j), ii) * (dz(ind(0, j))^2) + weights(1L, j) * deltadRdtauabstemp(*, ind(1L, j), ii) * (dz(ind(1L, j))^2)) * abscoeff(*, j) * dabscoeffdrm(*, j) + $
;              (weights(0, j) * deltadRdtauscatemp(*, indloc, ind(0, j), ii) * (dz(ind(0, j))^2) + weights(1L, j) * deltadRdtauscatemp(*, indloc, ind(1L, j), ii) * (dz(ind(1L, j))^2)) * scacoeff(*, j) * dscacoeffdrm(*, j))

            ; Aerosols variance radius
            FOR j = 0, nlayer - 1L DO Jactemp(*, j, ii, ngas_fitted + 2 * naerosols_fitted + k + 1L) = ((weights(0, j) * dRdtauabstemp(*, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * dRdtauabstemp(*, ind(1L, j), ii) * dz(ind(1L, j))) * dabscoeffds(*, j) + $
              (weights(0, j) * dRdtauscatemp(*, indloc, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * dRdtauscatemp(*, indloc, ind(1L, j), ii) * dz(ind(1L, j))) * dscacoeffds(*, j))
            FOR j = 0, nlayer - 1L DO dJactemp(*, j, ii, ngas_fitted + 2 * naerosols_fitted + k + 1L) = ((weights(0, j) * deltadRdtauabstemp(*, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * deltadRdtauabstemp(*, ind(1L, j), ii) * dz(ind(1L, j))) * dabscoeffds(*, j) + $
              (weights(0, j) * deltadRdtauscatemp(*, indloc, ind(0, j), ii) * dz(ind(0, j)) + weights(1L, j) * deltadRdtauscatemp(*, indloc, ind(1L, j), ii) * dz(ind(1L, j))) * dscacoeffds(*, j))
          ENDFOR
        ENDFOR

        Jacnoconv(indwn : indwn + nwntemp - 1L, *, *, *) = Jactemp(inds, *, *, *)
        dJacnoconv(indwn : indwn + nwntemp - 1L, *, *, *) = dJactemp(inds, *, *, *)

        IF nwn_k(i) GE N_ELEMENTS(ILS) THEN BEGIN
          FOR j = 0, nspectra - 1L DO BEGIN
            FOR k = 0, nlayer - 1L DO BEGIN
              FOR ii = 0, 1L + (ngas_fitted + 3 * naerosols_fitted) - 1L DO BEGIN
                temp = CONVOL(Jactemp(*, k, j, ii), ILS, /CENTER, /EDGE_WRAP) * (wn(1L) - wn(0))
                Jacconv(indwn : indwn + nwntemp - 1L, k, j, ii) = temp(inds)
              ENDFOR
            ENDFOR
          ENDFOR
        ENDIF ELSE BEGIN
          Jacconv(indwn : indwn + nwntemp - 1L, *, *, *) = Jactemp(inds, *, *, *)
        ENDELSE

        dRdtauabstemp = dRdtauabstemp(inds, *, *)
        dRdtauscatemp = dRdtauscatemp(inds, *, *, *)
      ENDIF

      ; Update the wavenumber indices
      indwn = indwn + nwntemp - 1L
    ENDFOR
  ENDELSE

  IF convolveallinone THEN BEGIN
    inds = WHERE(wn GE wnrange(0, 0) AND wn LE wnrange(1L, nwindows - 1L))
    IF N_ELEMENTS(Radnoconv(*, 0)) GE N_ELEMENTS(ILS) THEN BEGIN
      FOR i = 0, nspectra - 1L DO BEGIN
        temp = CONVOL(Radnoconv(*, i), ILS, /CENTER, /EDGE_WRAP) * (wn(1L) - wn(0))
        Radconv(*, i) = temp(inds)
        temp = CONVOL(dRadnoconv(*, i), ILS, /CENTER, /EDGE_WRAP) * (wn(1L) - wn(0))
        dRadconv(*, i) = temp(inds)
        IF ~nojacobian THEN BEGIN
          FOR j = 0, nlayer - 1L DO BEGIN
            FOR k = 0, 1L + (ngas_fitted + 3 * naerosols_fitted) - 1L DO BEGIN
              temp = CONVOL(Jacnoconv(*, j, i, k), ILS, /CENTER, /EDGE_WRAP) * (wn(1L) - wn(0))
              Jacconv(*, j, i, k) = temp(inds)
            ENDFOR
          ENDFOR
        ENDIF
      ENDFOR
    ENDIF ELSE BEGIN
      Radconv = Radnoconv
      dRadconv = dRadnoconv
      Jacconv = Jacnoconv
    ENDELSE
  ENDIF
  Radnoconv = Radnoconv(inds, *)
  dRadnoconv = dRadnoconv(inds, *)
  Jacnoconv = Jacnoconv(inds, *, *, *, *)
  dJacnoconv = dJacnoconv(inds, *, *, *, *)
  Jacconv = Jacconv(inds, *, *, *, *)

  RETURN, Radconv
END