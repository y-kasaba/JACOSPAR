FUNCTION UPDATE_JACOSPAR_INPUT_PARALLEL_Limb_ShoheiDatabase, seq, $
  outputfile, $
  fov, $
  npf, $
  nang, $
  precision, $
  seed, $
  ngroups, $
  nmembers, $
  naerosols, $
  ngas, $
  pathname, $
  altrange, $
  z, $
  NOJACOBIAN = nojacobian, $
  ONLYNAMES = onlynames, $
  NOFOV = nofov, $
  NOMULTSCAT = nomultscat, $
  GEOMETRY = geom
  path_phs = path_phs

  IF ~KEYWORD_SET(nojacobian) THEN nojacobian = 0
  IF ~KEYWORD_SET(onlynames) THEN onlynames = 0
  IF ~KEYWORD_SET(nofov) THEN nofov = 0
  IF ~KEYWORD_SET(nomultscat) THEN nomultscat = 0

  ; Function to upadte the jacospar input file for a given orbit
  ;
  ; Inputs:
  ;         - orbit: orbit number [INT or STRING]
  ;         - seq: the sequence of measurements that should be considered [INT]
  ;         - outpufile: the name of the outputfile [STRING]
  ;         - npf: number of tabulated phase functions in the database [INT]
  ;         - nang: number of tabulated angles in the database [INT]
  ;         - filenameOptProp: filename of the file containing the optical properties [STRING]
  ;         - filenamePhsProp: filename of the file containing the phase function properties [STRING]
  ;         - precision: mean rms in the data to be fitted [FLOAT]
  ;         - dwn_sc: the wavenumber step on which the scattering properties are given [FLOAT]
  ;         - wn_ac: the wavenumber grid on which the absorption properties are given [FLTARR]
  ;         - pathname: the path where the file should be written [STRING]
  ;
  ; Output: File is written
  ;
  pathnames = PATH_MANAGEMENT()

  ; Load the orbit geometry
  IF ~KEYWORD_SET(geom) THEN BEGIN
    PRINT, "Please define the geometry file."
  ENDIF
  nlayer = N_ELEMENTS(z)
  nspectra = N_ELEMENTS(geom(0, *))
  Rmars = 3.3962D3

  ; Write the input file
  filename = STRARR(nspectra)

  ; Calculation of the pointing and solar angles, based on note from Mayorov, 2006
  ; lat, lon, LT, Ls, inc, em, phase, altimetry, slant_dist, dphi, thesol, Mex alt
  ; (angles in radian, distances in km)
;    zobs= dblarr(n_elements(geom(0,*)));;kogure to calculate nadir 2020/5/28 ;;geom(0, *)
;    zobs(*) =  0d; + Rmars;kogure to calculate nadir 2020/5/28
  zobs= geom(0, *)
  hd = geom(1, *)
;    hd = dblarr(n_elements(geom(1, *)))
;    hd(*)=300d; + Rmars
;    thetaP = dblarr(n_elements(geom(0,*)));;kogure to calculate nadir 2020/5/28 
;    thetaP(*) = 180d;;kogure to calculate nadir 2020/5/28 theP
  thetaP = ACOS((zobs + Rmars) / (hd + Rmars)) / !DPI * 180D ; angle VOP;;2021/4/21 kogure
  dphi = ABS(geom(2, *))
;    dphi = dblarr(n_elements(geom(0,*)))
;    dphi(*) = 0d;相対位相角
  thesol = geom(3, *)
;    thesol = dblarr(n_elements(geom(3, *)))
;    thesol(*)=30d;zenith

  ;precision=0.3d

  ind = 0
  FOR i = 0, nspectra - 1 DO BEGIN
    filename(ind) = pathname + outputfile + "_" + STRING(i, FORMAT = "(I02)")
    filenameOptProp = pathname + outputfile + ".opt"
    filenamePhsProp =pathname + outputfile + '.phs';pathname + outputfile + ".phs";2021/2/12 kogure '/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Aerosols/jacospar_phs_test_4100_4200.phs'
    IF ~onlynames THEN BEGIN
      OPENW, lun, filename(ind) + ".inp", /get_lun

      PRINTF, lun, "&jacospar_env"
      PRINTF, lun, "  logfile = '" + filename(ind) + ".log'"
      PRINTF, lun, "  mverb = 3";2021/10/7 kogure risei change from 1 to 3
      PRINTF, lun, "/"
      PRINTF, lun, ""

      PRINTF, lun, "&jacospar_sys"
      PRINTF, lun, "  nsm_i = " + STRING(nmembers, FORMAT = "(I5)") + "    ! # of system members per GS in the input"
      PRINTF, lun, "  ngs_i = " + STRING(ngroups, FORMAT = "(I5)") + "    ! # of groups of systems (GSs) in the input"
      PRINTF, lun, "  ncom = " + STRING(naerosols, FORMAT = "(I2)") + "       ! # of scattering components (should be > 0)"
      PRINTF, lun, "  nlay = " + STRING(nlayer, FORMAT = "(I4)") + "     ! # of layers (should be > 0)"
      PRINTF, lun, "  nphs = " + STRING(npf, FORMAT = "(I6)") + "   ! # of tabulated phase functions in the database" ;2021/2/12 kogure change to 200
      PRINTF, lun, "  nang_i = 181   ! # of tabulated angles in the database"; + STRING(nang, FORMAT = "(I4)") + "   ! # of tabulated angles in the database";2021/2/12 kogure change to 181
      texte = ""
      FOR j = 0, naerosols - 2 DO texte = texte + '-1, '
      texte = texte + '-1'
      PRINTF, lun, "  mphs(1:2)=1,1! method for specifying phase function" ;+ STRING(naerosols, FORMAT = "(I2)") + ") = " + texte + "! method for specifying phase function";kogure change to :2)=1,1
      PRINTF, lun, "  msph = 1        ! flag for spherical shell geometry (0:no, 1:yes)"
      PRINTF, lun, "  rpla = " + STRING(Rmars * 1D3, FORMAT = "(F12.2)") + " ! radius (m) of the sphere with zero height"
      PRINTF, lun, "  hbot = " + STRING(z(0) * 1D3, "(F7.1)") + " ! altitude (m) of bottom of atmosphere (BOA)"
      stringalt = ''
      FOR j = 1, N_ELEMENTS(z) - 1 DO stringalt = stringalt + STRING(z(j) * 1D3, "(F7.1)") + ", "
      stringalt += " 100000.0"
      PRINTF, lun, "  hgrd = " + stringalt + " ! altitude (m) of the layer top boundary"
      PRINTF, lun, "  optfile = '" + filenameOptProp + "' ! file name for optical properties"
      PRINTF, lun, "  phsfile = '" + filenamePhsProp + "' ! file name for phase functions"
      PRINTF, lun, "  naps_i  = 0 ! # of atmosphere perturbations in the input"
      PRINTF, lun, "/"
      PRINTF, lun, ""

      PRINTF, lun, "&jacospar_expr"
      PRINTF, lun, "  jseed = " + STRING(seed, FORMAT = "(I2)") + " ! a seed for random number generator [1,2**31-1] (0 for automatic)"
      PRINTF, lun, "  njob = " + STRING(nmembers, FORMAT = "(I5)") + " ! # of computational jobs"
      PRINTF, lun, "  ngs = " + STRING(ngroups, FORMAT = "(I5)") + "  ! # of system groups with common scattering properties"
      PRINTF, lun, "  nsm   = 1   ! # of member systems per group"
      PRINTF, lun, "  nsrc  = 1   ! # of radiation sources"
      PRINTF, lun, "  mthm  = 0   ! flag of thermal sources (0 or 1)"
      PRINTF, lun, "  mdtmp = 0   ! flag for Jtmp (0:off, 1:on)"
      PRINTF, lun, "  mdalb = " + STRING(~nojacobian, FORMAT = "(I1)") + "   ! flag for Jalb (0:off, 1:on)"
      PRINTF, lun, "  mdabs = " + STRING(~nojacobian, FORMAT = "(I1)") + "   ! flag for Jabs (0:off, 1:on)"
      IF ~nojacobian THEN PRINTF, lun, "  ndsc  = " + STRING(naerosols, FORMAT = "(I1)") + "   ! # of components in [1, ncom] for Jsca_com" ELSE $
        PRINTF, lun, "  ndsc  = 0   ! # of components in [1, ncom] for Jsca_com"
      PRINTF, lun, "  nsps  = 0   ! # of surface perturbations"
      PRINTF, lun, "  naps  = 0   ! # of atmosphere perturbations"
      PRINTF, lun, "  mpplk = 0   ! flag [0/1] of Plank function perturbations"
      IF ~nojacobian THEN PRINTF, lun, "  idmin = 1   ! Minimum index for Jacobian calculation"
      IF ~nojacobian THEN PRINTF, lun, "  idmax = " + STRING(nlayer, FORMAT = "(I4)") + "   ! Maximum index for Jacobian calculation"
      PRINTF, lun, "  msolv = 1   ! -1:RTRN, 0:MC, 1:MC+analytical"
      PRINTF, lun, "  mtech = 2   ! flag [0,3] for automatic configuration of technical parameters" ; changed from 2 to 0 (AM, Shohei's email from 19/10/2017, 11/11/2017)
      PRINTF, lun, "  outfile = '" + filename(ind) + ".res'"
      PRINTF, lun, "/"
      PRINTF, lun, ""

      PRINTF, lun, "&jacospar_tech"
      IF nomultscat THEN PRINTF, lun, "  iocfor = 1        ! # of forced scattering events" ELSE PRINTF, lun, "  iocfor = 5        ! # of forced scattering events"
      IF nomultscat THEN PRINTF, lun, "  iocmax = 1        ! max order of scattering events" ELSE PRINTF, lun, "  iocmax = 100     ! max order of scattering events"
      PRINTF, lun, "  dtaufor = 0.1d   ! min layer optical thickness for collision-forcing of methodII";Add this parameter. kogure risei 2021/10/06 (Dr.Iwabuchi's and Shohei's email from 06/10/2021)
      PRINTF, lun, "frcsca0 = -0.8       ! parameter to set scattering coefficient in BS";岩渕先生メール 2021/11/15 このパラメータを加える。
      PRINTF, lun, "wmin_IS = 0.01 ! min. spectral weight";岩渕先生メール 2021/11/15 このパラメータを加える。
      PRINTF, lun, "/"
      PRINTF, lun, ""

      PRINTF, lun, "&jacospar_job"
      PRINTF, lun, "  ism_s = 1        ! location of the start line read"
      PRINTF, lun, "  nbat_min = 3     ! min # of batches of photons"; Change from 2 to 5. kogure risei 2021/10/06 (Dr.Iwabuchi's and Shohei's email from 06/10/2021)
      PRINTF, lun, "  nbat_max = 1000   ! max # of batches of photons"; Change from 100 to 1000. kogure risei 2021/10/06 (Dr.Iwabuchi's and Shohei's email from 06/10/2021)
      PRINTF, lun, "  npho = 10000      ! # of photons per batch"
      PRINTF, lun, "  macc = 1         ! flag for accuracy kind (0:mean, 1:rms, 2:max, 3:mean, 4:rms, 5:max)"
      PRINTF, lun, "  acc  = " + STRING(precision, FORMAT = "(F8.3)") + "      ! required relative RMSE (%)"
      PRINTF, lun, "  iblo = 1         ! data block number"
      PRINTF, lun, "  igs_s = 1        ! location of the start point read"
      PRINTF, lun, "  igs_w = 1        ! spacing in read points"
      PRINTF, lun, "  ism_s = 1        ! location of the start line read"
      PRINTF, lun, "  ism_w = 1        ! spacing in read lines"
      PRINTF, lun, "  hgtV = " + STRING(hd(i) * 1D3, FORMAT = "(F9.1)") + " ! height (m) at the view point V"
      PRINTF, lun, "  hgtP = " + STRING(zobs(i) * 1D3, FORMAT = "(F9.1)") + " ! height (m) at the target point P"
      PRINTF, lun, "  mlos = 2         ! method for defining plos (0:theV, 1:theP, 2:angPOV, 3:hgtH)";PRINTF, lun, "  mlos = 2         ! method for defining plos (0:theV, 1:theP, 2:angPOV, 3:hgtH)";;2021/4/19
      PRINTF, lun, "  plos = " + STRING(thetaP(i), FORMAT = "(F11.7)") + " ! theP"
      PRINTF, lun, "  mgeoH = 1        ! tangent point geometry, 1:forward(VPH/PVH), -1:backward(PHV/VHP)"
      PRINTF, lun, "  phiobs = " + STRING(180D, FORMAT = "(F11.7)") + " ! azimuth angle (deg.) of V-to-P vector";;;;;;;;;;;;;;;;;;;;;;change to nadir parameter  original =180d
      ;PRINTF, lun, "  phiobs = " + STRING(-35D, FORMAT = "(F11.7)") + " ! azimuth angle (deg.) of V-to-P vector";;;;;;;;;;;;;;;;;;;;;;change to nadir parameter  original =180d
      IF nofov THEN PRINTF, lun, "  conFOV = 0.0     ! half cone angle (deg.) of FOV" $
      ELSE PRINTF, lun, "  conFOV = " + STRING(fov, FORMAT = "(F11.7)")  + "! half cone angle (deg.) of FOV"
      PRINTF, lun, "  maa = 2          ! angular averaging of radiance (0:infinitesimal, 1:S-average, 2:P-average)"
      PRINTF, lun, "  msrck  = 1       ! radiation source kind (1=solar, 2=thermal)"
      PRINTF, lun, "  thesol = " + STRING(thesol(i), FORMAT = "(F11.7)") + "   ! zenith  angle (deg.) of solar direction"
      PRINTF, lun, "  phisol = " + STRING(dphi(i), FORMAT = "(F11.7)") + "   ! azimuth angle (deg.) of solar direction"
      PRINTF, lun, "/"
      PRINTF, lun, ""

      FOR j = 2, nmembers DO BEGIN
        PRINTF, lun, "&jacospar_job"
        PRINTF, lun, " ism_s = " + STRING(j, FORMAT = "(I4)") + "       ! location of the start line read"
        PRINTF, lun, "/"
      ENDFOR
      FREE_LUN, lun
    ENDIF
    ind = ind + 1
  ENDFOR
  RETURN, filename(0 : ind - 1)
END
