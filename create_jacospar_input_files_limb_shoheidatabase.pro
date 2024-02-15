FUNCTION CREATE_JACOSPAR_INPUT_FILES_Limb_ShoheiDatabase, orbit, $
  wnrange, $
  gas_species, $
  aerosols_species, $
  alt_range, $
  dwn, $
  FOV, $
  z, $
  ztpn, $
  Ht, $
  geometry, $
  Radiusm, $
  Radiuss, $
  wn_albedo, $
  S_albedo, $

  GAS_ACS = Gas_ACS, $
  AERO_K_ABS = Aero_k_abs, $
  AERO_K_SCA = Aero_k_sca, $
  DAERODRM_K_ABS = dAerodRm_k_abs, $
  DAERODRM_K_SCA = dAerodRm_k_sca, $
  DAERODS_K_ABS = dAerodS_k_abs, $
  DAERODS_K_SCA = dAerodS_k_sca, $

  SAVENAME = savename, $
  NOABS = noabs, $

  PATH_OUTPUT = path_output, $
  PRECISION = precision, $
  SEED = seed, $
  wn = wnaero, $

  GET_TAU_SCA = get_tau_sca, $
  GET_TAU_ABS = get_tau_abs, $
  TAU_AERO_SCA = tau_aero_sca, $
  TAU_ABS = tau_abs, $
  QEXT = Qext, $

  DELTATAUABS = deltatauabs, $
  DELTATAUSCA = deltatausca, $

  DWNAERO = dwnaero, $
  PLOTTOTALOPTICALDEPTH = plottotalopticaldepth, $

  _EXTRA = extra

  ngas = N_ELEMENTS(gas_species)
  naerosols = N_ELEMENTS(aerosols_species)

  IF ~KEYWORD_SET(noabs) THEN noabs = 0L
  IF ~KEYWORD_SET(notdotauintegration) THEN notdotauintegration = 0L
  IF ~KEYWORD_SET(precision) THEN precision = 1D
  IF ~KEYWORD_SET(seed) THEN seed = 1L
  IF ~KEYWORD_SET(get_tau_sca) THEN get_tau_sca = 0L
  IF ~KEYWORD_SET(get_tau_abs) THEN get_tau_abs = 0L
  IF ~KEYWORD_SET(dwnaero) THEN dwnaero = dwn
  IF ~KEYWORD_SET(plottotalopticaldepth) THEN plottotalopticaldepth = 0

  ; -------------------
  ; Path and file names
  ; -------------------
  pathnames = PATH_MANAGEMENT()
  CD, CURRENT = home
  IF ~KEYWORD_SET(path_output) THEN path_output = pathnames.results ; where the results will be written
  IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output

  IF ~KEYWORD_SET(savename) THEN savename = "Orbit_" + orbit

  ; ----------
  ; Parameters
  ; ----------
  n_angles = 181 ; Number of angles for the phase functions
  DEG2RAD = !pi / 180D
  kB = 1.38D-16 ; [cgs]
  gM = 371.1D ; [cgs]
  epsilon = 1D-2 ; To calculate the Jacobians to the extinction of the gas species
  nspectra = N_ELEMENTS(geometry(0, *))
  nlayer = N_ELEMENTS(z)
  dz = z(1 : nlayer - 1) - z(0 : nlayer - 2)
  nu_ref = 8500D

  ; -----------------------------------------------------
  ; Check format of entry variables, all should be double
  ; -----------------------------------------------------
  IF NOT (ISA(dwn, "DOUBLE") AND ISA(z, "DOUBLE") AND ISA(alt_range, "DOUBLE")) THEN $
    MESSAGE, "The variables 'z', 'dwn', 'dwn' and 'alt_range' must be double precision numbers"

  ; -----------------------------
  ; Defining the wavenumber grids
  ; -----------------------------
  wn = dwn * DINDGEN(ROUND((wnrange(1) - wnrange(0)) / dwn) + 1L) + wnrange(0)
  N_wn = N_ELEMENTS(wn)
  wnaero = dwnaero * DINDGEN(ROUND((wnrange(1) - wnrange(0)) / dwnaero) + 1L) + wnrange(0)
  N_wn_aero = N_ELEMENTS(wnaero)

  ; ------------------------------------------
  ; Absorption coefficients of the gas species
  ; ------------------------------------------
  ngroups = N_wn_aero - 1L;波数の数 number of wavenumber
  nmembers = ROUND((N_wn - 1D) / (N_wn_aero - 1D));波数とエアロゾルの波数の比 ratio of wavenumber
  IF nmembers * ngroups + 1 NE N_wn THEN BEGIN;比が１じゃないとき
    wn = wn(0 : nmembers * ngroups)
    N_wn = N_ELEMENTS(wn)
  ENDIF
  IF ~KEYWORD_SET(Gas_ACS) THEN BEGIN
    Gas_ACS = DBLARR(N_wn, nlayer, ngas)
    FOR i = 0, ngas - 1 DO BEGIN
      k_temp = INTERPOL_ACS(gas_species(i), [wnrange(0), wnrange(1)], wn_temp, ztpn(2, *), ztpn(1, *), FileName, dwn)
      wntemp = wn(0) + DINDGEN(N_wn) * dwn
      IF ROUND(ALOG10((wn_temp(1) - wn_temp(0)) / (wn(1) - wn(0)))) NE 0D THEN k_temp = RESAMPLE(wn_temp, k_temp, wntemp, dwn)
      IF N_ELEMENTS(k_temp) EQ 0 THEN BEGIN
        MESSAGE, "The pressure and/or temperature are out of range for the KDB files."
      ENDIF
      Gas_ACS(*, *, i) = k_temp
    ENDFOR
    DELVAR, i, inds, wn_temp, filename
  ENDIF

  ; -------------------------------
  ; Calculate tau for all the gases
  ; -------------------------------
  IF ngas GT 0 THEN BEGIN
    tau_gas_abs = DBLARR(nlayer, N_wn, ngas);gas吸収の配列（高度、波数、ガス種類）
    dzloc = 0.01D
    FOR i = 0, ngas - 1 DO BEGIN
      CASE STRLOWCASE(gas_species(i)) OF
        'co2': column = 5
        'co': column = 6
        'h2o': column = 7
      ENDCASE
      FOR j = 0, nlayer - 2L DO tau_gas_abs(j, *, i) = (Gas_ACS(*, j, i) * ztpn(column, j) + Gas_ACS(*, j + 1, i) * ztpn(column, j + 1)) * dz(j) * .5D3
      tau_gas_abs(nlayer - 1L, *, i) = Gas_ACS(*, nlayer - 1L, i) * Ht(nlayer - 1L) * ztpn(column, nlayer - 1L) * 1D3
    ENDFOR
    DELVAR, column
  ENDIF ELSE tau_gas_abs = DBLARR(nspectra - 1, N_wn)

  ; -----------------------------
  ; Angles for the phase function
  ; -----------------------------
  angle = DINDGEN(n_angles) * 180D / DOUBLE(n_angles - 1)
  cosangle = COS(angle * DEG2RAD)

  ; ------------------------
  ; Read aerosols parameters
  ; ------------------------
  angle = INDGEN(181, START = 0, /DOUBLE) ; angles from 0 to 180° by steps of 1°
  nangles = N_ELEMENTS(angle)
  Dqv = cos(angle / !RADEG)

  Aero_k_sca = DBLARR(N_wn_aero, naerosols, nlayer)
  Aero_k_abs = DBLARR(N_wn_aero, naerosols, nlayer)
  dAerodRm_k_sca = DBLARR(N_wn_aero, naerosols, nlayer)
  dAerodRm_k_abs = DBLARR(N_wn_aero, naerosols, nlayer)
  dAerodS_k_sca = DBLARR(N_wn_aero, naerosols, nlayer)
  dAerodS_k_abs = DBLARR(N_wn_aero, naerosols, nlayer)
  Phasefunc = DBLARR(N_wn_aero, n_angles, naerosols, nlayer, /NOZERO)
  Phasefunc(*) = 0;;2021/2/11 1→0
  ;g = DBLARR(N_wn_aero, naerosols, nlayer)
  radiusPath = strarr(naerosols, nlayer)

  FOR i = 0, naerosols - 1 DO BEGIN
    CASE STRLOWCASE(aerosols_species(i)) OF
      'dust': column = 8
      'waterice': column = 9
    ENDCASE
    FOR j = 0, nlayer - 1 DO BEGIN
      IF ztpn(column, j) GT 0 THEN BEGIN
        data = INTERPOL_AEROSOLS_OPTICAL_ShoheiDatabase(aerosols_species(i), Radiusm(i, j), Radiuss(i, j), wnrange = wnrange, wnstep = dwn, _EXTRA = extra);;kogure 20210810 To CREATE_JACOSPAR_INPUT_FILES_Limb_ShoheiDatabase 20211111 Add "wnrange=wnrange, wnstep=dwn" Risei kogure 
        Aero_k_abs(*, i, j) = INTERPOL(data(0:-2, 1), data(0:-2, 0), wnaero)
        Aero_k_sca(*, i, j) = INTERPOL(data(0:-2, 2), data(0:-2, 0), wnaero)
        dAerodRm_k_abs(*, i, j) = INTERPOL(data(0:-2, 3), data(0:-2, 0), wnaero) * 1D4;cm2/um *1d4 = cm2/cm  データベースにはdQ/dR=cm2/cmと書いてあるが、実際はcm2/umであるため、cm2/cmに直す
        dAerodRm_k_sca(*, i, j) = INTERPOL(data(0:-2, 4), data(0:-2, 0), wnaero) * 1D4;cm2/um *1d4 = cm2/cm   データベースにはdQ/dR=cm2/cmと書いてあるが、実際はcm2/umであるため、cm2/cmに直す
        dAerodS_k_abs(*, i, j) = INTERPOL(data(0:-2, 5), data(0:-2, 0), wnaero)
        dAerodS_k_sca(*, i, j) = INTERPOL(data(0:-2, 6), data(0:-2, 0), wnaero)
        ;g(*, i, j) = INTERPOL(data(*, 7), data(*, 0), wnaero)      ;;kogure 20210810 散乱位相関数を用いるので不要。青木さんデータベースにはない。 
        radiusPath(i, j) = STRLOWCASE(aerosols_species(i)) + '_' + STRING(data(-1,0), "(F5.3)") + '00_' + STRING(data(-1,1), "(F5.3)");;;;;;;;;ここまで20210811あとでチェック
      ENDIF
    ENDFOR
  ENDFOR
  
  ;;;;Risei . phs file 20210817 Shoheidatabase ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;N_wn_aero=n_elements(phs(0,*))
  phasefunc_dust=fltarr(181, nlayer,N_wn_aero)
  phasefunc_ice=fltarr(181, nlayer,N_wn_aero)
  phs_dust_fin=fltarr(181,N_wn_aero)
  phs_ice_fin=fltarr(181,N_wn_aero)
  
  
  for j = 0, nlayer - 1 DO BEGIN
    IF ztpn(8, j) GT 0 THEN BEGIN
      phs_dustpre=read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Aerosols_phs/' + radiuspath(0,j) + '_phs.dat',COMMENT_SYMBOL='%');,data_start=1)
      phs_dustpre=phs_dustpre.field001
      phs_dust=fltarr(182,10001)
      phs_ice=fltarr(182,10001)
      phs_icepre=read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Aerosols_phs/' + radiuspath(1,j) + '_phs.dat',COMMENT_SYMBOL='%')
      phs_icepre=phs_icepre.field001
      for i=0,10000 do begin
        phs_dust(0:180,i)=phs_dustpre(*,i*2)
        phs_dust(181,i)=phs_dustpre(0,i*2+1)
        phs_ice(0:180,i)=phs_icepre(*,i*2)
        phs_ice(181,i)=phs_icepre(0,i*2+1)
      endfor
      wn_phs=phs_dust(0,*)
      wnind=intarr(N_wn_aero)
      for i=0,N_wn_aero-1 do begin
        near=Min((1d/wn_phs*1d4 - wn(i))^2, index)
        wnind(i)=index
      endfor

      phs_dust_fin(*)=phs_dust(1:-1,wnind)
      phs_ice_fin(*)=phs_ice(1:-1,wnind)

      phasefunc_dust(*, j,*) = phs_dust_fin
      phasefunc_ice(*, j,*) = phs_ice_fin
    ENDIF
  endfor
  ;
  ;Grit of angles (0-180 deg. with 1 deg interval)
  np=n_elements(phasefunc_dust(*,0,0))
  
  ang = fltarr(np)
  ang0 = 180./double(np-1)
  ang(0) = 0.d
  for i = 1, np-1 do ang(i) = ang(i-1) + ang0
  openw,lun,path_output + savename + '.phs',/get_lun
  writeu,lun,ang
  ;angles def (0-180 deg. with 1 deg interval)
  for j = 0, nlayer-1 do begin
    for i = 0l, N_wn_aero-2 do begin ;nsc_sca = number of wn grid
      writeu,lun,Phasefunc_dust(*,j,i) ;phase function of dust (181 points at ith wavenubemr)
      writeu,lun,Phasefunc_ice(*,j,i);phase function of ice (181 points at ith wavenubemr)
    endfor
  endfor
  free_lun,lun
  ;;;;Risei;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  ;;;;;read DISORT aerosol file:: kogure 20200326;;;;;;;;;;;;;;;;;
;  IF ztpn(8, 0) GT 0 and ztpn(9, 0) GT 0 THEN BEGIN
;    path_DISORT=pathnames.aerosols + 'Dust_gamma_Kleinbohl.p48'
;    path_disort_ice=pathnames.aerosols + 'waterice_warren_modified_gamma.p48'
;
;    data_DISORT=read_ascii(path_DISORT)
;    data_DISORT=data_DISORT.field01
;    data_DISORT_ice=read_ascii(path_DISORT_ice)
;    data_DISORT_ice=data_DISORT_ice.field01
;
;    wvn_DISORT=data_DISORT(0,*)
;    ext_DISORT=data_DISORT(1,*)*1e-8
;    ext_DISORT_ice=data_DISORT_ice(1,*)*1e-8
;    dust_alb=data_DISORT(2,*)
;    ice_alb=data_DISORT_ice(2,*)
;
;    sca_dust=ext_DISORT*dust_alb
;    abs_dust=ext_DISORT-sca_dust
;    sca_ice=ext_DISORT_ice*ice_alb
;    abs_ice=ext_DISORT_ice-sca_ice
;    for j=0,nlayer - 1 DO BEGIN
;      Aero_k_abs(*, 0, j) = INTERPOL(abs_dust, wvn_DISORT, wnaero)
;      Aero_k_sca(*, 0, j) = INTERPOL(sca_dust, wvn_DISORT, wnaero)
;      Aero_k_abs(*, 1, j) = INTERPOL(abs_ice, wvn_DISORT, wnaero)
;      Aero_k_sca(*, 1, j) = INTERPOL(sca_ice, wvn_DISORT, wnaero)
;    endfor
;  endif
;
;  ;;;;;;;;;;;;;;;;;;;;;

  ; ----------------------------------------
  ; Calculate the aerosols vertical profiles
  ; ----------------------------------------
  tau_aero_sca = DBLARR(nlayer, N_wn_aero, naerosols)
  tau_aero_abs = DBLARR(nlayer, N_wn_aero, naerosols)
  FOR i = 0, naerosols - 1 DO BEGIN
    CASE STRLOWCASE(aerosols_species(i)) OF
      'dust': column = 8
      'waterice': column = 9
    ENDCASE
    FOR j = 0, nlayer - 2 DO BEGIN ; Aero_k_sca : cm2 ; ztpn(3,j): cm-3; z * 1e5: cm
      tau_aero_sca(j, *, i) = (Aero_k_sca(*, i, j) * ztpn(column, j) + Aero_k_sca(*, i, j + 1) * ztpn(column, j + 1)) * dz(j) * .5D5
      tau_aero_abs(j, *, i) = (Aero_k_abs(*, i, j) * ztpn(column, j) + Aero_k_abs(*, i, j + 1) * ztpn(column, j + 1)) * dz(j) * .5D5
    ENDFOR
    tau_aero_sca(nlayer - 1L, *, i) = Aero_k_sca(*, i, j) * ztpn(column, nlayer - 1L) * Ht(nlayer - 1L) * 1D5
    tau_aero_abs(nlayer - 1L, *, i) = Aero_k_abs(*, i, j) * ztpn(column, nlayer - 1L) * Ht(nlayer - 1L) * 1D5
  ENDFOR
  DELVAR, column
  ;-------------- kogure Risei 2021/1/7 write down gas & aerosol absorption---------------------------------------------
  openw,lun,path_output+'Gas_ACS',/get_lun;;;;kogure 20200505
  printf,lun,transpose(Gas_ACS(*,0,*))
  free_lun,lun

;  openw,lun,path_output+'Gas_ACS_40km',/get_lun;;;;kogure 20200505
;  printf,lun,transpose(Gas_ACS(*,40,*))
;  free_lun,lun
;
;  openw,lun,path_output+'Gas_ACS_60km',/get_lun;;;;kogure 20200505
;  printf,lun,transpose(Gas_ACS(*,60,*))
;  free_lun,lun

  openw,lun,path_output+'tau_gas_abs',/get_lun;;;;kogure 20200505
  printf,lun,transpose(tau_gas_abs(0,*,*))
  free_lun,lun

  openw,lun,path_output+'ztpn',/get_lun;;;;kogure 20200505
  printf,lun,ztpn
  free_lun,lun;;;;;;;;;;;;;;;;;kogure 20200505

  openw,lun,path_output + savename +'Aero_k_abs.txt',/get_lun
  data1=dblarr(3,N_wn_aero)
  data1(0,*)=wn
  data1(1,*)=Aero_k_abs(*,0,0)
  data1(2,*)=Aero_k_abs(*,1,0)
  printf,lun,['wn','dust_cross','ice_cross']
  printf,lun,data1
  close,lun
  free_lun,lun

  openw,lun,path_output + savename +'Aero_k_sca.txt',/get_lun
  data2=dblarr(3,N_wn_aero)
  data2(0,*)=wn
  data2(1,*)=Aero_k_sca(*,0,0)
  data2(2,*)=Aero_k_sca(*,1,0)
  printf,lun,['wn','dust_cross','ice_cross']
  printf,lun,data2
  close,lun
  free_lun,lun

  openw,lun,path_output + savename +'tau_aero_abs.txt',/get_lun
  data3=dblarr(3,N_wn_aero)
  data3(0,*)=wn
  data3(1,*)=tau_aero_abs(0,*,0)
  data3(2,*)=tau_aero_abs(0,*,1)
  printf,lun,['wn','dust_abs','ice_abs']
  printf,lun,data3
  close,lun
  free_lun,lun

  openw,lun,path_output + savename +'tau_aero_sca.txt',/get_lun
  data4=dblarr(nlayer,naerosols)
  if naerosols eq 2 then begin
    data4=dblarr(3,N_wn_aero)
    data4(0,*)=wn
    data4(1,*)=tau_aero_sca(0,*,0)
    data4(2,*)=tau_aero_sca(0,*,1)
    printf,lun,['wn','dust_sca','ice_sca']
    printf,lun,data4
    close,lun
    free_lun,lun
  endif
  if naerosols eq 1 then begin
    printf,lun,'wavenumber=', wnaero(0)
    printf,lun,'dust tau sca'
    printf,lun,tau_aero_sca(*,0,0)
  endif
;  data4(0,*)=wn
;  data4(1,*)=tau_aero_sca(*,0,0)
;  data4(2,*)=tau_aero_sca(*,0,1)
;  printf,lun,['tau abs','ice_abs']
;  printf,lun,data4
  close,lun
  free_lun,lun

;  openw,lun,path_output + savename +'g.txt',/get_lun
;  data5=dblarr(3,N_wn_aero)
;  data5(0,*)=wn
;  data5(1,*)=g(*,0,0)
;  data5(2,*)=g(*,1,0)
;  printf,lun,['wn','dust','ice']
;  printf,lun,data5
;  close,lun
;  free_lun,lun
  ;-------------- end   Risei---------------------------------------------------------

  ; ------------------------------
  ; Calculate the local extinction
  ; ------------------------------
  Qext = DBLARR(nlayer, N_wn_aero, naerosols)
  ind = WHERE(wnaero EQ nu_ref, /NULL)
  FOR i = 0, naerosols - 1 DO BEGIN
    CASE STRLOWCASE(aerosols_species(i)) OF
      'dust': column = 8
      'waterice': column = 9
    ENDCASE
    FOR j = 0, nlayer - 1 DO Qext(j, *, i) = (Aero_k_sca(*, i, j) + Aero_k_abs(*, i, j)) * (ztpn(column, j))
  ENDFOR

  ; ----------------------------------
  ; Calculate the total tau absorption
  ; ----------------------------------
  IF ngas EQ 0 THEN BEGIN
    IF naerosols EQ 1 THEN BEGIN
      tau_abs = INTERPOL(tau_aero_abs, wnaero, wn)
    ENDIF ELSE BEGIN
      tau_abs = DBLARR(nlayer, N_wn)
      FOR i = 0, nlayer - 1L DO tau_abs(i, *) = INTERPOL(TOTAL(tau_aero_abs(i, *, *), 3), wnaero, wn)
    ENDELSE
  ENDIF ELSE BEGIN
    IF ngas EQ 1 THEN BEGIN
      IF naerosols EQ 1 THEN BEGIN
        tau_abs = tau_gas_abs + INTERPOL(tau_aero_abs, wnaero, wn)
      ENDIF ELSE BEGIN
        tau_abs = tau_gas_abs
        FOR i = 0, nlayer - 1L DO tau_abs(i, *) += INTERPOL(TOTAL(tau_aero_abs(i, *, *), 3), wnaero, wn)
      ENDELSE
    ENDIF ELSE BEGIN
      IF naerosols EQ 1 THEN BEGIN
        tau_abs = TOTAL(tau_gas_abs, 3) + INTERPOL(tau_aero_abs, wnaero, wn)
        ;tau_abs = INTERPOL(tau_aero_abs, wnaero, wn) ;;;;;;;;;;;;;;;;2020/07/05 to calculate absorption by earosols. kogure
      ENDIF ELSE BEGIN
        tau_abs = TOTAL(tau_gas_abs, 3)
        FOR i = 0, nlayer - 1L DO tau_abs(i, *) += INTERPOL(TOTAL(tau_aero_abs(i, *, *), 3), wnaero, wn); A=A+100  eq  A+=100
      ENDELSE
    ENDELSE
  ENDELSE
  ;-------------- kogure Risei 2021/1/7 -----------------------------------
  openw,lun,path_output + savename +'tau_abs.txt',/get_lun
  printf,lun,tau_abs
  close,lun
  free_lun,lun
  ;-------------- end -----------------------------------------------
  tau_abs = TRANSPOSE(tau_abs)
  IF noabs THEN tau_abs(*) = 0D

  IF KEYWORD_SET(deltatauabs) THEN tau_abs(*, deltatauabs(0)) += tau_abs(*, deltatauabs(0)) * deltatauabs(1)
  IF KEYWORD_SET(deltatausca) THEN tau_aero_sca(deltatausca(0), *, deltatausca(1)) += tau_aero_sca(deltatausca(0), *, deltatausca(1)) * deltatausca(2)

  IF get_tau_sca EQ 1L THEN RETURN, tau_aero_sca
  IF get_tau_abs EQ 1L THEN RETURN, TRANSPOSE(tau_abs)

  IF plottotalopticaldepth THEN BEGIN
    p = PLOT(wnaero, TOTAL(tau_abs, 2), XTITLE = "Wavenumber [cm-1]", YTITLE = "Total absorption optical depth [1]", THICK = 2, /BUFFER, /YLOG)
    p.SAVE, path_output  +  '/Figures/' + savename + '_TotalOpticalDepthAbsorption.jpg', RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
    p.CLOSE
    p = PLOT(wnaero, TOTAL(TOTAL(tau_aero_sca, 3), 1), XTITLE = "Wavenumber [cm-1]", YTITLE = "Total scattering optical depth [1]", THICK = 2, /BUFFER, /YLOG)
    p.SAVE, path_output  +  '/Figures/' + savename + '_TotalOpticalDepthsScattering.jpg', RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
    p.CLOSE
  ENDIF

  ; -----------------------------------------------------------
  ; Write the optical properties of the atmosphere for JACOSPAR
  ; -----------------------------------------------------------
  IF N_ELEMENTS(wn_albedo) EQ 1 THEN BEGIN
    albedoloc = DBLARR(N_ELEMENTS(wn))
    albedoloc(*) = S_albedo
  ENDIF ELSE albedoloc = INTERPOL(S_albedo, wn_albedo, wn)
  dummy = FLTARR(naerosols * 2)
  OPENW, lun, path_output + savename + '.opt', /GET_LUN
  FOR j = 0, N_wn_aero - 2 DO WRITEU, lun, FLOAT(albedoloc(j * nmembers : (j + 1) * nmembers - 1)), dummy
  FOR i = 0, nlayer - 1 DO BEGIN
    FOR j = 0, N_wn_aero - 2 DO BEGIN
      ;WRITEU, lun, FLOAT(tau_abs(j * nmembers : (j + 1) * nmembers - 1, i));2021/2/21 kogure comment out
      ;      FOR k = 0, naerosols - 1 DO WRITEU, lun, FLOAT(tau_aero_sca(i, j, k)), FLOAT(g(j, k, i));2021/2/21 kogure comment out
      writeu,lun,FLOAT(tau_abs(j * nmembers : (j + 1) * nmembers - 1, i)), float(tau_aero_sca(i, j, 0)), float((j)*2.0+1.0), float(tau_aero_sca(i, j, 1)), float((j)*2.0+2.0);;;0.1cm-1で計算するためにj/10をした。Use phase function, not g-parameter. Risei
      ;writeu,lun,FLOAT(tau_abs(j * nmembers : (j + 1) * nmembers - 1, i)), float(tau_aero_sca(i, j, 0)), float(j)*2.0+1.0, tau_aero_sca(i, j, 1), float(j)*2.0+2.0
    ENDFOR
  ENDFOR
  FREE_LUN, lun

  ;;;;;;;;;;;;;;;;;;;;;; write down opt file , kogure 2020/7/1 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF N_ELEMENTS(wn_albedo) EQ 1 THEN BEGIN
    albedoloc = DBLARR(N_ELEMENTS(wn))
    albedoloc(*) = S_albedo
  ENDIF ELSE albedoloc = INTERPOL(S_albedo, wn_albedo, wn)
  dummy = FLTARR(naerosols * 2)
  
  OPENW, lun, path_output + savename + 'tau.txt', /GET_LUN
  pretausca = total(tau_aero_sca, 3)
  printf,lun,'alt(i), tau(alt_i)', wnaero(-1),'cm-1'
  for i = 0, nlayer - 1 do begin
    printf,lun, ztpn(0, i), tau_abs(-1, i) + pretausca(i, -1)
  endfor
  close,lun
  free_lun,lun
  
  OPENW, lun, path_output + savename + 'optcheck.txt', /GET_LUN
  ;FOR j = 0, N_wn_aero - 2 DO printf, lun, FLOAT(albedoloc(j * nmembers : (j + 1) * nmembers - 1));, dummy
  ;printf,lun,'albedo'
  ;printf,lun,albedoloc
  printf,lun,'wavenumber'
  printf,lun,wn
  printf,lun,'tau_absorption'
  printf,lun,tau_abs
  printf,lun,'tau_aero_sca'
  printf,lun,transpose(tau_aero_sca)
  ;  FOR i = 0, nlayer - 1 DO BEGIN
  ;    FOR j = 0, N_wn_aero - 2 DO BEGIN
  ;      printf, lun, FLOAT(tau_abs(j * nmembers : (j + 1) * nmembers - 1, i))
  ;      FOR k = 0, naerosols - 1 DO printf, lun, FLOAT(tau_aero_sca(i, j, k)), FLOAT(g(j, k, i))
  ;    ENDFOR
  ;  ENDFOR
  FREE_LUN, lun
  ;;;;;;;;;;;;;;;   end , kogure 2020/7/1   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; --------------
  ; Phase function;  2021/2/12 kogure comment out
  ; --------------
  ;  OPENW, lun, path_output + savename + '.phs', /GET_LUN
  ;  WRITEU, lun, FLOAT(angle)
  ;  FOR i = 0, N_wn_aero - 2 DO FOR j = 0, naerosols - 1 DO WRITEU, lun, FLOAT(Phasefunc(i, *, j, 0))
  ;  FREE_LUN, lun

  ; -----------------------------------------
  ; Write the configuration file for JACOSPAR
  ; -----------------------------------------
  FileNames = UPDATE_JACOSPAR_INPUT_PARALLEL_Limb_ShoheiDatabase(sequence, savename, FOV,naerosols * (N_wn_aero - 1), n_angles, precision, seed, ngroups, nmembers, naerosols, ngas, path_output, alt_range, z, GEOMETRY = geometry, _EXTRA = extra) ; dz

  RETURN, FileNames
END