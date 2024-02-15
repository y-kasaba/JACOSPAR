FUNCTION READ_JACOSPAR_OUTPUT_FILE, orbit, $
  wnrange, $
  gas_species, $
  aerosols_species, $
  alt_range, $
  dwn, $
  dwnaero, $
  z, $
  geometry, $
  
  wn_gas, $
  dRdalb, $
  dRdabs, $
  dRdsca, $
  R_std, $
  dRdalb_std, $
  dRdabs_std, $
  dRdsca_std, $

  SAVENAME = savename, $
  SOLAR = solar, $

  PATH_OUTPUT = path_output, $

  _EXTRA = extra
  
  wn_gas = wnrange(0) + DINDGEN(ROUND((wnrange(1) - wnrange(0)) / dwn) + 1L) * dwn
  wn_aero = wnrange(0) + dwnaero / 2D + DINDGEN(ROUND((wnrange(1) - wnrange(0)) / dwnaero) + 1L) * dwn
  N_wn_gas = N_ELEMENTS(wn_gas)
  N_wn_aero = N_ELEMENTS(wn_aero)
  ngas = N_ELEMENTS(gas_species)
  naerosols = N_ELEMENTS(aerosols_species)
  nspectra = N_ELEMENTS(geometry(0, *))
  nlayer = N_ELEMENTS(z)
  
  ; -------------------
  ; Path and file names
  ; -------------------
  pathnames = PATH_MANAGEMENT()
  CD, CURRENT = home
  IF ~KEYWORD_SET(path_output) THEN path_output = pathnames.results ; where the results will be written
  IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
  path_solar = pathnames.solar ; path of the solar radiance

  IF ~KEYWORD_SET(savename) THEN savename = "Orbit_" + orbit
  
  ; -----------------------
  ; Read the output spectra
  ; -----------------------
  ngroups = N_wn_aero - 1L
  nmembers = ROUND((N_wn_gas - 1D) / (N_wn_aero - 1D))  

  R_std = DBLARR(ngroups * nmembers + 1, nspectra)
  dRdalb = DBLARR(ngroups * nmembers + 1, nspectra)
  dRdalb_std = DBLARR(ngroups * nmembers + 1, nspectra)
  dRdabs = DBLARR(ngroups * nmembers + 1, nlayer, nspectra)
  dRdabs_std = DBLARR(ngroups * nmembers + 1, nlayer, nspectra)
  dRdsca = DBLARR(ngroups * nmembers + 1, naerosols, nlayer, nspectra)
  dRdsca_std = DBLARR(ngroups * nmembers + 1, naerosols, nlayer, nspectra)
  sun_distance =  DBLARR(nspectra) & sun_distance(*) = 1.6446D
  R = OUTPUT_JACOSPAR_PARALLEL(savename, path_output, nlayer, nspectra, ngroups, nmembers, naerosols, ngas, -1, -1, sun_distance, $
    R_std, dRdalb, dRdalb_std, dRdabs, dRdabs_std, dRdsca, dRdsca_std, _EXTRA = extra)
  
  R *= !dpi
  R_std *= !dpi
  dRdalb *= !dpi
  dRdalb_std *= !dpi
  dRdabs *= !dpi
  dRdabs_std *= !dpi
  dRdsca *= !dpi
  dRdsca_std *= !dpi
  
  RETURN, R  
END