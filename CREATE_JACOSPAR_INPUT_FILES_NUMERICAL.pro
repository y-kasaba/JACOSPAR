FUNCTION CREATE_JACOSPAR_INPUT_FILES_NUMERICAL, orbit, $
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
                                      savename, $
                                      epsilon, $
                                      fitradius, $
                                      fitsigma, $
                                      
                                      PATH_OUTPUT = path_output, $
                                      NOFITRADIUS = nofitradius, $
                                      ONLYFITSPECIES = onlyfitspecies, $
                                    
                                      _EXTRA = extra

  ngas = N_ELEMENTS(gas_species)
  naerosols = N_ELEMENTS(aerosols_species)
  nspectra = N_ELEMENTS(geometry(0, *))
  nlayer = N_ELEMENTS(ztpn(0, *))

  ; -------------------
  ; Path and file names
  ; -------------------
  pathnames = PATH_MANAGEMENT()
  CD, CURRENT = home
  IF ~KEYWORD_SET(path_output) THEN path_output = pathnames.results ; where the results will be written
  IF ~KEYWORD_SET(fitradius) THEN fitradius = 0
  IF ~KEYWORD_SET(fitsigma) THEN fitsigma = 0
  IF ~KEYWORD_SET(onlyfitspecies) THEN onlyfitspecies = INDGEN(5)
  
  ; -----------------
  ; Create file lists
  ; -----------------
  FileNames = STRARR(nspectra * (1L + nlayer * (TOTAL(onlyfitspecies LT 3) + TOTAL(onlyfitspecies GE 3) * (1L + fitradius + fitsigma))))
  indF = 0
  ; Unperturbed
  FileNames(0 : nspectra - 1L) = CREATE_JACOSPAR_INPUT_FILES(orbit, wnrange, gas_species, aerosols_species, alt_range, dwn, FOV, z, ztpn, Ht, geometry, Radiusm, Radiuss, wn_albedo, S_albedo, $
    PATH_OUTPUT = path_output, SAVENAME = savename, /NOJACOBIAN, _EXTRA = extra)
  indF += nspectra
  ; Perturbed gases and aerosols
  FOR i = 0, ngas - 1L DO BEGIN
    FOR j = 0, nlayer - 1L DO BEGIN
      IF ANY(onlyfitspecies EQ i) AND ztpn(i + 5, j) NE 0D THEN BEGIN
        ztpn_loc = ztpn
        ztpn_loc(i + 5, j) *= (1 + epsilon)
        FileNames(indF : indF + nspectra - 1L) = CREATE_JACOSPAR_INPUT_FILES(orbit, wnrange, gas_species, aerosols_species, alt_range, dwn, FOV, z, ztpn_loc, Ht, geometry, Radiusm, Radiuss, wn_albedo, S_albedo, $
          PATH_OUTPUT = path_output, SAVENAME = savename + '_' + gas_species(i) + '_L' + STRING(j, "(I03)"), /NOJACOBIAN, _EXTRA = extra)
        indF += nspectra
      ENDIF
    ENDFOR
  ENDFOR
  FOR i = 0, naerosols - 1L DO BEGIN
    FOR j = 0, nlayer - 1L DO BEGIN
      IF ANY(onlyfitspecies - 3L EQ i) AND ztpn(i + 5 + ngas, j) NE 0D THEN BEGIN
        ztpn_loc = ztpn
        ztpn_loc(i + 5 + ngas, j) *= (1 + epsilon)
        FileNames(indF : indF + nspectra - 1L) = CREATE_JACOSPAR_INPUT_FILES(orbit, wnrange, gas_species, aerosols_species, alt_range, dwn, FOV, z, ztpn_loc, Ht, geometry, Radiusm, Radiuss, wn_albedo, S_albedo, $
          PATH_OUTPUT = path_output, SAVENAME = savename + '_' + aerosols_species(i) + '_L' + STRING(j, "(I03)"), /NOJACOBIAN, _EXTRA = extra)
        indF += nspectra
      ENDIF
    ENDFOR
  ENDFOR  
  ; Aerosols distribution parameters
  IF fitradius THEN BEGIN
    FOR i = 0, naerosols - 1L DO BEGIN
      FOR j = 0, nlayer - 1L DO BEGIN
        IF ANY(onlyfitspecies - 3L EQ i) AND ztpn(i + 5 + ngas, j) NE 0D THEN BEGIN
          Radiusm_loc = Radiusm
          Radiusm_loc(i, j) += 1D-4
          FileNames(indF : indF + nspectra - 1L) = CREATE_JACOSPAR_INPUT_FILES(orbit, wnrange, gas_species, aerosols_species, alt_range, dwn, FOV, z, ztpn, Ht, geometry, Radiusm_loc, Radiuss, wn_albedo, S_albedo, $
            PATH_OUTPUT = path_output, SAVENAME = savename + '_' + aerosols_species(i) + '_M_L' + STRING(j, "(I03)"), /NOJACOBIAN, _EXTRA = extra)
          indF += nspectra
        ENDIF   
      ENDFOR
    ENDFOR
  ENDIF
  IF fitsigma THEN BEGIN
    FOR i = 0, naerosols - 1L DO BEGIN
      FOR j = 0, nlayer - 1L DO BEGIN
        IF ANY(onlyfitspecies - 3L EQ i) AND ztpn(i + 5 + ngas, j) NE 0D THEN BEGIN
          Radiuss_loc = Radiuss
          Radiuss_loc(i, j) *= (1 + epsilon)
          FileNames(indF : indF + nspectra - 1L) = CREATE_JACOSPAR_INPUT_FILES(orbit, wnrange, gas_species, aerosols_species, alt_range, dwn, FOV, z, ztpn, geometry, Radiusm, Radiuss_loc, wn_albedo, S_albedo, $
            PATH_OUTPUT = path_output, SAVENAME = savename + '_' + aerosols_species(i) + '_S_L' + STRING(j, "(I03)"), /NOJACOBIAN, _EXTRA = extra)
          indF += nspectra
        ENDIF
      ENDFOR
    ENDFOR
  ENDIF

  RETURN, FileNames
END