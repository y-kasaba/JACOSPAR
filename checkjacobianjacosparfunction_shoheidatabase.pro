PRO CHECKJACOBIANJACOSPARFUNCTION_shoheidatabase, wn_apriori, nlayer, dwn_gas, orbit, alt_range, z, ztpn, geometry, dwn, FOV, Radiusm, Radiuss, wn_albedo, S_albedo, solar, ILS, $
  radtemp, radtempnoconv, dRdabs, dRdsca, tau_abs, tau_aero_sca, PathSave, extra

  delta = .1D
  nspectra = N_ELEMENTS(geometry(0, *))

  radtemp_deltaabs = DBLARR(N_ELEMENTS(wn_apriori), nlayer - 1, nspectra)
  radtemp_deltasca = DBLARR(N_ELEMENTS(wn_apriori), nlayer - 1, nspectra)
  radtempnoconv_deltaabs = DBLARR(N_ELEMENTS(wn_apriori), nlayer - 1, nspectra)
  radtempnoconv_deltasca = DBLARR(N_ELEMENTS(wn_apriori), nlayer - 1, nspectra)
  FOR i = 0, nlayer - 2 DO BEGIN
    radtemp_deltaabs(*, i, *) = CREATE_SYNTHETIC_SPECTRA_limb_shoheidatabase(orbit, alt_range, wn_apriori, z, ztpn, geometry, dwn, FOV, Radiusm(*, *, 0), Radiuss, wn_albedo, S_albedo(*, 0), $
      SOLAR = solar, ILS = ILS, NAMEPLUSFILE = 'deltaabs' + STRING(i, "(I1)") + '_', $
      RADNOCONV = radnoconv, DELTATAUABS = [i, delta], TAU_ABS = tau_abs_loc, TAU_AERO_SCA = tau_aero_sca_loc, _EXTRA = extra)
    radtempnoconv_deltaabs(*, i, *) = radnoconv
    DELVAR, radnoconv
    radtemp_deltasca(*, i, *) = CREATE_SYNTHETIC_SPECTRA_limb_shoheidatabase(orbit, alt_range, wn_apriori, z, ztpn, geometry, dwn, FOV, Radiusm(*, *, 0), Radiuss, wn_albedo, S_albedo(*, 0), $
      SOLAR = solar, ILS = ILS, NAMEPLUSFILE = 'deltasca' + STRING(i, "(I1)") + '_', $
      RADNOCONV = radnoconv, DELTATAUSCA = [i, 0, delta], TAU_ABS = tau_abs_loc, TAU_AERO_SCA = tau_aero_sca_loc, _EXTRA = extra)
    radtempnoconv_deltasca(*, i, *) = radnoconv
    DELVAR, radnoconv
  ENDFOR

  FOR i = 0, nlayer - 2 DO BEGIN
    FOR j = 0, nspectra - 1 DO BEGIN
      p = PLOT(wn_apriori, radtempnoconv_deltaabs(*, i, j) - radtempnoconv(*, j), /BUFFER, COLOR = 'b', LAYOUT = [1, 2, 1], XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Numerical Jacobian", TITLE = "Jacobian absorption, layer " + STRING(i, "(1I)") + " relative to spectrum " + STRING(j, "(1I)"))
      p = PLOT(wn_apriori, dRdabs(*, i, j) * tau_abs(*, i) * delta, /OVERPLOT, /CURRENT, COLOR = 'r', XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Analytical Jacobian")
      p = PLOT(wn_apriori, (radtempnoconv_deltaabs(*, i, j) - radtempnoconv(*, j) - dRdabs(*, i, j) * tau_abs(*, i) * delta), COLOR = 'b', LAYOUT = [1, 2, 2], /CURRENT, XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Difference")
      p.SAVE, STRCOMPRESS(PathSave + 'Abs_noconv_' + STRING(i, "(1I)") + "_" + STRING(j, "(1I)") + '.jpg', /REMOVE_ALL) , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
      p.CLOSE

      p = PLOT(wn_apriori, radtempnoconv_deltasca(*, i, j) - radtempnoconv(*, j), /BUFFER, COLOR = 'b', LAYOUT = [1, 2, 1], XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Numerical Jacobian", TITLE = "Jacobian scattering, layer " + STRING(i, "(1I)") + " relative to spectrum " + STRING(j, "(1I)"))
      p = PLOT(wn_apriori, dRdsca(*, 0, i, j) * tau_aero_sca(i, *, 0) * delta, /OVERPLOT, /CURRENT, COLOR = 'r', XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Analytical Jacobian")
      p = PLOT(wn_apriori, radtempnoconv_deltasca(*, i, j) - radtempnoconv(*, j) - dRdsca(*, 0, i, j) * tau_aero_sca(i, *, 0) * delta, COLOR = 'b', LAYOUT = [1, 2, 2], /CURRENT, XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Difference")
      p.SAVE, STRCOMPRESS(PathSave + 'Sca_noconv' + STRING(i, "(1I)") + "_" + STRING(j, "(1I)") + '.jpg', /REMOVE_ALL) , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
      p.CLOSE
    ENDFOR
  ENDFOR

  FOR i = 0, nlayer - 2 DO BEGIN
    FOR j = 0, nspectra - 1 DO BEGIN
      p = PLOT(wn_apriori, radtemp_deltaabs(*, i, j) - radtemp(*, j), /BUFFER, COLOR = 'b', LAYOUT = [1, 2, 1], XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Numerical Jacobian", TITLE = "Jacobian absorption, layer " + STRING(i, "(1I)") + " relative to spectrum " + STRING(j, "(1I)"))
      p = PLOT(wn_apriori, CONVOL(dRdabs(*, i, j) * tau_abs(*, i) * delta, ILS, /CENTER, /EDGE_TRUNCATE) * dwn_gas, /OVERPLOT, /CURRENT, COLOR = 'r', XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Analytical Jacobian")
      p = PLOT(wn_apriori, radtemp_deltaabs(*, i, j) - radtemp(*, j) - CONVOL(dRdabs(*, i, j) * tau_abs(*, i) * delta, ILS, /CENTER, /EDGE_TRUNCATE) * dwn_gas, COLOR = 'b', LAYOUT = [1, 2, 2], /CURRENT, XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Difference")
      p.SAVE, STRCOMPRESS(PathSave + 'Abs_conv_' + STRING(i, "(1I)") + "_" + STRING(j, "(1I)") + '.jpg', /REMOVE_ALL) , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
      p.CLOSE

      p = PLOT(wn_apriori, radtemp_deltasca(*, i, j) - radtemp(*, j), /BUFFER, COLOR = 'b', LAYOUT = [1, 2, 1], XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Numerical Jacobian", TITLE = "Jacobian scattering, layer " + STRING(i, "(1I)") + " relative to spectrum " + STRING(j, "(1I)"))
      p = PLOT(wn_apriori, CONVOL(dRdsca(*, 0, i, j) * tau_aero_sca(i, *, 0) * delta, ILS, /CENTER, /EDGE_TRUNCATE) * dwn_gas, /OVERPLOT, /CURRENT, COLOR = 'r', XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Analytical Jacobian")
      p = PLOT(wn_apriori, radtemp_deltasca(*, i, j) - radtemp(*, j) - CONVOL(dRdsca(*, 0, i, j) * tau_aero_sca(i, *, 0) * delta, ILS, /CENTER, /EDGE_TRUNCATE) * dwn_gas, COLOR = 'b', LAYOUT = [1, 2, 2], /CURRENT, XTITLE = "Wavenumber [cm^{-1}]", YTITLE = "Difference")
      p.SAVE, STRCOMPRESS(PathSave + 'Sca_conv' + STRING(i, "(1I)") + "_" + STRING(j, "(1I)") + '.jpg', /REMOVE_ALL) , RESOLUTION = 300, /TRANSPARENT, /LANDSCAPE
      p.CLOSE
    ENDFOR
  ENDFOR
END