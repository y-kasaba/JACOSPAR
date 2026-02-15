FUNCTION READ_GEOMETRY_2025, orbit

  ; Structure of the geometry files
  ; lat, lon, LT, Ls, inc, em, phase, altimetry, slant_dist, dphi, ina, altmex

  pathnames = PATH_MANAGEMENT()
  pathname = pathnames.geometry
  sdir = FILE_SEARCH(pathname + 'Geom_' + STRMID(orbit, 0, 4) + "_kogure20260131.dat")

  IF STRLEN(sdir(0)) EQ 0 THEN BEGIN
    PRINT, 'Files ' + orbit + ' do not exist.'
    RETURN, !NULL
  ENDIF
  
  geom = READ_ASCII(sdir(0), DATA_START = 1)
  geom = geom.FIELD01

  RETURN, geom
END