FUNCTION READ_GEOMETRY_Emiliano, orbit

  ; Structure of the geometry files
  ; lat, lon, LT, Ls, inc, em, phase, altimetry, slant_dist, dphi, ina, altmex

  pathnames = PATH_MANAGEMENT()
  pathname = pathnames.geometry
  sdir = FILE_SEARCH(pathname + 'Geom_EmilianoAlt_' + STRMID(orbit, 0, 4) + "_v5.dat")

  IF STRLEN(sdir(0)) EQ 0 THEN BEGIN
    PRINT, 'Files ' + orbit + ' do not exist.'
    RETURN, !NULL
  ENDIF
  
  geom = READ_ASCII(sdir(0), DATA_START = 1)
  geom = geom.FIELD1

  RETURN, geom
END