FUNCTION READ_SPECTRA_LimbCase1, orbit, WN = wn, ERROR = error, DOERROR = doerror

  pathnames = PATH_MANAGEMENT()
  PathName = pathnames.spectra
  sdir = FILE_SEARCH(pathname + 'Spectra_' + STRMID(orbit, 0, 4) + "_ORIGINAL.dat")
  ;sdir = FILE_SEARCH(pathname + 'Spectra_' + STRMID(orbit, 0, 4) + ".dat")

  IF STRLEN(sdir(0)) EQ 0 THEN BEGIN
    PRINT, 'Files ' + orbit + ' do not exist.'
    RETURN, !NULL
  ENDIF

  spectra = READ_ASCII(sdir(0), DATA_START = 1)
  spectra = spectra.FIELD001
  ;  spectra = spectra.FIELD1
  ;spectra = spectra.FIELD01
  wn = spectra(*, 0)
  nspectra = N_ELEMENTS(spectra(0, *)) - 1L
  spectra = spectra(*, 1 : nspectra)
  spectra(where(spectra le 0)) = 1d-4; 小暮追加
  IF KEYWORD_SET(doerror) THEN BEGIN
    sdir = FILE_SEARCH(pathname + 'Sigma_Emiliano' + STRMID(orbit, 0, 4) + ".dat")
    error = READ_ASCII(sdir(0), DATA_START = 1)
    ;    error = error.FIELD001
    ;    error = error.FIELD1
    error = error.FIELD01
    error = error(*, 1 : nspectra - 1L)
  ENDIF ELSE BEGIN
    ; Assign 5% error for spectra values >= 0.04, and 10% error for spectra values < 0.04.
    error = spectra * 0.05
    error(where(spectra LT 0.04)) = spectra(where(spectra LT 0.04)) * 0.10
    ;error = spectra / 20D;100D  ;Kogure Risei. D'Aversa+2022より、OMEGA_SWIR_error＝2~5% at I/F>0.04, error=10% at I/F < 0.04.
  ENDELSE; BEGIN
  error(WHERE(error LT 2d-3)) = 2d-3;;;;;add this line. risei kogre 2021/10/8

  RETURN, spectra
END