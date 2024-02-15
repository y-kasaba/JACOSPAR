FUNCTION CREATE_AEROSOLS_OPTICAL, specie, RMean = RMean, INDICES = indices, FORCE = force, _EXTRA = extra
  
  IF ~KEYWORD_SET(force) THEN force = 0
  pathnames = PATH_MANAGEMENT()
  CASE STRLOWCASE(specie) OF
    'dust': data = READ_ASCII(pathnames.aerosols + 'Dust_rfi_Wolff_2009_in_wn.dat')
    'ice': data = READ_ASCII(pathnames.aerosols + 'waterice_warren.txt')
    ;'so2': data = TRANSPOSE([[11000D, 28600D],[1.34D, 1.34D],[0D, 0D]])
    ;'s2': data = TRANSPOSE([[11000D, 28600D],[1.96D, 1.96D],[0D, 0D]])
  ENDCASE
  IF STRCMP(specie, 'dust') OR STRCMP(specie, 'ice') THEN data = DOUBLE(data.FIELD1)

  IF ~KEYWORD_SET(Rmean) THEN Rmean = [DINDGEN(9) * 0.001D + 0.001D, DINDGEN(20) * 0.01D + 0.01D, DINDGEN(18) * 0.1D + 0.3D, DINDGEN(6) * 0.5D + 2.5D, DINDGEN(5) * 5D + 10D] * 1D-3
  IF KEYWORD_SET(indices) THEN Rmean = Rmean(indices)
  FOR i = 0, N_ELEMENTS(Rmean) - 1 DO BEGIN
    IF NOT(FILE_TEST(pathnames.aerosols + specie + '_' + STRING(Rmean(i), "(E8.2)") + '.aero')) OR force THEN BEGIN
      data = GET_AEROSOLS_OPTICAL(specie, Rmean(i), .3D, OPTICALFILE = data, _EXTRA = extra)
      SAVE, data, FILENAME = pathnames.aerosols + specie + '_' + STRING(Rmean(i), "(E8.2)") + '.aero'
    ENDIF
  ENDFOR
END