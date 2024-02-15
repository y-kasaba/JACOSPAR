FUNCTION CREATE_COMBINED_ACS, gas_species, wnstep, UPDATE = update, INDMIN = indmin, INDMAX = indmax, WN = wn

  IF ~KEYWORD_SET(update) THEN update = 0
  IF ~KEYWORD_SET(indmin) THEN indmin = 0
  IF ~KEYWORD_SET(indmax) THEN indmax = 100000L
  
  wnmin = 3000D
  wnmax = 10000D
  nwn = LONG((wnmax - wnmin) / wnstep)
  wn = wnmin + DINDGEN(nwn) * wnstep

  path = PATH_MANAGEMENT()
  pathacs = path.acs
  FileSave = pathacs + STRUPCASE(gas_species) + "_" + STRING(wnstep, "(E6.0)") + "_0"
  
  listfiles_wn = FILE_SEARCH(pathacs + STRUPCASE(gas_species) + '_*_wn.kasi')  
  listfiles_p = FILE_SEARCH(pathacs + STRUPCASE(gas_species) + '_*_p.kasi')
  listfiles_T = FILE_SEARCH(pathacs + STRUPCASE(gas_species) + '_*_T.kasi')
  
  prsgrid = READ_BINARY(listfiles_p[0], DATA_TYPE = 5)
  tmpgrid = READ_BINARY(listfiles_T[0], DATA_TYPE = 5)
  nprs = N_ELEMENTS(prsgrid)
  ntmp = N_ELEMENTS(tmpgrid)
  
  IF update THEN BEGIN
    acs = READ_BINARY(FileSave + ".kasi", DATA_TYPE = 5, DATA_DIMS = [nwn, ntmp, nprs])
  ENDIF ELSE acs = DBLARR(nwn, nprs, ntmp)
  
  FOR i = indmin, MIN([indmax, N_ELEMENTS(listfiles_wn) - 1]) DO BEGIN
    IF STRPOS(listfiles_wn(i), 'E+') GT 0  OR STRPOS(listfiles_wn(i), 'E-') GT 0 THEN CONTINUE
    listfiles_acs = STRJOIN(STRSPLIT(listfiles_wn(i), "_wn", /EXTRACT, /REGEX),'')
    PRINT, "Adding " + STRING(i, "(I04)") + "/" + STRING(N_ELEMENTS(listfiles_wn) - 1, "(I04)")
    wnloc = READ_BINARY(listfiles_wn(i), DATA_TYPE = 5)
    IF N_ELEMENTS(wnloc) NE 3 THEN BEGIN
      PRINT, "Problem with " + listfiles_acs + ", please correct the wavenumber grid using MATLAB"
      CONTINUE
    ENDIF
    IF wnloc(2) GT wnstep THEN BEGIN
      PRINT, "Error with " + listfiles_acs + ", the wavenumber step is too big!"
      CONTINUE
    ENDIF
    wnlocsave = wnloc
    wnlocall = wnloc(0) + DINDGEN(ROUND((wnloc(1) - wnloc(0)) / wnloc(2)) + 1) * wnloc(2)
    nwnloc = N_ELEMENTS(wnlocall)
    tokeep= WHERE(wnlocall GE wn(0) AND wnlocall LE wn(nwn - 1), /NULL)
    IF N_ELEMENTS(tokeep) GT 0 THEN BEGIN
      wnloc = DBLARR(nwnloc)
      wnloc = wnlocall(tokeep)
      dumb = MIN(ABS(wnloc(0) - wn), indsmin)
      dumb = MIN(ABS(wnloc(N_ELEMENTS(wnloc) - 1) - wn), indsmax)
      info = FILE_INFO(listfiles_acs)
      IF info.size / 8D - (nwnloc * ntmp * nprs) GT 1000D THEN BEGIN ;N_ELEMENTS(acsloc) NE N_ELEMENTS(wnloc) * ntmp * nprs THEN BEGIN
        PRINT, "Error with " + listfiles_acs + ", the length of the pressure and temperature is incorrect!"
        DELVAR, acsloc
        CONTINUE
      ENDIF
;      IF info.size GT 4D9 THEN BEGIN
;        pos = 0LL
;        cc = COLORTABLE(39, NCOLORS = nprs + 1, /TRANSPOSE)
;        FOR k = 0, nprs - 1 DO BEGIN
;          FOR j = 0, ntmp - 1 DO BEGIN
;            acsloc = READ_BINARY(listfiles_acs, DATA_TYPE = 5, DATA_DIMS = nwnloc, DATA_START = pos * 8LL)
;            acsloc = acsloc(tokeep)
;            acs(indsmin : indsmax, k, j) += RESAMPLE(wnloc, acsloc, wn, wnstep)
;            pos += nwnloc
;          ENDFOR
;        ENDFOR
;      ENDIF ELSE BEGIN
        acsloc = READ_BINARY(listfiles_acs, DATA_TYPE = 5, DATA_DIMS = [nwnloc, ntmp, nprs])
        acsloc = acsloc(tokeep, *, *)
        IF wnstep NE wnlocsave(2) THEN acs(indsmin : indsmax, *, *) += RESAMPLE(wnloc, acsloc, wn, wnstep) ELSE acs(indsmin : indsmax, *, *) += acsloc
      DELVAR, acsloc
;      ENDELSE
    ENDIF
  ENDFOR
   
  OPENW, lun, FileSave + "_wn.kasi", /get_lun
  WRITEU, lun, DOUBLE([wnmin, wnmax, wnstep])
  FREE_LUN, lun
  
  OPENW, lun, FileSave + "_p.kasi", /get_lun
  WRITEU, lun, DOUBLE(prsgrid)
  FREE_LUN, lun
  
  OPENW, lun, FileSave + "_T.kasi", /get_lun
  WRITEU, lun, DOUBLE(tmpgrid)
  FREE_LUN, lun
  
  OPENW, lun, FileSave + ".kasi", /get_lun
  WRITEU, lun, DOUBLE(acs)
  FREE_LUN, lun
  
  RETURN, acs
END