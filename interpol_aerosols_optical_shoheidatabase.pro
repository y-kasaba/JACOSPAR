FUNCTION INTERPOL_AEROSOLS_OPTICAL_ShoheiDatabase, specie, radiusM_eff, radiusS_eff, $
  FORCEAEROCALC = forceaerocalc, $
  WNVEC = wnvec, $
  WNRANGE = wnrange, $
  WNSTEP = wnstep, $
  DEN = den, $
  SAVECALC = savecalc, $
  _EXTRA = extra
  radiusPath = radiusPath

;;;Kogure Risei 2020/8/10 :  From Arnaud database to Shohei database. 

  IF ~KEYWORD_SET(forceaerocalc) THEN forceaerocalc = 0
  IF ~KEYWORD_SET(getmeanradius) THEN getmeanradius = 0
  IF ~KEYWORD_SET(geteffectiveradius) THEN geteffectiveradius = 0
  IF ~KEYWORD_SET(wnrange) THEN wnrange = [1300D, 25000D];[3000D, 20000D]
  IF ~KEYWORD_SET(wnstep) THEN wnstep = 50D
  IF ~KEYWORD_SET(den) THEN den = 1D
  IF ~KEYWORD_SET(writeascii) THEN writeascii = 0L
  IF ~KEYWORD_SET(noprint) THEN noprint = 0D
  IF ~KEYWORD_SET(savecalc) THEN savecalc = 0L

  ; ----------
  ; Parameters
  ; ----------
  pathnames = PATH_MANAGEMENT()
  PathInput = pathnames.aerosols
  ncol = 7L;kogure 20210810 ncol=8→ ncol=7(Gの列がない)

  radiusM_eff = DOUBLE(STRING(radiusM_eff, "(E9.3)"))
  radiusS_eff = DOUBLE(STRING(radiusS_eff, "(E9.3)"))

  ; -----------------
  ; Set un parameters
  ; -----------------
  IF KEYWORD_SET(wnvec) THEN BEGIN
    nwn = N_ELEMENTS(wnvec)
  ENDIF ELSE BEGIN
    nwn = LONG((wnrange(1) - wnrange(0)) / wnstep) + 1L
    wnvec = wnrange(0) + DINDGEN(nwn) * wnstep
  ENDELSE

  docalc = 1L
  filelist = FILE_SEARCH(pathnames.aerosols, specie + '_*', /MATCH_INITIAL_DOT);Data/aerosols/のリスト
  nfiles = N_ELEMENTS(filelist)
  radiuslist = DBLARR(nfiles)
  variancelist = DBLARR(nfiles)
  offset = STRLEN(pathnames.aerosols) + STRLEN(specie)
  FOR i = 0, nfiles - 1 DO BEGIN
    radiuslist(i) = DOUBLE(STRMID(filelist(i), offset + 1L, 5L)) * 1D-4;リストの中のエアロゾルの粒径を抜き出している
    variancelist(i) = DOUBLE(STRMID(filelist(i), offset + 10L, 5L));リストの中のエアロゾルの分散を抜き出している
    ;variancelist(i) = DOUBLE(STRMID(filelist(i), offset + 7L, 5L));リストの中のエアロゾルの分散を抜き出している
  ENDFOR
  value = MIN(radiuslist - radiusM_eff, ind, /ABSOLUTE)
  ;radiusPath = specie + '_' + STRING(radiuslist(ind) * 1D4, "(F5.3)") + '00_' + STRING(variancelist(ind), "(F5.3)")
  data = READ_ASCII(pathnames.aerosols + specie + '_' + STRING(radiuslist(ind) * 1D4, "(F5.3)") + '00_' + STRING(variancelist(ind), "(F5.3)") + '.dat', COMMENT_SYMBOL = '%');設定した粒径のエアロゾルデータを読み込む    ;kogure 20210810 F5.3→F5.5
  data = TRANSPOSE(data.FIELD1)
 ; data(*,0)=1d/data(*,0)*1d4 ;kogure 20210810 追加
  datafinal = DBLARR(nwn+1, ncol)
  datafinal(0:-2, 0) = wnvec
  FOR i = 1, ncol - 1 DO datafinal(0:-2, i) = INTERPOL(data(*, i), data(*, 0), [1D4/wnvec])
  radiuspath = [radiuslist(ind) * 1D4,variancelist(ind), 0,0,0,0,0 ]
  datafinal(-1,*) = radiuspath
  RETURN, datafinal;1列めは波数、他はエアロゾルファイルの列を波数でinterpolしたもの
END