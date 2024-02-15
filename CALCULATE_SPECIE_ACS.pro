; Calculating the file if not already existing
CASE molec OF
  "H2O": BEGIN
    file = pathnames.hitran + '01.par'
    vmr = 1D-6
  END
  "CO2": BEGIN
    file = pathnames.hitran + '02.par'
    vmr = 0.96D
  END
  "CO": BEGIN
    file = pathnames.hitran + '05.par'
    vmr = 1D-4
  END
ENDCASE

; Load the HITRAN file if not provided
IF ~KEYWORD_SET(hitran) THEN Hitran = READ_HITRAN2012(FILE = file)

; Set the isotopologue list if all isotopologues are requested
uniqueIso = Hitran[UNIQ(Hitran.isot, SORT(Hitran.isot))].isot
IF ~iso THEN isotodo = INDGEN(MAX(uniqueIso), START = 1) ELSE isotodo = iso

; Set the indices calculation if not specified
IF ~KEYWORD_SET(istart) THEN istart = 0L
IF ~KEYWORD_SET(iend) THEN iend = 10000L

; Prints info if requested
IF nbands THEN PRINT, molec + ":"

; To check the intensity with respect to HITRAN
kb = 1.38064852D-16
IF checkintensity THEN sumintensity = DBLARR(N_ELEMENTS(tmpval), 2)

; Loop on the isotopologues
FOR i = 0, N_ELEMENTS(isotodo) - 1 DO BEGIN
  ; Set the counter for missing bands if requested
  IF countmissing GE 0L THEN BEGIN
    countmissing = 0L
    somme = 0L
  ENDIF
  IF doprint THEN PRINT, "Preparing calculation of " + molec + " iso " + STRING(isotodo(i), "(I1)")
  ; Consider only the lines of the requested isotopologue
  HitranIso = Hitran[WHERE(Hitran.isot EQ isotodo(i))]
  ; List the unique vibrational bands
  inds = UNIQ(HitranIso.uppq + HitranIso.lowq, SORT(HitranIso.uppq + HitranIso.lowq))
  uniqueBand = [[HitranIso[inds].uppq], [HitranIso[inds].lowq]]
  IF nbands THEN BEGIN
    PRINT, "Iso " + STRING(isotodo(i), "(I1)") + ": " + STRING(N_ELEMENTS(uniqueBand(*, 0)), "(I4)")
    CONTINUE
  ENDIF
  ; Loop on the vibrational bands
  FOR j = MAX([0, istart]), MIN([N_ELEMENTS(uniqueBand(*, 0)) - 1, iend]) DO BEGIN
    ; Consider only the required ro-vibrational bands
    HitranIsoBand = HitranIso[WHERE(HitranIso.uppq EQ uniqueBand(j, 0) AND HitranIso.lowq EQ uniqueBand(j, 1))]
    ; List the wavenumbers of the band which are in the requested wavenumber region
    inds = WHERE(HitranIsoBand.wn GE wnmin - FBord AND HitranIsoBand.wn LE wnmax + FBord)
    IF inds(0) NE -1 THEN BEGIN
      IF doprint THEN BEGIN
        PRINT, STRING(j, "(I04)") + "/" + STRING(N_ELEMENTS(uniqueBand(*, 0)), "(I04)") + ": Computing " + molec + " iso " + STRING(isotodo(i), "(I1)") + " band " + STRTRIM(uniqueBand(j, 0), 1) + "-" + STRTRIM(uniqueBand(j, 1), 2)
        TIC
      ENDIF
      IF TOTAL(HitranIsoBand.INTENS) LT Imin THEN BEGIN
        IF doprint THEN PRINT, "Intensity too small"
        CONTINUE
      ENDIF
      ; Compute only if at least one line of the band is the wavenumber region
      IF countmissing GE 0 THEN BEGIN
        value = FILE_TEST(PathName + molec + "_" + STRING(isotodo(i), "(I1)") + "_" + STRTRIM(uniqueBand(j, 0), 2) + "_" + STRTRIM(uniqueBand(j, 1), 2) + ".kasi")
        countmissing = countmissing + value
        somme = somme + 1L
        IF ~value THEN PRINT, STRING(j, "(I04)") + ": " + STRTRIM(uniqueBand(j, 0), 2) + "_" + STRTRIM(uniqueBand(j, 1), 2)
        CONTINUE
      ENDIF
      ; Get the absorption coefficient of the band
      kloc = CREATE_ACS(molec, isotodo(i), STRTRIM(uniqueBand(j, 0), 2), STRTRIM(uniqueBand(j, 1), 2), vmr, wnloc, $
        HITRAN = HitranIsoBand, _EXTRA = extra)
      ; Skip if no absorption coefficient has been calculated (i.e. there are no lines in the wavenumber region)
      IF N_ELEMENTS(kloc) EQ 1L THEN CONTINUE
      IF correctisoratio THEN BEGIN
        kloc = kloc / MOLECULES_DATA(molec, iso, KEY = 'abundance')
        FileName = molec + "_" + STRING(iso, "(I1)") + "_" + STRTRIM(uniqueBand(j, 0), 2) + "_" + STRTRIM(uniqueBand(j, 1), 2)
        OPENW, lun, PathName + FileName + ".kasi", /get_lun
        WRITEU, lun, FLOAT(kloc)
        FREE_LUN, lun
      ENDIF
      IF checkintensity THEN BEGIN
        FOR k = 0, N_ELEMENTS(tmpval) - 1 DO BEGIN
          HitLine = CORRECT_LP(HitranIsoBand, tmpval(k), tmpval(k), prsval(k), prsval(k) * vmr, 1D, 0D)
          sumintensity(k, 0) = sumintensity(k, 0) + TOTAL(HitLine[WHERE(HitLine.wn GT wnmin AND HitLine.wn LT wnmax)].Sic)
        ENDFOR
      ENDIF
      ; Interpol on the pressure and temperature grid and add to the total absorption coefficient
      indswnloc = WHERE(wnloc GE wn(0) AND wnloc LE wn(nwn - 1))
      wnloc = wnloc(indswnloc)
      kloc = kloc(indswnloc, *, *)
      indswn = WHERE(wn GE wnloc(0) AND wn LE wnloc(N_ELEMENTS(wnloc) - 1))
      IF KEYWORD_SET(tmpval) THEN FOR k = 0, N_ELEMENTS(tmpval) - 1 DO kval(indswn, k) = kval(indswn, k) + $
        INTERPOL_K(kloc, prsval(k), tmpval(k), prs, tmp) $
      ELSE kval(indswn, *, *) = kval(indswn, *, *) + kloc
      IF doprint THEN TOC
    ENDIF
  ENDFOR
  IF countmissing GE 0 THEN PRINT, "For " + molec + " iso " + STRING(isotodo(i), "(I1)") + ": " + $
    STRING(countmissing, "(I04)") + "/" + STRING(somme, "(I04)")
ENDFOR

; Checkintensity: integrate the spectrum
IF checkintensity THEN BEGIN
  PRINT, "Intensity check: "
  FOR k = 0, N_ELEMENTS(tmpval) - 1 DO BEGIN
    sumintensity(k, 1) = INTEGRATE(wn, kval(*, k)) / (prsval(k) / kb / tmpval(k) * 1D3)
    diff = (sumintensity(k, 0) - sumintensity(k, 1)) / sumintensity(k, 0) * 100
    PRINT, "T = " + STRING(tmpval(k), "(I3)") + "K, p = " + STRING(prsval(k), "(E11.4)") + "mbar, HITRAN: " + $
      STRING(sumintensity(k, 0), "(E11.4)") + ", Integral: " + STRING(sumintensity(k, 1), "(E11.4)") + ", Percentage: " + $
      STRING(ABS(diff), "(E11.4)") + "%"
  ENDFOR
ENDIF

; Save if required
IF dosave THEN BEGIN
  OPENW, lun, PathName + FileName + "_wn.kasi", /get_lun
  WRITEU, lun, FLOAT(wn)
  FREE_LUN, lun
  
  OPENW, lun, PathName + FileName + ".kasi", /get_lun
  WRITEU, lun, FLOAT(kval)
  FREE_LUN, lun
ENDIF
IF countmissing GE 0 THEN RETURN, 0
ENDELSE

; Plot if required
IF doplot AND KEYWORD_SET(tmpval) AND KEYWORD_SET(prsval) THEN BEGIN
  p = OBJARR(N_ELEMENTS(tmpval))
  col = ['b', 'r', 'y', 'g', 'k', 'c']
  FOR i = 0, N_ELEMENTS(tmpval) - 1 DO $
    p(i) = PLOT(wn, kval(*, i), /YLOG, COLOR = col(i), XTITLE = "Wavenumber [$cm^{-1}$]", YTITLE = "Abs. coeff. [CGS]", $
    NAME = molec + "(" + STRING(iso, "(I1)") + ") - Tmp = " + STRING(tmpval(i), "(F5.1)") + " K - Prs = " + $
    STRING(prsval(i), "(E10.1)") + " mbar", /OVERPLOT, XRANGE = [wnmin, wnmax], YRANGE = [MAX([1D-28, MIN(kval)]), MAX(kval)])
  L = LEGEND(TARGET = p, /AUTO_TEXT_COLOR, /NORMAL)
ENDIF

; Return the absorption coefficient
RETURN, kval
