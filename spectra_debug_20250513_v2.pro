; spectra_debug_20250513.pro

PRO spectra_debug_20250513_v2
  ;--- パス設定 ---
  path_sav_file  = '/work1/data/JACOSPAR/Inversion/OMEGA/CUBE/'
  path_save = '/work1/data/JACOSPAR/Inversion/OMEGA/Results'
  dirsoft='/work1/data/JACOSPAR/Inversion/OMEGA/CUBE/SOFTWARE/SOFT05/'
  path_emiliano = '/work1/data/JACOSPAR/Inversion/OMEGA/SPECTRA/Converted_Spectra_0072.dat'
  path_matlab_spectra = '/work1/data/JACOSPAR/Inversion/OMEGA/SPECTRA/Spectra_0072_2.dat'

  xsize=1000
  ysize=600
  window,xsize=xsize,ysize=ysize
  device, retain=2, decomposed=0,SET_FONT='DejaVuSans', /TT_FONT
  loadct,39
  !p.background = 255
  !p.color = 0

  ;--- SAV ファイル一覧取得 ---
  files = FILE_SEARCH(path_sav_file + 'orb0072_2.sav', COUNT=count)
  FOR loop = 0, count-1 DO BEGIN


    sdir = FILE_SEARCH(path_save + STRUPCASE(FILE_BASENAME(files[loop],'.sav')) + '.jpg')


    PRINT, 'Processing ', loop+1, ' / ', count
    ;LOADCT, 39
    RESTORE, files[loop]

    ;--- Mars スペクトル読み込み ---
    openr,2,dirsoft+'specsol_0403.dat'
    specmars=0B
    specmars=fltarr(352)
    readf,2,specmars
    close,2
    specmars=specmars/dmars/dmars

    wvl=wvl[0:351]
    nwvl=n_elements(wvl)
    ;calculate altitude and set ind
    ;limb=where(alt gt 60)
    alt=alt-65.536
    maxaltitude=40
    if max(alt) le maxaltitude then maxaltitude=max(alt)
    if maxaltitude lt 10 then continue
    koudo=indgen(round(maxaltitude/5))*5+5

    ind=where(alt ge 0 and alt le maxaltitude+0.1)
    nind=n_elements(ind)
    if nind le 1 then continue

    alt=alt(ind)
    lati=lati(ind)
    longi=longi(ind)
    change=where(longi gt 180)
    if n_elements(change) gt 1 then longi(change)=longi(change)-360


    koudonew1=lindgen(n_elements(koudo))
    for i=0,n_elements(koudo)-1 do begin
      j=koudo[i]
      koudo1=where(alt ge j and alt le j+0.1)
      koudonew1(i)=koudo1(0)
    endfor


    altmax=max(alt(*,*))
    altmin=min(alt(*,*))
    longimax=max(longi(*,*))
    longimin=min(longi(*,*))

    ;normalized by specmars
    io=n_elements(jdat(*,1,1))
    ip=n_elements(jdat(1,1,*))
    RadPerSol0=dblarr(io,nwvl,ip)
    for i=0,io-1 do begin
      for o=0,ip-1 do begin
        RadPerSol0(i,*,o)=jdat(i,*,o)/specmars
      endfor
    endfor

    rad=dblarr(nind,nwvl)
    radpersol=dblarr(nind,nwvl)
    for i=0,nwvl-1 do begin
      radpersol1=reform(radpersol0(*,i,*))
      rad1=reform(jdat(*,i,*))
      radpersol(*,i)=double(radpersol1(ind))
      rad(*,i)=double(rad1(ind))
    endfor


    ;--- Emiliano データ読み込み ---
    sdir = FILE_SEARCH(path_emiliano)
    spectra = READ_ASCII(sdir(0), DATA_START = 1)
    spectra = spectra.FIELD001
    wn = spectra(*, 0)                                ; 波数 [cm^-1]
    nspectra = N_ELEMENTS(spectra(0, *)) - 1L
    spectra = spectra(*, 1:nspectra)                  ; 強度行列 [高度, 波数]
    spectra[WHERE(spectra LE 0)] = 1D-4               ; 小暮追加

    ; 波数→波長に変換してソートインデックスを作成
    lambda_emiliano = 1e4 / wn                        ; [μm] = 10^4(cm/μm) ÷ (cm^-1)
    order = SORT(lambda_emiliano)                    ; 昇順ソート用インデックス
    lambda_sorted = lambda_emiliano[order]            ; ソート済み波長ベクトル

    ;--- Emiliano 高度配列の初期化 ---
    altitudes_emiliano = FLTARR(8)
    altitudes_emiliano = [5.52,11.48,14.45,20.41,26.37,29.35,35.31,41.26]



    ; ---- Matlabのスペクトル読み込み ----
    sdir = FILE_SEARCH(path_matlab_spectra)  
    mat_spectra = READ_ASCII(sdir, DATA_START = 1)
    ;  spectra = spectra.FIELD001
    ;  spectra = spectra.FIELD1
    mat_spectra = mat_spectra.FIELD01
    mat_wn = mat_spectra(*, 0)
    mat_wavelen = 1/mat_wn * 1e4
    mat_nspectra = N_ELEMENTS(mat_spectra(0, *)) - 1L
    mat_spectra = mat_spectra(*, 1 : mat_nspectra)

    ;--- 各高度でプロット ---
    ;--- 各高度でプロット（点のみ） ---
    FOR k = 0, N_ELEMENTS(koudonew1)-1 DO BEGIN
      altitude_ours = alt[koudonew1[k]]
      rad_ours      = radpersol[koudonew1[k],*]

      ; 最も近い Emiliano 高度を探す
      diffs = ABS(altitudes_emiliano - altitude_ours)
      idx = WHERE(diffs EQ MIN(diffs))

      altitude_emiliano = altitudes_emiliano[idx[0]]
      observed_emiliano = spectra[*,idx[0]]
      observed_sorted = observed_emiliano[order]
      ; 波長変換・ソート済みなら observed_sorted, lambda_sorted を使っても可

      color = 255.0 * (k+1) / N_ELEMENTS(koudonew1)

      ; 我々のスペクトル：データ点だけ
      IF k EQ 0 THEN BEGIN
        PLOT, wvl, rad_ours, PSYM=4, SYMSIZE=0.7, COLOR=color, YRANGE=[0,0.3], $
             TITLE='Spectra Comparison', $
            XTITLE='Wavelength (μm)', YTITLE='Intensity', /XSTY, /YSTY
      ENDIF ELSE BEGIN
        OPLOT, wvl, rad_ours, PSYM=4, SYMSIZE=0.7, COLOR=color
      ENDELSE

      ; Emiliano のスペクトル：データ点だけ
      OPLOT, lambda_sorted, observed_sorted, LINESTYLE=0, THICK=2, COLOR=color

      ; Matlab のスペクトル:点線
      observed_matlab = mat_spectra[*, idx[0]]
      OPLOT, mat_wavelen, observed_matlab, LINESTYLE=1, THICK=2, COLOR=color

      ; 高度情報をグラフに表示
      text_position_x = MAX(wvl) * 0.9 ; X座標（波長の右端付近）
      text_position_y = 0.28 - 0.01 * k ; Y座標（上から順に表示）
      XYOUTS, text_position_x, text_position_y, $
         'Ours: ' + STRING(altitude_ours, FORMAT='(F6.1)') + ' km, ' + $
         'Emiliano: ' + STRING(altitude_emiliano, FORMAT='(F6.2)') + ' km', $
         COLOR=color, ALIGNMENT=1


      PRINT, 'Plotted ours=', STRING(altitude_ours,'(F6.1)'), 'km, Emiliano=', STRING(altitude_emiliano,'(F6.2)'), 'km'
    ENDFOR


  ENDFOR
END
