pro JACOSPAR_radiance_v2
  ;wn: wavenumber [cm^-1] of the calculation
  ;sp: output value from JACOSPAR

 
  dataName = 'OMEGA/Results/003720251010_dwn10_0-100km_limbCase2'
  data = read_ascii(dataName + '/rad')
  data = data.field1
  sp = data(1, *)         ; radiance
  wl = data(0, *)         ; wavelength [μm]
  wn = 1.0 / wl * 1e4     ; convert wavelength [μm] to wavenumber [cm^-1]

  ;; unit correction (optional)
  ;for i = 0L, 100 do sp(i) = sp(i) * (1d4 / (wn(i))) / wn(i) / 1d4 * 1d7

  ;--- Solar spectrum correction ---
  wn_pfsolspec_hr = fltarr(875482)
  rad_pfsolspec_hr = fltarr(875482)
  file_pfsolspec_hr = 'Data/Sun/pfsolspec_hr.dat'

  openr, 1, file_pfsolspec_hr
  for i = 0L, 875481 do begin
    readf, 1, a, b
    wn_pfsolspec_hr(i) = a
    rad_pfsolspec_hr(i) = b * !dpi / (1.5D^2D) ; convert to Mars [1.5 AU]
  endfor
  close, 1

  solar_flux = interpol(rad_pfsolspec_hr, wn_pfsolspec_hr, wn)
  sp = sp * solar_flux
  sp = sp / !pi

  ;--------------------------------
  ; Plot section
  ;--------------------------------
  device, decomposed=0  ; 旧式カラーインデックスモードに戻す
!p.multi = 0
!x.range = 0
!y.range = 0
  !p.background = 255   ; 白
  !p.color = 0           ; 黒
  erase, !p.background   ; 背景を白でクリア
  
  window, /free, title='Radiance Spectrum', xsize=800, ysize=500
  device, decomposed=0
  !p.background = 255
  !p.color = 0
  erase, !p.background

  plot, wn, sp, $
  xrange=[min(wn), max(wn)], $
  yrange=[0.0, 0.4], $
  xtitle='Wavenumber [cm!u-1!n]', $
  ytitle='Radiance (normalized)', $
  thick=2, $
  color=2     ; 青（インデックスカラー2）
  aspect_ratio=0.5


$stop
end
