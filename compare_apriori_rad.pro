; Code to read and compare OMEGA/Results/.../rad

; Get the list of folders
folder_list = FILE_SEARCH('/work1/data/JACOSPAR/Inversion/OMEGA/Results/007220250410_forward_daversa_0-45km_noice_AlbedoToSpectra_daversaSpectra_daversaApriori_EmilianoAlt_*km', /DIRECTORY)

; Loop through each folder and process the single 'rad' file
FOR i = 0, N_ELEMENTS(folder_list) - 1 DO BEGIN
    folder = folder_list[i]
    
    ; Search for the single 'rad' file in the folder
    rad_files = FILE_SEARCH(folder + '/*.rad')
    
    ; Ensure there is exactly one 'rad' file
    IF N_ELEMENTS(rad_files) EQ 1 THEN BEGIN
        rad_file = rad_files[0]
        
        ; Read the data from the rad file
        data = READ_ASCII(rad_file)
        
        ; Extract wavelength and intensity
        wavelength = data[*, 0]
        intensity = data[*, 1]
        
        ;グラフを描画（異なるファイルは異なる色で）
        PLOT, wavelength, intensity, TITLE='Intensity vs Wavelength', XTITLE='Wavelength (nm)', YTITLE='Intensity', COLOR=i + 1
        ;画像を保存
        filename = folder + '/plot_' + STRING(i) + '.png'



        ; Calculate and print statistics for the file
        mean_intensity = MEAN(intensity)
        PRINT, 'Folder: ', folder, ' - File: ', rad_file, ' - Mean Intensity: ', mean_intensity
    ENDIF ELSE BEGIN
        PRINT, 'Warning: Folder ', folder, ' does not contain exactly one .rad file.'
    ENDELSE
ENDFOR
