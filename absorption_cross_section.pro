pro absorption_cross_section
  
  
  aerosol='dust'
;  sla=double(indgen(100))/1d7
;  la=sla+1d-7; 波長
;  k=0.001681 ; 
  r=2d-6;粒径[m]
  pi=!pi
  
  if aerosol eq 'dust' then begin
    filename='absorption_cross_section_dust'
    file=read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Aerosols/Dust_rfi_Wolff_2009_in_wn.dat')
    file=file.field1
    wavenum=file(0,*);cm-1
    lamda=1d/wavenum*1d-2;[um]
  endif else begin
    filename='absorption_cross_section_ice'
    file=read_ascii('/Volumes/HD-PNFU3/JACOSPAR/Inversion/Data/Aerosols/waterice_warren.txt')
    file=file.field1
    lamda=file(0,*)*1d-6;波長[m]
  endelse
  
  k=file(2,*);消衰係数
  
  a=4*pi*k/lamda;吸収係数[m-1] 
  sigma=a*4/3*pi*r^3;吸収断面積(吸収係数[m-1]*分子一個の体積[m3/個]=分子一個の吸収断面積[m2/個])
  sigma=sigma*1d6;吸収係数[cm2/個]
  
  data=dblarr(2,n_elements(lamda))
  labe=['   lamda [m] ','  sigma [cm2] ']
  data(0,*)=lamda
  data(1,*)=sigma
  
  S=!pi*r^2*1d4;面積[cm2]
  
  
  openw,lun,'/Users/juice/Documents/Results_picture/'+filename+strmid(r,3,16)+'.txt',/get_lun
  printf,lun,'radius=',r,' [m]'
  printf,lun,'DiskArea*2',S*2, ' [cm2]'
  printf,lun,labe
  printf,lun,data
  free_lun,lun
  stop
  
  
end