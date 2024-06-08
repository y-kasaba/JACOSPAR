FUNCTION GETZTPN_LimbTest_new_20240220_invsaTest, lat, lon, lst, ls, z, cov, radius, appriori_files, cov_parameters
  nlayer = N_ELEMENTS(z)
  ztpn = DBLARR(10, nlayer)
  cov = DBLARR(10, nlayer)
  ztpn(0, *) = z
  cov(0, *) = z
  NAv = 6.02D23
  rhoice = .921 ; g/cm3 氷の密度　http://www.godac.jamstec.go.jp/catalog/data/doc_catalog/media/shiken38_02.pdf
  MM_H2O = 18D ; g/mol

  pathnames = PATH_MANAGEMENT()

  filenamePT = appriori_files[0]
  filenameCO2 = appriori_files[1]
  filenameCO = appriori_files[2]
  filenameH2O = appriori_files[3]
  filenamedust = appriori_files[4]
  filenameice = appriori_files[5]
  filename = appriori_files[6]
  
  cov_dust = cov_parameters[0]
  cov_ice = cov_parameters[1]
  cov_dust_r = cov_parameters[2]
  cov_ice_r = cov_parameters[3]

  data = READ_ASCII(filenamePT) & data = data.Field1;to compare with aoki-san
  datadust = READ_ASCII(pathnames.model + filenamedust, COMMENT_SYMBOL = '%') & datadust = datadust.Field1;to compare with aoki-san

  dataice = READ_ASCII(pathnames.model + filenameice, COMMENT_SYMBOL = '%') & dataice = dataice.Field1
  ztpn(1, *) = INTERPOL(data(1, *), data(0, *), z) ; temperature [K];to compare with aoki-san
  ztpn(2, *) = EXP(INTERPOL(ALOG(data(2, *)*1d2), data(0, *), z)) ; pressure [Pa]　          　;to compare with aoki-san

  R=8.314462618d;[J/K/mol];to compare with aoki-san

  ztpn(3, *) = ztpn(2,*)/R/ztpn(1,*)*NAv*1d-6;EXP(INTERPOL(ALOG(data(3, *)), data(0, *), z)) ; total number density [molec/cm3];to compare with aoki-san
  ztpn(5, *) =ztpn(3,*)*0.9532;10^(INTERPOL(dataCO2(1, *), dataCO2(0, *), z)) ; CO2 number density [molec/cm3];to compare with aoki-san
  ztpn(4, *) = 44.01D * ztpn(3, *) * ztpn(5, *) / NAv; total mass density [g/cm3];to compare with aoki-san
  ztpn(6, *) = ztpn(3,*)*8d-4;10^(INTERPOL(dataCO(1, *), dataCO(0, *), z)) ; CO number density [molec/cm3];to compare with aoki-san
  ztpn(7, *) = ztpn(3,*)*3d-4;10^(INTERPOL(dataH2O(1, *), dataH2O(0, *), z)) ; H2O number density [molec/cm3];to compare with aoki-san H2O (700ppm), CO(800ppm), CO2(0.9532) 
  ztpn(8, *) = (INTERPOL(datadust(1, *), datadust(0, *), z));*100d ; dust number density [molec/cm3];to compare with aoki-san
  ztpn(9, *) = (INTERPOL(dataice(1, *), dataice(0, *), z));*2d;*0.1d ; water ice number density [molec/cm3];to compare with aoki-san

  ztpn(8, where(ztpn(8, *) le 0)) = ztpn(8,min(where(ztpn(8, *) le 0))-1L)*0.1d; to calculate cov(log-dustdensity)  小暮　kogure 20210922
  ztpn(9, where(ztpn(9, *) le 0)) = ztpn(9,min(where(ztpn(9, *) le 0))-1L)*0.1d; to calculate cov(log-dustdensity)  小暮　kogure 20210922

  ;  ztpn(8, *) = ztpn(8, 28);;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;あとで消す！！2021/6/7
  ;  ztpn(9, *) = ztpn(9, 28) ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;あとで消す！！2021/6/7
  ;ztpn(9, *) /= 4 / 3 * !DPI * radius(1)^3 * rhoice / MM_H2O * Nav

  ;  ztpn(1, *) = INTERPOL(data(1, *), data(0, *), z) ; temperature [K]
  ;  ztpn(2, *) = EXP(INTERPOL(ALOG(data(2, *)), data(0, *), z)) ; pressure [Pa]
  ;  ztpn(3, *) = EXP(INTERPOL(ALOG(data(3, *)), data(0, *), z)) ; total number density [molec/cm3]
  ;  ztpn(5, *) =10^(INTERPOL(dataCO2(1, *), dataCO2(0, *), z)) ; CO2 number density [molec/cm3]
  ;  ztpn(4, *) = 44.01D * ztpn(3, *) * ztpn(5, *) / NAv; total mass density [g/cm3]
  ;  ztpn(6, *) = 10^(INTERPOL(dataCO(1, *), dataCO(0, *), z)) ; CO number density [molec/cm3]
  ;  ztpn(7, *) = 10^(INTERPOL(dataH2O(1, *), dataH2O(0, *), z)) ; H2O number density [molec/cm3]
  ;  ztpn(8, *) = (INTERPOL(datadust(1, *), datadust(0, *), z)) ; dust number density [molec/cm3]
  ;    ztpn(9, *) = (INTERPOL(dataice(1, *), dataice(0, *), z)) ; water ice number density [molec/cm3]
  ;  ztpn(9, *) /= 4 / 3 * !DPI * radius(1)^3 * rhoice / MM_H2O * Nav

  ;---------==========------======-------======------====-終わり--------========-----======-------===-
  ;  path_output = pathnames.results
  ;  ;SAVENAME = 'Orbit_' + 'step00_' + orbit + "_" + STRING(i, "(I02)")
  ;  openw,lun,path_output +'GetZTPN.txt',/get_lun
  ;  printf,lun,ztpn
  ;  close,lun
  ;  free_lun,lun


  filenameCO2 = filename + "_stdev_CO2.dat"
  filenameCO = filename + "_stdev_CO.dat"
  filenameH2O = filename + "_stdev_H2O.dat"
  filenamedust = filename + "_stdev_d1.5.dat"
  filenameice = filename + "_stdev_H2O_ice.dat"
  dataCO2 = READ_ASCII(pathnames.model + filenameCO2, COMMENT_SYMBOL = '%') & dataCO2 = dataCO2.Field1
  dataCO = READ_ASCII(pathnames.model + filenameCO, COMMENT_SYMBOL = '%') & dataCO = dataCO.Field1
  dataH2O = READ_ASCII(pathnames.model + filenameH2O, COMMENT_SYMBOL = '%') & dataH2O = dataH2O.Field1
  ;  datadust = READ_ASCII(pathnames.model + filenamedust, COMMENT_SYMBOL = '%') & datadust = datadust.Field1
  ;  dataice = READ_ASCII(pathnames.model + filenameice, COMMENT_SYMBOL = '%') & dataice = dataice.Field1
  cov(1, *) = INTERPOL(dataCO2(1, *), dataCO2(0, *), z)
  cov(2, *) = INTERPOL(dataCO(1, *), dataCO(0, *), z)
  cov(3, *) = INTERPOL(dataH2O(1, *), dataH2O(0, *), z)
  cov(4, *) = cov_dust;INTERPOL(datadust(1, *), datadust(0, *), z) * ztpn(8, *);2021/11/19 kogure risei. 高度方向で一定であるべき。そうでないと密度が小さい高高度のフィッティングがされにくくなる。
  cov(5, *) = cov_ice;INTERPOL(dataice(1, *), dataice(0, *), z) ;2021/11/19 kogure risei. 高度方向で一定であるべき。そうでないと密度が小さい高高度のフィッティングがされにくくなる。

  cov(6, *) = cov_dust_r;INTERPOL([.1D, .1D],[0D, 100D], cov(0, *))
  cov(7, *) = INTERPOL([.25D, .25D],[0D, 100D], cov(0, *))
  cov(8, *) = cov_ice_r;INTERPOL([.1D, .1D],[0D, 100D], cov(0, *))
  cov(9, *) = INTERPOL([.05D, .05D],[0D, 100D], cov(0, *))

  RETURN, ztpn
END
