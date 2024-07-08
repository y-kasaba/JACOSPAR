; 配列を文字列に変換するための関数
FUNCTION array_to_string, array
  result = '[' + STRJOIN(STRING(array, FORMAT='(F0.2)'), ', ') + ']'
  RETURN, result
END


PRO call_main_loop_dodependence0_invsaTest_20240601
  ;メイン関数（call_jaconspar）の感度実験を行うために、引数を変えながら何度もメイン関数を呼び出すプログラム。
  ;
  pathnames = PATH_MANAGEMENT()

  lat = 0
  lst = 12
  ls = 0


  IF lat GE 0D THEN $
    IF lst GE 6D AND lst LT 18D THEN $
    IF ls GE 0D AND ls LT 90D THEN filename = "gem-mars-a585_Spring_Northern_Day" $
  ELSE IF ls GE 90D AND ls LT 180D THEN filename = "gem-mars-a585_Summer_Northern_Day" $
  ELSE IF ls GE 180D AND ls LT 240D THEN filename = "gem-mars-a585_Autumn_Northern_Day" $
  ELSE filename = "gem-mars-a585_Winter_Northern_Day" $
  ELSE $
    IF ls GE 0D AND ls LT 90D THEN filename = "gem-mars-a585_Spring_Northern_Night" $
  ELSE IF ls GE 90D AND ls LT 180D THEN filename = "gem-mars-a585_Summer_Northern_Night" $
  ELSE IF ls GE 180D AND ls LT 240D THEN filename = "gem-mars-a585_Autumn_Northern_Night" $
  ELSE filename = "gem-mars-a585_Winter_Northern_Night" $
  ELSE $
    IF lst GT 6D AND lst LT 18D THEN $
    IF ls GE 0D AND ls LT 90D THEN filename = "gem-mars-a585_Spring_Southern_Day" $
  ELSE IF ls GE 90D AND ls LT 180D THEN filename = "gem-mars-a585_Summer_Southern_Day" $
  ELSE IF ls GE 180D AND ls LT 240D THEN filename = "gem-mars-a585_Autumn_Southern_Day" $
  ELSE filename = "gem-mars-a585_Winter_Southern_Day" $
  ELSE $
    IF ls GE 0D AND ls LT 90D THEN filename = "gem-mars-a585_Spring_Southern_Night" $
  ELSE IF ls GE 90D AND ls LT 180D THEN filename = "gem-mars-a585_Summer_Southern_Night" $
  ELSE IF ls GE 180D AND ls LT 240D THEN filename = "gem-mars-a585_Autumn_Southern_Night" $
  ELSE filename = "gem-mars-a585_Winter_Southern_Night"

  ;アプリオリファイルの定義
  ;ztpn_ini = GETZTPN_LimbTest_new(MEAN(geometry(4, *)), MEAN(geometry(5, *)), MEAN(geometry(6, *)), MEAN(geometry(7, *)), z, cov, StartRadius);T[K],P[pa],CO2ndens,,
  filenamePT = pathnames.acs + 'profile_atmos_for_ARSM.htp';filename + "_mean_atmo.dat";to compare with Aoki-san
  filenameCO2 = filename + "_mean_CO2.dat";to compare with aoki-san
  filenameCO = filename + "_mean_CO.dat";to compare with aoki-san
  filenameH2O = filename + "_mean_H2O.dat";to compare with aoki-san

  ;共通引数の設定
  WNRANGE = [2660D, 10600D]
  DWN = 10D
  ONLYFITSPECIES = [3L, 4L]
  FACTORINI = [0D, 0D, 0D, 1D, 1D]
  NCPUS = 6
  PRECISION = 1d
  altitude = [0D,  100D]

  ;--------------------- ここまでは共通 --------------------------;
 
      ;----------------------- 感度実験1 ---------------------------;
    ;保存ディレクトリ名
    NAMEPLUSDIR = '20240708_invsa_1000_1000_10_10_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

    ;orbit number
    orbit = '0034'

    ;エアロゾルアプリオリのファイルパス
    filenamedust = 'LimbTest_Dust.dat'
    filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
    appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

    ;エアロゾルアプリオリの粒径
    StartRadius = [1d-4, 1d-4]

    ;invsa
    cov_dust = 1000d
    cov_ice = 1000d
    cov_dust_r = 10d
    cov_ice_r = 10d
    cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

    ;Levenberg Marquardtパラメータ
    LM_initial_v = 1d3;vの初期値
    LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
    ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
    THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

    ; ログを記録
    path_output = pathnames.results + orbit + '/'
    path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
    IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
    write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

    ;関数の呼び出し
    c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
      WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
      StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
      /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
    
      ;----------------------- 感度実験2 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_1000_1000_40_40_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 1000d
      cov_ice = 1000d
      cov_dust_r = 40d
      cov_ice_r = 40d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験3 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_1000_1000_70_70_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 1000d
      cov_ice = 1000d
      cov_dust_r = 70d
      cov_ice_r = 70d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験4 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_4000_4000_40_40_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 4000d
      cov_ice = 4000d
      cov_dust_r = 40d
      cov_ice_r = 40d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験5 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_4000_4000_160_160_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 4000d
      cov_ice = 4000d
      cov_dust_r = 160d
      cov_ice_r = 160d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験6 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_4000_4000_280_280_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 4000d
      cov_ice = 4000d
      cov_dust_r = 280d
      cov_ice_r = 280d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験7 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_7000_7000_70_70_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 7000d
      cov_ice = 7000d
      cov_dust_r = 70d
      cov_ice_r = 70d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験8 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_7000_7000_280_280_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 7000d
      cov_ice = 7000d
      cov_dust_r = 280d
      cov_ice_r = 280d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験9 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_7000_7000_490_490_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 7000d
      cov_ice = 7000d
      cov_dust_r = 490d
      cov_ice_r = 490d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験10 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_10000_10000_100_100_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 10000d
      cov_ice = 10000d
      cov_dust_r = 100d
      cov_ice_r = 100d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験11 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_10000_10000_400_400_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 10000d
      cov_ice = 10000d
      cov_dust_r = 400d
      cov_ice_r = 400d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)
      ;----------------------- 感度実験12 ---------------------------;
      ;保存ディレクトリ名
      NAMEPLUSDIR = '20240708_invsa_10000_10000_700_700_r1d-4';NAMEPLUSDIR = '20240225_convergetest_1_LimbTestDust_1d-6_NadirTestH2Oice_1d-6'

      ;orbit number
      orbit = '0034'

      ;エアロゾルアプリオリのファイルパス
      filenamedust = 'LimbTest_Dust.dat'
      filenameice = 'NadirTest_H2O_ice_2times.dat';'NadirTest_H2O_ice.dat';;;;;phobos版のJACOSPAR同様、初めのリトリーバル感度をダストと同じにするために、ダストと同じアプリオリプロファイルを用いる。
      appriori_files = [filenamePT, filenameCO2, filenameCO, filenameH2O, filenamedust, filenameice, filename]

      ;エアロゾルアプリオリの粒径
      StartRadius = [1d-4, 1d-4]

      ;invsa
      cov_dust = 1000d
      cov_ice = 1000d
      cov_dust_r = 700d
      cov_ice_r = 700d
      cov_parameters = [cov_dust, cov_ice, cov_dust_r, cov_ice_r]

      ;Levenberg Marquardtパラメータ
      LM_initial_v = 1d3;vの初期値
      LM_update_rate = 10d;vの更新係数 　解から遠かった場合：v/2、 解に近付いた場合：v*2
      ;LM_v_limit = 1d10;vの最大側のlimitの値。この値よりもvが大きくならないようにする。
      THRESHOLD = 1d-4;収束の条件を規定するεに使用する値。小さいほど収束条件が厳しくなる。これまでのデフォルトは1d-4(phobos上のプログラムでは1d-3だったが、充足度合いから判断して1d-4の方が適切)

      ; ログを記録
      path_output = pathnames.results + orbit + '/'
      path_output = STRMID(path_output, 0, STRLEN(path_output) - 1) + nameplusdir + "/"
      IF ~FILE_TEST(path_output, /DIRECTORY) THEN FILE_MKDIR, path_output
      write_log,path_output + '/log.txt', WNRANGE, DWN, ONLYFITSPECIES, FACTORINI, PRECISION, altitude, orbit, appriori_files, StartRadius, LM_initial_v, LM_update_rate, THRESHOLD

      ;関数の呼び出し
      c = call_jacospar_limb_shoheidatabase_lmmethod_20240601_invsatest(orbit, altitude, NAMEPLUSDIR = NAMEPLUSDIR, cov_parameters = cov_parameters,$
        WNRANGE = WNRANGE, DWN = DWN, ONLYFITSPECIES = ONLYFITSPECIES, FACTORINI = FACTORINI, NCPUS = NCPUS, PRECISION = PRECISION,THRESHOLD = threshold, $
        StartRadius = StartRadius, appriori_files = appriori_files, LM_initial_v = LM_initial_v,LM_update_rate = LM_update_rate, $
        /plottotalopticaldepth, /nofov, /makeplot ,/MAKEPRINT,/fitradius)

  
    
END