FUNCTION CHECK_AEROSOLS

    radiusM_eff = 1D-4
    radiusS_eff = 0.3
    data = GET_AEROSOLS_OPTICAL('dust', radiusM_eff, radiusS_eff)
END