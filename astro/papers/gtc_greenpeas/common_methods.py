import pyneb as pn


def red_corr_HalphaHbeta_ratio(lines_df, default_cHbeta):

    # Normalizing flux
    if 'H1_6563A' in lines_df.index:
        flux_Halpha = lines_df.loc['H1_6563A', 'gauss_flux']
        flux_Hbeta = lines_df.loc['H1_4861A', 'intg_flux']
        halpha_norm = flux_Halpha / flux_Hbeta

        rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
        rc.setCorr(obs_over_theo=halpha_norm / 2.86, wave1=6563., wave2=4861.)
        cHbeta = rc.cHbeta
    else:
        flux_Hbeta = lines_df.loc['H1_4861A', 'intg_flux']
        rc = pn.RedCorr(R_V=3.4, law='G03 LMC', cHbeta=default_cHbeta)
        cHbeta = float(rc.cHbeta)

    return cHbeta, rc