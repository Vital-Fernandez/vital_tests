'''
Created on Feb 23, 2012

@author: william
'''

import numpy as np
import atpy
import os
import gzip
import bz2

def _getStarlightVersion():
#FIXME: Version 04 and 05 are scrambled. =(
#Example: File 1 -
    ###################################################################
    ## OUTPUT of StarlightChains_v04.for    [Cid@UFSC - 10/May/2008] ##
    ###################################################################
#         File 2 -
    ###################################################################
    ## OUTPUT of StarlightChains_v05.for    [Cid@UFSC - 10/May/2008] ##
    ###################################################################

    return np.array(['StarlightChains_v04.for', 'StarlightChains_v05.for'])

def read_set(self, filename, read_header=True, read_chains=False, read_spec=True, verbose=False):
    '''
        Returns an array w/ outputfile from StarlightChains_v05
                    (http://starlight.ufsc.br/)
    '''


    if not os.path.exists(filename):
        raise Exception('File not found: %s' % filename)

    if filename.endswith('.gz'):
        f = gzip.GzipFile(filename)
    elif filename.endswith('.bz2'):
        f = bz2.BZ2File(filename)
    else:
        f = open(filename)
    data = f.readlines()
    f.close()

    fileVersion = data[1].split()[3]
    #if fileVersion != _getStarlightVersion():
    if ( fileVersion != _getStarlightVersion() ).all():
        raise Exception('Incorrect starlight version. Expected: %s, found: %' (_getStarlightVersion(), fileVersion))

    self.reset()

    self.keywords['file_version'] = fileVersion
    self.keywords['arq_synt'] = os.path.basename(filename)

    if(read_header):
        if(verbose): print('read_header set to TRUE! Rea ding the spectra...')

        ## Some input info

        self.keywords['arq_spec']    = np.str(data[5].split()[0])
        self.keywords['arq_base']    = np.str(data[6].split()[0])
        self.keywords['arq_masks']    = np.str(data[7].split()[0])
        self.keywords['arq_config']  = np.str(data[8].split()[0])
        self.keywords['N_base']        =    np.int(data[9].split()[0])
        self.keywords['N_YAV_components']  =    np.int(data[10].split()[0])
        self.keywords['iFitPowerLaw']      =    np.int(data[11].split()[0])
        self.keywords['alpha_PowerLaw']    =  np.float32(data[12].split()[0])
        self.keywords['red_law_option']    = np.str(data[13].split()[0])
        self.keywords['q_norm']            =  np.float32(data[14].split()[0])

        ## (Re)Sampling Parameters

        self.keywords['l_ini']           = np.float32(data[17].split()[0])
        self.keywords['l_fin']           = np.float32(data[18].split()[0])
        self.keywords['dl']              = np.float32(data[19].split()[0])

        ## Normalization info

        self.keywords['l_norm']          = np.float32(data[22].split()[0])
        self.keywords['llow_norm']       = np.float32(data[23].split()[0])
        self.keywords['lupp_norm']       = np.float32(data[24].split()[0])
        self.keywords['fobs_norm']       = np.float32(data[25].split()[0])

        ## S/N

        self.keywords['llow_SN']         = np.float32(data[28].split()[0])
        self.keywords['lupp_SN']         = np.float32(data[29].split()[0])
        self.keywords['SN_snwin']        = np.float32(data[30].split()[0])
        self.keywords['SN_normwin']      = np.float32(data[31].split()[0])
        self.keywords['SNerr_snwin']     = np.float32(data[32].split()[0])
        self.keywords['SNerr_normwin']   = np.float32(data[33].split()[0])
        self.keywords['fscale_chi2']     = np.float32(data[34].split()[0])

        ## etc...

        self.keywords['idum_orig']       = np.int(data[37].split()[0])
        self.keywords['NOl_eff']         = np.int(data[38].split()[0])
        self.keywords['Nl_eff']          = np.int(data[39].split()[0])
        self.keywords['Ntot_clipped']    = np.int(data[40].split()[0])
        self.keywords['Nglobal_steps']   = np.int(data[41].split()[0])
        self.keywords['N_chains']        = np.int(data[42].split()[0])
        self.keywords['NEX0s_base']      = np.int(data[43].split()[0])

        ## Synthesis Results - Best model ##

        self.keywords['chi2']            = np.float32(data[49].split()[0])
        self.keywords['adev']            = np.float32(data[50].split()[0])

        self.keywords['sum_x']           = np.float32(data[52].split()[0])
        self.keywords['Flux_tot']        = np.float32(data[53].split()[0])
        self.keywords['Mini_tot']        = np.float32(data[54].split()[0])
        self.keywords['Mcor_tot']        = np.float32(data[55].split()[0])

        self.keywords['v_0']             = np.float32(data[57].split()[0])
        self.keywords['v_d']             = np.float32(data[58].split()[0])
        self.keywords['A_V']             = np.float32(data[59].split()[0])
        self.keywords['YA_V']            = np.float32(data[60].split()[0])


        # Read/define x, mu_ini, mu_cor, age_base, Z_base & YAV_flag arrays.
        _nlast = 62+self.keywords['N_base']
        # Reset populations lists
        #popx 2 popmu_ini 3 popmu_cor 4 popage_base 5 popZ_base 6  popYAV_flag 8 popMstars 9
        popx             = []    # column 2
        popmu_ini        = []    # column 3
        popmu_cor        = []    # column 4
        popage_base      = []    # column 5
        popZ_base        = []    # column 6
        popfbase_norm    = []    # column 7
        popYAV_flag      = []    # column 8
        popMstars        = []    # column 9


        for i in range(63,_nlast+1):
            tmp = data[i].split()
            popx.append(          np.float32(tmp[1]) )
            popmu_ini.append(     np.float32(tmp[2]) )
            popmu_cor.append(     np.float32(tmp[3]) )
            popage_base.append(   np.float32(tmp[4]) )
            popZ_base.append(     np.float32(tmp[5]) )
            popfbase_norm.append( np.float32(tmp[6]) )
            popYAV_flag.append(   np.float32(tmp[7]) )
            popMstars.append(     np.float32(tmp[8]) )

        t = atpy.Table()
        t.table_name = 'population'
        t.add_column('popx', np.array(popx), dtype='>f8')
        t.add_column('popmu_ini', np.array(popmu_ini), dtype='>f8')
        t.add_column('popmu_cor', np.array(popmu_cor), dtype='>f8')
        t.add_column('popage_base', np.array(popage_base), dtype='>f8')
        t.add_column('popZ_base', np.array(popZ_base), dtype='>f8')
        t.add_column('popfbase_norm', np.array(popfbase_norm), dtype='>f8')
        t.add_column('popYAV_flag', np.array(popYAV_flag), dtype='>f8')
        t.add_column('popMstars', np.array(popMstars), dtype='>f8')
        self.append(t)

        # Renormalize x to 100% sum!!
        #FIXME: Why not?
        #pop[0] = 100.*pop[0]/np.sum(pop[0])


        #FIXME: Skip this also...
        # OBS: PL have age = 0 in the Starlight output file:(
        #      Here I change it so that age_PL = 5e5 yr... & Z_PL = solar
        #      This is all obsolete anyway. The built-in PL is not used anymore.
        #if (int(self.keywords['iFitPowerLaw']) > 0):
        #    print '@@> [Warning!] ...Fixing PL age & Z ...????? CHECK THIS ?????'
        #    pop[3][self.keywords['N_base'] - 1] = 5e5 #popage_bae
        #    pop[4][self.keywords['N_base'] - 1]   = 0.02 #popZ_base

        #self.keywords['pop'] = pop


    if(read_spec == True):
        if(verbose): print('read_spec set to TRUE! Reading the spectra...')
        if(read_header == False):
                self.keywords['N_base']      =    np.int(data[9].split()[0])
                self.keywords['fobs_norm']   = np.float32(data[25].split()[0])

        # Read spectra (l_obs, f_obs, f_syn & f_wei)
        #l_obs 1 f_obs 2 f_syn 3 f_wei 4 Best_f_SSP 5
        iaux1 = 62 + self.keywords['N_base'] + 5 + self.keywords['N_base'] + 2 + self.keywords['N_base'] + 11
        self.keywords['Nl_obs'] = np.int(data[iaux1].split()[0])
        iaux2 = iaux1 + 1
        #iaux3 = iaux1 + StarlightOut['Nl_obs']

        try:
            dt = np.dtype([('wl', 'float32'), ('flux_obs', 'float32'), ('flux_syn', 'float32'), ('wei', 'float32'), ('Best_f_SSP', 'float32')])
            out_spec = np.loadtxt(filename, dtype=dt, skiprows=iaux2)
        except:
            if(verbose): print('Warning: Did not read Best_f_SSP!')
            dt = np.dtype([('wl', 'float32'), ('flux_obs', 'float32'), ('flux_syn', 'float32'), ('wei', 'float32')])
            out_spec = np.loadtxt(filename, dtype=dt, skiprows=iaux2)

        t = atpy.Table()
        t.table_name = 'spectra'
        t.add_column('l_obs', out_spec['wl'], dtype='>f8')
        t.add_column('f_obs', out_spec['flux_obs'], dtype='>f8')
        t.add_column('f_syn', out_spec['flux_syn'], dtype='>f8')
        t.add_column('f_wei', out_spec['wei'], dtype='>f8')
        if (len(out_spec.dtype.names) == 5):
            t.add_column('Best_f_SSP', np.array(out_spec['Best_f_SSP']), dtype='>f8')

        self.append(t)
