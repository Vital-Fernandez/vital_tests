import numpy as np
import pyCloudy as pc

#Files to exclude from the default cloudy output
components_remove = [#['radius', '.rad'],
                     #['continuum', '.cont'],
                     #['physical conditions', '.phy'],
                     ['overview', '.ovr'],
                     ['heating', '.heat'],
                     ['cooling', '.cool'],
                     ['optical depth', '.opd']]

for item in components_remove:
    pc.config.SAVE_LIST.remove(item)

elements_remove = [['hydrogen','.ele_H'],
                   ['helium','.ele_He'],
                   ['carbon', '.ele_C'],
                   ['nitrogen', '.ele_N'],
                   ['oxygen', '.ele_O'],
                   ['argon', '.ele_Ar'],
                   ['neon', '.ele_Ne'],
                   ['sulphur', '.ele_S'],
                   ['chlorin', '.ele_Cl'],
                   ['iron', '.ele_Fe'],
                   ['silicon', '.ele_Si']]

for item in elements_remove:
    pc.config.SAVE_LIST_ELEMS.remove(item)

#These are the commands common to all the models (here only one ...)
emis_tab = ['H  1 4861.33A',
            'H  1 6562.81A',
            'HE 1 3888.63A',
            'HE 1 4026.20A',
            'HE 1 4471.49A',
            'HE 1 5875.64A',
            'HE 1 6678.15A',
            'HE 1 7065.22A']


# Tell pyCloudy where your cloudy executable is:
pc.config.cloudy_exe = '/home/vital/cloudy/source/cloudy.exe'

# Define some parameters of the model:
dir_ = '/home/vital/Desktop/tests_cloudy/'
model_name = 'emissivity_grid'
full_model_name = '{0}{1}'.format(dir_, model_name)

# Physical parameters
dens    = 2.        #log cm-3
Teff    = 45000.    #K
qH      = 47.0      #s-1
dist    = 16.0      #kpc

def model_TempDen_variation(Te, ne, simu_address):
    logTe = np.log10(Te)
    cloudy_model = pc.CloudyInput('{}M_i_Te{}_ne{}'.format(simu_address, Te, ne))
    cloudy_model.set_other(('init file "hheonly.ini"'))
    cloudy_model.set_other(('element helium abundance -3'))
    cloudy_model.set_BB(Teff=Teff, lumi_unit='q(H)', lumi_value=qH)
    cloudy_model.set_radius(r_in=dist)
    cloudy_model.set_other(('constant temperature {}'.format(logTe)))
    cloudy_model.set_cste_density(dens)
    cloudy_model.set_other(('database H-like levels large element hydrogen'))
    cloudy_model.set_other(('database H-like levels large element helium'))
    cloudy_model.set_other(('set dr 0'))
    cloudy_model.set_stop(('zone 1'))
    cloudy_model.set_emis_tab(emis_tab)  # better use read_emis_file(file) for long list of lines, where file is an external file
    cloudy_model.print_input()
    return cloudy_model

# #Run cloudy
Te_interval = [9000.0, 9500.0, 10000.0, 10500.0]
for Te in Te_interval:
    model_i = model_TempDen_variation(Te, 100.0, dir_)
    pc.log_.timer('Starting Cloudy', quiet = True, calling = 'test1')
    model_i.run_cloudy()
    pc.log_.timer('Cloudy ended after seconds:', calling = 'test1')

#Loading simulation results
Ms = pc.load_models(dir_)

#Reading the physical parameters in vector form
Te_vector = np.empty(len(Ms))
ne_vector = np.empty(len(Ms))
emis_dict = dict.fromkeys(emis_tab, np.empty(len(Ms)))

for i in range(len(Ms)):
    ne_vector[i] = Ms[i].ne
    Te_vector[i] = Ms[i].T0
    for j in range(len(emis_dict)):
        line_label = emis_tab[j].replace(' ','_').replace('.','')
        emis_dict[emis_tab[j]][i] = Ms[i].get_emis(line_label)


