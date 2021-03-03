import numpy as np
import pandas as pd
from uncertainties import ufloat
from uncertainties import unumpy
from astro.ext_lib.bces_script import bces_functions as bc
from matplotlib import pyplot as plt, rcParams, spines

def load_excel_DF(frame_address):
    # File which stores all the original file sheets #WARNING: this will not work if more than one excel DF loaded
    ipExcel_sheetColumns = {}

    # Load excel file:
    with pd.ExcelFile(frame_address) as xlsx_file:

        # Load all sheets
        list_Df_sheet_i, sheets_names = [], xlsx_file.sheet_names
        for sheet in sheets_names:
            df_i = xlsx_file.parse(sheet, index_col=0)
            list_Df_sheet_i.append(df_i)
            ipExcel_sheetColumns[sheet] = list(df_i.columns.values)

    # Combine individual sheet-df into one
    df = pd.concat(list_Df_sheet_i, axis=1)

    # Combine nominal and error columns for the same variables
    df_columns = df.columns.values
    for column in df_columns:
        if column + '_err' in df_columns:
            # Scheme only to combine rows which only contain value
            idcs_nan = df[column + '_err'].isnull()

            # Empty error cells produce simple floats
            df.loc[idcs_nan, column] = df.loc[idcs_nan, column + '_err']

            # Otherwise they produce uarray
            df.loc[~idcs_nan, column] = unumpy.uarray(df.loc[~idcs_nan, column], df.loc[~idcs_nan, column + '_err'])

            # df[column] = unumpy.uarray(df[column].values, df[column + '_err'].values)

            # Remove error column from the dataframe
            df.drop(column + '_err', axis=1, inplace=True)

    return df

def quick_indexing(df):

    df['quick_index'] = np.nan

    counter = 1
    for obj in df.index:

        if df.loc[obj, 'Ignore_article'] != 'yes':

            if pd.notnull(df.loc[obj, 'Favoured_ref']):
                df.loc[obj, 'quick_index'] = df.loc[obj, 'Favoured_ref']
            else:
                df.loc[obj, 'quick_index'] = "FTDTR-" + str(counter)
                counter += 1

    idx_include = pd.notnull(df['quick_index'])

    return

size_dict = {'axes.labelsize':35, 'legend.fontsize':24, 'font.family':'Times New Roman', 'mathtext.default':'regular', 'xtick.labelsize':30, 'ytick.labelsize':30}
rcParams.update(size_dict)

catalogue_df = load_excel_DF('D:/Dropbox/Astrophysics/Data/PhD_HIIGalaxies/WHT_Galaxies_properties.xlsx')

quick_indexing(catalogue_df)

idcs = (~catalogue_df.TeOIII_emis2nd.isnull() & ~catalogue_df.TeSIII_emis2nd.isnull() & catalogue_df.quick_index.notnull())

TeOIII_array = catalogue_df.loc[idcs].TeOIII_emis2nd.values
TeSIII_array = catalogue_df.loc[idcs].TeSIII_emis2nd.values
objects = catalogue_df.loc[idcs].quick_index.values




#Make the plot
x_regression = np.linspace(0.8 * np.min(unumpy.nominal_values(TeOIII_array)), 1.20 * np.max(unumpy.nominal_values(TeOIII_array)), 10)


# for i in range(len(regr_dict['m'])):


y_regression_Garnet92 = (0.83 * x_regression/10000 + 0.17) * 10000
y_regression_EpmDiaz05 = (1.05 * x_regression/10000 - 0.08) * 10000
y_regression_Epm2014 = (0.92 * x_regression/10000 + 0.078) * 10000
y_regression_Hagele2006 = (1.19 * x_regression/10000 - 0.32) * 10000

regr_dict = bc.bces_regression(unumpy.nominal_values(TeOIII_array), unumpy.nominal_values(TeSIII_array), unumpy.std_devs(TeOIII_array), unumpy.std_devs(TeSIII_array))

reg_code = 3
y_fit = regr_dict['m'][reg_code] * x_regression + regr_dict['n'][reg_code]
print(f'm = {regr_dict["m"][reg_code]}, n = {regr_dict["n"][reg_code]/10000}')

fig, ax = plt.subplots(figsize=(14, 8))

TOIII_values = unumpy.nominal_values(TeOIII_array)
TSIII_values = unumpy.nominal_values(TeSIII_array)
ax.scatter(TOIII_values, TSIII_values, label='HII galaxies')
ax.plot(x_regression, y_regression_Garnet92, label='Garnett (1992)', linestyle=':')
ax.plot(x_regression, y_regression_EpmDiaz05, label=r'$P\'erez$-Montero et al (2005)', linestyle='--')
ax.plot(x_regression, y_regression_Epm2014, label=r'$P\'erez$-Montero (2014)', linestyle='-.')
ax.plot(x_regression, y_regression_Hagele2006, label=r'$Haegele (2006)', linestyle='--')
ax.plot(x_regression, y_fit, label='Linear fit', linestyle = '-')


ax.yaxis.set_ticks(np.arange(8000, 26000, 4000))
ax.update({'ylabel': r'$T_{e}[SIII]\,(K)$', 'xlabel': r'$T_{e}[OIII]\,(K)$', 'title': 'TSIII - TOIII model comparison'})
ax.legend()
plt.show()