import pandas as pd
print(f'pandas version {pd.__version__}')

idx = 'H1_6563A'
data = {'ion': 'H1',
        'wavelength': 6563.0,
        'latex_label': '$6563\AA\,HI$',
        'intgr_flux': 3.128572e-14,
        'dist': 2.8e20,
        'eqw': 1464.05371}

mySeries = pd.Series(index=data.keys())
for param, value in data.items():
    mySeries[param] = value
print(f'\nSeries: \n {mySeries}')

myDF = pd.DataFrame(columns=data.keys())
myDF.loc[idx] = mySeries
print(f'\nDataFrame:\n {myDF}')

# import numpy as np
# import pandas as pd
#
# myDF = pd.DataFrame(columns=['c0', 'c1', 'c2', 'c3'])
#
# myDF.loc['r0', 'c0'] = 'abc'
# myDF.loc['r0', 'c1'] = np.nan
# myDF.loc['r0', 'c2'] = 1234.0
# myDF.loc['r0', 'c3'] = 1.234e-18
#
# print('Output DF\n', myDF)

# # Save the table as a dataframe.
# with open(f'output_DF.txt', 'wb') as output_file:
#     string_DF = myDF.to_string()
#     output_file.write(string_DF.encode('UTF-8'))
#
# print('Output DF\n', myDF)