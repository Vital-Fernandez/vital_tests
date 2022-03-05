import numpy as np
import pandas as pd

print(f'Pandas {pd.__version__}')

columns = ['c0', 'c1', 'c2', 'c3']
mySeries = pd.Series(index=columns)

mySeries['c0'] = 'None'
mySeries['c1'] = np.nan
mySeries['c2'] = 1234.0
mySeries['c3'] = 1.234e-18

print(mySeries)
