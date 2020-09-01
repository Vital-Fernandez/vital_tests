import pandas as pd
import numpy as np

headers = ['caso', 'queso', 'quiso', 'coso']

df = pd.DataFrame(columns=headers)
row_i = [1, 2, '3', np.nan]

df.loc['uno'] = [1, 2, '3', np.nan]
df.loc['dos'] = [2.0, 2, 3, None]
df.loc['tres'] = [3.0, 2.0, 3, np.inf]

print(df)