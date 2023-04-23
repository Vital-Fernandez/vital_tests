
import pandas as pd
from pandas_indexing import modifying_df, XY

# Initialize data of lists
data = [{'a':4, 'b': 2, 'c': 3}, {'a': 10, 'b': 20, 'c': 30}]

# Creates pandas DataFrame by passing
df = pd.DataFrame(data, index=['first', 'second'])

print(df)
modifying_df(df)
print(df)

x, y = None, 1
my_values = XY(x, y)
print(my_values.x, my_values.y)
my_values.swap()
print(my_values.x, my_values.y)
print(x, y)


