import pandas as pd

# Create individual pandas DataFrame.
df1 = pd.DataFrame({'Col1': [1, 2, 3, 4], 'Col2': [99, 98, 95, 90]}, index=['A', 'B', 'C', 'D'])
df2 = pd.DataFrame({'Col1': [1, 2], 'Col2': [99, 98]}, index=['A', 'B'])
df3 = pd.DataFrame({'Col1': [3, 4], 'Col2': [95, 90]}, index=['C', 'D'])
df4 = pd.DataFrame({'Col1': [3, 4], 'Col2': [95, 90]}, index=['B', 'C'])

# Combine into one multi-index dataframe
df_dict = dict(obj1=df1, obj2=df2, obj3=df3, obj4=df4)

# Assign multi-index labels
mDF = pd.concat(list(df_dict.values()), keys=list(df_dict.keys()))
mDF.rename_axis(index=["ID", "property"], inplace=True)
print(mDF, '\n')

idcs_isin = mDF.index.get_level_values('property').isin(['A', 'B'])
idcs_query = mDF.query('property in ["A","B"]')
print(f'isin:\n{mDF.loc[idcs_isin]}\n')
print(f'Query:\n{idcs_query}')

bools = mDF.index.get_level_values('property').isin(['A','B'])
grouper = mDF.index.get_level_values('ID')
# there should be a minimum of two (`A`, `B`)
bools = pd.Series(bools).groupby(grouper).transform('sum').ge(2).array
print(mDF.loc[bools])