import pandas as pd
import timeit


import timeit
import_module = "import pandas as pd"
testcode = ''' 
def test(): 
    df = pd.DataFrame(columns=["A", "B", "C"])
    new_row = pd.Series({"A": 4, "B": 4, "C": 4})
    for i in range(16000):
        df[i] = new_row
        
'''
print(timeit.repeat(stmt=testcode, setup=import_module))
