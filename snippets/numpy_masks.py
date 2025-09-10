# import numpy as np
#
# # a = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]).astype(float)
# #
# # mask_word = '3-7,12.2,15-18'
# #
# # b = np.full(a.size, True)
# # for entry in mask_word.split(','):
# #     if '-' in entry:
# #         entry = np.array(entry.split('-')).astype(float)
# #         mask_entry = (a > entry[0]) & (a < entry[1])
# #         b = b * mask_entry
# #     else:
# #         b[np.searchsorted(a, float(entry))] = False
# #
# # print(a)
# # print(b)
# #
#
# import numpy as np
# x = np.array([1, 2, 3, -1, 5])
# mx = np.ma.masked_array(x, mask=[0, 0, 0, 1, 0])
# idcs_array = [False, False, True, True, False]
#
# print(mx)
# print(mx[idcs_array])
# print(np.ma.is_masked(x))
# print(np.ma.is_masked(mx))

# In[1]:


import lmfit
from lmfit.models import LinearModel
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


N = 100
x = np.arange(N)
a = np.random.randn(N)
a[:80] +=20
b = np.ma.masked_array(a, a < 10)

lm = LinearModel()
result = lm.fit(a, x=x) # wrong
result.plot()
plt.show()

mask = b.mask
w_masked = np.ma.masked_array(np.ones(x.size), mask)
x2 = np.ma.masked_array(x, mask)
result2 = lm.fit(b, x=x2, weights=w_masked) # correct
result2.plot()
plt.show()