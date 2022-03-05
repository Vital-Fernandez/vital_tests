import numpy as np

a = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]).astype(float)

mask_word = '3-7,12.2,15-18'

b = np.full(a.size, True)
for entry in mask_word.split(','):
    if '-' in entry:
        entry = np.array(entry.split('-')).astype(float)
        mask_entry = (a > entry[0]) & (a < entry[1])
        b = b * mask_entry
    else:
        b[np.searchsorted(a, float(entry))] = False

print(a)
print(b)

