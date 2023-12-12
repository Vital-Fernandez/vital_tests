import numpy as np
dist = 26.7
spaxel_size = 0.2

scale = 1 / 206265 * (dist * 1e6)
spaxel_size_pc = scale * spaxel_size
spaxel_area = (scale * spaxel_size)**2

print(f'Muse spaxel size for CGCG007 025: {spaxel_size_pc} pc')
print(f'Muse spaxel area for CGCG007 025: {spaxel_area} pc^2')