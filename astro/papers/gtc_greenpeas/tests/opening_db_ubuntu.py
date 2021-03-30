import pickle
import src.specsiser as sr
from pathlib import Path

file_address = Path('/home/vital/gp030321_BR_Direct-Teff-logU_it3.db')

# results_dict = sr.load_MC_fitting(file_address)

# Restore the pymc output file
with open('/home/vital/gp030321_BR_Direct-Teff-logU_it3.db', 'wb') as trace_restored:
    db = pickle.load(trace_restored)



exp_number = 3
exp_time = 1200
grism_b = R1000B
grism_r = R1000R
seeing = 0000
standard_blue = 0000
standard_red = 0000
night_type = 0000
slit_obj = 0.8"
slit_star = 2.5"

'''
Para el paper: 


Lámparas: HgAr (72957), Ne (image 72958), and Xe (image 72959). Confirmo entonces que es igual que en A+12

Para GP03030: 3x1200s (en cada grating) seeing~0.8"-0.9" dark night
Estándares: G158-100 with R1000 + slit 2.5" Paralactic. Tomada como baseline por el observatorio (fecha diferente al objeto)





  Para GP101157; 3x1200s, dark night seeing ~1"  
Estándar: Para el rojo, Ross640 with R1000R + s0.8" (notar que es la rendija que se usó en el objeto, y no la de 2.5" que se usa para las calibraciones...)
PAra el azul, G191-B2B with R1000B + s2.5" 





PAra GP121910: 3x1200s seeing ~0.9" dark night
Estándar: Azul Feige 34  with R1000B + s2.5" 
Rojo,  Ross640 with R1000R + s0.8" (notar que es la rendija que se usó en el objeto, y no la de 2.5" que se usa para las calibraciones...)



Al final logré encontrar lo que hice para combinar los dos brazos azul y rojo... usé scombine con mediana y sin reject, usando la estadística  entre ~5000-7000 (region de solapamiento). 
Seguramente para GP101157 y GP1219 la calibración en flujo tiene mayor error por combinar dos estrellas tomadas con diferente apertura? flux loss en el rojo? 

'''
