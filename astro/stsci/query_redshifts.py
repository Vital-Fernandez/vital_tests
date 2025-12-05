from astroquery.simbad import Simbad
from astroquery.ipac.ned import Ned
import toml

# List of galaxies
# galaxies = [
#     "Haro 11", "Haro 2", "He 2-10", "SBS 1159+545", "NGC 1705",
#     "Pox 186", "NGC 4861", "SBS 1415+437", "SBS 0335-052",
#     "I Zw 18", "NGC 2366", "UM 461", "UGC 4483", "UGCA 281", "VII Zw 403"
# ]

galaxies = ['MRK1450']

# Configure Simbad to return redshift (z)
Simbad.add_votable_fields("rvz_redshift")

# Dictionary to store results
redshifts = {}

for name in galaxies:
    result = Ned.query_object(name)
    if result is not None and "Redshift" in result.colnames:
        z = result["Redshift"][0]
        redshifts[name] = float(z) if z is not None else None
    else:
        redshifts[name] = None
    print(name, redshifts[name])

# # Save results to TOML
# with open("galaxy_redshifts.toml", "w") as f:
#     toml.dump(redshifts, f)
#
# print("Redshift values saved to galaxy_redshifts.toml")