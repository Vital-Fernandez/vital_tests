from astropy.cosmology import wCDM
import astropy.units as u

# Define the cosmology
cosmo = wCDM(H0=73.1 * u.km / u.s / u.Mpc,
             Om0=0.302,
             Ode0=1.0 - 0.302,  # Flat universe: ΩΛ = 1 - Ωm
             w0=-1.01)

# Define redshift
z = 9.23  # example redshift

# Compute lookback time
lookback = cosmo.lookback_time(z)

print(f"Lookback time at z = {z} is {lookback:.2f}")