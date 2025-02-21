import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define cities with (longitude, latitude)
cities = {
    "London": (-0.1278, 51.5074),
    "Toulouse": (1.4442, 43.6047),
    "Madrid": (-3.7038, 40.4168),
    "La Puebla": (-98.2062, 19.0414),
    "La Serena": (-71.2520, -29.9045),
    "Ann Arbor": (-83.7430, 42.2808)
}


# Create figure and axis with a Robinson projection
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson())

# Ensure the whole world is shown
ax.set_global()

# Set background color to black
ax.set_facecolor("#2B2B2B")

# Add coastlines and borders in white
ax.add_feature(cfeature.COASTLINE, edgecolor='#CCCCCC')
ax.add_feature(cfeature.BORDERS, edgecolor='#CCCCCC')
ax.add_feature(cfeature.LAKES, edgecolor="#CCCCCC", facecolor="#2B2B2B")

# Plot each city
for city, (lon, lat) in cities.items():
    ax.plot(lon, lat, marker="o", color="#d04350", markersize=5, transform=ccrs.PlateCarree())
    # ax.text(lon + 3, lat, city, color="#CCCCCC", fontsize=10, transform=ccrs.PlateCarree())

# # Show the plot
# plt.show()
# Save the figure with transparent outer edges
ax.spines['geo'].set_edgecolor('#CCCCCC')
plt.savefig("/home/vital/Desktop/world_map.png", dpi=300, bbox_inches='tight', transparent=True)
