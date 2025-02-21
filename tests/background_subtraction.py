from PIL import Image

# Open the image and ensure it has an alpha channel
img = Image.open("/home/vital/Desktop/ml_map_sc.png").convert("RGBA")
datas = img.getdata()

newData = []
for item in datas:
    # Check if the pixel is close to black (you can adjust the threshold as needed)
    if item[0] < 10 and item[1] < 10 and item[2] < 10:
        # Replace black with transparent
        newData.append((0, 0, 0, 0))
    else:
        newData.append(item)

img.putdata(newData)
img.save("/home/vital/Desktop/ml_map_sc_inv.png", "PNG")
