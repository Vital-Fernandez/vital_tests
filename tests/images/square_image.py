from PIL import Image


def make_square(image_path, output_path):
    # Open the image
    img = Image.open(image_path)

    # Calculate the size for the new square image
    size = max(img.size)

    # Create a new square image with a transparent background
    new_img = Image.new("RGBA", (size, size), (255, 255, 255, 0))

    # Calculate the position to paste the original image
    x = (size - img.width) // 2
    y = (size - img.height) // 2

    # Paste the original image onto the new square image
    new_img.paste(img, (x, y))

    # Save the new image
    new_img.save(output_path)


# Example usage
input_image = '/home/vital/PycharmProjects/Vital-Fernandez.github.io/assets/images/Profile_picture_avatar.png'
output_image = '/home/vital/PycharmProjects/Vital-Fernandez.github.io/assets/images/Profile_picture_avatar_squared.png'
make_square(input_image, output_image)
