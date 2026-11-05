from PIL import Image
im=Image.open("redock_overlay_labeled.png")
print(im.size, im.info.get("dpi"))
