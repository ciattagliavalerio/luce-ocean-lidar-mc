"""
qr_repo.py - QR code of the code repository for the poster footer.

Dotted modules in the frame blue, error correction H, centred lidar icon
(satellite, green beam, sea) covering 38% of the side. Decoding checked
with OpenCV on the full image and at 180 px.

Usage: python3 qr_repo.py   (needs qrcode[pil], numpy; opencv optional)
Writes: poster/figures/qr_repo.png
"""
import os
import numpy as np
import qrcode
from qrcode.image.styledpil import StyledPilImage
from qrcode.image.styles.moduledrawers.pil import CircleModuleDrawer, RoundedModuleDrawer
from qrcode.image.styles.colormasks import SolidFillColorMask
from PIL import Image, ImageDraw

URL = "https://github.com/ciattagliavalerio/luce-ocean-lidar-mc"
FRAME = (18, 53, 91)    # luceframe.
GREEN = (35, 164, 67)   # lucestep1.
BLUE = (0, 106, 255)    # lucestep2.
ICON_FRAC = 0.38        # Icon diameter / QR side.


def icon(n):
    """Satellite with a green beam onto the sea, on a white disc."""
    im = Image.new("RGBA", (n, n), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    cx = n / 2
    d.polygon([(cx - n * 0.03, n * 0.36), (cx + n * 0.03, n * 0.36),
               (cx + n * 0.20, n * 0.70), (cx - n * 0.20, n * 0.70)],
              fill=GREEN + (230,))
    d.rectangle([cx - n * 0.07, n * 0.20, cx + n * 0.07, n * 0.36], fill=FRAME)
    for s in (-1, 1):
        x0 = cx + s * n * 0.09
        x1 = cx + s * n * 0.30
        d.rectangle([min(x0, x1), n * 0.24, max(x0, x1), n * 0.32], fill=BLUE)
        d.line([(cx + s * n * 0.07, n * 0.28), (x0, n * 0.28)], fill=FRAME,
               width=n // 60)
    for k, y in enumerate((0.70, 0.79)):
        pts = [(x, n * y + n * 0.025 * np.sin((x / n) * 8 * np.pi + k))
               for x in np.linspace(n * 0.14, n * 0.86, 120)]
        d.line(pts, fill=BLUE if k == 0 else FRAME, width=n // 28)
    base = Image.new("RGBA", (n, n), (0, 0, 0, 0))
    ImageDraw.Draw(base).ellipse([0, 0, n - 1, n - 1], fill="white")
    return Image.alpha_composite(base, im)


def main():
    qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_H,
                       box_size=40, border=0)
    qr.add_data(URL)
    qr.make(fit=True)
    img = qr.make_image(image_factory=StyledPilImage,
                        module_drawer=CircleModuleDrawer(),
                        eye_drawer=RoundedModuleDrawer(radius_ratio=1.0),
                        color_mask=SolidFillColorMask(front_color=FRAME))
    img = img.convert("RGB")
    w = img.size[0]
    ic = icon(int(w * ICON_FRAC))
    img.paste(ic, ((w - ic.size[0]) // 2,) * 2, ic)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                       "poster", "figures", "qr_repo.png")
    img.save(out, dpi=(900, 900))
    print("written", out, img.size)
    try:
        import cv2
        a = cv2.copyMakeBorder(cv2.imread(out), 160, 160, 160, 160,
                               cv2.BORDER_CONSTANT, value=(255, 255, 255))
        det = cv2.QRCodeDetector()
        small = cv2.resize(a, (200, 200), interpolation=cv2.INTER_AREA)
        print("decode full:", det.detectAndDecode(a)[0] == URL,
              "decode 200 px:", det.detectAndDecode(small)[0] == URL)
    except ImportError:
        print("OpenCV not available: decoding not checked")


if __name__ == "__main__":
    main()
