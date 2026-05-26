#!/usr/bin/env python3
"""vscode_icons.py  –  Re-colour gmgpy toolbar icons to VS Code monochrome style.

Every visible pixel's RGB is replaced with #C5C5C5 (VS Code icon grey) while
the alpha channel is preserved, so anti-aliased edges remain correct on a dark
toolbar background.  The `off_F` disabled icon is dimmed to #646464.

Run from the repository root:
    conda run -n pygmt python src/gmgpy/icons/vscode_icons.py
"""

from PIL import Image
import os

ICONS_DIR = os.path.dirname(os.path.abspath(__file__))

# VS Code icon colours
GREY    = (197, 197, 197)   # #C5C5C5  — normal icons
DIM     = (100, 100, 100)   # #646464  — disabled / off-state icons

# (filename, colour)
ICONS = [
    # ── toolbar 24 px ──────────────────────────────────────────────────────
    ('save_24.png',            GREY),
    ('load_24.png',            GREY),
    ('G_24.png',               GREY),
    ('V_24.png',               GREY),
    ('M_24.png',               GREY),
    ('C_24.png',               GREY),
    ('large_up_24.png',        GREY),
    ('large_down_24.png',      GREY),
    ('small_up_24.png',        GREY),
    ('small_down_24.png',      GREY),
    ('zoom_in_24.png',         GREY),
    ('zoom_out_24.png',        GREY),
    ('full_extent_24.png',     GREY),
    ('pan_24.png',             GREY),
    ('left_small_24.png',      GREY),
    ('right_small_24.png',     GREY),
    ('large_left_24.png',      GREY),
    ('large_right_24.png',     GREY),
    ('well_24.png',            GREY),
    ('F_24.png',               GREY),
    ('off_F_24.png',           DIM),   # disabled state
    ('faultline_24.png',       GREY),
    ('undo_24.png',            GREY),
    ('pinch_24.png',           GREY),
    # ── statusbar 16 px ────────────────────────────────────────────────────
    ('G_16.png',               GREY),
    ('V_16.png',               GREY),
    ('M_16.png',               GREY),
    ('T_16.png',               GREY),
    ('large_up_16.png',        GREY),
    ('python_16.png',          GREY),
]


def _recolour(path: str, rgb: tuple) -> None:
    img = Image.open(path).convert('RGBA')
    r, g, b = rgb
    pixels = [(r, g, b, a) for _, _, _, a in img.getdata()]
    img.putdata(pixels)
    img.save(path)
    print(f'  {os.path.basename(path)}')


if __name__ == '__main__':
    print('Re-colouring icons to VS Code style...')
    missing = []
    for name, colour in ICONS:
        path = os.path.join(ICONS_DIR, name)
        if os.path.exists(path):
            _recolour(path, colour)
        else:
            missing.append(name)
    if missing:
        print(f'\nNot found (skipped): {missing}')
    print('Done.')
