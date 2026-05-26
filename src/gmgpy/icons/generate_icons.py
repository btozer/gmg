#!/usr/bin/env python3
"""
generate_icons.py — Draw every gmgpy toolbar / status-bar icon from scratch.

Design system
─────────────
• Final sizes : 24 × 24 px  (toolbar)   16 × 16 px (status-bar)
• Canvas size : final × 4   (draw at high-res, downscale → clean edges)
• Colour      : #C5C5C5  (197,197,197)  on transparent background
• Disabled    : #646464  (100,100,100)
• Stroke      : n // 12  (≈ 2 px at 24 px output)
• Margin      : n // 10  (≈ 2.4 px at 24 px output)

Usage
─────
    conda run -n pygmt python src/gmgpy/icons/generate_icons.py
"""

from PIL import Image, ImageDraw, ImageFont
import math, os

DIR  = os.path.dirname(os.path.abspath(__file__))
IC   = (197, 197, 197, 255)    # #C5C5C5  normal
DIM  = (100, 100, 100, 255)    # #646464  disabled
S    = 4                        # scale factor

FONT = '/System/Library/Fonts/Supplemental/Arial Bold.ttf'


# ── helpers ───────────────────────────────────────────────────────────────────

def _canvas(px: int):
    n = px * S
    return Image.new('RGBA', (n, n), (0, 0, 0, 0)), n


def _draw(img): return ImageDraw.Draw(img)


def _save(img, px, name):
    img.resize((px, px), Image.LANCZOS).save(os.path.join(DIR, name))
    print(f'  {name}')


def lw(n): return max(3, n // 12)   # stroke width
def mg(n): return max(6, n // 10)   # outer margin


# ── arrows ────────────────────────────────────────────────────────────────────

def _arrow(direction: str, thick: bool, px: int, name: str, colour=IC):
    """Solid filled arrow (chevron + shaft polygon)."""
    img, n = _canvas(px)
    d = _draw(img)

    pad = mg(n)
    ctr = n // 2
    tip  = pad
    base = n - pad
    span = base - tip

    if thick:
        hw = int(span * 0.42)   # arrowhead half-width
        sw = int(span * 0.19)   # shaft half-width
    else:
        hw = int(span * 0.32)
        sw = int(span * 0.13)

    join = tip + int(span * 0.46)   # where head meets shaft

    pts_up = [
        (ctr,       tip ),
        (ctr + hw,  join),
        (ctr + sw,  join),
        (ctr + sw,  base),
        (ctr - sw,  base),
        (ctr - sw,  join),
        (ctr - hw,  join),
    ]

    def rot(pts, deg):
        a = math.radians(deg)
        ca, sa = math.cos(a), math.sin(a)
        c0 = n / 2
        return [(c0 + (x-c0)*ca - (y-c0)*sa,
                 c0 + (x-c0)*sa + (y-c0)*ca) for x, y in pts]

    rots = {'up': 0, 'down': 180, 'left': 270, 'right': 90}
    d.polygon(rot(pts_up, rots[direction]), fill=colour)
    _save(img, px, name)


# ── magnifier ─────────────────────────────────────────────────────────────────

def _magnifier(sign: str, name: str):
    """sign = '+' or '-'"""
    img, n = _canvas(24)
    d = _draw(img)
    w  = lw(n)

    # glass – offset toward top-left so handle fits
    r  = int(n * 0.30)
    gx, gy = int(n * 0.38), int(n * 0.38)
    d.ellipse([gx-r, gy-r, gx+r, gy+r], outline=IC, width=w)

    # handle — exits bottom-right of glass toward icon corner
    angle = math.radians(45)              # 45° = lower-right in PIL coords
    hx0 = int(gx + r * math.cos(angle))
    hy0 = int(gy + r * math.sin(angle))
    hx1, hy1 = int(n * 0.88), int(n * 0.88)
    d.line([(hx0, hy0), (hx1, hy1)], fill=IC, width=w + 2)

    # symbol inside glass
    arm = int(r * 0.52)
    if sign == '+':
        d.line([(gx-arm, gy), (gx+arm, gy)], fill=IC, width=w)
        d.line([(gx, gy-arm), (gx, gy+arm)], fill=IC, width=w)
    else:  # minus
        d.line([(gx-arm, gy), (gx+arm, gy)], fill=IC, width=w)

    _save(img, 24, name)


# ── full extent (expand / fit-to-window) ──────────────────────────────────────

def _full_extent():
    img, n = _canvas(24)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)
    arm = int(n * 0.20)   # length of each bracket leg

    for bx, by, dx, dy in [
        (pad,   pad,   +1, +1),   # top-left
        (n-pad, pad,   -1, +1),   # top-right
        (pad,   n-pad, +1, -1),   # bottom-left
        (n-pad, n-pad, -1, -1),   # bottom-right
    ]:
        d.line([(bx, by), (bx + dx*arm, by)],        fill=IC, width=w)
        d.line([(bx, by), (bx,          by + dy*arm)], fill=IC, width=w)

    _save(img, 24, 'full_extent_24.png')


# ── pan (4-way move) ──────────────────────────────────────────────────────────

def _pan():
    img, n = _canvas(24)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)
    cx0 = cy0 = n // 2
    arr = int(n * 0.17)    # arrowhead height
    gap = int(n * 0.13)    # gap around centre

    def arrowhead(tx, ty, dx, dy):
        """Filled triangle arrowhead with tip at (tx,ty), direction (dx,dy)."""
        hw = int(n * 0.09)
        px0, py0 = tx - dx*arr, ty - dy*arr     # base centre
        perp_x, perp_y = -dy, dx
        d.polygon([
            (tx, ty),
            (px0 + perp_x*hw, py0 + perp_y*hw),
            (px0 - perp_x*hw, py0 - perp_y*hw),
        ], fill=IC)

    # Stems
    for (x1,y1),(x2,y2) in [
        ((cx0, cy0-gap), (cx0, pad+arr+w)),      # N stem
        ((cx0, cy0+gap), (cx0, n-pad-arr-w)),    # S stem
        ((cx0-gap, cy0), (pad+arr+w, cy0)),      # W stem
        ((cx0+gap, cy0), (n-pad-arr-w, cy0)),    # E stem
    ]:
        d.line([(x1,y1),(x2,y2)], fill=IC, width=w)

    # Arrowheads
    arrowhead(cx0, pad,    0, -1)   # N
    arrowhead(cx0, n-pad,  0, +1)   # S
    arrowhead(pad, cy0,   -1,  0)   # W
    arrowhead(n-pad, cy0,  1,  0)   # E

    _save(img, 24, 'pan_24.png')


# ── crosshair (capture coordinates) ──────────────────────────────────────────

def _crosshair():
    img, n = _canvas(24)
    d = _draw(img)
    w   = lw(n)
    cx0 = cy0 = n // 2
    r   = int(n * 0.24)    # circle radius
    gap = int(r * 0.6)     # gap between circle and crosshair lines
    pad = mg(n)

    d.ellipse([cx0-r, cy0-r, cx0+r, cy0+r], outline=IC, width=w)

    # 4 radial line segments (with gap around circle edge)
    d.line([(cx0, pad),    (cx0, cy0-r-gap)], fill=IC, width=w)  # N
    d.line([(cx0, cy0+r+gap), (cx0, n-pad)], fill=IC, width=w)  # S
    d.line([(pad, cy0),    (cx0-r-gap, cy0)], fill=IC, width=w)  # W
    d.line([(cx0+r+gap, cy0), (n-pad, cy0)], fill=IC, width=w)  # E

    # Centre dot
    rc = max(2, w-1)
    d.ellipse([cx0-rc, cy0-rc, cx0+rc, cy0+rc], fill=IC)

    _save(img, 24, 'C_24.png')


# ── letter badge ──────────────────────────────────────────────────────────────

def _badge(letter: str, px: int, name: str, colour=IC):
    img, n = _canvas(px)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)
    corner_r = max(4, n // 9)

    # Rounded-rect outline
    try:
        d.rounded_rectangle([pad, pad, n-pad, n-pad],
                             radius=corner_r, outline=colour, width=w)
    except TypeError:
        d.rectangle([pad, pad, n-pad, n-pad], outline=colour, width=w)

    # Centred letter
    try:
        font_size = int(n * 0.58)
        font = ImageFont.truetype(FONT, font_size)
        bb = d.textbbox((0, 0), letter, font=font)
        tw, th = bb[2]-bb[0], bb[3]-bb[1]
        tx = (n - tw) // 2 - bb[0]
        ty = (n - th) // 2 - bb[1]
        d.text((tx, ty), letter, font=font, fill=colour)
    except Exception as e:
        print(f'    font fallback for {letter!r}: {e}')
        d.text((n//2 - w*2, n//2 - w*3), letter, fill=colour)

    _save(img, px, name)


# ── geometric G badge ────────────────────────────────────────────────────────

def _G_badge(px: int, name: str, colour=IC):
    """Rounded-rect outline + thick C-arc + horizontal crossbar (= G)."""
    img, n = _canvas(px)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)
    corner_r = max(4, n // 9)

    # Rounded-rect outline
    try:
        d.rounded_rectangle([pad, pad, n-pad, n-pad],
                             radius=corner_r, outline=colour, width=w)
    except TypeError:
        d.rectangle([pad, pad, n-pad, n-pad], outline=colour, width=w)

    # G: thick arc (C-shape, opening on right) + crossbar
    inset = pad + max(3, int(n * 0.11))
    cx, cy = n // 2, n // 2
    sw = max(3, n // 9)                    # letter stroke width
    r  = cx - inset - sw // 2

    # Arc start=40→end=320 clockwise leaves a ~80° opening on the right side
    d.arc([cx-r, cy-r, cx+r, cy+r], start=40, end=320, fill=colour, width=sw)

    # Crossbar: from centre to right edge at mid-height
    bar_y = cy + int(r * 0.05)
    d.line([(cx, bar_y), (cx + r, bar_y)], fill=colour, width=sw)

    _save(img, px, name)


# ── geometric V badge ────────────────────────────────────────────────────────

def _V_badge(px: int, name: str, colour=IC):
    """Rounded-rect outline + two thick diagonal strokes meeting at a point."""
    img, n = _canvas(px)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)
    corner_r = max(4, n // 9)

    # Rounded-rect outline
    try:
        d.rounded_rectangle([pad, pad, n-pad, n-pad],
                             radius=corner_r, outline=colour, width=w)
    except TypeError:
        d.rectangle([pad, pad, n-pad, n-pad], outline=colour, width=w)

    # V: two thick lines meeting at the bottom centre
    inset = pad + max(3, int(n * 0.13))
    sw = max(3, n // 9)
    top_y = inset + sw // 2
    bot_y = n - inset
    cx    = n // 2
    d.line([(inset,     top_y), (cx, bot_y)], fill=colour, width=sw)
    d.line([(n - inset, top_y), (cx, bot_y)], fill=colour, width=sw)

    _save(img, px, name)


# ── save (floppy disk) ────────────────────────────────────────────────────────

def _save_icon():
    img, n = _canvas(24)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)

    notch_w = int(n * 0.25)     # top-right corner cut-in
    notch_h = int(n * 0.28)

    # Outer body (pentagon with notched corner)
    body = [
        (pad,              pad),
        (n - pad - notch_w, pad),
        (n - pad,          pad + notch_h),
        (n - pad,          n - pad),
        (pad,              n - pad),
    ]
    d.polygon(body, outline=IC, width=w)

    # Label window (upper inset rectangle)
    lx0 = pad + w
    lx1 = n - pad - notch_w - w
    ly0 = pad
    ly1 = pad + int(n * 0.32)
    d.rectangle([lx0, ly0, lx1, ly1], outline=IC, width=max(2, w // 2))

    # Write-protect notch (filled square on right inside the label window)
    tab = int(n * 0.13)
    d.rectangle([lx1 - tab, ly0 + w, lx1, ly1 - w], fill=IC)

    # Storage disk circle (lower centre)
    rd = max(w, int(n * 0.10))
    xc, yc = n // 2, int(n * 0.72)
    d.ellipse([xc-rd, yc-rd, xc+rd, yc+rd], outline=IC, width=max(2, w // 2))

    _save(img, 24, 'save_24.png')


# ── load (open folder) ────────────────────────────────────────────────────────

def _load_icon():
    img, n = _canvas(24)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)

    tab_h = int(n * 0.17)
    tab_w = int(n * 0.37)
    body_top = pad + tab_h

    # Folder tab (rounded bump, top-left)
    d.polygon([
        (pad,              body_top),
        (pad,              pad + w),
        (pad + tab_w,      pad + w),
        (pad + tab_w + tab_h, body_top),
    ], outline=IC, width=w)

    # Main folder body
    d.rectangle([pad, body_top, n - pad, n - pad], outline=IC, width=w)

    _save(img, 24, 'load_24.png')


# ── well (borehole) ───────────────────────────────────────────────────────────

def _well():
    img, n = _canvas(24)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)

    cx0 = n // 2
    # Drill-rig platform (small filled rectangle at top)
    pl_w = int(n * 0.30)
    pl_h = int(n * 0.16)
    d.rectangle([cx0 - pl_w, pad, cx0 + pl_w, pad + pl_h], fill=IC)

    # Borehole casing: two parallel vertical lines
    cas_x = int(n * 0.14)
    casing_top = pad + pl_h
    casing_bot = n - pad
    d.line([(cx0 - cas_x, casing_top), (cx0 - cas_x, casing_bot)], fill=IC, width=w)
    d.line([(cx0 + cas_x, casing_top), (cx0 + cas_x, casing_bot)], fill=IC, width=w)

    # Depth tick marks (right side)
    tick = int(n * 0.13)
    span = casing_bot - casing_top
    for i in range(1, 4):
        ty = casing_top + span * i // 4
        d.line([(cx0 + cas_x, ty), (cx0 + cas_x + tick, ty)], fill=IC, width=w)

    _save(img, 24, 'well_24.png')


# ── fault picker (lightning bolt) ────────────────────────────────────────────

def _faultline():
    """Lightning bolt — 8-point polygon with Z-kink at mid-height."""
    img, n = _canvas(24)
    d = _draw(img)
    pad = mg(n)

    # Bolt tilts upper-right → lower-left.  Two parallel bands share a Z-kink.
    # Points go clockwise:
    pts = [
        (int(n * 0.38), pad),               # top-left
        (int(n * 0.65), pad),               # top-right
        (int(n * 0.58), int(n * 0.50)),     # upper-right base
        (int(n * 0.74), int(n * 0.52)),     # notch — far right  (Z kink)
        (int(n * 0.60), n - pad),           # bottom-right
        (int(n * 0.35), n - pad),           # bottom-left
        (int(n * 0.44), int(n * 0.52)),     # lower-left apex
        (int(n * 0.28), int(n * 0.50)),     # notch — far left   (Z kink)
    ]
    d.polygon(pts, fill=IC)
    _save(img, 24, 'faultline_24.png')


# ── undo (curved left arrow) ──────────────────────────────────────────────────

def _undo():
    img, n = _canvas(24)
    d = _draw(img)
    w   = lw(n)

    # Centre and radius of the arc
    cx0, cy0 = int(n * 0.52), int(n * 0.46)
    r  = int(n * 0.28)

    # Draw top semicircle: PIL angles – 0=right, 90=bottom, 270=top
    # start=180 → end=360 traces left → top → right (the upper arc)
    d.arc([cx0-r, cy0-r, cx0+r, cy0+r], start=180, end=360, fill=IC, width=w)

    # Short descending stem at the RIGHT end of the arc (x=cx0+r, y=cy0)
    stem_len = int(r * 0.65)
    d.line([(cx0+r, cy0), (cx0+r, cy0+stem_len)], fill=IC, width=w)

    # Arrowhead at the LEFT end of the arc (x=cx0-r, y=cy0), pointing DOWN
    ah  = int(w * 2.2)    # arrowhead height
    ahw = int(w * 1.8)    # arrowhead half-width
    lx, ly = cx0 - r, cy0
    d.polygon([
        (lx,      ly + ah),        # tip (down)
        (lx - ahw, ly - ah // 2), # left
        (lx + ahw, ly - ah // 2), # right
    ], fill=IC)

    _save(img, 24, 'undo_24.png')


# ── pinch (two inward arrows + centre dot) ────────────────────────────────────

def _pinch():
    img, n = _canvas(24)
    d = _draw(img)
    w   = lw(n)
    pad = mg(n)
    cx0 = cy0 = n // 2
    arr = int(n * 0.16)     # arrowhead height
    ahw = int(n * 0.09)     # arrowhead half-width
    gap = int(n * 0.15)     # gap between arrowhead tip and centre

    # LEFT  arrow tip → pointing RIGHT toward centre
    tip_lx, tip_ly = cx0 - gap, cy0
    d.line([(pad + arr, cy0), (tip_lx, cy0)], fill=IC, width=w)
    d.polygon([
        (tip_lx,       tip_ly),
        (tip_lx - arr, tip_ly - ahw),
        (tip_lx - arr, tip_ly + ahw),
    ], fill=IC)

    # RIGHT arrow tip → pointing LEFT toward centre
    tip_rx, tip_ry = cx0 + gap, cy0
    d.line([(n - pad - arr, cy0), (tip_rx, cy0)], fill=IC, width=w)
    d.polygon([
        (tip_rx,       tip_ry),
        (tip_rx + arr, tip_ry - ahw),
        (tip_rx + arr, tip_ry + ahw),
    ], fill=IC)

    # Centre dot
    rc = max(2, w)
    d.ellipse([cx0-rc, cy0-rc, cx0+rc, cy0+rc], fill=IC)

    _save(img, 24, 'pinch_24.png')


# ── Python logo (two interlocked snake arcs) ─────────────────────────────────

def _python_16():
    """Two interlocked C-arcs + head dots — simplified Python snake logo."""
    img, n = _canvas(16)
    d = _draw(img)

    r   = int(n * 0.26)          # snake body arc radius
    sw  = max(5, n // 6)         # arc stroke width
    rh  = max(3, n // 12)        # head dot radius
    off = int(n * 0.09)          # centre offset for each snake's arc centre
    cx, cy = n // 2, n // 2

    # ── Upper snake ──────────────────────────────────────────────────────────
    # Arc centre: upper-left; 270° arc opens toward lower-right;
    # head dot placed at the upper-right endpoint (~310°).
    ux, uy = cx - off, cy - off
    d.arc([ux-r, uy-r, ux+r, uy+r], start=40, end=310, fill=IC, width=sw)
    ang1 = math.radians(310)
    hx1  = int(ux + r * math.cos(ang1))
    hy1  = int(uy + r * math.sin(ang1))
    d.ellipse([hx1-rh, hy1-rh, hx1+rh, hy1+rh], fill=IC)

    # ── Lower snake ──────────────────────────────────────────────────────────
    # Arc centre: lower-right (mirror of upper); 270° arc opens toward upper-left;
    # head dot placed at the lower-left endpoint (~130°).
    lx, ly = cx + off, cy + off
    d.arc([lx-r, ly-r, lx+r, ly+r], start=220, end=130, fill=IC, width=sw)
    ang2 = math.radians(130)
    hx2  = int(lx + r * math.cos(ang2))
    hy2  = int(ly + r * math.sin(ang2))
    d.ellipse([hx2-rh, hy2-rh, hx2+rh, hy2+rh], fill=IC)

    _save(img, 16, 'python_16.png')


# ── controls chevron (large_up_16) ────────────────────────────────────────────

def _chevron_up_16():
    img, n = _canvas(16)
    d = _draw(img)
    w   = lw(n)
    cx0 = n // 2
    pad = mg(n) + 2
    arm = int(n * 0.32)
    depth = int(n * 0.26)
    mid_y = n // 2 + depth // 2

    d.line([(cx0 - arm, mid_y + depth // 2),
            (cx0,        mid_y - depth)], fill=IC, width=w)
    d.line([(cx0,        mid_y - depth),
            (cx0 + arm,  mid_y + depth // 2)], fill=IC, width=w)

    _save(img, 16, 'large_up_16.png')


# ── generate all icons ────────────────────────────────────────────────────────

def main():
    print('Generating icons...')

    # File operations
    _save_icon()
    _load_icon()

    # Letter badges  (24 px)
    _G_badge(24, 'G_24.png')
    _V_badge(24, 'V_24.png')
    _badge('M', 24, 'M_24.png')
    _badge('F', 24, 'F_24.png')
    _badge('F', 24, 'off_F_24.png', colour=DIM)

    # Letter badges (16 px)
    _G_badge(16, 'G_16.png')
    _V_badge(16, 'V_16.png')
    _badge('M', 16, 'M_16.png')
    _badge('T', 16, 'T_16.png')

    # Aspect arrows – large (thick) & small (thin), up / down
    _arrow('up',    True,  24, 'large_up_24.png')
    _arrow('down',  True,  24, 'large_down_24.png')
    _arrow('up',    False, 24, 'small_up_24.png')
    _arrow('down',  False, 24, 'small_down_24.png')

    # Gain / transparency arrows – horizontal
    _arrow('left',  False, 24, 'left_small_24.png')
    _arrow('right', False, 24, 'right_small_24.png')
    _arrow('left',  True,  24, 'large_left_24.png')
    _arrow('right', True,  24, 'large_right_24.png')

    # Navigation
    _magnifier('+', 'zoom_in_24.png')
    _magnifier('-', 'zoom_out_24.png')
    _full_extent()
    _pan()

    # Domain tools
    _crosshair()
    _well()
    _faultline()
    _undo()
    _pinch()

    # Status-bar (16 px)
    _chevron_up_16()
    _python_16()

    print('Done.')


if __name__ == '__main__':
    main()
