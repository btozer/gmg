"""
theme.py — VS Code-inspired colour palette and style helpers for gmgpy.

All GUI modules should import colours and fonts from here so the entire
theme can be changed in one place.

Usage
-----
    from gmgpy.theme import THEME, get_font
    panel.SetBackgroundColour(THEME['bg_editor'])
    label.SetFont(get_font())
"""

import wx

# ---------------------------------------------------------------------------
# Colour palette (mirrors VS Code Dark+ defaults)
# ---------------------------------------------------------------------------

THEME = {
    # ---- Backgrounds -------------------------------------------------------
    "bg_base":          "#1e1e1e",  # Main frame / editor background
    "bg_panel":         "#252526",  # Side panel / activity bar
    "bg_sidebar":       "#2d2d2d",  # Scrolled windows, tree panels
    "bg_input":         "#3c3c3c",  # TextCtrl / SpinCtrl inputs
    "bg_toolbar":       "#333333",  # Toolbar strip
    "bg_statusbar":     "#007acc",  # Status bar (VS Code blue)
    "bg_statusbar_err": "#c72e0f",  # Status bar error state
    "bg_tooltip":       "#252526",  # Tooltip background
    "bg_canvas":        "#1e1e1e",  # Matplotlib figure background
    "bg_axes":          "#1e1e1e",  # Matplotlib axes face

    # ---- Foregrounds -------------------------------------------------------
    "fg_primary":       "#cccccc",  # Normal text
    "fg_secondary":     "#858585",  # Dimmed / secondary text
    "fg_statusbar":     "#ffffff",  # Status bar text
    "fg_caption":       "#bbbbbb",  # AUI pane captions
    "fg_input":         "#d4d4d4",  # Input widget text
    "fg_placeholder":   "#6a6a6a",  # Placeholder / hint text

    # ---- Accents -----------------------------------------------------------
    "accent":           "#007acc",  # Primary accent (VS Code blue)
    "accent_hover":     "#1a8ad4",  # Hover highlight
    "accent_active":    "#005f9e",  # Active / pressed

    # ---- Borders & separators ---------------------------------------------
    "border":           "#3c3c3c",  # Panel/pane borders
    "border_focus":     "#007acc",  # Focused input border
    "separator":        "#454545",  # Thin divider lines

    # ---- Matplotlib plot colours ------------------------------------------
    "plot_grid":        "#3c3c3c",  # Axes grid lines
    "plot_spine":       "#555555",  # Axes spine lines
    "plot_tick":        "#858585",  # Tick labels
    "plot_label":       "#cccccc",  # Axis labels
    "plot_title":       "#cccccc",  # Plot title text
    "plot_pred_grav":   "#e06c75",  # Predicted gravity line   (red)
    "plot_pred_vgg":    "#e5c07b",  # Predicted VGG line       (yellow)
    "plot_pred_mag":    "#98c379",  # Predicted magnetic line  (green)
    "plot_rms":         "#c678dd",  # RMS residual line        (purple)
}


def wx_colour(key: str) -> wx.Colour:
    """Return a wx.Colour for the given palette key."""
    return wx.Colour(THEME[key])


# ---------------------------------------------------------------------------
# Typography
# ---------------------------------------------------------------------------

def get_font(size: int = 10, bold: bool = False) -> wx.Font:
    """
    Return the UI font.  Uses the native macOS system font when available,
    falling back to the wx default UI font on other platforms.
    """
    weight = wx.FONTWEIGHT_BOLD if bold else wx.FONTWEIGHT_NORMAL
    font = wx.Font(
        size,
        wx.FONTFAMILY_DEFAULT,
        wx.FONTSTYLE_NORMAL,
        weight,
        faceName="SF Pro Text",   # macOS system font
    )
    if not font.IsOk():
        # Fallback: let wx choose the platform default
        font = wx.Font(
            size,
            wx.FONTFAMILY_DEFAULT,
            wx.FONTSTYLE_NORMAL,
            weight,
        )
    return font


def get_mono_font(size: int = 10) -> wx.Font:
    """Return a monospace font for console/code areas."""
    font = wx.Font(
        size,
        wx.FONTFAMILY_TELETYPE,
        wx.FONTSTYLE_NORMAL,
        wx.FONTWEIGHT_NORMAL,
        faceName="Menlo",   # macOS monospace
    )
    if not font.IsOk():
        font = wx.Font(
            size,
            wx.FONTFAMILY_TELETYPE,
            wx.FONTSTYLE_NORMAL,
            wx.FONTWEIGHT_NORMAL,
        )
    return font


# ---------------------------------------------------------------------------
# Matplotlib rcParams dict  (apply with matplotlib.rcParams.update())
# ---------------------------------------------------------------------------

MPL_DARK_RC = {
    "figure.facecolor":     THEME["bg_canvas"],
    "axes.facecolor":       THEME["bg_axes"],
    "axes.edgecolor":       THEME["plot_spine"],
    "axes.labelcolor":      THEME["plot_label"],
    "axes.grid":            True,
    "grid.color":           THEME["plot_grid"],
    "grid.linewidth":       0.5,
    "xtick.color":          THEME["plot_tick"],
    "ytick.color":          THEME["plot_tick"],
    "xtick.labelcolor":     THEME["plot_tick"],
    "ytick.labelcolor":     THEME["plot_tick"],
    "text.color":           THEME["fg_primary"],
    "legend.facecolor":     THEME["bg_panel"],
    "legend.edgecolor":     THEME["border"],
    "legend.labelcolor":    THEME["fg_primary"],
    "savefig.facecolor":    THEME["bg_canvas"],
    "savefig.edgecolor":    THEME["bg_canvas"],
}


# ---------------------------------------------------------------------------
# Icon loader  (dark-theme adapted bitmaps for the toolbar)
# ---------------------------------------------------------------------------

def load_icon(path: str) -> "wx.Bitmap":
    """
    Load a 24-px PNG toolbar icon adapted for a dark background:
      - near-white pixels  (R,G,B > 220)  → fully transparent
      - near-black pixels  (unsaturated, max channel < 80) → light grey #cccccc
      - all other (coloured) pixels are left untouched
    """
    img = wx.Image(path, wx.BITMAP_TYPE_PNG)
    if not img.HasAlpha():
        img.InitAlpha()
    for y in range(img.GetHeight()):
        for x in range(img.GetWidth()):
            r = img.GetRed(x, y)
            g = img.GetGreen(x, y)
            b = img.GetBlue(x, y)
            if r > 220 and g > 220 and b > 220:          # near-white → transparent
                img.SetAlpha(x, y, 0)
            elif max(r, g, b) < 80 and (max(r, g, b) - min(r, g, b)) < 30:
                img.SetRed(x, y, 204)                     # near-black → light grey
                img.SetGreen(x, y, 204)
                img.SetBlue(x, y, 204)
    return wx.Bitmap(img)
