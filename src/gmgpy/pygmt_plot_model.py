"""
Save the current GMG model as an external figure using PyGMT.

Called from gmg.plot_model().
"""

import numpy as np
import pygmt
import matplotlib.colors as mcolors


def _to_gmt_color(color):
    """Convert any matplotlib-compatible color to a GMT 'R/G/B' string."""
    r, g, b, *_ = mcolors.to_rgba(color)
    return f"{int(r * 255)}/{int(g * 255)}/{int(b * 255)}"


def plot_fig(file_path, area, xp, model_aspect,
             layer_list,
             topography_frame_visible, observed_topography_list, topo_ylim,
             gravity_frame_visible, observed_gravity_list, predicted_gravity, gravity_ylim,
             vgg_frame_visible, observed_vgg_list, predicted_vgg, vgg_ylim,
             magnetic_frame_visible, observed_magnetic_list, predicted_nt, mag_ylim,
             model_xlim, model_ylim):
    """
    Build a multi-panel PyGMT figure of the current GMG model and save to file.

    Panel layout (bottom to top):
        Model | Gravity | VGG | Magnetics | Topography
    Panels for frames that are not visible in the GUI are omitted.

    Parameters
    ----------
    file_path : str
        Output file path. The extension determines format (.pdf, .png, .eps, .svg).
    area : array-like
        Full model area in metres [x_min, x_max, z_min, z_max].
    xp : array-like
        Calculation x-positions in metres.
    model_aspect : float
        Vertical exaggeration (>1 = expanded vertically).
    layer_list : list of Layer
        All layer objects; index 0 is the background layer and is skipped.
    topography_frame_visible : bool
    observed_topography_list : list of ObservedData or None
    topo_ylim : (float, float)
    gravity_frame_visible : bool
    observed_gravity_list : list of ObservedData or None
    predicted_gravity : array or None
    gravity_ylim : (float, float)
    vgg_frame_visible : bool
    observed_vgg_list : list of ObservedData or None
    predicted_vgg : array or None
    vgg_ylim : (float, float)
    magnetic_frame_visible : bool
    observed_magnetic_list : list of ObservedData or None
    predicted_nt : array or None
    mag_ylim : (float, float)
    model_xlim : (float, float)   X-axis limits in km (from the matplotlib model axes).
    model_ylim : (float, float)   Z-axis limits in km (from the matplotlib model axes).
    """

    fig = pygmt.Figure()

    # ------------------------------------------------------------------
    # FIGURE DIMENSIONS
    # ------------------------------------------------------------------
    fig_width_cm = 15.0      # horizontal width shared by all panels
    data_panel_h_cm = 4.0    # fixed height for each data panel

    x_min, x_max = model_xlim

    # matplotlib returns ylim inverted for depth axes, e.g. (10.0, 0.0).
    # GMT requires y_lo < y_hi in -R; we flip via a negative projection height.
    z_lo = min(model_ylim)   # shallowest value (top of figure)
    z_hi = max(model_ylim)   # deepest value   (bottom of figure)

    x_range = x_max - x_min
    z_range = z_hi - z_lo

    # Model panel height respects the current vertical exaggeration.
    # A larger model_aspect means more vertical compression in the output.
    raw_h = (z_range / x_range) * fig_width_cm
    model_h_cm = max(raw_h / model_aspect, 4.0)

    def _proj(height_cm, invert_y=False):
        """Absolute Cartesian projection for a panel of the given height.
        Pass invert_y=True to flip the y-axis (depth-down convention).
        """
        sign = "-" if invert_y else ""
        return f"X{fig_width_cm}c/{sign}{height_cm}c"

    xp_km = np.asarray(xp, dtype=float) * 1e-3

    # ------------------------------------------------------------------
    # HELPER: draw one data panel at the current origin
    # ------------------------------------------------------------------
    def _draw_data_panel(ylim, y_label, obs_list, pred_array, pred_color):
        y_lo, y_hi = min(ylim), max(ylim)
        fig.basemap(
            region=[x_min, x_max, y_lo, y_hi],
            projection=_proj(data_panel_h_cm),
            frame=["WSen", f"yaf+l{y_label}", "xaf"],
        )
        # Observed data points
        if obs_list:
            for obs in obs_list:
                if obs is None:
                    continue
                if hasattr(obs, 'mpl_actor') and obs.mpl_actor is not None:
                    if not obs.mpl_actor.get_visible():
                        continue
                color_str = _to_gmt_color(obs.color if obs.color is not None else "black")
                fig.plot(
                    x=obs.data[:, 0],
                    y=obs.data[:, 1],
                    style="c0.08c",
                    fill=color_str,
                    pen="0.25p,black",
                )
        # Predicted curve
        if pred_array is not None:
            fig.plot(
                x=xp_km,
                y=np.asarray(pred_array, dtype=float),
                pen=f"1p,{pred_color}",
            )

    # ------------------------------------------------------------------
    # DECIDE WHICH DATA PANELS TO DRAW (bottom → top order above model)
    # ------------------------------------------------------------------
    data_panels = []
    if gravity_frame_visible:
        data_panels.append((gravity_ylim, "Gravity (mGal)",
                             observed_gravity_list, predicted_gravity, "red"))
    if vgg_frame_visible:
        data_panels.append((vgg_ylim, "VGG (Eotvos)",
                             observed_vgg_list, predicted_vgg, "blue"))
    if magnetic_frame_visible:
        data_panels.append((mag_ylim, "Magnetics (nT)",
                             observed_magnetic_list, predicted_nt, "green"))
    if topography_frame_visible:
        data_panels.append((topo_ylim, "Elevation (km)",
                             observed_topography_list, None, "black"))

    # ------------------------------------------------------------------
    # MODEL PANEL  (drawn first; sits at the bottom of the stacked figure)
    # ------------------------------------------------------------------
    # Use a negative projection height so depth increases downward in the figure.
    fig.basemap(
        region=[x_min, x_max, z_lo, z_hi],
        projection=_proj(model_h_cm, invert_y=True),
        frame=["WSen", "xaf+lDistance (km)", "yaf+lDepth (km)"],
    )

    # Layer polygons (skip index 0 = background layer)
    for layer in layer_list[1:]:
        if not layer.polygon_mpl_actor:
            continue
        actor = layer.polygon_mpl_actor[0]
        if not actor.get_visible():
            continue

        xy = actor.get_xy()
        x_poly = xy[:, 0]
        y_poly = xy[:, 1]

        fc = actor.get_facecolor()
        rgba = fc[0] if (hasattr(fc, 'ndim') and fc.ndim == 2) else fc
        alpha_pct = float(actor.get_alpha() or 1.0) * 100
        fill_str = _to_gmt_color(rgba)
        pen_color = _to_gmt_color(layer.color)

        fig.plot(
            x=x_poly,
            y=y_poly,
            fill=fill_str,
            transparency=alpha_pct,
            pen=f"0.5p,{pen_color}",
        )

    # ------------------------------------------------------------------
    # DATA PANELS  (each shifted upward from the previous panel)
    # ------------------------------------------------------------------
    prev_h = model_h_cm
    for ylim, ylabel, obs_list, pred, pred_color in data_panels:
        fig.shift_origin(yshift=f"{prev_h}c")
        _draw_data_panel(ylim, ylabel, obs_list, pred, pred_color)
        prev_h = data_panel_h_cm

    # ------------------------------------------------------------------
    # SAVE FIGURE
    # ------------------------------------------------------------------
    fig.savefig(file_path)
