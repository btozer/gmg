"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI application for Forward modelling 2D potential field profiles.
Written by Brook Tozer, University of Oxford 2015-17. SIO 2018-19.
Includes ability to import seismic reflection, well, surface outcrop and xy points into the model frame.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Dependencies**

SciPy
NumPy
Matplotlib
pylab
pickle
fatiando
obspy
wxpython
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**References**

***
polygons from Fatiando a Terra.

Uieda, L., V. C. Oliveira Jr and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings
of the 12th Python in Science Conference

www.fatiando.org/
***

***
Gravity algorithm written using NumPy by Brook Tozer (2015) (Modified from the Fatiando a Terra 2D gravity code).

CODE MODIFIED FROM: bott, M. H. P. (1969). GRAVN. Durham geophysical computer specification No. 1.
***

***
Magnetic algorithm written using NumPy by Brook Tozer (2015)

CODE MODIFIED FROM: Talwani, M., & Heirtzler, J. R. (1964). Computation of magnetic anomalies caused by two dimensional
structures of arbitrary shape, in Parks, G. A., Ed., Computers in the mineral industries, Part 1: Stanford Univ. Publ.,
Geological Sciences, 9, 464-480.
***

***
SEGY plotting is preformed using ObsPy.

ObsPy; a python toolbox for seismology Seismological Research Letters(May 2010), 81(3):530-533

obspy.org
***

***
Icons where designed using the Free icon Maker:

icons8.com
https://icons8.com/icon/set/circle/all

Pixel size: 24
Font: Roboto Slab
Font size: 200
Color: 3498db
***
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Documentation created using Sphinx.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NB. before launching gmg you must run:

    source activate py27-gmg

"""

# IMPORT MODULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import wx
import matplotlib

matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.colors as colors
from wx.lib.agw import floatspin as fs
import wx.grid as gridlib
import wx.lib.agw.customtreectrl as ct
import wx.py as py
import wx.lib.agw.aui as aui
from wx.lib.buttons import GenBitmapButton
import wx.lib.agw.foldpanelbar as fpb
import pylab as plt
import numpy as np
import csv
import math as m
import os
import sys
from sys import platform
from obspy import read
import cPickle as Pickle
from scipy import signal
from scipy import interpolate as ip
from fatiando.mesher import Polygon
import plot_model
import bott
import talwani_and_heirtzler
import model_stats
import struct
import gc
import webbrowser
import time
# FUTURE
# import wx.lib.agw.ribbon as RB
# import wx.EnhancedStatusBar as ESB
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class Gmg(wx.Frame):
    """
    Master Class for GMG GUI.
    Most functions are contained in this Class.
    Upon startup this sets the panels, sizer's and event bindings.
    Additional classes are used for handling "pop out" windows (Dialog boxes).
    Objects are passed between this Master Class and the Dialog Boxes as I/O.
    """

    def __init__(self, *args, **kwds):
        wx.Frame.__init__(self, None, wx.ID_ANY, 'gmg: 2D Geophysical Modelling GUI', size=(1800, 1050))

        # DIR CONTAINING PROGRAM ICONS
        self.gui_icons_dir = os.path.dirname(os.path.abspath(__file__))+"/icons/"

        # START AUI WINDOW MANAGER
        self.mgr = aui.AuiManager()

        # TELL AUI WHICH FRAME TO USE
        self.mgr.SetManagedWindow(self)

        # SET AUI ICON SIZING AND STYLES (ARROWS)
        images = wx.ImageList(16, 16)
        top = wx.ArtProvider.GetBitmap(wx.ART_GO_UP, wx.ART_MENU, (16, 16))
        bottom = wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN, wx.ART_MENU, (16, 16))
        images.Add(top)
        images.Add(bottom)

        # CREATE PANELS TO FILL WITH ATTRIBUTE CONTROLS, LAYER TREE CONTROL AND FAULT TREE CONTROL
        self.leftPanel = wx.SplitterWindow(self, wx.ID_ANY, size=(200, 1100), style=wx.SP_NOBORDER | wx.EXPAND)
        self.leftPanel_b = wx.SplitterWindow(self.leftPanel, wx.ID_ANY, size=(200, 705),
                                             style=wx.SP_NOBORDER | wx.EXPAND)
        self.leftPanel.SetMinimumPaneSize(1)
        self.leftPanel_b.SetMinimumPaneSize(1)

        # FIRST PANE; LEFT PANEL (=ATTRIBUTES)
        self.splitter_left_panel_one = wx.ScrolledWindow(self.leftPanel, wx.ID_ANY, size=(-1, -1),
                                                         style=wx.ALIGN_LEFT | wx.BORDER_RAISED | wx.EXPAND)
        self.controls_panel_bar_one = fpb.FoldPanelBar(self.splitter_left_panel_one, 1, size=(200, 395),
                                                       agwStyle=fpb.FPB_VERTICAL)
        self.fold_panel_one = self.controls_panel_bar_one.AddFoldPanel("Layer Attributes", collapsed=True,
                                                                       foldIcons=images)
        self.controls_panel_bar_one.Expand(self.fold_panel_one)  # ENSURES FOLD PANEL IS VISIBLE

        # SECOND PANE; LEFT PANEL (=LAYERS)
        self.splitter_left_panel_two = wx.ScrolledWindow(self.leftPanel_b, wx.ID_ANY, size=(-1, -1),
                                                         style=wx.ALIGN_LEFT | wx.BORDER_RAISED | wx.EXPAND)
        self.controls_panel_bar_two = fpb.FoldPanelBar(self.splitter_left_panel_two, 1, size=(200, 455),
                                                       agwStyle=fpb.FPB_VERTICAL)
        self.fold_panel_two = self.controls_panel_bar_two.AddFoldPanel("Layers", collapsed=True, foldIcons=images)
        self.controls_panel_bar_two.Expand(self.fold_panel_two)  # ENSURES FOLD PANEL IS VISIBLE

        # THIRD PANE; LEFT PANEL (=FAULTS)
        self.splitter_left_panel_three = wx.ScrolledWindow(self.leftPanel_b, wx.ID_ANY, size=(-1, -1),
                                                           style=wx.ALIGN_LEFT | wx.BORDER_RAISED | wx.EXPAND)
        self.controls_panel_bar_three = fpb.FoldPanelBar(self.splitter_left_panel_three, 1, size=(200, 250),
                                                         agwStyle=fpb.FPB_VERTICAL)
        self.fold_panel_three = self.controls_panel_bar_three.AddFoldPanel("Faults", collapsed=True, foldIcons=images)
        self.controls_panel_bar_three.Expand(self.fold_panel_three)  # ENSURES FOLD PANEL IS VISIBLE

        # SET SPLITTERS
        self.leftPanel_b.SplitHorizontally(self.splitter_left_panel_two, self.splitter_left_panel_three)
        self.leftPanel.SplitHorizontally(self.splitter_left_panel_one, self.leftPanel_b)

        self.splitter_left_panel_sizer = wx.BoxSizer(wx.VERTICAL)
        self.splitter_left_panel_sizer.Add(self.leftPanel, 1, wx.EXPAND)
        self.splitter_left_panel_one.SetScrollbar(1, 1, 10, 10)
        self.splitter_left_panel_two.SetScrollbar(1, 1, 10, 10)
        self.splitter_left_panel_three.SetScrollbar(1, 1, 10, 10)

        # CREATE PANEL TO FILL WITH MATPLOTLIB INTERACTIVE FIGURE (MAIN GUI MODELLING FRAME)
        self.rightPanel = wx.Panel(self, -1, size=(1700, 1100), style=wx.ALIGN_RIGHT | wx.BORDER_RAISED | wx.EXPAND)

        # CREATE PANEL FOR PYTHON CONSOLE (USED FOR DEBUGGING AND CUSTOM USAGES)
        self.ConsolePanel = wx.Panel(self, -1, size=(1700, 100), style=wx.ALIGN_LEFT | wx.BORDER_RAISED | wx.EXPAND)
        intro = "###############################################################\r" \
                "!USE import sys; then sys.Gmg.OBJECT TO ACCESS PROGRAM OBJECTS \r" \
                "ctrl+up FOR COMMAND HISTORY                                    \r" \
                "###############################################################"
        py_local = {'__app__': 'gmg Application'}
        sys.gmg = self
        self.win = py.shell.Shell(self.ConsolePanel, -1, size=(2200, 1100), locals=py_local, introText=intro)

        # ADD THE PANES TO THE AUI MANAGER
        self.mgr.AddPane(self.leftPanel, aui.AuiPaneInfo().Name('left').Left().Caption("Controls"))
        self.mgr.AddPane(self.rightPanel, aui.AuiPaneInfo().Name('right').CenterPane())
        self.mgr.AddPane(self.ConsolePanel, aui.AuiPaneInfo().Name('console').Bottom().Caption("Console"))
        self.mgr.GetPaneByName('console').Hide()  # HIDE PYTHON CONSOLE BY DEFAULT
        self.mgr.Update()

        # CREATE PROGRAM MENUBAR & TOOLBAR (PLACED AT TOP OF FRAME)
        self.create_menu()

        # CREATE STATUS BAR
        self.statusbar = self.CreateStatusBar(3)
        self.controls_button = GenBitmapButton(self.statusbar, -1, wx.Bitmap(self.gui_icons_dir + 'large_up_16.png'),
                                               pos=(0, -4), style=wx.NO_BORDER)
        self.Bind(wx.EVT_BUTTON, self.show_controls, self.controls_button)

        # PYTHON CONSOLE
        self.console_button = GenBitmapButton(self.statusbar, -1, wx.Bitmap(self.gui_icons_dir + 'python_16.png'),
                                              pos=(24, -4), style=wx.NO_BORDER)
        self.Bind(wx.EVT_BUTTON, self.show_console, self.console_button)

        # TOPOGRAPHY TOGGLE BUTTON
        self.topography_button = GenBitmapButton(self.statusbar, 601, wx.Bitmap(self.gui_icons_dir + 'T_16.png'),
                                                 pos=(48, -4), style=wx.NO_BORDER)
        self.Bind(wx.EVT_BUTTON, self.frame_adjustment, self.topography_button)

        # GRAVITY TOGGLE BUTTON
        self.gravity_button = GenBitmapButton(self.statusbar, 602, wx.Bitmap(self.gui_icons_dir + 'G_16.png'),
                                              pos=(72, -4), style=wx.NO_BORDER)
        self.Bind(wx.EVT_BUTTON, self.frame_adjustment, self.gravity_button)

        # MAGNETIC TOGGLE BUTTON
        self.magnetic_button = GenBitmapButton(self.statusbar, 603, wx.Bitmap(self.gui_icons_dir + 'M_16.png'),
                                               pos=(96, -4), style=wx.NO_BORDER)
        self.Bind(wx.EVT_BUTTON, self.frame_adjustment, self.magnetic_button)

        self.status_text = " "
        self.statusbar.SetStatusWidths([-1, -1, 1700])
        self.statusbar.SetStatusText(self.status_text, 2)
        self.statusbar.SetSize((1800, 24))

        # SET PROGRAM STATUS
        self.model_saved = False
        self.newmodel = False

        # BIND PROGRAM EXIT BUTTON WITH EXIT FUNCTION
        self.Bind(wx.EVT_CLOSE, self.on_close_button)

        # MAXIMIZE FRAME
        self.Maximize(True)

    def create_menu(self):
        """CREATE GUI MENUBAR"""
        self.menubar = wx.MenuBar()
        menu_file = wx.Menu()
        m_new_model = menu_file.Append(-1, "New Model...\tCtrl-N", "New Model...")
        self.Bind(wx.EVT_MENU, self.new_model, m_new_model)
        m_load_model = menu_file.Append(-1, "Load Model...\tCtrl-L", "Load Model...")
        self.Bind(wx.EVT_MENU, self.load_model, m_load_model)
        m_save_model = menu_file.Append(-1, "Save Model...\tCtrl-S", "Save Model...")
        self.Bind(wx.EVT_MENU, self.save_model, m_save_model)
        menu_file.AppendSeparator()
        m_save_xy = menu_file.Append(-1, "Save Layers As ASCII .xy File ...\tCtrl-Shift-L",
                                     "Save Layers as ASCII .xy...")
        self.Bind(wx.EVT_MENU, self.write_layers_xy, m_save_xy)
        m_save_c = menu_file.Append(-1, "Save Model As RayInvr c.in File...\tCtrl-Shift-C", "Save RAYINVR c.in File...")
        self.Bind(wx.EVT_MENU, self.write_c_xy, m_save_c)
        m_save_fig = menu_file.Append(-1, "Save Figure...\tCtrl-Shift-F", "Save model Figure...")
        self.Bind(wx.EVT_MENU, self.plot_model, m_save_fig)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "Exit...\tCtrl-X", "Exit...")
        self.Bind(wx.EVT_MENU, self.exit, m_exit)

        # DRAW MENU
        self.menubar.Append(menu_file, "&File")

        # MODEL VIEW MENU
        model_view_file = wx.Menu()
        m_modify_model_dimensions = model_view_file.Append(-1, "Modify Current Model Dimensions...\tCtrl-M",
                                                           "Modify Current Model Dimensions...")
        self.Bind(wx.EVT_MENU, self.modify_model_dimensions, m_modify_model_dimensions)
        model_view_file.AppendSeparator()

        # PROGRAM FRAME WINDOW SWITCHES
        self.topo_frame_switch = True
        self.gravity_frame_switch = True
        self.magnetic_frame_switch = True
        self.m_model_frames_submenu = wx.Menu()
        model_view_file.Append(-1, 'Toggle Model Frames', self.m_model_frames_submenu)
        self.m_model_frames_submenu.Append(601, 'Topography')
        self.Bind(wx.EVT_MENU, self.frame_adjustment, id=601)
        self.m_model_frames_submenu.Append(602, 'Gravity')
        self.Bind(wx.EVT_MENU, self.frame_adjustment, id=602)
        self.m_model_frames_submenu.Append(603, 'Magnetics')
        self.Bind(wx.EVT_MENU, self.frame_adjustment, id=603)
        model_view_file.AppendSeparator()

        m_aspect_increase = model_view_file.Append(-1, "&Increase Aspect Ratio...\tCtrl-shift-up",
                                                   "Increase the aspect ratio of the model plot...")
        self.Bind(wx.EVT_MENU, self.aspect_increase, m_aspect_increase)
        m_aspect_decrease = model_view_file.Append(-1, "&Decrease Aspect Ratio...\tCtrl-shift-down",
                                                   "Decrease the aspect ratio of the model plot...")
        self.Bind(wx.EVT_MENU, self.aspect_decrease, m_aspect_decrease)
        self.menubar.Append(model_view_file, "&Model View")

        # GRAVITY DATA MENU --------------------------------------------------------------------------------------------
        self.gravity_data = wx.Menu()
        # LOAD OBSERVED GRAVITY DATA
        m_load_obs_g = self.gravity_data.Append(-1, "&Load Gravity Anomaly...\tCtrl-L", "Load Observed Gravity Data...")
        self.Bind(wx.EVT_MENU, self.load_obs_g, m_load_obs_g)
        # EDIT
        self.m_obs_g_submenu = wx.Menu()
        self.gravity_data.Append(-1, 'Gravity Data...', self.m_obs_g_submenu)
        # FILTER MENU
        grav_m_filter_observed = self.gravity_data.Append(-1, "Median Filter...", "Filter Observed Anomaly")
        self.Bind(wx.EVT_MENU, self.filter_observed_gravity, grav_m_filter_observed)
        # HORIZONTAL DERIVATIVE
        grav_m_horizontal_derivative = self.gravity_data.Append(-1, "Take Horizontal Derivative...",
                                                                "Take Horizontal Derivative")
        self.Bind(wx.EVT_MENU, self.take_gravity_horizontal_derivative, grav_m_horizontal_derivative)
        # SET RMS OBS ARRAYS
        grav_m_set_rms = self.gravity_data.Append(-1, "Set RMS input...", "Set RMS input...")
        self.Bind(wx.EVT_MENU, self.set_obs_grav_rms, grav_m_set_rms)
        # SET BACKGROUND DENSITY
        m_set_background_density = self.gravity_data.Append(-1, "&Set Background Density...",
                                                            "Set Background Density...")
        self.Bind(wx.EVT_MENU, self.set_background_density, m_set_background_density)
        # SAVE PREDICTED ANOMALY TO DISC
        m_save_g_submenu = self.gravity_data.Append(-1, "&Save Predicted Anomaly...",
                                                    "Save Predicted Anomaly to Disc...")
        self.Bind(wx.EVT_MENU, self.save_modelled_grav, m_save_g_submenu)
        # DRAW MENU
        self.menubar.Append(self.gravity_data, "&Gravity")
        # --------------------------------------------------------------------------------------------------------------

        # MAGNETIC DATA MENU -------------------------------------------------------------------------------------------
        self.magnetic_data = wx.Menu()
        # LOAD OBSERVED MAGNETIC DATA
        m_load_obs_m = self.magnetic_data.Append(-1, "&Load Magnetic Anomaly...\tCtrl-L",
                                                 "Load Observed Magnetic Data...")
        self.Bind(wx.EVT_MENU, self.load_obs_m, m_load_obs_m)
        # EDIT
        self.m_obs_mag_submenu = wx.Menu()
        self.magnetic_data.Append(-1, 'Magnetic Anomalies...', self.m_obs_mag_submenu)
        # FILTER MENU
        mag_m_filter_observed = self.magnetic_data.Append(-1, "Median Filter...", "Filter Observed Anomaly")
        self.Bind(wx.EVT_MENU, self.filter_observed_magnetic, mag_m_filter_observed)
        # HORIZONTAL DERIVATIVE
        mag_m_horizontal_derivative = self.magnetic_data.Append(-1, "Take Horizontal Derivative...",
                                                                "Take Horizontal Derivative")
        self.Bind(wx.EVT_MENU, self.take_magnetic_horizontal_derivative, mag_m_horizontal_derivative)
        # SET RMS OBS ARRAYS
        mag_m_set_rms = self.magnetic_data.Append(-1, "Set RMS input..\tCtrl-L", "Set RMS input..")
        self.Bind(wx.EVT_MENU, self.set_obs_mag_rms, mag_m_set_rms)
        # SET MAG
        m_set_mag_variables = self.magnetic_data.Append(-1, "&Set Magnetic Field...\tCtrl-shift-up",
                                                        "Set Magnetic Feild...")
        self.Bind(wx.EVT_MENU, self.set_mag_variables, m_set_mag_variables)
        # SAVE PREDICTED ANOMALY TO DISC
        m_save_mag_submenu = self.magnetic_data.Append(-1, "&Save Predicted Anomaly...\tCtrl-shift-S",
                                                       "Save Predicted Anomaly to Disc...")
        self.Bind(wx.EVT_MENU, self.save_modelled_mag, m_save_mag_submenu)
        # DRAW MENU
        self.menubar.Append(self.magnetic_data, "&Magnetics")
        # --------------------------------------------------------------------------------------------------------------

        # TOPOGRAPHY DATA MENU -----------------------------------------------------------------------------------------
        self.topography_data = wx.Menu()
        # LOAD OBSERVED TOPO DATA
        m_load_topo = self.topography_data.Append(-1, "&Load Topography...\tCtrl-L", "Load Observed Topography...")
        self.Bind(wx.EVT_MENU, self.load_topo, m_load_topo)
        # EDIT
        self.m_topo_submenu = wx.Menu()
        self.topography_data.Append(-1, 'Topography Data', self.m_topo_submenu)
        # FILTER MENU
        topo_m_filter_observed = self.topography_data.Append(-1, "Median Filter...", "Filter Observed Anomaly")
        self.Bind(wx.EVT_MENU, self.filter_observed_topography, topo_m_filter_observed)
        # HORIZONTAL DERIVATIVE
        topo_m_horizontal_derivative = self.topography_data.Append(-1, "Take Horizontal Derivative...",
                                                                "Take Horizontal Derivative")
        self.Bind(wx.EVT_MENU, self.take_topography_horizontal_derivative, topo_m_horizontal_derivative)
        # DRAW MENU
        self.menubar.Append(self.topography_data, "&Topography")
        # --------------------------------------------------------------------------------------------------------------

        # XY DATA MENU -------------------------------------------------------------------------------------------------
        self.xy_data = wx.Menu()
        # LOAD XY DATA
        m_load_xy = self.xy_data.Append(-1, "&Load XY Points...\tCtrl-L", "Load XY Points...")
        self.Bind(wx.EVT_MENU, self.load_xy, m_load_xy)
        self.m_xy_submenu = wx.Menu()
        self.xy_data.Append(-1, 'XY Data...', self.m_xy_submenu)
        # DRAW MENU
        self.menubar.Append(self.xy_data, "&XY Data")
        # --------------------------------------------------------------------------------------------------------------

        # SEISMIC DATA -------------------------------------------------------------------------------------------------
        self.seismic_data = wx.Menu()
        # SEGY LOAD
        self.m_load_segy = self.seismic_data.Append(-1, "&Load Segy...\tCtrl-y", "Load Segy Data")
        self.Bind(wx.EVT_MENU, self.segy_input, self.m_load_segy)
        # SEGY NAME LIST
        self.m_segy_submenu = wx.Menu()
        self.seismic_data.Append(-1, 'SEGY Data...', self.m_segy_submenu)
        # GAIN
        self.m_gain = wx.Menu()
        self.seismic_data.Append(-1, 'Gain', self.m_gain)
        # COLOR PALETTE
        self.m_color_palette = wx.Menu()
        self.seismic_data.Append(-1, 'Color Palette', self.m_color_palette)
        self.m_color_palette.Append(901, 'Grey')
        self.Bind(wx.EVT_MENU, self.segy_color_adjustment, id=901)
        self.m_color_palette.Append(902, 'Seismic')
        self.Bind(wx.EVT_MENU, self.segy_color_adjustment, id=902)
        # GAIN INCREASE
        self.m_gain_increase = self.m_gain.Append(-1, "Increase...", "Increase...")
        self.Bind(wx.EVT_MENU, self.gain_increase, self.m_gain_increase)
        # GAIN DECREASE
        self.m_gain_decrease = self.m_gain.Append(-1, "Decrease...", "Decrease...")
        self.Bind(wx.EVT_MENU, self.gain_decrease, self.m_gain_decrease)
        # DRAW MENU
        self.menubar.Append(self.seismic_data, "&Seismic Data")
        # --------------------------------------------------------------------------------------------------------------

        # WELL DATA MENU -----------------------------------------------------------------------------------------------
        self.well_data = wx.Menu()
        self.m_load_well = self.well_data.Append(-1, "&Load well record...\tCtrl-Shift-w", "Load well record")
        self.Bind(wx.EVT_MENU, self.load_well, self.m_load_well)
        # WELL SUBMENU
        self.m_wells_submenu = wx.Menu()
        self.well_data.Append(-1, 'Wells...', self.m_wells_submenu)
        # DRAW MENU
        self.menubar.Append(self.well_data, "&Well Data")
        # --------------------------------------------------------------------------------------------------------------

        # OUTCROP MENU ------------------------------------------------------------------------------------------------
        self.outcrop_file = wx.Menu()
        self.m_load_outcrop_data = self.outcrop_file.Append(-1, "&Load Outcrop data...\tCtrl-Shift-w",
                                                            "Load outcrop data")
        self.Bind(wx.EVT_MENU, self.load_outcrop_data, self.m_load_outcrop_data)

        self.m_outcrop_submenu = wx.Menu()
        self.outcrop_file.Append(-1, 'Outcrop Data...', self.m_outcrop_submenu)

        self.menubar.Append(self.outcrop_file, "&Outcrop data")
        # --------------------------------------------------------------------------------------------------------------

        # MODEL LAYERS MENU --------------------------------------------------------------------------------------------
        self.layer_file = wx.Menu()
        # NEW LAYER
        self.m_new_layer = self.layer_file.Append(-1, "New Layer...", "New Layer...")
        self.Bind(wx.EVT_MENU, self.new_layer, self.m_new_layer)
        # LOAD LAYER
        self.m_load_layer = self.layer_file.Append(-1, "Load Layer...", "Load Layer...")
        self.Bind(wx.EVT_MENU, self.load_layer, self.m_load_layer)
        # TRANSPARENCY
        self.m_layer_transperency = wx.Menu()
        self.layer_file.Append(-1, 'Transparency', self.m_layer_transperency)
        # TRANSPARENCY INCREASE
        self.m_layer_transparency_increase = self.m_layer_transperency.Append(-1, "Increase...", "Increase...")
        self.Bind(wx.EVT_MENU, self.transparency_increase, self.m_layer_transparency_increase)
        # TRANSPARENCY DECREASE
        self.m_layer_transparency_decrease = self.m_layer_transperency.Append(-1, "Decrease...", "Decrease...")
        self.Bind(wx.EVT_MENU, self.transparency_decrease, self.m_layer_transparency_decrease)
        # BULK SHIFT
        self.m_bulk_shift = self.layer_file.Append(-1, "Bulk Shift...", "Bulk Shift...")
        self.Bind(wx.EVT_MENU, self.bulk_shift, self.m_bulk_shift)
        # DELETE LAYER
        self.m_delete_layer = self.layer_file.Append(-1, "Delete Current Layer...", "Delete Current Layer...")
        self.Bind(wx.EVT_MENU, self.delete_layer, self.m_delete_layer)
        # APPEND MENU
        self.menubar.Append(self.layer_file, "&Layers")
        # --------------------------------------------------------------------------------------------------------------

        # PINCH MENU ---------------------------------------------------------------------------------------------------
        layer_file = wx.Menu()
        m_pinch_layer = layer_file.Append(-1, "&Pinch Out Layer...\tCtrl-shift-up", "Pinch Out Layer...")
        self.Bind(wx.EVT_MENU, self.pinch_out_layer, m_pinch_layer)
        m_depinch_layer = layer_file.Append(-1, "&Depinch Layer...\tCtrl-shift-up", "Depinch Layer...")
        self.Bind(wx.EVT_MENU, self.depinch_layer, m_depinch_layer)
        self.menubar.Append(layer_file, "&Pinch Layer")
        # --------------------------------------------------------------------------------------------------------------

        # ATTRIBUTE TABLE MENU -----------------------------------------------------------------------------------------
        attribute_file = wx.Menu()
        m_attribute_table = attribute_file.Append(-1, "&Open Attribute Table...\tCtrl-shift-a",
                                                  "Open Attribute Table...")
        self.Bind(wx.EVT_MENU, self.open_attribute_table, m_attribute_table)
        self.menubar.Append(attribute_file, "&Attribute Table")
        # --------------------------------------------------------------------------------------------------------------

        # HELP MENU ----------------------------------------------------------------------------------------------------
        help_file = wx.Menu()
        m_help = help_file.Append(-1, "&Documentation...", "Open Documentation html...")
        self.Bind(wx.EVT_MENU, self.open_documentation, m_help)
        m_about = help_file.Append(-1, "&About...", "About gmg...")
        self.Bind(wx.EVT_MENU, self.about_gmg, m_about)
        m_legal = help_file.Append(-1, "&Legal...", "Legal...")
        self.Bind(wx.EVT_MENU, self.legal, m_legal)
        self.menubar.Append(help_file, "&Help")

        'SET MENUBAR'
        self.SetMenuBar(self.menubar)
        # --------------------------------------------------------------------------------------------------------------

        'TOOLBAR - (THIS IS THE ICON BAR BELOW THE MENU BAR)'
        self.toolbar = self.CreateToolBar()

        t_save_model = self.toolbar.AddTool(wx.ID_ANY, "Save model", wx.Bitmap(self.gui_icons_dir + 'save_24.png'),
                                            shortHelp="Save model")
        self.Bind(wx.EVT_TOOL, self.save_model, t_save_model)

        t_load_model = self.toolbar.AddTool(wx.ID_ANY, "Load model", wx.Bitmap(self.gui_icons_dir + 'load_24.png'),
                                            shortHelp="Load model")
        self.Bind(wx.EVT_TOOL, self.load_model, t_load_model)

        # t_calc_topo = self.toolbar.AddTool(wx.ID_ANY, "Calculate topography",
        #   wx.Bitmap(self.gui_icons_dir + 'T_24.png'), shortHelp="Calculate topography")
        # self.Bind(wx.EVT_TOOL, self.calc_topo_switch, t_calc_topo)  # FUTURE

        self.t_calc_grav = self.toolbar.AddCheckTool(toolId=wx.ID_ANY, label="Fault pick",
                                                     bitmap1=wx.Bitmap(self.gui_icons_dir + 'G_24.png'),
                                                     bmpDisabled=wx.Bitmap(self.gui_icons_dir + 'G_24.png'),
                                                     shortHelp="Calculate gravity anomaly",
                                                     longHelp="", clientData=None)
        self.Bind(wx.EVT_TOOL, self.calc_grav_switch, self.t_calc_grav)

        self.t_calc_mag = self.toolbar.AddCheckTool(toolId=wx.ID_ANY, label="Fault pick",
                                                    bitmap1=wx.Bitmap(self.gui_icons_dir + 'M_24.png'),
                                                    bmpDisabled=wx.Bitmap(self.gui_icons_dir + 'M_24.png'),
                                                    shortHelp="Calculate magnetic anomaly",
                                                    longHelp="", clientData=None)
        self.Bind(wx.EVT_TOOL, self.calc_mag_switch, self.t_calc_mag)


        self.t_capture_coordinates = self.toolbar.AddCheckTool(toolId=wx.ID_ANY, label="Capture coordinates",
                                                    bitmap1=wx.Bitmap(self.gui_icons_dir + 'C_24.png'),
                                                    bmpDisabled=wx.Bitmap(self.gui_icons_dir + 'C_24.png'),
                                                    shortHelp="Capture coordinates",
                                                    longHelp="", clientData=None)
        self.Bind(wx.EVT_TOOL, self.capture_coordinates, self.t_capture_coordinates)

        t_aspect_increase = self.toolbar.AddTool(wx.ID_ANY, "Aspect increase",
                                                 wx.Bitmap(self.gui_icons_dir + 'large_up_24.png'),
                                                 shortHelp="Aspect increase")
        self.Bind(wx.EVT_TOOL, self.aspect_increase, t_aspect_increase)

        t_aspect_decrease = self.toolbar.AddTool(wx.ID_ANY, "Aspect decrease",
                                                 wx.Bitmap(self.gui_icons_dir + 'large_down_24.png'),
                                                 shortHelp="Aspect decrease")
        self.Bind(wx.EVT_TOOL, self.aspect_decrease, t_aspect_decrease)

        t_aspect_increase2 = self.toolbar.AddTool(wx.ID_ANY, "Aspect increase x2",
                                                  wx.Bitmap(self.gui_icons_dir + 'small_up_24.png'),
                                                  shortHelp="Aspect increase x2")
        self.Bind(wx.EVT_TOOL, self.aspect_increase2, t_aspect_increase2)

        t_aspect_decrease2 = self.toolbar.AddTool(wx.ID_ANY, "Aspect decrease x2",
                                                  wx.Bitmap(self.gui_icons_dir + 'small_down_24.png'),
                                                  shortHelp="Aspect decrease x2")
        self.Bind(wx.EVT_TOOL, self.aspect_decrease2, t_aspect_decrease2)

        self.t_zoom = self.toolbar.AddCheckTool(toolId=wx.ID_ANY, label="Zoom in",
                                                    bitmap1=wx.Bitmap(self.gui_icons_dir + 'zoom_in_24.png'),
                                                    bmpDisabled=wx.Bitmap(self.gui_icons_dir + 'zoom_in_24.png'),
                                                    shortHelp="Zoom in",
                                                    longHelp="", clientData=None)
        self.Bind(wx.EVT_TOOL, self.zoom, self.t_zoom)

        t_zoom_out = self.toolbar.AddTool(wx.ID_ANY, "Zoom out",
                                          wx.Bitmap(self.gui_icons_dir + 'zoom_out_24.png'), shortHelp="Zoom out")
        self.Bind(wx.EVT_TOOL, self.zoom_out, t_zoom_out)

        t_full_extent = self.toolbar.AddTool(wx.ID_ANY, "Full extent",
                                             wx.Bitmap(self.gui_icons_dir + 'full_extent_24.png'),
                                             shortHelp="Full extent")
        self.Bind(wx.EVT_TOOL, self.full_extent, t_full_extent, id=604)

        self.t_pan = self.toolbar.AddCheckTool(toolId=wx.ID_ANY, label="Pan",
                                                    bitmap1=wx.Bitmap(self.gui_icons_dir + 'pan_24.png'),
                                                    bmpDisabled=wx.Bitmap(self.gui_icons_dir + 'pan_24.png'),
                                                    shortHelp="Pan",
                                                    longHelp="", clientData=None)
        self.Bind(wx.EVT_TOOL, self.pan, self.t_pan)

        t_gain_down = self.toolbar.AddTool(wx.ID_ANY, "Gain down",
                                           wx.Bitmap(self.gui_icons_dir + 'left_small_24.png'), shortHelp="Gain down")
        self.Bind(wx.EVT_TOOL, self.gain_decrease, t_gain_down)

        t_gain_up = self.toolbar.AddTool(wx.ID_ANY, "Gain up",
                                         wx.Bitmap(self.gui_icons_dir + 'right_small_24.png'), shortHelp="Gain up")
        self.Bind(wx.EVT_TOOL, self.gain_increase, t_gain_up)

        t_transparency_down = self.toolbar.AddTool(wx.ID_ANY, "Transparency down",
                                                   wx.Bitmap(self.gui_icons_dir + 'large_left_24.png'),
                                                   shortHelp="Transparency down")
        self.Bind(wx.EVT_TOOL, self.transparency_decrease, t_transparency_down)

        # INCREASE TRANSPARENCY ICON
        t_transparency_up = self.toolbar.AddTool(wx.ID_ANY, "Transparency up",
                                                 wx.Bitmap(self.gui_icons_dir + 'large_right_24.png'),
                                                 shortHelp="Transparency up")
        self.Bind(wx.EVT_TOOL, self.transparency_increase, t_transparency_up)

        # LOAD WELL ICON
        t_load_well = self.toolbar.AddTool(wx.ID_ANY, "Load well horizons",
                                           wx.Bitmap(self.gui_icons_dir + 'well_24.png'),
                                           shortHelp="Load well horizons")
        self.Bind(wx.EVT_TOOL, self.load_well, t_load_well)

        # TOOGLE FAULT PICKING MODE
        self.t_toogle_fault_mode = self.toolbar.AddCheckTool(toolId=10000, label="Fault pick",
                                                     bitmap1=wx.Bitmap(self.gui_icons_dir + 'F_24.png'),
                                                     bmpDisabled=wx.Bitmap(self.gui_icons_dir + 'off_F_24.png'),
                                                     shortHelp="Toogle fault picking")
        self.t_toogle_fault_mode.SetDisabledBitmap(wx.Bitmap(self.gui_icons_dir + 'off_F_24.png'))
        self.Bind(wx.EVT_TOOL, self.toogle_fault_mode, self.t_toogle_fault_mode)

        # FAULT PICKER ICON
        self.t_fault_pick = self.toolbar.AddTool(wx.ID_ANY, "Fault pick",
                                                 bitmap=wx.Bitmap(self.gui_icons_dir + 'faultline_24.png'),
                                                 shortHelp="Fault picker")
        self.Bind(wx.EVT_TOOL, self.pick_new_fault, self.t_fault_pick)

        # CREATE TOOLBAR
        self.toolbar.Realize()
        self.toolbar.SetSize((1790, 36))

    def start(self, area, xp, zp):
        """CREATE MPL FIGURE CANVAS"""

        self.fig = plt.figure()  # CREATE MPL FIGURE
        self.canvas = FigureCanvas(self.rightPanel, -1, self.fig)  # CREATE FIGURE CANVAS
        self.nav_toolbar = NavigationToolbar(self.canvas)  # CREATE DEFAULT NAVIGATION TOOLBAR
        self.nav_toolbar.Hide()  # HIDE DEFAULT NAVIGATION TOOLBAR

        # SET DRAW COMMAND WHICH CAN BE CALLED TO REDRAW THE FIGURE'
        self.draw = self.fig.canvas.draw

        # GET THE MODEL DIMENSIONS AND SAMPLE LOCATIONS'
        self.area = area
        self.x1, self.x2, self.z1, self.z2 = 0.001 * np.array(area)
        self.xp = np.array(xp, dtype='f')
        self.zp = np.array(zp, dtype='f')

        # DRAW MAIN PROGRAM WINDOW
        self.draw_main_frame()

        # CONNECT MPL FUNCTIONS
        self.connect()

        # UPDATE DISPLAY
        self.display_info()
        self.size_handler()
        if self.newmodel:
            self.update_layer_data()
        else:
            pass

        # REFRESH SIZER POSITIONS
        self.Hide()
        self.Show()

        # FINIALISE INIT PROCESS
        # self.run_algorithms()
        self.draw()

    def initalise_model(self):
        """INITIALISE OBSERVED DATA AND LAYERS"""

        # INITIALISE MODEL FRAME ATTRIBUTES
        self.showverts = True
        self.exsiting_node = True
        self.nodes = True
        self.zoom_on = False
        self.pan_on = False
        self.node_click_limit = 0.2  # CONTROLS HOW CLOSE A MOUSE CLICK MUST BE TO ACTIVATE A NODE
        self.layer_lock = True
        self.boundary_lock = True
        self.calc_padding = 5000.  # PADDING FOR POTENTIAL FIELD CALCULATION POINTS
        self.padding = 50000.  # PADDING FOR LAYERS
        self.currently_active_layer_id = 0  # LAYER COUNTER
        self.node_layer_reference = 0
        self.select_new_layer_nodes = False  # SWITCH TO TURN ON WHEN MOUSE CLICK IS TO BE CAPTURED FOR A NEW LAYER
        self.index_node = None  # THE ACTIVE NODE
        self.pred_topo = None  # FUTURE - PREDICTED TOPOGRAPHY FROM MOHO (ISOSTATIC FUNC)
        self.predgz = None  # THE CALCULATED GRAVITY RESPONSE
        self.prednt = None  # THE CALCULATED MAGNETIC RESPONSE
        self.gz = None
        self.obs_gz = None
        self.ntz = None
        self.obs_ntz = None
        self.colors = "bmycbmycbmyc"  # OBSERVED COLOR LIST
        self.colors_index = 0
        self.model_azimuth = 0.
        self.mag_observation_elv = 0.
        self.pinch_switch = False

        # INITALISE LAYER LIST
        self.layer_list = []  # LIST HOLDING ALL OF THE LAYER OBJECTS

        # INITIALISE POLYGON LISTS (USED AS MODEL LAYERS)
        self.mag_polygons = []
        self.gravity_polygons = []
        self.polyplots = []
        self.poly_fills = [[]]

        # INITIALISE LAYER LISTS (USED FOR STORING LAYER DATA)
        self.plotx_list = [[]]
        self.ploty_list = [[]]
        self.layer_colors = [[]]
        self.total_layer_count = 0
        self.boundary_lock_list = [0]
        self.boundary_lock_status = ['locked']
        self.layer_lock_list = [0]
        self.layer_lock_status = ['locked']
        self.layer_transparency = 0.4
        self.pinch_node_list = [[], []]
        self.pinch_count = 0

        # INITIALISE XY DATA ATTRIBUTES
        self.observed_xy_data_list = []
        self.xy_data_counter = 0

        # INITIALISE OBSERVED TOPOGRAPHY ATTRIBUTES
        self.observed_topography_list = []
        self.observed_topography_counter = 0
        self.observed_topography_switch = False

        # INITIALISE OBSERVED GRAVITY ATTRIBUTES
        self.observed_gravity_list = []
        self.observed_gravity_counter = 0
        self.observed_gravity_switch = False
        # INITIALISE MODELLING GRAVITY ATTRIBUTES
        self.background_density = 0
        self.absolute_densities = True
        self.calc_grav_switch = False
        self.obs_gravity_data_for_rms = []  # OBSERVED DATA LIST TO BE COMPARED TO CALCULATED
        self.grav_rms_value = 0.  # TOTAL RMS MISFIT VALUE
        self.grav_residuals = []  # CALCULATED RESIDUAL

        # INITIALISE OBSERVED MAGNETIC ATTRIBUTES
        self.observed_magnetic_list = []
        self.observed_magnetic_counter = 0
        self.observed_magnetic_switch = False
        # INITIALISE MODELLING MAGNETIC ATTRIBUTES
        self.earth_field = 0.
        self.profile_azimuth = 0.
        self.calc_mag_switch = False
        self.obs_mag_data_for_rms = []  # OBSERVED DATA LIST TO BE COMPARED TO CALCULATED
        self.mag_rms_value = 0.  # TOTAL RMS MISFIT VALUE (SINGLE INTEGER)
        self.mag_residuals = []  # CALCULATED RESIDUAL

        # INITIALISE GEOLOGICAL CONTACT ATTRIBUTES
        self.outcrop_data_list = []
        self.outcrop_data_count = 0

        # INITIALISE Well ATTRIBUTES
        self.well_data_list = []
        self.well_counter = 0

        # INITIALISE SEISMIC ATTRIBUTES
        self.segy_data_list = []
        self.segy_counter = 0
        self.segy_color_map = cm.gray
        self.segy_gain_neg = -4.0
        self.segy_gain_pos = 4.0

        # INITIALISE LAYER ATTRIBUTES
        self.densities = [0.]
        self.reference_densities = [2.67]
        self.susceptibilities = [0.]
        self.angle_a = [0.]
        self.angle_b = [0.]
        self.layers_calculation_switch = [1]

        # INITIALISE FAULT ATTRIBUTES
        self.fault_picking_switch = False
        self.faults = [[]]
        self.fault_names_list = []
        self.fault_counter = 0
        self.current_fault_index = 0
        self.fault_colors = ['k']
        self.fault_x_coords_list = [[]]
        self.fault_y_coords_list = [[]]
        self.selected_node = None

        # INITIALIZE COORDINATE CAPTURE
        self.capture = False
        self.linex = []
        self.liney = []

    def draw_main_frame(self):
        """
        DRAW THE PROGRAM CANVASES
        docs: https://matplotlib.org/api/axes_api.html
        """
        self.columns = 87  # NUMBER OF COLUMNS THE MODEL FRAMES WILL TAKE UP (89/100)
        self.x_orig = 10  # X ORIGIN OF MODEL FRAMES (RELATIVE TO 0 AT LEFT MARGIN)

        # TOPOGRAPHY CANVAS
        self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=2, colspan=self.columns)
        self.topo_frame.set_ylabel("Topo (km)")
        self.topo_frame.set_navigate(False)
        self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
        self.topo_frame.grid()
        self.topo_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        self.topo_d_frame = self.topo_frame.twinx()
        self.topo_d_frame.set_navigate(False)
        self.topo_d_frame.set_ylabel("dt/dx")
        self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
        self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        # GRAVITY CANVAS
        self.gravity_frame = plt.subplot2grid((26, 100), (2, self.x_orig), rowspan=3, colspan=self.columns)
        self.gravity_frame.set_navigate(False)
        self.gravity_frame.set_ylabel("Grav (mGal)")
        self.gravity_frame.xaxis.set_major_formatter(plt.NullFormatter())
        self.gravity_frame.grid()
        self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        self.gravity_d_frame = self.gravity_frame.twinx()
        self.gravity_d_frame.set_navigate(False)
        self.gravity_d_frame.set_ylabel("dg/dx")
        self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
        self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        # MAGNETIC CANVAS
        self.magnetic_frame = plt.subplot2grid((26, 100), (5, self.x_orig), rowspan=3, colspan=self.columns)
        self.magnetic_frame.set_ylabel("Mag (nT)")
        self.magnetic_frame.set_navigate(False)
        self.magnetic_frame.xaxis.set_major_formatter(plt.NullFormatter())
        self.magnetic_frame.grid()
        self.magnetic_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        self.magnetic_d_frame = self.magnetic_frame.twinx()
        self.magnetic_d_frame.set_ylabel("dnt/dx")
        self.magnetic_d_frame.set_navigate(False)
        self.magnetic_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
        self.magnetic_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        # MODEL CANVAS
        self.model_frame = plt.subplot2grid((26, 100), (8, self.x_orig), rowspan=17, colspan=self.columns)
        self.model_frame.set_ylabel("Depth (km)")
        self.model_frame.set_xlabel("x (km)")

        # CREATE DENSITY COLOUR BAR FOR COLORING LAYERS
        colormap = matplotlib.cm.coolwarm
        cnorm = colors.Normalize(vmin=-0.8, vmax=0.8)
        self.colormap = cm.ScalarMappable(norm=cnorm, cmap=colormap)

        # self.scale_canvas = plt.subplot2grid((24, 12), (23, 10), rowspan=1, colspan=2)
        # self.cb1 = matplotlib.colorbar.ColorbarBase(self.scale_canvas,
        #                                             cmap=colormap, norm=cnorm, orientation='horizontal')
        # self.cb1.ax.tick_params(labelsize=8)
        # self.cb1.set_label('Density contrast ($kg/m^{3}$)', fontsize=10)

        # SET CANVAS LIMITS
        self.model_frame.set_xlim(self.x1, self.x2)
        self.model_frame.set_ylim(self.z2, self.z1)
        self.model_frame.grid()
        self.topo_frame.set_xlim(self.model_frame.get_xlim())
        self.gravity_frame.set_xlim(self.model_frame.get_xlim())
        self.magnetic_frame.set_xlim(self.model_frame.get_xlim())
        self.fig.subplots_adjust(top=0.99, left=-0.045, right=0.99, bottom=0.02,
                                 hspace=1.5)

        # ADD FIRST LAYER
        if self.newmodel:
            # CREATE LAYER0 - PLACE HOLDER FOR THE TOP OF THE MODEL - NOT ACCESSIBLE BY USER
            layer0 = Layer()
            layer0.type = str('fixed')
            # CREATE THE XY NODES
            layer0.x_nodes = [-(float(self.padding)), 0., self.x2, self.x2 + (float(self.padding))]
            layer0.y_nodes = [0.001, 0.001, 0.001, 0.001]

            # SET CURRENT NODES
            self.current_x_nodes = layer0.x_nodes
            self.current_y_nodes = layer0.y_nodes

            # DRAW THE CURRENTLY ACTIVE LAYER (THE NODES THAT CAN BE INTERACTED WITH)
            self.currently_active_layer, = self.model_frame.plot(layer0.x_nodes, layer0.y_nodes, marker='o', color='k',
                                                                 linewidth=1.0, alpha=0.5, picker=True)
            # ADD THE LAYER LINE TO THE PLOT
            layer0.node_mpl_actor = self.model_frame.plot(layer0.x_nodes, layer0.y_nodes, color='black', linewidth=1.0,
                                                          alpha=1.0)
            # ADD THE LAYER POLYGON FILL TO THE PLOT
            layer0.polygon_mpl_actor = self.model_frame.fill(layer0.x_nodes, layer0.y_nodes, color='blue',
                                                             alpha=self.layer_transparency, closed=True, linewidth=None,
                                                             ec=None)

            # SET THE CURRENTLY ACTIVE RED NODE
            self.current_node = self.model_frame.scatter(-40000., 0, s=50, color='r', zorder=10)  # PALCE HOLDER ONLY

            self.layer_list.append(layer0)
            print self.layer_list
            print len(self.layer_list)

        # ADDITIONAL MAIN FRAME WIDGETS - PLACED ON LEFT HAND SIDE OF THE FRAME
        'Make Attribute Label'
        self.attr_text = wx.StaticText(self.fold_panel_one, -1, label="", style=wx.ALIGN_LEFT)
        self.attr_text2 = wx.StaticText(self.fold_panel_one, -1, label="", style=wx.ALIGN_LEFT)
        'Make NODE X Y Label'
        self.node_text = wx.StaticText(self.fold_panel_one, -1, label="Node Position:", style=wx.ALIGN_LEFT)
        'Make density spinner'
        self.density_text = wx.StaticText(self.fold_panel_one, -1, label="Density:         ", style=wx.ALIGN_LEFT)
        self.density_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-5, max_val=5, increment=0.001, value=0.00)
        self.density_input.SetFormat("%f")
        self.density_input.SetDigits(4)
        'Make refernece density spinner'
        self.ref_density_text = wx.StaticText(self.fold_panel_one, -1, label="Reference:     \nDensity",
                                              style=wx.ALIGN_LEFT)
        self.ref_density_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-5, max_val=5, increment=0.001,
                                              value=0.00)
        self.ref_density_input.SetFormat("%f")
        self.ref_density_input.SetDigits(4)
        'Make susceptibility spinner'
        self.susceptibility_text = wx.StaticText(self.fold_panel_one, -1, label="Susceptibility:", style=wx.ALIGN_LEFT)
        self.susceptibility_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-2.0, max_val=2.0, increment=0.00001,
                                                 value=0.00)
        self.susceptibility_input.SetFormat("%f")
        self.susceptibility_input.SetDigits(6)
        'MAKE ANGLE A SPINNER'
        self.angle_a_text = wx.StaticText(self.fold_panel_one, -1, label="Angle A (Inc): ", style=wx.ALIGN_LEFT)
        self.angle_a_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=0.0, max_val=90.0, increment=1.0, value=0.0)
        self.angle_a_input.SetFormat("%f")
        self.angle_a_input.SetDigits(1)
        'Make ANGLE B SPINNER'
        self.angle_b_text = wx.StaticText(self.fold_panel_one, -1, label="Angle B (Dec):", style=wx.ALIGN_LEFT)
        self.angle_b_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=0.0, max_val=180.0, increment=1.0, value=0.0)
        self.angle_b_input.SetFormat("%f")
        self.angle_b_input.SetDigits(1)
        'Make well text size slider'
        self.text_size_text = wx.StaticText(self.fold_panel_one, -1, label="Label Text Size:")
        self.text_size_input = wx.Slider(self.fold_panel_one, value=1, minValue=1, maxValue=20., size=(175, -1),
                                         style=wx.SL_HORIZONTAL)
        'Make Node XY spinners'
        self.x_text = wx.StaticText(self.fold_panel_one, -1, label="X:")
        self.x_input = fs.FloatSpin(self.fold_panel_one, -1, increment=0.001, value=0.00)
        self.x_input.SetDigits(4)
        self.y_text = wx.StaticText(self.fold_panel_one, -1, label="Y:")
        self.y_input = fs.FloatSpin(self.fold_panel_one, -1, increment=0.001, value=0.00)
        self.y_input.SetDigits(4)
        'Make Set button'
        self.node_set_button = wx.Button(self.fold_panel_one, -1, "Set layer attributes")

        'INITALISE CALCULATED P.F. LINES'
        self.predplot, = self.gravity_frame.plot([], [], '-r', linewidth=2, alpha=0.5)
        self.grav_rms_plot, = self.gravity_frame.plot([], [], color='purple', linewidth=1.5, alpha=0.5)
        self.prednt_plot, = self.magnetic_frame.plot([], [], '-g', linewidth=2, alpha=0.5)
        self.mag_rms_plot, = self.magnetic_frame.plot([], [], color='purple', linewidth=1.5, alpha=0.5)

        'MAKE LAYER TREE'
        self.tree = ct.CustomTreeCtrl(self.fold_panel_two, -1, size=(200, 280),
                                      agwStyle=wx.TR_DEFAULT_STYLE | wx.TR_EDIT_LABELS | wx.TR_HIDE_ROOT)
        self.tree.SetIndent(0.0)
        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.on_activated, self.tree)
        self.Bind(wx.EVT_TREE_BEGIN_LABEL_EDIT, self.on_begin_edit_label, self.tree)
        self.Bind(wx.EVT_TREE_END_LABEL_EDIT, self.on_end_edit_label, self.tree)

        'TREE ATTRIBUTES'
        self.root = self.tree.AddRoot("Layers:")
        self.tree.SetItemPyData(self.root, None)
        self.tree_items = ["Layer 0"]
        self.Bind(ct.EVT_TREE_ITEM_CHECKED, self.item_checked, self.tree)

        self.error = 0.
        self.last_layer = 0

        'MAKE FAULT TREE'
        self.fault_tree = ct.CustomTreeCtrl(self.fold_panel_three, -1, size=(200, 331),
                                            agwStyle=wx.TR_DEFAULT_STYLE | wx.TR_EDIT_LABELS | wx.TR_HIDE_ROOT)
        self.fault_tree.SetIndent(0.0)

        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.fault_activated, self.fault_tree)
        ### self.Bind(wx.EVT_TREE_BEGIN_LABEL_EDIT, self.on_begin_edit_label, self.fault_tree)
        ### self.Bind(wx.EVT_TREE_END_LABEL_EDIT, self.on_end_edit_label, self.fault_tree)

        'TREE ATTRIBUTES'
        self.fault_tree_root = self.fault_tree.AddRoot("Faults:")
        self.fault_tree.SetItemPyData(self.fault_tree_root, None)
        self.fault_tree_items = []
        self.Bind(ct.EVT_TREE_ITEM_CHECKED, self.fault_checked, self.fault_tree)

        'UPDATE INFO BAR'
        self.display_info()

        'DRAW MAIN'
        self.draw()

    def frame_adjustment(self, event):
        """FIND WHICH FRAME IS REFERENCED & CHANGE SWITCH"""

        self.current_xlim = self.model_frame.get_xlim()
        self.current_ylim = self.model_frame.get_ylim()

        if event.Id == 601:
            if self.topo_frame_switch is True:
                self.topo_frame.set_visible(False)
                self.topo_d_frame.set_visible(False)
                self.topo_frame_switch = False
            else:
                self.topo_frame_switch = True
                self.topo_frame.set_visible(True)
                self.topo_d_frame.set_visible(True)
        if event.Id == 602:
            if self.gravity_frame_switch is True:
                self.gravity_frame.set_visible(False)
                self.gravity_d_frame.set_visible(False)
                self.gravity_frame_switch = False
            else:
                self.gravity_frame_switch = True
                self.gravity_frame.set_visible(True)
                self.gravity_d_frame.set_visible(True)
        if event.Id == 603:
            if self.magnetic_frame_switch is True:
                self.magnetic_frame.set_visible(False)
                self.magnetic_d_frame.set_visible(False)
                self.magnetic_frame_switch = False
            else:
                self.magnetic_frame_switch = True
                self.magnetic_frame.set_visible(True)
                self.magnetic_d_frame.set_visible(True)

        # ADJUST FRAME SIZING AND SET PROGRAM WINDOW
        if self.topo_frame_switch is True and self.gravity_frame_switch is True and self.magnetic_frame_switch is True:
            'TRUE TRUE TRUE'
            'TOPO CANVAS'
            self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=2, colspan=self.columns)
            self.topo_frame.set_ylabel("(m)")
            self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_frame.grid()
            self.topo_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.topo_d_frame = self.topo_frame.twinx()
            self.topo_d_frame.set_ylabel("dt/dx")
            self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            'GRAV CANVAS'
            self.gravity_frame = plt.subplot2grid((26, 100), (2, self.x_orig), rowspan=3, colspan=self.columns)
            self.gravity_frame.set_ylabel("(mGal)")
            self.gravity_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_frame.grid()
            self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.gravity_d_frame = self.gravity_frame.twinx()
            self.gravity_d_frame.set_ylabel("dg/dx")
            self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            'MAG CANVAS'
            self.magnetic_frame = plt.subplot2grid((26, 100), (5, self.x_orig), rowspan=3, colspan=self.columns)
            self.magnetic_frame.set_ylabel("(nT)")
            self.magnetic_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.magnetic_frame.grid()
            self.magnetic_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.magnetic_d_frame = self.magnetic_frame.twinx()
            self.magnetic_d_frame.set_ylabel("dnt/dx")
            self.magnetic_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.magnetic_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        elif self.topo_frame_switch is False and self.gravity_frame_switch is True and \
                self.magnetic_frame_switch is True:
            'FALSE TRUE TRUE'
            'TOPO CANVAS'
            # HIDDEN
            'GRAV CANVAS'
            self.gravity_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=4, colspan=self.columns)
            self.gravity_frame.set_ylabel("(mGal)")
            self.gravity_frame.grid()
            self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.gravity_d_frame = self.gravity_frame.twinx()
            self.gravity_d_frame.set_ylabel("dg/dx")
            self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            'MAG CANVAS'
            self.magnetic_frame = plt.subplot2grid((26, 100), (4, self.x_orig), rowspan=4, colspan=self.columns)
            self.magnetic_frame.set_ylabel("(nT)")
            self.magnetic_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.magnetic_frame.grid()
            self.magnetic_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.magnetic_d_frame = self.magnetic_frame.twinx()
            self.magnetic_d_frame.set_ylabel("dnt/dx")
            self.magnetic_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.magnetic_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        elif self.topo_frame_switch is True and self.gravity_frame_switch is False and \
                self.magnetic_frame_switch is True:
            'TRUE FALSE TRUE'
            'TOPO CANVAS'
            self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=2, colspan=self.columns)
            self.topo_frame.set_ylabel("(m)")
            self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_frame.grid()

            self.topo_d_frame = self.topo_frame.twinx()
            self.topo_d_frame.set_ylabel("dt/dx")
            self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            'GRAV CANVAS'
            # HIDDEN
            'MAG CANVAS'
            self.magnetic_frame = plt.subplot2grid((26, 100), (2, self.x_orig), rowspan=6, colspan=self.columns)
            self.magnetic_frame.set_ylabel("(nT)")
            self.magnetic_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.magnetic_frame.grid()
            self.magnetic_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.magnetic_d_frame = self.magnetic_frame.twinx()
            self.magnetic_d_frame.set_ylabel("dnt/dx")
            self.magnetic_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.magnetic_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        elif self.topo_frame_switch is True and self.gravity_frame_switch is True and \
                self.magnetic_frame_switch is False:
            'TRUE TRUE FALSE'
            'TOPO CANVAS'
            self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=2, colspan=self.columns)
            self.topo_frame.set_ylabel("(m)")
            self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_frame.grid()

            self.topo_d_frame = self.topo_frame.twinx()
            self.topo_d_frame.set_ylabel("dt/dx")
            self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            'GRAV CANVAS'
            self.gravity_frame = plt.subplot2grid((26, 100), (2, self.x_orig), rowspan=6, colspan=self.columns)
            self.gravity_frame.set_ylabel("(mGal)")
            self.gravity_frame.grid()
            self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.gravity_d_frame = self.gravity_frame.twinx()
            self.gravity_d_frame.set_ylabel("dg/dx")
            self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            'MAG CANVAS'
            # HIDDEN

        elif self.topo_frame_switch is False and self.gravity_frame_switch is False and \
                self.magnetic_frame_switch is True:
            'FALSE FALSE TRUE'
            'TOPO CANVAS'
            # HIDDEN
            'GRAV CANVAS'
            # HIDDEN
            'MAG CANVAS'
            self.magnetic_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=8, colspan=self.columns)
            self.magnetic_frame.set_ylabel("(nT)")
            self.magnetic_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.magnetic_frame.grid()
            self.magnetic_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.magnetic_d_frame = self.magnetic_frame.twinx()
            self.magnetic_d_frame.set_ylabel("dnt/dx")
            self.magnetic_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.magnetic_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        elif self.topo_frame_switch is False and self.gravity_frame_switch is True and \
                self.magnetic_frame_switch is False:
            'FALSE TRUE FALSE'
            'TOPO CANVAS'
            # HIDDEN
            'GRAV CANVAS'
            self.gravity_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=8, colspan=self.columns)
            self.gravity_frame.set_ylabel("(mGal)")
            self.gravity_frame.grid()
            self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.gravity_d_frame = self.gravity_frame.twinx()
            self.gravity_d_frame.set_ylabel("dg/dx")
            self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            'MAG CANVAS'
            # HIDDEN

        elif self.topo_frame_switch is True and self.gravity_frame_switch is False and \
                self.magnetic_frame_switch is False:
            'TRUE FALSE FALSE'
            'TOPO CANVAS'
            self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=8, colspan=self.columns)
            self.topo_frame.set_ylabel("(m)")
            self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_frame.grid()

            self.topo_d_frame = self.topo_frame.twinx()
            self.topo_d_frame.set_ylabel("dt/dx")
            self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            'GRAV CANVAS'
            # HIDDEN
            'MAG CANVAS'
            # HIDDEN

        elif self.topo_frame_switch is False and self.gravity_frame_switch is False and \
                self.magnetic_frame_switch is False:
            pass
            'FALSE FALSE FALSE'
            'TOPO CANVAS'
            # HIDDEN
            'GRAV CANVAS'
            # HIDDEN
            'MAG CANVAS'
            # HIDDEN

        # SET CANVAS LIMITS
        self.model_frame.set_xlim(self.current_xlim)
        self.model_frame.set_ylim(self.current_ylim)
        self.model_frame.grid()
        if self.topo_frame is not None:
            self.topo_frame.set_xlim(self.model_frame.get_xlim())
            self.topo_d_frame.set_xlim(self.model_frame.get_xlim())
        if self.gravity_frame is not None:
            self.gravity_frame.set_xlim(self.model_frame.get_xlim())
            self.gravity_d_frame.set_xlim(self.model_frame.get_xlim())
        if self.magnetic_frame is not None:
            self.magnetic_frame.set_xlim(self.model_frame.get_xlim())
            self.magnetic_d_frame.set_xlim(self.model_frame.get_xlim())
        self.fig.subplots_adjust(top=0.99, left=-0.045, right=0.99, bottom=0.02, hspace=1.5)

        # INITALISE CALCULATED P.F. LINES
        if self.gravity_frame is not None:
            self.predplot, = self.gravity_frame.plot([], [], '-r', linewidth=2)
            self.grav_rms_plot, = self.gravity_frame.plot([], [], color='purple', linewidth=1.5)
        if self.magnetic_frame is not None:
            self.prednt_plot, = self.magnetic_frame.plot([], [], '-g', linewidth=2)
            self.mag_rms_plot, = self.magnetic_frame.plot([], [], color='purple', linewidth=1.5)

        # PLOT OBSERVED TOPO DATA
        if self.topo_frame is not None:
            # REPLOT OBSERVED TOPOGRAPHY DATA
            for x in range(len(self.observed_topography_list)):
                if self.observed_topography_list[x] is not None:
                    # DRAW DATA IN MODEL FRAME
                    if self.observed_topography_list[x].type is not "derivative":
                        self.observed_topography_list[x].mpl_actor = self.topo_frame.scatter(
                                                                           self.observed_topography_list[x].data[:, 0],
                                                                           self.observed_topography_list[x].data[:, 1],
                                                                           marker='o', s=5,
                                                                           color=self.observed_topography_list[x].color,
                                                                           gid=self.observed_topography_list[x].id)
                    else:
                        self.observed_topography_list[x].mpl_actor = self.topo_d_frame.scatter(
                                                                           self.observed_topography_list[x].data[:, 0],
                                                                           self.observed_topography_list[x].data[:, 1],
                                                                           marker='o', s=5,
                                                                           color=self.observed_topography_list[x].color,
                                                                           gid=self.observed_topography_list[x].id)
        # PLOT OBSERVED GRAVITY DATA
        if self.gravity_frame is not None:
            # REPLOT OBSERVED GRAVITY DATA
            for x in range(len(self.observed_gravity_list)):
                if self.observed_gravity_list[x] is not None:
                    # DRAW DATA IN MODEL FRAME
                    if self.observed_gravity_list[x].type is not "derivative":
                        self.observed_gravity_list[x].mpl_actor = self.gravity_frame.scatter(
                                                                            self.observed_gravity_list[x].data[:, 0],
                                                                            self.observed_gravity_list[x].data[:, 1],
                                                                            marker='o', s=5,
                                                                            color=self.observed_gravity_list[x].color,
                                                                            gid=11000+self.observed_gravity_list[x].id)
                    else:
                        self.observed_gravity_list[x].mpl_actor = self.gravity_d_frame.scatter(
                                                                            self.observed_gravity_list[x].data[:, 0],
                                                                            self.observed_gravity_list[x].data[:, 1],
                                                                            marker='o', s=5,
                                                                            color=self.observed_gravity_list[x].color,
                                                                            gid=11000+self.observed_gravity_list[x].id)

        # PLOT OBSERVED MAGNETIC DATA
        if self.magnetic_frame is not None:
            # REPLOT OBSERVED MAGNETIC DATA
            for x in range(len(self.observed_magnetic_list)):
                if self.observed_magnetic_list[x] is not None:
                    # DRAW DATA IN MODEL FRAME
                    if self.observed_magnetic_list[x].type is not "derivative":
                        self.observed_magnetic_list[x].mpl_actor = self.magnetic_frame.scatter(
                                                                            self.observed_magnetic_list[x].data[:, 0],
                                                                            self.observed_magnetic_list[x].data[:, 1],
                                                                            marker='o', s=5,
                                                                            color=self.observed_magnetic_list[x].color,
                                                                            gid=self.observed_magnetic_list[x].id)
                    else:
                        self.observed_magnetic_list[x].mpl_actor = self.magnetic_d_frame.scatter(
                                                                            self.observed_magnetic_list[x].data[:, 0],
                                                                            self.observed_magnetic_list[x].data[:, 1],
                                                                            marker='o', s=5,
                                                                            color=self.observed_magnetic_list[x].color,
                                                                            gid=self.observed_magnetic_list[x].id)

        # UPDATE FRAMES
        self.model_frame.grid(True)
        self.run_algorithms()
        self.draw()
        self.set_frame_limits()

    def size_handler(self):
        """ CREATE CANVAS BOX"""
        self.canvas_box = wx.BoxSizer(wx.HORIZONTAL)
        self.canvas_box.Add(self.canvas, 1, wx.ALL | wx.ALIGN_CENTER | wx.EXPAND, border=2)

        'LAYER ATTRIBUTES'
        self.attr_box = wx.BoxSizer(wx.VERTICAL)

        self.attr_text_box = wx.BoxSizer(wx.HORIZONTAL)
        self.attr_text_box.Add(self.attr_text, 0, wx.ALL | wx.ALIGN_CENTER | wx.EXPAND, 2)
        self.attr_text_box.Add(self.attr_text2, 0, wx.ALL | wx.ALIGN_CENTER | wx.EXPAND, 2)
        self.attr_box.Add(self.attr_text_box, 0, wx.ALL | wx.ALIGN_CENTER | wx.EXPAND)

        self.attr_box.Add(wx.StaticLine(self.fold_panel_one), 0, wx.ALL | wx.EXPAND, 5)

        'Density'
        self.density_box = wx.BoxSizer(wx.HORIZONTAL)
        self.density_box.Add(self.density_text, 0, wx.ALL | wx.LEFT | wx.EXPAND, 2)
        self.density_box.Add(self.density_input, 0, wx.ALL | wx.RIGHT | wx.EXPAND)
        self.attr_box.Add(self.density_box, 0, wx.LEFT | wx.EXPAND, 1)
        self.attr_box.Add(wx.StaticLine(self.fold_panel_one), 0, wx.ALL | wx.EXPAND, 5)
        'Reference Density'
        self.ref_density_box = wx.BoxSizer(wx.HORIZONTAL)
        self.ref_density_box.Add(self.ref_density_text, 0, wx.ALL | wx.LEFT | wx.EXPAND, 2)
        self.ref_density_box.Add(self.ref_density_input, 0, wx.ALL | wx.RIGHT | wx.EXPAND)
        self.attr_box.Add(self.ref_density_box, 0, wx.LEFT | wx.EXPAND, 1)
        self.attr_box.Add(wx.StaticLine(self.fold_panel_one), 0, wx.ALL | wx.EXPAND, 5)
        'Susceptibility'
        self.susept_box = wx.BoxSizer(wx.HORIZONTAL)
        self.susept_box.Add(self.susceptibility_text, 0, wx.ALL | wx.LEFT | wx.EXPAND, 3)
        self.susept_box.Add(self.susceptibility_input, 0, wx.ALL | wx.RIGHT | wx.EXPAND)
        self.attr_box.Add(self.susept_box, 0, wx.LEFT | wx.EXPAND)
        'Angle A'
        self.angle_a_box = wx.BoxSizer(wx.HORIZONTAL)
        self.angle_a_box.Add(self.angle_a_text, 0, wx.ALL | wx.LEFT | wx.EXPAND, 0)
        self.angle_a_box.Add(self.angle_a_input, 0, wx.ALL | wx.RIGHT | wx.EXPAND)
        self.attr_box.Add(self.angle_a_box, 0, wx.LEFT | wx.EXPAND, 5)
        'Angle B'
        self.angle_b_box = wx.BoxSizer(wx.HORIZONTAL)
        self.angle_b_box.Add(self.angle_b_text, 0, wx.ALL | wx.LEFT | wx.EXPAND, 0)
        self.angle_b_box.Add(self.angle_b_input, 0, wx.ALL | wx.RIGHT | wx.EXPAND)
        self.attr_box.Add(self.angle_b_box, 0, wx.LEFT | wx.EXPAND, 5)
        self.attr_box.Add(wx.StaticLine(self.fold_panel_one), 0, wx.ALL | wx.EXPAND, 5)
        'Well label text size'
        self.attr_box.Add(self.text_size_text, 0, wx.ALL | wx.ALIGN_CENTER | wx.ALIGN_TOP | wx.EXPAND)
        self.attr_box.Add(self.text_size_input, 0, wx.ALL | wx.ALIGN_CENTER | wx.ALIGN_TOP | wx.EXPAND)
        self.attr_box.Add(wx.StaticLine(self.fold_panel_one), 0, wx.ALL | wx.EXPAND, 5)
        'XY Nodes'
        self.node_box = wx.BoxSizer(wx.HORIZONTAL)
        self.node_box.Add(self.node_text, 0, wx.ALL | wx.ALIGN_CENTER | wx.ALIGN_TOP | wx.EXPAND, 2)
        self.attr_box.Add(self.node_box, 0, wx.ALL | wx.ALIGN_CENTER | wx.ALIGN_TOP | wx.EXPAND)
        self.attr_box.Add(wx.StaticLine(self.fold_panel_one), 0, wx.ALL | wx.EXPAND, 5)
        'X NODE'
        self.x_box = wx.BoxSizer(wx.HORIZONTAL)
        self.x_box.Add(self.x_text, 0, wx.LEFT | wx.EXPAND, 8)
        self.x_box.Add(self.x_input, 0, wx.LEFT | wx.EXPAND, 3)
        self.attr_box.Add(self.x_box, 0, wx.ALL | wx.LEFT | wx.ALIGN_TOP | wx.EXPAND)
        'Y NODE'
        self.y_box = wx.BoxSizer(wx.HORIZONTAL)
        self.y_box.Add(self.y_text, 0, wx.ALL | wx.EXPAND, 3)
        self.y_box.Add(self.y_input, 0, wx.ALL | wx.RIGHT | wx.EXPAND)
        self.attr_box.Add(self.y_box, 0, wx.ALL | wx.LEFT | wx.ALIGN_TOP | wx.EXPAND, 5)
        self.attr_box.Add(wx.StaticLine(self.fold_panel_one), 0, wx.ALL | wx.EXPAND, 5)
        'SET BUTTON'
        self.attr_box.Add(self.node_set_button, 0, wx.ALL | wx.LEFT | wx.EXPAND, 5)

        'CREATE LAYER TREE BOX'
        self.tree_box = wx.BoxSizer(wx.VERTICAL)
        self.tree_box.Add(self.tree, 1, wx.TOP | wx.ALIGN_CENTER | wx.EXPAND, border=20)

        'CREATE FAULT TREE BOX'
        self.fault_tree_box = wx.BoxSizer(wx.VERTICAL)
        self.fault_tree_box.Add(self.fault_tree, 1, wx.TOP | wx.ALIGN_CENTER | wx.EXPAND, border=20)

        '#PLACE BOX SIZERS IN CORRECT PANELS'
        self.fold_panel_one.SetSizerAndFit(self.attr_box)
        self.fold_panel_two.SetSizerAndFit(self.tree_box)
        self.fold_panel_three.SetSizerAndFit(self.fault_tree_box)
        self.leftPanel.SetSizer(self.splitter_left_panel_sizer)
        self.fold_panel_one.Collapse()
        self.fold_panel_one.Expand()
        self.fold_panel_two.Collapse()
        self.fold_panel_two.Expand()
        self.fold_panel_three.Collapse()
        self.fold_panel_three.Expand()

        self.rightPanel.SetSizerAndFit(self.canvas_box)
        self.rightPanel.SetSize(self.GetSize())

    def show_controls(self, event):
        """ SHOW CONTROL PANE"""
        self.mgr.GetPaneByName('left').Show()
        self.mgr.Update()

    def show_console(self, event):
        """ SHOW PYTHON CONSOLE"""
        self.mgr.GetPaneByName('console').Show()
        self.mgr.Update()

    def new_model(self, event):
        """NEW MODEL DIALOG BOX"""

        new_model_dialogbox = NewModelDialog(self, -1, 'Create New Model')
        new_model_dialogbox.ShowModal()
        if new_model_dialogbox.set_button is True:
            self.newmodel = True
            new_x1, new_x2 = float(new_model_dialogbox.x1) * 1000., float(new_model_dialogbox.x2) * 1000.
            new_z1, new_z2 = float(new_model_dialogbox.z1) * 1000., float(new_model_dialogbox.z2) * 1000.
            xp_inc = float(new_model_dialogbox.xp_inc) * 1000.

            self.layer_lines = [[]]
            self.polygon_fills = [[]]

            'INITAISE THE MODEL PARAMETERS'
            self.initalise_model()

            'ENSURES FOLD PANELS ARE VISIBLE'
            self.controls_panel_bar_one.Expand(self.fold_panel_one)
            self.controls_panel_bar_two.Expand(self.fold_panel_two)
            self.controls_panel_bar_three.Expand(self.fold_panel_three)

            'CREATE MODEL'
            self.model_aspect = 1.
            area = (new_x1, new_x2, new_z1, new_z2)
            xp = np.arange(new_x1 - self.calc_padding, new_x2 + self.calc_padding, xp_inc)
            zp = np.zeros_like(xp)
            self.start(area, xp, zp)
        else:
            self.newmodel = False
            return  # THE USER CHANGED THERE MIND

    def modify_model_dimensions(self, event):
        """MODIFY MODEL DIMENSIONS"""
        modify_model_dialogbox = NewModelDialog(self, -1, 'Modify Current Model', self.x1, self.x2, self.z1, self.z2)
        answer = modify_model_dialogbox.ShowModal()
        new_x1, new_x2 = float(modify_model_dialogbox.x1) * 1000., float(modify_model_dialogbox.x2) * 1000.
        new_z1, new_z2 = float(modify_model_dialogbox.z1) * 1000., float(modify_model_dialogbox.z2) * 1000.
        xp_inc = float(modify_model_dialogbox.xp_inc) * 1000.

        self.area = (new_x1, new_x2, new_z1, new_z2)
        self.model_frame.set_ylim(new_z2 / 1000., new_z1 / 1000.)
        self.model_frame.set_xlim(new_x1 / 1000., new_x2 / 1000.)
        xp = np.arange(new_x1 - self.calc_padding, new_x2 + self.calc_padding, xp_inc)
        self.xp = np.array(xp, dtype='f')
        zp = np.zeros_like(xp)
        self.zp = np.array(zp, dtype='f')
        self.update_layer_data()
        self.run_algorithms()

    def connect(self):
        """CONNECT MOUSE AND EVENT BINDINGS"""
        self.fig.canvas.mpl_connect('button_press_event', self.button_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self.move)
        self.fig.canvas.mpl_connect('button_release_event', self.button_release)
        self.fig.canvas.mpl_connect('key_press_event', self.key_press)
        # self.fig.canvas.mpl_connect('pick_event', self.on_pick)

        'Connect wx.widgets'
        self.density_input.Bind(fs.EVT_FLOATSPIN, self.set_density)
        self.ref_density_input.Bind(fs.EVT_FLOATSPIN, self.set_reference_density)
        self.susceptibility_input.Bind(fs.EVT_FLOATSPIN, self.set_susceptibility)
        self.angle_a_input.Bind(fs.EVT_FLOATSPIN, self.set_angle_a)
        self.angle_b_input.Bind(fs.EVT_FLOATSPIN, self.set_angle_b)
        self.text_size_input.Bind(wx.EVT_SLIDER, self.set_text_size)
        self.node_set_button.Bind(wx.EVT_BUTTON, self.b_node_set_button)

    # LAYER TREE FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def add_tree_nodes(self, parent_item, items):
        i = 0
        for item in items:
            if type(item) == str:
                new_item = self.tree.AppendItem(parent_item, item)
                self.tree.SetItemPyData(new_item, i)
                i += 1
            else:
                pass

    def add_new_tree_nodes(self, parent_item, item, data):
        new_item = self.tree.AppendItem(parent_item, item, ct_type=1)
        new_item.Check(checked=True)
        self.tree.SetItemPyData(new_item, data)

    def get_item_text(self, item):
        if item:
            return self.tree.GetItemText(item)
        else:
            return

    def on_activated(self, event):
        # FIRST CHECK IS FAULT PICKING MODE IS ON, IF IT IS, THEN TURN IT OFF
        if self.fault_picking_switch is True:
            self.fault_picking_switch = False

        # SET OBJECTS WITH THE CHOSEN LAYER
        self.currently_active_layer_id = self.tree.GetPyData(event.GetItem())
        self.nextdens = self.densities[self.currently_active_layer_id]
        self.density_input.SetValue(0.001 * self.densities[self.currently_active_layer_id])
        self.ref_density_input.SetValue(0.001 * self.reference_densities[self.currently_active_layer_id])
        self.susceptibility_input.SetValue(self.susceptibilities[self.currently_active_layer_id])
        self.angle_a_input.SetValue(self.angle_a[self.currently_active_layer_id])
        self.angle_b_input.SetValue(self.angle_b[self.currently_active_layer_id])
        self.current_x_nodes= self.plotx_list[self.currently_active_layer_id]
        self.current_y_nodes = self.ploty_list[self.currently_active_layer_id]
        self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])
        self.update_layer_data()
        self.run_algorithms()

    def on_begin_edit_label(self, event):
        # self.DoGetBestSize()
        pass

    def on_end_edit_label(self, event):
        self.currently_active_layer_id = self.tree.GetPyData(event.GetItem())
        new_label = self.get_item_text(event.GetItem())
        self.tree_items[self.currently_active_layer_id] = str(new_label)

    def delete_all_children(self, event):
        self.tree.DeleteChildren(event)

    def item_checked(self, event):
        """TOGGLE WHETHER OR NOT A LAYER IS INCLUDED IN THE CALCULATIONS"""

        layer = self.tree.GetPyData(event.GetItem())
        if self.layers_calculation_switch[layer] == 1:
            self.layers_calculation_switch[layer] = 0
        else:
            self.layers_calculation_switch[layer] = 1

    def display_info(self):
        self.statusbar.SetStatusText("                                                                                 "
                                     "                                                   "
                                     " || Currently Editing Layer: %s  || "
                                     " || Model Aspect Ratio = %s:1.0  || GRAV RMS = %s "
                                     " || MAG RMS = %s  ||" % (self.currently_active_layer_id, self.model_frame.get_aspect(), self.grav_rms_value,
                                                               self.mag_rms_value), 2)
        self.statusbar.Update()

    # FIGURE DISPLAY FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def zoom(self, event):
        if self.zoom_on is True:
            self.zoom_on = False
            self.nodes = True
            self.nav_toolbar.zoom()
        else:
            self.zoom_on = True
            self.nodes = False
            self.nav_toolbar.zoom()
        # REDRAW
        self.draw()

    def zoom_out(self, event):
        self.nav_toolbar.back()
        # REDRAW
        self.draw()

    def full_extent(self, event):
        """REDRAW MODEL FRAME WITH FULL EXTENT"""
        self.full_extent_adjustment()
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def full_extent_adjustment(self):
        """FIND WHICH FRAME IS REFERENCED & CHANGE SWITCH"""
        if not self.topo_frame.get_visible():
            self.topo_frame.set_visible(True)
        if not self.gravity_frame.get_visible():
            self.gravity_frame.set_visible(True)
        if not self.magnetic_frame.get_visible():
            self.magnetic_frame.set_visible(True)

        'SET CANVAS LIMITS'
        self.model_frame.set_xlim(self.x1, self.x2)
        self.model_frame.set_ylim(self.z2, self.z1)
        self.model_aspect = 1
        if self.topo_frame is not None:
            self.topo_frame.set_xlim(self.model_frame.get_xlim())
        if self.gravity_frame is not None:
            self.gravity_frame.set_xlim(self.model_frame.get_xlim())
        if self.magnetic_frame is not None:
            self.magnetic_frame.set_xlim(self.model_frame.get_xlim())
        self.fig.subplots_adjust(top=0.99, left=-0.045, right=0.99, bottom=0.02, hspace=1.5)

        'ADJUST FRAME SIZING AND SET PROGRAM WINDOW'
        if self.topo_frame is False:
            self.topo_frame.set_visible(False)
        if self.gravity_frame is False:
            self.gravity_frame.set_visible(False)
        if self.magnetic_frame is False:
            self.magnetic_frame.set_visible(False)

    def pan(self, event):
        """PAN MODEL VIEW USING MOUSE DRAG"""
        if self.nodes:
            self.nodes = False
            self.nav_toolbar.pan()
        else:
            self.nodes = True
            self.nav_toolbar.pan()
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def calc_grav_switch(self, event):
        """PREDICTED ANOMALY CALCULATION SWITCH: ON/OFF (SPEEDS UP PROGRAM WHEN OFF)"""
        if self.calc_grav_switch is True:
            self.calc_grav_switch = False
            if self.grav_residuals != [] and self.obs_gravity_data_for_rms != []:
                self.grav_residuals[:, 1] = np.zeros(len(self.obs_gravity_data_for_rms[:, 0]))
                self.grav_rms_plot.set_data(self.grav_residuals[:, 0], self.grav_residuals[:, 1])
        else:
            self.calc_grav_switch = True
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def calc_mag_switch(self, event):
        """PREDICTED ANOMALY CALCULATION SWITCH: ON/OFF (SPEEDS UP PROGRAM WHEN OFF)"""
        if self.calc_mag_switch is True:
            self.calc_mag_switch = False
            if self.mag_residuals != [] and self.obs_mag_data_for_rms != []:
                self.mag_residuals[:, 1] = np.zeros(len(self.obs_mag_data_for_rms[:, 0]))
                self.mag_rms_plot.set_data(self.mag_residuals[:, 0], self.mag_residuals[:, 1])
        else:
            self.calc_mag_switch = True
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    # MODEL WINDOW GRAPHICS OPTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def aspect_increase(self, event):
        if self.model_aspect >= 1:
            self.model_aspect = self.model_aspect + 1
            self.update_layer_data()
            self.draw()
        elif 1.0 > self.model_aspect >= 0.1:
            self.model_aspect = self.model_aspect + 0.1
            self.update_layer_data()
            self.draw()
        else:
            pass

    def aspect_decrease(self, event):
        if self.model_aspect >= 2:
            self.model_aspect = self.model_aspect - 1
            self.update_layer_data()
            self.draw()
        elif 1.0 >= self.model_aspect >= 0.2:
            self.model_aspect = self.model_aspect - 0.1
            self.update_layer_data()
            self.draw()
        else:
            pass

    def aspect_increase2(self, event):
        self.model_aspect = self.model_aspect + 2
        self.update_layer_data()
        self.draw()

    def aspect_decrease2(self, event):
        if self.model_aspect >= 3:
            self.model_aspect = self.model_aspect - 2
            self.update_layer_data()
            self.draw()
        else:
            pass

    def transparency_increase(self, event):
        self.layer_transparency = self.layer_transparency + 0.1
        for i in range(0, self.total_layer_count):
            self.polygon_fills[i][0].set_alpha(self.layer_transparency)
        self.draw()

    def transparency_decrease(self, event):
        if self.layer_transparency >= 0.1:
            self.layer_transparency = self.layer_transparency - 0.1
            for i in range(0, self.total_layer_count):
                self.polygon_fills[i][0].set_alpha(self.layer_transparency)
            self.draw()
        else:
            pass

    # SAVE/LOAD MODEL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def save_model(self, event):
        """
        SAVE MODEL TO DISC IN .Pickle FORMAT
        TO ADD NEW OBJECTS ADD THE OBJECT NAME TO BOTH THE header AND model_params.
        THEN MODIFY THE load_model FUNCTION TO ALSO INLCUDE THE NEW ITEMS
        """
        save_file_dialog = wx.FileDialog(self, "Save model file", "", "", "Model files (*.model)|*.model", wx.FD_SAVE
                                         | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # USER CHANGED THEIR MIND

        # CREATE SAVE DICTIONARY
        self.save_dict = {}

        header = ['model_aspect',  'area', 'xp', 'zp', 'tree_items',
                  'layer_colors', 'plotx_list', 'ploty_list', 'densities', 'reference_densities',
                  'boundary_lock_list', 'layer_lock_list',
                  'earth_field', 'model_azimuth',
                  'mag_observation_elv', 'susceptibilities', 'angle_a', 'angle_b',
                  'absolute_densities',
                  'obs_gravity_data_for_rms', 'obs_mag_data_for_rms',
                  'faults', 'fault_names_list', 'fault_counter', 'current_fault_index',
                  'fault_x_coords_list', 'fault_y_coords_list', 'fault_tree_items', 'fault_counter',
                  'observed_xy_data_list',
                  'observed_gravity_list',
                  'observed_magnetic_list',
                  'observed_topography_list',
                  'well_data_list',
                  'outcrop_data_list',
                  'segy_data_list']

        model_params = [self.model_aspect, self.area, self.xp, self.zp, self.tree_items,
                        self.layer_colors, self.plotx_list, self.ploty_list, self.densities, self.reference_densities,
                        self.boundary_lock_list, self.layer_lock_list,
                        self.earth_field, self.model_azimuth,
                        self.mag_observation_elv, self.susceptibilities, self.angle_a, self.angle_b,
                        self.absolute_densities,
                        self.obs_gravity_data_for_rms, self.obs_mag_data_for_rms,
                        self.faults, self.fault_names_list, self.fault_counter, self.current_fault_index,
                        self.fault_x_coords_list, self.fault_y_coords_list, self.fault_tree_items, self.fault_counter,
                        self.observed_xy_data_list,
                        self.observed_gravity_list,
                        self.observed_magnetic_list,
                        self.observed_topography_list,
                        self.well_data_list,
                        self.outcrop_data_list,
                        self.segy_data_list]

        for i in range(0, len(model_params)):
            try:
                self.save_dict[header[i]] = model_params[i]
            except IOError:
                print(header[i])
        try:
            output_stream = save_file_dialog.GetPath()
            out = open(output_stream, 'w')
            Pickle.dump(self.save_dict, out)
            self.model_saved = True
            self.update_layer_data()
            MessageDialog(self, -1, "Model saved successfully", "Save")
        except IOError:
            MessageDialog(self, -1, "Error in save process.\nModel not saved", "Save")

    def load_model(self, event):
        """LOAD MODEL FROM DISC IN .Pickle FORMAT"""
        try:
            if not self.model_saved:
                if wx.MessageBox("Current content has not been saved! Proceed?", "Please confirm",
                                 wx.ICON_QUESTION | wx.YES_NO, self) == wx.NO:
                    return
                open_file_dialog = wx.FileDialog(self, "Open model file", "", "", "Model files (*.model)|*.model",
                                                 wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
                if open_file_dialog.ShowModal() == wx.ID_CANCEL:
                    return  # USER CHANGED THEIR MIND

            # INITALISE MODEL PARAMETERS
            self.initalise_model()
            self.model_aspect = 1.

            # OPEN DATA STREAM
            file_in = open_file_dialog.GetPath()
            model_data = Pickle.load(open(file_in, "r+"))

            # CLEAR MEMORY
            gc.collect()
            del gc.garbage[:]

            # LOAD DATA INTO MODEL
            for x in range(len(model_data)):
                setattr(self, model_data.keys()[x], model_data.values()[x])

            # SAVE LOADED TREE ITEMS (WILL BE REMOVED BY self.start)
            self.loaded_tree_items = self.tree_items
            self.loaded_fault_tree_items = self.fault_tree_items

            # DRAW CANVAS
            self.start(self.area, self.xp, self.zp)

            # SET LAYERS self.boundary_lock_status
            self.boundary_lock_status = [[]]
            for x in range(len(self.boundary_lock_list)):
                if self.boundary_lock_list[x] == 0:
                    self.boundary_lock_status[x] = ['locked']
                    self.boundary_lock_status.append([])
                else:
                    self.boundary_lock_status[x] = ['unlocked']
                    self.boundary_lock_status.append([])

            # SET LAYERS self.layer_lock_status
            self.layer_lock_status = [[]]
            for x in range(len(self.layer_lock_list)):
                if self.layer_lock_list[x] == 0:
                    self.layer_lock_status[x] = ['locked']
                    self.layer_lock_status.append([])
                else:
                    self.layer_lock_status[x] = ['unlocked']
                    self.layer_lock_status.append([])

            # ----------------------------------------------------------------------------------------------------------
            # LOAD OBSERVED TOPOGRAPHY DATA
            if len(self.observed_topography_list) > 0:
                self.replot_observed_topography_data()

            # LOAD OBSERVED GRAVITY DATA
            if len(self.observed_gravity_list) > 0:
                self.replot_observed_gravity_data()

            # LOAD OBSERVED MAGNETIC DATA
            if len(self.observed_magnetic_list) > 0:
                self.replot_observed_magnetic_data()
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # Set VARIABLE VALUES FROM LOADED DATA
            self.currently_active_layer_id = (len(self.densities) - 1)
            self.total_layer_count = (len(self.densities) - 1)
            self.density_input.SetValue(0.001 * self.densities[self.currently_active_layer_id])
            self.ref_density_input.SetValue(0.001 * self.reference_densities[self.currently_active_layer_id])
            self.current_x_nodes= self.plotx_list[self.currently_active_layer_id]
            self.current_y_nodes = self.ploty_list[self.currently_active_layer_id]

            self.gravity_polygons = []
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            if self.layer_colors == [[]]:
                self.layer_colors = ['black' for _ in range(self.currently_active_layer_id + 1)]
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # LOAD LAYER TREE ITEMS
            self.tree.DeleteAllItems()  # DELETE CURRENT LAYER TREE
            self.root = self.tree.AddRoot("Layers:")  # CREATE NEW TREE
            self.tree.SetItemPyData(self.root, None)

            for i in range(1, len(self.loaded_tree_items)):
                tree_item = self.tree.AppendItem(self.root, "%s" % self.loaded_tree_items[i], ct_type=1)
                tree_item.Check(checked=True)
                self.tree.SetItemPyData(tree_item, i)

            self.layers_calculation_switch = [1] * len(self.loaded_tree_items)

            self.tree_items = self.loaded_tree_items
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # LOAD FAULT TREE ITEMS
            for i in range(0, len(self.loaded_fault_tree_items)):
                fault_tree_item = self.fault_tree.AppendItem(self.fault_tree_root, "%s" %
                                                             self.loaded_fault_tree_items[i], ct_type=1)
                fault_tree_item.Check(checked=True)
                self.fault_tree.SetItemPyData(fault_tree_item, i)

            self.fault_tree_items = self.loaded_fault_tree_items
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # MAKE LAYER LINES AND POLYGONS
            self.layer_lines = [[] for _ in range(self.currently_active_layer_id + 1)]
            self.polygon_fills = [[] for _ in range(self.currently_active_layer_id + 1)]
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # LOAD LAYERS
            self.load_layer_data()

            # SET CURRENT NODE AS A OFF STAGE (PLACE HOLDER)
            self.current_node = self.model_frame.scatter(-40000., 0., marker='o', color='r', zorder=10)
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            for i in range(0, self.fault_counter):
                # DRAW FAULTS
                self.faults[i] = self.model_frame.plot(self.fault_x_coords_list[i], self.fault_y_coords_list[i],
                                                   color='k', linewidth=0.5, zorder=1, marker='s',
                                                   alpha=1.0)
            # CREATE NEW CURRENT FAULT GRAPHIC
            self.faultline, = self.model_frame.plot([-100000, -100000], [-100000, -100000], marker='s', color='m',
                                                linewidth=0.75, alpha=1.0, zorder=2, picker=True)

            # Set Fault PICKING SWITCH OFF (DEFAULT TO LAYER MODE)
            self.fault_picking_swtich = False
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # LOAD OBSERVED WELL DATA
            if len(self.well_data_list) > 0:
                self.replot_well_data()

            # LOAD OBSERVED XY DATA
            if len(self.observed_xy_data_list) > 0:
                self.replot_observed_xy_data()

            # LOAD OUTCROP DATA
            if len(self.outcrop_data_list) > 0:
                self.replot_outcrop_data()

            # LOAD SEGY DATA
            if len(self.segy_data_list) > 0:
                self.replot_segy_data()
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # REFRESH SIZER POSITIONS
            self.Hide()
            self.Show()
            # UPDATE LAYER DATA AND PLOT
            self.draw()
            self.update_layer_data()
            self.run_algorithms()
            self.draw()
            self.Restore()  # FIX'S DISPLAY ISSUE
            # ----------------------------------------------------------------------------------------------------------

        # LOAD ERRORS
        except IOError:
            error_message = "IO ERROR IN LOADING PROCESS - MODEL NOT LOADED"
            MessageDialog(self, -1, error_message, "Load Error")
            raise
        except IndexError:
            error_message = "INDEX ERROR IN LOADING PROCESS - MODEL NOT LOADED"
            MessageDialog(self, -1, error_message, "Load Error")
            raise

        # MAXIMIZE FRAME
        self.Maximize(True)

    # MODEL LOADING PLOTTING FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def replot_observed_xy_data(self):
        """ADD LOADED OBSERVED XY DATA TO THE MODEL FRAME"""
        for x in range(len(self.observed_xy_data_list)):
            if self.observed_xy_data_list[x] is not None:
                # DRAW DATA IN MODEL FRAME
                self.self.observed_xy_data_list[x].mpl_actor = self.topo_frame.scatter(
                    self.self.observed_xy_data_list[x].data[:, 0], self.observed_xy_data_list[x].data[:, 1],
                    marker='o', color=self.observed_xy_data_list[x].color, s=5,
                    gid=self.observed_xy_data_list[x].id)

                # APPEND NEW DATA MENU TO 'XY data MENU'
                self.obs_submenu = wx.Menu()
                self.m_xy_submenu.Append(4000 + self.observed_xy_data_list[x].id, self.observed_xy_data_list[x].name,
                                                                                 self.obs_submenu)
                # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
                self.obs_submenu.Append(4000 + self.observed_xy_data_list[x].id, 'delete observed data')

                # BIND TO DEL XY FUNC
                self.Bind(wx.EVT_MENU, self.delete_xy, id=4000 + self.observed_xy_data_list[x].id)

        # SET GRAVITY COUNTER
        self.observed_xy_data_counter = len(self.observed_xy_data_list)

    def replot_observed_topography_data(self):
        """ADD LOADED OBSERVED TOPOGRAPHY TO THE MODEL FRAME"""
        for x in range(len(self.observed_topography_list)):
            if self.observed_topography_list[x] is not None:
                # DRAW DATA IN MODEL FRAME
                self.observed_topography_list[x].mpl_actor = self.topo_frame.scatter(
                    self.observed_topography_list[x].data[:, 0], self.observed_topography_list[x].data[:, 1],
                    marker='o', color=self.observed_topography_list[x].color, s=5,
                    gid=self.observed_topography_list[x].id)

                # ADD OBJECT TO MENUVAR
                self.obs_submenu = wx.Menu()
                self.m_topo_submenu.Append(10000 + self.observed_topography_list[x].id,
                                            self.observed_topography_list[x].name,
                                            self.obs_submenu)
                self.obs_submenu.Append(10000 + self.observed_topography_list[x].id, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_observed_topography, id=10000 + self.observed_topography_list[x].id)

                # TURN ON OBSERVED GRAVITY SWITCH
                self.observed_topography_switch = True

        # SET GRAVITY COUNTER
        self.observed_topography_counter = len(self.observed_topography_list)

    def replot_observed_gravity_data(self):
        """ADD LOADED OBSERVED GRAVITY TO THE MODEL FRAME"""
        for x in range(len(self.observed_gravity_list)):
            if self.observed_gravity_list[x] is not None:
                # DRAW DATA IN MODEL FRAME
                self.observed_gravity_list[x].mpl_actor = self.gravity_frame.scatter(
                    self.observed_gravity_list[x].data[:, 0], self.observed_gravity_list[x].data[:, 1], marker='o',
                    color=self.observed_gravity_list[x].color, s=5, gid=self.observed_gravity_list[x].id)

                # ADD OBJECT TO MENUVAR
                self.obs_submenu = wx.Menu()
                self.m_obs_g_submenu.Append(11000+self.observed_gravity_list[x].id, self.observed_gravity_list[x].name,
                                            self.obs_submenu)
                self.obs_submenu.Append(11000+self.observed_gravity_list[x].id, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000+self.observed_gravity_list[x].id)

                # TURN ON OBSERVED GRAVITY SWITCH
                self.observed_gravity_switch = True

        # SET GRAVITY COUNTER
        self.observed_gravity_counter = len(self.observed_gravity_list)

    def replot_observed_magnetic_data(self):
        """ADD LOADED OBSERVED MAGNETIC DATA TO THE MODEL FRAME"""
        for x in range(len(self.observed_magnetic_list)):
            if self.observed_magnetic_list[x] is not None:
                # DRAW DATA IN MODEL FRAME
                self.observed_magnetic_list[x].mpl_actor = self.magnetic_frame.scatter(
                                                                          self.observed_magnetic_list[x].data[:, 0],
                                                                          self.observed_magnetic_list[x].data[:, 1],
                                                                          marker='o',
                                                                          color=self.observed_magnetic_list[x].color,
                                                                          s=5, gid=self.observed_magnetic_list[x].id)

                # ADD OBJECT TO MENUVAR
                self.mag_submenu = wx.Menu()
                self.m_obs_mag_submenu.Append(12000+self.observed_magnetic_list[x].id,
                                              self.observed_magnetic_list[x].name,
                                              self.mag_submenu)
                self.mag_submenu.Append(12000+self.observed_magnetic_list[x].id, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12000+self.observed_magnetic_list[x].id)

                # TURN ON OBSERVED GRAVITY SWITCH
                self.observed_magnetic_switch = True

        # SET MAGNETIC COUNTER
        self.observed_magnetic_counter = len(self.observed_magnetic_list)

    def replot_well_data(self):
        """ADD LOADED WELL DATA TO THE MODEL FRAME"""
        for x in range(len(self.well_data_list)):
            if self.well_data_list[x] is not None:

                # SET CURRENT WELL
                self.loaded_well = self.well_data_list[x]
                well = self.well_data_list[x]

                # CREATE FILE MENU DATA
                self.well_name_submenu = wx.Menu()
                self.m_wells_submenu.Append(well.id + 3000, well.name, self.well_name_submenu)
                self.well_name_submenu.Append(well.id + 2000, 'Hide/Show')
                self.well_name_submenu.Append(well.id + 3000, 'Delete well')
                self.Bind(wx.EVT_MENU, self.show_hide_well, id=well.id + 2000)
                self.Bind(wx.EVT_MENU, self.delete_well, id=well.id + 3000)

                # DRAW WELL IN MODEL FRAME
                y1 = well.data[0][1].astype(float)
                y2 = well.data[-1][-1].astype(float)
                well_x_location = well.data[1][1].astype(float)
                wellx = (well_x_location, well_x_location)
                welly = (y1, y2)
                well.mpl_actor = self.model_frame.plot(wellx, welly, linestyle='-', linewidth='2', color='black')

                # PLOT WELL NAME
                well.mpl_actor_name = self.model_frame.annotate(well.name, xy=(well_x_location, -0.5),
                                                                xytext=(well_x_location, -0.5),
                                                                fontsize=well.textsize, weight='bold',
                                                                horizontalalignment='center', color='black',
                                                                bbox=dict(boxstyle="round,pad=.2", fc="0.8"),
                                                                clip_on=True)

                # PLOT WELL HORIZONS
                # SET EMPTY ARRAYS TO FILL WITH LABELS AND HORIZONS
                well.labels_list = [None] * (len(well.data) - 2)
                well.horizons_list = [None] * (len(well.data) - 2)
                for i in range(2, len(well.data)):
                    y = [well.data[i][1].astype(float), well.data[i][1].astype(float)]
                    x = [well.data[1][1].astype(float) - 1, well.data[1][1].astype(float) + 1]

                    # PLOT HORIZON LINE
                    well.horizons_list[i-2] = self.model_frame.plot(x, y, linestyle='-', linewidth='2', color='black')
                    horizon_y_pos = well.data[i][1].astype(float)
                    horizon = well.data[i][0].astype(str)

                    # ALTERNATE POSITION OF ODDs/EVENs TO TRY AND AVOID OVERLAP
                    if i % 2 == 0:
                        horizon_x_pos = well.data[1][1].astype(float) - 1.05
                        well.labels_list[i-2] = self.model_frame.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                        xytext=(horizon_x_pos, horizon_y_pos),
                                                                        fontsize=well.text_size, weight='bold',
                                                                        horizontalalignment='left',
                                                                        verticalalignment='top',
                                                                        color='black',
                                                                        bbox=dict(boxstyle="round,pad=.4",
                                                                                  fc="0.8", ec='None'), clip_on=True)
                    else:
                        horizon_x_pos = well.data[1][1].astype(float) + 1.05
                        well.labels_list[i-2] = self.model_frame.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                        xytext=(horizon_x_pos, horizon_y_pos),
                                                                        fontsize=well.text_size, weight='bold',
                                                                        horizontalalignment='right',
                                                                        verticalalignment='top',
                                                                        color='black',
                                                                        bbox=dict(boxstyle="round,pad=.4",
                                                                                  fc="0.8", ec='None'), clip_on=True)

        # SET WELL COUNTER
        self.well_counter = len(self.well_data_list)

    def replot_outcrop_data(self):
        """ADD LOADED OUTCROP DATA TO THE MODEL FRAME"""
        for x in range(len(self.outcrop_data_list)):
            if self.outcrop_data_list[x] is not None:

                # SET CURRENT OUTCROP DATA RECORD
                outcrop = self.outcrop_data_list[x]

                # PLOT MARKERS IN MODEL
                outcrop.lines = [None] * len(outcrop.data)
                for i in range(len(outcrop.data)):
                    x1 = outcrop.data[i, 0].astype(float)
                    y1 = outcrop.data[i, 1].astype(float)
                    y2 = outcrop.data[i, 2].astype(float)
                    x = (x1, x1)
                    y = (y1, y2)
                    outcrop.lines[i] = self.model_frame.plot(x, y, linestyle='-', linewidth='2', color=outcrop.color)

                # DRAW TEXT LABELS
                outcrop.labels = [None] * len(outcrop.data)

                # CREATE TEXT XYT
                text = zip(outcrop.data[:, 0].astype(float), outcrop.data[:, 1].astype(float),
                           outcrop.data[:, 3].astype(str))

                for i in range(len(outcrop.data)):

                    # ALTERNATE POSITION OF ODDs/EVENs To TRY AND AVOID OVERLAP
                    if i % 2 == 0:
                        outcrop.labels[i] = self.model_frame.annotate(text[i][2], xy=(text[i][0], text[i][1]),
                                                                      xytext=(text[i][0], text[i][1]),
                                                                      fontsize=outcrop.textsize,
                                                                      weight='regular', horizontalalignment='right',
                                                                      verticalalignment='bottom',
                                                                      color='black',
                                                                      bbox=dict(boxstyle="round,pad=.4", fc="0.8",
                                                                                ec='None'), clip_on=True)
                    else:
                        outcrop.labels[i] = self.model_frame.annotate(text[i][2], xy=(text[i][0], text[i][1]),
                                                                      xytext=(text[i][0], text[i][1]),
                                                                      fontsize=outcrop.textsize,
                                                                      weight='regular', horizontalalignment='left',
                                                                      verticalalignment='top',
                                                                      color='black',
                                                                      bbox=dict(boxstyle="round,pad=.4", fc="0.8",
                                                                                ec='None'), clip_on=True)

                #  APPEND NEW DATA MENU TO 'OUTCROP DATA MENU'
                self.outcrop_submenu = wx.Menu()
                self.m_outcrop_submenu.Append(self.outcrop_data_count, outcrop.name, self.outcrop_submenu)

                # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
                self.outcrop_submenu.Append(13000 + self.outcrop_data_count, 'delete observed data')

                # BIND TO DEL XY FUNC
                self.Bind(wx.EVT_MENU, self.delete_outcrop_data, id=13000 + self.outcrop_data_count)

                # ISET OUTCROP COUNTER
                self.outcrop_data_count = len(self.outcrop_data_list)

    def replot_segy_data(self):
        """ PLOT SEGY DATA"""
        for x in range(len(self.segy_data_list)):
            if self.segy_data_list[x] is not None:

                # SET CURRENT SEGY OBJECT
                segy = self.segy_data_list[x]

                # LOAD SEGY DATA
                try:
                    section = read(segy.file, unpack_trace_headers=False)
                except IOError:
                    load_error = MessageDialog(self, -1, "SEGY LOAD ERROR: FILE NOT FOOUND", "Segy load error")
                    return

                nsamples = len(section.traces[0].data)
                ntraces = len(section.traces)
                seis_data = plt.zeros((nsamples, ntraces))
                for i, tr in enumerate(section.traces):
                    seis_data[:, i] = tr.data
                del section

                # SET AXIS
                segy.axis = [segy.dimensions[0], segy.dimensions[1], segy.dimensions[3], segy.dimensions[2]]

                # PLOT SEGY DATA ON MODEL
                segy.mpl_actor = self.model_frame.imshow(seis_data, vmin=self.segy_gain_neg, vmax=self.segy_gain_pos,
                                                         aspect='auto', extent=segy.axis, cmap=self.segy_color_map,
                                                         alpha=0.75)
                # REMOVE SEGY DATA
                del seis_data

                # SET SEGY ON SWTICH
                self.segy_on = True

                # APPEND NEW SEGY TO THE SEGY DATA LIST
                self.segy_data_list.append(segy)

                # ADD SEGY_NAME TO SEGY MENU. NB: ID's START AT 1000
                self.segy_name_submenu = wx.Menu()
                self.m_segy_submenu.Append(segy.id + 1000, segy.name, self.segy_name_submenu)
                self.segy_name_submenu.Append(segy.id + 1000, 'delete segy')
                self.Bind(wx.EVT_MENU, self.remove_segy, id=segy.id + 1000)

            # INCREMENT COUNTER
            self.segy_counter = len(self.segy_data_list)


            # TOPOGRAPHY DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_topo(self, event):
        self.load_window = LoadObservedDataFrame(self, -1, 'Load observed data', 'topography')
        self.load_window.Show(True)

    def open_obs_t(self):
        """
        LOAD OBSERVE TOPOGRAPHY DATA.
        DATA ARE STORED IN gmg.observed_topography_list as a ObservedData object.
        object IDs start at 11000.
        """

        # PARSE USER INPUT FILE
        input_file = self.load_window.file_path

        # CREATE NEW OBSERVED GRAVITY OBJECT
        observed_topography = ObservedData()

        # SET ATTRIBUTES
        observed_topography.id = int(self.observed_topography_counter)
        observed_topography.type = str('observed')
        observed_topography.name = self.load_window.observed_name
        observed_topography.color = self.load_window.color_picked
        observed_topography.data = np.genfromtxt(input_file, delimiter=' ', dtype=float)
        observed_topography.mpl_actor = self.topo_frame.scatter(observed_topography.data[:, 0],
                                                             observed_topography.data[:, 1], marker='o',
                                                             color=observed_topography.color, s=5,
                                                             gid=observed_topography.id)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_topography_list.append(observed_topography)

        # TURN ON OBSERVED GRAVITY SWITCH
        self.observed_topography_switch = True

        # APPEND NEW DATA MENU TO 'TOPO data MENU'
        self.topo_submenu = wx.Menu()
        self.m_topo_submenu.Append(10000+observed_topography.id, observed_topography.name, self.topo_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.topo_submenu.Append(10000+observed_topography.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_observed_topography, id=10000+observed_topography.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_topography_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def delete_observed_topography(self, event):
        """DELETE AN OBSERVED TOPOGRAPHY DATA RECORD"""
        # DESTROY MENUBAR
        self.m_topo_submenu.DestroyItem(event.Id)

        # REMOVE OBJECT AND MPL ACTOR
        obj_id = event.Id-10000
        self.observed_topography_list[obj_id].mpl_actor.set_visible(False)
        self.observed_topography_list[obj_id] = None

        # UPDATE MODEL
        self.update_layer_data()
        self.draw()

    # GRAVITY DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_obs_g(self, event):
        self.load_window = LoadObservedDataFrame(self, -1, 'Load observed data', 'gravity')
        self.load_window.Show(True)

    def open_obs_g(self):
        """
        LOAD OBSERVE GRAVITY DATA.
        DATA ARE STORED IN gmg.observed_gravity_list as a ObservedData object.
        object IDs start at 11000.
        """

        # PARSE USER INPUT FILE
        input_file = self.load_window.file_path

        # CREATE NEW OBSERVED GRAVITY OBJECT
        observed_gravity = ObservedData()

        # SET ATTRIBUTES
        observed_gravity.id = int(self.observed_gravity_counter)
        observed_gravity.type = str('observed')
        observed_gravity.name = self.load_window.observed_name
        observed_gravity.color = self.load_window.color_picked
        observed_gravity.data = np.genfromtxt(input_file, delimiter=' ', dtype=float)
        observed_gravity.mpl_actor = self.gravity_frame.scatter(observed_gravity.data[:, 0],
                                                                observed_gravity.data[:, 1], marker='o',
                                                                color=observed_gravity.color, s=5,
                                                                gid=observed_gravity.id)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_gravity_list.append(observed_gravity)

        # TURN ON OBSERVED GRAVITY SWITCH
        self.observed_gravity_switch = True

        # APPEND NEW DATA MENU TO 'GRAV data MENU'
        self.grav_submenu = wx.Menu()
        self.m_obs_g_submenu.Append(11000+observed_gravity.id, observed_gravity.name, self.grav_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.grav_submenu.Append(11000+observed_gravity.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000+observed_gravity.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_gravity_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def delete_obs_grav(self, event):
        # DESTROY MENUBAR
        self.m_obs_g_submenu.DestroyItem(event.Id)

        # REMOVE OBJECT AND MPL ACTOR
        obj_id = event.Id-11000
        self.observed_gravity_list[obj_id].mpl_actor.set_visible(False)
        self.observed_gravity_list[obj_id] = None

        # UPDATE MODEL
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def save_modelled_grav(self, event):
        """SAVE PREDICTED GRAVITY TO EXTERNAL ASCII FILE"""
        save_file_dialog = wx.FileDialog(self, "Save Predicted Anomaly", "", "", "Predicted Anomaly (*.txt)|*.txt",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THEIR MIND
        # SAVE TO DISC
        outputfile = save_file_dialog.GetPath()
        np.savetxt(outputfile, zip((self.xp * 0.001), self.predgz), delimiter=' ', fmt='%.6f %.6f')

    # MAGNETIC DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_obs_m(self, event):
        self.load_window = LoadObservedDataFrame(self, -1, 'Load observed data', 'magnetics')
        self.load_window.Show(True)

    def open_obs_m(self):
        """
        LOAD OBSERVE GRAVITY DATA.
        DATA ARE STORED IN gmg.observed_gravity_list as a ObservedData object.
        object IDs start at 12000.
        """

        # PARSE USER INPUT FILE
        input_file = self.load_window.file_path

        # CREATE NEW OBSERVED GRAVITY OBJECT
        observed_magnetic = ObservedData()

        # SET ATTRIBUTES
        observed_magnetic.id = int(self.observed_magnetic_counter)
        observed_magnetic.type = str('observed')
        observed_magnetic.name = self.load_window.observed_name
        observed_magnetic.color = self.load_window.color_picked
        observed_magnetic.data = np.genfromtxt(input_file, delimiter=' ', dtype=float)
        observed_magnetic.mpl_actor = self.magnetic_frame.scatter(observed_magnetic.data[:, 0],
                                                                  observed_magnetic.data[:, 1], marker='o',
                                                                  color=observed_magnetic.color, s=5,
                                                                  gid=observed_magnetic.id)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_magnetic_list.append(observed_magnetic)

        # TURN ON OBSERVED GRAVITY SWITCH
        self.observed_magnetic_switch = True

        # APPEND NEW DATA MENU TO 'GRAV data MENU'
        self.mag_submenu = wx.Menu()
        self.m_obs_mag_submenu.Append(12000+observed_magnetic.id, observed_magnetic.name, self.mag_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.mag_submenu.Append(12000+observed_magnetic.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12000+observed_magnetic.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_magnetic_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def delete_obs_mag(self, event):
        # DESTROY MENUBAR
        self.m_obs_mag_submenu.DestroyItem(event.Id)

        # REMOVE OBJECT AND MPL ACTOR
        obj_id = event.Id-12000
        self.observed_magnetic_list[obj_id].mpl_actor.set_visible(False)
        self.observed_magnetic_list[obj_id] = None

        # UPDATE MODEL
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def set_mag_variables(self, event):
        mag_box = MagDialog(self, -1, 'Segy Dimensions', self.area)
        answer = mag_box.ShowModal()
        self.mag_observation_elv = mag_box.mag_observation_elv * 1000.  # CONVERT FROM (km) TO (m)
        self.model_azimuth = mag_box.model_azimuth
        self.earth_field = mag_box.earth_field
        self.draw()

    def save_modelled_mag(self, event):
        save_file_dialog = wx.FileDialog(self, "Save Predicted Anomaly", "", "", "Predicted Anomaly (*.txt)|*.txt",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THEIR MIND

        # SAVE TO DISC
        outputfile = save_file_dialog.GetPath()
        np.savetxt(outputfile, zip((self.xp * 0.001), self.prednt), delimiter=' ', fmt='%.6f %.6f')

    # SEGY DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def segy_input(self, event):
        seismic_data_box = SeisDialog(self, -1, 'Segy Dimensions', self.area)
        answer = seismic_data_box.ShowModal()
        self.d = seismic_data_box.dimensions
        self.segy_name = seismic_data_box.segy_name_input
        self.sx1, self.sx2, self.sz1, self.sz2 = self.d
        self.load_segy(self)

    def load_segy(self, event):
        try:
            open_file_dialog = wx.FileDialog(self, "Open Observed file", "", "", "All files (*.*)|*.*",
                                             wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
            if open_file_dialog.ShowModal() == wx.ID_CANCEL:
                return  # THE USER CHANGED THEIR MIND

            # PARSE INPUT FILE
            file_in = open_file_dialog.GetPath()

            # CREATE NEW SEGY OBJECT
            segy = SegyData()

            # ASSIGN ATTRIBUTES
            segy.id = self.segy_counter
            segy.file = file_in
            segy.name = self.segy_name
            segy.dimensions = self.d

            # LOAD SEGY DATA
            section = read(file_in, unpack_trace_headers=False)
            nsamples = len(section.traces[0].data)
            ntraces = len(section.traces)
            seis_data = plt.zeros((nsamples, ntraces))
            for i, tr in enumerate(section.traces):
                seis_data[:, i] = tr.data
            del section

            # SET AXIS
            segy.axis = [segy.dimensions[0], segy.dimensions[1], segy.dimensions[3], segy.dimensions[2]]

            # PLOT SEGY DATA ON MODEL
            segy.mpl_actor = self.model_frame.imshow(seis_data, vmin=self.segy_gain_neg, vmax=self.segy_gain_pos,
                                                           aspect='auto', extent=segy.axis, cmap=self.segy_color_map,
                                                           alpha=0.75)
            # REMOVE SEGY DATA
            del seis_data

            # SET SEGY ON SWTICH
            self.segy_on = True

            # APPEND NEW SEGY TO THE SEGY DATA LIST
            self.segy_data_list.append(segy)

            # ADD SEGY_NAME TO SEGY MENU. NB: ID's START AT 1000
            self.segy_name_submenu = wx.Menu()
            self.m_segy_submenu.Append(self.segy_counter + 1000, segy.name, self.segy_name_submenu)
            self.segy_name_submenu.Append(self.segy_counter + 1000, 'delete segy')
            self.Bind(wx.EVT_MENU, self.remove_segy, id=self.segy_counter + 1000)

            # INCREMENT COUNTER
            self.segy_counter += 1

        except AttributeEditor:
            load_error = MessageDialog(self, -1, "SEGY LOAD ERROR", "segy load error")

        # REPLOT MODEL
        self.update_layer_data()
        self.draw()

    def remove_segy(self, event):
        """DELETE SEGY DATA. NB 1000 IS TAKEN FROM EVENT.ID TO PREVENT OVERLAP WITH GRAV EVENT.IDS"""
        if self.segy_on:

            # DELETE SEGY OBJECT
            obj_id = event.Id - 1000
            self.segy_data_list[obj_id].mpl_actor.set_visible(False)
            self.segy_data_list[obj_id].mpl_actor.remove()
            self.segy_data_list[obj_id] = None

            # REMOVE MENUBAR
            self.m_segy_submenu.DestroyItem(event.Id)

            # UPDATE MODEL
            self.update_layer_data()
            self.set_frame_limits()
            self.draw()

    def segy_color_adjustment(self, event):
        if event.Id == 901:
            for s in range(0, len(self.segy_data_list)):
                if self.segy_data_list[s] is not None:
                    self.segy_data_list[s].mpl_actor.set_cmap(cm.gray)
        else:
            if event.Id == 902:
                for s in range(0, len(self.segy_data_list)):
                    if self.segy_data_list[s] is not None:
                        self.segy_data_list[s].mpl_actor.set_cmap(cm.seismic)
        # REDRAW MODEL
        self.draw()

    def gain_increase(self, event):
        # CHANGE THE GAIN VALUE
        gain_pos = self.segy_gain_pos - 1.0
        if gain_pos < 1.0:
            return
        else:
            self.segy_gain_pos = gain_pos
            self.segy_gain_neg = -gain_pos
            # REDRAW THE SEGY PML ACTOR
            for s in range(0, len(self.segy_data_list)):
                if self.segy_data_list[s] is not None:
                    self.segy_data_list[s].mpl_actor.set_clim(vmax=self.segy_gain_pos)
                    self.segy_data_list[s].mpl_actor.set_clim(vmin=self.segy_gain_neg)

        # REDRAW MODEL
        self.draw()

    def gain_decrease(self, event):
        # CHANGE THE GAIN VALUE
        gain_pos = self.segy_gain_pos + 1.0
        if gain_pos < 1.0:
            return
        else:
            self.segy_gain_pos = gain_pos
            self.segy_gain_neg = -gain_pos
            # REDRAW THE SEGY PML ACTOR
            for s in range(0, len(self.segy_data_list)):
                if self.segy_data_list[s] is not None:
                    self.segy_data_list[s].mpl_actor.set_clim(vmax=self.segy_gain_pos)
                    self.segy_data_list[s].mpl_actor.set_clim(vmin=self.segy_gain_neg)

        # REDRAW MODEL
        self.draw()

    # XY DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_xy(self, event):
        """ LOAD & PLOT XY DATA E.G. EQ HYPOCENTERS. NB: ID's start at 5000"""
        self.load_window = LoadObservedDataFrame(self, -1, 'Load observed data', 'XY')
        self.load_window.Show(True)

    def open_xy_data(self):
        """
        OPEN THE XY DATA SELECTED BY THE USER USING THE load_xy FUNC
        NB: IDs START AT 5000
        """
        xy_input_file = self.load_window.file_path
        self.xy_name = self.load_window.observed_name
        self.xy_color = self.load_window.color_picked

        # CREATE NEW OBSERVED GRAVITY OBJECT
        new_xy = ObservedData()

        # SET ATTRIBUTES
        new_xy.data = np.genfromtxt(xy_input_file, dtype=float, autostrip=True)
        new_xy.name = self.load_window.observed_name
        new_xy.color = self.load_window.color_picked
        new_xy.type = str('observed')
        new_xy.mpl_actor = self.model_frame.scatter(new_xy.data[:, 0], new_xy.data[:, 1], marker='o',
                                                    color=new_xy.color, s=3, gid=4000 + self.xy_data_counter)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_xy_data_list.append(new_xy)

        # APPEND NEW DATA MENU TO 'XY data MENU'
        self.obs_submenu = wx.Menu()
        self.m_xy_submenu.Append(4000 + self.xy_data_counter, new_xy.name, self.obs_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.obs_submenu.Append(4000 + self.xy_data_counter, 'delete observed data')

        # BIND TO DEL XY FUNC
        self.Bind(wx.EVT_MENU, self.delete_xy, id=4000 + self.xy_data_counter)

        # INCREMENT XY COUNTER
        self.xy_data_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.draw()

    def delete_xy(self, event):
        """"DELETE OBSERVED XY DATA NB: ID's start at 4000"""
        # DESTROY MENUBAR
        self.m_xy_submenu.DestroyItem(event.Id)

        # REMOVE OBJECT AND MPL ACTOR
        obj_id = event.Id - 4000
        self.observed_xy_data_list[obj_id].mpl_actor.set_visible(False)
        self.observed_xy_data_list[obj_id] = None

        # UPDATE MODEL
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    # WELL DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_well(self, event):
        """
        LOAD A WELL RECORD INTO THE MODEL FRAME.
        IDs BEGIN AT 3000.
        HIDE/SHOW TOGGLE IDs BEGIN AT 2000.    
        """

        # CREATE INSTANCE OF LOADING DIALOG BOX
        open_file_dialog = wx.FileDialog(self, "Open Observed file", "", "", "All files (*.*)|*.*",
                                         wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if open_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THEIR MIND
        else:
            well_in = open_file_dialog.GetPath()
            well_name_box = wx.TextEntryDialog(None, 'Please provide a name for the new well record:')
            answer = well_name_box.ShowModal()

        # CREATE A NEW WELL DATA OBJECT
        well = ObservedWellData()

        # SET WELL ATTRIBUTES
        well.id = self.well_counter

        # SET NAME
        well.name = well_name_box.GetValue()
        well.textsize = 2
        # SET DATA
        with open(well_in, 'r') as f:
            input_record = [line.strip().split(' ') for line in f]
        well.raw_record = input_record

        well.data = np.array(well.raw_record[1:])  # CREATE NP ARRAY WITHOUT HEADER INFO

        # CREATE FILE MENU DATA
        self.well_name_submenu = wx.Menu()
        self.m_wells_submenu.Append(self.well_counter + 3000, well.name, self.well_name_submenu)
        self.well_name_submenu.Append(self.well_counter + 2000, 'Hide/Show')
        self.well_name_submenu.Append(self.well_counter + 3000, 'Delete well')
        self.Bind(wx.EVT_MENU, self.show_hide_well, id=well.id + 2000)
        self.Bind(wx.EVT_MENU, self.delete_well, id=well.id + 3000)


        # DRAW WELL IN MODEL FRAME
        y1 = well.data[0][1].astype(float)
        y2 = well.data[-1][-1].astype(float)
        well_x_location = well.data[1][1].astype(float)
        wellx = (well_x_location, well_x_location)
        welly = (y1, y2)
        well.mpl_actor = self.model_frame.plot(wellx, welly, linestyle='-', linewidth='2', color='black')

        # PLOT WELL NAME
        well.mpl_actor_name = self.model_frame.annotate(well.name, xy=(well_x_location, -0.5),
                                                                         xytext=(well_x_location, -0.5),
                                                                         fontsize=well.textsize, weight='bold',
                                                                         horizontalalignment='center', color='black',
                                                                         bbox=dict(boxstyle="round,pad=.2", fc="0.8"),
                                                                         clip_on=True)

        # PLOT WELL HORIZONS
        # SET EMPTY ARRAYS TO FILL WITH LABELS AND HORIZONS
        well.labels_list = [None] * (len(well.data) - 2)
        well.horizons_list = [None] * (len(well.data) - 2)
        for i in range(2, len(well.data)):
            y = [well.data[i][1].astype(float), well.data[i][1].astype(float)]
            x = [well.data[1][1].astype(float) - 1, well.data[1][1].astype(float) + 1]

            # PLOT HORIZON LINE
            well.horizons_list[i-2] = self.model_frame.plot(x, y, linestyle='-', linewidth='2', color='black')
            horizon_y_pos = well.data[i][1].astype(float)
            horizon = well.data[i][0].astype(str)

            # ALTERNATE POSITION OF ODDs/EVENs TO TRY AND AVOID OVERLAP
            if i % 2 == 0:
                horizon_x_pos = well.data[1][1].astype(float) - 1.05
                well.labels_list[i-2] = self.model_frame.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                 xytext=(horizon_x_pos, horizon_y_pos),
                                                                 fontsize=well.text_size, weight='bold',
                                                                 horizontalalignment='left', verticalalignment='top',
                                                                 color='black', bbox=dict(boxstyle="round,pad=.4",
                                                                 fc="0.8", ec='None'), clip_on=True)
            else:
                horizon_x_pos = well.data[1][1].astype(float) + 1.05
                well.labels_list[i-2] = self.model_frame.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                xytext=(horizon_x_pos, horizon_y_pos),
                                                                fontsize=well.text_size, weight='bold',
                                                                 horizontalalignment='right', verticalalignment='top',
                                                                 color='black', bbox=dict(boxstyle="round,pad=.4",
                                                                 fc="0.8", ec='None'), clip_on=True)

        # APPEND WELL TO WELL DATA LIST
        self.well_data_list.append(well)

        # INCREAMENT WELL COUNTER
        self.well_counter += 1

        # UPDATE GMG
        self.update_layer_data()
        self.draw()

    def show_hide_well(self, event):
        print("Hide/Show well")
        print(event.Id)
        id = event.Id - 2000
        if self.well_data_list[id].mpl_actor[0].get_visible():
            # HIDE WELL
            self.well_data_list[id].mpl_actor[0].set_visible(False)
            self.well_data_list[id].mpl_actor_name.set_visible(False)
            # HIDE HORIZONS
            for h in range(len(self.well_data_list[id].horizons_list)):
                if self.well_data_list[id].horizons_list[h] is not None:
                    self.well_data_list[id].horizons_list[h][0].set_visible(False)
            for l in range(len(self.well_data_list[id].labels_list)):
                if self.well_data_list[id].labels_list[l] is not None:
                    self.well_data_list[id].labels_list[l].set_visible(False)

        else:
            # SHOW WELL
            self.well_data_list[id].mpl_actor[0].set_visible(True)
            self.well_data_list[id].mpl_actor_name.set_visible(True)
            # SHOW HORIZONS
            for h in range(len(self.well_data_list[id].horizons_list)):
                if self.well_data_list[id].horizons_list[h] is not None:
                    self.well_data_list[id].horizons_list[h][0].set_visible(True)
            for l in range(len(self.well_data_list[id].labels_list)):
                if self.well_data_list[id].labels_list[l] is not None:
                    self.well_data_list[id].labels_list[l].set_visible(True)

        # REDRAW
        self.draw()

    def delete_well(self, event):
        """"DELETE WELL DATA NB: ID's start at 2500"""

        # SET ID
        obj_id = event.Id - 3000

        # REMOVE OBJECT AND MPL ACTOR
        self.well_data_list[obj_id].mpl_actor[0].set_visible(False)

        # REMOVE HORIZON MPL ACTORS
        for h in range(len(self.well_data_list[obj_id].horizons_list)):
            if self.well_data_list[obj_id].horizons_list[h] is not None:
                self.well_data_list[obj_id].horizons_list[h][0].set_visible(False)

        # REMOVE LABEL MPL ACTORS
        for l in range(len(self.well_data_list[obj_id].labels_list)):
            if self.well_data_list[obj_id].labels_list[l] is not None:
                self.well_data_list[obj_id].labels_list[l].set_visible(False)

        # SET OBJECT AS NONE
        self.well_data_list[obj_id] = None

        # DELETE MENUBAR ENTRY
        self.m_wells_submenu.DestroyItem(event.Id)

        # UPDATE MODEL
        self.draw()

    # GEOLOGY OUTCROP DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_outcrop_data(self, event):
        """ LOAD & PLOT SURFACE CONTACT DATA. E.G. GEOLOGICAL CONTACTS  NB: ID's start at 13000"""
        self.load_window = LoadObservedDataFrame(self, -1, 'Load observed data', 'outcrop')
        self.load_window.Show(True)

    def open_outcrop_data(self):
        """OPEN THE OUTCROP DATA SELECTED BY THE USER USING THE load_outcrop_data FUNC"""

        # CREATE A NEW OUTCROP DATA OBJECT
        outcrop = ObservedOutcropData()

        # LOAD DATA FILE
        outcrop_input_file = self.load_window.file_path
        outcrop.id = self.outcrop_data_count
        outcrop.name = self.load_window.observed_name
        outcrop.color = self.load_window.color_picked
        outcrop.data = np.genfromtxt(outcrop_input_file, autostrip=True, dtype=str, comments='#')

        # PLOT MARKERS IN MODEL
        outcrop.lines = [None] * len(outcrop.data)
        for i in range(len(outcrop.data)):
            x1 = outcrop.data[i, 0].astype(float)
            y1 = outcrop.data[i, 1].astype(float)
            y2 = outcrop.data[i, 2].astype(float)
            x = (x1, x1)
            y = (y1, y2)
            outcrop.lines[i] = self.model_frame.plot(x, y, linestyle='-', linewidth='2', color=outcrop.color)

        # DRAW TEXT LABELS
        outcrop.labels = [None] * len(outcrop.data)

        # CREATE TEXT XYT
        text = zip(outcrop.data[:, 0].astype(float), outcrop.data[:, 1].astype(float),
                   outcrop.data[:, 3].astype(str))

        for i in range(len(outcrop.data)):

            # ALTERNATE POSITION OF ODDs/EVENs To TRY AND AVOID OVERLAP
            if i % 2 == 0:
                outcrop.labels[i] = self.model_frame.annotate(text[i][2], xy=(text[i][0], text[i][1]),
                                                            xytext=(text[i][0], text[i][1]),
                                                            fontsize=outcrop.textsize,
                                                            weight='regular', horizontalalignment='right',
                                                            verticalalignment='bottom',
                                                            color='black',
                                                            bbox=dict(boxstyle="round,pad=.4", fc="0.8", ec='None'),
                                                            clip_on=True)
            else:
                outcrop.labels[i] = self.model_frame.annotate(text[i][2], xy=(text[i][0], text[i][1]),
                                                            xytext=(text[i][0], text[i][1]),
                                                            fontsize=outcrop.textsize,
                                                            weight='regular', horizontalalignment='left',
                                                            verticalalignment='top',
                                                            color='black',
                                                            bbox=dict(boxstyle="round,pad=.4", fc="0.8", ec='None'),
                                                            clip_on=True)

        #  APPEND NEW DATA MENU TO 'OUTCROP DATA MENU'
        self.outcrop_submenu = wx.Menu()
        self.m_outcrop_submenu.Append(13000 + self.outcrop_data_count, outcrop.name, self.outcrop_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.outcrop_submenu.Append(13000 + self.outcrop_data_count, 'delete observed data')

        # BIND TO DEL XY FUNC
        self.Bind(wx.EVT_MENU, self.delete_outcrop_data, id=13000 + self.outcrop_data_count)

        # INCREMENT CONTACT COUNT
        self.outcrop_data_count += 1

        # APPEND NEW OUTCROP DATA OBJECT TO THE OUTCROP DATA LIST
        self.outcrop_data_list.append(outcrop)

        # UPDATE GMG GUI
        self.update_layer_data()
        self.draw()

    def delete_outcrop_data(self, event):
        """"DELETE OUTCROP DATA NB: ID's start at 13000"""

        # SET ID
        obj_id = event.Id - 13000

        # REMOVE LINE MPL ACTORS
        for i in range(len(self.outcrop_data_list[obj_id].lines)):
            self.outcrop_data_list[obj_id].lines[i][0].set_visible(False)

        # REMOVE LABEL MPL ACTORS
        for i in range(len(self.outcrop_data_list[obj_id].labels)):
            self.outcrop_data_list[obj_id].labels[i].set_visible(False)

        # SET OBJECT AS NONE
        self.outcrop_data_list[obj_id] = None

        # DELETE MENUBAR ENTRY
        self.m_outcrop_submenu.DestroyItem(event.Id)

        # UPDATE MODEL
        self.draw()

    # LAYER & NODE CONTROLS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def b_node_set_button(self, event):
        """
        CHECK IF A NODE FROM THE CURRENT LAYER IS SELECTED; IF NOT, THEN SKIP THIS PART AND ONLY UPDATE ATTRIBUTES
        """

        new_x = float(self.x_input.GetValue())
        new_y = float(self.y_input.GetValue())

        if self.currently_active_layer_id == self.node_layer_reference:
            xt = np.array(self.current_x_nodes)
            yt = np.array(self.current_y_nodes)

            if self.boundary_lock_list[self.currently_active_layer_id] == 0 and self.index_node is not None:
                if xt[self.index_node] == 0 and yt[self.index_node] != 0.001:
                    xt[self.index_node] = 0  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = new_y  # REPLACE OLD Y WITH NEW Y
                elif xt[self.index_node] == self.x2 and yt[self.index_node] != 0.001:
                    xt[self.index_node] = self.x2  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = new_y  # REPLACE OLD Y WITH NEW Y
                elif xt[self.index_node] == 0 and yt[self.index_node] == 0.001:
                    xt[self.index_node] = 0  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                elif xt[self.index_node] == self.x2 and yt[self.index_node] == 0.001:
                    xt[self.index_node] = self.x2  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                elif new_y <= 0:
                    xt[self.index_node] = new_x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                else:
                    xt[self.index_node] = new_x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = new_y  # REPLACE OLD Y WITH NEW Y
            elif self.boundary_lock_list[self.currently_active_layer_id] == 1:
                if new_y <= 0:
                    xt[self.index_node] = new_x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                else:
                    xt[self.index_node] = new_x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = new_y  # REPLACE OLD Y WITH NEW Y

            # DEAL WITH PINCHED NODE
            if self.pinch_switch is True:
                for k in range(0, len(self.index_arg2_list)):
                    if self.index_arg2_list[k] is not None:
                        next_x_list = self.plotx_list[k]
                        next_y_list = self.ploty_list[k]  # GET THE NODE LIST OF THE NEXT LAYER
                        next_x_list[self.index_arg2_list[k]] = new_x
                        next_y_list[self.index_arg2_list[k]] = new_y  # REPLACE THE PINCHED NODE WITH THE NEW NODE
                        self.plotx_list[k] = next_x_list
                        self.ploty_list[k] = next_y_list  # OVERWRITE THE NODE LIST WITH UPDATED LIST

            self.current_x_nodes= xt
            self.current_y_nodes = yt
            self.currently_active_layer.set_data(self.current_x_nodes, self.current_y_nodes)

            # COLOR CURRENTLY SELECTED NODE RED
            self.current_node.set_offsets([new_x, new_y])
        else:
            pass

        # UPDATE LAYER DATA
        self.set_density(self)
        self.set_susceptibility(self)
        self.set_angle_a(self)
        self.set_angle_b(self)

        # UPDATE GMG
        self.update_layer_data()
        self.run_algorithms()

    def get_node_under_point(self, event):
        """
        GET THE INDEX VALUE OF THE NODE UNDER POINT, AS LONG AS IT IS WITHIN NODE_CLICK_LIMIT TOLERANCE OF CLICK
        """
        self.node_layer_reference = self.currently_active_layer_id
        xyt = self.currently_active_layer.get_xydata()
        xt = xyt[:, 0]
        yt = xyt[:, 1]
        d = np.sqrt((xt - event.xdata) ** 2 + (yt - event.ydata) ** 2)
        self.index_arg = np.argmin(d)
        if d[self.index_arg] >= self.node_click_limit:
            return None, None
        else:
            # CHECK IF NODE IS A PINCHED POINT, IF YES FIND NODE OF ABOVE OR BELOW LAYER

            # RESET PINCH SWITCH
            self.pinch_switch = False

            # CREATE LIST OF NONES SAME LENGTH AS NUMBER OF LAYERS IN MODEL
            self.index_arg2_list = [None] * (self.total_layer_count + 1)

            for x in range(0, self.total_layer_count + 1):  # LOOP THROUGH ALL LAYERS TO CHECK FOR PINCHED NODES
                if x == self.currently_active_layer_id:  # SKIP CURRENT LAYER
                    continue
                else:
                    # CHECK FOR PINCHED NODES
                    node_list_x, node_list_y = self.plotx_list[x], self.ploty_list[x]
                    for i in range(0, len(node_list_x)):
                        if node_list_x[i] == xt[self.index_arg] and node_list_y[i] == yt[self.index_arg]:
                            # IF ONE OF THE NODES FROM LIST IS EQUAL TO A NODE FROM THE OTHER LAYER
                            # THEN RETURN THE INDEX
                            self.index_arg2_list[x] = i
                            print("self.index_arg2_list =")
                            print self.index_arg2_list
                            self.pinch_switch = True
                        else:
                            continue

            # RETURN NODE INDEX LISTS
            return self.index_arg, self.index_arg2_list

    def get_fault_node_under_point(self, event):
        """GET THE INDEX VALUE OF THE NODE UNDER POINT, AS LONG AS IT IS WITHIN NODE_CLICK_LIMIT TOLERANCE OF CLICK"""
        # GET FAULT NODE XY DATA
        xy_data = self.faultline.get_xydata()
        x = xy_data[:, 0]
        y = xy_data[:, 1]

        # FIND NODE CLOSEST TO EVENT CLICK POINT
        d = np.sqrt((x - event.xdata) ** 2 + (y - event.ydata) ** 2)
        self.index_arg = np.argmin(d)

        # RETURN RESULTING NODE OR NONE
        if d[self.index_arg] >= self.node_click_limit:
            return None
        else:
            return self.index_arg

    def button_press(self, event):
        """WHAT HAPPENS WHEN THE LEFT MOUSE BUTTON IS PRESSED"""

        if event.inaxes is None:
            return  # CLICK IS OUTSIDE MODEL FRAME SO RETURN
        if event.button != 1:
            return

        if self.fault_picking_switch is False and self.capture is False and self.select_new_layer_nodes is False:
            # THEN GMG IS IN LAYER MODE

            # GET THE NODE CLOSEST TO THE CLICK AND ANY PINCHED NODES
            self.index_node, self.index_arg2_list = self.get_node_under_point(event)
            if self.index_node is None:
                return

            xyt = self.currently_active_layer.get_xydata()
            xt, yt = xyt[:, 0], xyt[:, 1]
            self.x_input.SetValue(xt[self.index_node])
            self.y_input.SetValue(yt[self.index_node])

            # IF PINCH == TRUE, THEN PINCH THE NODE TO NEXT NODE
            if self.pinch_switch is True:
                print("nodes pinched")
                # GET THE NODE NUMBER AND LAYER NUMBER AND PLACE THEM IN "PINCH_NODE_LIST"
                if self.pinch_count == 0:
                    self.current_x_nodes= self.plotx_list[self.currently_active_layer_id]
                    self.current_y_nodes = self.ploty_list[self.currently_active_layer_id]
                    x1 = np.array(self.current_x_nodes)
                    y1 = np.array(self.current_y_nodes)
                    self.index_node = self.get_node_under_point(event)
                    self.pinch_node_list[self.pinch_count] = self.index_node[0]
                    self.pinch_node_list[self.pinch_count + 1] = self.currently_active_layer_id
                    self.pinch_count = + 1
                    # SET THE X AND Y OF THE FIRST NODE AS THAT OF THE SECOND NODE
                else:
                    # SET THE SECOND NODE X AND Y
                    self.plotx2 = self.plotx_list[self.currently_active_layer_id]
                    self.ploty2 = self.ploty_list[self.currently_active_layer_id]
                    x2 = np.array(self.plotx2)
                    y2 = np.array(self.ploty2)
                    self.index_node = self.get_node_under_point(event)
                    new_x = x2[int(self.index_node[0])]
                    new_y = y2[int(self.index_node[0])]

                    # SET THE FIRST NODE X AND Y
                    self.plotx1 = self.plotx_list[int(self.pinch_node_list[1])]
                    self.ploty1 = self.ploty_list[int(self.pinch_node_list[1])]
                    x1 = np.array(self.plotx1)
                    y1 = np.array(self.ploty1)

                    # REPLACE THE ORIGINAL NODE WITH THE NEW NODE
                    x1[int(self.pinch_node_list[0])] = new_x  # REPLACE OLD X WITH NEW X
                    y1[int(self.pinch_node_list[0])] = new_y  # REPLACE OLD Y WITH NEW Y
                    self.plotx_list[int(self.pinch_node_list[1])] = x1
                    self.ploty_list[int(self.pinch_node_list[1])] = y1
                    self.pinch_count = 0

            # COLOR CURRENTLY SELECTED NODE RED
            self.current_node.set_offsets([xt[self.index_node], yt[self.index_node]])

        elif self.fault_picking_switch is True and self.select_new_layer_nodes is False:
            # THEN GMG IS IN FAULT MODE

            # GET CURRENT NODE
            self.selected_node = self.get_fault_node_under_point(event)
            if self.selected_node is None:
                return

            # GET CURRENT X AND Y COORDS
            xyt = self.faultline.get_xydata()
            self.xt = xyt[:, 0]
            self.yt = xyt[:, 1]

            # COLOR CURRENTLY SELECTED NODE RED
            self.current_node.set_offsets([self.xt[self.selected_node], self.yt[self.selected_node]])

        elif self.capture is True or self.select_new_layer_nodes is True:
            # COORDINATE CAPTURE MODE OR NEW LAYER CREATION IS ON
            return

    def move(self, event):
        """WHAT HAPPEN WHEN THE LEFT MOUSE BUTTON IS HELD AND THE MOUSE IS MOVED"""

        if self.index_node is None and self.selected_node is None:
            # NO NODE WAS FOUND NEAR THE CLICK
            return
        if event.inaxes is None:
            # CLICK WAS OUTSIDE THE MODEL FRAME
            return
        if event.button != 1:
            return
        if self.pinch_switch is True:
            # PINCH MODE IS ON
            return
        if self.pan_on is True:
            # PAN MODE IS ON
            return
        if self.zoom_on is True:
            # ZOOM MODE IS ON
            return
        if self.select_new_layer_nodes is True:
            # CURRENTLY CREATING A NEW LAYER
            return
        if self.fault_picking_switch is True:
            # GMG IS IN FAULT MODE

            # ASSIGN NEW X AND Y POINTS
            self.new_x = event.xdata  # GET X OF NEW POINT
            self.new_y = event.ydata  # GET Y OF NEW POINT

            # UPDATE NODE ARRAY
            if self.xt[self.selected_node] == self.x1 and self.yt[self.selected_node] != 0.001:
                self.xt[self.selected_node] = self.x1  # REPLACE OLD X WITH NEW X
                self.yt[self.selected_node] = self.new_y  # REPLACE OLD Y WITH NEW Y
            elif self.xt[self.selected_node] == self.x2 and self.yt[self.selected_node] != 0.001:
                self.xt[self.selected_node] = self.x2  # REPLACE OLD X WITH NEW X
                self.yt[self.selected_node] = self.new_y  # REPLACE OLD Y WITH NEW Y
            elif self.xt[self.selected_node] == 0 and self.yt[self.selected_node] == 0.001:
                self.xt[self.selected_node] = 0  # REPLACE OLD X WITH NEW X
                self.yt[self.selected_node] = 0.001  # REPLACE OLD Y WITH NEW Y
            elif self.xt[self.selected_node] == self.x2 and self.yt[self.selected_node] == 0.001:
                self.xt[self.selected_node] = self.x2  # REPLACE OLD X WITH NEW X
                self.yt[self.selected_node] = 0.001  # REPLACE OLD Y WITH NEW Y
            elif self.new_y <= 0:
                self.xt[self.selected_node] = self.new_x  # REPLACE OLD X WITH NEW X
                self.yt[self.selected_node] = 0.001  # REPLACE OLD Y WITH NEW Y
            else:
                self.xt[self.selected_node] = self.new_x  # REPLACE OLD X WITH NEW X
                self.yt[self.selected_node] = self.new_y  # REPLACE OLD Y WITH NEW Y

            # UPDATE THE FAULT LIST RECORDS
            self.fault_x_coords_list[self.current_fault_index] = self.xt
            self.fault_y_coords_list[self.current_fault_index] = self.yt

            # UPDATE THE CURRENT VIEW OF THE FAULT
            self.faultline.set_data(self.xt, self.yt)

            # UPDATE "CURRENT NODE" RED DOT
            if self.xt[self.selected_node] == self.x1:
                self.current_node.set_offsets([self.x1, self.new_y])
            elif self.xt[self.selected_node] == self.x2:
                self.current_node.set_offsets([self.x2, self.new_y])
            else:
                self.current_node.set_offsets([self.new_x, self.new_y])

            self.update_layer_data()  # UPDATE LAYER DATA

        # GMG IS IN LAYER MODE
        if self.fault_picking_switch is False:
            if self.boundary_lock_list[self.currently_active_layer_id] == 0:
                x = event.xdata  # GET X OF NEW POINT
                y = event.ydata  # GET Y OF NEW POINT

                # GET CURRENT X AND Y ARRAYS
                xt = np.array(self.current_x_nodes)
                yt = np.array(self.current_y_nodes)

                # UPDATE NODE
                if xt[self.index_node] == self.x1 and yt[self.index_node] != 0.001:
                    xt[self.index_node] = self.x1  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = y  # REPLACE OLD Y WITH NEW Y
                elif xt[self.index_node] == self.x2 and yt[self.index_node] != 0.001:
                    xt[self.index_node] = self.x2  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = y  # REPLACE OLD Y WITH NEW Y
                elif xt[self.index_node] == 0 and yt[self.index_node] == 0.001:
                    xt[self.index_node] = 0  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                elif xt[self.index_node] == self.x2 and yt[self.index_node] == 0.001:
                    xt[self.index_node] = self.x2  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                elif y <= 0:
                    xt[self.index_node] = x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                else:
                    xt[self.index_node] = x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = y  # REPLACE OLD Y WITH NEW Y
            elif self.boundary_lock_list[self.currently_active_layer_id] == 1:
                x = event.xdata  # GET X OF NEW POINT
                y = event.ydata  # GET Y OF NEW POINT
                xt = np.array(self.current_x_nodes)
                yt = np.array(self.current_y_nodes)
                if y <= 0:
                    xt[self.index_node] = x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                else:
                    xt[self.index_node] = x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = y  # REPLACE OLD Y WITH NEW Y

            # DEAL WITH PINCHED MODE
            if self.pinch_switch is True:
                for k in range(0, len(self.index_arg2_list)):
                    if self.index_arg2_list[k] is not None:
                        # GET THE NODE LIST OF THE NEXT LAYER
                        next_x_list, next_y_list = self.plotx_list[k], self.ploty_list[k]

                        # REPLACE THE PINCHED NODE WITH THE NEW NODE
                        next_x_list[self.index_arg2_list[k]] = x
                        next_y_list[self.index_arg2_list[k]] = y
                        self.plotx_list[k] = next_x_list
                        self.ploty_list[k] = next_y_list  # OVERWRITE THE NODE LIST WITH UPDATED LIST

            self.current_x_nodes= xt
            self.current_y_nodes = yt
            self.currently_active_layer.set_data(self.current_x_nodes, self.current_y_nodes)

            # UPDATE "CURRENT NODE" RED DOT
            if xt[self.index_node] == self.x1:
                self.current_node.set_offsets([self.x1, y])
            elif xt[self.index_node] == self.x2:
                self.current_node.set_offsets([self.x2, y])
            else:
                self.current_node.set_offsets([x, y])

        # UPDATE LAYER DATA
        self.update_layer_data()

    def button_release(self, event):
        """WHAT HAPPENS WHEN THE LEFT MOUSE BUTTON IS RELEASED"""

        if event.inaxes is None:
            # CLICK WAS OUTSIDE THE MODEL FRAME
            return
        if event.button != 1:
            return

        if self.capture is True:
            # GMG IS IN COORDINATE CAPTURE MODE SO ADD THE CURRENT COORDINATES TO THE TABLE
            self.capture_window.table.Append((event.xdata, event.ydata))

        # NEW FLOATING LAYER CREATION SEQUENCE
        if self.select_new_layer_nodes is True:
            # APPEND NEW COORDINATES
            self.new_plotx.append(event.xdata)
            self.new_ploty.append(event.ydata)

            if self.click_count == 0:
                # PLOT NEW NODES
                self.new_layer_nodes = self.model_frame.plot(self.new_plotx, self.new_ploty, color='blue',marker='o')
                # FILL LAYER
                self.new_layer_fill = self.model_frame.fill(self.new_plotx, self.new_ploty, color='blue',
                                                        alpha=self.layer_transparency, closed=True, linewidth=None,
                                                        ec=None)
                # INCREMENT CLICK COUNTER
                self.click_count += 1

            elif self.click_count < 3:
                self.new_layer_nodes[0].set_xdata(self.new_plotx)
                self.new_layer_nodes[0].set_ydata(self.new_ploty)
                self.new_layer_fill[0].set_xy(zip(self.new_plotx, self.new_ploty))

                # INCREMENT CLICK COUNTER
                self.click_count += 1

            else:
                # REMOVE THE TEMP LAYER MPL ACTOR
                self.select_new_layer_nodes = False
                self.new_layer_nodes[0].set_visible(False)
                self.new_layer_fill[0].set_visible(False)
                self.new_layer_nodes[0].remove()
                self.new_layer_fill[0].remove()
                self.new_layer_nodes = None
                self.new_layer_fill = None
                # RUN FINAL PART OF LAYER LOADING
                self.create_new_floating_layer()

        # RUN MODELLING ALGORITHMS
        if self.fault_picking_switch is False:
            self.run_algorithms()
        else:
            # UPDATE GMG GRAPHICS
            self.draw()

    def key_press(self, event):
        """DEFINE KEY PRESS LINKS"""

        if self.fault_picking_switch is True:
            # GMG IS IN FAULT MODE SO USE FAULT MODE KEY FUNCTIONS
            self.fault_mode_key_press(event)
            return

        # f = ACTIVATE FAULT PICKING MODE
        if event.key == 'f':
            # TURN ON/OFF FAULT PICKING MODE
            if self.fault_picking_switch is True:
                self.fault_picking_switch = False
            else:
                self.fault_picking_switch = True

        # i = INSERT NEW NODE AT MOUSE POSITION
        if event.key == 'i':
            if event.inaxes is None:
                return

            # GET CURRENT LAYER XY
            xt = np.array(self.current_x_nodes)
            yt = np.array(self.current_y_nodes)

            # INSERT NEW NODES INTO LAYER X AND Y LISTS
            self.current_x_nodes= np.insert(xt, [self.index_arg + 1], event.xdata)
            self.current_y_nodes = np.insert(yt, [self.index_arg + 1], event.ydata)

            self.currently_active_layer.set_data(self.current_x_nodes, self.current_y_nodes)

            # UPDATE LAYER DATA AND PLOT
            self.update_layer_data()
            self.run_algorithms()
            self.draw()

        # d = DELETE NODE AT MOUSE POSITION
        if event.key == 'd':
            xt = np.array(self.current_x_nodes)
            yt = np.array(self.current_y_nodes)

            # FIND NODE CLOSEST TO CURSOR LOCATION
            d = np.sqrt((xt - event.xdata) ** 2 + (yt - event.ydata) ** 2)
            self.index_arg = np.argmin(d)
            ind = d[self.index_arg]

            if xt[self.index_arg] == 0:  # PREVENT END NODES BEING DELETED
                return 0
            if ind >= self.node_click_limit:
                return 0
            else:
                # DELETE NODE BY RECREATING XY DATA WITHOUT CURRENT NODE
                self.current_x_nodes= [tup for i, tup in enumerate(self.plotx) if i != self.index_arg]  # DELETE X
                self.current_y_nodes = [tup for i, tup in enumerate(self.ploty) if i != self.index_arg]  # DELETE Y
                self.currently_active_layer.set_data(self.current_x_nodes, self.current_y_nodes)

            # NOW CHECK FOR PINCHED NODES
            index_arg2 = None
            self.pinch_switch = False
            self.index_arg2_list = [None] * (self.total_layer_count + 1)  # CREATE LIST OF NONES = LENGTH AS NUMB OF LAYERS
            for x in range(0, self.total_layer_count + 1):  # LOOP THROUGH ALL LAYERS TO CHECK FOR PINCHED NODES
                if x == self.currently_active_layer_id:
                    pass
                x_node_list = self.plotx_list[x]
                y_node_list = self.ploty_list[x]
                for i in range(0, len(x_node_list)):
                    # NOW CHECK X AND Y VALUES ARE EQUAL
                    if x_node_list[i] == xt[self.index_arg] and y_node_list[i] == yt[self.index_arg]:
                        # IF ONE OF THE NODES FORM LIST IS EQUAL TO A NODE FROM THE OTHER LAYER THEN RETURN THE INDEX
                        self.index_arg2_list[x] = i
                        self.pinch_switch = True

            # REMOVE PINCHED NODES
            if self.pinch_switch is True:
                for k in range(len(self.index_arg2_list)):
                    if self.index_arg2_list[k] is not None:
                        next_x_list, next_y_list = self.plotx_list[k], self.ploty_list[
                            k]  # GET THE NODE LIST OF THE NEXT LAYER
                        next_x_list = [tup for i, tup in enumerate(next_x_list) if
                                       i != self.index_arg2_list[k]]  # DELETE X
                        next_y_list = [tup for i, tup in enumerate(next_y_list) if
                                       i != self.index_arg2_list[k]]  # DELETE Y
                        self.plotx_list[k], self.ploty_list[
                            k] = next_x_list, next_y_list  # OVERWRITE THE NODE LIST WITH UPDATED LIST

            # SHIFT CURRENT NODE COLORING TO PREVIOUS NODE
            self.current_node.set_offsets([xt[self.index_arg - 1], yt[self.index_arg - 1]])

            # UPDATE GMG
            self.update_layer_data()
            self.run_algorithms()
            self.draw()

        # q = BEAT THE ZOOM BUG
        if event.key == 'q':
            self.nodes = True

        # n = CREATE NEW LAYER AT MOUSE POINT
        if event.key == 'n':
            self.new_layer(event)

        # b = LOCK OR UNLOCK LAYER BOUNDARY LOCKED MODE
        if event.key == 'b':
            if self.boundary_lock:
                'unlock layer'
                self.boundary_lock = False
                self.boundary_lock_list[self.currently_active_layer_id] = 1
            elif not self.boundary_lock:
                'lock layer'
                self.boundary_lock = True
                self.boundary_lock_list[self.currently_active_layer_id] = 0

        # l = LOCK OR UNLOCK LAYER/POLYGON MODE
        if event.key == 'l':
            if self.layer_lock:
                'unlock layer'
                self.layer_lock = False
                self.layer_lock_list[self.currently_active_layer_id] = 1
            elif not self.layer_lock:
                'lock layer'
                self.layer_lock = True
                self.layer_lock_list[self.currently_active_layer_id] = 0

        # < = MOVE TO NEXT LAYER
        if event.key == '.':
            if self.currently_active_layer_id == self.total_layer_count:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.currently_active_layer_id = 0
                self.density_input.SetValue(0.001 * self.densities[self.currently_active_layer_id])
                self.ref_density_input.SetValue(0.001 * self.reference_densities[self.currently_active_layer_id])
                self.susceptibility_input.SetValue(self.susceptibilities[self.currently_active_layer_id])
                self.angle_a_input.SetValue(self.angle_a[self.currently_active_layer_id])
                self.current_x_nodes= self.plotx_list[self.currently_active_layer_id]
                self.current_y_nodes = self.ploty_list[self.currently_active_layer_id]
                self.x_input.SetValue(self.current_x_nodes[0])
                self.y_input.SetValue(self.current_y_nodes[0])
                self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])
                self.update_layer_data()
                self.run_algorithms()
            else:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.currently_active_layer_id = self.currently_active_layer_id + 1
                self.nextdens = self.densities[self.currently_active_layer_id]
                self.density_input.SetValue(0.001 * self.densities[self.currently_active_layer_id])
                self.ref_density_input.SetValue(0.001 * self.reference_densities[self.currently_active_layer_id])
                self.susceptibility_input.SetValue(self.susceptibilities[self.currently_active_layer_id])
                self.angle_a_input.SetValue(self.angle_a[self.currently_active_layer_id])
                self.current_x_nodes= self.plotx_list[self.currently_active_layer_id]
                self.current_y_nodes = self.ploty_list[self.currently_active_layer_id]
                self.x_input.SetValue(self.current_x_nodes[0])
                self.y_input.SetValue(self.current_y_nodes[0])
                self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])
                self.update_layer_data()
                self.run_algorithms()

        # < = MOVE TO NEXT LAYER
        if event.key == ',':
            if self.currently_active_layer_id == 0:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.currently_active_layer_id = self.total_layer_count
                self.density_input.SetValue(0.001 * self.densities[self.currently_active_layer_id])
                self.ref_density_input.SetValue(0.001 * self.reference_densities[self.currently_active_layer_id])
                self.susceptibility_input.SetValue(self.susceptibilities[self.currently_active_layer_id])
                self.angle_a_input.SetValue(self.angle_a[self.currently_active_layer_id])
                self.current_x_nodes= self.plotx_list[self.currently_active_layer_id]
                self.current_y_nodes = self.ploty_list[self.currently_active_layer_id]
                self.x_input.SetValue(self.current_x_nodes[0])
                self.y_input.SetValue(self.current_y_nodes[0])
                self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])
                self.update_layer_data()
                self.run_algorithms()
            else:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.currently_active_layer_id = self.currently_active_layer_id - 1
                self.nextdens = self.densities[self.currently_active_layer_id]
                self.density_input.SetValue(0.001 * self.densities[self.currently_active_layer_id])
                self.ref_density_input.SetValue(0.001 * self.reference_densities[self.currently_active_layer_id])
                self.susceptibility_input.SetValue(self.susceptibilities[self.currently_active_layer_id])
                self.angle_a_input.SetValue(self.angle_a[self.currently_active_layer_id])
                self.current_x_nodes= self.plotx_list[self.currently_active_layer_id]
                self.current_y_nodes = self.ploty_list[self.currently_active_layer_id]
                self.x_input.SetValue(self.current_x_nodes[0])
                self.y_input.SetValue(self.current_y_nodes[0])
                self.current_node.set_offsets([self.current_y_nodes[0], self.current_y_nodes[0]])
                self.update_layer_data()
                self.run_algorithms()
            self.draw()

        # z = ZOOM IN MODE
        if event.key == 'z':
            self.zoom(event)

        # ctrl+z = ZOOM OUT
        if event.key == 'ctrl+z':
            self.zoom_out(event)

        # shift = PAN MODE
        if event.key == 'ctrl+p':
            self.pan(event)

        # a = FULL EXTENT VIEW'
        if event.key == 'a':
            self.full_extent(event)

        # p = TURN ON PINCH NODE MODE
        if event.key == 'p':
            if self.pinch_switch is False:
                self.pinch_switch = True
            else:
                self.pinch_switch = False
                self.pinch_node_list = [[], []]
                self.pinch_count = 0

        # ctrl+i = INCREASE LAYER TRANSPARENCY
        if event.key == 'ctrl+i':
            self.transparency_increase(event)

        # ctrl+d = INCREASE ASPECT TRANSPARENCY
        if event.key == 'ctrl+d':
            self.transparency_decrease(event)

        # up arrow = INCREASE ASPECT RATIO
        if event.key == 'up':
            self.aspect_increase(event)

        # ctrl+up = INCREASE ASPECT RATIO X2
        if event.key == 'ctrl+up':
            self.aspect_increase2(event)

        # down arrow = DECREASE ASPECT RATIO
        if event.key == 'down':
            self.aspect_decrease(event)

        # ctrl+down = DECREASE ASPECT RATIO
        if event.key == 'ctrl+down':
            self.aspect_decrease2(event)

    def toogle_fault_mode(self, event):
        """SWITCH FAULT PICKING MODE ON AND OFF"""
        if self.fault_picking_switch is True:
            self.fault_picking_switch = False

        elif self.fault_picking_switch is False:
            self.fault_picking_switch = True

    def fault_mode_key_press(self, event):
        """KEY PRESS CALLBACKS WHEN FAULT MODE IS ACTIVATED"""

        'i = INSERT NEW NODE AT MOUSE POSITION'
        if event.key == 'i':
            if event.inaxes is None:
                return

            # INSERT NEW NODE INTO XY LIST
            self.xt = np.insert(self.xt, [self.index_arg + 1], event.xdata)
            self.yt = np.insert(self.yt, [self.index_arg + 1], event.ydata)

            # UPDATE THE FAULT LIST RECORDS
            self.fault_x_coords_list[self.current_fault_index] = self.xt
            self.fault_y_coords_list[self.current_fault_index] = self.yt

            # UPDATE FAULT GRAPHICS
            self.faults[self.current_fault_index][0].set_xdata(self.xt)
            self.faults[self.current_fault_index][0].set_ydata(self.yt)

            # UPDATE CURRENT FAULT OVERLAY GRAPHIC
            self.faultline.set_data(self.xt, self.yt)

        'd = DELETE NODE AT MOUSE POSITION'
        if event.key == 'd':
            if event.inaxes is None:
                return
            # FIND NODE CLOSEST TO CURSOR LOCATION
            d = np.sqrt((self.xt - event.xdata) ** 2 + (self.yt - event.ydata) ** 2)
            self.index_arg = np.argmin(d)
            self.distance = d[self.index_arg]

            if self.index_arg == 0 or self.index_arg == (len(self.fault_x_coords_list[self.current_fault_index]) - 1):
                # PREVENT END NODES BEING DELETED
                return 0
            if self.distance >= self.node_click_limit:
                # CLICK WAS TO FAR AWAY FROM A NODE TO DELETE IT
                return 0
            else:
                # DELETE NODE BY RECREATING XY DATA WITHOUT CURRENT NODE
                self.fault_x_coords_list[self.current_fault_index] = \
                    [tup for i, tup in enumerate(self.xt) if i != self.index_arg]  # DELETE X

                self.fault_y_coords_list[self.current_fault_index] = \
                    [tup for i, tup in enumerate(self.yt) if i != self.index_arg]  # DELETE Y

                # RESET CURRENT NOT POSITION TO FIRST NODE
                self.current_node.set_offsets([self.xt[0], self.yt[0]])

                # UPDATE CURRENT FAULT OVERLAY GRAPHIC
                self.faultline.set_data(self.xt, self.yt)

            # UPDATE GMG
            self.update_layer_data()

        '< = INCREMENT WHICH FAULT IS BEING EDITED'
        if event.key == ',':
            if self.current_fault_index <= self.fault_counter - 1 and self.current_fault_index > 0:
                # INCREMENT TO NEXT FAULT
                self.current_fault_index -= 1
            else:
                # GO TO NEWEST FAULT
                self.current_fault_index = self.fault_counter - 1

            # UPDATE CURRENT PLOT GRAPHICS
            self.faultline.set_data(self.fault_x_coords_list[self.current_fault_index],
                                    self.fault_y_coords_list[self.current_fault_index])

            self.xt = self.fault_x_coords_list[self.current_fault_index]
            self.yt = self.fault_y_coords_list[self.current_fault_index]

        '> = INCREMENT WHICH FAULT IS BEING EDITED'
        if event.key == '.':
            if self.current_fault_index < self.fault_counter - 1:
                # INCREMENT TO NEXT FAULT
                self.current_fault_index += 1
            elif self.current_fault_index == self.fault_counter - 1:
                # GO BACK TO FIRST FAULT
                self.current_fault_index = 0

            # UPDATE CURRENT PLOT GRAPHICS
            self.faultline.set_data(self.fault_x_coords_list[self.current_fault_index],
                                    self.fault_y_coords_list[self.current_fault_index])

            self.xt = self.fault_x_coords_list[self.current_fault_index]
            self.yt = self.fault_y_coords_list[self.current_fault_index]

        # UPDATE GMG
        self.update_layer_data()

    def pick_new_fault(self, event):
        """FAULT PICKING/LINE DRAWING MODE"""

        # CHECK IF FAULT PICKING MODE IS ON
        if self.fault_picking_switch is False:
            MessageDialog(self, -1, "Faulting picking mode is not activated.\nTurn on fault picking mode first.",
                          "Fault picker")
        else:
            # PROMPT NEW FAULT DIALOG BOX
            new_fault_dialogbox = NewFaultDialog(self, -1, 'Create New Fault')
            answer = new_fault_dialogbox.ShowModal()

            # GET NEW FAULT VALUES
            self.new_x1, self.new_y1 = new_fault_dialogbox.x1, new_fault_dialogbox.y1
            self.new_x2, self.new_y2 = new_fault_dialogbox.x2, new_fault_dialogbox.y2

            if self.fault_counter == 0:
                # CREATE NEW CURRENT FAULT GRAPHIC
                self.faultline, = self.model_frame.plot([-50000, 0], [-49999, 0], marker='s', color='m', linewidth=0.75,
                                                    alpha=1.0, zorder=2, picker=True)
            # INCREMENT THE TOTAL FAULT COUNT
            self.fault_counter += 1

            # SET THE CURRENT FAULT INDEX AS THE NEW FAULT
            self.current_fault_index = self.fault_counter - 1

            # APPEND BLANKS TO THE OBJECTS (USED FOR THE NEXT NEW FAULT)
            self.faults.append([])
            self.fault_names_list.append([])
            self.fault_x_coords_list.append([])
            self.fault_y_coords_list.append([])
            self.fault_colors.append('k')

            # LIST OF FAULT NAMES
            self.fault_tree_items.append('fault %s' % (int(self.current_fault_index)))

            # APPEND THE NEW FAULT TO THE FAULT TREE SIDE PANEL USING add_new_tree_nodes FUNC
            self.fault_item = 'fault %s' % (int(self.current_fault_index))
            self.add_new_tree_nodes(self.fault_tree_root, self.fault_item, self.current_fault_index)

            self.fold_panel_three.Collapse()
            self.fold_panel_three.Expand()

            # ADD NAME FOR NEW FAULT
            self.fault_names_list[self.current_fault_index] = 'Fault'

            # ADD NEW FAULT XY TO XY LISTS
            self.fault_x_coords_list[self.current_fault_index] = [self.new_x1, self.new_x2]
            self.fault_y_coords_list[self.current_fault_index] = [self.new_y1, self.new_y2]

            # SET CURRENT FAULT X AND Y ARRAYS
            self.xt = [self.new_x1, self.new_x2]
            self.yt = [self.new_y1, self.new_y2]

            # ADD FAULT LINE TO CANVAS
            self.faults[self.current_fault_index] = self.model_frame.plot(
                self.fault_x_coords_list[self.current_fault_index],
                self.fault_y_coords_list[self.current_fault_index],
                color='k', linewidth=0.5, zorder=1, alpha=1.0)

            # UPDATE CURRENT PLOT GRAPHICS
            self.faultline.set_data(self.fault_x_coords_list[self.current_fault_index],
                                    self.fault_y_coords_list[self.current_fault_index])

            # UPDATE CURRENT NODE RED DOT GRAPHIC
            self.current_node.set_offsets([self.new_x1, self.new_y1])

            # UPDATE GMG
            self.draw()

    def fault_activated(self, event):
        """RESPONSE WHEN A FAULT NAME IS SELECTED"""

        # GET THE SELECTED FAULT INDEX NUMBER
        self.current_fault_index = self.fault_tree.GetPyData(event.GetItem())

        if self.fault_picking_switch is False:
            self.fault_picking_switch = True

        # SET CHECKBOX AS CHECKED
        self.fault_tree.GetSelection().Check(checked=True)

        # UPDATE CURRENT PLOT GRAPHICS
        self.faultline.set_data(self.fault_x_coords_list[self.current_fault_index],
                                self.fault_y_coords_list[self.current_fault_index])

        # UPDATE GRAPHICS WITH CURRENT FAULT SELECTED
        self.update_layer_data()

    def fault_checked(self, event):
        """TOGGLE WHETHER OR NOT A FAULT WILL BE PLOTTED IN THE MODEL FIGURE"""
        i = self.fault_tree.GetPyData(event.GetItem())

        if self.faults[i][0].get_visible() == True:
            # HIDE FAULT
            self.faults[i][0].set_visible(False)
            self.faultline.set_visible(False)
        else:
            # SHOW FAULT
            self.faults[i][0].set_visible(True)
            self.faultline.set_visible(True)

        # UPDATE FIGURE
        self.draw()

    def on_begin_edit_fault_label(self, event):
        self.current_fault_index = self.fault_tree.GetPyData(event.GetItem())

    def on_end_edit_fault_label(self, event):
        new_label = self.fault_tree.GetItemText(event.GetItem())
        self.fault_tree_items[self.current_fault_index] = str(new_label)

    def write_layers_xy(self, event):
        """OUTPUT LAYER DATA TO FILE"""

        # CREATE OUTPUT FILE
        save_file_dialog = wx.FileDialog(self, "Save XY data", "", "", "xy files (*.xy)|*.xy",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THERE MIND

        # OUTPUT FILE
        all_layers_output_file = save_file_dialog.GetPath()
        # THE OUTPUT DIRECTORY
        output_dir = os.path.dirname(all_layers_output_file)

        # NOW WRITE OUT THE DATA
        with open(all_layers_output_file, 'wb') as f:
            try:
                for i in range(1, self.total_layer_count + 1):
                    out = csv.writer(f, delimiter=' ')
                    f.write('>\n')
                    data = [self.plotx_list[i], self.ploty_list[i]]
                    out.writerows(zip(*data))
                    layer_write = zip(self.plotx_list[i], self.ploty_list[i])
                    # WRITE INDIVIDUAL LAYER
                    np.savetxt(output_dir + '/' + self.loaded_tree_items[i] + '.xy', layer_write, delimiter=' ',
                               fmt='%f %f')
                    f.close()
            except IndexError:
                f.close()
                pass

    def write_c_xy(self, event):
        # CREATE OUTPUT FILE
        save_file_dialog = wx.FileDialog(self, "Save XY data", "", "", "xy files (*.xy)|*.xy",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THEIR MIND

        # NOW WRITE OUT THE DATA
        output_stream = save_file_dialog.GetPath()
        with open(output_stream, 'wb') as f:
            # LAYER NODES
            for i in range(0, self.total_layer_count + 1):
                f.write('B  {0}\n'.format(i + 1))
                data = zip(self.plotx_list[i], self.ploty_list[i], np.ones(len(self.ploty_list[i])))
                # print data
                np.savetxt(f, data, delimiter=' ', fmt='%6.02f %3.02f %1d')

            # VELOCITY NODES
            for i in range(0, self.total_layer_count):
                density = (self.background_density + self.densities[i])

                # CONVERT DENSITY TO VELOCITY USING GARNERS RULE
                velocity = round((m.pow((density / 1670.), (1. / .25))), 2)

                # CONVERT DENSITY TO VELOCITY USING NAFE-DRAKE EQUATION
                # velocity = (1.6612*density) - (0.4721*density)**2 + (0.0671*density)**3 -
                # (0.0043*density)**4 + (0.000106*density)**5

                # FORMAT c.in FILE
                f.write('B  {0}\n'.format(i))
                data = zip(self.plotx_list[i], np.linspace(velocity, velocity, len(self.ploty_list[i])),
                           np.ones(len(self.ploty_list[i])), np.linspace(velocity, velocity, len(self.ploty_list[i])),
                           np.ones(len(self.ploty_list[i])))

                # OUTPUT FILE
                np.savetxt(f, data, delimiter=' ', fmt='%6.02f %3.02f %1d %3.02f %1d')

    def capture_coordinates(self, event):
        if self.capture is False:
            self.capture = True

            # CREATE INSTANCE OF CAPTURE COORDINATES
            self.capture_window = CaptureCoordinates(self, -1, 'Capture Coordinates')
            self.capture_window.Show(True)

    def pinch_out_layer(self, event):
        pinch_box = PinchDialog(self, -1, 'Pinch Out Layer:', self.plotx_list, self.ploty_list, self.currently_active_layer_id)
        answer = pinch_box.ShowModal()
        self.current_x_nodes= pinch_box.pinched_x
        self.current_y_nodes = pinch_box.pinched_y
        self.update_layer_data()
        self.draw()

    def depinch_layer(self, event):
        depinch_box = DepinchDialog(self, -1, 'Depinch layer', self.plotx_list, self.ploty_list, self.currently_active_layer_id,
                                    self.total_layer_count)
        answer = depinch_box.ShowModal()
        self.current_x_nodes= depinch_box.depinched_x
        self.current_y_nodes = depinch_box.depinched_y
        self.update_layer_data()
        self.draw()

    def bulk_shift(self, event):
        bulk_shift_box = BulkShiftDialog(self, -1, 'Layer bulk shift', self.plotx_list, self.ploty_list, self.currently_active_layer_id)
        answer = bulk_shift_box.ShowModal()
        self.current_x_nodes= bulk_shift_box.new_x
        self.current_y_nodes = bulk_shift_box.new_y
        self.update_layer_data()
        self.draw()

    def filter_observed_topography(self, event):
        """FILTER OBSERVED TOPOGRAPHY USING MEDIAN FILTER - CALLS class MedianFilterDialog"""

        # RUN FILTER
        median_filter_box = MedianFilterDialog(self, -1, 'median filter', self.observed_topography_list)
        answer = median_filter_box.ShowModal()

        # CREATE NEW OBSERVED GRAVITY OBJECT
        observed = ObservedData()

        # SET ATTRIBUTES
        observed.id = int(self.observed_topography_counter)
        observed.type = str('filtered')
        observed.name = median_filter_box.output_name
        observed.color = median_filter_box.output_color
        observed.data = median_filter_box.filtered_output
        observed.mpl_actor = self.topo_frame.scatter(observed.data[:, 0], observed.data[:, 1], marker='o',
                                                     color=observed.color, s=5, gid=observed.id)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_topography_list.append(observed)

        # APPEND NEW DATA MENU TO 'GRAV data MENU'
        self.topo_submenu = wx.Menu()
        self.m_topo_submenu.Append(10000+observed.id, observed.name, self.topo_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.topo_submenu.Append(10000+observed.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_observed_topography, id=10000+observed.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_topography_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def filter_observed_gravity(self, event):
        """FILTER OBSERVED ANOMALY USING MEDIAN FILTER - CALLS class MedianFilterDialog"""

        # RUN FILTER
        median_filter_box = MedianFilterDialog(self, -1, 'median filter', self.observed_gravity_list)
        answer = median_filter_box.ShowModal()

        # CREATE NEW OBSERVED GRAVITY OBJECT
        observed = ObservedData()

        # SET ATTRIBUTES
        observed.id = int(self.observed_gravity_counter)
        observed.type = str('filtered')
        observed.name = median_filter_box.output_name
        observed.color = median_filter_box.output_color
        observed.data = median_filter_box.filtered_output
        observed.mpl_actor = self.gravity_frame.scatter(observed.data[:, 0],
                                                                observed.data[:, 1], marker='o',
                                                                color=observed.color, s=5,
                                                                gid=observed.id)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_gravity_list.append(observed)

        # TURN ON OBSERVED GRAVITY SWITCH
        self.observed_gravity_switch = True

        # APPEND NEW DATA MENU TO 'GRAV data MENU'
        self.grav_submenu = wx.Menu()
        self.m_obs_g_submenu.Append(11000+observed.id, observed.name, self.grav_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.grav_submenu.Append(11000+observed.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000+observed.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_gravity_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def filter_observed_magnetic(self, event):
        """FILTER OBSERVED MAGNETIC USING MEDIAN FILTER - CALLS class MedianFilterDialog"""

        # RUN FILTER
        median_filter_box = MedianFilterDialog(self, -1, 'median filter', self.observed_magnetic_list)
        answer = median_filter_box.ShowModal()

        # CREATE NEW OBSERVED GRAVITY OBJECT
        observed = ObservedData()

        # SET ATTRIBUTES
        observed.id = int(self.observed_magnetic_counter)
        observed.type = str('filtered')
        observed.name = median_filter_box.output_name
        observed.color = median_filter_box.output_color
        observed.data = median_filter_box.filtered_output
        observed.mpl_actor = self.magnetic_frame.scatter(observed.data[:, 0], observed.data[:, 1], marker='o',
                                                        color=observed.color, s=5, gid=observed.id)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_magnetic_list.append(observed)

        # APPEND NEW DATA MENU TO 'GRAV data MENU'
        self.mag_submenu = wx.Menu()
        self.m_obs_mag_submenu.Append(12000+observed.id, observed.name, self.mag_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.mag_submenu.Append(12000+observed.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12000+observed.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_magnetic_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def take_topography_horizontal_derivative(self, event):
        """
        TAKE HORIZONTAL DERIVATIVE OF OBSERVED DATA.
        CALLS class HorizontalDerivative
        """

        # OPEN THE HORIZONTAL DERIVATIVE INPUT WINDOW
        horizontal_derivative_box = HorizontalDerivative(self, -1, 'Horizontal derivative',
                                                         self.observed_topography_list)
        answer = horizontal_derivative_box.ShowModal()

        # CREATE NEW DATA OBJECT AND PARSE OUTPUT TO THE OBJECT
        new_derivative = ObservedData()
        new_derivative.data = horizontal_derivative_box.deriv
        new_derivative.name = horizontal_derivative_box.output_name
        new_derivative.color = horizontal_derivative_box.output_color
        new_derivative.id = 10000 + self.observed_topography_counter
        new_derivative.type = str('derivative')
        new_derivative.mpl_actor = self.topo_d_frame.scatter(new_derivative.data[:, 0], new_derivative.data[:, 1],
                                                         marker='o', color=new_derivative.color, s=5,
                                                         gid=10000 + self.observed_topography_counter)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_topography_list.append(new_derivative)

        #  APPEND NEW MENUBAR TO THE GRAVITY MENUBAR
        self.topo_submenu = wx.Menu()
        self.m_topo_submenu.Append(10000 + self.observed_topography_counter, new_derivative.name,
                                    self.topo_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.topo_submenu.Append(10000 + self.observed_topography_counter, 'delete observed data')

        # BIND TO DEL FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_topo, id=10000 + self.observed_topography_counter)

        # INCREMENT GRAV DERIV COUNTER
        self.observed_topography_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.draw()

    def take_gravity_horizontal_derivative(self, event):
        """
        TAKE HORIZONTAL DERIVATIVE OF OBSERVED DATA.
        CALLS class HorizontalDerivative
        """

        # OPEN THE HORIZONTAL DERIVATIVE INPUT WINDOW
        horizontal_derivative_box = HorizontalDerivative(self, -1, 'Horizontal derivative', self.observed_gravity_list)
        answer = horizontal_derivative_box.ShowModal()

        # CREATE NEW DATA OBJECT AND PARSE OUTPUT TO THE OBJECT
        new_derivative = ObservedData()
        new_derivative.data = horizontal_derivative_box.deriv
        new_derivative.name = horizontal_derivative_box.output_name
        new_derivative.color = horizontal_derivative_box.output_color
        new_derivative.id = 11000 + self.observed_gravity_counter
        new_derivative.type = str('derivative')
        new_derivative.mpl_actor = self.gravity_d_frame.scatter(new_derivative.data[:, 0], new_derivative.data[:, 1],
                                                         marker='o', color=new_derivative.color, s=5,
                                                         gid=11000 + self.observed_gravity_counter)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_gravity_list.append(new_derivative)

        #  APPEND NEW MENUBAR TO THE GRAVITY MENUBAR
        self.grav_submenu = wx.Menu()
        self.m_obs_g_submenu.Append(11000 + self.observed_gravity_counter, new_derivative.name, self.grav_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.grav_submenu.Append(11000 + self.observed_gravity_counter, 'delete observed data')

        # BIND TO DEL FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000 + self.observed_gravity_counter)

        # INCREMENT GRAV DERIV COUNTER
        self.observed_gravity_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.draw()

    def take_magnetic_horizontal_derivative(self, event):
        """
        TAKE HORIZONTAL DERIVATIVE OF OBSERVED DATA.
        CALLS class HorizontalDerivative
        """

        # OPEN THE HORIZONTAL DERIVATIVE INPUT WINDOW
        horizontal_derivative_box = HorizontalDerivative(self, -1, 'Horizontal derivative', self.observed_magnetic_list)
        answer = horizontal_derivative_box.ShowModal()

        # CREATE NEW DATA OBJECT AND PARSE OUTPUT TO THE OBJECT
        new_derivative = ObservedData()
        new_derivative.data = horizontal_derivative_box.deriv
        new_derivative.name = horizontal_derivative_box.output_name
        new_derivative.color = horizontal_derivative_box.output_color
        new_derivative.id = 12000 + self.observed_magnetic_counter
        new_derivative.type = str('derivative')
        new_derivative.mpl_actor = self.magnetic_d_frame.scatter(new_derivative.data[:, 0], new_derivative.data[:, 1],
                                                             marker='o', color=new_derivative.color, s=5,
                                                             gid=12000 + self.observed_magnetic_counter)

        # APPEND NEW DATA TO THE OBSERVED GRAVITY GMG LIST
        self.observed_magnetic_list.append(new_derivative)

        #  APPEND NEW MENUBAR TO THE GRAVITY MENUBAR
        self.mag_submenu = wx.Menu()
        self.m_obs_mag_submenu.Append(12500 + self.observed_magnetic_counter, new_derivative.name,
                                      self.mag_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.mag_submenu.Append(12500 + self.observed_magnetic_counter, 'delete observed data')

        # BIND TO DEL FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12500 + self.observed_magnetic_counter)

        # INCREMENT GRAV DERIV COUNTER
        self.observed_magnetic_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.draw()

    def new_layer(self, event):
        new_layer_dialogbox = NewLayerDialog(self, -1, 'Create New Layer')
        answer = new_layer_dialogbox.ShowModal()

        if new_layer_dialogbox.fixed:
            # CREATING A NEW FIXED LAYER

            # INCREMENT THE CURRENT LAYER INDEX VALUE (self.currently_active_layer_id)
            self.total_layer_count += 1

            # SET THE ACTIVE LAYER AS THE NEWLY CREATED LAYER
            self.currently_active_layer_id = self.total_layer_count
            print("self.currently_active_layer_id = %s") % self.currently_active_layer_id
            print("self.total_layer_count = %s") % self.total_layer_count
            # CREATE A NEW LAYER OBJECT
            new_layer = Layer()

            # SET SOME OF THE NEW LAYERS ATTRIBUTES
            new_layer.id = self.currently_active_layer_id
            new_layer.name = str('layer %s') % self.currently_active_layer_id
            new_layer.type = str('fixed')

            # ADD NEW LAYER TO THE LAYER TREE DISPLAY
            self.tree_items.append('layer %s' % (int(self.currently_active_layer_id)))
            self.item = 'layer %s' % (int(self.currently_active_layer_id))
            self.add_new_tree_nodes(self.root, self.item, self.currently_active_layer_id)

            # DETERMINE WHICH LAYER IS THE LAST PREVIOUS "FIXED LAYER"; SET THIS LAYER AS "previous_fixed_layer"
            if self.total_layer_count > 0:
                for i in range(0, self.total_layer_count):
                    if self.layer_list[i].type == 'fixed':
                        previous_fixed_layer = i
                    else:
                        continue
            else:
                previous_fixed_layer = 0

            # SET NEW LAYER NODES
            layer_above_x = np.array(self.layer_list[previous_fixed_layer].x_nodes)
            layer_above_y = np.array(self.layer_list[previous_fixed_layer].y_nodes)
            new_layer_thickness = new_layer_dialogbox.new_thickness
            new_layer.x_nodes = layer_above_x
            new_layer.y_nodes = layer_above_y + new_layer_thickness

            # CREATE LAYER LINE
            new_layer.node_mpl_actor = self.model_frame.plot(new_layer.x_nodes, new_layer.y_nodes, color='blue',
                                                                   linewidth=1.0, alpha=1.0)
            # CREATE LAYER POLYGON FILL
            new_layer.polygon_mpl_actor = self.model_frame.fill(new_layer.x_nodes, new_layer.y_nodes, color='blue',
                                                                         alpha=self.layer_transparency, closed=True,
                                                                         linewidth=None, ec=None)

            # SET CURRENTLY ACTIVE LAYER AS THE NEW LAYER
            # self.currently_active_layer.set_data(new_layer.x_nodes, new_layer.y_nodes)
            self.currently_active_layer.set_xdata(new_layer.x_nodes)
            self.currently_active_layer.set_ydata(new_layer.y_nodes)
            self.currently_active_layer.set_color(new_layer.color)

            # SET CURRENTLY ACTIVE LAYER NODE OBJECTS
            self.current_x_nodes = new_layer.x_nodes
            self.current_y_nodes = new_layer.y_nodes

            # SET THE CURRENT NODE
            self.current_node.set_offsets([new_layer.x_nodes[0], new_layer.y_nodes[0]])

            self.nextdens = new_layer.density

            # SET CURRENT ATTRIBUTE INPUTS IN LEFT PANEL
            self.density_input.SetValue(new_layer.density)
            self.ref_density_input.SetValue(new_layer.reference_density)
            self.susceptibility_input.SetValue(new_layer.susceptibility)
            self.angle_a_input.SetValue(new_layer.angle_a)
            self.angle_b_input.SetValue(new_layer.angle_b)
            
            # APPEND NEW LAYER TO THE LAYER LIST
            self.layer_list.append(new_layer)
            
            # UPDATE GMG FRAME
            print("y nodes =" )
            print self.layer_list[1].y_nodes
            self.update_layer_data()
            self.draw()

        elif not new_layer_dialogbox.fixed:
            # CREATEING A NEW FLOATING LAYER
            self.new_plotx = []
            self.new_ploty = []
            self.click_count = 0
            self.new_layer_nodes = None
            self.new_layer_fill = None
            # SWITCH ON MOUSE CLICK CAPTURE MODE (SEE button_release func for continuation of code)
            self.select_new_layer_nodes = True

        else:
            # USER CHANGED THEIR MIND - NO NEW LAYER ADDED'
            pass

    def create_new_floating_layer(self):
            """CREATE A NEW FLOATING LAYER USING FOUR USER INPUT MOUSE CLICKS"""
            # SOURCE NEW NODES FROM USER CLICKS
            self.current_x_nodes= self.new_plotx
            self.current_y_nodes = self.new_ploty

            # INCREMENT THE LAYER COUNT'
            self.currently_active_layer_id = self.total_layer_count
            self.currently_active_layer_id = self.currently_active_layer_id + 1

            # INCREMENT THE TOTAL LAYER COUNT
            self.total_layer_count += 1

            # ADD A NEW BLANK LAYER TO THE PLOT LISTS
            self.polygon_fills.append([])
            self.layer_lines.append([])
            self.plotx_list.append([])
            self.ploty_list.append([])
            self.poly_fills.append([])
            self.layer_colors.append('black')
            self.densities.append(0.)
            self.reference_densities.append(0.)
            self.susceptibilities.append(0.)
            self.angle_a.append(0.)
            self.angle_b.append(0.)
            self.layer_lock_list.append(1)
            self.boundary_lock_list.append(1)
            self.layer_lock_status.append('unlocked')
            self.boundary_lock_status.append('unlocked')
            self.tree_items.append('layer %s' % (int(self.currently_active_layer_id)))
            self.item = 'layer %s' % (int(self.currently_active_layer_id))
            self.add_new_tree_nodes(self.root, self.item, self.currently_active_layer_id)
            self.layers_calculation_switch.append(1)

            # CREATE LAYER LINE
            self.layer_lines[self.currently_active_layer_id] = self.model_frame.plot(self.plotx_list[self.currently_active_layer_id],
                                                                     self.ploty_list[self.currently_active_layer_id],
                                                         color='blue', linewidth=1.0, alpha=1.0)
            # CREATE LAYER POLYGON FILL
            self.plotx_polygon = np.array(self.current_x_nodes)
            self.ploty_polygon = np.array(self.current_y_nodes)
            self.polygon_fills[self.currently_active_layer_id] = self.model_frame.fill(self.plotx_polygon, self.ploty_polygon,
                                                                       color='blue', alpha=self.layer_transparency,
                                                                       closed=True, linewidth=None, ec=None)

            # UPDATE LAYER DATA
            self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])
            self.nextdens = self.densities[self.currently_active_layer_id]
            self.density_input.SetValue(0.001 * self.densities[self.currently_active_layer_id])
            self.ref_density_input.SetValue(0.001 * self.reference_densities[self.currently_active_layer_id])
            self.susceptibility_input.SetValue(self.susceptibilities[self.currently_active_layer_id])
            self.angle_a_input.SetValue(self.angle_a[self.currently_active_layer_id])
            self.update_layer_data()
            self.run_algorithms()
            self.draw()

    def load_layer(self, event):
        open_file_dialog = wx.FileDialog(self, "Open Layer", "", "", "Layer XY files (*.txt)|*.txt",
                                         wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if open_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THERE MIND

        file_in = open_file_dialog.GetPath()
        new_layer = np.genfromtxt(file_in, autostrip=True, delimiter=' ', dtype=float)

        # INCREMENT THE LAYER COUNT
        self.currently_active_layer_id = self.total_layer_count
        self.currently_active_layer_id += 1

        # INCREMENT THE TOTAL LAYER COUNT
        self.total_layer_count += 1

        # ADD A NEW BLANK LAYER TO THE PLOT LISTS

        self.polygon_fills.append([])
        self.layer_lines.append([])
        self.plotx_list.append([])
        self.ploty_list.append([])
        self.poly_fills.append([])
        self.layer_colors.append('black')
        self.densities.append(0.)
        self.reference_densities.append(0.)
        self.susceptibilities.append(0.)
        self.angle_a.append(0.)
        self.angle_b.append(0.)
        self.layer_lock_list.append(1)
        self.boundary_lock_list.append(1)
        self.layer_lock_status.append('unlocked')
        self.boundary_lock_status.append('unlocked')
        self.tree_items.append('layer %s' % (int(self.currently_active_layer_id + 1)))
        self.item = 'layer %s' % (int(self.currently_active_layer_id + 1))
        self.add_new_tree_nodes(self.root, self.item, self.currently_active_layer_id)
        self.layers_calculation_switch.append(0)
        self.current_x_nodes= new_layer[:, 0]
        self.current_y_nodes = new_layer[:, 1]

        # CREATE LAYER LINE
        self.layer_lines[self.currently_active_layer_id] = self.model_frame.plot(self.plotx_list[self.currently_active_layer_id],
                                                                     self.ploty_list[self.currently_active_layer_id],
                                                                     color='blue', linewidth=1.0, alpha=1.0)
        # CREATE LAYER POLYGON FILL
        self.plotx_polygon = np.array(self.current_x_nodes)
        self.ploty_polygon = np.array(self.current_y_nodes)
        self.polygon_fills[self.currently_active_layer_id] = self.model_frame.fill(self.plotx_polygon, self.ploty_polygon,
                                                                       color='blue', alpha=self.layer_transparency,
                                                                       closed=True, linewidth=None, ec=None)
        # UPDATE LAYER DATA
        self.nextdens = self.densities[self.currently_active_layer_id]
        self.density_input.SetValue(0.001 * self.densities[self.currently_active_layer_id])
        self.ref_density_input.SetValue(0.001 * self.reference_densities[self.currently_active_layer_id])
        self.susceptibility_input.SetValue(self.susceptibilities[self.currently_active_layer_id])
        self.angle_a_input.SetValue(self.angle_a[self.currently_active_layer_id])
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def delete_layer(self, event):
        """Delete LAYER DATA"""

        # SET CURRENT NODE AS A OFF STAGE (PLACE HOLDER)
        self.current_node = self.model_frame.scatter(-40000., 0., marker='o', color='r', zorder=10)

        # REMOVE POLYGON
        self.polygon_fills[self.currently_active_layer_id][0].remove()
        del self.polygon_fills[self.currently_active_layer_id]
        # REMOVE POLYGON
        self.layer_lines[self.currently_active_layer_id][0].remove()
        del self.layer_lines[self.currently_active_layer_id]
        # SET CURRENT AS THE NEXT LAYER'
        if self.currently_active_layer_id == self.total_layer_count:
            self.currently_active_layer.set_xdata(self.plotx_list[0])
            self.currently_active_layer.set_ydata(self.ploty_list[0])
            self.currently_active_layer.set_color(self.layer_colors[0])
            self.current_x_nodes= self.plotx_list[0]
            self.current_y_nodes = self.ploty_list[0]
        else:
            self.currently_active_layer.set_xdata(self.plotx_list[self.currently_active_layer_id + 1])
            self.currently_active_layer.set_ydata(self.ploty_list[self.currently_active_layer_id + 1])
            self.currently_active_layer.set_color(self.layer_colors[self.currently_active_layer_id + 1])
            self.current_x_nodes= self.plotx_list[self.currently_active_layer_id + 1]
            self.current_y_nodes = self.ploty_list[self.currently_active_layer_id + 1]
        self.draw()

        # REMOVE META DATA
        del self.plotx_list[self.currently_active_layer_id]
        del self.ploty_list[self.currently_active_layer_id]
        del self.densities[self.currently_active_layer_id]
        del self.reference_densities[self.currently_active_layer_id]
        del self.susceptibilities[self.currently_active_layer_id]
        del self.angle_a[self.currently_active_layer_id]
        del self.angle_b[self.currently_active_layer_id]
        del self.layer_lock_list[self.currently_active_layer_id]
        del self.boundary_lock_list[self.currently_active_layer_id]
        del self.layer_lock_status[self.currently_active_layer_id]
        del self.boundary_lock_status[self.currently_active_layer_id]
        del self.layer_colors[self.currently_active_layer_id]
        del self.layers_calculation_switch[self.currently_active_layer_id]

        #  REMOVE TREE ITEMS
        del self.tree_items[self.currently_active_layer_id-1]
        layers = self.tree.GetRootItem().GetChildren()
        self.tree.Delete(layers[self.currently_active_layer_id-1])
        # RESET TREE ITEM ID'S
        layers = self.tree.GetRootItem().GetChildren()
        for i in range(len(layers)):
            self.tree.SetPyData(layers[i], i)

        # UPDATE COUNTERS
        self.total_layer_count -= 1

        # SET CURRENT LAYER AS PREVIOUS LAYER TO THAT JUST DELETED
        self.currently_active_layer_id -= 1

        # UPDATE
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    # LAYER AND MODEL ATTRIBUTE CONTROLS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_background_density(self, event):
        grav_box = SetBackgroundDensityDialog(self, -1, 'Set background density')
        answer = grav_box.ShowModal()
        self.background_density_upper = float(grav_box.background_density_upper)

        for i in range(0, self.total_layer_count):
            self.reference_densities[i] = float(self.background_density_upper) * 1000.
        # self.background_density_upper = float((grav_box.background_density_lower))
        # self.background_density_upper = float((grav_box.background_density_lid))
        self.absolute_densities = True
        self.draw()

    def set_density(self, value):
        if self.density_input.GetValue() == 0:
            self.densities[self.currently_active_layer_id] = 0
        else:
            self.densities[self.currently_active_layer_id] = float(self.density_input.GetValue() * 1000.)

    def set_reference_density(self, value):
        if self.ref_density_input.GetValue() == 0:
            self.reference_densities[self.currently_active_layer_id] = 0
        else:
            self.reference_densities[self.currently_active_layer_id] = float(self.ref_density_input.GetValue() * 1000.)

    def set_susceptibility(self, value):
        self.susceptibilities[self.currently_active_layer_id] = float(self.susceptibility_input.GetValue())

    def set_angle_a(self, value):
        self.angle_a[self.currently_active_layer_id] = float(self.angle_a_input.GetValue())

    def set_angle_b(self, value):
        self.angle_b[self.currently_active_layer_id] = float(self.angle_b_input.GetValue())

    def set_text_size(self, value):
        """GET NEW TEXT SIZE"""

        self.textsize = float(self.text_size_input.GetValue())

        # WELL DATA
        # LOOP THROUGH ALL WELL NAMES
        for i in range(len(self.well_data_list)):
            self.well_data_list[i].text_size = self.textsize
            self.well_data_list[i].mpl_actor_name.set_size(self.textsize)

            # LOOP THROUGH ALL WELL HORIZON LABELS
            for l in range(len(self.well_data_list[i].labels_list)):
                if self.well_data_list[i].labels_list[l] is not None:
                    self.well_data_list[i].labels_list[l].set_size(self.textsize)

        # # LOOP THROUGH OUTCROP DATA LABELS
        if self.outcrop_data_count > 0:
            for i in range(self.outcrop_data_count):
                if self.outcrop_data_list[i] is not None:
                    for t in range(len(self.outcrop_data_list[i].labels)):
                        self.outcrop_data_list[i].labels[t].set_fontsize(self.textsize)

        # REDRAW ANNOTATIONS WITH NEW TEXT SIZE
        self.draw()

    def set_obs_grav_rms(self, value):
        """SET THE DATA TO BE USED FOR CALCULATING THE RMS MISTFIT"""
        selection = SetObsRmsDialog(self, -1, 'Set RMS Input', self.observed_gravity_list)
        answer = selection.ShowModal()
        for i in range(0, len(self.observed_gravity_list)):
                if self.observed_gravity_list[i].name == selection.obs_name:
                    self.obs_gravity_data_for_rms = self.observed_gravity_list[i].data

    def set_obs_mag_rms(self, value):
        """SET THE DATA TO BE USED FOR CALCULATING THE RMS MISTFIT"""
        selection = SetObsRmsDialog(self, -1, 'Set RMS Input', self.observed_magnetic_list)
        answer = selection.ShowModal()
        for i in range(0, len(self.observed_magnetic_list)):
                if self.observed_magnetic_list[i].name == selection.obs_name:
                    self.obs_mag_data_for_rms = self.observed_magnetic_list[i].data

    def model_rms(self, xp):
        """CALCULATE RMS MISFIT OF OBSERVED VS CALCULATED"""
        if self.obs_gravity_data_for_rms != [] and self.calc_grav_switch is True:
            x = xp * 0.001
            y = self.predgz
            self.grav_rms_value, self.grav_residuals = model_stats.rms(self.obs_gravity_data_for_rms[:, 0],
                                                                       self.obs_gravity_data_for_rms[:, 1], x, y)
        else:
            pass

        if self.obs_mag_data_for_rms != [] and self.calc_mag_switch is True:
            x = self.xp * 0.001
            y = self.prednt
            self.mag_rms_value, self.mag_residuals = model_stats.rms(self.obs_mag_data_for_rms[:, 0],
                                                                     self.obs_mag_data_for_rms[:, 1], x, y)
        else:
            pass

    # LAYER ATTRIBUTE TABLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def open_attribute_table(self, event):
        attribute_table = AttributeEditor(self, -1, 'attribute editor', self.tree_items, self.densities,
                                          self.reference_densities, self.susceptibilities, self.angle_a,
                                          self.angle_b, self.layer_colors)
        attribute_table.Show(True)

    def attribute_set(self, new_tree_items, densities, reference_densities, susceptibilities, angle_a, angle_b,
                      layer_colors):
        """UPDATE GMG ATTRIBUTES WITH NEW ATTRIBUTES FROM THE ATTRIBUTE TABLE"""

        # UPDATE MAIN FRAME TREE LIST
        current_tree_items = self.tree.GetRootItem().GetChildren()
        for i in range(0, len(self.tree_items) - 1):
            new_label = new_tree_items[i]
            self.tree.SetItemText(current_tree_items[i], new_tree_items[i + 1])

        # UPDATE MAIN FRAME ATTRIBUTES
        self.densities = densities
        self.reference_densities = reference_densities
        self.susceptibilities = susceptibilities
        self.angle_a = angle_a
        self.angle_b = angle_b
        self.layer_colors = layer_colors

        # UPDATE GMG STATE
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    # LIVE GRAPHICS UPDATES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_layer_data(self):
        """LOAD LAYER DATA INTO GMG"""

        # UPDATE MODEL CANVAS (mcanvas) LIMITS
        xmin, xmax = self.model_frame.get_xlim()
        self.topo_frame.set_xlim(xmin, xmax)
        self.gravity_frame.set_xlim(xmin, xmax)
        self.magnetic_frame.set_xlim(xmin, xmax)

        # SET LISTS WITH UPDATED LAYER DATA
        self.plotx_list[self.currently_active_layer_id] = self.current_x_nodes
        self.ploty_list[self.currently_active_layer_id] = self.ploty

        # RESET LISTS (UPDATED BY THE CALCULATE_GRAVITY FUNC)
        self.gravity_polygons = []
        self.mag_polygons = []

        # REMOVE EXISTING DATA
        if len(self.model_frame.lines) >= 1:
            while len(self.model_frame.lines) > 0:
                self.model_frame.lines[0].remove()
            while len(self.model_frame.patches) > 0:
                self.model_frame.patches[0].remove()
            while len(self.model_frame.texts) > 0:
                self.model_frame.texts[0].remove()

        # UPDATE DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # FIRST CREATE THE POLYLINE DATA (THE BOTTOM LINE OF THE LAYER POLYGON - DONE FIRST SO THE WHOLE POLYGON
        # ISN'T PASSED TO SELF.POLYPLOTS)
        for i in range(0, self.total_layer_count + 1):
            # CREATE THE LAYER POLYGONS TO PASS TO SELF.POLYGONS AND ONTO THE TALWANI ALGORITHM

            # FIRST SET UP XY DATA; IF THE LAYER IS BELOW LAYER 1 THEN ATTACH THE ABOVE LAYER TO COMPLETE POLYGON;
            # ELSE USE TOP LAYER CHECK FOR LAYER MODE AND FIND LAST LAYER TO MAKE POLYGON
            if i >= 1 and self.layer_lock_list[i] == 0:
                for layers in range(i, 0, -1):
                    if self.layer_lock_list[layers - 1] == 0:
                        self.last_layer = layers - 1

                        # NOW APPEND NODES FOR BOUNDARY CONDITIONS (CONTINUOUS SLAB)
                        plotx = np.array(self.plotx_list[i])
                        ploty = np.array(self.ploty_list[i])

                        # SET PADDING NODES TO DEPTH EQUAL TO MODEL LIMIT NODES TO CREATE SLAB'
                        ploty[0] = ploty[1]
                        ploty[-1] = ploty[-2]
                        self.plotx_list[i], self.ploty_list[i] = plotx, ploty

                        # ADD NODES FROM ABOVE LAYER TO COMPETE POLYGON'
                        layer_above_x = np.array(self.plotx_list[self.last_layer])[::-1]
                        layer_above_y = np.array(self.ploty_list[self.last_layer])[::-1]

                        # CREATE POLYGON
                        self.plotx_polygon = np.append(np.array(plotx), np.array(layer_above_x))
                        self.ploty_polygon = np.append(np.array(ploty), np.array(layer_above_y))
                        break
            else:
                # IF THE LAYER IS A SIMPLE 'FLOATING LAYER'
                self.plotx_polygon = np.array(self.plotx_list[i])
                self.ploty_polygon = np.array(self.ploty_list[i])

            if self.densities[i] != 0 and self.absolute_densities is True:
                next_color = self.colormap.to_rgba(0.001 * self.densities[i] - 0.001 * self.reference_densities[i])
            elif self.densities[i] != 0:
                next_color = self.colormap.to_rgba(0.001 * self.densities[i])
            else:
                next_color = self.colormap.to_rgba(0)

            # CREATE LAYER POLYGON FILLS
            self.polygon_fills[i] = self.model_frame.fill(self.plotx_polygon, self.ploty_polygon, color=next_color,
                                                      alpha=self.layer_transparency, closed=True, linewidth=None,
                                                      ec=None)

            self.gravity_polygons.append(zip(self.plotx_polygon, self.ploty_polygon))

        # MODEL LAYERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CREATE LAYER LINES
        for i in range(0, self.total_layer_count + 1):
            self.layer_lines[i] = self.model_frame.plot(self.plotx_list[i], self.ploty_list[i],
                                                        color=self.layer_colors[i], linewidth=1.0, alpha=0.5)
        # CURRENT LAYER LINE
        self.currently_active_layer, = self.model_frame.plot(self.current_x_nodes, self.current_y_nodes, marker='o',
                                               color=self.layer_colors[self.currently_active_layer_id], linewidth=1.0,
                                                             alpha=0.5)

        # # FAULTS LAYERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # # CREATE FAULT LINES
        # for i in range(0, self.fault_counter + 1):
        #     self.faults[i] = self.model_frame.plot(self.fault_plotx_list[i], self.fault_ploty_list[i],
        #                                             color=self.fault_colors[i], linewidth=1.0, alpha=1.0)
        # # CURRENT FAULT LINE
        # self.fault_polyline, = self.model_frame.plot(self.fault_plotx, self.fault_ploty, marker='o',
        #                                          color=self.fault_colors[self.current_fault_index],
        #                                          linewidth=1.0, alpha=0.5)

        # UPDATE GMG
        self.model_frame.set_aspect(self.model_aspect)
        self.display_info()
        self.run_algorithms()
        self.draw()

    def update_layer_data(self):
        """UPDATE PROGRAM GRAPHICS AFTER A CHANGE IS MADE - A.K.A REDRAW EVERYTHING"""

        # UPDATE FRAME LIMITS
        xmin, xmax = self.model_frame.get_xlim()
        if self.topo_frame:
            self.topo_frame.set_xlim(xmin, xmax)
            self.topo_d_frame.set_xlim(xmin, xmax)
        if self.gravity_frame:
            self.gravity_frame.set_xlim(xmin, xmax)
            self.gravity_d_frame.set_xlim(xmin, xmax)
        if self.magnetic_frame:
            self.magnetic_frame.set_xlim(xmin, xmax)
            self.magnetic_d_frame.set_xlim(xmin, xmax)

        if self.fault_picking_switch is True:
            # GMG IS IN FAULT MODE

            # UPDATE FAULT LINE GRAPHICS
            for i in range(0, self.fault_counter):
                try:
                    self.faults[i][0].set_xdata(self.fault_x_coords_list[i])
                    self.faults[i][0].set_ydata(self.fault_y_coords_list[i])
                except IndexError:
                    # USED FOR THE CASE WHEN ONLY 1 FAULT EXISTS
                    pass

            # UPDATE CURRENT PLOT GRAPHICS
            self.faultline.set_data(self.fault_x_coords_list[self.current_fault_index],
                                    self.fault_y_coords_list[self.current_fault_index])
        else:
            # GMG IS IN LAYER MODE
            # UPDATE PLOT LISTS WITH LATEST EDIT
            self.layer_list[self.currently_active_layer_id].x_nodes = self.current_x_nodes
            self.layer_list[self.currently_active_layer_id].y_nodes = self.current_y_nodes
            print("total layer count = %s") % self.total_layer_count
            print("")

            # RESET LISTS (USED AS INPUT FOR THE GRAV/MAG ALGORITHMS)
            self.gravity_polygons = []
            self.mag_polygons = []

            # CREATE UPDATED POLYGON XYs -------------------------------------------------------------------------------
            # FIRST CREATE THE POLYLINE DATA (THE BOTTOM LINE OF THE LAYER POLYGON - (THIS DONE FIRST SO THE WHOLE
            # POLYGON ISN'T PASSED TO SELF.POLYPLOTS)
            for i in range(0, self.total_layer_count + 1):
                print("i = %s") % i

                # CREATE THE LAYER POLYGONS TO PASS TO SELF.POLYGONS AND ONTO THE GRAV/MAG ALGORITHMS
                # FIRST SET UP XY DATA; IF LAYER IS BELOW LAYER 0 THEN ATTACH THE ABOVE LAYER TO COMPLETE THE POLYGON;
                # ELSE USE TOP LAYER CHECK FOR 'FIXED' LAYER MODE AND FIND LAST LAYER TO MAKE POLYGON
                print self.layer_list[i].name
                print self.layer_list[i].id
                print self.layer_list[i].type
                print self.layer_list[i].color
                print self.layer_list[i].x_nodes
                print self.layer_list[i].y_nodes
                print self.layer_list[i].node_mpl_actor
                print self.layer_list[i].polygon_mpl_actor
                print("")
                if i >= 1 and self.layer_list[i].type == 'fixed':
                    for layers in range(i, 0, -1):
                        if self.layer_list[i-1].type == 'fixed':
                            self.last_layer = layers - 1

                            # NOW APPEND NODES FOR BOUNDARY CONDITIONS (CONTINUOUS SLAB)
                            plotx = np.array(self.layer_list[i].x_nodes)
                            ploty = np.array(self.layer_list[i].y_nodes)

                            # SET PADDING NODES TO DEPTH EQUAL TO MODEL LIMIT NODES TO CREATE SLAB
                            ploty[0], ploty[-1] = ploty[1], ploty[-2]
                            self.layer_list[i].x_nodes = plotx
                            self.layer_list[i].y_nodes = ploty

                            # ADD NODES FROM ABOVE LAYER TO COMPETE POLYGON
                            layer_above_x = np.array(self.layer_list[i-1].x_nodes)[::-1]
                            layer_above_y = np.array(self.layer_list[i-1].y_nodes)[::-1]

                            # CREATE POLYGON
                            self.plotx_polygon = np.append(np.array(plotx), np.array(layer_above_x))
                            self.ploty_polygon = np.append(np.array(ploty), np.array(layer_above_y))
                            break
                else:
                    # IF THE LAYER IS A SIMPLE 'FLOATING LAYER'
                    self.plotx_polygon = np.array(self.layer_list[i].x_nodes)
                    self.ploty_polygon = np.array(self.layer_list[i].y_nodes)

                # APPEND LAYER TO THE LIST OF LAYERS SUPPLIED TO THE GRAVITY ALGORITHM
                self.gravity_polygons.append(zip(self.plotx_polygon, self.ploty_polygon))
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # UPDATE LAYER POLYGONS AND LINES
            for i in range(0, self.total_layer_count + 1):

                # SET POLYGON FILL COLOR BASED ON DENSITY
                if self.layer_list[i].density != 0.0 and self.layer_list[i].reference_density is True:
                    # DETERMINE DENSITY CONTRAST FROM (DENSITY - REF DENSITY)
                    next_color = self.colormap.to_rgba(0.001 * self.layer_list[i].density -
                                                       0.001 * self.layer_list[i].reference_density)
                elif self.layer_list[i].density != 0.0:
                    # NO REF DENSITY, SO JUST USE DENSITY VALUE
                    next_color = self.colormap.to_rgba(0.001 * self.layer_list[i].density)
                else:
                    # NO DENSITY HAS BEEN SET SO LEAVE BLANK
                    next_color = self.colormap.to_rgba(0.)

                # # UPDATE POLYGON XY AND COLOR FILL
                self.layer_list[i].polygon_mpl_actor[0].set_xy(self.gravity_polygons[i])
                self.layer_list[i].polygon_mpl_actor[0].set_color(next_color)
                #
                # # UPDATE LAYER LINES
                self.layer_list[i].node_mpl_actor[0].set_xdata(self.layer_list[i].x_nodes)
                self.layer_list[i].node_mpl_actor[0].set_ydata(self.layer_list[i].y_nodes)
            # ----------------------------------------------------------------------------------------------------------

            # UPDATE CURRENTLY ACTIVE LAYER LINE AND NODES
            self.currently_active_layer.set_xdata(self.layer_list[self.currently_active_layer_id].x_nodes)
            self.currently_active_layer.set_ydata(self.layer_list[self.currently_active_layer_id].y_nodes)
            self.currently_active_layer.set_color(self.layer_list[self.currently_active_layer_id].color)

        # DRAW CANVAS FEATURES
        self.model_frame.set_aspect(self.model_aspect)
        self.grav_frame_aspect = ((self.gravity_frame.get_xlim()[1] - self.gravity_frame.get_xlim()[0]) /
                                  (self.gravity_frame.get_ylim()[1] - self.gravity_frame.get_ylim()[0]))
        # UPDATE INFO
        self.display_info()

        # CONTENT HAS NOT BEEN SAVED SINCE LAST MODIFICATION
        self.model_saved = False

        # UPDATE GMG GRAPHICS
        self.draw()

    def run_algorithms(self):
        """
        RUN POTENTIAL FIELD CALCULATION ALGORITHMS
        """
        # --------------------------------------------------------------------------------------------------------------
        # CALCULATE TOPOGRAPHY - :FUTURE: PREDICTED TOPOGRAPHY FROM ISOSTATIC FUNC
        self.pred_topo = np.zeros_like(self.xp)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # CALCULATE GRAVITY
        # CREATE ARRAY OF DENSITY CONTRASTS TO BE PASSED TO THE BOTT ALGORITHM
        self.density_contrasts = np.zeros(len(self.densities))
        for i in range(len(self.densities)):
            self.density_contrasts[i] = (self.densities[i] - self.reference_densities[i])

        # ZIP POLYGONS WITH DENSITY CONTRASTS AND PASS TO BOTT
        if self.gravity_polygons and self.calc_grav_switch is True:
            # SELECT ONLY THOSE LAYERS THAT ARE CHECKED
            polygons_to_use = []
            densities_to_use = []
            for layer in range(len(self.densities)):
                if self.layers_calculation_switch[layer] == 1:
                    polygons_to_use.append(self.gravity_polygons[layer])
                    # NB: NODES ARE INPUT LEFT TO RIGHT SO WE MUST MULTIPLY BY -1 TO PRODUCE THE CORRECT SIGN AT OUTPUT
                    densities_to_use.append(self.density_contrasts[layer] * -1)

            # PASS TO BOTT ALGORITHM
            polys = []
            for p, d in zip(polygons_to_use, densities_to_use):
                polys.append(Polygon(1000 * np.array(p), {'density': d}))
            self.predgz = bott.Gz(self.xp, self.zp, polys)
        else:
            self.predgz = np.zeros_like(self.xp)
        self.predplot.set_data(self.xp * 0.001, self.predgz)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # CALCULATE MAGNETICS
        # ZIP POLYGONS WITH SUSCEPT/STRIKES AND PASS TO TALWANI AND HEIRTZLER ALGORITHM
        if self.gravity_polygons and self.earth_field != 0.0 and self.calc_mag_switch is True:
            # SELECT ONLY THOSE LAYERS THAT ARE CHECKED
            polygons_to_use, susceptibilities_to_use = [], []
            for layer in range(0, len(self.gravity_polygons)):
                if self.layers_calculation_switch[layer] == 1:
                    polygons_to_use.append(self.gravity_polygons[layer])
                    # NB: NODES ARE INPUT LEFT TO RIGHT SO WE MUST MULTIPLY BY -1 TO PRODUCE THE CORRECT SIGN AT OUTPUT
                    susceptibilities_to_use.append(self.susceptibilities[layer] * -1)

            # PASS TO TALWANI&HEIRTZLER ALGORITHM
            self.polys = []
            for p, s, in zip(polygons_to_use, susceptibilities_to_use):
                self.polys.append(Polygon(1000. * np.array(p), {'susceptibility': s}))
            self.prednt = talwani_and_heirtzler.ntz(self.xp, self.zp, self.polys, self.model_azimuth, self.earth_field,
                                                    self.angle_a, self.angle_b, self.mag_observation_elv)
        else:
            self.prednt = np.zeros_like(self.xp)

        # PLOT MAGNETIC PROFILE
        self.prednt_plot.set_data(self.xp * 0.001, self.prednt)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # UPDATE RMS VALUES

        # RUN THE RMS CALC CODE
        self.model_rms(self.xp)

        # 2SET GRAVITY RMS
        if self.obs_gravity_data_for_rms != [] and self.observed_gravity_switch is True and \
                self.calc_grav_switch is True and self.predgz != []:
            self.grav_rms_plot.set_data(self.grav_residuals[:, 0], self.grav_residuals[:, 1])
        else:
            pass

        # SET MAGNETIC RMS
        if self.obs_mag_data_for_rms != [] and self.mag_obs_switch is True and self.calc_mag_switch is True \
                and self.prednt != []:
            self.mag_rms_plot.set_data(self.mag_residuals[:, 0], self.mag_residuals[:, 1])
        else:
            pass
        # --------------------------------------------------------------------------------------------------------------

        # SET FRAME X AND Y LIMITS
        self.set_frame_limits()

        # AFTER RUNNING ALGORITHMS, SET MODEL AS UNSAVED
        self.model_saved = False

        # UPDATE GMG GRAPHICS
        self.draw()

    def set_frame_limits(self):
        """SET FRAME X AND Y LIMITS"""
        # --------------------------------------------------------------------------------------------------------------
        # SET GRAVITY DISPLAY BOX LIMITS
        if self.observed_gravity_switch is True and self.grav_residuals != []:
            # CREATE EMPTY LIST
            ymin_list = []
            ymax_list = []

            # APPEND OBSERVED MIN AND MAX
            for i in range(len(self.observed_gravity_list)):
                if self.observed_gravity_list[i] is not None:
                    ymin_list.append(self.observed_gravity_list[i].data[:, 1].min() - 2.0)
                    ymax_list.append(self.observed_gravity_list[i].data[:, 1].max() + 2.0)

            # APPEND PREDICTED GRAVITY ANOMALY
            ymin_list.append(self.predgz.min())
            ymax_list.append(self.predgz.max())

            # APPEND RMS GRAVITY ANOMALY
            ymin_list.append(self.grav_residuals.min() - 2.0)
            ymax_list.append(self.grav_residuals.max() + 2.0)

            # SET YMIN AND YMAX
            ymin = min(ymin_list)
            ymax = max(ymax_list)

        elif self.observed_gravity_switch is True:
            # CREATE EMPTY LIST
            ymin_list = []
            ymax_list = []

            # APPEND OBSERVED MIN AND MAX
            for i in range(len(self.observed_gravity_list)):
                if self.observed_gravity_list[i] is not None:
                    ymin_list.append(self.observed_gravity_list[i].data[:, 1].min() - 2.0)
                    ymax_list.append(self.observed_gravity_list[i].data[:, 1].max() + 2.0)

            # APPEND PREDICTED GRAVITY ANOMALY
            if self.predgz != []:
                ymin_list.append(self.predgz.min() - 2.0)
                ymax_list.append(self.predgz.max() + 2.0)

            # SET YMIN AND YMAX
            ymin = min(ymin_list)
            ymax = max(ymax_list)

        else:
            ymin = self.predgz.min() - 2.0
            ymax = self.predgz.max() + 2.0


        if self.gravity_frame is not None:
            self.gravity_frame.set_ylim(ymin, ymax)
        # --------------------------------------------------------------------------------------------------------------
        # SET DERIVATIVE Y-AXIS LIMITS

        # CREATE EMPTY LIST
        ymin_list = [-1]
        ymax_list = [1]
        for i in range(len(self.observed_gravity_list)):
            if self.observed_gravity_list[i].type == str('derivative'):
                ymin_list.append(self.observed_gravity_list[i].data[:, 1].min() - 0.1)
                ymax_list.append(self.observed_gravity_list[i].data[:, 1].max() + 0.1)
        if self.gravity_frame is not None:
            self.gravity_d_frame.set_ylim(min(ymin_list), max(ymax_list))

        # --------------------------------------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------------------------------------
        # SET MAGNETIC DISPLAY BOX LIMITS
        if self.observed_magnetic_switch is True and self.mag_residuals != []:
            # CREATE EMPTY LIST
            ymin_list = []
            ymax_list = []

            # APPEND OBSERVED MIN AND MAX
            for i in range(len(self.observed_magnetic_list)):
                if self.observed_magnetic_list[i] is not None:
                    ymin_list.append(self.observed_magnetic_list[i].data[:, 1].min() - 2.0)
                    ymax_list.append(self.observed_magnetic_list[i].data[:, 1].max() + 2.0)

            # APPEND PREDICTED GRAVITY ANOMALY
            ymin_list.append(self.prednt.min())
            ymax_list.append(self.prednt.max())

            # APPEND RMS GRAVITY ANOMALY
            ymin_list.append(self.mag_residuals.min() - 2.0)
            ymax_list.append(self.mag_residuals.max() + 2.0)

            # SET YMIN AND YMAX
            ymin = min(ymin_list)
            ymax = max(ymax_list)

        elif self.observed_magnetic_switch is True:
            # CREATE EMPTY LIST
            ymin_list = []
            ymax_list = []

            # APPEND OBSERVED MIN AND MAX
            for i in range(len(self.observed_magnetic_list)):
                if self.observed_magnetic_list[i] is not None:
                    ymin_list.append(self.observed_magnetic_list[i].data[:, 1].min() - 2.0)
                    ymax_list.append(self.observed_magnetic_list[i].data[:, 1].max() + 2.0)

            # APPEND PREDICTED GRAVITY ANOMALY
            ymin_list.append(self.prednt.min() - 2.0)
            ymax_list.append(self.prednt.max() + 2.0)

            # SET YMIN AND YMAX
            ymin = min(ymin_list)
            ymax = max(ymax_list)

        else:
            # APPEND PREDICTED GRAVITY ANOMALY
            ymin = self.prednt.min() - 2.0
            ymax = self.prednt.max() + 2.0

        if self.magnetic_frame is not None:
            self.magnetic_frame.set_ylim(ymin, ymax)

        # SET DERIVATIVE Y-AXIS LIMITS
        # --------------------------------------------------------------------------------------------------------------
        # CREATE EMPTY LIST
        ymin_list = []
        ymax_list = []
        for i in range(len(self.observed_magnetic_list)):
            if self.observed_magnetic_list[i].type == 'derivative':
                ymin_list.append(self.observed_magnetic_list[i].data[:, 1].min() - 0.1)
                ymax_list.append(self.observed_magnetic_list[i].data[:, 1].max() + 0.1)
        self.magnetic_d_frame.set_ylim(ymin, ymax)
        # --------------------------------------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------------------------------------

        # UPDATE GMG GRAPHICS
        self.draw()

    # EXTERNAL FIGURE CONSTRUCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def plot_model(self, event):
        """CREATE EXTERNAL FIGURE OF MODEL USING INBUILT FIGURE CONSTRUCTION TOOL"""

        # GET PLOTTING PARAMETERS FROM DIALOG BOX
        self.set_values = PlotSettingsDialog(self, -1, 'Set figure parameters', self.model_aspect, self.grav_frame_aspect)
        self.set_values.Show(True)

    def draw_model(self):
        self.file_path = self.set_values.file_path
        self.file_type = self.set_values.file_type
        self.fs = self.set_values.fs
        self.ms = self.set_values.ms
        self.lw = self.set_values.lw
        self.font_type = self.set_values.font_type
        self.aspect_ratio = self.set_values.aspect_ratio
        self.use_tight_layout = self.set_values.use_tight_layout
        self.poly_alpha = self.set_values.poly_alpha
        self.draw_polygons = self.set_values.draw_polygons
        self.draw_layers = self.set_values.draw_layers
        self.floating_layers = self.set_values.draw_floating_layers
        self.draw_colorbar = self.set_values.draw_colorbar
        self.draw_wells = self.set_values.draw_wells
        self.draw_xy_data = self.set_values.draw_xy_data
        self.well_fs = self.set_values.well_fs
        self.well_line_width = self.set_values.well_line_width
        self.draw_faults = self.set_values.draw_faults
        self.xy_size = self.set_values.xy_size
        self.xy_color = self.set_values.xy_color
        self.colorbar_x = self.set_values.colorbar_x
        self.colorbar_y = self.set_values.colorbar_y
        self.colorbar_size_x = self.set_values.colorbar_size_x
        self.colorbar_size_y = self.set_values.colorbar_size_y
        self.layer_line_width = self.set_values.layer_line_width
        self.layer_alpha = self.set_values.layer_alpha
        self.grav_y_min = self.set_values.grav_frame_min
        self.grav_y_max = self.set_values.grav_frame_max

        # GET FIGURE DIMENSIONS
        xmin, xmax = self.model_frame.get_xlim()
        ymin, ymax = self.model_frame.get_ylim()
        area = np.array([xmin, xmax, ymin, ymax])

        # RUN PLOT MODEL CODE
        fig_plot = plot_model.plot_fig(self.file_path, self.file_type, area, self.xp, self.obs_topo, self.pred_topo,
                                       self.obs_grav, self.predgz, self.obs_mag, self.prednt, self.total_layer_count,
                                       self.layer_lock_list, self.plotx_list, self.ploty_list, self.densities,
                                       self.absolute_densities, self.reference_densities, self.segy_plot_list,
                                       self.well_list, self.well_name_list, self.topo_frame, self.gravity_frame,
                                       self.magnetic_frame, self.aspect_ratio, self.use_tight_layout, self.poly_alpha,
                                       self.fs, self.ms, self.lw, self.font_type, self.layer_colors, self.draw_polygons,
                                       self.draw_layers, self.floating_layers, self.draw_colorbar, self.draw_xy_data,
                                       self.xy_size, self.xy_color, self.colorbar_x, self.colorbar_y,
                                       self.colorbar_size_x, self.colorbar_size_y, self.layer_line_width,
                                       self.layer_alpha, self.grav_rms_value, self.mag_rms_value, self.grav_y_min,
                                       self.grav_y_max, self.xy_data_list_save, self.draw_wells, self.wells, self.well_fs,
                                       self.well_line_width, self.draw_faults, self.faults)
        del fig_plot

        # # IF ON A LINUX SYSTEM OPEN THE FIGURE WITH PDF VIEWER
        # if sys.platform == 'linux2':
        #     subprocess.call(["xdg-open", self.file_path])
        # # IF ON A macOS SYSTEM OPEN THE FIGURE WITH PDF VIEWER
        # elif sys.platform == 'darwin':
        #     os.open(self.file_path)

        # UPDATE GMG
        self.update_layer_data()
        self.draw()

        return

    # DOCUMENTATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def open_documentation(self, event):
        """OPEN DOCUMENTATION HTML"""

        self.doc_dir = os.path.dirname(os.path.abspath(__file__)).split('/')
        doc_url = self.doc_dir[0] + '/' + self.doc_dir[1] + '/' + self.doc_dir[2] + '/' + self.doc_dir[3] + \
                  '/docs/html/gmg_documentation.html'

        if platform == "linux" or platform == "linux2":
            # LINUX
            webbrowser.open_new(doc_url)
        elif platform == "darwin":
            # OS X
            client = webbrowser.get("open -a /Applications/Safari.app %s")
            client.open(doc_url)
        elif platform == "win32":
            # WINDOWS
            webbrowser.open_new(doc_url)

    def about_gmg(self, event):
        """ SHOW SOFTWARE INFORMATION"""
        about = [
            "GMG is an Open Source Graphical User Interface (GUI) designed principally for modelling 2D potential "
            "field (gravity and magnetic) profiles. The software also includes functions for loading XY data, "
            "seismic reflection SEGY data and exploration well horizons. The software therefore provides an "
            "integrated geological/geophysical interpretation package. It is anticipated that GMG will also be "
            "useful for teaching purposes. \n \n"
            "Data I/O is made as simple as possible using space delimited ASCII text files. \n \n"
            "The project was instigated after failing to find an adequate open source option (in which the source "
            "code can be viewed and modified by the user) for performing 2D geophysical modeling tasks.           "
            "Inspiration came from fatiando a terra and GMT. \n \n"
            "GMG was initially developed at the University of Oxford 2014-2017. \n \n"
            "B. Tozer"]

        dlg = wx.MessageDialog(self, about[0], "About", wx.OK | wx.ICON_INFORMATION)
        result = dlg.ShowModal()
        dlg.Destroy()

    def legal(self, event):
        """ SHOW LICENCE"""
        licence = ["Copyright 2015-2017 Brook Tozer \n\nRedistribution and use in source and binary forms, with or "
                   "without modification, are permitted provided that the following conditions are met: \n \n"
                   "1. Redistributions of source code must retain the above copyright notice, this list of conditions "
                   "and the following disclaimer. \n\n2. Redistributions in binary form must reproduce the above "
                   "copyright notice, this list of conditions and the following disclaimer in the documentation and/or "
                   "other materials provided with the distribution. \n\n3. Neither the name of the copyright holder "
                   "nor the names of its contributors may be used to endorse or promote products  derived from this "
                   "software without specific prior written permission. \n\nTHIS SOFTWARE IS PROVIDED BY THE "
                   "COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT "
                   "NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE "
                   "DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, "
                   "INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, "
                   "PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS "
                   "INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,"
                   " OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, "
                   "EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."]

        dlg = wx.MessageDialog(self, licence[0], "BSD-3-Clause Licence", wx.OK | wx.ICON_INFORMATION)
        result = dlg.ShowModal()
        dlg.Destroy()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXIT FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def exit(self, event):
        """ SHUTDOWN APP (FROM FILE MENU)"""
        dlg = wx.MessageDialog(self, "Do you really want to exit", "Confirm Exit", wx.OK | wx.CANCEL | wx.ICON_QUESTION)
        result = dlg.ShowModal()
        if result == wx.ID_OK:
            wx.GetApp().ExitMainLoop()

    def on_close_button(self, event):
        """ SHUTDOWN APP (X BUTTON)"""
        dlg = wx.MessageDialog(self, "Do you really want to exit", "Confirm Exit", wx.OK | wx.CANCEL | wx.ICON_QUESTION)
        result = dlg.ShowModal()
        if result == wx.ID_OK:
            wx.GetApp().ExitMainLoop()

    # FUTURE MODULES (IN PROCESS)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_error(self, value):
        pass
        # self.error = value
        # self.update_layer_data()
        # self.run_algorithms()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GMG OBJECT CLASSES GO HERE


class Layer:
    """
    CREATE A MODEL LAYER OBJECT.
    THE LAYER WILL BE STORED IN THE gmg.layers_list LIST
    """
    def __init__(self):
        self.id = None  # THE LAYER NUMBER
        self.name = None  # THE LAYER NAME
        self.layer_type = None  # EITHER str('fixed') OR str('floating')
        self.node_pml_actor = None
        self.polygon_mpl_actor = None
        self.poly_plot = None
        self.x_nodes = []
        self.y_nodes = []
        self.density = 0.
        self.reference_density = 0.
        self.susceptibility = 0.
        self.angle_a = 0.
        self.angle_b = 0.
        self.color = 'k'
        self.layer_transparency = 0.4
        self.pinch_node_list = [[], []]
        self.pinch_count = 0
        self.layers_calculation_switch = True  # SWITCH TO DICTATE IF THE LAYER IS INCLUDED IN THE CURRENT CALC


class ObservedData:
    def __init__(self):
        """GENERIC CLASS FOR AN OBSERVATIONAL DATA OBJECT"""
        self.id = None  # OBJECT ID VALUE FOR wx
        self.name = None  # THE NAME ASSIGNED TO THE DATA (str)
        self.color = None  # THE COLOR USED FOR PLOTTING THE DATA (str)
        self.type = None  # THE KIND OF DATA (str: 'observed', 'filtered', 'derivative')
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.mpl_actor = None  # THE PML ACTOR ELEMENT (PLOTTING OBJECT)


class ObservedOutcropData:
    """GENERIC CLASS FOR AN OBSERVATIONAL OUTCROP DATA OBJECT"""
    def __init__(self):
        self.id = None
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.name = None  # THE NAME ASSIGNED TO THE DATA (str)
        self.color = None  # THE COLOR USED FOR PLOTTING THE DATA (str)
        self.textsize = 2  # THE SIZE OF TEXT LABELS

class SegyData:
    def __init__(self):
        self.file = None
        self.name = None  # THE NAME ASSIGNED TO THE DATA (str)
        self.dimensions = None
        self.mpl_actor = None
        self.axis = None
        self.color_map = cm.gray
        
        
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.gain_positive = 4.0
        self.gain_neg = -self.gain_positive
        self.segy_show = False
        self.plot_list = None


class ObservedWellData:
    def __init__(self):
        self.name = None  # NAME OF WELL RECORD
        self.raw_record = None  # RAW TEXT RECORD
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.mpl_actor = None  # MPL ACTOR USED TO STORE MPL LINE
        self.mpl_actor_name = None  # MPL ACTOR USED TO STORE WELL NAME LABEL
        self.labels_list = []
        self.horizons_list = []
        self.text_size = 2


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GMG EXTERNAL DIALOG BOXES GO HERE


class NewModelDialog(wx.Dialog):
    """CREATE A NEW MODEL. RETURNS MODEL PARAMETERS AND CALCULATION INCREMENT"""
    def __init__(self, parent, id, title, m_x1=None, m_x2=None, m_z1=None, m_z2=None):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.MAXIMIZE_BOX | wx.OK | wx.CANCEL
                                                          | wx.BORDER_RAISED)

        self.floating_panel = wx.Panel(self, -1)
        self.set_button = False

        self.main_box = wx.BoxSizer(wx.HORIZONTAL)

        self.title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.title_text = wx.StaticText(self.floating_panel, -1, "Enter the region boundary for the new model\n"
                                                                 "and the x-increment for anomaly calculations:",
                                        style=wx.ALIGN_CENTRE_HORIZONTAL)
        self.title_sizer.Add(self.title_text)

        self.line1 = (wx.StaticLine(self.floating_panel), 0, wx.ALL | wx.EXPAND, 5)

        self.x_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.x_Text = wx.StaticText(self.floating_panel, -1, "X1, X2 (km)      ")
        self.new_x1 = wx.TextCtrl(self.floating_panel, -1, "0", size=(100, -1))
        self.new_x2 = wx.TextCtrl(self.floating_panel, -1, "0", size=(100, -1))
        self.x_sizer.AddMany([self.x_Text, self.new_x1, self.new_x2])
        self.z_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.z_Text = wx.StaticText(self.floating_panel, -1, "Z1, Z2 (km)      ")
        self.new_z1 = wx.TextCtrl(self.floating_panel, -1, "0", size=(100, -1))
        self.new_z2 = wx.TextCtrl(self.floating_panel, -1, "0", size=(100, -1))
        self.z_sizer.AddMany([self.z_Text, self.new_z1, self.new_z2])
        self.xp1_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.xp1_text = wx.StaticText(self.floating_panel, -1, "Calculation\nIncrement (km)")
        self.xp1_inc = wx.TextCtrl(self.floating_panel, -1, "0", size=(100, -1))
        self.xp1_sizer.AddMany([self.xp1_text, self.xp1_inc])

        self.line2 = (wx.StaticLine(self.floating_panel), 0, wx.ALL | wx.EXPAND, 5)

        self.create_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.spacing = wx.StaticText(self.floating_panel, -1, "                           ")
        self.b_create_button = wx.Button(self.floating_panel, -1, "Create Model")
        self.Bind(wx.EVT_BUTTON, self.create_button, self.b_create_button)
        self.create_button_sizer.AddMany([self.spacing, self.b_create_button])

        self.line3 = (wx.StaticLine(self.floating_panel), 0, wx.ALL | wx.EXPAND, 5)

        self.footer_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.footer_text = wx.StaticText(self.floating_panel, -1,
                                         "NB. these values can be modified during modelling,\n"
                                         "it may be beneficial to begin with a coarse spacing.")
        self.footer_sizer.Add(self.footer_text)

        #  POPULATE BOX WITH EXISTING VALUES
        if m_x1:
            self.new_x1.SetValue(str(m_x1))
        if m_x2:
            self.new_x2.SetValue(str(m_x2))
        if m_z1:
            self.new_z1.SetValue(str(m_z1))
        if m_z2:
            self.new_z2.SetValue(str(m_z2))

        self.sizer = wx.FlexGridSizer(cols=1, hgap=10, vgap=10)
        self.sizer.AddMany([self.title_text, self.line1, self.x_sizer, self.z_sizer, self.xp1_sizer,
                            self.line2, self.create_button_sizer, self.line3, self.footer_text])

        self.main_box.Add(self.sizer, proportion=1, flag=wx.ALL | wx.EXPAND, border=10)
        self.floating_panel.SetSizerAndFit(self.main_box)
        self.main_box.Fit(self)

    def create_button(self, event):
        self.x1 = float(self.new_x1.GetValue())
        self.x2 = float(self.new_x2.GetValue())
        self.z1 = float(self.new_z1.GetValue())
        self.z2 = float(self.new_z2.GetValue())
        self.xp_inc = float(self.xp1_inc.GetValue())
        self.set_button = True
        self.EndModal(1)


class LoadObservedDataFrame(wx.Frame):
    """LOAD OBSERVED DATA"""

    def __init__(self, parent, id, title, type):
        wx.Frame.__init__(self, None, wx.ID_ANY, 'Load observed data')

        floating_panel = wx.Panel(self, -1)

        # CREATE INSTANCE OF MAIN FRAME CLASS TO RECEIVE NEW ATTRIBUTES
        self.parent = parent

        # SET TYPE OF DATA BEING LOADED
        self.type = type

        # FILE PATH
        self.file_path_text = wx.TextCtrl(floating_panel, -1, value="input_file", size=(150, -1))
        self.b_file_path = wx.Button(floating_panel, -1, "File...")
        self.Bind(wx.EVT_BUTTON, self.file_path, self.b_file_path)

        # OBSERVED NAME
        self.observed_name_text = wx.StaticText(floating_panel, -1, "Name for observed data:")
        self.observed_name = wx.TextCtrl(floating_panel, -1, value="Observed_name", size=(150, -1))

        # COLOR
        self.colors = ['red', 'orange', 'green', 'blue', 'black', 'purple', 'pink']
        self.color_text = wx.StaticText(floating_panel, -1, "Display color:")
        self.color_picked = wx.ComboBox(floating_panel, -1, value='blue', choices=self.colors, size=(75, -1),
                                        style=wx.CB_DROPDOWN)

        # LOAD BUTTON
        self.b_load_button = wx.Button(floating_panel, -1, "Load")
        self.Bind(wx.EVT_BUTTON, self.load_button, self.b_load_button)

        # CANCEL BUTTON
        self.b_cancel_button = wx.Button(floating_panel, -1, "Cancel")
        self.Bind(wx.EVT_BUTTON, self.cancel_button, self.b_cancel_button)

        # SIZER
        sizer = wx.FlexGridSizer(cols=2, hgap=6, vgap=6)
        sizer.AddMany([self.file_path_text, self.b_file_path,
                       self.observed_name_text, self.observed_name,
                       self.color_text, self.color_picked,
                       self.b_load_button, self.b_cancel_button])
        floating_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def file_path(self, event):
        """RETURN FILE PATH TO LOAD"""
        self.open_file_dialog = wx.FileDialog(self, "Open Observed file", "", "", "All files (*.*)|*.*",
                                              wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if self.open_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # the user changed idea...

        #SET FILE PATH IN CLASS PANEL
        self.chosen_path = self.open_file_dialog.GetPath()
        self.file_path_text.SetValue(str(self.chosen_path))

    def load_button(self, event):
        """WHEN THE LOAD BUTTON IS PRESSED -> SET VALUES TO PASS TO MOULDER, THEN CLOSE WINDOW"""
        self.file_path = str(self.file_path_text.GetValue())
        self.observed_name = str(self.observed_name.GetValue())
        self.color_picked = str(self.color_picked.GetValue())

        if self.type == 'gravity':
            self.parent.open_obs_g()
        elif self.type == 'magnetics':
            self.parent.open_obs_m()
        elif self.type == 'topography':
            self.parent.open_obs_t()
        elif self.type == 'XY':
            self.parent.open_xy_data()
        elif self.type == 'outcrop':
            self.parent.open_outcrop_data()
        else:
            pass
        self.Destroy()

    def cancel_button(self, event):
        """USER CHANGED THEIR MIND"""
        self.Destroy()


class SeisDialog(wx.Dialog):
    """LOAD SEGY DATA"""

    def __init__(self, parent, id, title, area):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.MAXIMIZE_BOX
                                                          | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)
        self.area = area
        x1, x2, z1, z2 = 0.001 * np.array(self.area)
        self.x2 = str(x2)
        self.z2 = str(z2)

        self.x_start = wx.StaticText(input_panel, -1, "Model X start coordinate:")
        self.x_start_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.x_start_Text.SetInsertionPoint(0)
        self.x_end = wx.StaticText(input_panel, -1, "Model X end coordinate:")
        self.x_end_Text = wx.TextCtrl(input_panel, -1, value=self.x2, size=(100, -1))
        self.z_start = wx.StaticText(input_panel, -1, "Model Z start coordinate:")
        self.z_start_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.z_end = wx.StaticText(input_panel, -1, "Model Z end coordinate:")
        self.z_end_Text = wx.TextCtrl(input_panel, -1, value=self.z2, size=(100, -1))

        self.segy_name = wx.StaticText(input_panel, -1, "Segy Name:")
        self.segy_name_Text = wx.TextCtrl(input_panel, -1, size=(100, -1))

        sizer = wx.FlexGridSizer(cols=2, hgap=6, vgap=6)
        self.b_set_button = wx.Button(input_panel, -1, "Set Coordinates")
        self.Bind(wx.EVT_BUTTON, self.set_button, self.b_set_button)
        sizer.AddMany(
            [self.x_start, self.x_start_Text, self.x_end, self.x_end_Text, self.z_start, self.z_start_Text, self.z_end,
             self.z_end_Text, self.segy_name, self.segy_name_Text, self.b_set_button])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_button(self, event):
        self.sx1 = float(self.x_start_Text.GetValue())
        self.sx2 = float(self.x_end_Text.GetValue())
        self.sz1 = float(self.z_start_Text.GetValue())
        self.sz2 = float(self.z_end_Text.GetValue())
        self.dimensions = [self.sx1, self.sx2, self.sz1, self.sz2]
        self.segy_name_input = str(self.segy_name_Text.GetValue())
        self.EndModal(1)


class MagDialog(wx.Dialog):
    """SET MAGNETIC FIELD PARAMETERS"""

    def __init__(self, parent, id, title, area):
        wx.Dialog.__init__(self, parent, id, 'Set Magnetic Field', style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                                         | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        # OBSERVATION ELEVATION
        self.set_mag_observation_elv = wx.StaticText(input_panel, -1, "Observation\nElevation:")
        self.mag_observation_elv_text = wx.TextCtrl(input_panel, -1, "0", size=(75, -1))
        self.mag_observation_elv_units = wx.StaticText(input_panel, -1, "km")
        self.mag_observation_elv_text.SetInsertionPoint(0)

        # MODEL AZIMUTH
        self.set_model_azimuth = wx.StaticText(input_panel, -1, "Profile Azimuth\n(angle C):")
        self.model_azimuth_text = wx.TextCtrl(input_panel, -1, "0", size=(75, -1))
        self.model_azimuth_units = wx.StaticText(input_panel, -1, "Degrees")
        self.model_azimuth_text.SetInsertionPoint(0)

        # EARTHS FIELD
        self.set_earth_field = wx.StaticText(input_panel, -1, "Earth's Regional\nField Intensity:")
        self.earth_field_text = wx.TextCtrl(input_panel, -1, "0", size=(75, -1))
        self.earth_field_units = wx.StaticText(input_panel, -1, "nT")
        sizer = wx.FlexGridSizer(cols=3, hgap=7, vgap=7)
        self.b_set_button_mag = wx.Button(input_panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.set_button_mag, self.b_set_button_mag)
        sizer.AddMany([self.set_mag_observation_elv, self.mag_observation_elv_text, self.mag_observation_elv_units,
                       self.set_model_azimuth, self.model_azimuth_text, self.model_azimuth_units,
                       self.set_earth_field, self.earth_field_text, self.earth_field_units,
                       self.b_set_button_mag])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_button_mag(self, event):
        self.mag_observation_elv = float(self.mag_observation_elv_text.GetValue())
        self.model_azimuth = float(self.model_azimuth_text.GetValue())
        self.earth_field = float(self.earth_field_text.GetValue())
        self.EndModal(1)


class PinchDialog(wx.Dialog):
    """PINCH A LAYER TO ANOTHER LAYER"""

    def __init__(self, parent, id, title, plotx_list, ploty_list, i):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)
        self.plotx_list = plotx_list
        self.ploty_list = ploty_list
        self.currently_active_layer_id = i
        self.p_start = wx.StaticText(input_panel, -1, "Pinch From (km):")
        self.p_start_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.p_start_Text.SetInsertionPoint(0)
        self.p_end = wx.StaticText(input_panel, -1, "Pinch To (km):")
        self.p_end_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        sizer = wx.FlexGridSizer(cols=2, hgap=6, vgap=6)
        self.b_up_set_button = wx.Button(input_panel, -1, "Pinch up")
        self.b_down_set_button = wx.Button(input_panel, -1, "Pinch down")
        self.Bind(wx.EVT_BUTTON, self.up_set_button, self.b_up_set_button)
        self.Bind(wx.EVT_BUTTON, self.down_set_button, self.b_down_set_button)
        sizer.AddMany([self.p_start, self.p_start_Text, self.p_end, self.p_end_Text, self.b_up_set_button,
                       self.b_down_set_button])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def up_set_button(self, event):
        p_start = float(self.p_start_Text.GetValue())
        p_end = float(self.p_end_Text.GetValue())

        current_x = self.plotx_list[self.currently_active_layer_id]
        current_y = self.ploty_list[self.currently_active_layer_id]
        above_x = self.plotx_list[self.currently_active_layer_id - 1]
        above_y = self.ploty_list[self.currently_active_layer_id - 1]

        # REMOVE PINCHED NODES
        self.pinched_x = [tup for i, tup in enumerate(current_x) if current_x[i] < p_start or current_x[
            i] > p_end]  # PASS X VALUES IF THEY ARE OUTSIDE THE PINCH range
        self.pinched_y = [tup for i, tup in enumerate(current_y) if current_x[i] < p_start or current_x[
            i] > p_end]  # PASS Y VALUES IF THEY ARE OUTSIDE THE PINCH range

        # FIND NODES OF ABOVE LAYER WITHIN PINCH ZONE
        new_xs = [tup for i, tup in enumerate(above_x) if
                  p_start < above_x[i] < p_end]  # PASS X VALUES IF THEY ARE INSIDE THE PINCH range
        new_ys = [tup for i, tup in enumerate(above_y) if
                  p_start < above_x[i] < p_end]  # PASS Y VALUES IF THEY ARE INSIDE THE PINCH range

        # DETERMINE INSERT POSITION
        insert_point = 0
        for i in range(0, len(self.pinched_x)):
            if self.pinched_x[i] < p_start:
                insert_point = i + 1
        # INSERT NEW NODES
        self.pinched_x[insert_point:1] = new_xs
        self.pinched_y[insert_point:1] = new_ys

        self.EndModal(1)

    def down_set_button(self, event):
        p_start = float(self.p_start_Text.GetValue())
        p_end = float(self.p_end_Text.GetValue())

        current_x = self.plotx_list[self.currently_active_layer_id]
        current_y = self.ploty_list[self.currently_active_layer_id]
        below_x = self.plotx_list[self.currently_active_layer_id + 1]
        below_y = self.ploty_list[self.currently_active_layer_id + 1]

        # REMOVE PINCHED NODES
        self.pinched_x = [tup for i, tup in enumerate(current_x) if current_x[i] < p_start or current_x[
            i] > p_end]  # PASS X VALUES IF THEY ARE OUTSIDE THE PINCH range
        self.pinched_y = [tup for i, tup in enumerate(current_y) if current_x[i] < p_start or current_x[
            i] > p_end]  # PASS Y VALUES IF THEY ARE OUTSIDE THE PINCH range

        #  FIND NODES OF ABOVE LAYER WITHIN PINCH ZONE
        new_xs = [tup for i, tup in enumerate(below_x) if
                  p_start < below_x[i] < p_end]  # PASS X VALUES IF THEY ARE INSIDE THE PINCH range
        new_ys = [tup for i, tup in enumerate(below_y) if
                  p_start < below_x[i] < p_end]  # PASS Y VALUES IF THEY ARE INSIDE THE PINCH range

        # DETERMINE INSERT POSITION
        insert_point = 0
        for i in range(0, len(self.pinched_x)):
            if self.pinched_x[i] < p_start:
                insert_point = i + 1

        # INSERT NEW NODES
        self.pinched_x[insert_point:1] = new_xs
        self.pinched_y[insert_point:1] = new_ys

        self.EndModal(1)


class DepinchDialog(wx.Dialog):
    """DEPINCH A LAYER"""

    def __init__(self, parent, id, title, plotx_list, ploty_list, i, layer_count):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)
        self.plotx_list = plotx_list
        self.ploty_list = ploty_list
        self.currently_active_layer_id = i
        self.total_layer_count = layer_count
        self.p_start = wx.StaticText(input_panel, -1, "Pinch From:")
        self.p_start_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.p_start_Text.SetInsertionPoint(0)
        self.p_end = wx.StaticText(input_panel, -1, "Pinch Too")
        self.p_end_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        sizer = wx.FlexGridSizer(cols=2, hgap=6, vgap=6)
        self.b_depinch_button = wx.Button(input_panel, -1, "depinch")
        self.Bind(wx.EVT_BUTTON, self.depinch_button, self.b_depinch_button)
        sizer.AddMany([self.p_start, self.p_start_Text, self.p_end, self.p_end_Text, self.b_depinch_button])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def depinch_button(self, event):
        p_start = float(self.p_start_Text.GetValue())
        p_end = float(self.p_end_Text.GetValue())

        current_x = self.plotx_list[self.currently_active_layer_id]
        current_y = self.ploty_list[self.currently_active_layer_id]

        if self.currently_active_layer_id != 0:
            above_x = self.plotx_list[self.currently_active_layer_id - 1]
            above_y = self.ploty_list[self.currently_active_layer_id - 1]
        if self.currently_active_layer_id != self.total_layer_count:
            below_x = self.plotx_list[self.currently_active_layer_id + 1]
            below_y = self.ploty_list[self.currently_active_layer_id + 1]

        # REMOVE PINCHED NODES
        self.depinched_x = [tup for i, tup in enumerate(current_x) if current_x[i] < p_start or current_x[
            i] > p_end]  # PASS X VALUES IF THEY ARE OUTSIDE THE PINCH range
        self.depinched_y = [tup for i, tup in enumerate(current_y) if current_x[i] < p_start or current_x[
            i] > p_end]  # PASS Y VALUES IF THEY ARE OUTSIDE THE PINCH range

        # FIND NODES WITHIN PINCH ZONE - TRY ABOVE LAYER FIRST
        # PASS X VALUES IF THEY ARE INSIDE THE PINCH range
        new_xs = [tup for i, tup in enumerate(above_x) if p_start < above_x[i] < p_end]
        # PASS Y VALUES IF THEY ARE INSIDE THE PINCH range
        new_ys = [tup for i, tup in enumerate(above_y) if p_start < above_x[i] < p_end]
        if not new_xs:
            # PASS X VALUES IF THEY ARE INSIDE THE PINCH range
            new_xs = [tup for i, tup in enumerate(below_x) if p_start < below_x[i] < p_end]
            # PASS Y VALUES IF THEY ARE INSIDE THE PINCH range
            new_ys = [tup for i, tup in enumerate(below_y) if p_start < below_x[i] < p_end]

        new_ys_list = [y + 0.1 for y in new_ys]

        # DETERMINE INSERT POSITION
        insert_point = 0
        for i in range(0, len(self.depinched_x)):
            if self.depinched_x[i] < p_start:
                insert_point = i + 1

        # INSERT NEW NODES
        self.depinched_x[insert_point:1] = new_xs
        self.depinched_y[insert_point:1] = new_ys_list

        self.EndModal(1)


class MedianFilterDialog(wx.Dialog):
    """APPLY A MEDIAN FILTER TO OBSERVED DATA"""
    def __init__(self, parent, id, title, observed_list):
        wx.Dialog.__init__(self, parent, id, "Apply Median Filter", style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        # CREATE DATA LIST MENU
        self.observed_list = observed_list
        self.obs_combo_list_text = wx.StaticText(input_panel, -1, "Input Data:")
        self.obs_combo_list = wx.ComboBox(input_panel, id=-1, value="", choices=[])

        # POPULATE OBSERVED DATA LIST
        for i in range(len(observed_list)):
            if observed_list[i] is not None:
                self.obs_combo_list.Append(observed_list[i].name)

        # DEFINE FILTER LENGTH
        self.filter_window = wx.StaticText(input_panel, -1, "Filter Length (odd):")
        self.filter_window_text = wx.TextCtrl(input_panel, -1, "9")

        # DEFINE OUTPUT NAME
        self.output_name = wx.StaticText(input_panel, -1, "OUTPUT DATA NAME:")
        self.output_name_text = wx.TextCtrl(input_panel, -1)

        # DEFINE OUTPUT COLOR
        self.colors = ['red', 'orange', 'yellow', 'green', 'blue', 'grey', 'white', 'black']
        self.output_color = wx.StaticText(input_panel, -1, "OUTPUT DATA COLOR:")
        self.output_color_text = wx.ComboBox(input_panel, -1, value='green', choices=self.colors, size=(75, -1),
                                             style=wx.CB_DROPDOWN)

        # DEFINE SET BUTTON
        self.b_apply_filter = wx.Button(input_panel, -1, "Apply")
        self.Bind(wx.EVT_BUTTON, self.median_pass, self.b_apply_filter)

        # DEFINE SIZER'
        sizer = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        sizer.AddMany([self.obs_combo_list_text, self.obs_combo_list, self.output_name, self.output_name_text,
                       self.output_color, self.output_color_text, self.filter_window, self.filter_window_text,
                       self.b_apply_filter])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def median_pass(self, event):
        """APPLY A MEDIAN FILTER TO OBSERVED DATA"""
        self.obs_to_filter_name = str(self.obs_combo_list.GetValue())
        self.filter_length = int(self.filter_window_text.GetValue())
        self.output_name = str(self.output_name_text.GetValue())
        self.output_color = str(self.output_color_text.GetValue())

        for i in range(len(self.observed_list)):
            if self.observed_list[i].name == self.obs_to_filter_name:
                self.filter_input = self.observed_list[i].data

                self.filtered_output = np.zeros(shape=(len(self.filter_input), 2))
                self.filtered_output[:, 0] = self.filter_input[:, 0]
                self.filtered_output[:, 1] = signal.medfilt(self.filter_input[:, 1], self.filter_length)
                self.EndModal(1)


class HorizontalDerivative(wx.Dialog):
    """ESTIMATE THE HORIZONTAL DERIVATIVE OF THE OBSERVED DATA"""
    def __init__(self, parent, id, title, observed_list):
        wx.Dialog.__init__(self, parent, id, "Take Horizontal Derivative", style=wx.DEFAULT_DIALOG_STYLE |
                                                                                 wx.RESIZE_BORDER | wx.MAXIMIZE_BOX
                                                                                 | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        # CREATE DATA LIST MENU
        self.observed_list = observed_list
        self.obs_combo_list_text = wx.StaticText(input_panel, -1, "Input Data:")
        self.obs_combo_list = wx.ComboBox(input_panel, id=-1, value="", choices=[])

        # POPULATE OBSERVED DATA LIST
        for i in range(len(observed_list)):
            if observed_list[i] is not None:
                self.obs_combo_list.Append(observed_list[i].name)

        # DEFINE PREPROCESSING MEDIAN FILTER LENGTH
        self.filter_window = wx.StaticText(input_panel, -1, "Filter Length (odd):")
        self.filter_window_text = wx.TextCtrl(input_panel, -1, "3")

        # DEFINE PREPROCESSING MEDIAN FILTER LENGTH
        self.x_increment = wx.StaticText(input_panel, -1, "X increment (km):")
        self.x_increment_text = wx.TextCtrl(input_panel, -1, "1")

        # DEFINE OUTPUT NAME
        self.output_name = wx.StaticText(input_panel, -1, "Output data name:")
        self.output_name_text = wx.TextCtrl(input_panel, -1)

        # DEFINE OUTPUT COLOR
        self.colors = ['red', 'orange', 'yellow', 'green', 'blue', 'grey', 'white', 'black']
        self.output_color = wx.StaticText(input_panel, -1, "Output data color:")
        self.output_color_text = wx.ComboBox(input_panel, -1, value='red', choices=self.colors, size=(75, -1),
                                             style=wx.CB_DROPDOWN)

        # DEFINE SET BUTTON
        self.b_calculate = wx.Button(input_panel, -1, "Calculate")
        self.Bind(wx.EVT_BUTTON, self.calculate_derivative, self.b_calculate)

        # DEFINE SIZER
        sizer = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        sizer.AddMany([self.obs_combo_list_text, self.obs_combo_list, self.output_name, self.output_name_text,
                       self.output_color, self.output_color_text, self.x_increment, self.x_increment_text,
                       self.filter_window, self.filter_window_text, self.b_calculate])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def calculate_derivative(self, event):
        """TAKE HORIZONTAL DERIVATIVE OF SELECTED DATA"""

        # 0. PARSE THE USER DEFINED PARAMETERS
        self.obs_to_filter_name = str(self.obs_combo_list.GetValue())

        for i in range(len(self.observed_list)):
            if self.observed_list[i].name == self.obs_to_filter_name:
                self.obs_to_filter = self.observed_list[i].data

        self.filter_length = int(self.filter_window_text.GetValue())
        self.x_inc = float(self.x_increment_text.GetValue())
        self.output_name = str(self.output_name_text.GetValue())
        self.output_color = str(self.output_color_text.GetValue())

        # 1. INTERPOLATE THE DATA SO THE X INCREMENT IS EVENLY SPACED
        interpolate_func = ip.interp1d(self.obs_to_filter[:, 0],
                                       self.obs_to_filter[:, 1], kind='slinear')  # CREATE INTERPOLATION FUNCTION

        x_min = np.round(self.obs_to_filter[:, 0].min(), 0) + float(self.x_inc)
        x_max = np.round(self.obs_to_filter[:, 0].max(), 0) - float(self.x_inc)
        x_interp_values = np.arange(int(x_min), int(x_max) + 1, float(self.x_inc))  # X VALUES FOR INTERPOLATED ARRAY
        y_interp_values = np.array(interpolate_func(x_interp_values))  # Y VALUES FOR INTERPOLATED ARRAY

        input_interpolated = np.column_stack((x_interp_values, y_interp_values))

        # 2. APPLY A MEDIAN FILTER TO THE INTERPOLATED DATA
        self.filtered_output = np.zeros(shape=(len(input_interpolated), 2))
        self.filtered_output[:, 0] = input_interpolated[:, 0]
        self.filtered_output[:, 1] = signal.medfilt(input_interpolated[:, 1], self.filter_length)

        # 3. TAKE THE HORIZONTAL DERIVATIVE
        N = len(input_interpolated)
        self.deriv = np.zeros(shape=((len(input_interpolated) - 1), 2))

        for i in range(1, len(input_interpolated) - 1):
            self.deriv[i, 0] = input_interpolated[i, 0]

            # CALC USING FINITE DIFFERENCE METHOD
            self.deriv[i, 1] = abs((float(input_interpolated[i + 1, 1]) - float(input_interpolated[i - 1, 1]))) \
                          / (2. * float(self.x_inc))

        # SET FIRST AND LAST VALUES (WHICH ARE NOT CALCULATED BY FD METHOD) EQUAL TO SECOND AND SECOND-TO-LAST VALUES
        self.deriv[0, 0] = self.deriv[1, 0]-self.x_inc
        self.deriv[0, 1] = self.deriv[1, 1]
        self.deriv[-1, 0] = self.deriv[-2, 0]+self.x_inc
        self.deriv[-1, 1] = self.deriv[-2, 1]

        # CLOSE
        self.EndModal(1)


class SetBackgroundDensityDialog(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        # CHOOSE BACKGROUND DENSITY
        self.background_den_label_upper = wx.StaticText(input_panel, -1, "Reference density (g/cm^3):")
        self.background_density_text_input_upper = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.background_density_text_input_upper.SetInsertionPoint(0)
        # self.background_den_label_lower = wx.StaticText(input_panel, -1, "Lower Crust (g/cm^3):")
        # self.background_density_text_input_lower = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        # self.background_den_label_lid = wx.StaticText(input_panel, -1, "Mantle Lid (g/cm^3):")
        # self.background_density_text_input_lid = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))

        # DEFINE SET BUTTON
        self.b_apply_background_den = wx.Button(input_panel, -1, "Apply")
        self.Bind(wx.EVT_BUTTON, self.set_background_density, self.b_apply_background_den)

        # DEFINE SIZER
        sizer = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        sizer.AddMany([self.background_den_label_upper, self.background_density_text_input_upper,
                       # self.background_den_label_lower, self.background_density_text_input_lower,
                       # self.background_den_label_lid, self.background_density_text_input_lid,
                       self.b_apply_background_den])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_background_density(self, event):
        self.background_density_upper = float(self.background_density_text_input_upper.GetValue())
        # self.background_density_lower = float(self.background_density_text_input_lower.GetValue())
        # self.background_density_lid   = float(self.background_density_text_input_lid.GetValue())

        # CLOSE
        self.EndModal(1)


class BulkShiftDialog(wx.Dialog):
    """APPLY A BULK SHIFT OF A LAYER IN THE Z DIRECTION"""

    def __init__(self, parent, id, title, plotx_list, ploty_list, i):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)
        self.ploty_list = ploty_list
        self.plotx_list = plotx_list
        self.currently_active_layer_id = i
        self.bulk_panel_x = wx.StaticText(input_panel, -1, "Set x bulk shift value:")
        self.bulk_panel_x_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.bulk_panel_x_Text.SetInsertionPoint(0)
        self.bulk_panel_y = wx.StaticText(input_panel, -1, "Set y bulk shift value:")
        self.bulk_panel_y_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.bulk_panel_y_Text.SetInsertionPoint(0)
        sizer = wx.FlexGridSizer(cols=2, hgap=6, vgap=6)
        self.b_bulk_shift_button = wx.Button(input_panel, -1, "Bulk shift")
        self.Bind(wx.EVT_BUTTON, self.bulk_shift_button, self.b_bulk_shift_button)
        sizer.AddMany([self.bulk_panel_x, self.bulk_panel_y, self.bulk_panel_x_Text, self.bulk_panel_y_Text,
                       self.b_bulk_shift_button])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def bulk_shift_button(self, event):
        self.x_shift_value = float(self.bulk_panel_x_Text.GetValue())
        self.y_shift_value = float(self.bulk_panel_y_Text.GetValue())
        current_y = self.ploty_list[self.currently_active_layer_id]
        current_x = self.plotx_list[self.currently_active_layer_id]

        # BULKSHIFT Y NODES, SET TO ZERO IF NEW VALUE IS < 0
        self.new_x = [x + self.x_shift_value if x + self.x_shift_value > 0. else 0.01 for x in current_x]
        self.new_y = [y + self.y_shift_value if y + self.y_shift_value > 0. else 0.01 for y in current_y]

        # CLOSE
        self.EndModal(1)


class SetObsRmsDialog(wx.Dialog):
    """SET THE DATASET TO USE FOR CALCULATING THE RMS MISFIT"""

    def __init__(self, parent, id, title, observed_data_list):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        '# CHOOSE INPUT OBSERVED DATA'
        self.obs_combo_list_text = wx.StaticText(input_panel, -1, "Observed data file:")
        self.obs_combo_list = wx.ComboBox(input_panel, id=-1, value="", choices=[])

        for i in range(0, len(observed_data_list)):
            try:
                self.obs_combo_list.Append(str(observed_data_list[i].name))
            except:
                pass

        '# DEFINE SET BUTTON'
        self.b_set = wx.Button(input_panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.set_obs, self.b_set)

        '# DEFINE SIZER'
        sizer = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        sizer.AddMany([self.obs_combo_list_text, self.obs_combo_list,
                       self.b_set])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_obs(self, event):
        self.obs_name = str(self.obs_combo_list.GetValue())
        self.EndModal(1)


class NewLayerDialog(wx.Dialog):
    """CREATE A NEW MODEL LAYER"""

    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        # CREATE FIXED LAYER BUTTON
        self.b_fixed = wx.Button(input_panel, -1, "New fixed layer")
        self.Bind(wx.EVT_BUTTON, self.set_fixed, self.b_fixed)
        # CREATE FLOATING LAYER BUTTON
        self.b_floating = wx.Button(input_panel, -1, "New floating layer")
        self.Bind(wx.EVT_BUTTON, self.set_floating, self.b_floating)
        #  DEFINE SIZER
        sizer = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        sizer.AddMany([self.b_fixed, self.b_floating])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_fixed(self, event):
        """ APPEND NEW LAYER BELOW LATEST FIXED LAYER"""
        self.fixed = True
        floating_dialogbox = self.SetNewThickness(self, -1, "Set new layer thickness")
        answer = floating_dialogbox.ShowModal()
        self.new_thickness = floating_dialogbox.new_thickness
        self.EndModal(1)

    def set_floating(self, event):
        """ APPEND NEW LAYER BELOW LATEST FIXED LAYER"""
        self.fixed = False
        self.EndModal(1)

        # floating_dialogbox = self.SetFloatingLayer(self, -1, "Set new layer nodes")
        # answer = floating_dialogbox.ShowModal()
        # self.x1, self.y1 = floating_dialogbox.x1, floating_dialogbox.y1
        # self.x2, self.y2 = floating_dialogbox.x2, floating_dialogbox.y2
        # self.x3, self.y3 = floating_dialogbox.x3, floating_dialogbox.y3
        # self.x4, self.y4 = floating_dialogbox.x4, floating_dialogbox.y4


    class SetNewThickness(wx.Dialog):
        """APPEND NEW FLOATING LAYER AT USER SPECIFIED POSITION. RETURNS XY NODE POSITIONS"""

        def __init__(self, parent, id, title):
            wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                              | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
            floating_panel = wx.Panel(self, -1)

            self.new_Text = wx.StaticText(floating_panel, -1, "Set New layer thickness (km)")
            self.n1_sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.new_thickness_ctrl = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.b_make_button = wx.Button(floating_panel, -1, "Make")
            self.Bind(wx.EVT_BUTTON, self.make_button, self.b_make_button)
            sizer = wx.FlexGridSizer(cols=1, hgap=6, vgap=6)
            sizer.AddMany([self.new_Text, self.new_thickness_ctrl, self.b_make_button])
            floating_panel.SetSizerAndFit(sizer)
            sizer.Fit(self)

        def make_button(self, event):
            self.new_thickness = float(self.new_thickness_ctrl.GetValue())
            self.EndModal(1)

    class SetFloatingLayer(wx.Dialog):
        """APPEND NEW FLOATING LAYER AT USER SPECIFIED POSITION. RETURNS XY NODE POSITIONS"""

        def __init__(self, parent, id, title):
            wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                              | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
            floating_panel = wx.Panel(self, -1)

            self.n1_sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.X1_Text = wx.StaticText(floating_panel, -1, "X1:")
            self.new_x1 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.Y1_Text = wx.StaticText(floating_panel, -1, "Y1:")
            self.new_y1 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.n1_sizer.AddMany([self.X1_Text, self.new_x1, self.Y1_Text, self.new_y1])
            self.n2_sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.X2_Text = wx.StaticText(floating_panel, -1, "X2:")
            self.new_x2 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.Y2_Text = wx.StaticText(floating_panel, -1, "Y2:")
            self.new_y2 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.n2_sizer.AddMany([self.X2_Text, self.new_x2, self.Y2_Text, self.new_y2])
            self.n3_sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.X3_Text = wx.StaticText(floating_panel, -1, "X3:")
            self.new_x3 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.Y3_Text = wx.StaticText(floating_panel, -1, "Y3:")
            self.new_y3 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.n3_sizer.AddMany([self.X3_Text, self.new_x3, self.Y3_Text, self.new_y3])
            self.n4_sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.X4_Text = wx.StaticText(floating_panel, -1, "X4:")
            self.new_x4 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.Y4_Text = wx.StaticText(floating_panel, -1, "Y4:")
            self.new_y4 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
            self.n4_sizer.AddMany([self.X4_Text, self.new_x4, self.Y4_Text, self.new_y4])
            sizer = wx.FlexGridSizer(cols=1, hgap=6, vgap=6)
            self.b_make_button = wx.Button(floating_panel, -1, "Make")
            self.Bind(wx.EVT_BUTTON, self.make_button, self.b_make_button)
            sizer.AddMany([self.n1_sizer, self.n2_sizer, self.n3_sizer, self.n4_sizer, self.b_make_button])
            floating_panel.SetSizerAndFit(sizer)
            sizer.Fit(self)

        def make_button(self, event):
            self.x1, self.y1 = float(self.new_x1.GetValue()), float(self.new_y1.GetValue())
            self.x2, self.y2 = float(self.new_x2.GetValue()), float(self.new_y2.GetValue())
            self.x3, self.y3 = float(self.new_x3.GetValue()), float(self.new_y3.GetValue())
            self.x4, self.y4 = float(self.new_x4.GetValue()), float(self.new_y4.GetValue())
            self.EndModal(1)


class NewFaultDialog(wx.Dialog):
    """CREATE A NEW FAULT. RETURNS XY NODE POSITIONS"""

    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)

        floating_panel = wx.Panel(self, -1)

        self.n1_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.X1_Text = wx.StaticText(floating_panel, -1, "X1:")
        self.new_x1 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
        self.Y1_Text = wx.StaticText(floating_panel, -1, "Y1:")
        self.new_y1 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
        self.n1_sizer.AddMany([self.X1_Text, self.new_x1, self.Y1_Text, self.new_y1])
        self.n2_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.X2_Text = wx.StaticText(floating_panel, -1, "X2:")
        self.new_x2 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
        self.Y2_Text = wx.StaticText(floating_panel, -1, "Y2:")
        self.new_y2 = wx.TextCtrl(floating_panel, -1, "0", size=(100, -1))
        self.n2_sizer.AddMany([self.X2_Text, self.new_x2, self.Y2_Text, self.new_y2])
        sizer = wx.FlexGridSizer(cols=1, hgap=6, vgap=6)
        self.b_make_button = wx.Button(floating_panel, -1, "Make")
        self.Bind(wx.EVT_BUTTON, self.make_button, self.b_make_button)
        sizer.AddMany([self.n1_sizer, self.n2_sizer, self.b_make_button])
        floating_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def make_button(self, event):
        self.x1, self.y1 = float(self.new_x1.GetValue()), float(self.new_y1.GetValue())
        self.x2, self.y2 = float(self.new_x2.GetValue()), float(self.new_y2.GetValue())
        self.EndModal(1)


class MessageDialog(wx.MessageDialog):
    """GENERIC MESSAGE DIALOG BOX"""

    def __init__(self, parent, id, message_text, title):
        wx.MessageDialog.__init__(self, parent, message_text, title)
        dlg = wx.MessageDialog(self, message_text, title, wx.OK)
        answer = dlg.ShowModal()
        dlg.Destroy()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# POP OUT FRAMES ARE PLACED HERE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class CaptureCoordinates(wx.Frame):
    """CAPTURE MOUSE CLICK COORDINATES AND WRITE TO DISK FILE. RETURNS ASCII TEXT FILE"""
    def __init__(coordinate_list, parent, id, title):
        wx.Frame.__init__(coordinate_list, None, wx.ID_ANY, 'Capture coordinates', size=(350, 500))
        coordinate_list.input_panel = wx.Panel(coordinate_list)

        # SET INSTANCE OF GMG CLASS TO RECEIVE NEW ATTRIBUTES
        coordinate_list.parent = parent

        # BIND PROGRAM EXIT BUTTON WITH EXIT FUNCTION
        coordinate_list.Bind(wx.EVT_CLOSE, coordinate_list.on_close_button)

        # CREATE LIST CONTROL
        coordinate_list.table = wx.ListCtrl(coordinate_list.input_panel, size=(350, 500), style=wx.LC_REPORT)
        coordinate_list.table.InsertColumn(0, 'X')
        coordinate_list.table.InsertColumn(1, 'Y')

        # CREATE SAVE BUTTON
        coordinate_list.save_btn = wx.Button(coordinate_list.input_panel, label="Save coordinates")
        coordinate_list.save_btn.Bind(wx.EVT_BUTTON, coordinate_list.on_save)

        # ADD FEATURES TO SIZER
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(coordinate_list.table, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(coordinate_list.save_btn, 0, wx.ALL | wx.CENTER, 5)
        coordinate_list.input_panel.SetSizer(sizer)
        sizer.Fit(coordinate_list)

    def on_save(coordinate_list, event):
        # CREATE OUTPUT FILE
        save_file_dialog = wx.FileDialog(coordinate_list, "Save XY data", "", "", "xy files (*.xy)|*.xy",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # USER CHANGED THEIR MIND

        # NOW WRITE OUT THE DATA
        # GET OUTPUT FILENAME
        output_stream = save_file_dialog.GetPath()

        # GET DATA
        x = []
        y = []
        for i in range(coordinate_list.table.GetItemCount()):
            x.append(float(coordinate_list.table.GetItem(itemIdx=i, col=0).GetText()))
            y.append(float(coordinate_list.table.GetItem(itemIdx=i, col=1).GetText()))

        # OUTPUT DATA
        np.savetxt(output_stream, zip(x, y), delimiter=' ', fmt='%0.6f %0.6f')

    def on_close_button(coordinate_list, event):
        coordinate_list.parent.capture = False
        coordinate_list.Destroy()


class PlotSettingsDialog(wx.Frame):
    """CREATE AN EXTERNAL FIGURE PLOT FROM THE MODEL. RETURNS DRAW PARAMETERS"""

    def __init__(self, parent, id, title, model_aspect, dcanvas_aspect):
        wx.Frame.__init__(self, None, wx.ID_ANY, 'Figure settings')
        input_panel = wx.Panel(self, -1)

        # CREATE INSTANCE OF MAIN FRAME CLASS TO RECEIVE NEW ATTRIBUTES
        self.parent = parent

        self.model_aspect = model_aspect
        self.grav_frame_aspect = dcanvas_aspect

        # CREATE MAIN WINDOW FRAME
        self.main_box = wx.BoxSizer(wx.HORIZONTAL)

        # FIGURE FILE
        self.file_path_text = wx.TextCtrl(input_panel, -1, value="model_figure", size=(150, -1))
        self.file_types = ['pdf', 'png', 'eps', 'ps']
        self.file_type_text = wx.ComboBox(input_panel, -1, value='pdf', choices=self.file_types, size=(75, -1),
                                          style=wx.CB_DROPDOWN)
        self.b_file_path = wx.Button(input_panel, -1, "File...")
        self.Bind(wx.EVT_BUTTON, self.file_path, self.b_file_path)

        # FIGURE FONT SIZE
        self.set_fs = wx.StaticText(input_panel, -1, "Text size:")
        self.fs_text = fs.FloatSpin(input_panel, -1, min_val=0.01, max_val=20.0, increment=0.1, value=8.,
                                    size=(75, -1))
        self.fs_text.SetFormat("%f")
        self.fs_text.SetDigits(3)
        self.place_holder_text_02 = wx.StaticText(input_panel, -1, "")

        # MODEL ASPECT RATIO
        self.set_aspect_ratio = wx.StaticText(input_panel, -1, "Model Aspect ratio:")
        self.aspect_ratio_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=2000.0, increment=0.1,
                                              value=self.model_aspect, size=(75, -1))

        # POLYGON ALPHA VALUE
        self.set_poly_alpha = wx.StaticText(input_panel, -1, " Polygon \n alpha val:")
        self.poly_alpha_text = fs.FloatSpin(input_panel, -1, min_val=0.01, max_val=1.0, increment=0.01, value=0.5,
                                            size=(75, -1))
        self.poly_alpha_text.SetFormat("%f")
        self.poly_alpha_text.SetDigits(3)
        self.place_holder_text_03 = wx.StaticText(input_panel, -1, "")

        # MARKER SIZE
        self.set_ms = wx.StaticText(input_panel, -1, "Observed point size:")
        self.ms_text = fs.FloatSpin(input_panel, -1, min_val=0.01, max_val=20.0, increment=0.1, value=0.5,
                                    size=(75, -1))
        self.ms_text.SetFormat("%f")
        self.ms_text.SetDigits(3)
        self.place_holder_text_04 = wx.StaticText(input_panel, -1, "")

        # LINE WIDTH
        self.set_lw = wx.StaticText(input_panel, -1, " Calculated anomaly \n line width:")
        self.lw_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=1.0,
                                    size=(75, -1))
        self.lw_text.SetFormat("%f")
        self.lw_text.SetDigits(3)
        self.place_holder_text_05 = wx.StaticText(input_panel, -1, "")

        # FONT TYPE
        self.fonts = ['Times New Roman', 'Times', 'Courier', 'Courier New', 'Helvetica', 'Sans', 'verdana', 'Arial']
        self.set_font_type = wx.StaticText(input_panel, -1, "Font type:")
        self.font_type_text = wx.ComboBox(input_panel, -1, value='Times New Roman', choices=self.fonts, size=(75, -1),
                                          style=wx.CB_DROPDOWN)
        self.place_holder_text_06 = wx.StaticText(input_panel, -1, "")

        # LAYER POLYGONS
        self.draw_polygons_checkbox = wx.CheckBox(input_panel, -1, " Draw layer \n polygons")
        self.draw_polygons_checkbox.SetValue(True)

        # LAYER LINES
        self.draw_layer_lines_checkbox = wx.CheckBox(input_panel, -1, " Draw layer \n lines")
        self.draw_layer_lines_checkbox.SetValue(True)
        self.place_holder_text_07 = wx.StaticText(input_panel, -1, "")

        # FLOATING LAYER LINES
        self.draw_floating_layer_lines_checkbox = wx.CheckBox(input_panel, -1, " Draw floating \n layer outlines")
        self.draw_floating_layer_lines_checkbox.SetValue(True)

        # DRAW COLORBAR LINES
        self.draw_colorbar_checkbox = wx.CheckBox(input_panel, -1, "Draw colorbar")
        self.draw_colorbar_checkbox.SetValue(True)
        self.place_holder_text_08 = wx.StaticText(input_panel, -1, "")

        # DRAW XY DATA
        self.draw_xy_checkbox = wx.CheckBox(input_panel, -1, "Draw XY data")
        self.draw_xy_checkbox.SetValue(True)

        # DRAW WELLS
        self.draw_wells_checkbox = wx.CheckBox(input_panel, -1, "Draw Wells")
        self.draw_wells_checkbox.SetValue(True)
        self.place_holder_text_09 = wx.StaticText(input_panel, -1, "")

        # DRAW FAULTS
        self.draw_faults_checkbox = wx.CheckBox(input_panel, -1, "Draw Faults")
        self.draw_faults_checkbox.SetValue(True)
        self.place_holder_text_10 = wx.StaticText(input_panel, -1, "")

        # TIGHT LAYOUT
        self.use_tight_layout_checkbox = wx.CheckBox(input_panel, -1, " Tight \n layout?")
        self.draw_wells_checkbox.SetValue(True)
        self.place_holder_text_11 = wx.StaticText(input_panel, -1, "")

        # WELL FONT SIZE
        self.set_well_font_size = wx.StaticText(input_panel, -1, "Well font size:")
        self.well_font_size_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=1.5,
                                                size=(75, -1))
        self.well_font_size_text.SetFormat("%f")
        self.well_font_size_text.SetDigits(2)
        self.place_holder_text_12 = wx.StaticText(input_panel, -1, "")

        # WELL LINE WIDTH
        self.set_well_lw = wx.StaticText(input_panel, -1, "Well line width:")
        self.well_lw_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=0.1,
                                         size=(75, -1))
        self.well_lw_text.SetFormat("%f")
        self.well_lw_text.SetDigits(3)
        self.place_holder_text_13 = wx.StaticText(input_panel, -1, "")

        # XY POINT SIZE
        self.set_xy_size = wx.StaticText(input_panel, -1, "XY point size:")
        self.xy_size_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=1.5,
                                         size=(75, -1))
        self.xy_size_text.SetFormat("%f")
        self.xy_size_text.SetDigits(2)
        self.place_holder_text_14 = wx.StaticText(input_panel, -1, "")

        # XY COLOR
        self.colors = ['red', 'orange', 'yellow', 'green', 'blue', 'grey', 'white', 'black']
        self.set_xy_color = wx.StaticText(input_panel, -1, "XY color:")
        self.xy_color_text = wx.ComboBox(input_panel, -1, value='black', choices=self.colors, size=(75, -1),
                                         style=wx.CB_DROPDOWN)
        self.place_holder_text_15 = wx.StaticText(input_panel, -1, "")

        # LAYER LINE WIDTH
        self.set_layer_lw = wx.StaticText(input_panel, -1, "Layer line width:")
        self.layer_lw_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=0.1,
                                          size=(75, -1))
        self.layer_lw_text.SetFormat("%f")
        self.layer_lw_text.SetDigits(3)
        self.place_holder_text_16 = wx.StaticText(input_panel, -1, "")

        # LAYER TRANSPARENCY
        self.set_layer_alpha = wx.StaticText(input_panel, -1, "Layer transparency:")
        self.layer_alpha_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.7,
                                             size=(75, -1))
        self.layer_alpha_text.SetFormat("%f")
        self.layer_alpha_text.SetDigits(2)
        self.place_holder_text_17 = wx.StaticText(input_panel, -1, "")

        # COLORBAR XY
        self.colorbar_xy = wx.StaticText(input_panel, -1, "Colorbar xy position:")
        self.colorbar_x_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.7,
                                            size=(75, -1))
        self.colorbar_x_text.SetFormat("%f")
        self.colorbar_x_text.SetDigits(2)
        self.colorbar_y_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.55,
                                            size=(75, -1))
        self.colorbar_y_text.SetFormat("%f")
        self.colorbar_y_text.SetDigits(2)

        # COLORBAR SIZE
        self.colorbar_size = wx.StaticText(input_panel, -1, "Colorbar size:")
        self.colorbar_size_x_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.15,
                                                 size=(75, -1))
        self.colorbar_size_x_text.SetFormat("%f")
        self.colorbar_size_x_text.SetDigits(2)
        self.colorbar_size_y_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.005,
                                                 size=(75, -1))
        self.colorbar_size_y_text.SetFormat("%f")
        self.colorbar_size_y_text.SetDigits(3)

        # GRAV FRAME Y VALUES
        self.grav_frame_text = wx.StaticText(input_panel, -1, "Gravity frame \nY-axis max/min:")
        self.grav_min_text = fs.FloatSpin(input_panel, -1, min_val=-2000., max_val=2000.0, increment=1, value=-100,
                                          size=(75, -1))
        self.grav_min_text.SetFormat("%f")
        self.grav_min_text.SetDigits(2)
        self.grav_max_text = fs.FloatSpin(input_panel, -1, min_val=-2000., max_val=2000.0, increment=1, value=100,
                                          size=(75, -1))
        self.grav_max_text.SetFormat("%f")
        self.grav_max_text.SetDigits(2)

        # DRAW BUTTON
        self.b_draw_button = wx.Button(input_panel, -1, "Draw figure")
        self.Bind(wx.EVT_BUTTON, self.draw_button, self.b_draw_button)

        # EXIT BUTTON
        self.b_exit = wx.Button(input_panel, -1, "Exit")
        self.Bind(wx.EVT_BUTTON, self.exit, self.b_exit)

        # MAKE SIZER
        sizer = wx.FlexGridSizer(cols=3, hgap=8, vgap=8)
        sizer.AddMany([self.file_path_text, self.file_type_text, self.b_file_path,
                       self.set_fs, self.fs_text, self.place_holder_text_02,
                       self.set_aspect_ratio, self.aspect_ratio_text, self.place_holder_text_03,
                       self.set_poly_alpha, self.poly_alpha_text, self.place_holder_text_04,
                       self.set_ms, self.ms_text, self.place_holder_text_05,
                       self.set_lw, self.lw_text, self.place_holder_text_06,
                       self.set_font_type, self.font_type_text, self.place_holder_text_07,
                       self.draw_polygons_checkbox, self.draw_layer_lines_checkbox, self.place_holder_text_08,
                       self.draw_floating_layer_lines_checkbox, self.draw_colorbar_checkbox, self.place_holder_text_09,
                       self.draw_xy_checkbox, self.draw_wells_checkbox, self.place_holder_text_10,
                       self.draw_faults_checkbox, self.use_tight_layout_checkbox, self.place_holder_text_11,
                       self.set_well_font_size, self.well_font_size_text, self.place_holder_text_12,
                       self.set_well_lw, self.well_lw_text, self.place_holder_text_13,
                       self.set_xy_size, self.xy_size_text, self.place_holder_text_14,
                       self.set_xy_color, self.xy_color_text, self.place_holder_text_15,
                       self.set_layer_lw, self.layer_lw_text, self.place_holder_text_16,
                       self.set_layer_alpha, self.layer_alpha_text, self.place_holder_text_17,
                       self.colorbar_xy, self.colorbar_x_text, self.colorbar_y_text,
                       self.colorbar_size, self.colorbar_size_x_text, self.colorbar_size_y_text,
                       self.grav_frame_text, self.grav_min_text, self.grav_max_text,
                       self.b_draw_button, self.b_exit])

        self.main_box.Add(sizer, proportion=1, flag=wx.ALL | wx.EXPAND, border=10)
        input_panel.SetSizerAndFit(self.main_box)
        self.main_box.Fit(self)

    def draw_button(self, event):
        self.file_path = str(self.file_path_text.GetValue())
        self.file_type = str(self.file_type_text.GetValue())
        self.fs = float(self.fs_text.GetValue())
        self.ms = float(self.ms_text.GetValue())
        self.lw = float(self.lw_text.GetValue())
        self.font_type = str(self.font_type_text.GetValue())
        self.poly_alpha = float(self.poly_alpha_text.GetValue())
        self.aspect_ratio = float(self.aspect_ratio_text.GetValue())
        self.use_tight_layout = self.use_tight_layout_checkbox.GetValue()
        self.draw_polygons = self.draw_polygons_checkbox.GetValue()
        self.draw_layers = self.draw_layer_lines_checkbox.GetValue()
        self.draw_floating_layers = self.draw_floating_layer_lines_checkbox.GetValue()
        self.draw_colorbar = self.draw_colorbar_checkbox.GetValue()
        self.draw_wells = self.draw_wells_checkbox.GetValue()
        self.well_fs = float(self.well_font_size_text.GetValue())
        self.well_line_width = float(self.well_lw_text.GetValue())
        self.draw_faults = self.draw_faults_checkbox.GetValue()
        self.draw_xy_data = self.draw_xy_checkbox.GetValue()
        self.xy_size = self.xy_size_text.GetValue()
        self.xy_color = str(self.xy_color_text.GetValue())
        self.colorbar_x = float(self.colorbar_x_text.GetValue())
        self.colorbar_y = float(self.colorbar_y_text.GetValue())
        self.colorbar_size_x = float(self.colorbar_size_x_text.GetValue())
        self.colorbar_size_y = float(self.colorbar_size_y_text.GetValue())
        self.layer_line_width = float(self.layer_lw_text.GetValue())
        self.layer_alpha = float(self.layer_alpha_text.GetValue())
        self.grav_frame_min = float(self.grav_min_text.GetValue())
        self.grav_frame_max = float(self.grav_max_text.GetValue())

        self.parent.draw_model()

    def exit(self, event):
        self.Destroy()

    def file_path(self, event):
        self.save_file_dialog = wx.FileDialog(self, "Save As", "", "", "Figure files (*.*)|*.*",
                                              wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if self.save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # the user changed idea...
        self.chosen_path = self.save_file_dialog.GetPath()
        self.file_path_text.SetValue(str(self.chosen_path))
        self.save_file_dialog.Destroy()


class AttributeEditor(wx.Frame):
    """OPENS A TABLE FOR VIEWING AND EDITING LABEL ATTRIBUTES"""

    def __init__(attribute_edit, parent, id, title, tree_items, densities, reference_densities, susceptibilities,
                 angle_a, angle_b, layer_colors):
        wx.Frame.__init__(attribute_edit, None, wx.ID_ANY, 'Attribute editor', size=(650, 1000))
        attribute_edit.input_panel = wx.Panel(attribute_edit)

        # CREATE INSTANCE OF MAIN FRAME CLASS TO RECEIVE NEW ATTRIBUTES
        attribute_edit.parent = parent

        # SET VARIABLES
        attribute_edit.tree_items = tree_items
        attribute_edit.densities = densities
        attribute_edit.reference_densities = reference_densities
        attribute_edit.susceptibilities = susceptibilities
        attribute_edit.angle_a = angle_a
        attribute_edit.angle_b = angle_b
        attribute_edit.layer_colors = layer_colors

        # DEFINE ATTRIBUTE GRID
        attribute_edit.attr_grid = gridlib.Grid(attribute_edit.input_panel, -1, size=(650, 980))
        attribute_edit.attr_grid.CreateGrid(len(attribute_edit.tree_items) - 1, 7)
        attribute_edit.attr_grid.SetColLabelValue(0, 'Layer Name')
        attribute_edit.attr_grid.SetColLabelValue(1, 'Density')
        attribute_edit.attr_grid.SetColLabelValue(2, 'Density Reference')
        attribute_edit.attr_grid.SetColLabelValue(3, 'Susceptibility')
        attribute_edit.attr_grid.SetColLabelValue(4, 'Angle A')
        attribute_edit.attr_grid.SetColLabelValue(5, 'Angle B')
        attribute_edit.attr_grid.SetColLabelValue(6, 'Layer color')

        # SET COLUMN FORMATS
        attribute_edit.attr_grid.SetColFormatFloat(1, 3, 2)
        attribute_edit.attr_grid.SetColFormatFloat(2, 3, 2)
        attribute_edit.attr_grid.SetColFormatFloat(3, 9, 7)
        attribute_edit.attr_grid.SetColFormatFloat(4, 3, 2)
        attribute_edit.attr_grid.SetColFormatFloat(5, 3, 2)

        # CREATE SET BUTTON
        attribute_edit.b_set_attr_button = wx.Button(attribute_edit.input_panel, -1, "Set")
        attribute_edit.Bind(wx.EVT_BUTTON, attribute_edit.set_attr_button, attribute_edit.b_set_attr_button)

        # POPULATE ATTRIBUTE TABLE
        attribute_edit.length = len(attribute_edit.tree_items)

        for i in range(attribute_edit.length - 1):
            attribute_edit.attr_grid.SetCellValue(i, 0, attribute_edit.tree_items[i + 1])
            attribute_edit.attr_grid.SetCellValue(i, 1, str(float(attribute_edit.densities[i + 1]) / 1000.))
            attribute_edit.attr_grid.SetCellValue(i, 2, str(float(attribute_edit.reference_densities[i + 1]) / 1000.))
            attribute_edit.attr_grid.SetCellValue(i, 3, str(attribute_edit.susceptibilities[i]))
            attribute_edit.attr_grid.SetCellValue(i, 4, str(attribute_edit.angle_a[i + 1]))
            attribute_edit.attr_grid.SetCellValue(i, 5, str(attribute_edit.angle_b[i + 1]))
            attribute_edit.attr_grid.SetCellValue(i, 6, str(attribute_edit.layer_colors[i + 1]))

        # SET SIZER
        for col in range(6):
            attribute_edit.attr_grid.SetColSize(col, 12)

        attribute_edit.attr_grid.AutoSize()
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(attribute_edit.b_set_attr_button)
        sizer.Add(attribute_edit.attr_grid)
        attribute_edit.input_panel.SetSizer(sizer)
        sizer.Fit(attribute_edit)

        # ACTION BINDINGS
        attribute_edit.Bind(wx.EVT_CHAR_HOOK, attribute_edit.on_key)
        attribute_edit.attr_grid.Bind(wx.grid.EVT_GRID_CELL_LEFT_DCLICK, attribute_edit.open_colour_box)
        attribute_edit.attr_grid.Bind(wx.EVT_SIZE, attribute_edit.on_size)

    def on_size(attribute_edit, event):
        width, height = attribute_edit.GetClientSize()
        for col in range(6):
            attribute_edit.attr_grid.SetColSize(col, width / (10 + 1))

    def open_colour_box(attribute_edit, event):
        if event.GetCol() == 5:
            row = event.GetRow()
            attribute_edit.on_color_dlg(event, row)
        else:
            pass

    def on_color_dlg(attribute_edit, event, row):
        """SET COLOUR FOR LAYER"""
        dlg = wx.ColourDialog(attribute_edit)

        # ENSURE THE FULL COLOUR DIALOG IS DISPLAYED, NOT THE ABBREVIATED VERSION
        dlg.GetColourData().SetChooseFull(True)

        if dlg.ShowModal() == wx.ID_OK:
            rgb = dlg.GetColourData().GetColour().Get()
            html = struct.pack('BBB', *rgb).encode('hex')
            attribute_edit.attr_grid.SetCellValue(row, 5, '#' + str(html))
        dlg.Destroy()

    def on_key(attribute_edit, event):
        if event.ControlDown() and event.GetKeyCode() == 67:
            attribute_edit.selection()
            attribute_edit.copy()  # CALL COPY METHOD
        # SKIP OTHER KEY EVENTS
        if event.GetKeyCode():
            event.Skip()
            return

    def selection(attribute_edit):
        # SHOW CELL SELECTION
        # IF SELECTION IS CELL
        if attribute_edit.attr_grid.GetSelectedCells():
            print("Selected cells " + str(attribute_edit.GetSelectedCells()))
        # IF SELECTION IS BLOCK
        if attribute_edit.attr_grid.GetSelectionBlockTopLeft():
            print("Selection block top left " + str(attribute_edit.attr_grid.GetSelectionBlockTopLeft()))
        if attribute_edit.attr_grid.GetSelectionBlockbottomRight():
            print("Selection block bottom right " + str(attribute_edit.attr_grid.GetSelectionBlockbottomRight()))
        # IF SELECTION IS COL
        if attribute_edit.attr_grid.GetSelectedCols():
            print("Selected cols " + str(attribute_edit.attr_grid.GetSelectedCols()))
        # IF SELECTION IS ROW
        if attribute_edit.attr_grid.GetSelectedRows():
            print("Selected rows " + str(attribute_edit.attr_grid.GetSelectedRows()))

    def currentcell(attribute_edit):
        # SHOW CURSOR POSITION
        row = attribute_edit.attr_grid.GetGridCursorRow()
        col = attribute_edit.attr_grid.GetGridCursorCol()
        cell = (row, col)
        print("Current cell " + str(cell))

    def copy(attribute_edit):
        # NUMBER OF ROWS AND COLS
        rows = attribute_edit.attr_grid.GetSelectionBlockbottomRight()[0][0] - \
               attribute_edit.attr_grid.GetSelectionBlockTopLeft()[0][0] + 1
        cols = attribute_edit.attr_grid.GetSelectionBlockbottomRight()[0][1] - \
               attribute_edit.attr_grid.GetSelectionBlockTopLeft()[0][1] + 1

        # DATA VARIABLE CONTAIN TEXT THAT MUST BE SET IN THE CLIPBOARD
        data = ''

        # FOR EACH CELL IN SELECTED range APPEND THE CELL VALUE IN THE DATA VARIABLE
        # TABS '\t' FOR COLS AND '\r' FOR ROWS
        for r in range(rows):
            for c in range(cols):
                data = data + str(
                    attribute_edit.attr_grid.GetCellValue(attribute_edit.attr_grid.GetSelectionBlockTopLeft()[0][0] + r,
                                                          attribute_edit.attr_grid.GetSelectionBlockTopLeft()[0][
                                                              1] + c))
                if c < cols - 1:
                    data = data + '\t'
            data = data + '\n'

        # CREATE TEXT DATA OBJECT
        clipboard = wx.TextDataObject()

        # SET DATA OBJECT VALUE
        clipboard.SetText(data)

        # PUT THE DATA IN THE CLIPBOARD
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Error")

    def set_attr_button(attribute_edit, event):
        # RECREATE ARRAYS (INCLUDE VALUES FOR "LAYER 0)"
        attribute_edit.tree_items = ['Layer 1']
        attribute_edit.densities = [0.0]
        attribute_edit.reference_densities = [0.0]
        attribute_edit.susceptibilities = [0.0]
        attribute_edit.angle_a = [0.0]
        attribute_edit.angle_b = [0.0]
        attribute_edit.layer_colors = ['b']

        # SET NEW DATA FOR LAYERS 1 to i
        for i in range(attribute_edit.length - 1):
            attribute_edit.tree_items.append(str(attribute_edit.attr_grid.GetCellValue(i, 0)))
            attribute_edit.densities.append(float(attribute_edit.attr_grid.GetCellValue(i, 1)) * 1000.)
            attribute_edit.reference_densities.append(float(attribute_edit.attr_grid.GetCellValue(i, 2)) * 1000.)
            attribute_edit.susceptibilities.append(float(attribute_edit.attr_grid.GetCellValue(i, 3)))
            attribute_edit.angle_a.append(float(attribute_edit.attr_grid.GetCellValue(i, 4)))
            attribute_edit.angle_b.append(float(attribute_edit.attr_grid.GetCellValue(i, 5)))
            attribute_edit.layer_colors.append(str(attribute_edit.attr_grid.GetCellValue(i, 6)))

        # UPDATE MAIN FRAME
        attribute_edit.parent.attribute_set(attribute_edit.tree_items, attribute_edit.densities,
                                            attribute_edit.reference_densities, attribute_edit.susceptibilities,
                                            attribute_edit.angle_a, attribute_edit.angle_b, attribute_edit.layer_colors)

        attribute_edit.parent.update_layer_data()
        attribute_edit.parent.run_algorithms()
        attribute_edit.parent.draw()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# START SOFTWARE
if __name__ == "__main__":
    app = wx.App(False)
    fr = wx.Frame(None, title='GMG: Geophysical Modelling GUI')
    app.frame = Gmg()
    app.frame.CenterOnScreen()
    app.frame.Show()
    app.MainLoop()
