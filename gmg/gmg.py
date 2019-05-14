"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI application for Forward modelling 2D potential field profiles.
Written by Brook Tozer, University of Oxford 2015-17.
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
Gravity algorithm written using NumPy by Brook Tozer (2015) (Modified from the Fatiando a Terra grav code).

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

NB. before launching gmg, run:

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
from fatiando.mesher import Polygon
import plot_model
import bott
import talwani_and_heirtzler
import model_stats
import struct
import gc
import webbrowser


# FUTURE
# import wx.lib.agw.ribbon as RB
# import wx.EnhancedStatusBar as ESB
# from scipy import interpolate as ip

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        self.dir = os.path.dirname(os.path.abspath(__file__)).split("/")
        self.gui_icons_dir = "/"+self.dir[1]+"/"+self.dir[2]+"/"+self.dir[3]+"/docs/icons/"

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
        sys.Gmg = self
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
        self.t_canvas = True
        self.d_canvas = True
        self.nt_canvas = True
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
        grav_m_observed_filter = self.gravity_data.Append(-1, "Filter Anomaly\tCtrl-L", "Filter Observed Anomaly")
        self.Bind(wx.EVT_MENU, self.observed_filter, grav_m_observed_filter)
        # SET RMS OBS ARRAYS
        grav_m_set_rms_arrays = self.gravity_data.Append(-1, "Set RMS\tCtrl-L", "Set RMS")
        self.Bind(wx.EVT_MENU, self.set_obs_rms, grav_m_set_rms_arrays)
        # SET BACKGROUND DENSITY
        m_set_background_density = self.gravity_data.Append(-1, "&Set Background Density...\tCtrl-shift-down",
                                                            "Set Background Density...")
        self.Bind(wx.EVT_MENU, self.set_background_density, m_set_background_density)
        # REMOVE
        self.m_obs_g_submenu.Append(100, 'Remove All Gravity Anomalies')
        self.Bind(wx.EVT_MENU, self.delete_all_obs_grav, id=100)
        # SAVE PREDICTED ANOMALY TO DISC
        m_save_g_submenu = self.gravity_data.Append(-1, "&Save Predicted Anomaly...\tCtrl-shift-S",
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
        mag_m_observed_filter = self.magnetic_data.Append(-1, "Filter Anomaly\tCtrl-L", "Filter Observed Anomaly")
        self.Bind(wx.EVT_MENU, self.observed_filter, mag_m_observed_filter)
        # SET RMS OBS ARRAYS
        mag_m_set_rms_arrays = self.magnetic_data.Append(-1, "Set RMS\tCtrl-L", "Set RMS")
        self.Bind(wx.EVT_MENU, self.set_obs_rms, mag_m_set_rms_arrays)
        # SET MAG
        m_set_mag_variables = self.magnetic_data.Append(-1, "&Set Magnetic Field...\tCtrl-shift-up",
                                                        "Set Magnetic Feild...")
        self.Bind(wx.EVT_MENU, self.set_mag_variables, m_set_mag_variables)
        # REMOVE
        self.m_obs_mag_submenu.Append(101, 'Remove all Magnetic Anomalies')
        self.Bind(wx.EVT_MENU, self.delete_all_obs_mag, id=101)
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
        # REMOVE
        self.m_xy_submenu.Append(102, 'Remove all XY Data')
        self.Bind(wx.EVT_MENU, self.delete_all_xy, id=102)
        # DRAW MENU
        self.menubar.Append(self.xy_data, "&XY Data")
        # --------------------------------------------------------------------------------------------------------------

        # SEISMIC DATA -------------------------------------------------------------------------------------------------
        self.seismic_data = wx.Menu()
        # SEGY LOAD
        self.m_load_segy = self.seismic_data.Append(-1, "&Load Segy...\tCtrl-y", "Load Segy Data")
        self.Bind(wx.EVT_MENU, self.segy_input, self.m_load_segy)
        # SEGY REMOVE_ALL
        self.m_remove_all_segy = self.seismic_data.Append(-1, "&Remove All Segy...\tCtrl-y", "Remove All Segy Data")
        self.Bind(wx.EVT_MENU, self.remove_all_segy, self.m_remove_all_segy)
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
        self.m_color_palette.Append(902, 'Sesimic')
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
        self.layer_file.Append(-1, 'Transperency', self.m_layer_transperency)
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
        # --------------------------------------------------------------------------------------------------------------

        'SET MENUBAR'
        self.SetMenuBar(self.menubar)

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
        ### self.Bind(wx.EVT_TOOL, self.calc_topo_switch, t_calc_topo)  # FUTURE

        t_calc_model_bott = self.toolbar.AddTool(wx.ID_ANY, "Calculate gravity",
                                                 wx.Bitmap(self.gui_icons_dir + 'G_24.png'),
                                                 shortHelp="Calculate gravity")
        self.Bind(wx.EVT_TOOL, self.calc_grav_switch, t_calc_model_bott)

        t_calc_mag = self.toolbar.AddTool(wx.ID_ANY, "Calculate magnetic",
                                          wx.Bitmap(self.gui_icons_dir + 'M_24.png'), shortHelp="Calculate magnetic")
        self.Bind(wx.EVT_TOOL, self.calc_mag_switch, t_calc_mag)

        t_capture_coordinates = self.toolbar.AddTool(wx.ID_ANY, "Capture coordinates",
                                                     wx.Bitmap(self.gui_icons_dir + 'C_24.png'),
                                                     shortHelp="Capture coordinates")
        self.Bind(wx.EVT_TOOL, self.capture_coordinates, t_capture_coordinates)

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

        t_zoom = self.toolbar.AddTool(wx.ID_ANY, "Zoom in",
                                      wx.Bitmap(self.gui_icons_dir + 'zoom_in_24.png'), shortHelp="Zoom in")
        self.Bind(wx.EVT_TOOL, self.zoom, t_zoom)

        t_zoom_out = self.toolbar.AddTool(wx.ID_ANY, "Zoom out",
                                          wx.Bitmap(self.gui_icons_dir + 'zoom_out_24.png'), shortHelp="Zoom out")
        self.Bind(wx.EVT_TOOL, self.zoom_out, t_zoom_out)

        t_full_extent = self.toolbar.AddTool(wx.ID_ANY, "Full extent",
                                             wx.Bitmap(self.gui_icons_dir + 'full_extent_24.png'),
                                             shortHelp="Full extent")
        self.Bind(wx.EVT_TOOL, self.full_extent, t_full_extent, id=604)

        t_pan = self.toolbar.AddTool(wx.ID_ANY, "Pan",
                                     wx.Bitmap(self.gui_icons_dir + 'pan_24.png'), shortHelp="Pan")
        self.Bind(wx.EVT_TOOL, self.pan, t_pan)

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

        # INCREASE TRANSPERANCY ICON
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
        t_toogle_fault_mode = self.toolbar.AddTool(wx.ID_ANY, "Fault pick",
                                                   wx.Bitmap(self.gui_icons_dir + 'F_24.png'),
                                                   shortHelp="Toogle fault picking")
        self.Bind(wx.EVT_TOOL, self.toogle_fault_mode, t_toogle_fault_mode)

        # FAULT PICKER ICON
        t_fault_pick = self.toolbar.AddTool(wx.ID_ANY, "Fault pick",
                                            wx.Bitmap(self.gui_icons_dir + 'faultline_24.png'),
                                            shortHelp="Fault picker")
        self.Bind(wx.EVT_TOOL, self.pick_new_fault, t_fault_pick)

        # CREATE TOOLBAR
        self.toolbar.Realize()
        self.toolbar.SetSize((1790, 36))

    def start(self, area, xp, zp):
        """CREATE MPL FIGURE CANVAS"""

        self.fig = plt.figure()  # CREATE MPL FIGURE
        self.canvas = FigureCanvas(self.rightPanel, -1, self.fig)  # CREATE FIGURE CANVAS
        self.nav_toolbar = NavigationToolbar(self.canvas)  # CREATE DEFAULT NAVIGATION TOOLBAR
        self.nav_toolbar.Hide()  # HIDE DEFAULT NAVIGATION TOOLBAR

        'SET DRAW COMMAND WHICH CAN BE CALLED TO REDRAW THE FIGURE'
        self.draw = self.fig.canvas.draw

        'GET THE MODEL DIMENSIONS AND SAMPLE LOCATIONS'
        self.area = area
        self.x1, self.x2, self.z1, self.z2 = 0.001 * np.array(area)
        self.xp = np.array(xp, dtype='f')
        self.zp = np.array(zp, dtype='f')

        'DRAW MAIN PROGRAM WINDOW'
        self.draw_main_frame()

        'CONNECT MPL FUNCTIONS'
        self.connect()

        'UPDATE DISPLAY'
        self.display_info()
        self.size_handler()

        if self.newmodel:
            self.update_layer_data()
        else:
            pass

        'REFRESH SIZER POSITIONS'
        self.Hide()
        self.Show()

    def initalise_model(self):
        """INITIALISE OBSERVED DATA AND LAYERS"""

        self.pinch = False
        self.showverts = True
        self.exsiting_node = True
        self.nodes = True
        self.zoom_on = False
        self.pan_on = False
        self.node_click_limit = 0.1
        self.layer_lock = True
        self.boundary_lock = True
        self.calc_padding = 5000.  # PADDING FOR POTENTIAL FIELD CALCULATION POINTS
        self.padding = 50000.  # PADDING FOR LAYERS
        self.i = 0  # LAYER COUNTER
        self.node_layer_reference = 0
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

        'INITIALISE POLYGON LISTS (USED AS MODEL LAYERS)'
        self.mag_polygons = []
        self.polygons = []
        self.polyplots = []
        self.poly_fills = [[]]

        'INITIALISE LAYER LISTS (USED FOR STORING LAYER DATA)'
        self.plotx_list = [[]]
        self.ploty_list = [[]]
        self.layer_colors = [[]]
        self.layer_count = 0
        self.boundary_lock_list = [0]
        self.boundary_lock_status = ['locked']
        self.layer_lock_list = [0]
        self.layer_lock_status = ['locked']
        self.layer_transparency = 0.4
        self.pinch_node_list = [[], []]
        self.pinch_count = 0

        'INITIALISE XY DATA ATTRIBUTES'
        self.xy = []
        self.xy_list = [[]]
        self.xy_list_save = [[]]
        self.xy_name_list = []
        self.xy_color_list =[[]]
        self.xy_count = 0

        'INITIALISE GEOLOGICAL CONTACT ATTRIBUTES'
        self.outcrop_data = [[]]
        self.outcrop_data_list = [[]]
        self.outcrop_data_list_save = [[]]
        self.outcrop_data_name_list = []
        self.outcrop_data_color_list = [[]]
        self.outcrop_text_list = [[]]
        self.outcrop_data_count = 0

        'INITIALISE TOPOGRAPHY ATTRIBUTES'
        self.obs_topo = []
        self.obs_topo_list = [[]]
        self.obs_topo_list_save = [[]]
        self.obs_topo_name_list = []
        self.obs_topo_colors = [[]]
        self.obs_topo_count = 0

        'INITIALISE GRAVITY ATTRIBUTES'
        self.grav_obs_switch = False
        self.obs_grav = []
        self.obs_grav_list = [[]]
        self.obs_grav_list_save = [[]]
        self.obs_grav_name_list = []
        self.obs_grav_colors = [[]]
        self.obs_grav_count = 0
        self.background_density = 0
        self.absolute_densities = True
        self.calc_grav_switch = False
        self.obs_gravity_data_for_rms = []  # OBSERVED DATA LIST TO BE COMPARED TO CALCULATED
        self.grav_rms_value = 0.  # TOTAL RMS MISFIT VALUE
        self.grav_residuals = []  # CALCULATED RESIDUAL

        'INITIALISE MAGNETIC ATTRIBUTES'
        self.mag_obs_switch = False
        self.obs_mag = []
        self.obs_mag_list = [[]]
        self.obs_mag_list_save = [[]]
        self.obs_mag_name_list = []
        self.obs_mag_colors = [[]]
        self.obs_mag_count = 0
        self.earth_field = 0.
        self.profile_azimuth = 0.
        self.calc_mag_switch = False
        self.obs_mag_data_for_rms = []  # OBSERVED DATA LIST TO BE COMPARED TO CALCULATED
        self.mag_rms_value = 0.  # TOTAL RMS MISFIT VALUE (SINGLE INTEGER)
        self.mag_residuals = []  # CALCULATED RESIDUAL

        'INITIALISE SEISMIC ATTRIBUTES'
        self.gain = 4.0
        self.gain_neg = -self.gain
        self.segy_on = False
        self.segy_name_list = []
        self.segy_file_list = []
        self.segy_plot_list = []
        self.segy_dimension_list = []
        self.segy_count = 0

        'INITIALISE Well ATTRIBUTES'
        self.well_list = []
        self.well_list_hidden = []
        self.well_name_list = []
        self.well_name_text = [[]]
        self.well_labels_list = [[]]
        self.well_horizons_list = [[]]
        self.well_count = 0
        self.well_list_switch = [[]]
        self.well_textsize = 2

        'INITIALISE LAYER ATTRIBUTES'
        self.densities = [0.]
        self.reference_densities = [2.67]
        self.susceptibilities = [0.]
        self.angle_a = [0.]
        self.angle_b = [0.]
        self.layers_calculation_switch = [1]

        'INITIALISE FAULTS'
        self.fault_picking_switch = False
        self.faults = [[]]
        self.fault_names_list = []
        self.fault_counter = 0
        self.current_fault_index = 0
        self.fault_colors = ['k']
        self.fault_x_coords_list = [[]]
        self.fault_y_coords_list = [[]]
        self.selected_node = None

        'INITIALIZE COORDINATE CAPTURE'
        self.capture = False
        self.linex = []
        self.liney = []

    def draw_main_frame(self):
        """DRAW THE PROGRAM CANVASES"""
        columns=93
        x_orig=9
        'TOPO CANVAS'
        self.tcanvas = plt.subplot2grid((26, 100), (0, x_orig), rowspan=2, colspan=columns)
        self.tcanvas.set_ylabel("Topo (m)")
        self.tcanvas.xaxis.set_major_formatter(plt.NullFormatter())
        self.tcanvas.grid()
        'GRAV CANVAS'
        self.dcanvas = plt.subplot2grid((26, 100), (2, x_orig), rowspan=3, colspan=columns)
        self.dcanvas.set_ylabel("Grav (mGal)")
        self.dcanvas.xaxis.set_major_formatter(plt.NullFormatter())
        self.dcanvas.grid()
        self.dcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        'MAG CANVAS'
        self.ntcanvas = plt.subplot2grid((26, 100), (5, x_orig), rowspan=3, colspan=columns)
        self.ntcanvas.set_ylabel("Mag (nT)")
        self.ntcanvas.xaxis.set_major_formatter(plt.NullFormatter())
        self.ntcanvas.grid()
        self.ntcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        'MODEL CANVAS'
        self.mcanvas = plt.subplot2grid((26, 100), (8, x_orig), rowspan=17, colspan=columns)
        self.mcanvas.set_ylabel("Depth (km)")
        self.mcanvas.set_xlabel("x (km)")

        'DENSITY COLOUR BAR CANVAS'
        colormap = matplotlib.cm.coolwarm
        cnorm = colors.Normalize(vmin=-0.8, vmax=0.8)
        self.colormap = cm.ScalarMappable(norm=cnorm, cmap=colormap)
        self.scale_canvas = plt.subplot2grid((24, 12), (23, 10), rowspan=1, colspan=2)
        self.cb1 = matplotlib.colorbar.ColorbarBase(self.scale_canvas,
                                                    cmap=colormap, norm=cnorm, orientation='horizontal')
        self.cb1.ax.tick_params(labelsize=8)
        self.cb1.set_label('Density contrast ($kg/m^{3}$)', fontsize=10)

        'SET CANVAS LIMITS'
        self.mcanvas.set_xlim(self.x1, self.x2)
        self.mcanvas.set_ylim(self.z2, self.z1)
        self.mcanvas.grid()
        self.tcanvas.set_xlim(self.mcanvas.get_xlim())
        self.dcanvas.set_xlim(self.mcanvas.get_xlim())
        self.ntcanvas.set_xlim(self.mcanvas.get_xlim())
        self.fig.subplots_adjust(top=0.99, left=-0.045, right=0.99, bottom=0.02,
                                 hspace=1.5)

        'ADD FIRST LAYER'
        if self.newmodel:
            self.plotx = [-(float(self.padding)), 0., self.x2, self.x2 + (float(self.padding))]
            self.ploty = [0.001, 0.001, 0.001, 0.001]
            self.plotx_list[0] = self.plotx
            self.ploty_list[0] = self.ploty
            self.nextpoly = zip(self.plotx, self.ploty)
            self.polyline, = self.mcanvas.plot(self.plotx, self.ploty, marker='o',
                                               color='k', linewidth=1.0, alpha=0.5, picker=True)
            self.plotx_polygon = np.array(self.plotx)
            self.ploty_polygon = np.array(self.ploty)
            self.layer_lines[0] = self.mcanvas.plot(self.plotx, self.ploty, color='black', linewidth=1.0, alpha=1.0)
            self.layer_colors[0] = 'black'
            self.polygon_fills[0] = self.mcanvas.fill(self.plotx_polygon, self.ploty_polygon, color='blue',
                                                      alpha=self.layer_transparency, closed=True, linewidth=None,
                                                      ec=None)
            self.current_node = self.mcanvas.scatter(-40000., 0, s=50, color='r', zorder=10)
        else:
            self.nextpoly = []

        '''ADDITIONAL MAIN FRAME WIDGETS - PLACED ON LEFT HAND SIDE OF THE FRAME'''
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
        self.predplot, = self.dcanvas.plot([], [], '-r', linewidth=2, alpha=0.5)
        self.grav_rms_plot, = self.dcanvas.plot([], [], color='purple', linewidth=1.5, alpha=0.5)
        self.prednt_plot, = self.ntcanvas.plot([], [], '-g', linewidth=2, alpha=0.5)
        self.mag_rms_plot, = self.ntcanvas.plot([], [], color='purple', linewidth=1.5, alpha=0.5)

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

        self.current_xlim = self.mcanvas.get_xlim()
        self.current_ylim = self.mcanvas.get_ylim()

        if event.Id == 601:
            if self.t_canvas is True:
                self.tcanvas.set_visible(False)
                self.t_canvas = False
            else:
                self.t_canvas = True
                self.tcanvas.set_visible(True)
        if event.Id == 602:
            if self.d_canvas is True:
                self.dcanvas.set_visible(False)
                self.d_canvas = False
            else:
                self.d_canvas = True
                self.dcanvas.set_visible(True)
        if event.Id == 603:
            if self.nt_canvas is True:
                self.ntcanvas.set_visible(False)
                self.nt_canvas = False
            else:
                self.nt_canvas = True
                self.ntcanvas.set_visible(True)

        'ADJUST FRAME SIZING AND SET PROGRAM WINDOW'
        if self.t_canvas is True and self.d_canvas is True and self.nt_canvas is True:
            'TRUE TRUE TRUE'
            'TOPO CANVAS'
            self.tcanvas = plt.subplot2grid((26, 12), (0, 1), rowspan=2, colspan=12)
            self.tcanvas.set_ylabel("Topo (m)")
            self.tcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.tcanvas.grid()
            self.tcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            'GRAV CANVAS'
            self.dcanvas = plt.subplot2grid((26, 12), (2, 1), rowspan=3, colspan=12)
            self.dcanvas.set_ylabel("Grav An. (mGal)")
            self.dcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.dcanvas.grid()
            self.dcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            'MAG CANVAS'
            self.ntcanvas = plt.subplot2grid((26, 12), (5, 1), rowspan=3, colspan=12)
            self.ntcanvas.set_ylabel("Mag An. (nT)")
            self.ntcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.ntcanvas.grid()
            self.ntcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        elif self.t_canvas is False and self.d_canvas is True and self.nt_canvas is True:
            'FALSE TRUE TRUE'
            'TOPO CANVAS'
            # HIDDEN
            'GRAV CANVAS'
            self.dcanvas = plt.subplot2grid((26, 12), (0, 1), rowspan=4, colspan=12)
            self.dcanvas.set_ylabel("Grav An. (mGal)")
            self.dcanvas.grid()
            self.dcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            'MAG CANVAS'
            self.ntcanvas = plt.subplot2grid((26, 12), (4, 1), rowspan=4, colspan=12)
            self.ntcanvas.set_ylabel("Mag An. (nT)")
            self.ntcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.ntcanvas.grid()
            self.ntcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        elif self.t_canvas is True and self.d_canvas is False and self.nt_canvas is True:
            'TRUE FALSE TRUE'
            'TOPO CANVAS'
            self.tcanvas = plt.subplot2grid((26, 12), (0, 1), rowspan=2, colspan=12)
            self.tcanvas.set_ylabel("Topo (m)")
            self.tcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.tcanvas.grid()
            'GRAV CANVAS'
            # HIDDEN
            'MAG CANVAS'
            self.ntcanvas = plt.subplot2grid((26, 12), (2, 1), rowspan=6, colspan=12)
            self.ntcanvas.set_ylabel("Mag An. (nT)")
            self.ntcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.ntcanvas.grid()
            self.ntcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        elif self.t_canvas is True and self.d_canvas is True and self.nt_canvas is False:
            'TRUE TRUE FALSE'
            'TOPO CANVAS'
            self.tcanvas = plt.subplot2grid((26, 12), (0, 1), rowspan=2, colspan=12)
            self.tcanvas.set_ylabel("Topo (m)")
            self.tcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.tcanvas.grid()
            'GRAV CANVAS'
            self.dcanvas = plt.subplot2grid((26, 12), (2, 1), rowspan=6, colspan=12)
            self.dcanvas.set_ylabel("Grav An. (mGal)")
            self.dcanvas.grid()
            self.dcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            'MAG CANVAS'
            # HIDDEN

        elif self.t_canvas is False and self.d_canvas is False and self.nt_canvas is True:
            'FALSE FALSE TRUE'
            'TOPO CANVAS'
            # HIDDEN
            'GRAV CANVAS'
            # HIDDEN
            'MAG CANVAS'
            self.ntcanvas = plt.subplot2grid((26, 12), (0, 1), rowspan=8, colspan=12)
            self.ntcanvas.set_ylabel("Mag An. (nT)")
            self.ntcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.ntcanvas.grid()
            self.ntcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        elif self.t_canvas is False and self.d_canvas is True and self.nt_canvas is False:
            'FALSE TRUE FALSE'
            'TOPO CANVAS'
            # HIDDEN
            'GRAV CANVAS'
            self.dcanvas = plt.subplot2grid((26, 12), (0, 1), rowspan=8, colspan=12)
            self.dcanvas.set_ylabel("Grav An. (mGal)")
            self.dcanvas.grid()
            self.dcanvas.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            'MAG CANVAS'
            # HIDDEN

        elif self.t_canvas is True and self.d_canvas is False and self.nt_canvas is False:
            'TRUE FALSE FALSE'
            'TOPO CANVAS'
            self.tcanvas = plt.subplot2grid((26, 12), (0, 1), rowspan=8, colspan=12)
            self.tcanvas.set_ylabel("Topo (m)")
            self.tcanvas.xaxis.set_major_formatter(plt.NullFormatter())
            self.tcanvas.grid()
            'GRAV CANVAS'
            self.dcanvas.set_visible(False)
            'MAG CANVAS'
            self.ntcanvas.set_visible(False)

        elif self.t_canvas is False and self.d_canvas is False and self.nt_canvas is False:
            pass
            'FALSE FALSE FALSE'
            'TOPO CANVAS'
            # HIDDEN
            'GRAV CANVAS'
            # HIDDEN
            'MAG CANVAS'
            # HIDDEN
        self.draw()

        'SET CANVAS LIMITS'
        self.mcanvas.set_xlim(self.current_xlim)
        self.mcanvas.set_ylim(self.current_ylim)
        self.mcanvas.grid()
        if self.tcanvas is not None:
            self.tcanvas.set_xlim(self.mcanvas.get_xlim())
        if self.dcanvas is not None:
            self.dcanvas.set_xlim(self.mcanvas.get_xlim())
        if self.ntcanvas is not None:
            self.ntcanvas.set_xlim(self.mcanvas.get_xlim())
        self.fig.subplots_adjust(top=0.99, left=-0.045, right=0.99, bottom=0.02, hspace=1.5)

        'INITALISE CALCULATED P.F. LINES'
        if self.dcanvas is not None:
            self.predplot, = self.dcanvas.plot([], [], '-r', linewidth=2)
            self.grav_rms_plot, = self.dcanvas.plot([], [], color='purple', linewidth=1.5)
        if self.ntcanvas is not None:
            self.prednt_plot, = self.ntcanvas.plot([], [], '-g', linewidth=2)
            self.mag_rms_plot, = self.ntcanvas.plot([], [], color='purple', linewidth=1.5)

        'PLOT OBSERVED TOPO DATA'
        if self.tcanvas is not None:
            pass
            # for x in range(len(self.obs_topo_list_save)):
            #     if self.obs_topo_list_save[x] is not None:
            #         topo = self.obs_topo_list_save[x]
            #         self.obs_topo_list[x] = self.dcanvas.scatter(topo[:, 0], topo[:, 1], marker='o',
            #                                                      color=self.obs_topo_colors[x], s=5, gid=x)
        'PLOT OBSERVED GRAVITY DATA'
        if self.dcanvas is not None:
            for x in range(len(self.obs_grav_list_save)):
                if len(self.obs_grav_list_save[x]) >= 1:
                    grav = self.obs_grav_list_save[x]
                    grav_x = grav[:, 0]
                    grav_y = grav[:, 1]
                    self.obs_grav_list[x] = self.dcanvas.scatter(grav_x, grav_y, marker='o',
                                                                 color=self.obs_grav_colors[x], s=5, gid=x)
        'PLOT OBSERVED MAG DATA'
        if self.ntcanvas is not None:
            for x in range(len(self.obs_mag_list_save)):
                if len(self.obs_mag_list_save[x]) >= 1:
                    mag = self.obs_mag_list_save[x]
                    mag_x = mag[:, 0]
                    mag_y = mag[:, 1]
                    self.obs_mag_list[x] = self.ntcanvas.scatter(mag_x, mag_y, marker='o',
                                                                color=self.obs_mag_colors[x], s=5, gid=x)

        'UPDATE FRAMES'
        self.update_layer_data()
        self.run_algorithms()
        self.set_frame_limits()
        self.mcanvas.grid(True)
        self.draw()

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
        self.mcanvas.set_ylim(new_z2 / 1000., new_z1 / 1000.)
        self.mcanvas.set_xlim(new_x1 / 1000., new_x2 / 1000.)
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

    # %LAYER TREE FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        self.i = self.tree.GetPyData(event.GetItem())
        self.nextdens = self.densities[self.i]
        self.density_input.SetValue(0.001 * self.densities[self.i])
        self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
        self.susceptibility_input.SetValue(self.susceptibilities[self.i])
        self.angle_a_input.SetValue(self.angle_a[self.i])
        self.angle_b_input.SetValue(self.angle_b[self.i])
        self.plotx = self.plotx_list[self.i]
        self.ploty = self.ploty_list[self.i]
        self.current_node.set_offsets([self.plotx[0], self.ploty[0]])
        self.update_layer_data()
        self.run_algorithms()

    def on_begin_edit_label(self, event):
        # self.DoGetBestSize()
        pass

    def on_end_edit_label(self, event):
        self.i = self.tree.GetPyData(event.GetItem())
        new_label = self.get_item_text(event.GetItem())
        self.tree_items[self.i] = str(new_label)

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
                                     " || MAG RMS = %s  ||" % (self.i, self.mcanvas.get_aspect(), self.grav_rms_value,
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
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def zoom_out(self, event):
        self.nav_toolbar.back()
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def full_extent(self, event):
        """REDRAW MODEL FRAME WITH FULL EXTENT"""
        self.full_extent_adjustment()
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def full_extent_adjustment(self):
        """FIND WHICH FRAME IS REFERENCED & CHANGE SWITCH"""
        if not self.tcanvas.get_visible():
            self.tcanvas.set_visible(True)
        if not self.dcanvas.get_visible():
            self.dcanvas.set_visible(True)
        if not self.ntcanvas.get_visible():
            self.ntcanvas.set_visible(True)

        'SET CANVAS LIMITS'
        self.mcanvas.set_xlim(self.x1, self.x2)
        self.mcanvas.set_ylim(self.z2, self.z1)
        self.model_aspect = 1
        if self.tcanvas is not None:
            self.tcanvas.set_xlim(self.mcanvas.get_xlim())
        if self.dcanvas is not None:
            self.dcanvas.set_xlim(self.mcanvas.get_xlim())
        if self.ntcanvas is not None:
            self.ntcanvas.set_xlim(self.mcanvas.get_xlim())
        self.fig.subplots_adjust(top=0.99, left=-0.045, right=0.99, bottom=0.02, hspace=1.5)

        'ADJUST FRAME SIZING AND SET PROGRAM WINDOW'
        if self.t_canvas is False:
            self.tcanvas.set_visible(False)
        if self.d_canvas is False:
            self.dcanvas.set_visible(False)
        if self.nt_canvas is False:
            self.ntcanvas.set_visible(False)

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
        for i in range(0, self.layer_count):
            self.polygon_fills[i][0].set_alpha(self.layer_transparency)
        self.draw()

    def transparency_decrease(self, event):
        if self.layer_transparency >= 0.1:
            self.layer_transparency = self.layer_transparency - 0.1
            for i in range(0, self.layer_count):
                self.polygon_fills[i][0].set_alpha(self.layer_transparency)
            self.draw()
        else:
            pass

    # SAVE/LOAD MODEL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def save_model(self, event):
        """
        SAVE MODEL TO DISC IN .Pickle FORMAT
        TO ADD NEW OBJECTS ADD THE OBJECT NAME TO BOTH THE header AND model_params. THEN MODIFY THE load_model FUNCTION
        """
        save_file_dialog = wx.FileDialog(self, "Save model file", "", "", "Model files (*.model)|*.model", wx.FD_SAVE
                                         | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # USER CHANGED THEIR MIND

        # CREATE SAVE DICTIONARY
        self.save_dict = {}

        header = ['layer_colors', 'plotx_list', 'ploty_list', 'densities', 'reference_densities',
                  'boundary_lock_list', 'layer_lock_list', 'well_list', 'well_name_list',
                  'well_list_hidden', 'well_list_switch', 'segy_name_list', 'segy_dimension_list',
                  'obs_grav_list_save', 'obs_grav_name_list',
                  'obs_mag_list_save', 'obs_mag_name_list', 'earth_field', 'model_azimuth',
                  'mag_observation_elv', 'susceptibilities', 'angle_a', 'angle_b',
                  'tree_items', 'segy_file_list', 'absolute_densities', 'area', 'xp', 'zp',
                  'obs_gravity_data_for_rms', 'obs_mag_data_for_rms',
                  'obs_grav_colors', 'obs_mag_colors', 'obs_topo_colors', 'model_aspect',
                  'faults', 'fault_names_list',
                  'fault_counter', 'current_fault_index', 'fault_x_coords_list',
                  'fault_y_coords_list', 'fault_tree_items', 'fault_counter',
                  'outcrop_data_list', 'outcrop_data_list_save', 'outcrop_data_name_list', 'outcrop_data_color_list',
                  'outcrop_text_list',
                  'xy_name_list', 'xy_list_save', 'xy_color_list']

        model_params = [self.layer_colors, self.plotx_list, self.ploty_list, self.densities, self.reference_densities,
                        self.boundary_lock_list, self.layer_lock_list, self.well_list, self.well_name_list,
                        self.well_list_hidden, self.well_list_switch, self.segy_name_list, self.segy_dimension_list,
                        self.obs_grav_list_save, self.obs_grav_name_list,
                        self.obs_mag_list_save, self.obs_mag_name_list, self.earth_field, self.model_azimuth,
                        self.mag_observation_elv, self.susceptibilities, self.angle_a, self.angle_b,
                        self.tree_items, self.segy_file_list, self.absolute_densities, self.area, self.xp, self.zp,
                        self.obs_gravity_data_for_rms, self.obs_mag_data_for_rms,
                        self.obs_grav_colors, self.obs_mag_colors, self.obs_topo_colors, self.model_aspect,
                        self.faults, self.fault_names_list,
                        self.fault_counter, self.current_fault_index, self.fault_x_coords_list,
                        self.fault_y_coords_list, self.fault_tree_items, self.fault_counter,
                        self.outcrop_data_list, self.outcrop_data_list_save, self.outcrop_data_name_list,
                        self.outcrop_data_color_list, self.outcrop_text_list,
                        self.xy_name_list, self.xy_list_save, self.xy_color_list]

        for i in range(0, len(model_params)):
            try:
                self.save_dict[header[i]] = model_params[i]
            except IOError:
                print header[i]
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
            # PLOT OBSERVED DATA
            if self.xy_list_save != [[]]:
                self.plot_obs_xy()
            if self.obs_grav_list_save != [[]]:
                self.plot_obs_grav()
            if self.obs_mag_list_save != [[]]:
                self.plot_obs_mag()
            if self.obs_topo_list_save != [[]]:
                self.plot_obs_topo()
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # LOAD & PLOT SEGY
            self.plot_segy()
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # Set VARIABLE VALUES FROM LOADED DATA
            self.i = (len(self.densities) - 1)
            self.layer_count = (len(self.densities) - 1)
            self.density_input.SetValue(0.001 * self.densities[self.i])
            self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
            self.plotx = self.plotx_list[self.i]
            self.ploty = self.ploty_list[self.i]
            self.obs_grav_count = len(self.obs_grav_list_save) - 1
            self.obs_mag_list_count = len(self.obs_mag_list_save) - 1
            self.polygons = []
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            if self.layer_colors == [[]]:
                self.layer_colors = ['black' for _ in range(self.i + 1)]
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
            self.layer_lines = [[] for _ in range(self.i + 1)]
            self.polygon_fills = [[] for _ in range(self.i + 1)]
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # LOAD LAYERS
            self.load_layer_data()
            # ----------------------------------------------------------------------------------------------------------
            
            # ----------------------------------------------------------------------------------------------------------
            for i in range(0, self.fault_counter):
                # DRAW FAULTS
                self.faults[i] = self.mcanvas.plot(self.fault_x_coords_list[i], self.fault_y_coords_list[i],
                                                   color='k', linewidth=0.5, zorder=1, marker='s',
                                                   alpha=1.0)
            # CREATE NEW CURRENT FAULT GRAPHIC
            self.faultline, = self.mcanvas.plot([-50000, 0], [-49999, 0], marker='s', color='m', linewidth=0.75,
                                                alpha=1.0, zorder=2, picker=True)
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # DRAW OUTCROP DATA
            self.outcrop_data_count = 0
            for self.x in range(len(self.outcrop_data_name_list)):
                print self.outcrop_data_count
                if self.outcrop_data_name_list[self.x] == "None" or self.outcrop_data_name_list[self.x] == []:
                    self.outcrop_data_count += 1
                else:
                    # LOAD DATA INTO MODEL
                    print self.outcrop_data_count
                    self.open_outcrop_data('open')
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # LOAD WELL MENU --- SET WELL MENU - IDs start at 2000
            self.count = 0
            for x in range(len(self.well_name_list)):
                if self.well_name_list[x] == "None" or self.well_name_list[x] == []:
                    self.count += 1
                    continue
                else:
                    self.well_name_submenu = wx.Menu()
                    self.m_wells_submenu.Append(self.count + 3000, self.well_name_list[x], self.well_name_submenu)
                    self.well_name_submenu.Append(self.count + 2000, 'Hide/Show')
                    self.well_name_submenu.Append(self.count + 3000, 'Delete well')
                    self.Bind(wx.EVT_MENU, self.show_hide_well, id=self.count + 2000)
                    self.Bind(wx.EVT_MENU, self.delete_well, id=self.count + 3000)
                    self.count += 1
            self.well_count = len(self.well_name_list)
            # ----------------------------------------------------------------------------------------------------------
            
            # ----------------------------------------------------------------------------------------------------------
            # SET CURRENT NODE AS A OFF STAGE (PLACE HOLDER)
            self.current_node = self.mcanvas.scatter(-40000., 0., marker='o', color='r', zorder=10)
            # ----------------------------------------------------------------------------------------------------------
            
            # ----------------------------------------------------------------------------------------------------------
            # Set Fault PICKIGN SWTICH OFF (DEFAULT TO LAYER MODE)
            self.fault_picking_swtich = False
            # ----------------------------------------------------------------------------------------------------------
            
            # ----------------------------------------------------------------------------------------------------------
            # REFRESH SIZER POSITIONS
            self.Hide()
            self.Show()
            # ----------------------------------------------------------------------------------------------------------
            
            # ----------------------------------------------------------------------------------------------------------
            # UPDATE LAYER DATA AND PLOT
            self.draw()
            self.update_layer_data()
            self.run_algorithms()
            self.draw()
            self.draw_well()  # DRAW WELLS
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

    # PLOTTING FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def plot_segy(self):
        """ PLOT SEGY DATA"""
        self.segy_count = len(self.segy_file_list)

        if self.segy_count >= 1:
            for s in range(0, self.segy_count):
                if self.segy_name_list[s] != "NaN":
                    file_in = self.segy_file_list[s]
                    self.sx1, self.sx2, self.sz1, self.sz2 = self.segy_dimension_list[s]
                    section = read(file_in, unpack_trace_headers=False)
                    nsamples = len(section.traces[0].data)
                    ntraces = len(section.traces)
                    self.seis_data = plt.zeros((nsamples, ntraces))
                    for i, tr in enumerate(section.traces):
                        self.seis_data[:, i] = tr.data
                    del section
                    gc.collect()
                    del gc.garbage[:]

                    self.seis_axis = [self.sx1, self.sx2, self.sz2, self.sz1]
                    self.segy_plot_list.append(self.mcanvas.imshow(self.seis_data, vmin=self.gain_neg,
                                                                   vmax=self.gain, aspect='auto',
                                                                   extent=self.seis_axis, cmap=cm.gray,
                                                                   alpha=0.75))
                    del self.seis_data
                    gc.collect()
                    del gc.garbage[:]

                    self.segy_on = True
                    '''ADD SEGY_NAME TO SEGY MENU'''
                    # 1000 IS ADDED TO EVENT.ID TO PREVENT OVERLAP WITH GRAV EVENT.IDS
                    self.segy_name_submenu = wx.Menu()
                    self.m_segy_submenu.Append(s+1000, self.segy_name_list[s], self.segy_name_submenu)
                    self.segy_name_submenu.Append(s+1000, 'delete segy')
                    self.Bind(wx.EVT_MENU, self.remove_segy, id=s+1000)  # ID+1000 TO AVOID GRAV SUBMENU ID CONFLICT

                else:
                    continue

    def plot_obs_xy(self):
        """PLOT OBSERVED XY"""
        self.xy_list = [[]]
        for i in range(0, len(self.xy_list_save)-1):
            print("xy count =")
            print self.xy_count
            self.xy_list.append([])
            if self.xy_list_save[i] != []:
                self.xy = self.xy_list_save[i]
                self.xy_list[i] = self.mcanvas.scatter(self.xy[:, 0], self.xy[:, 1], marker='o',
                                                       color=self.xy_color_list[i], s=5, gid=self.xy_count)
                self.xy_list.append([])
                self.colors_index += 1
                self.xy_name = self.xy_name_list[i]
                self.xy_submenu = wx.Menu()
                self.m_xy_submenu.Append(self.xy_count, self.xy_name, self.xy_submenu)
                self.xy_submenu.Append(self.xy_count+4000, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_xy, id=self.xy_count+4000)
                self.xy_count += 1
            else:
                self.xy_count += 1

    def plot_obs_topo(self):
        """PLOT OBSERVED TOPOGRAPHY"""
        self.obs_topo_list = [[]]
        for x in range(0, len(self.obs_topo_name_list)):
            self.obs_topo_list.append([])
            if self.obs_topo_name_list[x]:
                self.topo = self.obs_topo_list_save[x]
                self.obs_topo_list[x] = self.dcanvas.scatter(self.topo[:, 0], self.topo[:, 1], marker='o',
                                                             color=self.obs_topo_colors[x], s=5,
                                                             gid=self.obs_topo_count)
                self.obs_topo_list.append([])
                self.obs_name = self.obs_topo_name_list[x]
                self.obs_submenu = wx.Menu()
                self.m_topo_submenu.Append(self.obs_topo_count, self.obs_name, self.obs_submenu)
                self.obs_submenu.Append(self.obs_topo_count+10000, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_topo, id=self.obs_topo_count+10000)
                self.obs_topo_count += 1
                self.obs_topo_switch = True
            else:
                self.obs_topo_count += 1

    def plot_obs_grav(self):
        """PLOT OBSERVED GRAVITY"""
        self.obs_grav_list = [[]]
        for x in range(0, len(self.obs_grav_name_list)):
            self.obs_grav_list.append([])
            if self.obs_grav_name_list[x]:
                self.obs_grav = self.obs_grav_list_save[x]
                self.obs_grav_list[x] = self.dcanvas.scatter(self.obs_grav[:, 0], self.obs_grav[:, 1], marker='o',
                                                             color=self.obs_grav_colors[x], s=5,
                                                             gid=self.obs_grav_count)
                self.obs_grav_list.append([])
                self.obs_name = self.obs_grav_name_list[x]
                self.obs_submenu = wx.Menu()
                self.m_obs_g_submenu.Append(self.obs_grav_count, self.obs_name, self.obs_submenu)
                self.obs_submenu.Append(self.obs_grav_count+11000, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=self.obs_grav_count+11000)
                self.obs_grav_count += 1
                self.grav_obs_switch = True
            else:
                self.obs_grav_count += 1

    def plot_obs_mag(self):
        """PLOT OBSERVED MAG"""
        self.obs_mag_list = [[]]
        for x in range(0, len(self.obs_mag_name_list)):
            self.obs_mag_list.append([])
            if self.obs_mag_name_list[x]:
                self.obs_mag = self.obs_mag_list_save[x]
                self.obs_mag_list[x] = self.ntcanvas.scatter(self.obs_mag[:, 0], self.obs_mag[:, 1], marker='o',
                                                             color=self.obs_mag_colors[x], s=5,
                                                             gid=self.obs_mag_count)
                self.obs_mag_list.append([])
                self.obs_mag_name = self.obs_mag_name_list[x]
                self.obs_submenu = wx.Menu()
                self.m_obs_mag_submenu.Append(self.obs_mag_count, self.obs_mag_name, self.obs_submenu)
                self.obs_submenu.Append(self.obs_mag_count+12000, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=self.obs_mag_count+12000)
                self.obs_mag_count += 1
                self.mag_obs_switch = True
            else:
                self.obs_mag_count += 1

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

        # GET XY DATA NAME
        self.xy_name_list.append([])
        self.xy_name_list[self.xy_count] = str(self.xy_name)

        #GET XY_COLOR
        self.xy_color_list.append([])
        self.xy_color_list[self.xy_count] = str(self.xy_color)

        # LOAD DATA
        self.xy = np.genfromtxt(xy_input_file, delimiter=' ', dtype=float)
        self.xy_list_save.append([])
        self.xy_list_save[self.xy_count] = self.xy
        self.xy_list.append([])

        # PLOT DATA IN MODEL
        self.xy_list[self.xy_count] = self.mcanvas.scatter(self.xy[:, 0], self.xy[:, 1], marker='o',
                                                           color=self.xy_color_list[self.xy_count], s=3,
                                                           gid=self.xy_count)

        #  APPEND NEW DATA MENU TO 'XY data MENU'
        self.obs_submenu = wx.Menu()
        self.m_xy_submenu.Append(self.xy_count, self.xy_name, self.obs_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.obs_submenu.Append(4000+self.xy_count, 'delete observed data')

        # BIND TO DEL XY FUNC
        self.Bind(wx.EVT_MENU, self.delete_xy, id=4000+self.xy_count)

        # INCREMENT XY COUNTER
        self.xy_count += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.draw()

    def delete_xy(self, event):
        """"DELETE OBSERVED XY DATA NB: ID's start at 4000"""
        id = event.Id-4000
        self.m_xy_submenu.DestroyItem(id)
        self.xy_list[id].set_visible(False)
        self.xy_list[id] = []
        self.xy_name_list[id] = []
        self.xy_color_list[id] = []
        self.xy_list_save[id] = []

        # UPDATE GMG
        self.update_layer_data()
        self.draw()

    def delete_all_xy(self, event):
        for i in range(0, self.xy_count):
            try:
                self.m_xy_submenu.DestroyId(i)
                self.xy_list[i].set_visible(False)
            except:
                continue
        # RESET DATA LISTS
        self.xy = []
        self.xy_list = [[]]
        self.xy_list_save = [[]]
        self.xy_name_list = []
        self.xy_color_list = []
        self.xy_count = 0

        # UPDATE GMG
        self.update_layer_data()
        self.draw()

    # TOPOGRAPHY DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_topo(self, event):
        self.load_window = LoadObservedDataFrame(self, -1, 'Load observed data', 'topography')
        self.load_window.Show(True)

    def open_obs_t(self):
        """LOAD OBSERVE TOPO DATA: IDs start at 10000"""
        obs_t_input = self.load_window.file_path
        self.obs_name = self.load_window.observed_name
        self.color = self.load_window.color_picked

        self.obs_topo_name_list.append([])
        self.obs_topo_name_list[self.obs_topo_count] = str(self.obs_name)
        self.obs_topo_colors.append([])
        self.obs_topo_colors[self.obs_topo_count] = str(self.color)

        self.topo = np.genfromtxt(obs_t_input, delimiter=' ', dtype=float)
        self.obs_topo_list_save[self.obs_topo_count] = self.topo
        self.obs_topo_list_save.append([])
        self.obs_topo_list[self.obs_topo_count] = self.tcanvas.scatter(self.topo[:, 0], self.topo[:, 1], marker='o',
                                                                       color=self.colors[self.colors_index],
                                                                       s=5, gid=self.obs_grav_count)
        self.colors_index += 1
        self.obs_topo_list.append([])

        #  APPEND NEW DATA MENU TO 'Topo data MENU'
        self.topo_submenu = wx.Menu()
        self.m_topo_submenu.Append(self.obs_topo_count, self.obs_name, self.topo_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.topo_submenu.Append(10000+self.obs_topo_count, 'delete observed data')

        # BIND TO DEL TOPO FUNC
        self.Bind(wx.EVT_MENU, self.delete_topo, id=10000+self.obs_topo_count)

        self.obs_topo_count += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def delete_topo(self, event):
        """DELETE AN OBSERVED TOPO DATA RECORD"""
        self.m_topo_submenu.DestroyItem(event.Id-10000)
        self.obs_topo_name_list[event.Id-10000] = []
        self.obs_topo_list[event.Id-10000].set_visible(False)
        self.obs_topo_list[event.Id-10000] = []
        self.obs_topo_list_save[event.Id-10000] = []

        # UPDATE MODEL
        self.update_layer_data()
        self.draw()

    # GRAVITY DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_obs_g(self, event):
        self.load_window = LoadObservedDataFrame(self, -1, 'Load observed data', 'gravity')
        self.load_window.Show(True)

    def open_obs_g(self):
        """LOAD OBSERVE GRAVITY DATA: IDs start at 11000"""
        obs_g_input = self.load_window.file_path
        self.obs_name = self.load_window.observed_name
        self.color = self.load_window.color_picked

        self.obs_grav_name_list.append([])
        self.obs_grav_name_list[self.obs_grav_count] = str(self.obs_name)
        self.obs_grav_colors.append([])
        self.obs_grav_colors[self.obs_grav_count] = str(self.color)

        # DRAW LOADED DATA IN GRAVITY FRAME
        self.obs_grav = np.genfromtxt(obs_g_input, delimiter=' ', dtype=float)
        self.obs_gz = self.obs_grav[:, 1]
        self.obs_grav_list_save.append([])
        self.obs_grav_name_list.append([])
        self.obs_grav_list_save[self.obs_grav_count] = self.obs_grav
        self.obs_grav_list[self.obs_grav_count] = self.dcanvas.scatter(self.obs_grav[:, 0], self.obs_grav[:, 1],
                                                                       marker='o', color=self.color,
                                                                       s=5, gid=self.obs_grav_count)
        self.colors_index += 1
        self.obs_grav_list.append([])
        self.grav_obs_switch = True

        #  APPEND NEW DATA MENU TO 'GRAV data MENU'
        self.grav_submenu = wx.Menu()
        self.m_obs_g_submenu.Append(self.obs_grav_count, self.obs_name, self.grav_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.grav_submenu.Append(11000+self.obs_grav_count, 'delete observed data')

        # BIND TO DEL GRAV FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000+self.obs_grav_count)

        # INCREMENT GRAV COUNTER
        self.obs_grav_count += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def delete_obs_grav(self, event):
        self.m_obs_g_submenu.DestroyItem(event.Id-11000)
        self.obs_grav_name_list[event.Id-11000] = []
        self.obs_grav_list[event.Id-11000].set_visible(False)
        self.obs_grav_list[event.Id-11000] = []
        self.obs_grav_list_save[event.Id-11000] = []

        # UPDATE MODEL
        self.update_layer_data()
        self.draw()

    def delete_all_obs_grav(self, event):
        for x in range(0, self.obs_grav_count):
            try:
                self.m_obs_g_submenu.DestroyItem(x)
                self.obs_grav_list[x].set_visible(False)
            except IndexError:
                continue
        self.grav_obs_switch = False
        self.obs_grav = []
        self.obs_grav_list = [[]]
        self.obs_grav_list_save = [[]]
        self.obs_grav_name_list = []
        self.obs_grav_count = 0
        self.update_layer_data()
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
        """LOAD OBSERVE MAGNETIC DATA: IDs start at 12000"""
        # GET VALUES FROM LOAD WINDOW
        obs_m_input = self.load_window.file_path
        self.obs_mag_name = self.load_window.observed_name
        self.color = self.load_window.color_picked

        self.obs_mag_name_list.append([])
        self.obs_mag_name_list[self.obs_mag_count] = str(self.obs_mag_name)
        self.obs_mag_colors.append([])
        self.obs_mag_colors[self.obs_mag_count] = str(self.color)

        self.obs_mag = np.genfromtxt(obs_m_input, delimiter=' ', dtype=float)
        self.obs_mz = self.obs_mag[:, 1]
        self.obs_mag_list_save[self.obs_mag_count] = self.obs_mag
        self.obs_mag_list_save.append([])
        self.obs_mag_list[self.obs_mag_count] = self.ntcanvas.scatter(self.obs_mag[:, 0], self.obs_mag[:, 1],
                                                                      marker='o', color=self.colors[self.colors_index],
                                                                      s=5, gid=self.obs_grav_count)
        self.colors_index += 1
        self.obs_mag_list.append([])
        self.mag_obs_switch = True

        #  APPEND NEW DATA MENU TO 'Mag data MENU'
        self.obs_mag_submenu = wx.Menu()
        self.m_obs_mag_submenu.Append(self.obs_mag_count, self.obs_mag_name, self.obs_mag_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.obs_mag_submenu.Append(12000+self.obs_mag_count, 'delete observed data')

        # BIND TO DEL MAG FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12000+self.obs_mag_count)

        self.obs_mag_count += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def delete_obs_mag(self, event):
        self.obs_mag_submenu.DestroyItem(event.Id-12000)
        self.obs_mag_list[event.Id].set_visible(False)
        self.obs_mag_list[event.Id-12000] = []
        self.obs_mag_name_list[event.Id-12000] = []
        self.obs_mag_list_save[event.Id-12000] = []

        # UPDATE MODEL
        self.update_layer_data()
        self.draw()

    def delete_all_obs_mag(self, event):
        for x in range(0, self.obs_mag_count):
            try:
                self.m_obs_mag_submenu.DestroyItem(x)
                self.obs_mag_list[x].set_visible(False)
            except:
                continue
        self.mag_obs_switch = False
        self.obs_mag = []
        self.obs_mag_list = [[]]
        self.obs_mag_list_save = [[]]
        self.obs_mag_name_list = []
        self.obs_mag_count = 0
        self.update_layer_data()
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
                return  # the user changed idea...
            file_in = open_file_dialog.GetPath()
            self.segy_dimension_list.append(self.d)
            self.segy_file_list.append(file_in)
            self.segy_name_list.append(self.segy_name)
            section = read(file_in, unpack_trace_headers=False)
            nsamples = len(section.traces[0].data)
            ntraces = len(section.traces)
            self.seis_data = plt.zeros((nsamples, ntraces))
            for i, tr in enumerate(section.traces):
                self.seis_data[:, i] = tr.data
            del section
            self.seis_axis = [self.sx1, self.sx2, self.sz2, self.sz1]
            self.segy_plot_list.append(self.mcanvas.imshow(self.seis_data, vmin=self.gain_neg, vmax=self.gain,
                                                           aspect='auto', extent=self.seis_axis, cmap=cm.gray,
                                                           alpha=0.75))
            del self.seis_data
            self.segy_on = True
            self.update_layer_data()

            # ADD SEGY_NAME TO SEGY MENU
            # 1000 IS ADDED TO EVENT.ID TO PREVENT OVERLAP WITH GRAV EVENT.IDS
            self.segy_name_submenu = wx.Menu()
            self.m_segy_submenu.Append(self.segy_count + 1000, self.segy_name_list[self.segy_count],
                                       self.segy_name_submenu)
            self.segy_name_submenu.Append(self.segy_count + 1000, 'delete segy')
            self.Bind(wx.EVT_MENU, self.remove_segy,
                      id=self.segy_count + 1000)  # ID IS SET TO +1000 TO AVOID CONFLICT WITH GRAV SUB MENU IDS

            self.segy_count = self.segy_count + 1
        except AttributeEditor:
            load_error = MessageDialog(self, -1, "SEGY LOAD ERROR", "segy load error")

    def remove_all_segy(self, event):
        self.segy_plot_list = []
        self.segy_name_list = []
        self.segy_dimension_list = []
        self.segy_dimension_list = []
        self.segy_file_list = []
        self.segy_count = 0
        self.draw()

    def remove_segy(self, event):
        # 1000 IS TAKEN FROM EVENT.ID TO PREVENT OVERLAP WITH GRAV EVENT.IDS
        if self.segy_on:
            img = self.segy_plot_list[event.Id - 1000]
            img.remove()
            self.segy_name_list[event.Id - 1000] = "NaN"
            self.segy_plot_list[event.Id - 1000] = 0
            self.m_segy_submenu.DestroyId(event.Id)
            self.draw()

    def segy_color_adjustment(self, event):
        if event.Id == 901:
            if self.segy_on is False:
                return 0
            else:
                for s in range(0, len(self.segy_plot_list)):
                    if self.segy_plot_list[s]:
                        self.segy_plot_list[s].set_cmap(cm.gray)
                self.draw()
        else:
            if event.Id == 902:
                if self.segy_on is False:
                    return 0
                else:
                    for s in range(0, len(self.segy_plot_list)):
                        if self.segy_plot_list[s]:
                            self.segy_plot_list[s].set_cmap(cm.seismic)
                    self.draw()

    def gain_increase(self, event):
        if self.segy_on is False:
            return 0
        elif self.gain >= 0.5:
            self.gain = self.gain - 0.5
            self.gain_neg = -self.gain
            for s in range(0, len(self.segy_plot_list)):
                if self.segy_plot_list[s]:
                    self.segy_plot_list[s].set_clim(vmax=self.gain)
                    self.segy_plot_list[s].set_clim(vmin=self.gain_neg)
            self.draw()
        else:
            return 0

    def gain_decrease(self, event):
        if self.segy_on is False:
            return 0
        elif self.gain >= 0.5:
            self.gain = self.gain + 0.5
            self.gain_neg = -self.gain
            for s in range(0, len(self.segy_plot_list)):
                if self.segy_plot_list[s]:
                    self.segy_plot_list[s].set_clim(vmax=self.gain)
                    self.segy_plot_list[s].set_clim(vmin=self.gain_neg)
            self.draw()
        else:
            return 0

    # WELL DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_well(self, event):
        open_file_dialog = wx.FileDialog(self, "Open Observed file", "", "", "All files (*.*)|*.*",
                                         wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if open_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THEIR MIND
        else:
            well_in = open_file_dialog.GetPath()
            well_name_box = wx.TextEntryDialog(None, 'Please provide a name for the new well record:')
            answer = well_name_box.ShowModal()

        self.well_name = well_name_box.GetValue()
        self.well_name_list.append(str(self.well_name))

        with open(well_in, 'r') as f:
            input_well_data = [line.strip().split(' ') for line in f]
        self.well_list.append(input_well_data)

        # CREATE FILE MENU DATA
        self.well_name_submenu = wx.Menu()
        self.m_wells_submenu.Append(int(self.well_count) + 3000, self.well_name, self.well_name_submenu)
        self.well_name_submenu.Append(int(self.well_count) + 2000, 'Hide/Show')
        self.well_name_submenu.Append(int(self.well_count) + 3000, 'Delete well')
        self.Bind(wx.EVT_MENU, self.show_hide_well, id=int(self.well_count) + 2000)
        self.Bind(wx.EVT_MENU, self.delete_well, id=int(self.well_count) + 3000)
        self.well_count += 1

        # DRAW WELL
        well_data = np.array(self.well_list[self.well_count - 1][1:])  # CREATE NP ARRAY WITHOUT HEADER INFO
        y1 = well_data[0][1].astype(float)
        y2 = well_data[-1][-1].astype(float)
        well_x_location = well_data[1][1].astype(float)
        wellx = (well_x_location, well_x_location)
        welly = (y1, y2)
        self.wells.append([])
        self.wells[self.well_count - 1] = self.mcanvas.plot(wellx, welly, linestyle='-', linewidth='2', color='black')

        self.well_name_text.append([])

        # PLOT WELL NAME
        self.well_name_text[self.well_count - 1] = self.mcanvas.annotate(self.well_name, xy=(well_x_location, -0.5),
                                                                         xytext=(well_x_location, -0.5),
                                                                         fontsize=self.well_textsize, weight='bold',
                                                                         horizontalalignment='center', color='black',
                                                                         bbox=dict(boxstyle="round,pad=.2", fc="0.8"),
                                                                         clip_on=True)
        self.well_name_text.append([])

        # PLOT WELL HORIZONS
        # SET EMPTY ARRAYS TO FILL WITH LABELS AND HORIZONS
        self.well_labels.append([])
        self.horizons.append([])
        self.well_labels[self.well_count - 1] = [None] * len(well_data)
        self.horizons[self.well_count - 1] = [None] * len(well_data)
        for i in range(2, len(well_data)):
            y = [well_data[i][1].astype(float), well_data[i][1].astype(float)]
            x = [well_data[1][1].astype(float) - 1, well_data[1][1].astype(float) + 1]

            # PLOT HORIZON LINE
            self.horizons[self.well_count - 1][i] = self.mcanvas.plot(x, y, linestyle='-', linewidth='2', color='black')
            horizon_y_pos = well_data[i][1].astype(float)
            horizon = well_data[i][0].astype(str)

            # ALTERNATE POSITION OF ODDs/EVENs TO TRY AND AVOID OVERLAP
            if i % 2 == 0:
                horizon_x_pos = well_data[1][1].astype(float) - 1.05
                self.well_labels[self.well_count - 1][i] = self.mcanvas.annotate(horizon,
                                                                                 xy=(horizon_x_pos, horizon_y_pos),
                                                                                 xytext=(horizon_x_pos, horizon_y_pos),
                                                                                 fontsize=self.well_textsize,
                                                                                 weight='bold',
                                                                                 horizontalalignment='left',
                                                                                 verticalalignment='top',
                                                                                 color='black',
                                                                                 bbox=dict(boxstyle="round,pad=.4",
                                                                                           fc="0.8", ec='None'),
                                                                                 clip_on=True)
            else:
                horizon_x_pos = well_data[1][1].astype(float) + 1.05
                self.well_labels[self.well_count - 1][i] = self.mcanvas.annotate(horizon,
                                                                                 xy=(horizon_x_pos, horizon_y_pos),
                                                                                 xytext=(horizon_x_pos, horizon_y_pos),
                                                                                 fontsize=self.well_textsize,
                                                                                 weight='bold',
                                                                                 horizontalalignment='right',
                                                                                 verticalalignment='bottom',
                                                                                 color='black',
                                                                                 bbox=dict(boxstyle="round,pad=.4",
                                                                                           fc="0.8", ec='None'),
                                                                                 clip_on=True)
        self.well_labels_list[self.well_count - 1] = self.well_labels
        self.well_horizons_list[self.well_count - 1] = self.horizons
        self.well_horizons_list.append([])
        self.well_labels_list.append([])

        # UPDATE GMG
        self.update_layer_data()

    def draw_well(self):
        """DRAW WELLS ON MODEL CANVAS"""

        # CREATE EMPTY ARRAYS TO FILL WITH WELL DATA
        self.wells = [[]] * len(self.well_list)
        self.well_name_text = [[]] * len(self.well_list)
        self.well_labels = [[]] * len(self.well_list)
        self.horizons = [[]] * len(self.well_list)

        for w in range(0, len(self.well_list)):
            if self.well_list[w] == "None" or self.well_list[w] == []:
                self.well_horizons_list.append([])
                self.well_labels_list.append([])
                continue
            else:

                self.wd = np.array(self.well_list[w][1:])  # CREATE NP ARRAY WITHOUT HEADER INFO
                y1 = float(self.wd[0][1])
                y2 = float(self.wd[-1][-1])
                well_x_location = float(self.wd[1][1])
                wellx = (well_x_location, well_x_location)
                welly = (y1, y2)

                self.wells[w] = self.mcanvas.plot(wellx, welly, linestyle='-', linewidth='2', color='black')

                # PLOT WELL NAME
                well_name = self.well_name_list[w]
                self.well_name_text[w] = self.mcanvas.annotate(well_name, xy=(well_x_location, -0.5),
                                                               xytext=(well_x_location, -0.5),
                                                               fontsize=self.well_textsize, weight='bold',
                                                               horizontalalignment='center', color='black',
                                                               bbox=dict(boxstyle="round,pad=.2", fc="0.8"),
                                                               clip_on=True)

                # PLOT WELL HORIZONS
                # SET EMPTY ARRAYS TO FILL WITH LABELS AND HORIZONS
                self.well_labels[w] = [None] * len(self.wd)
                self.horizons[w] = [None] * len(self.wd)
                for i in range(2, len(self.wd)):
                    y = [self.wd[i][1].astype(float), self.wd[i][1].astype(float)]
                    x = [self.wd[1][1].astype(float) - 1, self.wd[1][1].astype(float) + 1]

                    # PLOT HORIZON LINE
                    self.horizons[w][i] = self.mcanvas.plot(x, y, linestyle='-', linewidth='2', color='black')
                    horizon_y_pos = self.wd[i, 1].astype(float)
                    horizon = self.wd[i, 0].astype(str)

                    # ALTERNATE POSITION OF ODDs/EVENs TO TRY AND AVOID OVERLAP
                    if i % 2 == 0:
                        horizon_x_pos = self.wd[1, 1].astype(float) - 1.05
                        self.well_labels[w][i] = self.mcanvas.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                       xytext=(horizon_x_pos, horizon_y_pos),
                                                                       fontsize=self.well_textsize,
                                                                       weight='bold', horizontalalignment='left',
                                                                       verticalalignment='top',
                                                                       color='black',
                                                                       bbox=dict(boxstyle="round,pad=.4", fc="0.8",
                                                                                 ec='None'),
                                                                       clip_on=True)
                    else:
                        horizon_x_pos = self.wd[1][1].astype(float) + 1.05
                        self.well_labels[w][i] = self.mcanvas.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                       xytext=(horizon_x_pos, horizon_y_pos),
                                                                       fontsize=self.well_textsize,
                                                                       weight='bold', horizontalalignment='right',
                                                                       verticalalignment='bottom',
                                                                       color='black',
                                                                       bbox=dict(boxstyle="round,pad=.4", fc="0.8",
                                                                                 ec='None'),
                                                                       clip_on=True)
                self.well_labels_list[w] = self.well_labels
                self.well_horizons_list[w] = self.horizons
                self.well_horizons_list.append([])
                self.well_labels_list.append([])
        self.draw()

    def show_hide_well(self, event):

        if self.wells[event.Id - 2000][0].get_visible():
            # HIDE WELL
            self.wells[event.Id - 2000][0].set_visible(False)
            self.well_name_text[event.Id - 2000].set_visible(False)
            horizons = self.horizons[event.Id - 2000]
            for i in range(2, len(horizons)):
                horizons[i][0].set_visible(False)
            labels = self.well_labels[event.Id - 2000]
            for i in range(2, len(labels)):
                labels[i].set_visible(False)

            self.update_layer_data()
        else:
            # SHOW WELL
            self.wells[event.Id - 2000][0].set_visible(True)
            self.well_name_text[event.Id - 2000].set_visible(True)
            horizons = self.horizons[event.Id - 2000]
            for i in range(2, len(horizons)):
                horizons[i][0].set_visible(True)
            labels = self.well_labels[event.Id - 2000]
            for i in range(2, len(labels)):
                labels[i].set_visible(True)

            self.update_layer_data()

    def delete_well(self, event):
        """REMOVE PLOT GRAPHICS"""
        self.wells[event.Id - 3000][0].set_visible(False)
        self.well_name_text[event.Id - 3000].set_visible(False)
        horizons = self.horizons[event.Id - 3000]
        for i in range(2, len(horizons)):
            horizons[i][0].set_visible(False)
            horizons[i] = "None"
        labels = self.well_labels[event.Id - 3000]
        for i in range(2, len(labels)):
            labels[i].set_visible(False)
            labels[i] = "None"

        # SET DATA TO NONE
        self.well_list[event.Id - 3000] = "None"
        self.well_name_list[event.Id - 3000] = "None"
        self.wells[event.Id - 3000][0] = "None"
        self.well_name_text[event.Id - 3000] = "None"

        # REMOVE SUBMENU
        self.m_wells_submenu.DestroyItem(event.Id)
        self.update_layer_data()

    # GEOLOGY OUTCROP DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_outcrop_data(self, event):
        """ LOAD & PLOT SURFACE CONTACT DATA. E.G. GEOLOGICAL CONTACTS  NB: ID's start at 13000"""
        self.load_window = LoadObservedDataFrame(self, -1, 'Load observed data', 'outcrop')
        self.load_window.Show(True)

    def open_outcrop_data(self, type):
        """OPEN THE OUTCROP DATA SELECTED BY THE USER USING THE load_outcrop_data FUNC"""
        if type == 'load':
            outcrop_input_file = self.load_window.file_path
            self.outcrop_data_name = self.load_window.observed_name
            self.outcrop_data_color = self.load_window.color_picked
            self.outcrop_data = np.genfromtxt(outcrop_input_file, autostrip=True, delimiter=' ', dtype=str,
                                              comments='#')
            self.outcrop_data_list.append([])
        else:
            self.outcrop_data = self.outcrop_data_list_save[self.x]
            self.outcrop_data_name =  self.outcrop_data_name_list[self.x]
            self.outcrop_data_color = self.outcrop_data_color_list[self.x]

        # GET XY DATA NAME
        self.outcrop_data_name_list.append([])
        self.outcrop_data_name_list[self.outcrop_data_count] = str(self.outcrop_data_name)

        #GET XY_COLOR
        self.outcrop_data_color_list.append([])
        self.outcrop_data_color_list[self.outcrop_data_count] = str(self.outcrop_data_color)

        # LOAD DATA
        self.outcrop_data_list_save.append([])
        self.outcrop_data_list_save[self.outcrop_data_count] = self.outcrop_data

        # PLOT MARKERS IN MODEL
        self.contact_labels = [None] * len(self.outcrop_data)
        self.contacts = [None] * len(self.outcrop_data)
        lines = [None] * len(self.outcrop_data)
        for i in range(len(self.outcrop_data)):
            x1 = self.outcrop_data[i, 0].astype(float)
            y1 = self.outcrop_data[i, 1].astype(float)
            y2 = self.outcrop_data[i, 2].astype(float)
            x = (x1, x1)
            y = (y1, y2)
            lines[i] = self.mcanvas.plot(x, y, linestyle='-', linewidth='2', color='black')

        self.outcrop_data_list.append([])
        self.outcrop_data_list[self.outcrop_data_count] = lines

        # DRAW TEXT
        self.text_labels = [None] * len(self.outcrop_data)
        text = zip(self.outcrop_data[:, 0].astype(float), self.outcrop_data[:, 1].astype(float),
                   self.outcrop_data[:, 3].astype(str))
        for i in range(len(self.outcrop_data)):

            # ALTERNATE POSITION OF ODDs/EVENs To TRY AND AVOID OVERLAP
            if i % 2 == 0:
                self.text_labels[i] = self.mcanvas.annotate(text[i][2], xy=(text[i][0], text[i][1]),
                                                            xytext=(text[i][0], text[i][1]),
                                                            fontsize=self.well_textsize,
                                                            weight='regular', horizontalalignment='right',
                                                            verticalalignment='bottom',
                                                            color='black',
                                                            bbox=dict(boxstyle="round,pad=.4", fc="0.8", ec='None'),
                                                            clip_on=True)
            else:
                self.text_labels[i] = self.mcanvas.annotate(text[i][2], xy=(text[i][0], text[i][1]),
                                                            xytext=(text[i][0], text[i][1]),
                                                            fontsize=self.well_textsize,
                                                            weight='regular', horizontalalignment='left',
                                                            verticalalignment='top',
                                                            color='black',
                                                            bbox=dict(boxstyle="round,pad=.4", fc="0.8", ec='None'),
                                                            clip_on=True)

        self.outcrop_text_list.append([])
        self.outcrop_text_list[self.outcrop_data_count] = self.text_labels


        #  APPEND NEW DATA MENU TO 'OUTCROP DATA MENU'
        self.outcrop_submenu = wx.Menu()
        self.m_outcrop_submenu.Append(self.outcrop_data_count, self.outcrop_data_name, self.outcrop_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.outcrop_submenu.Append(13000+self.outcrop_data_count, 'delete observed data')

        # BIND TO DEL XY FUNC
        self.Bind(wx.EVT_MENU, self.delete_outcrop_data, id=13000+self.outcrop_data_count)

        # INCREMENT CONTACT COUNT
        self.outcrop_data_count += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.draw()

    def delete_outcrop_data(self, event):
        """"DELETE OUTCROP DATA NB: ID's start at 13000"""
        id = event.Id-13000
        self.m_outcrop_submenu.DestroyItem(id)

        for i in range(0, len(self.outcrop_data_list[id])):
            self.outcrop_data_list[id][i][0].set_visible(False)

        for i in range (0, len(self.outcrop_text_list[id])):
            self.outcrop_text_list[id][i].set_visible(False)

        self.outcrop_data_list[id] = []
        self.outcrop_data_name_list[id] = []
        self.outcrop_data_color_list[id] = []
        self.outcrop_data_list_save[id] = []
        self.outcrop_text_list[id] = []

        self.update_layer_data()
        self.draw()

    # LAYER & NODE CONTROLS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def b_node_set_button(self, event):
        """
        CHECK IF A NODE FROM THE CURRENT LAYER IS SELECTED; IF NOT, THEN SKIP THIS PART AND ONLY UPDATE ATTRIBUTES
        """

        new_x = float(self.x_input.GetValue())
        new_y = float(self.y_input.GetValue())

        if self.i == self.node_layer_reference:
            xt = np.array(self.plotx)
            yt = np.array(self.ploty)

            if self.boundary_lock_list[self.i] == 0 and self.index_node is not None:
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
            elif self.boundary_lock_list[self.i] == 1:
                if new_y <= 0:
                    xt[self.index_node] = new_x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                else:
                    xt[self.index_node] = new_x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = new_y  # REPLACE OLD Y WITH NEW Y

            # DEAL WITH PINCHED NODE
            if self.pinch_switch != 0:
                for k in range(0, len(self.index_arg2_list)):
                    if self.index_arg2_list[k] is not None:
                        next_x_list = self.plotx_list[k]
                        next_y_list = self.ploty_list[k]  # GET THE NODE LIST OF THE NEXT LAYER
                        next_x_list[self.index_arg2_list[k]] = new_x
                        next_y_list[self.index_arg2_list[k]] = new_y  # REPLACE THE PINCHED NODE WITH THE NEW NODE
                        self.plotx_list[k] = next_x_list
                        self.ploty_list[k] = next_y_list  # OVERWRITE THE NODE LIST WITH UPDATED LIST

            self.plotx = xt
            self.ploty = yt
            self.polyline.set_data(self.plotx, self.ploty)

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
        GET THE INDEX VALUE OF THE NODE UNDER POINT IF IT IS WITHIN NODE_CLICK_LIMIT TOLERANCE OF CLICK
        """

        self.node_layer_reference = self.i
        xyt = self.polyline.get_xydata()
        xt = xyt[:, 0]
        yt = xyt[:, 1]
        d = np.sqrt((xt - event.xdata) ** 2 + (yt - event.ydata) ** 2)
        self.index_arg = np.argmin(d)
        if d[self.index_arg] >= self.node_click_limit:
            return None, None
        else:
            # CHECK IF NODE IS A PINCHED POINT, IF YES FIND NODE OF ABOVE OR BELOW LAYER
            self.pinch_switch = 0

            # CREATE LIST OF NONES SAME LENGTH AS NUMBER OF LAYERS IN MODEL
            self.index_arg2_list = [None] * (self.layer_count + 1)

            for x in range(0, self.layer_count + 1):  # LOOP THROUGH ALL LAYERS TO CHECK FOR PINCHED NODES
                if x == self.i:  # SKIP CURRENT LAYER
                    pass
                node_list_x, node_list_y = self.plotx_list[x], self.ploty_list[x]

                for i in range(0, len(node_list_x)):
                    if node_list_x[i] == xt[self.index_arg] and node_list_y[i] == yt[self.index_arg]:
                        # IF ONE OF THE NODES FROM LIST IS EQUAL TO A NODE FROM THE OTHER LAYER, THEN RETURN THE INDEX
                        self.index_arg2_list[x] = i
                        self.pinch_switch = 1

            'return node index'
            return self.index_arg, self.index_arg2_list

    def get_fault_node_under_point(self, event):
        """GET THE INDEX VALUE OF THE NODE UNDER POINT IF IT IS WITHIN NODE_CLICK_LIMIT TOLERANCE OF CLICK"""
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
        """ WHEN THE LEFT MOUSE BUTTON IS PRESSED"""

        if event.inaxes is None:
            return  # CLICK IS OUTSIDE MODEL FRAME
        if event.button != 1:
            return

        if self.fault_picking_switch is False and self.capture is False:
            # THEN GMG IS IN LAYER MODE

            # GET THE NODE CLOSEST TO THE CLICK AND ANY PINCHED NODES
            self.index_node, self.index_arg2_list = self.get_node_under_point(event)
            if self.index_node is None:
                return

            xyt = self.polyline.get_xydata()
            xt, yt = xyt[:, 0], xyt[:, 1]
            self.x_input.SetValue(xt[self.index_node])
            self.y_input.SetValue(yt[self.index_node])

            # IF PINCH == TRUE, THEN PINCH THE NODE TO NEXT NODE
            if self.pinch:
                # GET THE NODE NUMBER AND LAYER NUMBER AND PLACE THEM IN "PINCH_NODE_LIST"
                if self.pinch_count == 0:
                    self.plotx = self.plotx_list[self.i]
                    self.ploty = self.ploty_list[self.i]
                    x1 = np.array(self.plotx)
                    y1 = np.array(self.ploty)
                    self.index_node = self.get_node_under_point(event)
                    self.pinch_node_list[self.pinch_count] = self.index_node[0]
                    self.pinch_node_list[self.pinch_count + 1] = self.i
                    self.pinch_count = + 1
                    # SET THE X AND Y OF THE FIRST NODE AS THAT OF THE SECOND NODE
                else:
                    # SET THE SECOND NODE X AND Y
                    self.plotx2 = self.plotx_list[self.i]
                    self.ploty2 = self.ploty_list[self.i]
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

        elif self.fault_picking_switch is True:
            # GMG IS IN FAULT MODE

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

        elif self.capture is True:
            # COORDINATE CAPTURE MODE
            pass

    def toogle_fault_mode(self, event):
        """SWITCH FAULT PICKING MODE ON AND OFF"""
        if self.fault_picking_switch is True:
            self.fault_picking_switch = False
        elif self.fault_picking_switch is False:
            self.fault_picking_switch = True

    def move(self, event):
        """WHEN THE MOUSE IS MOVED"""

        if self.index_node is None and self.selected_node is None:
            # NO NODE WAS FOUND NEAR THE CLICK
            return
        if event.inaxes is None:
            # CLICK WAS OUTSIDE THE MODEL FRAME
            return
        if event.button != 1:
            return
        if self.pinch is True:
            # PINCH MODE IS ON
            return
        if self.pan_on is True:
            # PAN MODE IS ON
            return
        if self.zoom_on is True:
            # ZOOM MODE IS ON
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
            if self.boundary_lock_list[self.i] == 0:
                x = event.xdata  # GET X OF NEW POINT
                y = event.ydata  # GET Y OF NEW POINT

                # GET CURRENT X AND Y ARRAYS
                xt = np.array(self.plotx)
                yt = np.array(self.ploty)

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
            elif self.boundary_lock_list[self.i] == 1:
                x = event.xdata  # GET X OF NEW POINT
                y = event.ydata  # GET Y OF NEW POINT
                xt = np.array(self.plotx)
                yt = np.array(self.ploty)
                if y <= 0:
                    xt[self.index_node] = x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = 0.001  # REPLACE OLD Y WITH NEW Y
                else:
                    xt[self.index_node] = x  # REPLACE OLD X WITH NEW X
                    yt[self.index_node] = y  # REPLACE OLD Y WITH NEW Y

            # DEAL WITH PINCHED MODE
            if self.pinch_switch != 0:
                for k in range(0, len(self.index_arg2_list)):
                    if self.index_arg2_list[k] is not None:
                        # GET THE NODE LIST OF THE NEXT LAYER
                        next_x_list, next_y_list = self.plotx_list[k], self.ploty_list[k]

                        # REPLACE THE PINCHED NODE WITH THE NEW NODE
                        next_x_list[self.index_arg2_list[k]] = x
                        next_y_list[self.index_arg2_list[k]] = y
                        self.plotx_list[k] = next_x_list
                        self.ploty_list[k] = next_y_list  # OVERWRITE THE NODE LIST WITH UPDATED LIST

            self.plotx = xt
            self.ploty = yt
            self.polyline.set_data(self.plotx, self.ploty)

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
        """WHEN MOUSE BUTTON IS RELEASED"""

        if event.inaxes is None:
            # CLICK WAS OUTSIDE THE MODEL FRAME
            return
        if event.button != 1:
            return

        # UPDATE LAYER ATTRIBUTES TABLE WITH CURRENT COORDINATES
        # self.x_input.SetValue(event.xdata)
        # self.y_input.SetValue(event.ydata)

        if self.capture is True:
            # GMG IS IN COORDINATE CAPTURE MODE SO ADD THE CURRENT COORDINATES TO THE TABLE
            self.capture_window.table.Append((event.xdata, event.ydata))

        # UPDATE LAYERS
        self.update_layer_data()

        # RUN MODELLING ALGORITHMS
        if self.fault_picking_switch is False:
            self.run_algorithms()

    def key_press(self, event):
        """DEFINE KEY PRESS LINKS"""

        if self.fault_picking_switch is True:
            # GMG IS IN FAULT MODE SO USE FAULT MODE KEY FUNCTIONS
            self.fault_mode_key_press(event)
            return

        'f = ACTIVATE FAULT PICKING MODE'
        if event.key == 'f':
            # TURN ON/OFF FAULT PICKING MODE
            if self.fault_picking_switch is True:
                self.fault_picking_switch = False
            else:
                self.fault_picking_switch = True

        'i = INSERT NEW NODE AT MOUSE POSITION'
        if event.key == 'i':
            if event.inaxes is None:
                return

            # GET CURRENT LAYER XY
            xt = np.array(self.plotx)
            yt = np.array(self.ploty)

            # INSERT NEW NODES INTO LAYER X AND Y LISTS
            self.plotx = np.insert(xt, [self.index_arg + 1], event.xdata)
            self.ploty = np.insert(yt, [self.index_arg + 1], event.ydata)

            self.nextpoly.append([event.xdata, event.ydata])
            self.polyline.set_data(self.plotx, self.ploty)

            # UPDATE LAYER DATA AND PLOT
            self.update_layer_data()
            self.run_algorithms()
            self.draw()

        'd = DELETE NODE AT MOUSE POSITION'
        if event.key == 'd':
            xt = np.array(self.plotx)
            yt = np.array(self.ploty)

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
                self.plotx = [tup for i, tup in enumerate(self.plotx) if i != self.index_arg]  # DELETE X
                self.ploty = [tup for i, tup in enumerate(self.ploty) if i != self.index_arg]  # DELETE Y
                self.polyline.set_data(self.plotx, self.ploty)

            # NOW CHECK FOR PINCHED NODES
            index_arg2 = None
            self.pinch_switch = 0
            self.index_arg2_list = [None] * (self.layer_count + 1)  # CREATE LIST OF NONES = LENGTH AS NUMB OF LAYERS
            for x in range(0, self.layer_count + 1):  # LOOP THROUGH ALL LAYERS TO CHECK FOR PINCHED NODES
                if x == self.i:
                    pass
                x_node_list = self.plotx_list[x]
                y_node_list = self.ploty_list[x]
                for i in range(0, len(x_node_list)):
                    # NOW CHECK X AND Y VALUES ARE EQUAL
                    if x_node_list[i] == xt[self.index_arg] and y_node_list[i] == yt[self.index_arg]:
                        # IF ONE OF THE NODES FORM LIST IS EQUAL TO A NODE FROM THE OTHER LAYER THEN RETURN THE INDEX
                        self.index_arg2_list[x] = i
                        self.pinch_switch = 1

            # REMOVE PINCHED NODES
            if self.pinch_switch != 0:
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

        'q = BEAT THE ZOOM BUG'
        if event.key == 'q':
            self.nodes = True

        'n = CREATE NEW LAYER AT MOUSE POINT'
        if event.key == 'n':
            self.new_layer(event)

        'b = LOCK OR UNLOCK LAYER BOUNDARY LOCKED MODE'
        if event.key == 'b':
            if self.boundary_lock:
                'unlock layer'
                self.boundary_lock = False
                self.boundary_lock_list[self.i] = 1
            elif not self.boundary_lock:
                'lock layer'
                self.boundary_lock = True
                self.boundary_lock_list[self.i] = 0

        'l = LOCK OR UNLOCK LAYER/POLYGON MODE'
        if event.key == 'l':
            if self.layer_lock:
                'unlock layer'
                self.layer_lock = False
                self.layer_lock_list[self.i] = 1
            elif not self.layer_lock:
                'lock layer'
                self.layer_lock = True
                self.layer_lock_list[self.i] = 0

        '< = MOVE TO NEXT LAYER'
        if event.key == '.':
            if self.i == self.layer_count:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.i = 0
                self.density_input.SetValue(0.001 * self.densities[self.i])
                self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
                self.susceptibility_input.SetValue(self.susceptibilities[self.i])
                self.angle_a_input.SetValue(self.angle_a[self.i])
                self.plotx = self.plotx_list[self.i]
                self.ploty = self.ploty_list[self.i]
                self.x_input.SetValue(self.plotx[0])
                self.y_input.SetValue(self.ploty[0])
                self.current_node.set_offsets([self.plotx[0], self.ploty[0]])
                self.update_layer_data()
                self.run_algorithms()
            else:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.i = self.i + 1
                self.nextdens = self.densities[self.i]
                self.density_input.SetValue(0.001 * self.densities[self.i])
                self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
                self.susceptibility_input.SetValue(self.susceptibilities[self.i])
                self.angle_a_input.SetValue(self.angle_a[self.i])
                self.plotx = self.plotx_list[self.i]
                self.ploty = self.ploty_list[self.i]
                self.x_input.SetValue(self.plotx[0])
                self.y_input.SetValue(self.ploty[0])
                self.current_node.set_offsets([self.plotx[0], self.ploty[0]])
                self.update_layer_data()
                self.run_algorithms()

        '< = MOVE TO NEXT LAYER'
        if event.key == ',':
            if self.i == 0:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.i = self.layer_count
                self.density_input.SetValue(0.001 * self.densities[self.i])
                self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
                self.susceptibility_input.SetValue(self.susceptibilities[self.i])
                self.angle_a_input.SetValue(self.angle_a[self.i])
                self.plotx = self.plotx_list[self.i]
                self.ploty = self.ploty_list[self.i]
                self.x_input.SetValue(self.plotx[0])
                self.y_input.SetValue(self.ploty[0])
                self.current_node.set_offsets([self.plotx[0], self.ploty[0]])
                self.update_layer_data()
                self.run_algorithms()
            else:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.i = self.i - 1
                self.nextdens = self.densities[self.i]
                self.density_input.SetValue(0.001 * self.densities[self.i])
                self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
                self.susceptibility_input.SetValue(self.susceptibilities[self.i])
                self.angle_a_input.SetValue(self.angle_a[self.i])
                self.plotx = self.plotx_list[self.i]
                self.ploty = self.ploty_list[self.i]
                self.x_input.SetValue(self.plotx[0])
                self.y_input.SetValue(self.ploty[0])
                self.current_node.set_offsets([self.plotx[0], self.ploty[0]])
                self.update_layer_data()
                self.run_algorithms()
            self.draw()

        'z = ZOOM IN MODE'
        if event.key == 'z':
            self.zoom(event)

        'ctrl+z = ZOOM OUT'
        if event.key == 'ctrl+z':
            self.zoom_out(event)

        'shift = PAN MODE'
        if event.key == 'ctrl+p':
            self.pan(event)

        'a = FULL EXTENT VIEW'
        if event.key == 'a':
            self.full_extent(event)

        'p = TURN ON PINCH NODE MODE'
        if event.key == 'p':
            if not self.pinch:
                self.pinch = True
            else:
                self.pinch = False
                self.pinch_node_list = [[], []]
                self.pinch_count = 0

        'ctrl+i = INCREASE LAYER TRANSPARENCY'
        if event.key == 'ctrl+i':
            self.transparency_increase(event)

        'ctrl+d = INCREASE ASPECT TRANSPARENCY'
        if event.key == 'ctrl+d':
            self.transparency_decrease(event)

        'up arrow = INCREASE ASPECT RATIO'
        if event.key == 'up':
            self.aspect_increase(event)

        'ctrl+up = INCREASE ASPECT RATIO X2'
        if event.key == 'ctrl+up':
            self.aspect_increase2(event)

        'down arrow = DECREASE ASPECT RATIO'
        if event.key == 'down':
            self.aspect_decrease(event)

        'ctrl+down = DECREASE ASPECT RATIO'
        if event.key == 'ctrl+down':
            self.aspect_decrease2(event)

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

            if self.index_arg == 0 or self.index_arg == (len(self.fault_x_coords_list[self.current_fault_index])-1):
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
                self.faultline, = self.mcanvas.plot([-50000, 0], [-49999, 0], marker='s', color='m', linewidth=0.75,
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
            self.faults[self.current_fault_index] = self.mcanvas.plot(
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
            for i in range(0, self.layer_count + 1):
                out = csv.writer(f, delimiter=' ')
                f.write('>\n')
                data = [self.plotx_list[i], self.ploty_list[i]]
                out.writerows(zip(*data))
                layer_write = zip(self.plotx_list[i], self.ploty_list[i])
                # WRITE INDIVIDUAL LAYER
                np.savetxt(output_dir + '/' + self.loaded_tree_items[i] + '.xy', layer_write, delimiter=' ',
                           fmt='%f %f')

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
            for i in range(0, self.layer_count + 1):
                f.write('B  {0}\n'.format(i + 1))
                data = zip(self.plotx_list[i], self.ploty_list[i], np.ones(len(self.ploty_list[i])))
                # print data
                np.savetxt(f, data, delimiter=' ', fmt='%6.02f %3.02f %1d')

            # VELOCITY NODES
            for i in range(0, self.layer_count):
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
        pinch_box = PinchDialog(self, -1, 'Pinch Out Layer:', self.plotx_list, self.ploty_list, self.i)
        answer = pinch_box.ShowModal()
        self.plotx = pinch_box.pinched_x
        self.ploty = pinch_box.pinched_y
        self.update_layer_data()
        self.draw()

    def depinch_layer(self, event):
        depinch_box = DepinchDialog(self, -1, 'Depinch layer', self.plotx_list, self.ploty_list, self.i,
                                    self.layer_count)
        answer = depinch_box.ShowModal()
        self.plotx = depinch_box.depinched_x
        self.ploty = depinch_box.depinched_y
        self.update_layer_data()
        self.draw()

    def bulk_shift(self, event):
        bulk_shift_box = BulkShiftDialog(self, -1, 'Layer bulk shift', self.plotx_list, self.ploty_list, self.i)
        answer = bulk_shift_box.ShowModal()
        self.plotx = bulk_shift_box.new_x
        self.ploty = bulk_shift_box.new_y
        self.update_layer_data()
        self.draw()

    def observed_filter(self, event):
        """FILTER OBSERVED ANOMALY USING MEDIAN FILTER - CALLS class MedianFilterDialog"""

        # RUN FILTER
        median_filter_box = MedianFilterDialog(self, -1, 'median filter', self.obs_grav_name_list,
                                               self.obs_grav_list_save, self.obs_mag_name_list, self.obs_mag_list_save)
        answer = median_filter_box.ShowModal()

        # GET FILTERED OUTPUT
        filtered_data = median_filter_box.filtered_output
        filtered_name = median_filter_box.output_name
        filtered_color = median_filter_box.output_color
        filter_type = median_filter_box.filter_type

        # LOAD FILTERED DATA
        if filter_type == "gravity":
            self.obs_grav_name_list.append(filtered_name)
            self.obs_grav_list_save.append([])
            self.obs_grav_list_save[self.obs_grav_count] = filtered_data
            self.obs_grav_list.append([])
            self.obs_grav_list[self.obs_grav_count] = self.dcanvas.scatter(filtered_data[:, 0], filtered_data[:, 1],
                                                                           marker='o', color=filtered_color, s=5,
                                                                           gid=self.obs_grav_count)

            #  APPEND NEW DATA MENU TO 'Topo data MENU'
            self.grav_submenu = wx.Menu()
            self.m_obs_g_submenu.Append(self.obs_grav_count,filtered_name, self.grav_submenu)

            # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
            self.grav_submenu.Append(11000 + self.obs_grav_count, 'delete observed data')

            # BIND TO DEL TOPO FUNC
            self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000 + self.obs_grav_count)

            # INCREMENT GRAV COUNTER
            self.obs_grav_count += 1

            # UPDATE GMG GUI
            self.update_layer_data()
            self.draw()

        if filter_type == "magnetics":
            self.obs_mag_name_list.append(filtered_name)
            self.obs_mag_list_save.append([])
            self.obs_mag_list_save[self.obs_mag_count] = filtered_data
            self.obs_mag_list.append([])
            self.obs_mag_list[self.obs_mag_count] = self.dcanvas.scatter(filtered_data[:, 0], filtered_data[:, 1],
                                                                         marker='o', color=filtered_color, s=5,
                                                                         gid=self.obs_mag_count)

        #  APPEND NEW DATA MENU TO 'Mag data MENU'
        self.obs_mag_submenu = wx.Menu()
        self.m_obs_mag_submenu.Append(self.obs_mag_count, filtered_name, self.obs_mag_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.obs_mag_submenu.Append(12000+self.obs_mag_count, 'delete observed data')

        # BIND TO DEL MAG FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12000+self.obs_mag_count)

        self.obs_mag_count += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()


    def new_layer(self, event):
        new_layer_dialogbox = NewLayerDialog(self, -1, 'Create New Layer')
        answer = new_layer_dialogbox.ShowModal()
        
        if new_layer_dialogbox.fixed:
            # DETERMINE THE LAST "FIXED LAYER" X AND Y VALUES; SET THIS LAYER AS "k"
            if self.layer_count > 0:
                for i in range(0, self.layer_count):
                    if self.layer_lock_list[i] == 0:
                        k = i
                    else:
                        continue
            else:
                k = 0

            # NOW APPEND NODES FOR BOUNDARY CONDITIONS (CONTINUOUS SLAB)
            layer_above_x = np.array(self.plotx_list[k])
            layer_above_y = np.array(self.ploty_list[k])

            # INCREMENT THE CURRENT LAYER INDEX VALUE (self.i)
            self.i = self.layer_count
            self.i += 1

            # INCREMENT THE TOTAL LAYER COUNT
            self.layer_count += 1

            # ADD A NEW BLANK LAYER TO THE PLOT LISTS
            self.polygon_fills.append([])
            self.layer_lines.append([])
            self.plotx_list.append([])
            self.ploty_list.append([])
            self.poly_fills.append([])
            self.densities.append(0.)
            self.reference_densities.append(0.)
            self.susceptibilities.append(0.)
            self.angle_a.append(0.)
            self.angle_b.append(0.)
            self.layer_lock_list.append(0)
            self.boundary_lock_list.append(0)
            self.layer_lock_status.append('locked')
            self.boundary_lock_status.append('locked')
            self.layer_colors.append('black')
            self.tree_items.append('layer %s' % (int(self.i)))
            self.item = 'layer %s' % (int(self.i ))
            self.layers_calculation_switch.append(1)
            self.add_new_tree_nodes(self.root, self.item, self.i)

            # SET NEW LAYER NODES
            self.new_layer_thickness = new_layer_dialogbox.new_thickness
            self.plotx = layer_above_x
            self.ploty = layer_above_y + self.new_layer_thickness
            self.nextpoly = zip(self.plotx, self.ploty)

            # CREATE LAYER LINE
            self.layer_lines[self.i] = self.mcanvas.plot(self.plotx_list[self.i], self.ploty_list[self.i],
                                                         color='blue', linewidth=1.0, alpha=1.0)
            # CREATE LAYER POLYGON FILL
            self.plotx_polygon = np.array(self.plotx)
            self.ploty_polygon = np.array(self.ploty)
            self.polygon_fills[self.i] = self.mcanvas.fill(self.plotx_polygon, self.ploty_polygon, color='blue',
                                                           alpha=self.layer_transparency, closed=True, linewidth=None,
                                                           ec=None)
            # UPDATE LAYER DATA
            self.current_node.set_offsets([self.plotx[0], self.ploty[0]])
            self.nextdens = self.densities[self.i]
            self.density_input.SetValue(0.001 * self.densities[self.i])
            self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
            self.susceptibility_input.SetValue(self.susceptibilities[self.i])
            self.angle_a_input.SetValue(self.angle_a[self.i])
            self.update_layer_data()
            self.draw()

        elif not new_layer_dialogbox.fixed:
            # GET NEW NODE VALUES
            self.new_x1, self.new_y1 = new_layer_dialogbox.x1, new_layer_dialogbox.y1
            self.new_x2, self.new_y2 = new_layer_dialogbox.x2, new_layer_dialogbox.y2
            self.new_x3, self.new_y3 = new_layer_dialogbox.x3, new_layer_dialogbox.y3
            self.new_x4, self.new_y4 = new_layer_dialogbox.x4, new_layer_dialogbox.y4

            # INCREMENT THE LAYER COUNT'
            self.i = self.layer_count
            self.i = self.i + 1

            # INCREMENT THE TOTAL LAYER COUNT
            self.layer_count += 1

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
            self.tree_items.append('layer %s' % (int(self.i)))
            self.item = 'layer %s' % (int(self.i))
            self.add_new_tree_nodes(self.root, self.item, self.i)
            self.layers_calculation_switch.append(1)
            self.plotx = [self.new_x1, self.new_x2, self.new_x3, self.new_x4]
            self.ploty = [self.new_y1, self.new_y2, self.new_y3, self.new_y4]
            self.nextpoly = zip(self.plotx, self.ploty)

            # CREATE LAYER LINE
            self.layer_lines[self.i] = self.mcanvas.plot(self.plotx_list[self.i], self.ploty_list[self.i],
                                                         color='blue', linewidth=1.0, alpha=1.0)
            # CREATE LAYER POLYGON FILL
            self.plotx_polygon = np.array(self.plotx)
            self.ploty_polygon = np.array(self.ploty)
            self.polygon_fills[self.i] = self.mcanvas.fill(self.plotx_polygon, self.ploty_polygon, color='blue',
                                                           alpha=self.layer_transparency, closed=True, linewidth=None,
                                                           ec=None)

            # UPDATE LAYER DATA
            self.current_node.set_offsets([self.plotx[0], self.ploty[0]])
            self.nextdens = self.densities[self.i]
            self.density_input.SetValue(0.001 * self.densities[self.i])
            self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
            self.susceptibility_input.SetValue(self.susceptibilities[self.i])
            self.angle_a_input.SetValue(self.angle_a[self.i])
            self.update_layer_data()
            self.run_algorithms()
            self.draw()
        else:
            # USER CHANGED THEIR MIND - NO NEW LAYER ADDED'
            pass

    def load_layer(self, event):
        open_file_dialog = wx.FileDialog(self, "Open Layer", "", "", "Layer XY files (*.txt)|*.txt",
                                         wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if open_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THERE MIND

        file_in = open_file_dialog.GetPath()
        new_layer = np.genfromtxt(file_in, autostrip=True, delimiter=' ', dtype=float)

        # INCREMENT THE LAYER COUNT
        self.i = self.layer_count
        self.i += 1

        # INCREMENT THE TOTAL LAYER COUNT
        self.layer_count += 1

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
        self.tree_items.append('layer %s' % (int(self.i + 1)))
        self.item = 'layer %s' % (int(self.i + 1))
        self.add_new_tree_nodes(self.root, self.item, self.i)
        self.layers_calculation_switch.append(0)
        self.plotx = new_layer[:, 0]
        self.ploty = new_layer[:, 1]
        self.nextpoly = zip(self.plotx, self.ploty)

        'CREATE LAYER LINE'
        self.layer_lines[self.i] = self.mcanvas.plot(self.plotx_list[self.i], self.ploty_list[self.i],
                                                     color='blue', linewidth=1.0, alpha=1.0)
        'CREATE LAYER POLYGON FILL'
        self.plotx_polygon = np.array(self.plotx)
        self.ploty_polygon = np.array(self.ploty)
        self.polygon_fills[self.i] = self.mcanvas.fill(self.plotx_polygon, self.ploty_polygon, color='blue',
                                                       alpha=self.layer_transparency, closed=True, linewidth=None,
                                                       ec=None)
        'Update layer data'
        self.nextdens = self.densities[self.i]
        self.density_input.SetValue(0.001 * self.densities[self.i])
        self.ref_density_input.SetValue(0.001 * self.reference_densities[self.i])
        self.susceptibility_input.SetValue(self.susceptibilities[self.i])
        self.angle_a_input.SetValue(self.angle_a[self.i])
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def delete_layer(self, event):
        """Delete LAYER DATA"""

        # REMOVE POLYGON
        self.polygon_fills[self.i][0].remove()
        del self.polygon_fills[self.i]
        # REMOVE POLYGON
        self.layer_lines[self.i][0].remove()
        del self.layer_lines[self.i]
        # SET CURRENT AS THE NEXT LAYER'
        if self.i == self.layer_count:
            self.polyline.set_xdata(self.plotx_list[0])
            self.polyline.set_ydata(self.ploty_list[0])
            self.polyline.set_color(self.layer_colors[0])
            self.plotx = self.plotx_list[0]
            self.ploty = self.ploty_list[0]
        else:
            self.polyline.set_xdata(self.plotx_list[self.i + 1])
            self.polyline.set_ydata(self.ploty_list[self.i + 1])
            self.polyline.set_color(self.layer_colors[self.i + 1])
            self.plotx = self.plotx_list[self.i + 1]
            self.ploty = self.ploty_list[self.i + 1]
        self.draw()

        # REMOVE META DATA'
        del self.plotx_list[self.i]
        del self.ploty_list[self.i]
        del self.densities[self.i]
        del self.reference_densities[self.i]
        del self.susceptibilities[self.i]
        del self.angle_a[self.i]
        del self.angle_b[self.i]
        del self.layer_lock_list[self.i]
        del self.boundary_lock_list[self.i]
        del self.layer_lock_status[self.i]
        del self.boundary_lock_status[self.i]
        del self.layer_colors[self.i]
        del self.layers_calculation_switch[self.i]

        #  REMOVE TREE ITEMS
        del self.tree_items[self.i]
        layers = self.tree.GetRootItem().GetChildren()
        self.tree.Delete(layers[self.i])
        # RESET TREE ITEM ID'S
        layers = self.tree.GetRootItem().GetChildren()
        for i in range(len(layers)):
            self.tree.SetPyData(layers[i], i)

        # UPDATE COUNTERS
        self.layer_count -= 1

        # UPDATE
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    # LAYER AND MODEL ATTRIBUTE CONTROLS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_background_density(self, event):
        grav_box = SetBackgroundDensityDialog(self, -1, 'Set background density')
        answer = grav_box.ShowModal()
        self.background_density_upper = float(grav_box.background_density_upper)

        for i in range(0, self.layer_count):
            self.reference_densities[i] = float(self.background_density_upper) * 1000.
        # self.background_density_upper = float((grav_box.background_density_lower))
        # self.background_density_upper = float((grav_box.background_density_lid))
        self.absolute_densities = True
        self.draw()

    def set_density(self, value):
        if self.density_input.GetValue() == 0:
            self.densities[self.i] = 0
        else:
            self.densities[self.i] = float(self.density_input.GetValue() * 1000.)

    def set_reference_density(self, value):
        if self.ref_density_input.GetValue() == 0:
            self.reference_densities[self.i] = 0
        else:
            self.reference_densities[self.i] = float(self.ref_density_input.GetValue() * 1000.)

    def set_susceptibility(self, value):
        self.susceptibilities[self.i] = float(self.susceptibility_input.GetValue())

    def set_angle_a(self, value):
        self.angle_a[self.i] = float(self.angle_a_input.GetValue())

    def set_angle_b(self, value):
        self.angle_b[self.i] = float(self.angle_b_input.GetValue())

    def set_text_size(self, value):
        """GET NEW TEXT SIZE"""

        self.textsize = float(self.text_size_input.GetValue())

        # WELL DATA
        # LOOP THROUGH ALL WELL NAMES
        if len(self.well_name_text) > 0:
            for i in range(len(self.well_name_text)):
                if self.well_name_text[i] != "None" and self.well_name_text[i] != []:
                    self.well_name_text[i].set_size(self.textsize)
        # LOOP THROUGH ALL WELL HORIZON LABELS
        if len(self.well_labels) > 0:
            for i in range(len(self.well_labels)):
                labels = self.well_labels[i]
                for l in range(2, len(labels)):
                    if labels[l] != "None" and labels[l] != []:
                        labels[l].set_size(self.textsize)

        # LOOP THROUGH OUTCROP DATA LABELS
        if self.outcrop_data_count > 0:
            for i in range(self.outcrop_data_count):
                outcrop_text = self.outcrop_text_list[i]
                for k in range(len(outcrop_text)):
                    if outcrop_text[k] is not None:
                        outcrop_text[k].set_size(self.textsize)
                    else:
                        pass

        # REDRAW ANNOTATIONS WITH NEW TEXT SIZE
        self.draw()

    def set_obs_rms(self, value):
        set_values = SetObsRmsDialog(self, -1, 'Set RMS Input', self.obs_grav_name_list, self.obs_grav_list_save,
                                     self.obs_mag_name_list, self.obs_mag_list_save)
        answer = set_values.ShowModal()
        self.obs_gravity_data_for_rms = set_values.obs_gravity_data_for_rms
        self.obs_mag_data_for_rms = set_values.obs_mag_data_for_rms

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
        self.draw()

    # LIVE GRAPHICS UPDATES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def load_layer_data(self):
        """LOAD LAYER DATA INTO GMG"""

        # UPDATE MODEL CANVAS (mcanvas) LIMITS
        xmin, xmax = self.mcanvas.get_xlim()
        self.dcanvas.set_xlim(xmin, xmax)
        self.ntcanvas.set_xlim(xmin, xmax)

        # SET LISTS WITH UPDATED LAYER DATA
        self.plotx_list[self.i] = self.plotx
        self.ploty_list[self.i] = self.ploty

        # RESET LISTS (UPDATED BY THE CALCULATE_GRAVITY FUNC)
        self.polygons = []
        self.mag_polygons = []

        # REMOVE EXISTING DATA
        if len(self.mcanvas.lines) >= 1:
            while len(self.mcanvas.lines) > 0:
                self.mcanvas.lines[0].remove()
            while len(self.mcanvas.patches) > 0:
                self.mcanvas.patches[0].remove()
            while len(self.mcanvas.texts) > 0:
                self.mcanvas.texts[0].remove()

        # UPDATE DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # FIRST CREATE THE POLYLINE DATA (THE BOTTOM LINE OF THE LAYER POLYGON - DONE FIRST SO THE WHOLE POLYGON
        # ISN'T PASSED TO SELF.POLYPLOTS)
        for i in range(0, self.layer_count + 1):
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
            self.polygon_fills[i] = self.mcanvas.fill(self.plotx_polygon, self.ploty_polygon, color=next_color,
                                                      alpha=self.layer_transparency, closed=True, linewidth=None,
                                                      ec=None)

            self.polygons.append(zip(self.plotx_polygon, self.ploty_polygon))

        # MODEL LAYERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CREATE LAYER LINES
        for i in range(0, self.layer_count + 1):
            self.layer_lines[i] = self.mcanvas.plot(self.plotx_list[i], self.ploty_list[i], color=self.layer_colors[i],
                                                    linewidth=1.0, alpha=0.5)
        # CURRENT LAYER LINE
        self.polyline, = self.mcanvas.plot(self.plotx, self.ploty, marker='o', color=self.layer_colors[self.i],
                                           linewidth=1.0, alpha=0.5)

        # # FAULTS LAYERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # # CREATE FAULT LINES
        # for i in range(0, self.fault_counter + 1):
        #     self.faults[i] = self.mcanvas.plot(self.fault_plotx_list[i], self.fault_ploty_list[i],
        #                                             color=self.fault_colors[i], linewidth=1.0, alpha=1.0)
        # # CURRENT FAULT LINE
        # self.fault_polyline, = self.mcanvas.plot(self.fault_plotx, self.fault_ploty, marker='o',
        #                                          color=self.fault_colors[self.current_fault_index], linewidth=1.0, alpha=0.5)

        # UPDATE GMG
        self.mcanvas.set_aspect(self.model_aspect)
        self.display_info()
        self.draw()

    def update_layer_data(self):
        """UPDATE PROGRAM GRAPHICS AFTER A CHANGE IS MADE (REDRAWS EVERYTHING)"""

        # UPDATE MODEL CANVAS (mcanvas) LIMITS
        xmin, xmax = self.mcanvas.get_xlim()
        if self.t_canvas:
            self.dcanvas.set_xlim(xmin, xmax)
        if self.d_canvas:
            self.dcanvas.set_xlim(xmin, xmax)
        if self.nt_canvas:
            self.ntcanvas.set_xlim(xmin, xmax)

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
            self.plotx_list[self.i] = self.plotx
            self.ploty_list[self.i] = self.ploty

            # RESET LISTS (USED AS INPUT FOR THE GRAV/MAG ALGORITHMS)
            self.polygons = []
            self.mag_polygons = []

            # CREATE UPDATED POLYGON XYs -------------------------------------------------------------------------------
            # FIRST CREATE THE POLYLINE DATA (THE BOTTOM LINE OF THE LAYER POLYGON - (THIS DONE FIRST SO THE WHOLE
            # POLYGON ISN'T PASSED TO SELF.POLYPLOTS)
            for i in range(0, self.layer_count + 1):
                # CREATE THE LAYER POLYGONS TO PASS TO SELF.POLYGONS AND ONTO THE BOTT ALGORITHM
                # FIRST SET UP XY DATA; IF THE LAYER IS BELOW LAYER 1 THEN ATTACH THE ABOVE LAYER TO COMPLETE POLYGON;
                # ELSE USE TOP LAYER CHECK FOR LAYER MODE AND FIND LAST LAYER TO MAKE POLYGON
                if i >= 1 and self.layer_lock_list[i] == 0:
                    for layers in range(i, 0, -1):
                        if self.layer_lock_list[layers - 1] == 0:
                            self.last_layer = layers - 1

                            # NOW APPEND NODES FOR BOUNDARY CONDITIONS (CONTINUOUS SLAB)
                            plotx, ploty = np.array(self.plotx_list[i]), np.array(self.ploty_list[i])

                            # SET PADDING NODES TO DEPTH EQUAL TO MODEL LIMIT NODES TO CREATE SLAB
                            ploty[0], ploty[-1] = ploty[1], ploty[-2]
                            self.plotx_list[i], self.ploty_list[i] = plotx, ploty

                            # ADD NODES FROM ABOVE LAYER TO COMPETE POLYGON
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

                self.polygons.append(zip(self.plotx_polygon, self.ploty_polygon))
            # -----------------------------------------------------------------------------------------------------------

            # UPDATE LAYER POLYGONS AND LINES
            for i in range(0, self.layer_count + 1):

                # SET POLYGON FILL COLOR
                if self.densities[i] != 0 and self.absolute_densities is True:
                    # DETERMINE DENSITY CONTRAST FROM (DENSITY - REF DENSITY)
                    next_color = self.colormap.to_rgba(0.001 * self.densities[i] - 0.001 * self.reference_densities[i])
                elif self.densities[i] != 0:
                    # NO REF DENSITY, SO JUST USE DENSITY VALUE
                    next_color = self.colormap.to_rgba(0.001 * self.densities[i])
                else:
                    # NO DENSITY HAS BEEN SET SO LEAVE BLANK
                    next_color = self.colormap.to_rgba(0.)

                # UPDATE POLYGON XY AND COLOR FILL
                self.polygon_fills[i][0].set_xy(self.polygons[i])
                self.polygon_fills[i][0].set_color(next_color)

                # LINES
                self.layer_lines[i][0].set_xdata(self.plotx_list[i])
                self.layer_lines[i][0].set_ydata(self.ploty_list[i])

            # SET LINE FOR CURRENT LAYER BEING EDITED
            self.polyline.set_xdata(self.plotx)
            self.polyline.set_ydata(self.ploty)
            self.polyline.set_color(self.layer_colors[self.i])

        # DRAW CANVAS FEATURES
        self.mcanvas.set_aspect(self.model_aspect)
        self.dcanvas_aspect = ((self.dcanvas.get_xlim()[1] - self.dcanvas.get_xlim()[0]) /
                               (self.dcanvas.get_ylim()[1] - self.dcanvas.get_ylim()[0]))
        # UPDATE INFO
        self.display_info()

        # CONTENT HAS NOT BEEN SAVED SINCE LAST MODIFICATION
        self.model_saved = False

        # UPDATE GMG GRAPHICS
        self.draw()

    def run_algorithms(self):
        """RUN POTENTIAL FIELD CALCULATION ALGORITHMS"""
        # REVERSE ORDER OF POLYGONS SO NODES ARE INPUT TO ALGORITHMS IN CLOCKWISE DIRECTION
        # (AS LAYER NODES ARE STORED IN LEFT TO RIGHT ORDER IN THE NUMPY ARRAYS
        if self.calc_grav_switch is True or self.calc_mag_switch is True:
            if len(self.polygons) > 1:
                for i in range(0, len(self.polygons)):
                    self.polygons[i] = self.polygons[i][::-1]
            else:
                self.polygons[0] = self.polygons[0][::-1]

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # TOPOGRAPHY - :FUTURE: PREDICTED TOPOGRAPHY FROM ISOSTATIC FUNC
        self.pred_topo = np.zeros_like(self.xp)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # GRAVITY
        # ARRAY OF DENSITY CONTRASTS TO BE PASSED TO BOTT ALGORITHM
        self.density_contrasts = np.zeros(len(self.densities))
        for i in range(len(self.densities)):
            self.density_contrasts[i] = (self.densities[i] - self.reference_densities[i])

        # ZIP POLYGONS WITH DENSITY CONTRASTS AND PASS TO BOTT
        if self.polygons and self.calc_grav_switch is True:
            # SELECT ONLY THOSE LAYERS THAT ARE CHECKED
            polygons_to_use, densities_to_use = [], []
            for layer in range(len(self.densities)):
                if self.layers_calculation_switch[layer] == 1:

                    polygons_to_use.append(self.polygons[layer])
                    densities_to_use.append(self.density_contrasts[layer])

            # PASS TO BOTT ALGORITHM
            polys = []
            for p, d in zip(polygons_to_use, densities_to_use):
                polys.append(Polygon(1000 * np.array(p), {'density': d}))
            self.predgz = bott.Gz(self.xp, self.zp, polys)
        else:
            self.predgz = np.zeros_like(self.xp)
        self.predplot.set_data(self.xp * 0.001, self.predgz)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # MAGNETICS
        # ZIP POLYGONS WITH SUSCEPT/STRIKES AND PASS TO TALWANI AND HEIRTZLER ALGORITHM
        if self.polygons and self.earth_field != 0.0 and self.calc_mag_switch is True:
            # SELECT ONLY THOSE LAYERS THAT ARE CHECKED
            polygons_to_use, susceptibilities_to_use = [], []
            for layer in range(0, len(self.polygons)):
                if self.layers_calculation_switch[layer] == 1:
                    polygons_to_use.append(self.polygons[layer])
                    susceptibilities_to_use.append(self.susceptibilities[layer])

            # PASS TO TALWANI AND HEIRTZLER ALGORITHM
            self.polys = []
            for p, s, in zip(polygons_to_use, susceptibilities_to_use):
                self.polys.append(Polygon(1000. * np.array(p), {'susceptibility': s}))
            self.prednt = talwani_and_heirtzler.ntz(self.xp, self.zp, self.polys, self.model_azimuth, self.earth_field,
                                                    self.angle_a, self.angle_b, self.mag_observation_elv)
        else:
            self.prednt = np.zeros_like(self.xp)

        # PLOT MAGNETIC PROFILE
        self.prednt_plot.set_data(self.xp * 0.001, self.prednt)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # UPDATE RMS VALUES

        #  RUN RMS CALC CODE
        self.model_rms(self.xp)

        # SET GRAV RMS VALUES
        if self.obs_gravity_data_for_rms != [] and self.grav_obs_switch == True and self.calc_grav_switch == True \
                and self.predgz != []:
            self.grav_rms_plot.set_data(self.grav_residuals[:, 0], self.grav_residuals[:, 1])
        elif self.obs_gravity_data_for_rms is None:
            self.grav_rms_plot.set_data(self.obs_gravity_data_for_rms[:, 0],
                                        np.zeros_like(self.obs_gravity_data_for_rms[:, 0]))
        else:
            pass

        # SET MAG RMS VALUES
        if self.obs_mag_data_for_rms != [] and self.mag_obs_switch == True and self.calc_mag_switch == True \
                and self.prednt != []:
            self.mag_rms_plot.set_data(self.mag_residuals[:, 0], self.mag_residuals[:, 1])
        elif self.obs_mag_data_for_rms is None:
            self.mag_rms_plot.set_data(self.obs_mag_data_for_rms[:, 0],
                                       np.zeros_like(self.obs_mag_data_for_rms[:, 0]))
        else:
            pass
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # SET FRAME X AND Y LIMITS
        self.set_frame_limits()

        # AFTER RUNNING ALGORITHMS, SET MODEL AS UNSAVED
        self.model_saved = False

        # UPDATE GMG GRAPHICS
        self.draw()

    def set_frame_limits(self):
        """SET FRAME X AND Y LIMITS"""

        # SET GRAVITY DISPLAY BOX LIMITS
        if self.grav_obs_switch is True and self.grav_residuals != []:
            ymin = (np.hstack((self.obs_grav[:, 1], self.predgz, self.grav_residuals[:, 1])).min()) - 2.0
            ymax = (np.hstack((self.obs_grav[:, 1], self.predgz, self.grav_residuals[:, 1])).max()) + 2.0
        elif self.grav_obs_switch:
            ymin = (np.append(self.obs_grav[:, 1], self.predgz).min()) - 2.
            ymax = (np.append(self.obs_grav[:, 1], self.predgz).max()) + 2.0
        else:
            ymin = (self.predgz.min()) - 2.
            ymax = (self.predgz.max()) + 2.
        if self.dcanvas is not None:
            self.dcanvas.set_ylim(ymin, ymax)

        # SET MAGNETICS DISPLAY BOX LIMITS
        if self.mag_obs_switch is True and self.mag_residuals != []:
            ymin_mag = (np.hstack((self.obs_mag[:, 1], self.prednt, self.mag_residuals[:, 1])).min()) - 2.
            ymax_mag = (np.hstack((self.obs_mag[:, 1], self.prednt, self.mag_residuals[:, 1])).max()) + 2.

        elif self.mag_obs_switch:
            ymin_mag = (np.append(self.obs_mag[:, 1], self.prednt).min()) - 2.
            ymax_mag = (np.append(self.obs_mag[:, 1], self.prednt).max()) + 2.
        else:
            ymin_mag = (self.prednt.min()) - 2.
            ymax_mag = (self.prednt.max()) + 2.
        if self.ntcanvas is not None:
            self.ntcanvas.set_ylim(ymin_mag, ymax_mag)

        # UPDATE GMG GRAPHICS
        self.draw()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXTERNAL FIGURE CONSTRUCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def plot_model(self, event):
        """CREATE EXTERNAL FIGURE OF MODEL USING INBUILT FIGURE CONSTRUCTION TOOL"""

        # GET PLOTTING PARAMETERS FROM DIALOG BOX
        self.set_values = PlotSettingsDialog(self, -1, 'Set figure parameters', self.model_aspect, self.dcanvas_aspect)
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
        xmin, xmax = self.mcanvas.get_xlim()
        ymin, ymax = self.mcanvas.get_ylim()
        area = np.array([xmin, xmax, ymin, ymax])

        # RUN PLOT MODEL CODE
        fig_plot = plot_model.plot_fig(self.file_path, self.file_type, area, self.xp, self.obs_topo, self.pred_topo,
                                       self.obs_grav, self.predgz, self.obs_mag, self.prednt, self.layer_count,
                                       self.layer_lock_list, self.plotx_list, self.ploty_list, self.densities,
                                       self.absolute_densities, self.reference_densities, self.segy_plot_list,
                                       self.well_list, self.well_name_list, self.t_canvas, self.d_canvas,
                                       self.nt_canvas, self.aspect_ratio, self.use_tight_layout, self.poly_alpha,
                                       self.fs, self.ms, self.lw, self.font_type, self.layer_colors, self.draw_polygons,
                                       self.draw_layers, self.floating_layers, self.draw_colorbar, self.draw_xy_data,
                                       self.xy_size, self.xy_color, self.colorbar_x, self.colorbar_y,
                                       self.colorbar_size_x, self.colorbar_size_y, self.layer_line_width,
                                       self.layer_alpha, self.grav_rms_value, self.mag_rms_value, self.grav_y_min,
                                       self.grav_y_max, self.xy_list_save, self.draw_wells, self.wells, self.well_fs,
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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # DOCUMENTATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FUTURE MODULES (IN PROCESS)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def horizontal_derivative(self, event):
        pass
        #      """CLACULATE HORIZONTAL DERIVATIVE OF PROFILE"""
        #
        #     '# First interpolate the data onto a even x spacing as required by the
        #     '# finite difference derivative approximation'
        #     N = len(in_data)
        #     x_inc = 1.
        #     interpolat_func = ip.interp1d(in_data[:, 0], in_data[:, 1], kind='slinear')  #create interpolation func
        #     x_interp_values = np.arange(np.round(in_data[:, 0].min(), 0),
        #                                    np.round(in_data[:, 0].max(), 0), x_inc)  # xvalues for interpolated array
        #     y_interp_values = np.array(interpolat_func(x_interp_values))             # yvalues for interpolated array
        #     data            = np.column_stack((x_interp_values, y_interp_values))
        #
        #     'NOW CALCULATE THE DERIVATIVE USING THE INTERPOLATED (EVENLY SPACED) INPUT DATA'
        #     for i in range(1, N-1):
        #         deriv = (data[i+1] - data[i-1])/(2*x_inc)
        #
        #     output = np.column_stack((x_interp_values, deriv))
        #     return output

    def set_error(self, value):
        pass
        # self.error = value
        # self.update_layer_data()
        # self.run_algorithms()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''DIALOG BOXES GO HERE'''


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
        'SET FILE PATH IN CLASS PANEL'
        self.chosen_path = self.open_file_dialog.GetPath()
        self.file_path_text.SetValue(str(self.chosen_path))
        self.open_file_dialog.Destroy()

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
            self.parent.open_outcrop_data('load')
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
        self.i = i
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

        current_x = self.plotx_list[self.i]
        current_y = self.ploty_list[self.i]
        above_x = self.plotx_list[self.i - 1]
        above_y = self.ploty_list[self.i - 1]

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

        current_x = self.plotx_list[self.i]
        current_y = self.ploty_list[self.i]
        below_x = self.plotx_list[self.i + 1]
        below_y = self.ploty_list[self.i + 1]

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
        self.i = i
        self.layer_count = layer_count
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

        current_x = self.plotx_list[self.i]
        current_y = self.ploty_list[self.i]

        if self.i != 0:
            above_x = self.plotx_list[self.i - 1]
            above_y = self.ploty_list[self.i - 1]
        if self.i != self.layer_count:
            below_x = self.plotx_list[self.i + 1]
            below_y = self.ploty_list[self.i + 1]

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

    def __init__(self, parent, id, title, obs_grav_name_list, obs_grav_list_save, obs_mag_name_list, obs_mag_list_save):
        wx.Dialog.__init__(self, parent, id, "Apply Median Filter", style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX
                           )
        input_panel = wx.Panel(self, -1)

        # GET OBSERVED DATA FROM gmg
        self.obs_grav_name_list = obs_grav_name_list
        self.obs_grav_list_save = obs_grav_list_save
        self.obs_mag_name_list = obs_mag_name_list
        self.obs_mag_list_save = obs_mag_list_save

        # CHOOSE INPUT OBSERVED DATA
        self.obs_combo_list_text = wx.StaticText(input_panel, -1, "Input Data:")
        self.obs_combo_list = wx.ComboBox(input_panel, id=-1, value="", choices=[])
        for i in self.obs_grav_name_list:
            try:
                self.obs_combo_list.Append(i)
            except:
                pass
        for i in self.obs_mag_name_list:
            try:
                self.obs_combo_list.Append(i)
            except:
                pass
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
        self.obs_to_filter = str(self.obs_combo_list.GetValue())
        self.filter_length = int(self.filter_window_text.GetValue())
        self.output_name = str(self.output_name_text.GetValue())
        self.output_color = str(self.output_color_text.GetValue())

        for i in range(len(self.obs_grav_name_list)):
            # CHECK GRAV
            if self.obs_grav_name_list[i] == self.obs_to_filter:
                self.filter_type = str('gravity')
                self.obs_input = self.obs_grav_list_save[i]
                self.obs_input_name = str(self.obs_grav_name_list[i])
                self.filtered_output = np.zeros(shape=(len(self.obs_input), 2))
                self.filtered_output[:, 0] = self.obs_input[:, 0]
                self.filtered_output[:, 1] = signal.medfilt(self.obs_input[:, 1], self.filter_length)
                self.EndModal(1)

        for i in range(len(self.obs_mag_name_list)):
            # CHECK MAG
            if self.obs_mag_name_list[i] == self.obs_to_filter:
                self.filter_type = str('magnetics')
                self.obs_input = self.obs_mag_list_save[i]
                self.obs_input_name = str(self.obs_mag_name_list[i])
                self.filtered_output = np.zeros(shape=(len(self.obs_input), 2))
                self.filtered_output[:, 0] = self.obs_input[:, 0]
                self.filtered_output[:, 1] = signal.medfilt(self.obs_input[:, 1], self.filter_length)
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
        self.EndModal(1)


class BulkShiftDialog(wx.Dialog):
    """APPLY A BULK SHIFT OF A LAYER IN THE Z DIRECTION"""

    def __init__(self, parent, id, title, plotx_list, ploty_list, i):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)
        self.ploty_list = ploty_list
        self.plotx_list = plotx_list
        self.i = i
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
        current_y = self.ploty_list[self.i]
        current_x = self.plotx_list[self.i]

        # BULKSHIFT Y NODES, SET TO ZERO IF NEW VALUE IS < 0
        self.new_x = [x + self.x_shift_value if x + self.x_shift_value > 0. else 0.01 for x in current_x]
        self.new_y = [y + self.y_shift_value if y + self.y_shift_value > 0. else 0.01 for y in current_y]
        self.EndModal(1)


class SetObsRmsDialog(wx.Dialog):
    """SET THE DATASET TO USE FOR CALCULATING THE RMS MISFIT"""

    def __init__(self, parent, id, title, obs_grav_name_list, obs_grav_list_save, obs_mag_name_list, obs_mag_list_save):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        self.obs_grav_name_list = obs_grav_name_list
        self.obs_grav_list_save = obs_grav_list_save
        self.obs_mag_name_list = obs_mag_name_list
        self.obs_mag_list_save = obs_mag_list_save

        '# CHOOSE INPUT OBSERVED DATA'
        self.obs_combo_list_grav_text = wx.StaticText(input_panel, -1, "Gravity Input:")
        self.obs_combo_list_grav = wx.ComboBox(input_panel, id=-1, value="", choices=[])
        self.obs_combo_list_mag_text = wx.StaticText(input_panel, -1, "Magnetic Input:")
        self.obs_combo_list_mag = wx.ComboBox(input_panel, id=-1, value="", choices=[])
        for i in self.obs_grav_name_list:
            try:
                self.obs_combo_list_grav.Append(i)
            except:
                pass
        for i in self.obs_mag_name_list:
            try:
                self.obs_combo_list_mag.Append(i)
            except:
                pass

        '# DEFINE SET BUTTON'
        self.b_set = wx.Button(input_panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.set_obs, self.b_set)

        '# DEFINE SIZER'
        sizer = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        sizer.AddMany([self.obs_combo_list_grav_text, self.obs_combo_list_grav,
                       self.obs_combo_list_mag_text, self.obs_combo_list_mag, self.b_set])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_obs(self, event):
        obs_grav = str(self.obs_combo_list_grav.GetValue())
        obs_mag = str(self.obs_combo_list_mag.GetValue())

        if obs_grav is not None:
            for i in range(len(self.obs_grav_name_list)):
                '# CHECK GRAV AND SET VALUES'
                if self.obs_grav_name_list[i] == obs_grav:
                    self.obs_gravity_data_for_rms = self.obs_grav_list_save[i]
        else:
            self.obs_gravity_data_for_rms = None

        if obs_mag is not None:
            '# CHECK MAG AND SET VALUES'
            for i in range(len(self.obs_mag_name_list)):
                if self.obs_mag_name_list[i] == obs_mag:
                    self.obs_mag_data_for_rms = self.obs_mag_list_save[i]
        else:
            self.obs_mag_data_for_rms = None
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
        floating_dialogbox = self.SetFloatingLayer(self, -1, "Set new layer nodes")
        answer = floating_dialogbox.ShowModal()
        self.x1, self.y1 = floating_dialogbox.x1, floating_dialogbox.y1
        self.x2, self.y2 = floating_dialogbox.x2, floating_dialogbox.y2
        self.x3, self.y3 = floating_dialogbox.x3, floating_dialogbox.y3
        self.x4, self.y4 = floating_dialogbox.x4, floating_dialogbox.y4
        self.EndModal(1)

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

'''# FRAMES ARE STORED HERE'''


class CaptureCoordinates(wx.Frame):
    """CAPTURE MOUSE CLICK COORDINATES AND WRITE TO DISK FILE. RETURNS ASCII TEXT FILE"""

    def __init__(coordinate_list, parent, id, title):
        wx.Frame.__init__(coordinate_list, None, wx.ID_ANY, 'Capture coordinates', size=(350, 500))
        coordinate_list.input_panel = wx.Panel(coordinate_list)

        # CREATE INSTANCE OF MAIN FRAME CLASS TO RECEIVE NEW ATTRIBUTES
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
        self.dcanvas_aspect = dcanvas_aspect

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

"""START SOFTWARE"""
if __name__ == "__main__":
    app = wx.App(False)
    fr = wx.Frame(None, title='GMG: Geophysical Modelling GUI')
    app.frame = Gmg()
    app.frame.CenterOnScreen()
    app.frame.Show()
    app.MainLoop()
