"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI application for forward modelling 2D potential field profiles with adidtional geophysical data.
Written by Brook Tozer, University of Oxford 2015-17; SIO 2018-19; GNS Science 2020-Present.
Includes ability to import seismic reflection, well, surface outcrop and xy points into the model frame.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Dependencies**

SciPy
NumPy
Matplotlib
pylab
pickle
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
Icons where designed using the free icon maker:

icons8.com: https://icons8.com/icon/set/circle/all

Pixel size: 24
Font: Roboto Slab
Font size: 200
Color: 3498db
***

***
Documentation is created using Sphinx.
***

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NB. If you created a seperate conda env for gmg you must activate it before launching gmg. e.g.::

    source activate py3-gmg

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
import wx.adv
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
import pickle as Pickle
from polygon import Polygon
import plot_model
import bott
import talwani_and_heirtzler
import kim_and_wessel
from frames import *
from dialogs import *
from objects import *
import model_stats
import struct
import gc
import webbrowser

# FUTURE
# from wx.lib.agw import floatspin as fs
# import wx.grid as gridlib
# import time
# from scipy import signal
# from scipy import interpolate as ip
# from colormap import rgb2hex
# import wx.lib.agw.ribbon as RB
# import wx.EnhancedStatusBar as ESB
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class Gmg(wx.Frame):
    """
    Master Class for GMG GUI. Most functions are contained in this Class.
    Upon App launch this sets the panels, sizer's and event bindings.
    Additional classes are used for handling "pop out" windows (Dialog boxes).
    Objects are passed between this "master" GUI class and the Dialog Boxes.
    """

    def __init__(self, *args, **kwds):
        wx.Frame.__init__(self, None, wx.ID_ANY, 'gmg: 2D Geophysical Modelling GUI', size=(1800, 1050))

        # DEFIND ICONS DIRECTORY
        self.gui_icons_dir = os.path.dirname(os.path.abspath(__file__)) + "/icons/"

        # SET AND SHOW SPLASH SCREEN
        bitmap = wx.Bitmap(self.gui_icons_dir + "gmg_logo_scaled.png")
        splash = wx.adv.SplashScreen(bitmap, wx.adv.SPLASH_CENTER_ON_SCREEN | wx.adv.SPLASH_TIMEOUT, 3000,
                                     self, id=wx.ID_ANY, size=(1, 1), style=wx.BORDER_SIMPLE | wx.FRAME_NO_TASKBAR)
        splash.Show()
        self.Show()

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
        self.leftPanel = wx.SplitterWindow(self, wx.ID_ANY, size=(200, 500), style=wx.SP_NOBORDER)
        self.leftPanel_b = wx.SplitterWindow(self.leftPanel, wx.ID_ANY, size=(200, 500),
                                             style=wx.SP_NOBORDER)
        # self.leftPanel.SetMinimumPaneSize(1)
        # self.leftPanel_b.SetMinimumPaneSize(1)

        # FIRST PANE; LEFT PANEL (=ATTRIBUTES)
        self.splitter_left_panel_one = wx.ScrolledWindow(self.leftPanel, wx.ID_ANY, size=(200, 500),
                                                         style=wx.ALIGN_LEFT | wx.BORDER_RAISED)
        self.controls_panel_bar_one = fpb.FoldPanelBar(self.splitter_left_panel_one, 1, size=(200, 500),
                                                       agwStyle=fpb.FPB_VERTICAL)
        self.fold_panel_one = self.controls_panel_bar_one.AddFoldPanel("Layer Attributes", collapsed=True,
                                                                       foldIcons=images)
        self.controls_panel_bar_one.Expand(self.fold_panel_one)  # ENSURES FOLD PANEL IS VISIBLE

        # SECOND PANE; LEFT PANEL (=LAYERS)
        # GREY wx PANEL
        self.splitter_left_panel_two = wx.ScrolledWindow(self.leftPanel_b, wx.ID_ANY, size=(200, 250),
                                                         style=wx.ALIGN_LEFT | wx.BORDER_RAISED)
        # THE LAYER TREE SCROLL BAR GOES IN HERE
        self.controls_panel_bar_two = fpb.FoldPanelBar(self.splitter_left_panel_two, 1, size=(200, 250),
                                                       agwStyle=fpb.FPB_VERTICAL)
        self.fold_panel_two = self.controls_panel_bar_two.AddFoldPanel("Layers", collapsed=True, foldIcons=images)
        self.controls_panel_bar_two.Expand(self.fold_panel_two)  # ENSURES FOLD PANEL IS VISIBLE

        # THIRD PANE; LEFT PANEL (=FAULTS)
        # GREY wx PANEL
        self.splitter_left_panel_three = wx.ScrolledWindow(self.leftPanel_b, wx.ID_ANY, size=(200, 250),
                                                           style=wx.ALIGN_LEFT | wx.BORDER_RAISED)
        # THE FAULT TREE SCROLL BAR GOES IN HERE
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
        m_new_model = menu_file.Append(-1, "New Model...", "New Model...")
        self.Bind(wx.EVT_MENU, self.new_model, m_new_model)
        m_load_model = menu_file.Append(-1, "Load Model...", "Load Model...")
        self.Bind(wx.EVT_MENU, self.load_model, m_load_model)
        m_save_model = menu_file.Append(-1, "Save Model...", "Save Model...")
        self.Bind(wx.EVT_MENU, self.save_model, m_save_model)
        menu_file.AppendSeparator()
        m_save_xy = menu_file.Append(-1, "Save Layers As ASCII .xy File ...",
                                     "Save Layers as ASCII .xy...")
        self.Bind(wx.EVT_MENU, self.write_layers_xy, m_save_xy)
        m_save_c = menu_file.Append(-1, "Save Model As RayInvr c.in File...", "Save RayInvr c.in File...")
        self.Bind(wx.EVT_MENU, self.write_c_xy, m_save_c)
        m_save_fig = menu_file.Append(-1, "Save Figure...", "Save model Figure...")
        self.Bind(wx.EVT_MENU, self.plot_model, m_save_fig)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "Exit...", "Exit...")
        self.Bind(wx.EVT_MENU, self.exit, m_exit)

        # DRAW MENU
        self.menubar.Append(menu_file, "&File")

        # MODEL VIEW MENU
        model_view_file = wx.Menu()
        m_modify_model_dimensions = model_view_file.Append(-1, "Modify Current Model Dimensions...",
                                                           "Modify Current Model Dimensions...")
        self.Bind(wx.EVT_MENU, self.modify_model_dimensions, m_modify_model_dimensions)
        model_view_file.AppendSeparator()

        # PROGRAM FRAME WINDOW SWITCHES
        self.topo_frame_switch = True
        self.gravity_frame_switch = True
        self.magnetic_frame_switch = True
        self.m_model_frames_submenu = wx.Menu()
        model_view_file.AppendSubMenu(self.m_model_frames_submenu, "Toggle Model Frames")
        self.m_model_frames_submenu.Append(601, 'Topography')
        self.Bind(wx.EVT_MENU, self.frame_adjustment, id=601)
        self.m_model_frames_submenu.Append(602, 'Gravity')
        self.Bind(wx.EVT_MENU, self.frame_adjustment, id=602)
        self.m_model_frames_submenu.Append(603, 'Magnetics')
        self.Bind(wx.EVT_MENU, self.frame_adjustment, id=603)
        model_view_file.AppendSeparator()

        m_aspect_increase = model_view_file.Append(-1, "&Increase Aspect Ratio...",
                                                   "Increase the aspect ratio of the model plot...")
        self.Bind(wx.EVT_MENU, self.aspect_increase, m_aspect_increase)
        m_aspect_decrease = model_view_file.Append(-1, "&Decrease Aspect Ratio...",
                                                   "Decrease the aspect ratio of the model plot...")
        self.Bind(wx.EVT_MENU, self.aspect_decrease, m_aspect_decrease)
        self.menubar.Append(model_view_file, "&Model View")

        # GRAVITY DATA MENU --------------------------------------------------------------------------------------------
        self.gravity_data = wx.Menu()
        # LOAD OBSERVED GRAVITY DATA
        m_load_obs_g = self.gravity_data.Append(-1, "&Load Gravity Anomaly...", "Load Observed Gravity Data...")
        self.Bind(wx.EVT_MENU, self.load_obs_g, m_load_obs_g)
        # EDIT
        self.m_obs_g_submenu = wx.Menu()
        self.gravity_data.AppendSubMenu(self.m_obs_g_submenu, "Gravity Data...")
        # FILTER MENU
        grav_m_filter_observed = self.gravity_data.Append(-1, "Median Filter...", "Filter Observed Anomaly")
        self.Bind(wx.EVT_MENU, self.filter_observed_gravity, grav_m_filter_observed)
        # HORIZONTAL DERIVATIVE
        grav_m_horizontal_derivative = self.gravity_data.Append(-1, "Take Horizontal Derivative...",
                                                                "Take Horizontal Derivative")
        self.Bind(wx.EVT_MENU, self.take_gravity_horizontal_derivative, grav_m_horizontal_derivative)
        # SET RMS OBS ARRAYS
        grav_m_set_rms = self.gravity_data.Append(-1, "Set RMS Input...", "Set RMS Input...")
        self.Bind(wx.EVT_MENU, self.set_obs_grav_rms, grav_m_set_rms)
        # SET ELEVATION FOR CALCULATIONS
        m_set_grav_elv = self.gravity_data.Append(-1, "&Set Calculation Elevation...",
                                                  "Set Calculation Elevation...")
        self.Bind(wx.EVT_MENU, self.set_gravity_elv, m_set_grav_elv)
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
        m_load_obs_m = self.magnetic_data.Append(-1, "&Load Magnetic Anomaly...",
                                                 "Load Observed Magnetic Data...")
        self.Bind(wx.EVT_MENU, self.load_obs_m, m_load_obs_m)
        # EDIT
        self.m_obs_mag_submenu = wx.Menu()
        self.magnetic_data.AppendSubMenu(self.m_obs_mag_submenu, "Magnetic Anomalies...")
        # FILTER MENU
        mag_m_filter_observed = self.magnetic_data.Append(-1, "Median Filter...", "Filter Observed Anomaly")
        self.Bind(wx.EVT_MENU, self.filter_observed_magnetic, mag_m_filter_observed)
        # HORIZONTAL DERIVATIVE
        mag_m_horizontal_derivative = self.magnetic_data.Append(-1, "Take Horizontal Derivative...",
                                                                "Take Horizontal Derivative")
        self.Bind(wx.EVT_MENU, self.take_magnetic_horizontal_derivative, mag_m_horizontal_derivative)
        # SET RMS OBS ARRAYS
        mag_m_set_rms = self.magnetic_data.Append(-1, "Set RMS Input..", "Set RMS input..")
        self.Bind(wx.EVT_MENU, self.set_obs_mag_rms, mag_m_set_rms)
        # SET MAG
        m_set_mag_variables = self.magnetic_data.Append(-1, "&Set Calculation Elevation...",
                                                        "Set Calculation Elevation...")
        self.Bind(wx.EVT_MENU, self.set_mag_variables, m_set_mag_variables)
        # SAVE PREDICTED ANOMALY TO DISC
        m_save_mag_submenu = self.magnetic_data.Append(-1, "&Save Predicted Anomaly...",
                                                       "Save Predicted Anomaly to Disc...")
        self.Bind(wx.EVT_MENU, self.save_modelled_mag, m_save_mag_submenu)
        # DRAW MENU
        self.menubar.Append(self.magnetic_data, "&Magnetics")
        # --------------------------------------------------------------------------------------------------------------

        # TOPOGRAPHY DATA MENU -----------------------------------------------------------------------------------------
        self.topography_data = wx.Menu()
        # LOAD OBSERVED TOPO DATA
        m_load_topo = self.topography_data.Append(-1, "&Load Topography...", "Load Observed Topography...")
        self.Bind(wx.EVT_MENU, self.load_topo, m_load_topo)
        # EDIT
        self.m_topo_submenu = wx.Menu()
        self.topography_data.AppendSubMenu(self.m_topo_submenu, "Topography Data")
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
        m_load_xy = self.xy_data.Append(-1, "&Load XY Points...", "Load XY Points...")
        self.Bind(wx.EVT_MENU, self.load_xy, m_load_xy)
        self.m_xy_submenu = wx.Menu()
        self.xy_data.AppendSubMenu(self.m_xy_submenu, "XY Data...")
        # DRAW MENU
        self.menubar.Append(self.xy_data, "&XY Data")
        # --------------------------------------------------------------------------------------------------------------

        # SEISMIC DATA -------------------------------------------------------------------------------------------------
        self.seismic_data = wx.Menu()
        # SEGY LOAD
        self.m_load_segy = self.seismic_data.Append(-1, "&Load Segy...", "Load Segy Data")
        self.Bind(wx.EVT_MENU, self.segy_input, self.m_load_segy)
        # SEGY NAME LIST
        self.m_segy_submenu = wx.Menu()
        self.seismic_data.AppendSubMenu(self.m_segy_submenu, "SEGY Data...")
        # GAIN
        self.m_gain = wx.Menu()
        self.seismic_data.AppendSubMenu(self.m_gain, "Gain")
        # COLOR PALETTE
        self.m_color_palette = wx.Menu()
        self.seismic_data.AppendSubMenu(self.m_color_palette, "Color Palette")
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
        self.well_data.AppendSubMenu(self.m_wells_submenu, "Wells...")
        # DRAW MENU
        self.menubar.Append(self.well_data, "&Well Data")
        # --------------------------------------------------------------------------------------------------------------

        # OUTCROP MENU ------------------------------------------------------------------------------------------------
        self.outcrop_file = wx.Menu()
        self.m_load_outcrop_data = self.outcrop_file.Append(-1, "&Load Outcrop Data...",
                                                            "Load Outcrop Data")
        self.Bind(wx.EVT_MENU, self.load_outcrop_data, self.m_load_outcrop_data)

        self.m_outcrop_submenu = wx.Menu()
        self.outcrop_file.AppendSubMenu(self.m_outcrop_submenu, "Outcrop Data...")

        self.menubar.Append(self.outcrop_file, "&Outcrop Data")
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
        self.layer_file.AppendSubMenu(self.m_layer_transperency, "Transparency")
        # TRANSPARENCY INCREASE
        self.m_layer_transparency_increase = self.m_layer_transperency.Append(-1, "Increase...", "Increase...")
        self.Bind(wx.EVT_MENU, self.transparency_increase, self.m_layer_transparency_increase)
        # TRANSPARENCY DECREASE
        self.m_layer_transparency_decrease = self.m_layer_transperency.Append(-1, "Decrease...", "Decrease...")
        self.Bind(wx.EVT_MENU, self.transparency_decrease, self.m_layer_transparency_decrease)
        # BULK SHIFT
        self.m_bulk_shift = self.layer_file.Append(-1, "Bulk Shift...", "Bulk Shift...")
        self.Bind(wx.EVT_MENU, self.bulk_shift, self.m_bulk_shift)
        # PINCH/DEPINCH LAYER
        self.pinch_submenu = wx.Menu()
        self.layer_file.AppendSubMenu(self.pinch_submenu, "Pinch")
        self.pinch_submenu.Append(701, "&Pinch Out Layer...", "Pinch Out Layer...")
        self.Bind(wx.EVT_MENU, self.pinch_out_layer, id=701)
        self.pinch_submenu.Append(702, "&Depinch Layer...", "Depinch Layer...")
        self.Bind(wx.EVT_MENU, self.depinch_layer, id=702)
        # SEPARATOR
        self.layer_file.AppendSeparator()
        # DELETE LAYER
        self.m_delete_layer = self.layer_file.Append(-1, "Delete Current Layer...", "Delete Current Layer...")
        self.Bind(wx.EVT_MENU, self.delete_layer, self.m_delete_layer)
        # APPEND MENU
        self.menubar.Append(self.layer_file, "&Layers")
        # --------------------------------------------------------------------------------------------------------------

        # ATTRIBUTE TABLE MENU -----------------------------------------------------------------------------------------
        attribute_file = wx.Menu()
        m_attribute_table = attribute_file.Append(-1, "&Open Attribute Table...",
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

        # SET MENUBAR
        self.SetMenuBar(self.menubar)
        # --------------------------------------------------------------------------------------------------------------

        # TOOLBAR - (THIS IS THE ICON BAR BELOW THE MENU BAR)
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

        # GET THE MODEL DIMENSIONS AND SAMPLE LOCATIONS
        self.area = area
        self.x1, self.x2, self.z1, self.z2 = 0.001 * np.array(area)
        self.xp = np.array(xp, dtype='f')

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
        self.run_algorithms()
        self.draw()

    def initalize_model(self):
        """INITIALISE OBSERVED DATA AND LAYERS"""

        # INITIALISE MODEL FRAME ATTRIBUTES
        self.nodes = True  # SWITCH FOR NODE EDITING MODE
        self.zoom_on = False  # SWITCH FOR ZOOM MODE
        self.pan_on = False  # SWITCH PANNING MODE
        self.pinch_switch = False  # SWITCH FOR NODE PINCHING MODE
        self.pinch_count = 0
        self.didnt_get_node = False  # SWITCH FOR MOUSE CLICK (EITHER CLICKED A NODE OR DIDN'T)
        self.node_click_limit = 1  # CONTROLS HOW CLOSE A MOUSE CLICK MUST BE TO ACTIVATE A NODE
        self.index_node = None  # THE ACTIVE NODE
        self.select_new_layer_nodes = False  # SWITCH TO TURN ON WHEN MOUSE CLICK IS TO BE CAPTURED FOR A NEW LAYER
        self.currently_active_layer_id = 0  # LAYER COUNTER
        self.pred_topo = None  # FUTURE - PREDICTED TOPOGRAPHY FROM MOHO (ISOSTATIC FUNC)
        self.predicted_gravity = None  # THE CALCULATED GRAVITY RESPONSE
        self.predicted_nt = None  # THE CALCULATED MAGNETIC RESPONSE
        self.calc_padding = 50000.  # PADDING FOR POTENTIAL FIELD CALCULATION POINTS
        self.padding = 50000.  # PADDING FOR LAYERS
        self.mag_observation_elv = 0.  # OBSERVATION LEVEL FOR MAGNETIC DATA
        self.gravity_observation_elv = 0.  # OBSERVATION LEVEL FOR GRAVITY DATA

        # INITIALISE LAYER LIST
        self.layer_list = []  # LIST HOLDING ALL OF THE LAYER OBJECTS
        self.total_layer_count = 0
        self.layer_transparency = 0.4

        # INITIALISE POLYGON LISTS (USED AS MODEL LAYERS)
        self.gravity_polygons = []
        self.mag_polygons = []
        self.polyplots = []
        self.poly_fills = [[]]

        # INITIALISE FAULT ATTRIBUTES
        self.fault_list = []
        self.currently_actives_fault_id = 0
        self.total_fault_count = 0
        self.fault_picking_switch = False
        self.select_new_fault_nodes = False
        self.selected_node = None

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
        self.gravity_rms_value = None  # TOTAL RMS MISFIT VALUE
        self.grav_residuals = []  # CALCULATED RESIDUAL

        # INITIALISE OBSERVED MAGNETIC ATTRIBUTES
        self.observed_magnetic_list = []
        self.observed_magnetic_counter = 0
        self.observed_magnetic_switch = False
        # INITIALISE MODELLING MAGNETIC ATTRIBUTES
        self.earth_field = 0.

        self.calc_mag_switch = False
        self.obs_mag_data_for_rms = []  # OBSERVED DATA LIST TO BE COMPARED TO CALCULATED
        self.magnetic_rms_value = None  # TOTAL RMS MISFIT VALUE (SINGLE INTEGER)
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

        # INITIALIZE COORDINATE CAPTURE
        self.capture = False
        self.linex = []
        self.liney = []

    def on_wx_key_press(self, evt):
        evt.Skip(False)

    def on_wx_key_release(self, evt):
        evt.Skip(False)

    def draw_main_frame(self):
        """
        DRAW THE GUI FRAMES ( 1. TOPO; 2. GRAVITY; 3. MAGNETICS; 4. MODEL)
        docs: https://matplotlib.org/api/axes_api.html
        """
        self.columns = 87  # NUMBER OF COLUMNS THE MODEL FRAMES WILL TAKE UP (I.E. 87/100)
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
            # CREATE LAYER 0 - THIS IS A PLACE HOLDER FOR THE TOP OF THE MODEL AND NOT ACCESSIBLE BY USER
            layer0 = Layer()
            layer0.type = str('fixed')
            # CREATE THE XY NODES
            layer0.x_nodes = [self.x1 - (float(self.padding)), self.x1, self.x2, self.x2 + (float(self.padding))]
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
            self.current_node = self.model_frame.scatter(-40000., 0, s=50, color='r', zorder=10)  # PLACE HOLDER ONLY

            self.layer_list.append(layer0)

        # ADDITIONAL MAIN FRAME WIDGETS - PLACED ON LEFT HAND SIDE OF THE FRAME

        # MAKE NODE X Y LABEL
        self.node_text = wx.StaticText(self.fold_panel_one, -1, label="Node Position:", style=wx.ALIGN_LEFT)

        # MAKE DENSITY SPINNER
        self.density_text = wx.StaticText(self.fold_panel_one, -1, label="Density:         ", style=wx.ALIGN_LEFT)
        self.density_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-5, max_val=5, increment=0.001, value=0.00)
        self.density_input.SetFormat("%f")
        self.density_input.SetDigits(4)

        # MAKE REFERNECE DENSITY SPINNER
        self.ref_density_text = wx.StaticText(self.fold_panel_one, -1, label="Reference:     \nDensity",
                                              style=wx.ALIGN_LEFT)
        self.ref_density_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-5, max_val=5, increment=0.001,
                                              value=0.00)
        self.ref_density_input.SetFormat("%f")
        self.ref_density_input.SetDigits(4)

        # MAKE SUSCEPTIBILITY SPINNER
        self.susceptibility_text = wx.StaticText(self.fold_panel_one, -1, label="Susceptibility:", style=wx.ALIGN_LEFT)
        self.susceptibility_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-2.0, max_val=1000000.0,
                                                 increment=0.00001, value=0.00, style=wx.ALIGN_RIGHT)
        self.susceptibility_input.SetFormat("%f")
        self.susceptibility_input.SetDigits(6)

        # MAKE ANGLE A SPINNER
        self.angle_a_text = wx.StaticText(self.fold_panel_one, -1, label="Angle A (Inc): ", style=wx.ALIGN_LEFT)
        self.angle_a_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-90.0, max_val=90.0, increment=1.0,
                                          value=0.0, style=wx.ALIGN_RIGHT)
        self.angle_a_input.SetFormat("%f")
        self.angle_a_input.SetDigits(1)

        # MAKE ANGLE B SPINNER
        self.angle_b_text = wx.StaticText(self.fold_panel_one, -1, label="Angle B (Dec):", style=wx.ALIGN_LEFT)
        self.angle_b_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-360.0, max_val=360.0, increment=1.0, value=0.0,
                                          style=wx.ALIGN_RIGHT)
        self.angle_b_input.SetFormat("%f")
        self.angle_b_input.SetDigits(1)

        # MAKE ANGLE C SPINNER
        self.angle_c_text = wx.StaticText(self.fold_panel_one, -1, label="Angle C (Azm):", style=wx.ALIGN_LEFT)
        self.angle_c_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=-360.0, max_val=360.0, increment=1.0,
                                          value=0.0, style=wx.ALIGN_RIGHT)
        self.angle_c_input.SetFormat("%f")
        self.angle_c_input.SetDigits(1)

        # MAKE EARTH FIELD SPINNER
        self.earth_field_text = wx.StaticText(self.fold_panel_one, -1, label="F:", style=wx.ALIGN_LEFT)
        self.earth_field_input = fs.FloatSpin(self.fold_panel_one, -1, min_val=0.0, max_val=1000000.0, increment=1.0,
                                              value=0.0, style=wx.ALIGN_RIGHT)
        self.earth_field_input.SetFormat("%f")
        self.earth_field_input.SetDigits(1)

        # MAKE WELL TEXT SIZE SLIDER
        self.text_size_text = wx.StaticText(self.fold_panel_one, -1, label="Label Text Size:")
        self.text_size_input = wx.Slider(self.fold_panel_one, value=1, minValue=1, maxValue=20, size=(175, -1),
                                         style=wx.SL_HORIZONTAL)
        # MAKE NODE XY SPINNERS
        self.x_text = wx.StaticText(self.fold_panel_one, -1, label="X Value:")
        self.x_input = fs.FloatSpin(self.fold_panel_one, -1, increment=0.001, value=0.00)
        self.x_input.SetDigits(4)
        self.y_text = wx.StaticText(self.fold_panel_one, -1, label="Y Value:")
        self.y_input = fs.FloatSpin(self.fold_panel_one, -1, increment=0.001, value=0.00)
        self.y_input.SetDigits(4)
        # Make Set button
        self.node_set_button = wx.Button(self.fold_panel_one, -1, "Set")

        # MAKE DENSITY CONTRAST SCALEBAR
        # colormap = matplotlib.cm.coolwarm
        # cnorm = colors.Normalize(vmin=-0.8, vmax=0.8)
        # self.cb1 = matplotlib.colorbar.ColorbarBase(self.fold_panel_one, cmap=colormap, norm=cnorm,
        #                                             orientation='horizontal')
        # self.cb1.ax.tick_params(labelsize=6)
        # self.cb1.set_label('Density contrast ($kg/m^{3}$)', fontsize=6)

        # INITALISE CALCULATED ANOMALY LINES
        self.pred_gravity_plot, = self.gravity_frame.plot([], [], '-r', linewidth=2, alpha=0.5)
        self.gravity_rms_plot, = self.gravity_frame.plot([], [], color='purple', linewidth=1.5, alpha=0.5)
        self.predicted_nt_plot, = self.magnetic_frame.plot([], [], '-g', linewidth=2, alpha=0.5)
        self.mag_rms_plot, = self.magnetic_frame.plot([], [], color='purple', linewidth=1.5, alpha=0.5)

        # MAKE LAYER TREE
        self.tree = ct.CustomTreeCtrl(self.fold_panel_two, -1,
                                      agwStyle=wx.TR_DEFAULT_STYLE | wx.TR_ROW_LINES | wx.TR_HIDE_ROOT)
        self.tree.SetIndent(0.0)
        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.on_layer_activated, self.tree)
        self.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.on_tree_right_click_down, self.tree)

        # TREE ATTRIBUTES
        self.root = self.tree.AddRoot("Layers:")
        self.tree.SetItemPyData(self.root, None)
        self.tree_items = ["Layer 0"]
        self.Bind(ct.EVT_TREE_ITEM_CHECKED, self.item_checked, self.tree)

        # MAKE FAULT TREE
        self.fault_tree = ct.CustomTreeCtrl(self.fold_panel_three, -1,
                                            agwStyle=wx.TR_DEFAULT_STYLE | wx.TR_ROW_LINES | wx.TR_HIDE_ROOT)
        self.fault_tree.SetIndent(0.0)
        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.on_fault_activated, self.fault_tree)
        self.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.on_fault_tree_right_click_down, self.fault_tree)

        # TREE ATTRIBUTES
        self.fault_tree_root = self.fault_tree.AddRoot("Faults:")
        self.fault_tree.SetItemPyData(self.fault_tree_root, None)
        self.fault_tree_items = []
        self.Bind(ct.EVT_TREE_ITEM_CHECKED, self.fault_checked, self.fault_tree)

        # UPDATE INFO BAR
        self.display_info()

        # self.blocker= wx.EventBlocker(self.rightPanel, wx.Bell())

        # REDRAW MAIN
        self.draw()

    def size_handler(self):
        """PLACE THE GUI FRAMES IN THE wxSIZER WINDOWS"""

        # --------------------------------------------------------------------------------------------------------------
        # POPULATE ATTRIBUTES PANEL (LEFT PANEL OF GUI)
        self.attributes_box = wx.GridBagSizer(hgap=2, vgap=3)
        r = 1  # CURRENT ROW
        c = 0  # CURRENT COLUMN

        # LINE SEPx
        line = wx.StaticLine(self.fold_panel_one)
        self.attributes_box.Add(line, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND | wx.EXPAND, border=5)

        # DENSITY
        r += 1
        self.attributes_box.Add(self.density_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=5)
        c += 1
        self.attributes_box.Add(self.density_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=5)

        # LINE SEP
        r += 1
        c += 0
        line = wx.StaticLine(self.fold_panel_one)
        self.attributes_box.Add(line, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND, border=5)

        # REFERENCE DENSITY
        r += 1
        c = 0
        self.attributes_box.Add(self.ref_density_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=5)
        c += 1
        self.attributes_box.Add(self.ref_density_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=5)

        # LINE SEP
        r += 1
        c = 0
        line = wx.StaticLine(self.fold_panel_one)
        self.attributes_box.Add(line, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND, border=5)

        # SUSCEPTIBILITY
        r += 1
        c = 0
        self.attributes_box.Add(self.susceptibility_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)
        c += 1
        self.attributes_box.Add(self.susceptibility_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)

        # ANGLE A
        r += 1
        c = 0
        self.attributes_box.Add(self.angle_a_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)
        c += 1
        self.attributes_box.Add(self.angle_a_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)

        # ANGLE B
        r += 1
        c = 0
        self.attributes_box.Add(self.angle_b_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)
        c += 1
        self.attributes_box.Add(self.angle_b_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)

        # Angle C
        r += 1
        c = 0
        self.attributes_box.Add(self.angle_c_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)
        c += 1
        self.attributes_box.Add(self.angle_c_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)

        # Earth Field
        r += 1
        c = 0
        self.attributes_box.Add(self.earth_field_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)
        c += 1
        self.attributes_box.Add(self.earth_field_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)

        # LINE SEP
        r += 1
        c = 0
        line = wx.StaticLine(self.fold_panel_one)
        self.attributes_box.Add(line, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND, border=5)

        # XY NODES
        r += 1
        c = 0
        self.attributes_box.Add(self.node_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)

        # LINE SEP
        r += 1
        c = 0
        line = wx.StaticLine(self.fold_panel_one)
        self.attributes_box.Add(line, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND, border=5)

        # X NODE
        r += 1
        c = 0
        self.attributes_box.Add(self.x_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)
        c = + 1
        self.attributes_box.Add(self.x_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)

        # Y NODE
        r += 1
        c = 0
        self.attributes_box.Add(self.y_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)
        c = + 1
        self.attributes_box.Add(self.y_input, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                                border=5)

        # LINE SEP
        r += 1
        c = 0
        line = wx.StaticLine(self.fold_panel_one)
        self.attributes_box.Add(line, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND, border=5)

        # SET BUTTON
        r += 1
        c = 0
        self.attributes_box.Add(self.node_set_button, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND,
                                border=5)

        # LINE SEP
        r += 1
        c = 0
        line = wx.StaticLine(self.fold_panel_one)
        self.attributes_box.Add(line, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND, border=5)

        # LABEL TEXT SIZE
        r += 1
        c = 0
        self.attributes_box.Add(self.text_size_text, pos=(r, c), span=(1, 2), flag=wx.ALIGN_CENTER,
                                border=5)
        r += 1
        c = 0
        self.attributes_box.Add(self.text_size_input, pos=(r, c), span=(1, 2), flag=wx.ALIGN_CENTER,
                                border=5)

        # LINE SEP
        r += 1
        c = 0
        line = wx.StaticLine(self.fold_panel_one)
        self.attributes_box.Add(line, pos=(r, c), span=(1, 2), flag=wx.ALIGN_LEFT | wx.EXPAND, border=5)

        # DENSITY SCALE BAR
        # self.attr_box.Add(self.cb1, 0, wx.ALL | wx.LEFT | wx.EXPAND, 5)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # CREATE LAYER TREE BOX
        self.tree_box = wx.BoxSizer(wx.VERTICAL)
        self.tree_box.Add(self.tree, 1, wx.TOP | wx.ALIGN_CENTER | wx.EXPAND, border=20)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # CREATE FAULT TREE BOX
        self.fault_tree_box = wx.BoxSizer(wx.VERTICAL)
        self.fault_tree_box.Add(self.fault_tree, 1, wx.TOP | wx.ALIGN_CENTER | wx.EXPAND, border=20)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # CREATE A BOX SIZER FOR THE MAIN MODELLING FRAME
        self.canvas_box = wx.BoxSizer(wx.HORIZONTAL)
        # ADD THE MAIN MODELLING FRAME TO IT'S A BOX SIZER
        self.canvas_box.Add(self.canvas, 1, wx.ALL | wx.ALIGN_CENTER | wx.EXPAND, border=2)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # PLACE BOX SIZERS IN CORRECT PANELS
        self.fold_panel_one.SetSizerAndFit(self.attributes_box)
        self.fold_panel_two.SetSizerAndFit(self.tree_box)
        self.fold_panel_three.SetSizerAndFit(self.fault_tree_box)
        self.leftPanel.SetSizer(self.splitter_left_panel_sizer)
        self.fold_panel_two.Collapse()
        self.fold_panel_two.Expand()
        self.fold_panel_three.Collapse()
        self.fold_panel_three.Expand()
        self.fold_panel_one.Collapse()
        self.fold_panel_one.Expand()

        self.rightPanel.SetSizerAndFit(self.canvas_box)
        self.rightPanel.SetSize(self.GetSize())
        # --------------------------------------------------------------------------------------------------------------

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
            # TRUE TRUE TRUE

            # TOPO CANVAS
            self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=2, colspan=self.columns)
            self.topo_frame.set_ylabel("(m)")
            self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_frame.grid()
            self.topo_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.topo_d_frame = self.topo_frame.twinx()
            self.topo_d_frame.set_ylabel("dt/dx")
            self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            # GRAV CANVAS
            self.gravity_frame = plt.subplot2grid((26, 100), (2, self.x_orig), rowspan=3, colspan=self.columns)
            self.gravity_frame.set_ylabel("(mGal)")
            self.gravity_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_frame.grid()
            self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.gravity_d_frame = self.gravity_frame.twinx()
            self.gravity_d_frame.set_ylabel("dg/dx")
            self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            # MAG CANVAS
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
            # FALSE TRUE TRUE

            # TOPO CANVAS - HIDDEN

            # GRAV CANVAS
            self.gravity_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=4, colspan=self.columns)
            self.gravity_frame.set_ylabel("(mGal)")
            self.gravity_frame.grid()
            self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.gravity_d_frame = self.gravity_frame.twinx()
            self.gravity_d_frame.set_ylabel("dg/dx")
            self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            # MAG CANVAS
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
            # TRUE FALSE TRUE

            # TOPO CANVAS
            self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=2, colspan=self.columns)
            self.topo_frame.set_ylabel("(m)")
            self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_frame.grid()

            self.topo_d_frame = self.topo_frame.twinx()
            self.topo_d_frame.set_ylabel("dt/dx")
            self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            # GRAV CANVAS - HIDDEN

            # MAG CANVAS
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
            # TRUE TRUE FALSE
            # TOPO CANVAS
            self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=2, colspan=self.columns)
            self.topo_frame.set_ylabel("(m)")
            self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_frame.grid()

            self.topo_d_frame = self.topo_frame.twinx()
            self.topo_d_frame.set_ylabel("dt/dx")
            self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            # GRAV CANVAS
            self.gravity_frame = plt.subplot2grid((26, 100), (2, self.x_orig), rowspan=6, colspan=self.columns)
            self.gravity_frame.set_ylabel("(mGal)")
            self.gravity_frame.grid()
            self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.gravity_d_frame = self.gravity_frame.twinx()
            self.gravity_d_frame.set_ylabel("dg/dx")
            self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            # MAG CANVAS - HIDDEN

        elif self.topo_frame_switch is False and self.gravity_frame_switch is False and \
                self.magnetic_frame_switch is True:
            # FALSE FALSE TRUE'
            # TOPO CANVAS - HIDDEN

            # GRAV CANVAS - HIDDEN

            # MAG CANVAS'
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
            # FALSE TRUE FALSE

            # TOPO CANVAS - HIDDEN

            # GRAV CANVAS
            self.gravity_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=8, colspan=self.columns)
            self.gravity_frame.set_ylabel("(mGal)")
            self.gravity_frame.grid()
            self.gravity_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            self.gravity_d_frame = self.gravity_frame.twinx()
            self.gravity_d_frame.set_ylabel("dg/dx")
            self.gravity_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.gravity_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            # MAG CANVAS -HIDDEN

        elif self.topo_frame_switch is True and self.gravity_frame_switch is False and \
                self.magnetic_frame_switch is False:
            # TRUE FALSE FALSE

            # TOPO CANVAS
            self.topo_frame = plt.subplot2grid((26, 100), (0, self.x_orig), rowspan=8, colspan=self.columns)
            self.topo_frame.set_ylabel("(m)")
            self.topo_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_frame.grid()

            self.topo_d_frame = self.topo_frame.twinx()
            self.topo_d_frame.set_ylabel("dt/dx")
            self.topo_d_frame.xaxis.set_major_formatter(plt.NullFormatter())
            self.topo_d_frame.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            # GRAV CANVAS - HIDDEN

            # MAG CANVAS - HIDDEN

        elif self.topo_frame_switch is False and self.gravity_frame_switch is False and \
                self.magnetic_frame_switch is False:
            pass
            # FALSE FALSE FALSE
            # TOPO CANVAS - HIDDEN
            # GRAV CANVAS - HIDDEN
            # MAG CANVAS - HIDDEN

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
            self.pred_gravity_plot, = self.gravity_frame.plot([], [], '-r', linewidth=2, alpha=0.5)
            self.gravity_rms_plot, = self.gravity_frame.plot([], [], color='purple', linewidth=1.5, alpha=0.5)
        if self.magnetic_frame is not None:
            self.predicted_nt_plot, = self.magnetic_frame.plot([], [], '-g', linewidth=2, alpha=0.5)
            self.mag_rms_plot, = self.magnetic_frame.plot([], [], color='purple', linewidth=1.5, alpha=0.5)

        # PLOT OBSERVED TOPO DATA
        if self.topo_frame is not None:
            # REPLOT OBSERVED TOPOGRAPHY DATA
            for x in range(len(self.observed_topography_list)):
                if self.observed_topography_list[x] is not None:
                    # DRAW DATA IN MODEL FRAME
                    if self.observed_topography_list[x].type != "derivative":
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
                    if self.observed_gravity_list[x].type != "derivative":
                        self.observed_gravity_list[x].mpl_actor = self.gravity_frame.scatter(
                            self.observed_gravity_list[x].data[:, 0],
                            self.observed_gravity_list[x].data[:, 1],
                            marker='o', s=5,
                            color=self.observed_gravity_list[x].color,
                            gid=11000 + self.observed_gravity_list[x].id)
                    else:
                        self.observed_gravity_list[x].mpl_actor = self.gravity_d_frame.scatter(
                            self.observed_gravity_list[x].data[:, 0],
                            self.observed_gravity_list[x].data[:, 1],
                            marker='o', s=5,
                            color=self.observed_gravity_list[x].color,
                            gid=11000 + self.observed_gravity_list[x].id)

        # PLOT OBSERVED MAGNETIC DATA
        if self.magnetic_frame is not None:
            # REPLOT OBSERVED MAGNETIC DATA
            for x in range(len(self.observed_magnetic_list)):
                if self.observed_magnetic_list[x] is not None:
                    # DRAW DATA IN MODEL FRAME
                    if self.observed_magnetic_list[x].type != "derivative":
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
            self.initalize_model()

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
        self.gravity_observation_elv = np.array(zp, dtype='f')
        self.update_layer_data()
        self.run_algorithms()

    def connect(self):
        """CONNECT MOUSE AND EVENT BINDINGS"""

        # CONNECT MPL INTERACTIONS
        self.fig.canvas.mpl_connect('button_press_event', self.button_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self.move)
        self.fig.canvas.mpl_connect('button_release_event', self.button_release)
        self.fig.canvas.mpl_connect('key_press_event', self.key_press)
        # self.fig.canvas.mpl_connect('pick_event', self.on_pick)

        # CONNECT wx.widgets
        self.density_input.Bind(fs.EVT_FLOATSPIN, self.set_density)
        self.ref_density_input.Bind(fs.EVT_FLOATSPIN, self.set_reference_density)
        self.susceptibility_input.Bind(fs.EVT_FLOATSPIN, self.set_susceptibility)
        self.angle_a_input.Bind(fs.EVT_FLOATSPIN, self.set_angle_a)
        self.angle_b_input.Bind(fs.EVT_FLOATSPIN, self.set_angle_b)
        self.angle_c_input.Bind(fs.EVT_FLOATSPIN, self.set_angle_c)
        self.earth_field_input.Bind(fs.EVT_FLOATSPIN, self.set_earth_field)
        self.text_size_input.Bind(wx.EVT_SLIDER, self.set_text_size)
        self.node_set_button.Bind(wx.EVT_BUTTON, self.on_menu_set_button_press)

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

    def on_layer_activated(self, event):
        """WHEN A LAYER IS SELECTED IN THE TREE"""
        # FIRST CHECK IS FAULT PICKING MODE IS ON, IF IT IS, THEN TURN IT OFF
        if self.fault_picking_switch is True:
            self.fault_picking_switch = False

        # GET HTE TREE ID AND USE IT TO SET THE CURRENT LAYER
        self.currently_active_layer_id = self.tree.GetPyData(event.GetItem())

        # SET OBJECTS WITH THE CHOSEN LAYER
        self.density_input.SetValue(0.001 * self.layer_list[self.currently_active_layer_id].density)
        self.ref_density_input.SetValue(0.001 * self.layer_list[self.currently_active_layer_id].reference_density)
        self.susceptibility_input.SetValue(self.layer_list[self.currently_active_layer_id].susceptibility)
        self.angle_a_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_a)
        self.angle_b_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_b)
        self.angle_c_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_c)
        self.earth_field_input.SetValue(self.layer_list[self.currently_active_layer_id].earth_field)

        # GET THE XY NODES FROM THE ACTIVE LAYER AND SET THE CUURENTLY ACTIVE NODES (I.E. MAKE THEM INTERACTIVE)
        self.current_x_nodes = self.layer_list[self.currently_active_layer_id].x_nodes
        self.current_y_nodes = self.layer_list[self.currently_active_layer_id].y_nodes

        # SET CURRENTLY ACTIVE (RED) NODE
        self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])

        # UPDATE MODEL
        self.update_layer_data()
        self.run_algorithms()

    def on_tree_right_click_down(self, event):
        """WHEN A LAYER IN THE LAYER TREE MENU IS RIGHT CLICKED"""

        # FIRST RUN THE on_layer_activated FUNC
        self.on_layer_activated(event)

        # CREATE POPOUT MENU WITH OPTIONS AND BIND OPTIONS TO ACTIONS
        menu = wx.Menu()

        item1 = menu.Append(wx.ID_ANY, "Change layer colour")
        item2 = menu.Append(wx.ID_ANY, "Rename layer")

        self.Bind(wx.EVT_MENU, self.change_color, item1)
        self.Bind(wx.EVT_MENU, self.rename_layer, item2)

        self.PopupMenu(menu)
        menu.Destroy()

    def change_color(self, event):
        """SET COLOUR FOR LAYER"""

        # CREATE DIALOG
        starting_color = wx.ColourData()
        starting_color.SetColour('white')
        dlg = wx.ColourDialog(self, starting_color)

        # ENSURE THE FULL COLOUR DIALOG IS DISPLAYED, NOT THE ABBREVIATED VERSION
        dlg.GetColourData().SetChooseFull(True)
    
        # WHAT TO DO WHEN THE DIALOG IS CLOSED
        if dlg.ShowModal() == wx.ID_OK:
            rgb = dlg.GetColourData().GetColour().Get()
            rgb = rgb[0:3]
            html = bytes.hex(struct.pack('BBB',*rgb))

            # SET FAULT OR LAYER COLOR
            if self.fault_picking_switch == True:
                self.fault_list[self.currently_active_fault_id].color = str('#' + str(html))
            else:
                self.layer_list[self.currently_active_layer_id].color = str('#' + str(html))

        # CLOSE DIALOG
        dlg.Destroy()

        # REDRAW
        self.update_layer_data()
        self.draw()

    def rename_layer(self, event):
        """USE A POPUP MENU TO RENAME THE LAYER"""

        # CREATE POP OUT MENU AND SHOW
        layer_name_box = LayerNameDialog(self, -1, 'Rename layer', self.tree_items[self.currently_active_layer_id])
        new = layer_name_box.ShowModal()

        # WAIT FOR USER TO CLOSE POP OUT

        # GET THE NEW LAYER NAME FROM POP OUT
        new_name = layer_name_box.name

        # SET THE TREE AND LAYER OBJECT WITH THE NEW NAME
        current_tree_items = self.tree.GetRootItem().GetChildren()
        self.tree.SetItemText(current_tree_items[self.currently_active_layer_id - 1], str(new_name))
        self.tree_items[self.currently_active_layer_id] = str(new_name)
        self.layer_list[self.currently_active_layer_id].name = str(new_name)

    def delete_all_children(self, event):
        self.tree.DeleteChildren(event)

    def item_checked(self, event):
        """TOGGLE WHETHER OR NOT THE SELECTED LAYER IS INCLUDED IN THE MODELLED ANOMALY CALCULATIONS"""
        layer = self.tree.GetPyData(event.GetItem())

        checked_value = self.tree.GetRootItem().GetChildren()[layer - 1].GetValue()

        if checked_value is True:
            self.layer_list[layer].include_in_calculations_switch = True
        else:
            self.layer_list[layer].include_in_calculations_switch = False

        # UPDATE MODEL
        self.run_algorithms()
        self.draw()

    def display_info(self):
        self.statusbar.SetStatusText("                                                                                 "
                                     "                                                   "
                                     " || Currently Editing Layer: %s  || "
                                     " || Model Aspect Ratio = %s:1.0  || GRAV RMS = %s "
                                     " || MAG RMS = %s  ||" % (self.currently_active_layer_id,
                                                               self.model_frame.get_aspect(), self.gravity_rms_value,
                                                               self.magnetic_rms_value), 2)
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
                self.gravity_rms_plot.set_data(self.grav_residuals[:, 0], self.grav_residuals[:, 1])
            self.gravity_rms_plot.set_visible(False)
        else:
            self.calc_grav_switch = True
            self.gravity_rms_plot.set_visible(True)
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def calc_mag_switch(self, event):
        """PREDICTED ANOMALY CALCULATION SWITCH: ON/OFF (SPEEDS UP PROGRAM WHEN OFF)"""
        if self.calc_mag_switch is True:
            self.calc_mag_switch = False
            if self.mag_residuals != [] and self.obs_mag_data_for_rms != []:
                self.mag_residuals[:, 1] = np.zeros(len(self.obs_mag_data_for_rms[:, 0]))
                self.mag_rms_plot.set_visible(False)
        else:
            self.calc_mag_switch = True
            self.mag_rms_plot.set_visible(True)
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
        for l in range(0, len(self.layer_list)):
            if self.layer_list[l].polygon_mpl_actor[0] is not None \
                    and self.layer_list[l].polygon_mpl_actor[0].get_alpha() <= 0.9:
                new_alpha = self.layer_list[l].polygon_mpl_actor[0].get_alpha() + 0.1
                self.layer_list[l].polygon_mpl_actor[0].set_alpha(new_alpha)
        self.draw()

    def transparency_decrease(self, event):
        for l in range(0, len(self.layer_list)):
            if self.layer_list[l].polygon_mpl_actor[0] is not None \
                    and self.layer_list[l].polygon_mpl_actor[0].get_alpha() >= 0.1:
                new_alpha = self.layer_list[l].polygon_mpl_actor[0].get_alpha() - 0.1
                self.layer_list[l].polygon_mpl_actor[0].set_alpha(new_alpha)
        self.draw()

    # SAVE/LOAD MODEL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def save_model(self, event):
        """
        SAVE MODEL TO DISC IN .Pickle FORMAT
        TO ADD NEW OBJECTS ADD THE OBJECT NAME TO BOTH THE header AND model_params.
        THEN MODIFY THE load_model FUNCTION TO ALSO INLCUDE THE NEW ITEMS
        """
        save_file_dialog = wx.FileDialog(self, "Save model file", "", "", "Model files (*.model)|*.model", 
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # USER CHANGED THEIR MIND

        # CREATE SAVE DICTIONARY
        self.save_dict = {}

        header = ['model_aspect', 'area', 'xp',
                  'tree_items', 'fault_tree_items',
                  'gravity_observation_elv', 'mag_observation_elv',
                  'obs_gravity_data_for_rms', 'obs_mag_data_for_rms',
                  'layer_list',
                  'fault_list',
                  'observed_xy_data_list',
                  'observed_gravity_list',
                  'observed_magnetic_list',
                  'observed_topography_list',
                  'well_data_list',
                  'outcrop_data_list',
                  'segy_data_list']

        model_params = [self.model_aspect, self.area, self.xp,
                        self.tree_items, self.fault_tree_items,
                        self.gravity_observation_elv, self.mag_observation_elv,
                        self.obs_gravity_data_for_rms, self.obs_mag_data_for_rms,
                        self.layer_list,
                        self.fault_list,
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
                print((header[i]))
        try:
            output_stream = save_file_dialog.GetPath()
            with open(output_stream, 'wb') as output_file:
                Pickle.dump(self.save_dict, output_file, protocol=Pickle.HIGHEST_PROTOCOL)
            self.model_saved = True
            self.update_layer_data()

            # DISPLAY MESSAGE
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
            self.initalize_model()
            self.model_aspect = 1.

            # OPEN DATA STREAM
            with open(open_file_dialog.GetPath(), 'rb') as input_file:
                model_data = Pickle.load(input_file)

            # CLEAR MEMORY
            gc.collect()
            del gc.garbage[:]

            # LOAD DATA INTO MODEL
            for x in range(len(model_data)):
                setattr(self, list(model_data.keys())[x], list(model_data.values())[x])

            # SAVE LOADED TREE ITEMS (WILL BE REMOVED BY self.start)
            self.loaded_tree_items = self.tree_items
            self.loaded_fault_tree_items = self.fault_tree_items

            # DRAW CANVAS
            self.start(self.area, self.xp, self.gravity_observation_elv)

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
            # SET VARIABLE VALUES FROM LOADED DATA - USING LAYER 1
            # SET LAYER ATTRIBUTES SIDE BAR
            if len(self.layer_list) > 1:
                self.currently_active_layer_id = 1
                self.total_layer_count = len(self.layer_list) - 1
                self.density_input.SetValue(0.001 * self.layer_list[1].density)
                self.ref_density_input.SetValue(0.001 * self.layer_list[1].reference_density)
                self.current_x_nodes = self.layer_list[1].x_nodes
                self.current_y_nodes = self.layer_list[1].y_nodes
            else:
                self.currently_active_layer_id = 0
                self.total_layer_count = 0
                self.density_input.SetValue(0.0)
                self.ref_density_input.SetValue(0.0)
                self.current_x_nodes = self.layer_list[0].x_nodes
                self.current_y_nodes = self.layer_list[0].y_nodes   
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            # CREATE LAYERS
            if len(self.layer_list) > 1:
                for l in range(0, self.total_layer_count + 1):
                    self.layer_list[l].node_mpl_actor = self.model_frame.plot(self.layer_list[l].x_nodes,
                                                                            self.layer_list[l].y_nodes, color='blue',
                                                                            linewidth=1.0, alpha=1.0)
                    # CREATE LAYER POLYGON FILL
                    self.layer_list[l].polygon_mpl_actor = self.model_frame.fill(self.layer_list[l].x_nodes,
                                                                                self.layer_list[l].y_nodes, color='blue',
                                                                                alpha=self.layer_transparency,
                                                                                closed=True, linewidth=None, ec=None)

                self.currently_active_layer, = self.model_frame.plot(self.current_x_nodes, self.current_y_nodes, marker='o',
                                                                    color=self.layer_list[
                                                                        self.currently_active_layer_id].color,
                                                                    linewidth=1.0,
                                                                    alpha=0.5)
            else:
                self.currently_active_layer, = self.model_frame.plot(self.current_x_nodes, self.current_y_nodes, marker='o',
                                                                    color=self.layer_list[0].color,
                                                                    linewidth=1.0,
                                                                    alpha=0.5)
                
            # SET CURRENT NODE AS A OFF STAGE (PLACE HOLDER)
            self.current_node = self.model_frame.scatter(-40000., 0., marker='o', color='r', zorder=10)
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

            self.tree_items = self.loaded_tree_items
            # ----------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------
            self.total_fault_count = len(self.fault_list)

            for i in range(0, self.total_fault_count):
                # DRAW FAULTS
                self.fault_list[i].mpl_actor = self.model_frame.plot(self.fault_list[i].x_nodes,
                                                                     self.fault_list[i].y_nodes,
                                                                     color=self.fault_list[i].color,
                                                                     linewidth=0.5, zorder=1, marker='s', alpha=1.0)
            # CREATE NEW CURRENT FAULT GRAPHIC
            self.currently_active_fault, = self.model_frame.plot([-100000, -100000], [-100000, -100000], marker='s',
                                                                 color='g', linewidth=0.75, alpha=1.0, zorder=2,
                                                                 picker=True)

            # SET FAULT PICKING SWITCH OFF (DEFAULT TO LAYER MODE)
            self.fault_picking_swtich = False
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
            if 1 in self.layer_list:
                self.update_layer_data()  # IF THE MODEL INCLUDES LAYERS THEN UPDATE
            self.run_algorithms()
            self.draw()
            self.Restore()  # FIX'S DISPLAY ISSUE
            self.fold_panel_two.SetSize(200, 300)
            self.fold_panel_three.SetSize(200, 300)


            # ----------------------------------------------------------------------------------------------------------

        # HANDLE MODEL LOAD ERRORS
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
        self.update_layer_data()  # UPDATE LAYER DATA
    
    def replot_observed_xy_data(self):
        """ADD LOADED OBSERVED XY DATA TO THE MODEL FRAME"""
        for x in range(len(self.observed_xy_data_list)):
            if self.observed_xy_data_list[x] is not None:
                # DRAW DATA IN MODEL FRAME

                self.observed_xy_data_list[x].mpl_actor = self.model_frame.scatter(
                    self.observed_xy_data_list[x].data[:, 0], self.observed_xy_data_list[x].data[:, 1],
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
                self.m_obs_g_submenu.Append(11000 + self.observed_gravity_list[x].id,
                                            self.observed_gravity_list[x].name,
                                            self.obs_submenu)
                self.obs_submenu.Append(11000 + self.observed_gravity_list[x].id, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000 + self.observed_gravity_list[x].id)

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
                self.m_obs_mag_submenu.Append(12000 + self.observed_magnetic_list[x].id,
                                              self.observed_magnetic_list[x].name,
                                              self.mag_submenu)
                self.mag_submenu.Append(12000 + self.observed_magnetic_list[x].id, 'delete observed data')
                self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12000 + self.observed_magnetic_list[x].id)

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
                    well.horizons_list[i - 2] = self.model_frame.plot(x, y, linestyle='-', linewidth='2', color='black')
                    horizon_y_pos = well.data[i][1].astype(float)
                    horizon = well.data[i][0].astype(str)

                    # ALTERNATE POSITION OF ODDs/EVENs TO TRY AND AVOID OVERLAP
                    if i % 2 == 0:
                        horizon_x_pos = well.data[1][1].astype(float) - 1.05
                        well.labels_list[i - 2] = self.model_frame.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                            xytext=(horizon_x_pos, horizon_y_pos),
                                                                            fontsize=well.text_size, weight='bold',
                                                                            horizontalalignment='left',
                                                                            verticalalignment='top',
                                                                            color='black',
                                                                            bbox=dict(boxstyle="round,pad=.4",
                                                                                      fc="0.8", ec='None'),
                                                                            clip_on=True)
                    else:
                        horizon_x_pos = well.data[1][1].astype(float) + 1.05
                        well.labels_list[i - 2] = self.model_frame.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                            xytext=(horizon_x_pos, horizon_y_pos),
                                                                            fontsize=well.text_size, weight='bold',
                                                                            horizontalalignment='right',
                                                                            verticalalignment='top',
                                                                            color='black',
                                                                            bbox=dict(boxstyle="round,pad=.4",
                                                                                      fc="0.8", ec='None'),
                                                                            clip_on=True)

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
                text = list(zip(outcrop.data[:, 0].astype(float), outcrop.data[:, 1].astype(float),
                                outcrop.data[:, 3].astype(str)))

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

    # TOPOGRAPHY DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
        self.m_topo_submenu.Append(10000 + observed_topography.id, observed_topography.name, self.topo_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.topo_submenu.Append(10000 + observed_topography.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_observed_topography, id=10000 + observed_topography.id)

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
        obj_id = event.Id - 10000
        self.observed_topography_list[obj_id].mpl_actor.set_visible(False)
        self.observed_topography_list[obj_id] = None

        # UPDATE MODEL
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

        # APPEND NEW DATA MENU TO 'TOPO data MENU'
        self.topo_submenu = wx.Menu()
        self.m_topo_submenu.Append(10000 + observed.id, observed.name, self.topo_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.topo_submenu.Append(10000 + observed.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_observed_topography, id=10000 + observed.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_topography_counter += 1

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

        # INCREMENT TOPO DERIV COUNTER
        self.observed_topography_counter += 1

        # UPDATE GMG GUI
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
        self.m_obs_g_submenu.Append(11000 + observed_gravity.id, observed_gravity.name, self.grav_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.grav_submenu.Append(11000 + observed_gravity.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000 + observed_gravity.id)

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
        obj_id = event.Id - 11000
        self.observed_gravity_list[obj_id].mpl_actor.set_visible(False)
        self.observed_gravity_list[obj_id] = None

        # UPDATE MODEL
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def set_gravity_elv(self, event):
        """POPOUT BOX TO LET USER DEFINE THE ELEVATION AT WHICH TO CALCULATE THE GRAVITY ANOMALY"""

        # CREATE THE POPOUT BOX FOR USER UNPUT
        grav_box = GravDialog(self, -1, 'Gravity elevation', self.gravity_observation_elv)
        answer = grav_box.ShowModal()

        # SET THE NEW CALCULATION ELEVATION
        self.gravity_observation_elv = grav_box.grav_observation_elv * 1000.  # CONVERT FROM (km) TO (m)

        # UPDATE GMG
        self.run_algorithms()
        self.draw()

    def save_modelled_grav(self, event):
        """SAVE PREDICTED GRAVITY TO EXTERNAL ASCII FILE"""
        save_file_dialog = wx.FileDialog(self, "Save Predicted Anomaly", "", "", "Predicted Anomaly (*.txt)|*.txt",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THEIR MIND
        # SAVE TO DISC
        outputfile = save_file_dialog.GetPath()
        np.savetxt(outputfile, list(zip((self.xp * 0.001), self.predicted_gravity)), delimiter=' ', fmt='%.6f %.6f')

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
        self.m_obs_g_submenu.Append(11000 + observed.id, observed.name, self.grav_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.grav_submenu.Append(11000 + observed.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_grav, id=11000 + observed.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_gravity_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
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
        self.m_obs_mag_submenu.Append(12000 + observed_magnetic.id, observed_magnetic.name, self.mag_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.mag_submenu.Append(12000 + observed_magnetic.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12000 + observed_magnetic.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_magnetic_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def delete_obs_mag(self, event):
        """DELETE A MAGNETIC ANOMALY FROM THE MAGNETIC DATA PANEL"""

        # DESTROY MENUBAR
        self.m_obs_mag_submenu.DestroyItem(event.Id)

        # REMOVE OBJECT AND MPL ACTOR
        obj_id = event.Id - 12000
        self.observed_magnetic_list[obj_id].mpl_actor.set_visible(False)
        self.observed_magnetic_list[obj_id] = None

        # UPDATE MODEL
        self.update_layer_data()
        self.set_frame_limits()
        self.draw()

    def set_mag_variables(self, event):
        """POPOUT BOX TO LET USER DEFINE MAGNETIC FIELD VALUES"""

        # CREATE POPOUT MENU
        mag_box = MagDialog(self, -1, 'Magnetic parameters', self.mag_observation_elv)
        answer = mag_box.ShowModal()

        # UPDATE MAGNETIC OBSERVATION ELEVATION
        self.mag_observation_elv = mag_box.mag_observation_elv * 1000.  # CONVERT FROM (km) TO (m)

        # UPDATE GMG
        self.run_algorithms()
        self.draw()

    def save_modelled_mag(self, event):
        """SAVE THE MODELLED MAGNETIC ANOMALY TO DISC"""

        save_file_dialog = wx.FileDialog(self, "Save Predicted Anomaly", "", "", "Predicted Anomaly (*.txt)|*.txt",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THEIR MIND

        # SAVE TO DISC
        outputfile = save_file_dialog.GetPath()
        np.savetxt(outputfile, list(zip((self.xp * 0.001), self.predicted_nt)), delimiter=' ', fmt='%.6f %.6f')

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

        # APPEND NEW DATA MENU TO 'MAG data MENU'
        self.mag_submenu = wx.Menu()
        self.m_obs_mag_submenu.Append(12000 + observed.id, observed.name, self.mag_submenu)

        # APPEND DELETE DATA OPTION TO THE NEW DATA MENU
        self.mag_submenu.Append(12000 + observed.id, 'delete observed data')

        # BIND TO DELETE OBSERVED GRAVITY FUNC
        self.Bind(wx.EVT_MENU, self.delete_obs_mag, id=12000 + observed.id)

        # INCREMENT OBSERVED GRAVITY COUNTER
        self.observed_magnetic_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.set_frame_limits()
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

        # INCREMENT MAG DERIV COUNTER
        self.observed_magnetic_counter += 1

        # UPDATE GMG GUI
        self.update_layer_data()
        self.draw()

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
        new_xy.id = self.xy_data_counter
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
            well.horizons_list[i - 2] = self.model_frame.plot(x, y, linestyle='-', linewidth='2', color='black')
            horizon_y_pos = well.data[i][1].astype(float)
            horizon = well.data[i][0].astype(str)

            # ALTERNATE POSITION OF ODDs/EVENs TO TRY AND AVOID OVERLAP
            if i % 2 == 0:
                horizon_x_pos = well.data[1][1].astype(float) - 1.05
                well.labels_list[i - 2] = self.model_frame.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                    xytext=(horizon_x_pos, horizon_y_pos),
                                                                    fontsize=well.text_size, weight='bold',
                                                                    horizontalalignment='left', verticalalignment='top',
                                                                    color='black', bbox=dict(boxstyle="round,pad=.4",
                                                                                             fc="0.8", ec='None'),
                                                                    clip_on=True)
            else:
                horizon_x_pos = well.data[1][1].astype(float) + 1.05
                well.labels_list[i - 2] = self.model_frame.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                                    xytext=(horizon_x_pos, horizon_y_pos),
                                                                    fontsize=well.text_size, weight='bold',
                                                                    horizontalalignment='right',
                                                                    verticalalignment='top',
                                                                    color='black', bbox=dict(boxstyle="round,pad=.4",
                                                                                             fc="0.8", ec='None'),
                                                                    clip_on=True)

        # APPEND WELL TO WELL DATA LIST
        self.well_data_list.append(well)

        # INCREAMENT WELL COUNTER
        self.well_counter += 1

        # UPDATE GMG
        self.update_layer_data()
        self.draw()

    def show_hide_well(self, event):
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
        text = list(zip(outcrop.data[:, 0].astype(float), outcrop.data[:, 1].astype(float),
                        outcrop.data[:, 3].astype(str)))

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

    def on_menu_set_button_press(self, event):
        """
        CHECK IF A NODE FROM THE CURRENT LAYER IS SELECTED; IF NOT, THEN SKIP THIS PART AND ONLY UPDATE ATTRIBUTES
        """

        # print()
        # self.index_node
        # print()
        # self.layer_list[self.currently_active_layer_id].type

        # GET NEW XY POINT
        new_x = float(self.x_input.GetValue())
        new_y = float(self.y_input.GetValue())

        # GET CURRENTLY ACTIVE LAYER NODES AND LABEL THEM xt AND yt
        xt = self.layer_list[self.currently_active_layer_id].x_nodes
        yt = self.layer_list[self.currently_active_layer_id].y_nodes

        # MODIFY xt and yt DEPENDING ON WHAT CONDITIONS ARE MET
        if self.layer_list[self.currently_active_layer_id].type == 'fixed' and self.index_node is not None:
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
        elif self.layer_list[self.currently_active_layer_id].type == 'floating' and self.index_node is not None:
            print("MOVING NODE ON FLOATING LAYER")
            xt[self.index_node] = new_x  # REPLACE OLD X WITH NEW X
            yt[self.index_node] = new_y  # REPLACE OLD Y WITH NEW Y

        # UPDATE THE CURRENTLY ACTIVE LAYER NODE LIST
        self.current_x_nodes = xt
        self.current_y_nodes = yt
        self.currently_active_layer.set_data(self.current_x_nodes, self.current_y_nodes)

        # UPDATE THE CURRENTLY ACTIVE LAYER layer_list ENTRY
        self.layer_list[self.currently_active_layer_id].x_nodes = xt
        self.layer_list[self.currently_active_layer_id].y_nodes = yt

        # COLOR CURRENTLY SELECTED NODE RED
        self.current_node.set_offsets([new_x, new_y])

        # UPDATE LAYER DATA
        self.set_density(self)
        self.set_susceptibility(self)
        self.set_angle_a(self)
        self.set_angle_b(self)
        self.set_angle_c(self)
        self.set_earth_field(self)

        # UPDATE GMG
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def get_node_under_point(self, event):
        """
        GET THE INDEX VALUE OF THE NODE UNDER POINT, AS LONG AS IT IS WITHIN NODE_CLICK_LIMIT TOLERANCE OF CLICK
        """
        # RESET NODE SWITCH
        self.didnt_get_node = False

        if self.pinch_switch is False:
            # PINCH MODE ISN'T ON, SO DO THE NORMAL ROUTINE

            # GET THE CURRENT LAYERS NODE VALUES
            xyt = self.currently_active_layer.get_xydata()
            xt = xyt[:, 0]
            yt = xyt[:, 1]

            # CALCULATE DISTANCES FROM THE MOUSE CLICK TO THE LAYER NODES
            d = np.sqrt((xt - event.xdata) ** 2 + (yt - event.ydata) ** 2)

            # FIND THE NODE INDEX VALUE FOR THE NODE CLOSEST TO THE CLICK
            self.index_arg = np.argmin(d)

            # CHECK IF THE NODE IS WITHIN THE "MEANT TO CLICK" DISTANCE
            if d[self.index_arg] >= self.node_click_limit:
                self.didnt_get_node = True
                return None, None
            else:
                # CHECK IF NODE IS A PINCHED POINT, IF YES FIND NODE OF ABOVE OR BELOW LAYER
                if self.layer_list[self.currently_active_layer_id].pinched is True:

                    # CREATE EMPTY LIST THE FILL WIH THE NODE INDEX VALUES
                    self.pinched_index_arg_list = []

                    for l in self.layer_list[self.currently_active_layer_id].pinched_list:
                        # GET LAYER NODES
                        xt = self.layer_list[l].x_nodes
                        yt = self.layer_list[l].y_nodes

                        # CALCULATE DISTANCES FROM THE MOUSE CLICK TO THE LAYER NODES
                        d = np.sqrt((xt - event.xdata) ** 2 + (yt - event.ydata) ** 2)

                        # FIND THE NODE INDEX VALUE FOR THE NODE CLOSEST TO THE CLICK AND APPEND IT
                        # TO THE PINCHED INDEX LIST
                        self.pinched_index_arg_list.append(np.argmin(d))
                else:
                    self.pinched_index_arg_list = None

                # NOW RETURN:
                # 1) THE INDEX VALUE OF THE NODE CLICKED.
                # 2) A LIST OF INDEX VALUES FOR NODES PINCHED TO INDEX_ARG (NONE IF THE NODE ISN'T PINCHED).
                return self.index_arg, self.pinched_index_arg_list
        else:
            # GMG IS IN PINCH MODE - SO JUST RETURN THE INDEX OF THE NODE

            # GET THE CURRENT LAYERS NODE VALUES
            xyt = self.currently_active_layer.get_xydata()
            xt = xyt[:, 0]
            yt = xyt[:, 1]

            # CALCULATE DISTANCES FROM THE MOUSE CLICK TO THE LAYER NODES
            d = np.sqrt((xt - event.xdata) ** 2 + (yt - event.ydata) ** 2)

            # FIND THE NODE INDEX VALUE FOR THE NODE CLOSEST TO THE CLICK
            self.index_arg = np.argmin(d)

            # CHECK IF THE NODE IS WITHIN THE "MEANT TO CLICK" DISTANCE
            if d[self.index_arg] >= self.node_click_limit:
                return None, None
            else:
                return self.index_arg, None

    def get_fault_node_under_point(self, event):
        """GET THE INDEX VALUE OF THE NODE UNDER POINT, AS LONG AS IT IS WITHIN NODE_CLICK_LIMIT TOLERANCE OF CLICK"""
        # RESET NODE SWITCH
        self.didnt_get_node = False

        # GET FAULT NODE XY DATA
        xy_data = self.currently_active_fault.get_xydata()
        x = xy_data[:, 0]
        y = xy_data[:, 1]

        # FIND NODE CLOSEST TO EVENT CLICK POINT
        d = np.sqrt((x - event.xdata) ** 2 + (y - event.ydata) ** 2)
        self.index_arg = np.argmin(d)

        # RETURN RESULTING NODE OR NONE
        if d[self.index_arg] >= self.node_click_limit:
            self.didnt_get_node = True
            return None
        else:
            return self.index_arg

    def button_press(self, event):
        """WHAT HAPPENS WHEN THE LEFT MOUSE BUTTON IS PRESSED"""
        # print("clicked")
        if event.inaxes is None:
            return  # CLICK IS OUTSIDE MODEL FRAME SO RETURN
        if event.button != 1:
            return

        if self.fault_picking_switch is False and self.capture is False and self.select_new_layer_nodes is False:
            # THEN GMG IS IN LAYER MODE
            # print("layer_mode")
            # GET THE NODE CLOSEST TO THE CLICK AND ANY PINCHED NODES
            self.index_node, self.pinched_index_arg_list = self.get_node_under_point(event)
            if self.index_node is None:
                return

            # CHECK if the 'p' KEY IS ON (I.E. IN PINCH MODE)
            if self.pinch_switch is True:
                # SET THE FIRST NODE CLICKED AS NODE 1
                if self.pinch_count == 0:
                    self.layer_getting_pinched = self.currently_active_layer_id
                    self.node_to_pinch_index = self.index_node
                    self.pinch_count += 1

                elif self.pinch_count == 1 and self.currently_active_layer_id == self.layer_getting_pinched:
                    # USER HASN'T CHANGED TO A DIFFERENT LAYER YET
                    return

                else:
                    # USER IS UP TO SECOND MOUSE CLICK. SET THE NODE TO BE PINCHED AS THE SECOND NODE CLICKED

                    # GET CURRENT LAYER NODES
                    xyt = self.currently_active_layer.get_xydata()
                    xt, yt = xyt[:, 0], xyt[:, 1]

                    # SET THE PINCHED NODE
                    self.layer_list[self.layer_getting_pinched].x_nodes[self.node_to_pinch_index] = xt[self.index_node]
                    self.layer_list[self.layer_getting_pinched].y_nodes[self.node_to_pinch_index] = yt[self.index_node]

                    # UPDATE LAYER LINE
                    self.layer_list[self.layer_getting_pinched].node_mpl_actor[0].set_visible(False)
                    self.layer_list[self.layer_getting_pinched].node_mpl_actor[0].remove()
                    self.layer_list[self.layer_getting_pinched].node_mpl_actor = self.model_frame.plot(
                        self.layer_list[self.layer_getting_pinched].x_nodes,
                        self.layer_list[self.layer_getting_pinched].y_nodes,
                        color='blue', linewidth=1.0, alpha=1.0)

                    # UPDATE LAYER POLYGON FILL
                    current_color = self.layer_list[self.layer_getting_pinched].polygon_mpl_actor[0].get_fc()
                    self.layer_list[self.layer_getting_pinched].polygon_mpl_actor[0].set_visible(False)
                    self.layer_list[self.layer_getting_pinched].polygon_mpl_actor[0].remove()

                    self.layer_list[self.layer_getting_pinched].polygon_mpl_actor = self.model_frame.fill(
                        self.layer_list[self.layer_getting_pinched].x_nodes,
                        self.layer_list[self.layer_getting_pinched].y_nodes, color=current_color, alpha=0.4,
                        closed=True, linewidth=None, ec=None)

                    # RESET PINCH COUNT
                    self.pinch_count = 0

                    # SET THE PINCH ATTRIBUTES IN THE LAYER OBJECTS
                    self.layer_list[self.layer_getting_pinched].pinched = True
                    if self.currently_active_layer_id in self.layer_list[self.layer_getting_pinched].pinched_list:
                        pass
                    else:
                        self.layer_list[self.layer_getting_pinched].pinched_list.append(self.currently_active_layer_id)

                    # SET THE PINCH ATTRIBUTES IN THE LAYER OBJECTS
                    self.layer_list[self.currently_active_layer_id].pinched = True
                    if self.layer_getting_pinched in self.layer_list[self.currently_active_layer_id].pinched_list:
                        pass
                    else:
                        self.layer_list[self.currently_active_layer_id].pinched_list.append(self.layer_getting_pinched)

                    # REDRAW MODEL
                    self.update_layer_data()
                    self.draw()
            else:
                # GMG IS IN SIMPLE LAYER MODE - SO JUST SET THE NEW NODE LOCATION
                xyt = self.currently_active_layer.get_xydata()
                xt, yt = xyt[:, 0], xyt[:, 1]
                self.x_input.SetValue(xt[self.index_node])
                self.y_input.SetValue(yt[self.index_node])

                # COLOR CURRENTLY SELECTED NODE RED
                self.current_node.set_offsets([xt[self.index_node], yt[self.index_node]])

        elif self.fault_picking_switch is True and self.select_new_layer_nodes is False \
                and self.select_new_fault_nodes is False:
            # THEN GMG IS IN FAULT MODE

            # GET CURRENT NODE
            self.selected_node = self.get_fault_node_under_point(event)
            if self.selected_node is None:
                return

            # GET CURRENT X AND Y COORDS
            xyt = self.currently_active_fault.get_xydata()
            self.xt = xyt[:, 0]
            self.yt = xyt[:, 1]

            # COLOR CURRENTLY SELECTED NODE RED
            self.current_node.set_offsets([self.xt[self.selected_node], self.yt[self.selected_node]])

        elif self.capture is True or self.select_new_layer_nodes is True or self.select_new_fault_nodes is True:
            # COORDINATE CAPTURE MODE OR NEW LAYER CREATION OR NEW FAULT CREATION IS ON. SO PASS
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
        if self.didnt_get_node is True:
            # NO NODE WAS SELECTED WHEN CLICKING
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
            else:
                self.xt[self.selected_node] = self.new_x  # REPLACE OLD X WITH NEW X
                self.yt[self.selected_node] = self.new_y  # REPLACE OLD Y WITH NEW Y

            # UPDATE THE FAULT LIST RECORDS
            self.fault_list[self.currently_active_fault_id].x_nodes = self.xt
            self.fault_list[self.currently_active_fault_id].y_nodes = self.yt

            # UPDATE THE FAULT MPL ACTOR
            self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_xdata(self.xt)
            self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_ydata(self.yt)

            # UPDATE THE CURRENT VIEW OF THE FAULT
            self.currently_active_fault.set_data(self.xt, self.yt)

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
            if self.layer_list[self.currently_active_layer_id].type == str('fixed'):
                # GET X AND Y VALUES
                x = event.xdata  # GET X OF NEW POINT
                y = event.ydata  # GET Y OF NEW POINT

                # GET CURRENT X AND Y ARRAYS
                xt = self.layer_list[self.currently_active_layer_id].x_nodes
                yt = self.layer_list[self.currently_active_layer_id].y_nodes
                current_x_value = xt[self.index_node]
                current_y_value = yt[self.index_node]

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
            elif self.layer_list[self.currently_active_layer_id].type == str('floating'):
                # GET THE X AND Y VALUES
                x = event.xdata  # GET X OF NEW POINT
                y = event.ydata  # GET Y OF NEW POINT
                xt = self.layer_list[self.currently_active_layer_id].x_nodes
                yt = self.layer_list[self.currently_active_layer_id].y_nodes
                current_x_value = xt[self.index_node]
                current_y_value = yt[self.index_node]

                xt[self.index_node] = x  # REPLACE OLD X WITH NEW X
                yt[self.index_node] = y  # REPLACE OLD Y WITH NEW Y

            # RESET THE LAYER WITH THE NEW NODE POSITION
            self.layer_list[self.currently_active_layer_id].x_nodes = xt
            self.layer_list[self.currently_active_layer_id].y_nodes = yt

            # DEAL WITH PINCHED NODE
            if self.layer_list[self.currently_active_layer_id].pinched is True:
                self.layer_list[self.currently_active_layer_id].x_nodes
                i = 0
                for l in self.layer_list[self.currently_active_layer_id].pinched_list:
                    if self.layer_list[l].x_nodes[self.pinched_index_arg_list[i]] == current_x_value and \
                            self.layer_list[l].y_nodes[self.pinched_index_arg_list[i]] == current_y_value:
                        # SET PINCHED NODE AS NEW VALUE
                        self.layer_list[l].x_nodes[self.pinched_index_arg_list[i]] = xt[self.index_node]
                        self.layer_list[l].y_nodes[self.pinched_index_arg_list[i]] = yt[self.index_node]
                    i += 1

            # SET THE CURRENTLY ACTIVE LAYER WITH THE NEW NODE
            self.current_x_nodes = xt
            self.current_y_nodes = yt
            self.currently_active_layer.set_data(self.current_x_nodes, self.current_y_nodes)

            # UPDATE "CURRENT NODE" RED DOT
            if xt[self.index_node] == self.x1:
                self.current_node.set_offsets([self.x1, y])
                # UPDATE NODE POSITION TEXT BOXES
            elif xt[self.index_node] == self.x2:
                self.current_node.set_offsets([self.x2, y])
            else:
                self.current_node.set_offsets([x, y])

            # UPDATE NODE POSITION TEXT BOXES
            self.x_input.SetValue(x)
            self.y_input.SetValue(y)

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
                self.new_layer_nodes = self.model_frame.plot(self.new_plotx, self.new_ploty, color='blue', marker='o')
                # FILL LAYER
                self.new_layer_fill = self.model_frame.fill(self.new_plotx, self.new_ploty, color='blue',
                                                            alpha=self.layer_transparency, closed=True, linewidth=None,
                                                            ec=None)
                # INCREMENT CLICK COUNTER
                self.click_count += 1

            elif self.click_count < 3:
                self.new_layer_nodes[0].set_xdata(self.new_plotx)
                self.new_layer_nodes[0].set_ydata(self.new_ploty)
                self.new_layer_fill[0].set_xy(list(zip(self.new_plotx, self.new_ploty)))

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

        # NEW FAULT CREATION SEQUENCE
        if self.select_new_fault_nodes is True:
            # APPEND NEW COORDINATES
            self.new_plotx.append(event.xdata)
            self.new_ploty.append(event.ydata)

            if self.click_count == 0:
                # PLOT NEW NODES
                self.new_fault_nodes = self.model_frame.plot(self.new_plotx, self.new_ploty, color='green', marker='o')
                # INCREMENT CLICK COUNTER
                self.click_count += 1

            elif self.click_count < 2:
                self.new_fault_nodes[0].set_xdata(self.new_plotx)
                self.new_fault_nodes[0].set_ydata(self.new_ploty)
                # INCREMENT CLICK COUNTER
                self.click_count += 1
            else:
                # REMOVE THE TEMP FAULT MPL ACTOR
                self.new_fault_nodes[0].set_visible(False)
                self.new_fault_nodes[0].remove()
                self.new_fault_nodes = None

                # SWITCH OFF NEW FAULT MODE
                self.select_new_fault_nodes = False
                self.click_count = 0

                # RUN FINAL PART OF FAULT LOADING
                self.create_new_fault()

        # RUN MODELLING ALGORITHMS
        if self.fault_picking_switch is False:
            self.run_algorithms()
        else:
            # UPDATE GMG GRAPHICS
            self.draw()

    def key_release(self, event):
        pass

    def key_press(self, event):

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
            self.current_x_nodes = np.insert(xt, [self.index_arg + 1], event.xdata)
            self.current_y_nodes = np.insert(yt, [self.index_arg + 1], event.ydata)

            # SET THE CURRENT LAYER WITH THE UPDATED NODE LIST
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
                self.current_x_nodes = [tup for i, tup in enumerate(xt) if i != self.index_arg]  # DELETE X
                self.current_y_nodes = [tup for i, tup in enumerate(yt) if i != self.index_arg]  # DELETE Y

            # SET THE CURRENT LAYER WITH THE UPDATED NODE LIST
            self.currently_active_layer.set_data(self.current_x_nodes, self.current_y_nodes)

            # # NOW CHECK FOR PINCHED NODES
            # index_arg2 = None
            # self.pinch_switch = False
            # self.pinched_index_arg_list = [None] * (self.total_layer_count + 1)  # CREATE LIST OF NONES = LENGTH AS NUMB OF LAYERS
            # for x in range(0, self.total_layer_count + 1):  # LOOP THROUGH ALL LAYERS TO CHECK FOR PINCHED NODES
            #     if x == self.currently_active_layer_id:
            #         pass
            #     x_node_list = self.layer_list[x].x_nodes
            #     y_node_list = self.layer_list[x].y_nodes
            #
            #     for i in range(0, len(x_node_list)):
            #         # NOW CHECK X AND Y VALUES ARE EQUAL
            #         if x_node_list[i] == xt[self.index_arg] and y_node_list[i] == yt[self.index_arg]:
            #             # IF ONE OF THE NODES FORM LIST IS EQUAL TO A NODE FROM THE OTHER LAYER THEN RETURN THE INDEX
            #             self.pinched_index_arg_list[x] = i
            #             self.pinch_switch = True
            #
            # # REMOVE PINCHED NODES
            # if self.pinch_switch is True:
            #     for k in range(len(self.pinched_index_arg_list)):
            #         if self.pinched_index_arg_list[k] is not None:
            #             next_x_list = self.plotx_list[k]
            #             next_y_list = self.ploty_list[k]
            #             # GET THE NODE LIST OF THE NEXT LAYER
            #             next_x_list = [tup for i, tup in enumerate(next_x_list) if
            #                            i != self.pinched_index_arg_list[k]]  # DELETE X
            #             next_y_list = [tup for i, tup in enumerate(next_y_list) if
            #                            i != self.pinched_index_arg_list[k]]  # DELETE Y
            #             # OVERWRITE THE NODE LIST WITH UPDATED LIST
            #             self.plotx_list[k], self.ploty_list[k] = next_x_list, next_y_list

            # SHIFT CURRENT NODE COLORING TO PREVIOUS NODE
            self.current_node.set_offsets([xt[self.index_arg - 1], yt[self.index_arg - 1]])

            # UPDATE LAYER DATA AND PLOT
            self.update_layer_data()
            self.run_algorithms()
            self.draw()

        # q = BEAT THE ZOOM BUG
        if event.key == 'q':
            self.nodes = True

        # n = CREATE NEW LAYER AT MOUSE POINT
        if event.key == 'n':
            self.new_layer(event)

        # < = MOVE TO NEXT LAYER
        if event.key == '.':
            if self.currently_active_layer_id == self.total_layer_count:
                # UPDATE LAYER DATA
                self.update_layer_data()

                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.currently_active_layer_id = 0

                # SET CURRENT INPUT VALUES IN MENU
                self.density_input.SetValue(self.layer_list[self.currently_active_layer_id].density)
                self.ref_density_input.SetValue(self.layer_list[self.currently_active_layer_id].reference_density)
                self.susceptibility_input.SetValue(self.layer_list[self.currently_active_layer_id].susceptibility)
                self.angle_a_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_a)
                self.angle_b_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_b)
                self.angle_c_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_c)
                self.earth_field_input.SetValue(self.layer_list[self.currently_active_layer_id].earth_field)

                # SET CURRENT NODE VALUES
                self.x_input.SetValue(self.layer_list[self.currently_active_layer_id].x_nodes[0])
                self.y_input.SetValue(self.layer_list[self.currently_active_layer_id].y_nodes[0])
                self.current_x_nodes = self.layer_list[self.currently_active_layer_id].x_nodes
                self.current_y_nodes = self.layer_list[self.currently_active_layer_id].y_nodes
                self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])

                # UPDATE MODEL
                self.update_layer_data()
                self.run_algorithms()
            else:
                # UPDATE LAYER DATA
                self.update_layer_data()

                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.currently_active_layer_id += 1

                # SET CURRENT INPUT VALUES IN MENU
                self.density_input.SetValue(self.layer_list[self.currently_active_layer_id].density)
                self.ref_density_input.SetValue(self.layer_list[self.currently_active_layer_id].reference_density)
                self.susceptibility_input.SetValue(self.layer_list[self.currently_active_layer_id].susceptibility)
                self.angle_a_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_a)
                self.angle_b_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_b)
                self.angle_c_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_c)
                self.earth_field_input.SetValue(self.layer_list[self.currently_active_layer_id].earth_field)

                # SET CURRENT NODE VALUES
                self.x_input.SetValue(self.layer_list[self.currently_active_layer_id].x_nodes[0])
                self.y_input.SetValue(self.layer_list[self.currently_active_layer_id].y_nodes[0])
                self.current_x_nodes = self.layer_list[self.currently_active_layer_id].x_nodes
                self.current_y_nodes = self.layer_list[self.currently_active_layer_id].y_nodes
                self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])

                # UPDATE MODEL
                self.update_layer_data()
                self.run_algorithms()

        # < = MOVE TO NEXT LAYER
        if event.key == ',':
            if self.currently_active_layer_id == 0:
                # UPDATE LAYER DATA
                self.update_layer_data()

                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.currently_active_layer_id = self.total_layer_count

                # SET CURRENT INPUT VALUES IN MENU
                self.density_input.SetValue(self.layer_list[self.currently_active_layer_id].density)
                self.ref_density_input.SetValue(self.layer_list[self.currently_active_layer_id].reference_density)
                self.susceptibility_input.SetValue(self.layer_list[self.currently_active_layer_id].susceptibility)
                self.angle_a_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_a)
                self.angle_b_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_b)
                self.angle_c_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_c)
                self.earth_field_input.SetValue(self.layer_list[self.currently_active_layer_id].earth_field)

                # SET CURRENT NODE VALUES
                self.x_input.SetValue(self.layer_list[self.currently_active_layer_id].x_nodes[0])
                self.y_input.SetValue(self.layer_list[self.currently_active_layer_id].y_nodes[0])
                self.current_x_nodes = self.layer_list[self.currently_active_layer_id].x_nodes
                self.current_y_nodes = self.layer_list[self.currently_active_layer_id].y_nodes
                self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])

                # UPDATE MODEL
                self.update_layer_data()
                self.run_algorithms()
            else:
                # UPDATE LAYER DATA
                self.update_layer_data()
                # MOVE TO LAYER ABOVE INPUT LAYER & ASSIGN CURRENT LAYER XY DATA
                self.currently_active_layer_id -= 1

                # SET CURRENT INPUT VALUES IN MENU
                self.density_input.SetValue(self.layer_list[self.currently_active_layer_id].density)
                self.ref_density_input.SetValue(self.layer_list[self.currently_active_layer_id].reference_density)
                self.susceptibility_input.SetValue(self.layer_list[self.currently_active_layer_id].susceptibility)
                self.angle_a_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_a)
                self.angle_b_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_b)
                self.angle_c_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_c)
                self.earth_field_input.SetValue(self.layer_list[self.currently_active_layer_id].earth_field)

                # SET CURRENT NODE VALUES
                self.x_input.SetValue(self.layer_list[self.currently_active_layer_id].x_nodes[0])
                self.y_input.SetValue(self.layer_list[self.currently_active_layer_id].y_nodes[0])
                self.current_x_nodes = self.layer_list[self.currently_active_layer_id].x_nodes
                self.current_y_nodes = self.layer_list[self.currently_active_layer_id].y_nodes
                self.current_node.set_offsets([self.current_x_nodes[0], self.current_y_nodes[0]])

                # UPDATE MODEL
                self.update_layer_data()
                self.run_algorithms()

            # UPDATE MODEL FRAME
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

        # if event.key == 'return':
        #     pass

    def pinch_out_layer(self, event):
        """PINCH OUT A FIXED LAYER OVER A GIVEN X RANGE"""
        if self.layer_list[self.currently_active_layer_id].type == 'floating':
            error_message = "Only fixed layers can use the bulk pinch function! \n" \
                            "Use the p key to pinch individual floating layer nodes"
            MessageDialog(self, -1, error_message, "Error")
            return
        else:
            # CREATE AND SHOW POP OUT BOX
            pinch_box = PinchDialog(self, -1, 'Pinch Out Layer:', self.layer_list, self.currently_active_layer_id)
            answer = pinch_box.ShowModal()

            # SET NEW NODE VALUES
            self.current_x_nodes = pinch_box.pinched_x
            self.current_y_nodes = pinch_box.pinched_y
            self.layer_list[self.currently_active_layer_id].x_nodes = pinch_box.pinched_x
            self.layer_list[self.currently_active_layer_id].y_nodes = pinch_box.pinched_y

            self.layer_getting_pinched = pinch_box.layer_getting_pinched

            # SET THE PINCH ATTRIBUTES IN THE LAYER OBJECTS
            self.layer_list[self.layer_getting_pinched].pinched = True
            if self.currently_active_layer_id in self.layer_list[self.layer_getting_pinched].pinched_list:
                pass
            else:
                self.layer_list[self.layer_getting_pinched].pinched_list.append(self.currently_active_layer_id)

            # SET THE PINCH ATTRIBUTES IN THE LAYER OBJECTS
            self.layer_list[self.currently_active_layer_id].pinched = True
            if self.layer_getting_pinched in self.layer_list[self.currently_active_layer_id].pinched_list:
                pass
            else:
                self.layer_list[self.currently_active_layer_id].pinched_list.append(self.layer_getting_pinched)

            # SET THE CURRENTLY SELECTED (RED) ACTIVE NODE
            self.current_node.set_offsets([self.layer_list[self.currently_active_layer_id].x_nodes[0],
                                           self.layer_list[self.currently_active_layer_id].y_nodes[0]])

            # REDRAW MODEL
            self.update_layer_data()
            self.draw()

    def depinch_layer(self, event):
        """PINCH OUT A FIXED LAYER OVER A GIVEN X RANGE"""
        # CREATE AND SHOW POP OUT BOX
        depinch_box = DepinchDialog(self, -1, 'Depinch layer', self.layer_list, self.currently_active_layer_id,
                                    self.total_layer_count)
        answer = depinch_box.ShowModal()

        # SET NEW NODE VALUES
        self.current_x_nodes = depinch_box.depinched_x
        self.current_y_nodes = depinch_box.depinched_y

        self.layer_list[self.currently_active_layer_id].x_nodes = depinch_box.depinched_x
        self.layer_list[self.currently_active_layer_id].y_nodes = depinch_box.depinched_y

        # CHECK IF THE PINCHED LAYER CONNECTIONS ATE STILL ACTIVE. IF NOT, THEN REMOVE THEM FROM THE LIST
        # SET THE CURRENT LAYER NODES (USE X+Y FOR COMPARISON)
        current_nodes = self.layer_list[self.currently_active_layer_id].x_nodes \
                        + self.layer_list[self.currently_active_layer_id].y_nodes

        # SET THE PINCHED LAYERS NODES (USE X+Y FOR COMPARISON)
        for l in self.layer_list[self.currently_active_layer_id].pinched_list:
            pinched_nodes = self.layer_list[l].x_nodes + self.layer_list[l].y_nodes

            # LOOP THROUGH NODES AND COMPARE
            still_pinch_connected = False
            for node in current_nodes:
                for node2 in pinched_nodes:
                    if node == node2:
                        still_pinch_connected = True
            # IF THE LAYERS ARE NO LONGER CONNECTED, THEN REMOVE THE CONNECTION IN THE PINCH LIST
            if still_pinch_connected == False:
                self.layer_list[self.currently_active_layer_id].pinched_list.remove(l)
                self.layer_list[l].pinched_list.remove(self.currently_active_layer_id)
            # IF THE LAYER NO LONGER HAS ANY CONNECTIONS, THEN SET PINCHED SWITCH AS FALSE
            if not self.layer_list[l].pinched_list:
                self.layer_list[l].pinched = False
            # IF THE LAYER NO LONGER HAS ANY CONNECTIONS, THEN SET PINCHED SWITCH AS FALSE
            if not self.layer_list[self.currently_active_layer_id].pinched_list:
                self.layer_list[self.currently_active_layer_id].pinched = False

        # SET THE CURRENTLY SELECTED (RED) ACTIVE NODE
        self.current_node.set_offsets([self.layer_list[self.currently_active_layer_id].x_nodes[0],
                                       self.layer_list[self.currently_active_layer_id].y_nodes[0]])

        # REDRAW MODEL
        self.update_layer_data()
        self.draw()

    def bulk_shift(self, event):
        """APPLY A BULK X AND/OR Y SHIFT TO A GIVEN LAYERS NODES"""

        # CREATE AND SHOW POP OUT BOX
        bulk_shift_box = BulkShiftDialog(self, -1, 'Layer bulk shift', self.layer_list, self.currently_active_layer_id)
        answer = bulk_shift_box.ShowModal()

        # SET NEW NODE VALUES
        self.current_x_nodes = bulk_shift_box.new_x
        self.current_y_nodes = bulk_shift_box.new_y

        self.layer_list[self.currently_active_layer_id].x_nodes = bulk_shift_box.new_x
        self.layer_list[self.currently_active_layer_id].y_nodes = bulk_shift_box.new_y

        # SET THE CURRENTLY SELECTED (RED) ACTIVE NODE
        self.current_node.set_offsets([self.layer_list[self.currently_active_layer_id].x_nodes[0],
                                       self.layer_list[self.currently_active_layer_id].y_nodes[0]])

        # REDRAW MODEL
        self.update_layer_data()
        self.draw()

    def new_layer(self, event):
        new_layer_dialogbox = NewLayerDialog(self, -1, 'Create New Layer', 'new')
        answer = new_layer_dialogbox.ShowModal()

        if new_layer_dialogbox.fixed:
            # CREATING A NEW FIXED LAYER

            # INCREMENT THE CURRENT LAYER INDEX VALUE (self.currently_active_layer_id)
            self.total_layer_count += 1

            # SET THE ACTIVE LAYER AS THE NEWLY CREATED LAYER
            self.currently_active_layer_id = self.total_layer_count

            # CREATE A NEW LAYER OBJECT
            new_layer = Layer()

            # SET SOME OF THE NEW LAYERS ATTRIBUTES
            new_layer.id = self.currently_active_layer_id
            new_layer.name = str('layer %s') % self.currently_active_layer_id
            new_layer.type = str('fixed')
            new_layer.include_in_calculations_switch = True

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
            self.currently_active_layer.set_xdata(new_layer.x_nodes)
            self.currently_active_layer.set_ydata(new_layer.y_nodes)
            self.currently_active_layer.set_color(new_layer.color)

            # SET CURRENTLY ACTIVE LAYER NODE OBJECTS
            self.current_x_nodes = new_layer.x_nodes
            self.current_y_nodes = new_layer.y_nodes

            # SET THE CURRENTLY SELECTED (RED) ACTIVE NODE
            self.current_node.set_offsets([new_layer.x_nodes[0], new_layer.y_nodes[0]])

            # SET CURRENT ATTRIBUTE INPUTS IN LEFT PANEL
            self.density_input.SetValue(new_layer.density)
            self.ref_density_input.SetValue(new_layer.reference_density)
            self.susceptibility_input.SetValue(new_layer.susceptibility)
            self.angle_a_input.SetValue(new_layer.angle_a)
            self.angle_b_input.SetValue(new_layer.angle_b)
            self.angle_c_input.SetValue(new_layer.angle_c)
            self.earth_field_input.SetValue(new_layer.earth_field)

            # APPEND NEW LAYER TO THE LAYER LIST
            self.layer_list.append(new_layer)

            # UPDATE GMG FRAME
            self.update_layer_data()
            self.draw()

        elif not new_layer_dialogbox.fixed:
            print("not")
            # CREATING A NEW FLOATING LAYER
            self.new_plotx = []
            self.new_ploty = []
            self.click_count = 0
            self.new_layer_nodes = None
            self.new_layer_fill = None

            # SWITCH ON MOUSE CLICK CAPTURE MODE TO CREATE NEW LAYER (SEE button_release FUNC FOR CONTINUATION OF CODE)
            self.select_new_layer_nodes = True

        else:
            # USER CHANGED THEIR MIND - NO NEW LAYER ADDED
            pass

    def create_new_floating_layer(self):
        """CREATE A NEW FLOATING LAYER USING FOUR USER INPUT MOUSE CLICKS"""
        # INCREMENT THE TOTAL LAYER COUNT
        self.total_layer_count += 1

        # SET CURRENTLY ACTIVE LAYER AS THE NEWLY CREATED LAYER
        self.currently_active_layer_id = self.total_layer_count

        # CREATE NEW LAYER OBJECT
        new_layer = Layer()

        # SOURCE NEW NODES FROM USER CLICKS
        new_layer.x_nodes = self.new_plotx
        new_layer.y_nodes = self.new_ploty

        # SET CURRENTLY ACTIVE LAYER NODE OBJECTS
        self.current_x_nodes = new_layer.x_nodes
        self.current_y_nodes = new_layer.y_nodes

        # SET SOME OF THE NEW LAYERS ATTRIBUTES
        new_layer.id = self.currently_active_layer_id
        new_layer.name = str('layer %s') % self.currently_active_layer_id
        new_layer.type = str('floating')
        new_layer.include_in_calculations_switch = True

        # ADD NEW LAYER TO THE LAYER TREE DISPLAY
        self.tree_items.append('layer %s' % (int(self.currently_active_layer_id)))
        self.item = 'layer %s' % (int(self.currently_active_layer_id))
        self.add_new_tree_nodes(self.root, self.item, self.currently_active_layer_id)

        # CREATE LAYER LINE
        new_layer.node_mpl_actor = self.model_frame.plot(new_layer.x_nodes, new_layer.y_nodes, color='blue',
                                                         linewidth=1.0, alpha=1.0)
        # CREATE LAYER POLYGON FILL
        new_layer.polygon_mpl_actor = self.model_frame.fill(new_layer.x_nodes, new_layer.y_nodes, color='blue',
                                                            alpha=self.layer_transparency, closed=True,
                                                            linewidth=None, ec=None)

        # SET THE CURRENTLY SELECTED (RED) ACTIVE NODE
        self.current_node.set_offsets([new_layer.x_nodes[0], new_layer.y_nodes[0]])

        # SET CURRENT ATTRIBUTE INPUTS IN LEFT PANEL
        self.density_input.SetValue(new_layer.density)
        self.ref_density_input.SetValue(new_layer.reference_density)
        self.susceptibility_input.SetValue(new_layer.susceptibility)
        self.angle_a_input.SetValue(new_layer.angle_a)
        self.angle_b_input.SetValue(new_layer.angle_b)
        self.angle_c_input.SetValue(new_layer.angle_c)
        self.earth_field_input.SetValue(new_layer.earth_field)

        # APPEND NEW LAYER TO THE LAYER LIST
        self.layer_list.append(new_layer)

        # UPDATE MODEL
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def load_layer(self, event):
        """LOAD A NEW FLOATING LAYER FROM A SPACE DELIMITED XY TEXT FILE"""

        # OPEN NEW LAYER DIALOG AND LET USER CHOSE BETWEEN FIXED AND FLOATING LAYER
        new_layer_dialogbox = NewLayerDialog(self, -1, 'Create New Layer', 'load')
        answer = new_layer_dialogbox.ShowModal()
        # IF THE USER CHANGES THERE MIND, THEN CLOSE THE DIALOG BOX AND RETURN TO MODELLNG
        if new_layer_dialogbox.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THERE MIND

        # 1. IF THE USER CHOOSES TO LOAD A NEW FIXED LAYER
        if new_layer_dialogbox.fixed:
            print("got to here")

            # BEGIN LOADING LAYER
            open_file_dialog = wx.FileDialog(self, "Open Layer", "", "", "Layer XY files (*.txt)|*.txt",
                                             wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
            answer = open_file_dialog.ShowModal()

            file_in = open_file_dialog.GetPath()
            new_layer_nodes = np.genfromtxt(file_in, autostrip=True, delimiter=' ', dtype=float)

            # 1.1 DETERMINE WHICH LAYER IS THE LAST PREVIOUS "FIXED LAYER"; SET THIS LAYER AS "previous_fixed_layer"
            if self.total_layer_count > 0:
                for i in range(0, self.total_layer_count):
                    if self.layer_list[i].type == 'fixed':
                        previous_fixed_layer = i
                    else:
                        continue
            else:
                previous_fixed_layer = 0

            # INCREMENT THE LAYER COUNT
            self.currently_active_layer_id = self.total_layer_count
            self.currently_active_layer_id += 1

            # INCREMENT THE TOTAL LAYER COUNT
            self.total_layer_count += 1

            # CREATE NEW LAYER OBJECT
            new_layer = Layer()

            # 1.2 SET NEW NODES FROM LOADED FILE
            # SET FIRST PADDING NODE
            new_layer.x_nodes = np.array(-(float(self.padding)) + new_layer_nodes[0, 0])
            new_layer.y_nodes = np.array(new_layer_nodes[0, 1])
            # APPEND LOADED NODES
            new_layer.x_nodes = np.append(new_layer.x_nodes, new_layer_nodes[:, 0])
            new_layer.y_nodes = np.append(new_layer.y_nodes, new_layer_nodes[:, 1])
            # SET LAST PADDING NODE
            new_layer.x_nodes = np.append(new_layer.x_nodes, (float(self.padding)+new_layer_nodes[0, 0]))
            new_layer.y_nodes = np.append(new_layer.y_nodes, new_layer_nodes[0, 1])

            new_layer.type = str('fixed')

        # 2. IF THE USER CHOOSES TO LOAD A NEW FLOATING LAYER
        elif not new_layer_dialogbox.fixed:

            # # BEGIN LOADING LAYER
            open_file_dialog = wx.FileDialog(self, "Open Layer", "", "", "Layer XY files (*.*)|*.*",
                                             wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
            answer = open_file_dialog.ShowModal()
            
            file_in = open_file_dialog.GetPath()
            new_layer_nodes = np.genfromtxt(file_in, autostrip=True, delimiter=' ', dtype=float)

            # INCREMENT THE LAYER COUNT
            self.currently_active_layer_id = self.total_layer_count
            self.currently_active_layer_id += 1

            # INCREMENT THE TOTAL LAYER COUNT
            self.total_layer_count += 1

            # CREATE NEW LAYER OBJECT
            new_layer = Layer()

            # 2.1 SET NEW NODES FROM LOADED FILE
            new_layer.x_nodes = new_layer_nodes[:, 0]
            new_layer.y_nodes = new_layer_nodes[:, 1]
            new_layer.type = str('floating')

        # SET CURRENTLY ACTIVE LAYER NODE OBJECTS
        self.current_x_nodes = new_layer.x_nodes
        self.current_y_nodes = new_layer.y_nodes

        # SET SOME OF THE NEW LAYERS ATTRIBUTES
        new_layer.id = self.currently_active_layer_id
        new_layer.name = str('layer %s') % self.currently_active_layer_id
        new_layer.include_in_calculations_switch = True

        # ADD NEW LAYER TO THE LAYER TREE DISPLAY
        self.tree_items.append('layer %s' % (int(self.currently_active_layer_id)))
        self.item = 'layer %s' % (int(self.currently_active_layer_id))
        self.add_new_tree_nodes(self.root, self.item, self.currently_active_layer_id)

        # CREATE LAYER LINE
        new_layer.node_mpl_actor = self.model_frame.plot(new_layer.x_nodes, new_layer.y_nodes, color='blue',
                                                         linewidth=1.0, alpha=1.0)
        # CREATE LAYER POLYGON FILL
        new_layer.polygon_mpl_actor = self.model_frame.fill(new_layer.x_nodes, new_layer.y_nodes, color='blue',
                                                            alpha=self.layer_transparency, closed=True,
                                                            linewidth=None, ec=None)

        # SET THE CURRENTLY SELECTED (RED) ACTIVE NODE
        self.current_node.set_offsets([new_layer.x_nodes[0], new_layer.y_nodes[0]])

        # SET CURRENT ATTRIBUTE INPUTS IN LEFT PANEL
        self.density_input.SetValue(new_layer.density)
        self.ref_density_input.SetValue(new_layer.reference_density)
        self.susceptibility_input.SetValue(new_layer.susceptibility)
        self.angle_a_input.SetValue(new_layer.angle_a)
        self.angle_b_input.SetValue(new_layer.angle_b)
        self.angle_c_input.SetValue(new_layer.angle_c)
        self.earth_field_input.SetValue(new_layer.earth_field)

        # APPEND NEW LAYER TO THE LAYER LIST
        self.layer_list.append(new_layer)

        # UPDATE MODEL
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    def delete_layer(self, event):
        """Delete LAYER DATA"""
        # CANNOT DELETE THE FIRST LAYER
        if self.total_layer_count == 1:
            msg = "Sorry - layer one cannot be deleted!"
            dlg = wx.MessageDialog(self, msg, "Warning", wx.OK | wx.ICON_INFORMATION)
            result = dlg.ShowModal()
            dlg.Destroy()
        else:
            # HIDE THE LAYER MPL ACTORS
            self.layer_list[self.currently_active_layer_id].node_mpl_actor[0].set_visible(False)
            self.layer_list[self.currently_active_layer_id].polygon_mpl_actor[0].set_visible(False)

            # DELETE THE LAYER OBJECT
            del self.layer_list[self.currently_active_layer_id]

            #  REMOVE THE LAYER TREE ITEMS
            del self.tree_items[self.currently_active_layer_id]
            layers = self.tree.GetRootItem().GetChildren()
            self.tree.Delete(layers[self.currently_active_layer_id - 1])

            # RESET TREE ITEM ID'S
            layers = self.tree.GetRootItem().GetChildren()
            for i in range(len(layers)):
                self.tree.SetPyData(layers[i], i + 1)

            # DECREASE LAYER ID BY 1 FOR EACH LAYER THAT COMES AFTER THE ONE BEING DELETED
            # (I.E. IF LAYER 3 IS DEL; THEN LAYER 4 ID BECOMES 3 etc)
            try:
                for i in range(self.currently_active_layer_id, len(self.layer_list)):
                    self.layer_list[i].id -= 1
            except IndexError:
                pass

            # INCREMENT THE TOTAL LAYER COUNT
            self.total_layer_count -= 1

            # SET CURRENTLY ACTIVE LAYER TO LAYER 1
            self.currently_active_layer_id = 1

            # SET OBJECTS WITH THE CHOSEN LAYER
            self.density_input.SetValue(0.001 * self.layer_list[self.currently_active_layer_id].density)
            self.ref_density_input.SetValue(0.001 * self.layer_list[self.currently_active_layer_id].reference_density)
            self.susceptibility_input.SetValue(self.layer_list[self.currently_active_layer_id].susceptibility)
            self.angle_a_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_a)
            self.angle_b_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_b)
            self.angle_c_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_c)
            self.earth_field_input.SetValue(self.layer_list[self.currently_active_layer_id].earth_field)

            # GET THE XY NODES FROM THE ACTIVE LAYER AND SET THE CURRENTLY ACTIVE NODES (I.E. MAKE THEM INTERACTIVE)
            self.current_x_nodes = self.layer_list[self.currently_active_layer_id].x_nodes
            self.current_y_nodes = self.layer_list[self.currently_active_layer_id].y_nodes

            # UPDATE MODEL
            self.draw()
            self.update_layer_data()

            # SET THE CURRENTLY SELECTED (RED) ACTIVE NODE
            self.current_node.set_offsets([self.layer_list[self.currently_active_layer_id].x_nodes[0],
                                           self.layer_list[self.currently_active_layer_id].y_nodes[0]])
            self.draw()

    def write_layers_xy(self, event):
        """OUTPUT ALL LAYERS XY DATA TO INDIVIDUAL TEXT FILES"""

        # CREATE OUTPUT FILE
        save_file_dialog = wx.FileDialog(self, "Save LAYER XY", "", "", "xy files (*.xy)|*.xy",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THERE MIND

        # OUTPUT FILE
        all_layers_output_file = save_file_dialog.GetPath()
        # THE OUTPUT DIRECTORY
        output_dir = os.path.dirname(all_layers_output_file)

        # NOW WRITE OUT THE DATA
        with open(all_layers_output_file, 'w') as f:
            try:

                # OPEN "ALL LAYERS" OUTPUT FILE
                out = csv.writer(f, delimiter=' ')

                # LOOP THROUGH THE LAYERS
                for i in range(1, self.total_layer_count + 1):
                    # DEFINE THE LAYER NODES (REMOVE FIXED LAYER PADDING NODES)
                    if self.layer_list[i].type == 'fixed':
                        data = [self.layer_list[i].x_nodes[1:-1], self.layer_list[i].y_nodes[1:-1]]
                        layer_write = list(zip(self.layer_list[i].x_nodes[1:-1], self.layer_list[i].y_nodes[1:-1]))
                    else:
                        data = [self.layer_list[i].x_nodes, self.layer_list[i].y_nodes]
                        layer_write = list(zip(self.layer_list[i].x_nodes, self.layer_list[i].y_nodes))

                    # OUTPUT THE LAYER NAME TO THE ALL LAYERS FILE
                    f.write(">" + self.loaded_tree_items[i] + "\n")

                    # WRITE LAYER TO "ALL LAYERS" FILE
                    out.writerows(list(zip(*data)))

                    # SAVE THE LAYER AS AN INDIVIDUAL FILE
                    np.savetxt(output_dir + '/' + self.loaded_tree_items[i] + '.xy', layer_write, delimiter=' ',
                               fmt='%f %f')

                # CLOSE THE "ALL LAYERS" OUTPUT FILE
                f.close()
            except IndexError:
                f.close()
                pass

    def write_c_xy(self, event):
        # CREATE OUTPUT FILE
        save_file_dialog = wx.FileDialog(self, "Save c.in RayInvr file", "", "", "in files (*.in)|*.in",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED THEIR MIND

        # NOW WRITE OUT THE DATA
        output_stream = save_file_dialog.GetPath()
        with open(output_stream, 'wb') as f:
            # LAYER NODES
            for i in range(0, self.total_layer_count + 1):
                f.write('B  {0}\n'.format(i + 1))

                # DEFINE THE LAYER NODES (REMOVE FIXED LAYER PADDING NODES)
                if self.layer_list[i].type == 'fixed':
                    x_nodes = self.layer_list[i].x_nodes[1:-1]
                    y_nodes = self.layer_list[i].y_nodes[1:-1]
                else:
                    x_nodes = self.layer_list[i].x_nodes
                    y_nodes = self.layer_list[i].y_nodes

                # SAVE FILE
                data = list(zip(x_nodes, y_nodes, np.ones(len(y_nodes))))
                np.savetxt(f, data, delimiter=' ', fmt='%6.02f %3.02f %1d')

            # VELOCITY NODES
            for i in range(0, self.total_layer_count):
                density = self.layer_list[i].density

                # CONVERT DENSITY TO VELOCITY USING GARNERS RULE
                velocity = round((m.pow((density / 1670.), (1. / .25))), 2)

                # CONVERT DENSITY TO VELOCITY USING NAFE-DRAKE EQUATION
                # velocity = (1.6612*density) - (0.4721*density)**2 + (0.0671*density)**3 -
                # (0.0043*density)**4 + (0.000106*density)**5

                # DEFINE THE LAYER NODES (REMOVE FIXED LAYER PADDING NODES)
                if self.layer_list[i].type == 'fixed':
                    x_nodes = self.layer_list[i].x_nodes[1:-1]
                    y_nodes = self.layer_list[i].y_nodes[1:-1]
                else:
                    x_nodes = self.layer_list[i].x_nodes
                    y_nodes = self.layer_list[i].y_nodes

                # FORMAT c.in FILE
                f.write('B  {0}\n'.format(i))
                data = list(zip(x_nodes, np.linspace(velocity, velocity, len(y_nodes)), np.ones(len(x_nodes)),
                                np.linspace(velocity, velocity, len(x_nodes)), np.ones(len(x_nodes))))

                # OUTPUT FILE
                np.savetxt(f, data, delimiter=' ', fmt='%6.02f %3.02f %1d %3.02f %1d')

    def capture_coordinates(self, event):
        if self.capture is False:
            self.capture = True

            # CREATE INSTANCE OF CAPTURE COORDINATES
            self.capture_window = CaptureCoordinates(self, -1, 'Capture Coordinates')
            self.capture_window.Show(True)

    # FAULT MODE CONTROLS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def toogle_fault_mode(self, event):
        """SWITCH FAULT PICKING MODE ON AND OFF"""
        if self.fault_picking_switch is True:
            self.fault_picking_switch = False

        elif self.fault_picking_switch is False:
            self.fault_picking_switch = True

    def pick_new_fault(self, event):
        """FAULT PICKING/LINE DRAWING MODE"""

        # CHECK IF FAULT PICKING MODE IS ON
        if self.fault_picking_switch is False:
            MessageDialog(self, -1, "Faulting picking mode is not activated.\nTurn on fault picking mode first.",
                          "Fault picker")
        else:
            # PROMPT NEW FAULT DIALOG BOX
            if self.total_fault_count == 0:
                # CREATE NEW CURRENT FAULT GRAPHIC
                self.currently_active_fault, = self.model_frame.plot([-100000, -100000], [-100000, -100000], marker='s',
                                                                     color='green', linewidth=0.75, alpha=1.0, zorder=2,
                                                                     picker=True)
            # START CREATE NEW FAULT PROCESS
            MessageDialog(self, -1, "Select three nodes to create a new fault", "Select New Fault")
            self.select_new_fault_nodes = True
            self.new_plotx = []
            self.new_ploty = []
            self.click_count = 0
            self.new_layer_nodes = None
            # NOW WAIT FOR USER NODE CLICKS - ON THE THIRD CLICK THE button_release FUNC
            # WILL ACTIVATE create_new_fault()

    def create_new_fault(self):
        """CREATE A NEW FAULT USING THREE USER INPUT MOUSE CLICKS"""

        # SET CURRENTLY ACTIVE LAYER AS THE NEWLY CREATED LAYER
        self.currently_active_fault_id = self.total_fault_count

        # CREATE NEW LAYER OBJECT
        new_fault = Fault()

        # SOURCE NEW NODES FROM USER CLICKS
        new_fault.id = self.currently_active_fault_id
        new_fault.name = str('Fault')
        new_fault.x_nodes = self.new_plotx
        new_fault.y_nodes = self.new_ploty

        # SET CURRENTLY ACTIVE LAYER NODE OBJECTS
        self.current_x_nodes = new_fault.x_nodes
        self.current_y_nodes = new_fault.y_nodes

        # SET SOME OF THE NEW LAYERS ATTRIBUTES
        new_fault.id = self.currently_active_layer_id
        new_fault.name = str('Fault %s') % self.currently_active_fault_id
        # CREATE LAYER LINE
        new_fault.mpl_actor = self.model_frame.plot(new_fault.x_nodes, new_fault.y_nodes, color='green',
                                                    marker='o', linewidth=0.5, zorder=1, alpha=1.0)

        # APPEND THE NEW FAULT TO THE FAULT TREE SIDE PANEL USING add_new_tree_nodes FUNC
        # LIST OF FAULT NAMES
        self.fault_tree_items.append('fault %s' % (int(self.currently_active_fault_id)))
        self.fault_item = 'fault %s' % (int(self.currently_active_fault_id))
        self.add_new_tree_nodes(self.fault_tree_root, self.fault_item, self.currently_active_fault_id)

        self.fault_tree.SetSpacing(40)

        # self.fold_panel_three.Collapse()
        # self.fold_panel_three.Expand()

        # UPDATE CURRENT PLOT GRAPHICS
        self.currently_active_fault.set_data(new_fault.x_nodes, new_fault.y_nodes)

        # UPDATE CURRENT NODE RED DOT GRAPHIC
        self.current_node.set_offsets([new_fault.x_nodes[0], new_fault.y_nodes[0]])

        # APPEND NEW LAYER TO THE LAYER LIST
        self.fault_list.append(new_fault)

        # INCREMENT THE TOTAL LAYER COUNT
        self.total_fault_count += 1

        # UPDATE MODEL
        self.draw()

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
            self.fault_list[self.currently_active_fault_id].x_nodes = self.xt
            self.fault_list[self.currently_active_fault_id].y_nodes = self.yt

            # UPDATE FAULT GRAPHICS
            self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_xdata(self.xt)
            self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_ydata(self.yt)

            # UPDATE CURRENT FAULT OVERLAY GRAPHIC
            self.currently_active_fault.set_data(self.xt, self.yt)

        'd = DELETE NODE AT MOUSE POSITION'
        if event.key == 'd':
            if event.inaxes is None:
                return
            # FIND NODE CLOSEST TO CURSOR LOCATION
            d = np.sqrt((self.xt - event.xdata) ** 2 + (self.yt - event.ydata) ** 2)
            self.index_arg = np.argmin(d)
            self.distance = d[self.index_arg]

            if self.index_arg == 0 or \
                    self.index_arg == (len(self.fault_list[self.currently_active_fault_id].x_nodes) - 1):
                # PREVENT END NODES BEING DELETED
                return 0
            if self.distance >= self.node_click_limit:
                # CLICK WAS TO FAR AWAY FROM A NODE TO DELETE IT
                return 0
            else:
                # DELETE NODE BY RECREATING XY DATA WITHOUT CURRENT NODE
                self.xt = [tup for i, tup in enumerate(self.xt) if i != self.index_arg]  # DELETE X

                self.yt = [tup for i, tup in enumerate(self.yt) if i != self.index_arg]  # DELETE Y

                # UPDATE THE FAULT LIST RECORDS
                self.fault_list[self.currently_active_fault_id].x_nodes = self.xt
                self.fault_list[self.currently_active_fault_id].y_nodes = self.yt

                # UPDATE FAULT GRAPHICS
                self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_xdata(self.xt)
                self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_ydata(self.yt)

                # UPDATE CURRENT FAULT OVERLAY GRAPHIC
                self.currently_active_fault.set_data(self.xt, self.yt)

                # RESET CURRENT NOT POSITION TO FIRST NODE
                self.current_node.set_offsets([self.xt[0], self.yt[0]])

            # UPDATE GMG
            self.update_layer_data()

        '< = INCREMENT WHICH FAULT IS BEING EDITED'
        if event.key == ',':
            if self.currently_active_fault_id <= self.total_fault_count - 1 and self.currently_active_fault_id > 0:
                # INCREMENT TO NEXT FAULT
                self.currently_active_fault_id -= 1
            else:
                # GO TO NEWEST FAULT
                self.currently_active_fault_id = self.total_fault_count - 1

            # UPDATE CURRENT PLOT GRAPHICS
            self.currently_active_fault.set_data(self.fault_list[self.currently_active_fault_id].x_nodes,
                                                 self.fault_list[self.currently_active_fault_id].y_nodes)

            self.xt = self.fault_list[self.currently_active_fault_id].x_nodes
            self.yt = self.fault_list[self.currently_active_fault_id].y_nodes

        '> = INCREMENT WHICH FAULT IS BEING EDITED'
        if event.key == '.':
            if self.currently_active_fault_id < self.total_fault_count - 1:
                # INCREMENT TO NEXT FAULT
                self.currently_active_fault_id += 1
            elif self.currently_active_fault_id == self.total_fault_count - 1:
                # GO BACK TO FIRST FAULT
                self.currently_active_fault_id = 0

            # UPDATE CURRENT PLOT GRAPHICS
            self.currently_active_fault.set_data(self.fault_list[self.currently_active_fault_id].x_nodes,
                                                 self.fault_list[self.currently_active_fault_id].y_nodes)

            self.xt = self.fault_list[self.currently_active_fault_id].x_nodes
            self.yt = self.fault_list[self.currently_active_fault_id].y_nodes

        # UPDATE GMG
        self.update_layer_data()

    def on_fault_activated(self, event):
        """RESPONSE WHEN A FAULT NAME IS SELECTED"""

        # GET THE SELECTED FAULT INDEX NUMBER
        self.currently_active_fault_id = self.fault_tree.GetPyData(event.GetItem())

        if self.fault_picking_switch is False:
            self.fault_picking_switch = True

        # SET CHECKBOX AS CHECKED
        self.fault_tree.GetSelection().Check(checked=True)

        # UPDATE CURRENT PLOT GRAPHICS
        self.currently_active_fault.set_data(self.fault_list[self.currently_active_fault_id].x_nodes,
                                             self.fault_list[self.currently_active_fault_id].y_nodes)
        self.xt = self.fault_list[self.currently_active_fault_id].x_nodes
        self.yt = self.fault_list[self.currently_active_fault_id].y_nodes

        # UPDATE GRAPHICS WITH CURRENT FAULT SELECTED
        self.update_layer_data()

    def fault_checked(self, event):
        """TOGGLE WHETHER OR NOT A FAULT WILL BE PLOTTED IN THE MODEL FIGURE"""
        i = self.fault_tree.GetPyData(event.GetItem())

        if self.faults[i][0].get_visible() == True:
            # HIDE FAULT
            self.faults[i][0].set_visible(False)
            self.currently_active_fault.set_visible(False)
        else:
            # SHOW FAULT
            self.faults[i][0].set_visible(True)
            self.currently_active_fault.set_visible(True)

        # UPDATE FIGURE
        self.draw()

    def on_fault_tree_right_click_down(self, event):
        """WHEN A FAULT IN THE FAULT TREE MENU IS RIGHT CLICKED"""

        # FIRST RUN on_fault_activated
        self.on_fault_activated(event)

        # CREATE POPOUT MENU WITH OPTIONS AND BIND OPTIONS TO ACTIONS
        menu = wx.Menu()

        item1 = menu.Append(wx.ID_ANY, "Change fault colour")
        item2 = menu.Append(wx.ID_ANY, "Rename fault")

        self.Bind(wx.EVT_MENU, self.change_color, item1)
        self.Bind(wx.EVT_MENU, self.rename_fault, item2)

        self.PopupMenu(menu)
        menu.Destroy()

    def rename_fault(self, event):
        """USE A POPUP MENU TO RENAME THE FAULT"""

        # CREATE POP OUT MENU AND SHOW
        fault_name_box = LayerNameDialog(self, -1, 'Rename fault',
                                         self.fault_tree_items[self.currently_active_fault_id])
        new = fault_name_box.ShowModal()

        # WAIT FOR USER TO CLOSE POP OUT

        # GET THE NEW LAYER NAME FROM POP OUT
        new_name = fault_name_box.name

        # SET THE TREE AND LAYER OBJECT WITH THE NEW NAME
        current_tree_items = self.fault_tree.GetRootItem().GetChildren()
        self.fault_tree.SetItemText(current_tree_items[self.currently_active_fault_id], str(new_name))
        self.fault_tree_items[self.currently_active_fault_id] = str(new_name)
        self.fault_list[self.currently_active_fault_id].name = str(new_name)

    # LAYER AND MODEL ATTRIBUTE CONTROLS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_density(self, value):
        self.layer_list[self.currently_active_layer_id].density = float(self.density_input.GetValue() * 1000.)

    def set_reference_density(self, value):
        self.layer_list[self.currently_active_layer_id].reference_density = \
            float(self.ref_density_input.GetValue() * 1000.)

    def set_background_density(self, event):
        grav_box = SetBackgroundDensityDialog(self, -1, 'Set background density')
        answer = grav_box.ShowModal()
        self.background_density_upper = float(grav_box.background_density_upper)

        for i in range(0, self.total_layer_count):
            self.layer_list[i].reference_density = float(self.background_density_upper) * 1000.
        # self.background_density_upper = float((grav_box.background_density_lower))
        # self.background_density_upper = float((grav_box.background_density_lid))
        self.absolute_densities = True
        self.draw()

    def set_susceptibility(self, value):
        self.layer_list[self.currently_active_layer_id].susceptibility = float(self.susceptibility_input.GetValue())

    def set_angle_a(self, value):
        self.layer_list[self.currently_active_layer_id].angle_a = float(self.angle_a_input.GetValue())

    def set_angle_b(self, value):
        self.layer_list[self.currently_active_layer_id].angle_b = float(self.angle_b_input.GetValue())

    def set_angle_c(self, value):
        self.layer_list[self.currently_active_layer_id].angle_c = float(self.angle_c_input.GetValue())

    def set_earth_field(self, value):
        self.layer_list[self.currently_active_layer_id].earth_field = float(self.earth_field_input.GetValue())

    def set_text_size(self, value):
        """GET NEW TEXT SIZE"""

        self.textsize = float(self.text_size_input.GetValue())

        # WELL DATA
        # LOOP THROUGH ALL WELL NAMES
        for i in range(len(self.well_data_list)):
            if self.well_data_list[i] is not None:
                self.well_data_list[i].text_size = self.textsize
                self.well_data_list[i].mpl_actor_name.set_size(self.textsize)
                # LOOP THROUGH ALL WELL HORIZON LABELS
                for l in range(len(self.well_data_list[i].labels_list)):
                    if self.well_data_list[i] is not None:
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

            if self.observed_gravity_list[i] is not None:
                if self.observed_gravity_list[i].name == selection.obs_name:
                    self.obs_gravity_data_for_rms = self.observed_gravity_list[i].data

    def set_obs_mag_rms(self, value):
        """SET THE DATA TO BE USED FOR CALCULATING THE RMS MISTFIT"""
        selection = SetObsRmsDialog(self, -1, 'Set RMS Input', self.observed_magnetic_list)
        answer = selection.ShowModal()
        for i in range(0, len(self.observed_magnetic_list)):
            if self.observed_magnetic_list[i] is not None:
                if self.observed_magnetic_list[i].name == selection.obs_name:
                    self.obs_mag_data_for_rms = self.observed_magnetic_list[i].data

    def model_rms(self, xp):
        """CALCULATE RMS MISFIT OF OBSERVED VS CALCULATED"""
        # GRAVITY RMS
        if self.obs_gravity_data_for_rms != [] and self.calc_grav_switch is True:
            x = xp * 0.001
            y = self.predicted_gravity
            self.gravity_rms_value, self.grav_residuals = model_stats.rms(self.obs_gravity_data_for_rms[:, 0],
                                                                          self.obs_gravity_data_for_rms[:, 1], x, y)
        else:
            pass
        # MAGNETICS RMS
        if self.obs_mag_data_for_rms != [] and self.calc_mag_switch is True:
            x = self.xp * 0.001
            y = self.predicted_nt
            self.magnetic_rms_value, self.mag_residuals = model_stats.rms(self.obs_mag_data_for_rms[:, 0],
                                                                          self.obs_mag_data_for_rms[:, 1], x, y)
        else:
            pass

    # LAYER ATTRIBUTE TABLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def open_attribute_table(self, event):
        self.attribute_table = AttributeEditor(self, -1, 'Attribute editor', self.tree_items, self.layer_list)
        self.attribute_table.Show(True)

    def attribute_set(self, new_tree_items, new_layer_list):
        """UPDATE GMG ATTRIBUTES WITH NEW ATTRIBUTES FROM THE ATTRIBUTE TABLE"""

        # UPDATE MAIN FRAME TREE ITEMS (RENAME THE ITEMS)
        current_tree_items = self.tree.GetRootItem().GetChildren()
        for i in range(0, len(self.tree_items) - 1):
            new_label = new_tree_items[i]
            self.tree.SetItemText(current_tree_items[i], new_tree_items[i + 1])

        # UPDATE MAIN FRAME ATTRIBUTES
        for l in range(0, len(self.layer_list)):
            self.layer_list[l].density = new_layer_list[l].density
            self.layer_list[l].reference_density = new_layer_list[l].reference_density
            self.layer_list[l].susceptibility = new_layer_list[l].susceptibility
            self.layer_list[l].angle_a = new_layer_list[l].angle_a
            self.layer_list[l].angle_b = new_layer_list[l].angle_b
            self.layer_list[l].angle_c = new_layer_list[l].angle_c
            self.layer_list[l].earth_field = new_layer_list[l].earth_field
            self.layer_list[l].color = new_layer_list[l].color

        # UPDATE LAYER ATTRIBUTE INPUTS
        self.density_input.SetValue(0.001 * self.layer_list[self.currently_active_layer_id].density)
        self.ref_density_input.SetValue(0.001 * self.layer_list[self.currently_active_layer_id].reference_density)
        self.susceptibility_input.SetValue(self.layer_list[self.currently_active_layer_id].susceptibility)
        self.angle_a_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_a)
        self.angle_b_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_b)
        self.angle_c_input.SetValue(self.layer_list[self.currently_active_layer_id].angle_c)
        self.earth_field_input.SetValue(self.layer_list[self.currently_active_layer_id].earth_field)

        # UPDATE GMG STATE
        self.update_layer_data()
        self.run_algorithms()
        self.draw()

    # LIVE GRAPHICS UPDATES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

            # UPDATE FAULT NODES
            self.fault_list[self.currently_active_fault_id].x_nodes = self.xt
            self.fault_list[self.currently_active_fault_id].y_nodes = self.yt

            # UPDATE FAULT MPL ACTOR
            self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_xdata(self.xt)
            self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_ydata(self.yt)
            self.fault_list[self.currently_active_fault_id].mpl_actor[0].set_color(
                self.fault_list[self.currently_active_fault_id].color)

            # UPDATE CURRENT PLOT GRAPHICS
            self.currently_active_fault.set_data(self.xt, self.yt)
            self.currently_active_fault.set_color(self.fault_list[self.currently_active_fault_id].color)
        else:
            # GMG IS IN LAYER MODE
            # UPDATE PLOT LISTS WITH LATEST EDIT
            self.layer_list[self.currently_active_layer_id].x_nodes = self.current_x_nodes
            self.layer_list[self.currently_active_layer_id].y_nodes = self.current_y_nodes

            # CREATE UPDATED POLYGON XYs -------------------------------------------------------------------------------
            # FIRST CREATE THE POLYLINE DATA (THE BOTTOM LINE OF THE LAYER POLYGON - (THIS DONE FIRST SO THE WHOLE
            # POLYGON ISN'T PASSED TO SELF.POLYPLOTS)

            for i in range(0, self.total_layer_count + 1):
                # CREATE THE LAYER POLYGONS TO PASS TO SELF.POLYGONS AND ONTO THE GRAV/MAG ALGORITHMS
                # FIRST SET UP XY DATA; IF LAYER IS BELOW LAYER 0 THEN ATTACH THE ABOVE LAYER TO COMPLETE THE POLYGON;
                # ELSE USE TOP LAYER CHECK FOR 'FIXED' LAYER MODE AND FIND LAST LAYER TO MAKE POLYGON

                if i >= 1 and self.layer_list[i].type == 'fixed':
                    # CHECK FOR LAST PREVIOUS FIXED LAYER AND USE ITS BASE TO COMPLETE THE POLYGON
                    for layer in range(i, 0, -1):
                        if self.layer_list[layer - 1].type == 'fixed':
                            # ASSIGN THE LAST FIXED LAYER INDEX
                            last_layer_index = layer - 1

                            # NOW APPEND NODES FOR BOUNDARY CONDITIONS (CONTINUOUS SLAB)
                            plotx = np.array(self.layer_list[i].x_nodes)
                            ploty = np.array(self.layer_list[i].y_nodes)

                            # SET THE PADDING NODES TO THE SAME DEPTH AS THE MODEL LIMIT NODES TO CREATE FLAT SLAB
                            ploty[0] = ploty[1]
                            ploty[-1] = ploty[-2]
                            self.layer_list[i].x_nodes = plotx
                            self.layer_list[i].y_nodes = ploty

                            # ADD NODES FROM ABOVE LAYER TO COMPETE POLYGON
                            layer_above_x = np.array(self.layer_list[last_layer_index].x_nodes)[::-1]
                            layer_above_y = np.array(self.layer_list[last_layer_index].y_nodes)[::-1]

                            polygon_x = np.append(np.array(layer_above_x), np.array(plotx))
                            polygon_y = np.append(np.array(layer_above_y), np.array(ploty))

                            # UPDATE LAYER POLYGON ATTRIBUTE
                            self.layer_list[i].polygon = list(zip(polygon_x, polygon_y))
                            break
                        else:
                            continue
                else:
                    # IF THE LAYER IS A SIMPLE 'FLOATING LAYER'
                    polygon_x = np.array(self.layer_list[i].x_nodes)
                    polygon_y = np.array(self.layer_list[i].y_nodes)

                    # UPDATE LAYER POLYGON ATTRIBUTE
                    self.layer_list[i].polygon = list(zip(polygon_x, polygon_y))
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

                # UPDATE POLYGON XY AND COLOR FILL
                self.layer_list[i].polygon_mpl_actor[0].set_xy(self.layer_list[i].polygon)
                self.layer_list[i].polygon_mpl_actor[0].set_color(next_color)

                # UPDATE LAYER LINES
                self.layer_list[i].node_mpl_actor[0].set_xdata(self.layer_list[i].x_nodes)
                self.layer_list[i].node_mpl_actor[0].set_ydata(self.layer_list[i].y_nodes)
                self.layer_list[i].node_mpl_actor[0].set_color(self.layer_list[i].color)
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
        """RUN POTENTIAL FIELD CALCULATION ALGORITHMS"""
        # --------------------------------------------------------------------------------------------------------------
        # :FUTURE: CALCULATE PREDICTED TOPOGRAPHY FROM ISOSTATIC FUNC
        self.pred_topo = np.zeros_like(self.xp)
        # --------------------------------------------------------------------------------------------------------------
        # tree.GetRootItem().GetChildren()[i].GetValue()
        # --------------------------------------------------------------------------------------------------------------

        # CALCULATE GRAVITY
        polygons_to_use = []
        densities_to_use = []
        if self.calc_grav_switch is True:
            # SELECT ONLY THOSE LAYERS THAT ARE CHECKED
            for layer in range(0, self.total_layer_count + 1):
                if self.layer_list[layer].include_in_calculations_switch is True:
                    # CHOSE POLYGONS
                    polygons_to_use.append(self.layer_list[layer].polygon)
                    # DETERMINE DENSITY CONTRASTS
                    densities_to_use.append((self.layer_list[layer].density -
                                             self.layer_list[layer].reference_density))

            # PASS POLYGONS TO BOTT ALGORITHM AND RETURN THE PREDICTED VALUES
            bott_input_polygons = []
            for p, d in zip(polygons_to_use, densities_to_use):
                bott_input_polygons.append(Polygon(1000 * np.array(p), {'density': d}))

            # SET THE PREDICTED VALUES AS THE BOTT OUTPUT
            # NB: NODES ARE INPUT LEFT TO RIGHT SO WE MUST MULTIPLY BY -1 TO PRODUCE THE CORRECT SIGN AT OUTPUT
            self.predicted_gravity = bott.gz(self.xp, self.gravity_observation_elv, bott_input_polygons) * -1
            # self.predicted_gravity = kim_and_wessel.gz(self.xp, self.gravity_observation_elv, bott_input_polygons)
        else:
            # SET THE PREDICTED VALUES AS ZEROS
            self.predicted_gravity = np.zeros_like(self.xp)

        # SET THE PREDICTED PLOT LINE WITH THE NEWLY CALCULATED VALUES
        self.pred_gravity_plot.set_data(self.xp * 0.001, self.predicted_gravity)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # CALCULATE MAGNETICS
        # ZIP POLYGONS WITH SUSCEPTIBILITIES AND PASS TO TALWANI AND HEIRTZLER ALGORITHM
        if self.calc_mag_switch is True:
            # SELECT ONLY THOSE LAYERS THAT ARE CHECKED
            polygons_to_use = []
            susceptibilities_to_use = []
            angle_a_to_use = []
            angle_b_to_use = []
            angle_c_to_use = []
            earth_field_to_use = []
            for layer in range(0, self.total_layer_count + 1):
                if self.layer_list[layer].include_in_calculations_switch is True:
                    polygons_to_use.append(self.layer_list[layer].polygon)
                    susceptibilities_to_use.append(self.layer_list[layer].susceptibility)
                    angle_a_to_use.append(self.layer_list[layer].angle_a)
                    angle_b_to_use.append(self.layer_list[layer].angle_b)
                    angle_c_to_use.append(self.layer_list[layer].angle_c)
                    earth_field_to_use.append(self.layer_list[layer].earth_field)

            # PASS TO TALWANI & HEIRTZLER ALGORITHM
            mag_input_polygons = []
            for p, s, a, b, c, f, in zip(polygons_to_use, susceptibilities_to_use, angle_a_to_use, angle_b_to_use,
                                         angle_c_to_use, earth_field_to_use):
                mag_input_polygons.append(Polygon(1000. * np.array(p), {'susceptibility': s, 'angle_a': a,
                                                                        'angle_b': b, 'angle_c': c, 'f': f}))
            # SET THE PREDICTED VALUES AS THE TALWANI & HEIRTZLER OUTPUT
            # NB: NODES ARE INPUT LEFT TO RIGHT SO WE MUST MULTIPLY BY -1 TO PRODUCE THE CORRECT SIGN AT OUTPUT
            self.predicted_nt = talwani_and_heirtzler.nt(self.xp, self.mag_observation_elv, mag_input_polygons) * -1
        else:
            # SET THE PREDICTED VALUES AS ZEROS
            self.predicted_nt = np.zeros_like(self.xp)

        # SET THE PREDICTED PLOT LINE WITH THE NEWLY CALCULATED VALUES
        self.predicted_nt_plot.set_data(self.xp * 0.001, self.predicted_nt)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # UPDATE RMS VALUES

        # RUN THE RMS CALC CODE
        self.model_rms(self.xp)

        # SET GRAVITY RMS
        if self.obs_gravity_data_for_rms != [] and self.calc_grav_switch is True and self.predicted_gravity != []:
            self.gravity_rms_plot.set_data(self.grav_residuals[:, 0], self.grav_residuals[:, 1])
        else:
            pass

        # SET MAGNETIC RMS
        if self.obs_mag_data_for_rms != [] and self.calc_mag_switch is True and self.predicted_nt != []:
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
            ymin_list.append(self.predicted_gravity.min())
            ymax_list.append(self.predicted_gravity.max())

            # # APPEND RMS GRAVITY ANOMALY
            # ymin_list.append(self.grav_residuals.min() - 2.0)
            # ymax_list.append(self.grav_residuals.max() + 2.0)

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
            if self.predicted_gravity is not None:
                ymin_list.append(self.predicted_gravity.min() - 2.0)
                ymax_list.append(self.predicted_gravity.max() + 2.0)

            # SET YMIN AND YMAX
            ymin = min(ymin_list)
            ymax = max(ymax_list)

        elif self.predicted_gravity is not None:
            ymin = self.predicted_gravity.min() - 2.0
            ymax = self.predicted_gravity.max() + 2.0
        else:
            pass

        if self.gravity_frame is not None:
            # print('')
            # print(ymin)
            # print(ymax)
            # print('')
            self.gravity_frame.set_ylim(ymin, ymax)
        # --------------------------------------------------------------------------------------------------------------
        # SET DERIVATIVE Y-AXIS LIMITS

        # CREATE EMPTY LIST
        ymin_list = [-1]
        ymax_list = [1]
        for i in range(len(self.observed_gravity_list)):
            if self.observed_gravity_list[i] == str('derivative'):
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
            ymin_list.append(self.predicted_nt.min())
            ymax_list.append(self.predicted_nt.max())

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
            ymin_list.append(self.predicted_nt.min() - 2.0)
            ymax_list.append(self.predicted_nt.max() + 2.0)

            # SET YMIN AND YMAX
            ymin = min(ymin_list)
            ymax = max(ymax_list)

        elif self.predicted_nt is not None:
            # APPEND PREDICTED GRAVITY ANOMALY
            ymin = self.predicted_nt.min() - 2.0
            ymax = self.predicted_nt.max() + 2.0
        else:
            pass

        if self.magnetic_frame is not None:
            self.magnetic_frame.set_ylim(ymin, ymax)

        # SET DERIVATIVE Y-AXIS LIMITS
        # --------------------------------------------------------------------------------------------------------------
        # CREATE EMPTY LIST
        ymin_list = []
        ymax_list = []
        for i in range(len(self.observed_magnetic_list)):
            if self.observed_magnetic_list[i].type == str('derivative'):
                ymin_list.append(self.observed_magnetic_list[i].data[:, 1].min() - 0.1)
                ymax_list.append(self.observed_magnetic_list[i].data[:, 1].max() + 0.1)
        self.magnetic_d_frame.set_ylim(ymin, ymax)
        # --------------------------------------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------------------------------------

        # UPDATE GMG GRAPHICS
        self.draw()

    # EXTERNAL FIGURE CONSTRUCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def plot_model(self, event):
        """CREATE EXTERNAL FIGURE OF MODEL USING INBUILT FIGURE CONSTRUCTION TOOL"""

        # GET PLOTTING PARAMETERS FROM DIALOG BOX
        self.set_values = PlotSettingsDialog(self, -1, 'Set figure parameters', self.model_aspect,
                                             self.grav_frame_aspect)
        self.set_values.Show(True)

    def draw_model(self):
        # GET USER INPUT FROM POPOUT BOX
        self.file_path = self.set_values.file_path
        self.file_type = self.set_values.file_type
        self.use_tight_layout = self.set_values.use_tight_layout

        self.fs = self.set_values.fs  # FONT SIZE
        self.aspect_ratio = self.set_values.aspect_ratio  # MODEL ASPECT RATIO
        self.ps = self.set_values.ps  # OBSERVED POINT SIZE
        self.calc_line_width = self.set_values.lw  # CALCUALTED LINE WIDTH

        self.font_type = self.set_values.font_type_text.GetValue()

        self.topo_frame_min = self.set_values.topo_min_text.GetValue()
        self.topo_frame_max = self.set_values.topo_max_text.GetValue()
        self.grav_frame_min = self.set_values.grav_min_text.GetValue()
        self.grav_frame_max = self.set_values.grav_max_text.GetValue()
        self.mag_frame_min = self.set_values.mag_min_text.GetValue()
        self.mag_frame_max = self.set_values.mag_max_text.GetValue()

        self.draw_polygons = self.set_values.draw_polygons
        self.polygon_alpha = self.set_values.polygon_alpha

        self.draw_fixed_layers = self.set_values.draw_fixed_layers
        self.layer_line_width = self.set_values.layer_line_width

        self.draw_floating_layers = self.set_values.draw_floating_layers
        self.layer_line_alpha = self.set_values.layer_line_alpha

        self.draw_colorbar = self.set_values.draw_colorbar
        self.colorbar_x = self.set_values.colorbar_x
        self.colorbar_y = self.set_values.colorbar_y
        self.colorbar_size_x = self.set_values.colorbar_size_x
        self.colorbar_size_y = self.set_values.colorbar_size_y

        self.draw_xy_data = self.set_values.draw_xy_data
        self.xy_size = self.set_values.xy_size
        self.xy_color = self.set_values.xy_color

        self.draw_wells = self.set_values.draw_wells
        self.well_fs = self.set_values.well_fs
        self.well_line_width = self.set_values.well_line_width

        self.draw_faults = self.set_values.draw_faults
        self.faults_lw = self.set_values.faults_lw

        # GET FIGURE DIMENSIONS
        xmin, xmax = self.model_frame.get_xlim()
        ymin, ymax = self.model_frame.get_ylim()
        area = np.array([xmin, xmax, ymin, ymax])

        # # RUN PLOT MODEL CODE
        fig_plot = plot_model.plot_fig(self.file_path, self.file_type, self.use_tight_layout, self.fs,
                                       self.aspect_ratio, self.ps, self.calc_line_width, self.font_type,
                                       self.topo_frame_min, self.topo_frame_max, self.grav_frame_min,
                                       self.grav_frame_max, self.mag_frame_min, self.mag_frame_max,
                                       self.draw_polygons, self.polygon_alpha, self.draw_fixed_layers,
                                       self.layer_line_width, self.draw_floating_layers, self.layer_line_alpha,
                                       self.draw_colorbar, self.colorbar_x, self.colorbar_y, self.colorbar_size_x,
                                       self.colorbar_size_y, self.draw_xy_data, self.xy_size, self.xy_color,
                                       self.draw_wells, self.well_fs, self.well_line_width, self.draw_faults,
                                       self.faults_lw,
                                       self.layer_list, self.fault_list, self.observed_topography_list,
                                       self.observed_gravity_list, self.observed_magnetic_list,
                                       self.outcrop_data_list, self.well_data_list, self.segy_data_list,
                                       self.topo_frame, self.gravity_frame, self.magnetic_frame, self.predicted_gravity,
                                       self.gravity_rms_value, self.predicted_nt, self.magnetic_rms_value, self.area,
                                       self.xp)
        del fig_plot
        #
        # # IF ON A LINUX SYSTEM OPEN THE FIGURE WITH PDF VIEWER
        # try:
        #     if sys.platform == 'linux2':
        #         subprocess.call(["xdg-open", self.file_path])
        #     # IF ON A macOS SYSTEM OPEN THE FIGURE WITH PDF VIEWER
        #     elif sys.platform == 'darwin':
        #         os.open(self.file_path)
        # except IOError:
        #     pass

        # UPDATE GMG
        self.update_layer_data()
        self.draw()

        return

    # DOCUMENTATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def open_documentation(self, event):
        """OPEN DOCUMENTATION HTML"""

        self.doc_dir = os.path.dirname(os.path.abspath(__file__)).split('/')
        doc_url = self.doc_dir[0] + '/' + self.doc_dir[1] + '/' + self.doc_dir[2] + '/' + self.doc_dir[3] + \
                  '/' + self.doc_dir[4] + '/docs/html/gmg_documentation.html'

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
        licence = ["Copyright 2015-2020 Brook Tozer \n\nRedistribution and use in source and binary forms, with or "
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

    # EXIT FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# START SOFTWARE
if __name__ == "__main__":
    app = wx.App(False)
    fr = wx.Frame(None, title='GMG: Geophysical Modelling GUI')
    app.frame = Gmg()
    app.frame.CenterOnScreen()
    app.frame.Show()
    app.MainLoop()
