"""
GMG OBJECT CLASSES GO HERE
"""

import matplotlib
matplotlib.use('WXAgg')
import matplotlib.cm as cm
from frames import *
from dialogs import *


class Layer:
    """
    CREATE A MODEL LAYER OBJECT.
    THE LAYER WILL BE STORED IN THE *gmg.layers_list'* LIST
    """
    def __init__(self):
        self.mpl_actor = None
        self.poly_plot = None
        self.poly_fill = None
        self.x_nodes = []
        self.y_nodes = []
        self.density = 0.
        self.reference_density = 0.
        self.susceptibility = 0.
        self.angle_a = 0.
        self.angle_b = 0.
        self.color = 0.
        self.boundary_lock_list = [0]
        self.boundary_lock_status = ['locked']
        self.layer_lock_list = [0]
        self.layer_lock_status = ['locked']
        self.layer_transparency = 0.4
        self.pinch_node_list = [[], []]
        self.pinch_count = 0
        self.layers_calculation_switch = True  # SWITCH TO DICTATE IF THE LAYER IS INCLUDED IN THE CURRENT CALCULATIONS


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
        self.id = None  # OBJECT ID VALUE FOR wx
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