"""
GMG OBJECT CLASSES GO HERE
"""

import matplotlib
matplotlib.use('WXAgg')
import matplotlib.cm as cm
from frames import *
from dialogs import *


class Layer:
    """GENERIC MODEL LAYER (POLYGON) OBJECT. THE LAYER WILL BE STORED IN THE gmg.layer_list LIST"""

    def __init__(self):
        self.id = None  # THE LAYER NUMBER
        self.name = None  # THE LAYER NAME
        self.type = None  # EITHER (1) str('fixed') OR (2) str('floating')
        self.node_mpl_actor = None
        self.polygon_mpl_actor = None
        self.poly_plot = None
        self.x_nodes = []
        self.y_nodes = []
        self.polygon = None  # x_nodes AND y_nodes ZIPPED. TO BE PASSED TO THE POTENTIAL FIELD ALGORITHMS
        self.density = 0.
        self.reference_density = 0.
        self.susceptibility = 0.
        self.angle_a = 0.
        self.angle_b = 0.
        self.angle_c = 0.
        self.earth_field = 0.
        self.color = 'k'
        self.layer_transparency = 0.4
        self.pinched = False  # SWITCH DICTATING IF THE LAYER HSA ANY NODES PINCHED TO ANOTHER
        self.pinched_list = []  # LIST OF LAYER ID's FOR LAYERS THAT HAVE NODES PINCHED TO THE CURRENT LAYER
        self.include_in_calculations_switch = True  # SWITCH THAT DICTATES IF THE LAYER IS INCLUDED IN THE CURRENT CALC


class Fault:
    """GENERIC MODEL FAULT OBJECT. THE FAULT WILL BE STORED IN THE gmg.fault_list LIST"""

    def __init__(self):
        self.id = None  # THE FAULTS NUMBER
        self.name = None  # THE FAULTS NAME
        self.color = 'k'
        self.mpl_actor = None  # THE MPL ACTOR ELEMENT (PLOTTING OBJECT)
        self.x_nodes = []  # X NODES OF THE FAULT
        self.y_nodes = []  # Y NODES OF THE FAULT


class ObservedData:
    def __init__(self):
        """GENERIC OBSERVATIONAL DATA OBJECT"""
        self.id = None  # OBJECT ID VALUE FOR wx
        self.name = None  # THE NAME ASSIGNED TO THE DATA (str)
        self.color = None  # THE COLOR USED FOR PLOTTING THE DATA (str)
        self.type = None  # THE KIND OF DATA (str: 'observed', 'filtered', 'derivative')
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.mpl_actor = None  # THE MPL ACTOR ELEMENT (PLOTTING OBJECT)


class ObservedOutcropData:
    """GENERIC OBSERVATIONAL OUTCROP DATA OBJECT"""

    def __init__(self):
        self.id = None
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.name = None  # THE NAME ASSIGNED TO THE DATA (str)
        self.color = None  # THE COLOR USED FOR PLOTTING THE DATA (str)
        self.textsize = 2  # THE SIZE OF TEXT LABELS


class SegyData:
    """GENERIC SEGY DATA OBJECT"""

    def __init__(self):
        self.file = None  # THE FULL FILE PATH TO THE SEGY DATA
        self.name = None  # THE NAME ASSIGNED TO THE DATA (str)
        self.dimensions = None
        self.mpl_actor = None  # THE MPL ACTOR ELEMENT (PLOTTING OBJECT)
        self.axis = None  # THE X AND Z AXIS LIMITS = X1, X2, Z1, Z2
        self.color_map = cm.gray  # THE COLOR SCALE USED TO PLOT THE DATA
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.gain_positive = 4.0
        self.gain_neg = -self.gain_positive
        self.segy_show = False
        self.plot_list = None


class ObservedWellData:
    """GENERIC WELL DATA OBJECT"""

    def __init__(self):
        self.name = None  # NAME OF WELL RECORD
        self.raw_record = None  # RAW TEXT RECORD
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.mpl_actor = None  # MPL ACTOR USED TO STORE MPL LINE
        self.mpl_actor_name = None  # MPL ACTOR USED TO STORE WELL NAME LABEL
        self.labels_list = []
        self.horizons_list = []
        self.text_size = 2
