"""
GMG OBJECT CLASSES GO HERE
"""

import matplotlib
matplotlib.use('WXAgg')
import matplotlib.cm as cm
from matplotlib.artist import Artist


def holds_mpl_artist(value):
    """
    TRUE IF value IS A MATPLOTLIB ARTIST, OR A LIST OR TUPLE HOLDING ONE.

    NB: THE OBJECTS BELOW STORE ARTISTS BOTH DIRECTLY (e.g. AN Annotation) AND AS THE LISTS RETURNED BY
    Axes.plot AND Axes.fill, SO BOTH FORMS HAVE TO BE RECOGNISED
    """
    if isinstance(value, Artist):
        return True
    elif isinstance(value, (list, tuple)):
        return any(holds_mpl_artist(item) for item in value)
    else:
        return False


def walk_objects(value):
    """
    YIELD EVERY gmg OBJECT REACHABLE WITHIN value, DESCENDING THROUGH LISTS AND TUPLES.

    :param value: an object, or a (possibly nested) list or tuple of them
    """
    if isinstance(value, (list, tuple)):
        for item in value:
            for found in walk_objects(item):
                yield found
    elif hasattr(value, '__dict__') and not isinstance(value, Artist):
        yield value


def strip_mpl_artists(obj):
    """
    SET EVERY ATTRIBUTE OF obj THAT HOLDS A MATPLOTLIB ARTIST TO None.

    Matplotlib artists must not be written into a .model file. They are recreated from the saved data whenever a
    model is loaded, so they carry nothing, and pickling them ties the file to the matplotlib version that wrote
    it. Which attributes hold artists is decided by type rather than by name, so an artist added to one of the
    objects in this module later cannot be missed.

    :param obj: the object to strip

    :returns: dict : the attributes that were removed, so that the caller can put them back
    """
    removed = {name: value for name, value in vars(obj).items() if holds_mpl_artist(value)}

    for name in removed:
        setattr(obj, name, None)

    return removed


class Layer:
    """GENERIC LAYER (POLYGON) OBJECT. THE LAYER WILL BE STORED IN THE gmg.layer_list LIST"""

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
        # OUT OF PLANE (STRIKE DIRECTION) EXTENTS IN km, MEASURED FROM THE MODEL SECTION. None MEANS THE LAYER IS OF
        # INFINITE EXTENT, I.E. IT IS MODELLED IN 2D. SET BOTH TO MAKE THE LAYER 2.75D
        self.y1 = None
        self.y2 = None
        self.color = 'k'
        self.layer_transparency = 0.4
        self.pinched = False  # SWITCH DICTATING IF THE LAYER HSA ANY NODES PINCHED TO ANOTHER
        self.pinched_list = []  # LIST OF LAYER ID's FOR LAYERS THAT HAVE NODES PINCHED TO THE CURRENT LAYER
        self.include_in_calculations_switch = True  # SWITCH THAT DICTATES IF THE LAYER IS INCLUDED IN THE CURRENT CALC


class Fault:
    """GENERIC FAULT OBJECT. THE FAULT WILL BE STORED IN THE gmg.fault_list LIST"""

    def __init__(self):
        self.id = None  # THE FAULTS NUMBER
        self.name = None  # THE FAULTS NAME
        self.color = 'k'
        self.mpl_actor = None  # THE MPL ACTOR ELEMENT (PLOTTING OBJECT)
        self.x_nodes = []  # X NODES OF THE FAULT
        self.y_nodes = []  # Y NODES OF THE FAULT


class ObservedData:
    def __init__(self):
        """GENERIC OBSERVATION DATA OBJECT"""

        self.id = None  # OBJECT ID VALUE FOR wx
        self.name = None  # THE NAME ASSIGNED TO THE DATA (str)
        self.color = None  # THE COLOR USED FOR PLOTTING THE DATA (str)
        self.type = None  # THE KIND OF DATA (str: 'observed', 'filtered', 'derivative')
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.mpl_actor = None  # THE MPL ACTOR ELEMENT (PLOTTING OBJECT)


class ObservedOutcropData:
    """GENERIC GEOLOGICAL OUTCROP DATA OBJECT"""

    def __init__(self):
        self.id = None
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE (numpy array)
        self.name = None  # THE NAME ASSIGNED TO THE DATA (str)
        self.color = None  # THE COLOR USED FOR PLOTTING THE DATA (str)
        self.textsize = 2  # THE SIZE OF TEXT LABELS


class SegyData:
    """GENERIC SEGY DATA OBJECT"""

    def __init__(self):
        self.file = None  # THE FULL FILE PATH TO THE SEGY DATA : str
        self.name = None  # THE NAME ASSIGNED TO THE DATA : str
        self.dimensions = None # LIST CONTAINING x1 x2 z1 z2
        self.mpl_actor = None  # THE MPL ACTOR ELEMENT (PLOTTING OBJECT)
        self.axis = None  # THE X AND Z AXIS LIMITS = X1, X2, Z1, Z2 : list : float
        self.color_map = cm.gray  # THE COLOR SCALE USED TO PLOT THE DATA
        self.data = None  # THE XY DATA LOADED FROM INPUT FILE : numpy array : float
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
