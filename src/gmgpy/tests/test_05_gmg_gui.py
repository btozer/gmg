"""
Drive the gmg application itself, rather than the forward modelling algorithms in isolation.

These tests cover the parts of gmg that only exist once there is a running application: saving and loading a
model, migrating a model written by an older gmg, and the layer attribute controls in the sidebar.

They need wxPython and a display, and gmg pulls in pygmt, so every one of them skips rather than fails when any
of that is missing. The rest of the suite is pure numpy and is unaffected.

.. note:: gmg.py imports its sibling modules flat rather than as a package, so its own directory has to be put on
          ``sys.path`` before it can be imported. It is then loaded from its path rather than by name, because the
          repository root holds an empty ``__init__.py`` that makes it importable as a package called ``gmg`` too,
          and that one shadows the real module whenever the root's parent is on ``sys.path``.
"""
import importlib.util
import os
import pickle
import sys

import numpy as np
import pytest

wx = pytest.importorskip('wx', reason='the gmg application tests need wxPython')

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PACKAGE_DIR not in sys.path:
    sys.path.insert(0, PACKAGE_DIR)

# THE MODEL BUILT FOR EVERY TEST: ONE FLOATING LAYER WITH BOTH A DENSITY AND A SUSCEPTIBILITY CONTRAST
MODEL_X1, MODEL_X2 = 0.0, 100.0
MODEL_Z1, MODEL_Z2 = 0.0, 40.0
LAYER_X_NODES = [20.0, 60.0, 60.0, 20.0]
LAYER_Y_NODES = [5.0, 5.0, 15.0, 15.0]
DENSITY = 2900.0
REFERENCE_DENSITY = 2700.0
SUSCEPTIBILITY = 0.02
FIELD_INTENSITY = 50000.0
ANGLE_A, ANGLE_B, ANGLE_C = 60.0, 10.0, 45.0

# DELIBERATELY ASYMMETRIC, SO THAT A SIGN ERROR IN THE STRIKE DIRECTION WOULD SHOW
EXTENT_Y1, EXTENT_Y2 = -7.5, 32.5

# THE wx APPLICATION, KEPT ALIVE FOR THE WHOLE SESSION BY THE gmg_module FIXTURE
_APP = None


class _AcceptingDialog(object):
    """A DIALOG THAT ACCEPTS ITSELF IMMEDIATELY, SO THE APPLICATION CAN BE DRIVEN WITHOUT A USER"""

    def __init__(self, *args, **kwargs):
        pass

    def ShowModal(self):
        return wx.ID_OK

    def Destroy(self):
        pass


def _load_gmg():
    """
    Import gmg.py from its path, registering it under the name its own modules expect.

    Importing it by name is not safe here: the repository root holds an empty ``__init__.py``, so it is also
    importable as a package called ``gmg``, and that one wins whenever the root's parent is on ``sys.path``,
    which is the case under pytest.
    """
    spec = importlib.util.spec_from_file_location('gmg', os.path.join(PACKAGE_DIR, 'gmg.py'))
    module = importlib.util.module_from_spec(spec)
    sys.modules['gmg'] = module
    spec.loader.exec_module(module)

    return module


def _new_model_dialog(x1, x2, z1, z2):
    """BUILD A STAND IN FOR NewModelDialog THAT ANSWERS WITH THE GIVEN MODEL DIMENSIONS"""

    class _NewModelDialog(_AcceptingDialog):
        def __init__(self, *args, **kwargs):
            self.set_button = True
            self.x1, self.x2, self.z1, self.z2 = x1, x2, z1, z2
            self.xp_inc = 1.0

    return _NewModelDialog


@pytest.fixture(scope='session')
def gmg_module():
    """
    Import gmg once, with its modal dialogs replaced so that nothing waits for a user.

    Skips the whole module when there is no display, no wxPython or no working pygmt, since gmg cannot start.
    """
    # NB: THE APPLICATION HAS TO BE KEPT ALIVE FOR THE WHOLE SESSION. LET IT BE COLLECTED AND EVERY LATER wx CALL
    # RAISES PyNoAppError
    global _APP
    try:
        _APP = wx.App(False)
    except Exception as exc:  # NO DISPLAY, OR wx CANNOT TALK TO THE WINDOW SERVER
        pytest.skip('a wx application could not be started: %s' % exc)

    try:
        gmg = _load_gmg()
    except Exception as exc:  # MOST LIKELY pygmt OR ONE OF ITS DEPENDENCIES
        pytest.skip('gmg could not be imported: %s' % exc)

    # wx 4.2 RAISES ON DEBUG ASSERTIONS THAT THE APPLICATION SURVIVES WHEN RUN NORMALLY
    wx.DisableAsserts()

    # THE "MODEL SAVED" AND SIMILAR ACKNOWLEDGEMENTS WOULD OTHERWISE BLOCK
    gmg.MessageDialog = _AcceptingDialog
    gmg.NewModelDialog = _new_model_dialog(MODEL_X1, MODEL_X2, MODEL_Z1, MODEL_Z2)

    return gmg


@pytest.fixture
def model_path(tmp_path):
    """Where a test may write a .model file, and what the file dialogs will answer with."""
    return tmp_path / 'test.model'


@pytest.fixture
def gmg_frame(gmg_module, model_path, monkeypatch):
    """A running gmg holding a one layer model, with the file dialogs pointed at model_path."""
    chosen = {'path': str(model_path)}

    class _FileDialog(_AcceptingDialog):
        def GetPath(self):
            return chosen['path']

    monkeypatch.setattr(wx, 'FileDialog', _FileDialog)
    monkeypatch.setattr(wx, 'MessageBox', lambda *args, **kwargs: wx.YES)

    frame = gmg_module.Gmg()
    frame.new_model(None)

    # ALL THREE CALCULATIONS ON, AS THE USER WOULD SWITCH THEM FROM THE TOOLBAR
    frame.calc_grav_switch = True
    frame.calc_vgg_switch = True
    frame.calc_mag_switch = True

    # new_layer GATHERS ITS NODES FROM FOUR MOUSE CLICKS, SO THEY ARE SUPPLIED DIRECTLY HERE
    frame.new_plotx = list(LAYER_X_NODES)
    frame.new_ploty = list(LAYER_Y_NODES)
    frame.create_new_floating_layer()

    layer = frame.layer_list[frame.currently_active_layer_id]
    layer.density = DENSITY
    layer.reference_density = REFERENCE_DENSITY
    layer.susceptibility = SUSCEPTIBILITY
    layer.earth_field = FIELD_INTENSITY
    layer.angle_a, layer.angle_b, layer.angle_c = ANGLE_A, ANGLE_B, ANGLE_C

    frame.update_layer_data()
    frame.run_algorithms()

    return frame


def test_model_round_trips(gmg_frame):
    """
    Test a saved model reloads with its layers, its attributes and its anomalies intact

    Loading used to append a spurious layer with no polygon, which then crashed the forward calculation, and
    loading with nothing outstanding to save raised UnboundLocalError before the file dialog was ever built.
    """
    layer_id = gmg_frame.currently_active_layer_id
    gmg_frame.layer_list[layer_id].y1 = EXTENT_Y1
    gmg_frame.layer_list[layer_id].y2 = EXTENT_Y2
    gmg_frame.run_algorithms()

    layer_count = len(gmg_frame.layer_list)
    expected_gravity = gmg_frame.predicted_gravity.copy()
    expected_magnetics = gmg_frame.predicted_nt.copy()

    gmg_frame.save_model(None)
    gmg_frame.model_saved = True  # THE PATH THAT USED TO RAISE UnboundLocalError
    gmg_frame.load_model(None)

    assert len(gmg_frame.layer_list) == layer_count
    assert all(layer.polygon is not None for layer in gmg_frame.layer_list)
    assert (gmg_frame.layer_list[layer_id].y1, gmg_frame.layer_list[layer_id].y2) == (EXTENT_Y1, EXTENT_Y2)

    gmg_frame.run_algorithms()
    np.testing.assert_allclose(gmg_frame.predicted_gravity, expected_gravity)
    np.testing.assert_allclose(gmg_frame.predicted_nt, expected_magnetics)


def test_loading_repeatedly_does_not_accumulate_layers(gmg_frame):
    """
    Test loading the same model several times leaves the layer count unchanged
    """
    layer_count = len(gmg_frame.layer_list)
    gmg_frame.save_model(None)

    for _ in range(3):
        gmg_frame.model_saved = True
        gmg_frame.load_model(None)
        assert len(gmg_frame.layer_list) == layer_count


def test_saved_model_holds_no_matplotlib_artists(gmg_frame, model_path):
    """
    Test a saved model carries no matplotlib artists, and that the live model keeps its own

    Pickling the artists tied a saved file to one matplotlib version, which is what stopped older models from
    opening at all.
    """
    artists_before = [layer.node_mpl_actor for layer in gmg_frame.layer_list]

    gmg_frame.save_model(None)

    assert b'matplotlib' not in model_path.read_bytes()
    assert [layer.node_mpl_actor for layer in gmg_frame.layer_list] == artists_before
    assert all(artist is not None for artist in artists_before)


def test_model_written_before_2_75d_still_loads(gmg_frame, model_path):
    """
    Test a model saved before layers had out of plane extents loads, with every layer of infinite extent
    """
    expected = gmg_frame.predicted_gravity.copy()
    gmg_frame.save_model(None)

    with open(str(model_path), 'rb') as saved_file:
        saved = pickle.load(saved_file)
    for layer in saved['layer_list']:
        del layer.y1
        del layer.y2
    with open(str(model_path), 'wb') as saved_file:
        pickle.dump(saved, saved_file, protocol=pickle.HIGHEST_PROTOCOL)

    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    assert all(layer.y1 is None and layer.y2 is None for layer in gmg_frame.layer_list)
    gmg_frame.run_algorithms()
    np.testing.assert_allclose(gmg_frame.predicted_gravity, expected)


def test_finite_extent_control_changes_the_anomaly(gmg_frame):
    """
    Test the sidebar finite extent checkbox switches a layer between 2D and 2.75D
    """
    infinite_gravity = gmg_frame.predicted_gravity.copy()
    infinite_vgg = gmg_frame.predicted_vgg.copy()

    gmg_frame.y1_input.SetValue(-20.0)
    gmg_frame.y2_input.SetValue(20.0)
    gmg_frame.finite_extent_input.SetValue(True)
    gmg_frame.set_finite_extent(None)

    layer = gmg_frame.layer_list[gmg_frame.currently_active_layer_id]
    assert (layer.y1, layer.y2) == (-20.0, 20.0)
    assert gmg_frame.y1_input.IsEnabled() and gmg_frame.y2_input.IsEnabled()

    # A FINITE BODY MUST PRODUCE A SMALLER ANOMALY THAN AN INFINITE ONE
    assert np.max(np.abs(gmg_frame.predicted_gravity)) < np.max(np.abs(infinite_gravity))

    # THERE IS NO 2.75D VGG, SO IT MUST BE UNTOUCHED AND SAID TO BE 2D ONLY
    np.testing.assert_array_equal(gmg_frame.predicted_vgg, infinite_vgg)
    assert '2D only' in gmg_frame.vertical_gg_frame.get_ylabel()

    # AND SWITCHING BACK MUST RESTORE THE ORIGINAL EXACTLY
    gmg_frame.finite_extent_input.SetValue(False)
    gmg_frame.set_finite_extent(None)
    assert layer.y1 is None and layer.y2 is None
    np.testing.assert_array_equal(gmg_frame.predicted_gravity, infinite_gravity)


def test_inducing_field_direction_changes_the_magnetic_anomaly(gmg_frame, gmg_module, monkeypatch):
    """
    Test stating an inducing field direction separates it from the magnetisation direction

    With magnetisation induced only the two coincide, so leaving the direction unset must reproduce the
    original anomaly exactly.
    """
    induced = gmg_frame.predicted_nt.copy()

    class _InducingField(_AcceptingDialog):
        inclination, declination = 25.0, -5.0

        def ShowModal(self):
            return 1

    monkeypatch.setattr(gmg_module, 'SetInducingFieldDialog', _InducingField)
    gmg_frame.set_inducing_field(None)

    assert (gmg_frame.inducing_field_inclination, gmg_frame.inducing_field_declination) == (25.0, -5.0)
    assert not np.allclose(gmg_frame.predicted_nt, induced)

    _InducingField.inclination, _InducingField.declination = None, None
    gmg_frame.set_inducing_field(None)

    assert gmg_frame.inducing_field_inclination is None
    np.testing.assert_array_equal(gmg_frame.predicted_nt, induced)


def test_interior_observation_points_are_marked(gmg_frame):
    """
    Test observation points falling inside a magnetised layer are detected and marked

    Both magnetic algorithms are wrong at such a point, which arises whenever a magnetised layer crosses the
    calculation elevation. They are marked rather than corrected, so that a wrong value is never read as a good
    one.
    """
    layer = gmg_frame.layer_list[gmg_frame.currently_active_layer_id]

    # AS BUILT THE LAYER SITS WHOLLY BELOW THE CALCULATION ELEVATION, SO NOTHING IS INSIDE IT
    assert not gmg_frame.observation_points_inside_magnetic_layers().any()
    assert '2D only' not in gmg_frame.magnetic_frame.get_ylabel()
    assert len(gmg_frame.untrusted_nt_plot.get_xdata()) == 0

    # RAISE THE LAYER'S TOP ABOVE THE CALCULATION ELEVATION SO THAT IT STRADDLES IT.
    # NB: THE NODES ARE EDITED IN PLACE BECAUSE gmg.current_y_nodes ALIASES THE SAME LIST AND update_layer_data
    # WRITES IT BACK OVER THE LAYER, SO REBINDING layer.y_nodes ALONE WOULD BE UNDONE
    layer.y_nodes[:] = [-2.0, -2.0, 15.0, 15.0]
    gmg_frame.update_layer_data()
    gmg_frame.run_algorithms()

    inside = gmg_frame.observation_points_inside_magnetic_layers()
    assert inside.any(), 'a layer crossing the calculation elevation must contain observation points'

    # THE MARKED SAMPLES MUST BE EXACTLY THOSE DETECTED.
    # NB: A LOOSE TOLERANCE, BECAUSE matplotlib STORES LINE DATA AS float32
    marked = gmg_frame.untrusted_nt_plot.get_xdata()
    detected = gmg_frame.xp[inside] * 0.001
    np.testing.assert_allclose(marked, detected, rtol=1.0e-6)

    # AND THEY MUST ALL LIE WITHIN THE LAYER'S HORIZONTAL SPAN.
    # NB: COMPARED IN METRES, SINCE SCALING BY 0.001 IS NOT EXACT IN BINARY AND PUTS 60000 m JUST OVER 60 km
    assert gmg_frame.xp[inside].min() >= 1000.0 * min(layer.x_nodes)
    assert gmg_frame.xp[inside].max() <= 1000.0 * max(layer.x_nodes)
    assert 'inside a layer' in gmg_frame.magnetic_frame.get_ylabel()

    # A LAYER WITH NO SUSCEPTIBILITY CANNOT MAKE A MAGNETIC ANOMALY UNTRUSTWORTHY
    layer.susceptibility = 0.0
    gmg_frame.run_algorithms()
    assert not gmg_frame.observation_points_inside_magnetic_layers().any()
    assert 'inside a layer' not in gmg_frame.magnetic_frame.get_ylabel()


def test_bulk_setters_reach_every_layer(gmg_frame, gmg_module, monkeypatch):
    """
    Test the bulk background density and out of plane extent setters apply to every layer

    The background density loop used to stop one short, leaving the last layer on its old reference density.
    """
    class _BackgroundDensity(_AcceptingDialog):
        background_density_upper = 2.5

    class _Extent(_AcceptingDialog):
        y1, y2 = -12.0, 18.0

        def ShowModal(self):
            return 1

    monkeypatch.setattr(gmg_module, 'SetBackgroundDensityDialog', _BackgroundDensity)
    monkeypatch.setattr(gmg_module, 'SetOutOfPlaneExtentDialog', _Extent)

    gmg_frame.set_background_density(None)
    assert all(layer.reference_density == 2500.0 for layer in gmg_frame.layer_list)

    gmg_frame.set_out_of_plane_extent(None)
    assert all((layer.y1, layer.y2) == (-12.0, 18.0) for layer in gmg_frame.layer_list)
    assert gmg_frame.finite_extent_input.GetValue()

    _Extent.y1, _Extent.y2 = None, None
    gmg_frame.set_out_of_plane_extent(None)
    assert all(layer.y1 is None and layer.y2 is None for layer in gmg_frame.layer_list)
