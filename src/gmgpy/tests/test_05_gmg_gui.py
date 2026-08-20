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


class _RecordingMessageDialog(object):
    """
    STANDS IN FOR THE MESSAGE DIALOG, RECORDING WHAT WAS REPORTED INSTEAD OF BLOCKING ON A MODAL.

    The saves report both success and failure this way, so recording the text is how a test can tell which
    happened, and lets a reported failure be asserted on rather than merely not crashing.
    """
    messages = []

    def __init__(self, parent, id, message_text, title):
        _RecordingMessageDialog.messages.append(message_text)

    @classmethod
    def reset(cls):
        cls.messages = []

    @classmethod
    def said(cls, fragment):
        """True if any recorded message contains the given text."""
        return any(fragment in message for message in cls.messages)


def _file_dialog_returning(path):
    """BUILD A STAND IN FOR THE FILE DIALOG THAT ANSWERS WITH THE GIVEN PATH"""

    class _FileDialog(_AcceptingDialog):
        def GetPath(self):
            return str(path)

    return _FileDialog


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
def reported(gmg_module, monkeypatch):
    """Capture what the saves report, instead of letting a modal block the run."""
    _RecordingMessageDialog.reset()

    # save_xy_file LIVES IN dialogs, SO ITS MESSAGES GO THROUGH THAT MODULE'S NAME WHICHEVER CALLER USED IT.
    # write_c_xy REPORTS THROUGH gmg'S OWN NAME, SO BOTH ARE REPLACED
    monkeypatch.setattr(sys.modules['dialogs'], 'MessageDialog', _RecordingMessageDialog)
    monkeypatch.setattr(gmg_module, 'MessageDialog', _RecordingMessageDialog)

    return _RecordingMessageDialog


def _capture_points(gmg_frame, points):
    """Turn on coordinate capture and click the given (x, y) points through the real handler."""
    from matplotlib.backend_bases import MouseEvent, MouseButton

    gmg_frame.capture_coordinates(None)

    for x, y in points:
        event = MouseEvent('button_release_event', gmg_frame.fig.canvas, 100, 100, button=MouseButton.LEFT)
        event.inaxes = gmg_frame.model_frame
        event.xdata, event.ydata = x, y
        gmg_frame.button_release(event)

    return gmg_frame.capture_window


def _press(button):
    """Press a wx button the way a user would, rather than calling its handler directly."""
    event = wx.CommandEvent(wx.EVT_BUTTON.typeId, button.GetId())
    event.SetEventObject(button)
    button.GetEventHandler().ProcessEvent(event)


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
    assert not gmg_frame.unreliable_magnetic_observation_points().any()
    assert '2D only' not in gmg_frame.magnetic_frame.get_ylabel()
    assert len(gmg_frame.untrusted_nt_plot.get_xdata()) == 0

    # RAISE THE LAYER'S TOP ABOVE THE CALCULATION ELEVATION SO THAT IT STRADDLES IT.
    # NB: THE NODES ARE EDITED IN PLACE BECAUSE gmg.current_y_nodes ALIASES THE SAME LIST AND update_layer_data
    # WRITES IT BACK OVER THE LAYER, SO REBINDING layer.y_nodes ALONE WOULD BE UNDONE
    layer.y_nodes[:] = [-2.0, -2.0, 15.0, 15.0]
    gmg_frame.update_layer_data()
    gmg_frame.run_algorithms()

    inside = gmg_frame.unreliable_magnetic_observation_points()
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
    assert 'unreliable' in gmg_frame.magnetic_frame.get_ylabel()

    # A LAYER WITH NO SUSCEPTIBILITY CANNOT MAKE A MAGNETIC ANOMALY UNTRUSTWORTHY
    layer.susceptibility = 0.0
    gmg_frame.run_algorithms()
    assert not gmg_frame.unreliable_magnetic_observation_points().any()
    assert 'unreliable' not in gmg_frame.magnetic_frame.get_ylabel()


def test_a_node_on_an_observation_point_is_marked(gmg_frame):
    """
    Test an observation point coinciding with a layer node is marked as unreliable

    The magnetic field at a corner of a uniformly magnetised body diverges logarithmically from every direction
    of approach, so no finite value exists there. The kernels nudge the node clear so that the calculation
    returns a number rather than nan, but that number is set by the size of the nudge and means nothing. A node
    is on the boundary rather than inside it, so a point-in-polygon test alone does not catch this.
    """
    layer = gmg_frame.layer_list[gmg_frame.currently_active_layer_id]
    assert not gmg_frame.unreliable_magnetic_observation_points().any()

    # PUT THE LAYER'S TOP EDGE ON THE CALCULATION ELEVATION, WHICH LANDS ITS TWO TOP NODES ON OBSERVATION POINTS
    # WITHOUT PUTTING ANY OBSERVATION POINT INSIDE THE LAYER
    layer.y_nodes[:] = [0.0, 0.0, 15.0, 15.0]
    gmg_frame.update_layer_data()
    gmg_frame.run_algorithms()

    unreliable = gmg_frame.unreliable_magnetic_observation_points()
    assert unreliable.sum() == 2, 'exactly the two top nodes should be flagged'
    np.testing.assert_allclose(np.sort(gmg_frame.xp[unreliable]),
                               [1000.0 * min(layer.x_nodes), 1000.0 * max(layer.x_nodes)])
    assert 'unreliable' in gmg_frame.magnetic_frame.get_ylabel()


def test_self_intersecting_layers_are_reported(gmg_frame):
    """
    Test a layer whose outline crosses itself is detected and reported in the status bar

    Such a layer is ambiguous input and every anomaly it contributes to is affected, not only the magnetic one,
    so it is reported against the model rather than against any one panel.
    """
    layer = gmg_frame.layer_list[gmg_frame.currently_active_layer_id]
    assert gmg_frame.self_intersecting_layers() == []

    # DRAG ONE NODE PAST ANOTHER TO MAKE A BOW TIE, WHICH IS ALL IT TAKES
    layer.x_nodes[:] = [20.0, 60.0, 20.0, 60.0]
    gmg_frame.update_layer_data()
    gmg_frame.run_algorithms()

    assert gmg_frame.self_intersecting_layers() == [gmg_frame.currently_active_layer_id]
    gmg_frame.display_info()
    assert 'SELF-INTERSECTING' in gmg_frame.statusbar.GetStatusText(1)

    # UNTANGLING IT CLEARS THE REPORT
    layer.x_nodes[:] = list(LAYER_X_NODES)
    gmg_frame.update_layer_data()
    gmg_frame.run_algorithms()

    assert gmg_frame.self_intersecting_layers() == []
    gmg_frame.display_info()
    assert 'SELF-INTERSECTING' not in gmg_frame.statusbar.GetStatusText(1)


def test_captured_coordinates_are_saved(gmg_frame, reported, tmp_path, monkeypatch):
    """
    Test clicking points in capture mode and pressing Save writes them to a file
    """
    target = tmp_path / 'captured.xy'
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(target))

    window = _capture_points(gmg_frame, [(10.5, 3.25), (20.0, 7.5), (33.125, 12.0)])
    assert window.table.GetItemCount() == 3

    _press(window.save_btn)

    assert target.exists()
    assert target.read_text().split() == ['10.500000', '3.250000',
                                          '20.000000', '7.500000',
                                          '33.125000', '12.000000']
    assert reported.said('Saved successfully')


def test_saving_with_nothing_captured_writes_no_file(gmg_frame, reported, tmp_path, monkeypatch):
    """
    Test pressing Save before any point is captured says so, rather than leaving an empty file

    numpy validates the format against the column count only after opening the file, so an empty dataset used to
    raise ValueError and leave a 0 byte file behind, with the exception swallowed by wx.
    """
    target = tmp_path / 'nothing.xy'
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(target))

    window = _capture_points(gmg_frame, [])
    assert window.table.GetItemCount() == 0

    _press(window.save_btn)

    assert not target.exists(), 'an empty save must not leave a file behind'
    assert reported.said('nothing to save')
    assert not reported.said('Saved successfully')


def test_saving_to_an_unwritable_path_is_reported(gmg_frame, reported, monkeypatch):
    """
    Test a path that cannot be written is reported rather than silently doing nothing

    This is the case that produced no file at all, and no message, because wx swallows an exception raised in an
    event handler and pythonw discards stderr.
    """
    unwritable = '/System/no/such/directory/coordinates.xy'
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(unwritable))

    window = _capture_points(gmg_frame, [(1.0, 2.0)])
    _press(window.save_btn)

    assert not os.path.exists(unwritable)
    assert reported.said('could not be written')
    assert not reported.said('Saved successfully')


def test_a_filename_without_a_suffix_gains_one(gmg_frame, reported, tmp_path, monkeypatch):
    """
    Test a name typed without an extension is written with the one the file dialog filtered on
    """
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(tmp_path / 'no_suffix'))

    window = _capture_points(gmg_frame, [(1.0, 2.0)])
    _press(window.save_btn)

    assert (tmp_path / 'no_suffix.xy').exists()
    assert not (tmp_path / 'no_suffix').exists()


def test_predicted_anomaly_save_reports_and_writes(gmg_frame, reported, tmp_path, monkeypatch):
    """
    Test the shared save helper is reached from gmg's own save routines, not only the capture window
    """
    target = tmp_path / 'predicted_gravity.txt'
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(target))

    gmg_frame.save_modelled_grav(None)

    assert target.exists()
    assert len(target.read_text().splitlines()) == len(gmg_frame.xp)
    assert reported.said('Saved successfully')


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
