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

# ALL FOUR BEHAVE IDENTICALLY ON SAVE AND LOAD, SO THE TESTS BELOW RUN OVER EVERY ONE OF THEM
OBSERVED_PANELS = {
    'topography': ('observed_topography_list', 'topography_frame', 'topography_d_frame',
                   'topography_frame_switch'),
    'gravity': ('observed_gravity_list', 'gravity_frame', 'gravity_d_frame',
                'gravity_frame_switch'),
    'vgg': ('observed_vgg_list', 'vertical_gg_frame', 'vertical_gg_d_frame',
            'vertical_gg_frame_switch'),
    'magnetic': ('observed_magnetic_list', 'magnetic_frame', 'magnetic_d_frame',
                 'magnetic_frame_switch'),
}

# THE wx APPLICATION, KEPT ALIVE FOR THE WHOLE SESSION BY THE gmg_module FIXTURE
_APP = None


class _SkippableEvent(object):
    """
    STANDS IN FOR A wx EVENT IN A HANDLER CALLED DIRECTLY.

    The real event classes take constructor arguments that vary between wx versions, and none of the
    handlers under test read anything from the event beyond skipping it.
    """

    def Skip(self):
        pass


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

    # wx 4.2 RAISES ON DEBUG ASSERTIONS THAT THE APPLICATION SURVIVES WHEN RUN NORMALLY, SUCH AS THE SIZER
    # COMPLAINT FROM size_handler.
    # NB: THIS ALSO HIDES REAL wx FAULTS. A wxWidgets 3.2.9 REGRESSION THAT MADE EVERY SAVE DIALOG RAISE WENT
    # UNSEEN HERE FOR EXACTLY THIS REASON, COMPOUNDED BY wx.FileDialog BEING REPLACED IN THESE TESTS. WHAT CAN BE
    # CHECKED WITHOUT A REAL DIALOG IS IN test_07_file_dialogs.py; DO NOT ASSUME THIS SUITE WOULD CATCH SUCH A BUG
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


def _attach_observed_data(gmg_module, gmg_frame, field, data_type='observed'):
    """
    Give the named panel one observed dataset, the way loading a data file does.

    The real handler reads a file through a modal dialog, so the object is built here directly and then
    registered through the same replot_* function the load path uses.
    """
    list_name, _frame_name, _d_frame_name, _switch_name = OBSERVED_PANELS[field]

    observed = gmg_module.ObservedData()
    observed.id = 0
    observed.type = str(data_type)
    observed.name = '%s test data' % field
    observed.color = 'black'
    observed.data = np.column_stack((np.linspace(MODEL_X1, MODEL_X2, 25), np.arange(25, dtype=float)))
    observed.mpl_actor = None

    getattr(gmg_frame, list_name).append(observed)
    getattr(gmg_frame, 'replot_observed_%s_data' % field)()
    gmg_frame.rebuild_data_frames()

    return observed


def _scatters_in(axes):
    """The scatter artists drawn on the given axes, ignoring the plotted anomaly and RMS lines."""
    from matplotlib.collections import PathCollection

    return [artist for artist in axes.collections if isinstance(artist, PathCollection)]


def _write_png(tmp_path, name='figure.png'):
    """Write a small PNG for the image loaders to read."""
    from matplotlib.image import imsave

    target = tmp_path / name
    imsave(str(target), np.linspace(0.0, 1.0, 64).reshape(8, 8), cmap='gray')

    return target


def _write_grid(tmp_path, name='grid.txt'):
    """Write a small whitespace delimited grid, the format load_grid_data reads."""
    target = tmp_path / name
    np.savetxt(str(target), np.arange(64, dtype=float).reshape(8, 8), delimiter=' ')

    return target


def _write_segy(tmp_path, name='line.segy'):
    """Write a small SEG-Y, the format load_segy reads. gmg imports obspy already, so this is free."""
    from obspy.core import Stream, Trace

    traces = []
    for _ in range(3):
        trace = Trace(data=np.linspace(0.0, 1.0, 40).astype(np.float32))
        trace.stats.delta = 0.004  # SEG-Y REJECTS THE 1.0 s DEFAULT
        traces.append(trace)

    target = tmp_path / name
    Stream(traces).write(str(target), format='SEGY', data_encoding=5)

    return target


def _filter_dialog(name, data, color='green'):
    """
    BUILD A STAND IN FOR GaussianFilterDialog THAT ANSWERS AS IF 'Filter' HAD BEEN PRESSED.

    The handlers read filtered_output, output_name and output_color off the dialog after ShowModal.
    """

    class _GaussianFilterDialog(_AcceptingDialog):
        def __init__(self, *args, **kwargs):
            self.filtered_output = np.array(data, dtype=float)
            self.output_name = name
            self.output_color = color

    return _GaussianFilterDialog


def _derivative_dialog(name, data, color='red'):
    """
    BUILD A STAND IN FOR HorizontalDerivative THAT ANSWERS AS IF 'Calculate' HAD BEEN PRESSED.

    The handlers check hasattr(box, 'deriv') first, so a dialog without it stands for a cancelled one.
    """

    class _HorizontalDerivative(_AcceptingDialog):
        def __init__(self, *args, **kwargs):
            self.deriv = np.array(data, dtype=float)
            self.output_name = name
            self.output_color = color

    return _HorizontalDerivative


def _load_segy(gmg_module, gmg_frame, monkeypatch, segy_path, name='test segy'):
    """Load a segy the way the menu item does, with both of its dialogs answered."""
    monkeypatch.setattr(gmg_module, 'SeisDialog', _seis_dialog(name))
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(segy_path))

    gmg_frame.segy_input(None)


def _seis_dialog(name, dimensions=(0.0, 100.0, 0.0, 40.0)):
    """
    BUILD A STAND IN FOR SeisDialog THAT ANSWERS AS IF 'Set Coordinates' HAD BEEN PRESSED.

    The real dialog assigns these two attributes only in that handler, which is why the loaders have
    to check for them before reading.
    """

    class _SeisDialog(_AcceptingDialog):
        def __init__(self, *args, **kwargs):
            self.dimensions = list(dimensions)
            self.segy_name_input = name

    return _SeisDialog


def _load_png(gmg_module, gmg_frame, monkeypatch, png_path, name='test png'):
    """Load a PNG the way the menu item does, with both of its dialogs answered."""
    monkeypatch.setattr(gmg_module, 'SeisDialog', _seis_dialog(name))
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(png_path))

    gmg_frame.png_input(None)


def _load_grid(gmg_module, gmg_frame, monkeypatch, grid_path, name='test grid'):
    """Load a grid the way the menu item does, with both of its dialogs answered."""
    monkeypatch.setattr(gmg_module, 'SeisDialog', _seis_dialog(name))
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(grid_path))

    gmg_frame.grid_input(None)


def _menu_labels(menu):
    """The item labels of a wx menu, so a test can say what is listed rather than only how many."""
    return [item.GetItemLabelText() for item in menu.GetMenuItems()]


def _fire_menu(gmg_frame, menu_id):
    """
    Fire a menu item by ID, through the frame's own bindings.

    NB: THIS IS DELIBERATELY NOT A DIRECT HANDLER CALL. THE MENU ID CLASH BEING GUARDED AGAINST LIVES IN
    THE BINDING TABLE, AND CALLING THE HANDLER WOULD STEP STRAIGHT PAST IT
    """
    event = wx.CommandEvent(wx.EVT_MENU.typeId, menu_id)
    gmg_frame.GetEventHandler().ProcessEvent(event)


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


@pytest.mark.parametrize('field', sorted(OBSERVED_PANELS))
def test_model_with_a_hidden_panel_holding_data_loads(gmg_module, gmg_frame, field):
    """
    Test a model whose data panel is saved hidden still loads when that panel holds observed data

    Hiding a panel leaves its axes as None. load_model replotted the observed data unconditionally,
    gated only on the data list being non-empty, so the load died on None.scatter() and left gmg with
    a wiped model and no currently_active_layer, at which point clicking the model frame raised too.
    """
    list_name, frame_name, _d_frame_name, switch_name = OBSERVED_PANELS[field]
    _attach_observed_data(gmg_module, gmg_frame, field)

    # THE USER CLOSES THE PANEL, THEN SAVES
    setattr(gmg_frame, switch_name, False)
    gmg_frame.save_model(None)

    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    # THE SAVED LAYOUT IS HONOURED: THE PANEL COMES BACK CLOSED, HOLDING ITS DATA
    assert getattr(gmg_frame, switch_name) is False
    assert getattr(gmg_frame, frame_name) is None
    assert len(getattr(gmg_frame, list_name)) == 1

    # AND THE DATA IS DRAWN THE MOMENT THE PANEL IS REOPENED
    gmg_frame.ensure_data_frame_shown(field)
    assert getattr(gmg_frame, frame_name) is not None
    assert len(_scatters_in(getattr(gmg_frame, frame_name))) == 1


@pytest.mark.parametrize('field', sorted(OBSERVED_PANELS))
def test_loading_draws_each_observed_dataset_once(gmg_module, gmg_frame, field):
    """
    Test loading draws one artist per observed dataset, not two

    rebuild_data_frames already replots the observed data into the panels it recreates. load_model then
    called replot_observed_*_data, which drew a second set over the top and overwrote mpl_actor, so the
    first set was orphaned but still on the axes and 'delete observed data' left it behind.
    """
    _list_name, frame_name, _d_frame_name, _switch_name = OBSERVED_PANELS[field]
    _attach_observed_data(gmg_module, gmg_frame, field)

    gmg_frame.save_model(None)
    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    frame = getattr(gmg_frame, frame_name)
    assert len(_scatters_in(frame)) == 1


@pytest.mark.parametrize('field', sorted(OBSERVED_PANELS))
def test_loading_puts_derivative_data_in_the_d_dx_axis(gmg_module, gmg_frame, field):
    """
    Test a horizontal derivative loads into the d/dx twin axis, not over the observed data

    replot_observed_*_data had no type check and always used the main axes, and ran last, so every
    saved derivative came back plotted against the wrong y scale.
    """
    _list_name, frame_name, d_frame_name, _switch_name = OBSERVED_PANELS[field]
    _attach_observed_data(gmg_module, gmg_frame, field, data_type='derivative')

    gmg_frame.save_model(None)
    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    assert len(_scatters_in(getattr(gmg_frame, d_frame_name))) == 1
    assert len(_scatters_in(getattr(gmg_frame, frame_name))) == 0


def test_modify_model_dimensions_carries_the_observation_profiles_over(gmg_module, gmg_frame, monkeypatch):
    """
    Test changing the model dimensions resamples the observation elevations onto the new sample points

    An observation level read from a profile holds one value per calculation point. The three lines that
    rebuilt them were commented out, so changing the dimensions left them at their old length and every
    later calculation raised 'operands could not be broadcast together'. A model saved in that state
    could not be reopened at all.
    """
    previous_xp = gmg_frame.xp.copy()
    gmg_frame.mag_observation_elv = np.full_like(previous_xp, -4000.0)
    gmg_frame.gravity_observation_elv = (previous_xp * 0.01).astype('f')  # VARYING, SO INTERPOLATION IS TESTED

    # THE USER WIDENS THE MODEL, WHICH LENGTHENS xp
    monkeypatch.setattr(gmg_module, 'NewModelDialog',
                        _new_model_dialog(MODEL_X1, MODEL_X2 * 2, MODEL_Z1, MODEL_Z2))
    gmg_frame.modify_model_dimensions(None)

    assert len(gmg_frame.xp) > len(previous_xp)
    assert gmg_frame.mag_observation_elv.shape == gmg_frame.xp.shape
    assert gmg_frame.gravity_observation_elv.shape == gmg_frame.xp.shape

    # A CONSTANT LEVEL IS UNCHANGED, AND A VARYING PROFILE IS CARRIED OVER RATHER THAN ZEROED
    np.testing.assert_allclose(gmg_frame.mag_observation_elv, -4000.0)
    overlap = gmg_frame.xp <= previous_xp[-1]
    np.testing.assert_allclose(gmg_frame.gravity_observation_elv[overlap],
                               gmg_frame.xp[overlap] * 0.01, atol=1.0e-3)


def test_model_with_observation_elevations_of_the_wrong_length_loads(gmg_frame, model_path, reported):
    """
    Test a model whose observation elevations do not match its calculation points is repaired on load

    gmg wrote files in this state for as long as modify_model_dimensions left the profiles behind, and
    every attempt to reopen one raised ValueError from the forward calculation.
    """
    gmg_frame.mag_observation_elv = np.full_like(gmg_frame.xp, -4000.0)
    gmg_frame.save_model(None)

    # IMITATE A FILE WRITTEN AFTER THE DIMENSIONS WERE CHANGED BUT BEFORE THE PROFILES WERE RESAMPLED
    with open(str(model_path), 'rb') as saved_file:
        saved = pickle.load(saved_file)
    saved['mag_observation_elv'] = saved['mag_observation_elv'][::2]
    with open(str(model_path), 'wb') as saved_file:
        pickle.dump(saved, saved_file, protocol=pickle.HIGHEST_PROTOCOL)

    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    assert gmg_frame.mag_observation_elv.shape == gmg_frame.xp.shape
    np.testing.assert_allclose(gmg_frame.mag_observation_elv, -4000.0)

    # THE REPAIR IS REPORTED RATHER THAN MADE QUIETLY, SINCE IT CHANGES WHAT THE MODEL CALCULATES FROM
    assert reported.said('do not match its')


def test_a_failed_load_reports_the_fault_and_leaves_gmg_usable(gmg_module, gmg_frame, reported, monkeypatch):
    """
    Test a load that fails part way reports it and leaves a model that can still be clicked on

    load_model wipes the current model before it reads the file and used to catch only IOError and
    IndexError. Any other fault escaped into the wx event handler, which swallows it, leaving a frame
    with no currently_active_layer - so the next click raised AttributeError as well.
    """
    from matplotlib.backend_bases import MouseButton, MouseEvent

    gmg_frame.save_model(None)
    gmg_frame.model_saved = True

    def _fail_to_read(*args, **kwargs):
        raise RuntimeError('deliberate read failure')

    monkeypatch.setattr(gmg_module.pickle, 'load', _fail_to_read)

    with pytest.raises(RuntimeError):
        gmg_frame.load_model(None)

    assert reported.said('MODEL NOT LOADED')
    assert reported.said('deliberate read failure')

    # AN EMPTY BUT WORKING MODEL IS LEFT BEHIND, AND A CLICK IN THE MODEL FRAME NO LONGER RAISES
    assert len(gmg_frame.layer_list) == 1
    event = MouseEvent('button_press_event', gmg_frame.fig.canvas, 100, 100, button=MouseButton.LEFT)
    event.inaxes = gmg_frame.model_frame
    event.xdata, event.ydata = 50.0, 20.0
    gmg_frame.button_press(event)


@pytest.mark.parametrize('field', sorted(OBSERVED_PANELS))
def test_repeated_loads_do_not_pile_up_menu_entries(gmg_module, gmg_frame, field):
    """
    Test loading a model several times leaves one menu entry per dataset, not one more each time

    The menubar is built in Gmg.__init__ and is not rebuilt by start(), so the data submenus outlive the
    model they were filled for. Every replot_* function appended without clearing first.
    """
    menus = {'topography': 'm_topo_submenu', 'gravity': 'm_obs_g_submenu',
             'vgg': 'm_obs_vgg_submenu', 'magnetic': 'm_obs_mag_submenu'}

    _attach_observed_data(gmg_module, gmg_frame, field)
    expected = _menu_labels(getattr(gmg_frame, menus[field]))
    assert len(expected) == 1

    gmg_frame.save_model(None)
    for _ in range(3):
        gmg_frame.model_saved = True
        gmg_frame.load_model(None)
        assert _menu_labels(getattr(gmg_frame, menus[field])) == expected


def test_loading_a_png_draws_it_immediately_and_lists_it(gmg_module, gmg_frame, tmp_path, monkeypatch, reported):
    """
    Test loading a PNG draws it there and then, and puts its name in the PNG menu

    load_png plotted the image but never repainted the canvas, so it only appeared once something else
    redrew - in practice when the user clicked the model. Everything that registered the PNG was
    commented out, so the 'PNGs...' menu stayed empty and there was no way to hide or delete one.
    """
    drawn = []
    monkeypatch.setattr(gmg_frame, 'draw', lambda *args, **kwargs: drawn.append(1))

    _load_png(gmg_module, gmg_frame, monkeypatch, _write_png(tmp_path), name='seismic line')

    assert not reported.said('PNG LOAD ERROR')

    # THE IMAGE IS ON THE MODEL FRAME, AND THE CANVAS WAS REPAINTED WITHOUT WAITING FOR A CLICK
    assert len(gmg_frame.model_frame.images) == 1
    assert drawn

    # AND IT IS LISTED, WITH THE NAME THE DIALOG WAS GIVEN
    assert len(gmg_frame.png_data_list) == 1
    assert gmg_frame.png_data_list[0].name == 'seismic line'
    assert _menu_labels(gmg_frame.m_png_submenu) == ['seismic line']


def test_an_unnamed_png_still_gets_a_menu_entry(gmg_module, gmg_frame, tmp_path, monkeypatch):
    """
    Test a PNG loaded with the name box left empty is listed under a name of gmg's own
    """
    _load_png(gmg_module, gmg_frame, monkeypatch, _write_png(tmp_path), name='')

    assert _menu_labels(gmg_frame.m_png_submenu) == ['PNG 0']


def test_png_hide_show_toggles_the_image(gmg_module, gmg_frame, tmp_path, monkeypatch):
    """
    Test the PNG Hide/Show menu item hides the image and brings it back
    """
    _load_png(gmg_module, gmg_frame, monkeypatch, _write_png(tmp_path))
    png = gmg_frame.png_data_list[0]
    assert png.mpl_actor.get_visible()

    _fire_menu(gmg_frame, png.id + 14500)
    assert not png.mpl_actor.get_visible()

    _fire_menu(gmg_frame, png.id + 14500)
    assert png.mpl_actor.get_visible()


def test_deleting_a_png_clears_the_image_and_the_menu(gmg_module, gmg_frame, tmp_path, monkeypatch):
    """
    Test 'delete png' takes the image off the model frame and the entry out of the menu
    """
    _load_png(gmg_module, gmg_frame, monkeypatch, _write_png(tmp_path))
    png_id = gmg_frame.png_data_list[0].id

    _fire_menu(gmg_frame, png_id + 14000)

    assert len(gmg_frame.model_frame.images) == 0
    assert gmg_frame.png_data_list[png_id] is None
    assert _menu_labels(gmg_frame.m_png_submenu) == []


def test_png_survives_a_save_and_load(gmg_module, gmg_frame, tmp_path, model_path, monkeypatch):
    """
    Test a PNG saved in a model is drawn and listed again on load, and does not double on a second load

    load_model held 'pass  # self.replot_png_data()', so a saved PNG simply vanished. The equivalent
    segy code appends to the list it is looping over, which doubles every entry; replot_png_data must
    not do the same.
    """
    _load_png(gmg_module, gmg_frame, monkeypatch, _write_png(tmp_path), name='saved figure')

    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(model_path))
    gmg_frame.save_model(None)

    for _ in range(2):
        gmg_frame.model_saved = True
        gmg_frame.load_model(None)

        assert len(gmg_frame.png_data_list) == 1
        assert len(gmg_frame.model_frame.images) == 1
        assert _menu_labels(gmg_frame.m_png_submenu) == ['saved figure']


def test_png_input_returns_quietly_when_the_dialog_is_dismissed(gmg_module, gmg_frame, monkeypatch):
    """
    Test dismissing the dimensions dialog without setting coordinates does nothing, rather than raising

    SeisDialog only assigns dimensions and segy_name_input inside its 'Set Coordinates' handler, so
    png_input raised AttributeError on the next line for any other way of closing it.
    """
    class _DismissedDialog(_AcceptingDialog):
        """CLOSED WITHOUT PRESSING 'Set Coordinates', SO IT CARRIES NO dimensions"""

    monkeypatch.setattr(gmg_module, 'SeisDialog', _DismissedDialog)

    gmg_frame.png_input(None)

    assert gmg_frame.png_data_list == []


def test_segy_survives_a_save_and_load_without_doubling(gmg_module, gmg_frame, tmp_path, model_path, monkeypatch):
    """
    Test a segy saved in a model comes back once per load, however many times it is loaded

    replot_segy_data appended to self.segy_data_list while looping over range(len()) of that same list,
    so every load duplicated each segy - and its menu entry with it.
    """
    _load_segy(gmg_module, gmg_frame, monkeypatch, _write_segy(tmp_path), name='line 2')
    assert len(gmg_frame.segy_data_list) == 1

    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(model_path))
    gmg_frame.save_model(None)

    for _ in range(3):
        gmg_frame.model_saved = True
        gmg_frame.load_model(None)

        assert len(gmg_frame.segy_data_list) == 1
        assert len(gmg_frame.model_frame.images) == 1
        assert _menu_labels(gmg_frame.m_segy_submenu) == ['line 2']


def test_a_segy_that_will_not_read_does_not_stop_the_others(gmg_module, gmg_frame, tmp_path, model_path,
                                                            monkeypatch, reported):
    """
    Test one unreadable segy is reported and skipped, leaving the rest of the model's segys plotted

    replot_segy_data used to return out of the loop, so every segy after the bad one was abandoned.
    """
    _load_segy(gmg_module, gmg_frame, monkeypatch, _write_segy(tmp_path, 'first.segy'), name='first')
    second = _write_segy(tmp_path, 'second.segy')
    _load_segy(gmg_module, gmg_frame, monkeypatch, second, name='second')
    assert len(gmg_frame.model_frame.images) == 2

    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(model_path))
    gmg_frame.save_model(None)

    # THE FIRST FILE GOES MISSING, AS IT WOULD IF THE MODEL WERE MOVED TO ANOTHER MACHINE
    os.remove(str(tmp_path / 'first.segy'))

    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    assert reported.said('SEGY LOAD ERROR')
    assert _menu_labels(gmg_frame.m_segy_submenu) == ['second']
    assert len(gmg_frame.model_frame.images) == 1


def test_grid_survives_a_save_and_load(gmg_module, gmg_frame, tmp_path, model_path, monkeypatch):
    """
    Test a grid saved in a model is drawn and listed again, and can then be deleted

    load_model held 'pass  # self.replot_grid_data()' and the function did not exist, so a saved grid
    simply vanished even though grid_data_list was being written to the file.
    """
    _load_grid(gmg_module, gmg_frame, monkeypatch, _write_grid(tmp_path), name='velocity model')

    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(model_path))
    gmg_frame.save_model(None)

    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    assert len(gmg_frame.model_frame.images) == 1
    assert _menu_labels(gmg_frame.m_grid_submenu) == ['velocity model']

    # AND THE MENU ENTRY STILL WORKS AFTER THE LOAD
    _fire_menu(gmg_frame, gmg_frame.grid_data_list[0].id + 5000)
    assert gmg_frame.grid_data_list[0] is None
    assert len(gmg_frame.model_frame.images) == 0
    assert _menu_labels(gmg_frame.m_grid_submenu) == []


def test_vgg_filter_and_derivative_are_listed_under_vgg(gmg_module, gmg_frame, monkeypatch):
    """
    Test VGG datasets go in the VGG menu, not the gravity one

    filter_observed_vgg and take_vgg_horizontal_derivative appended to m_obs_g_submenu while binding
    delete_obs_vgg, which calls DestroyItem on the VGG menu - so the entry showed up under Gravity and
    could never be removed.
    """
    _attach_observed_data(gmg_module, gmg_frame, 'vgg')
    gravity_before = _menu_labels(gmg_frame.m_obs_g_submenu)

    monkeypatch.setattr(gmg_module, 'GaussianFilterDialog',
                        _filter_dialog('vgg filtered', gmg_frame.observed_vgg_list[0].data))
    gmg_frame.filter_observed_vgg(None)

    monkeypatch.setattr(gmg_module, 'HorizontalDerivative',
                        _derivative_dialog('vgg derivative', gmg_frame.observed_vgg_list[0].data))
    gmg_frame.take_vgg_horizontal_derivative(None)

    assert _menu_labels(gmg_frame.m_obs_vgg_submenu) == ['vgg test data', 'vgg filtered', 'vgg derivative']
    assert _menu_labels(gmg_frame.m_obs_g_submenu) == gravity_before

    # AND BOTH DELETE CLEANLY, WHICH THEY COULD NOT WHILE THEY SAT IN THE WRONG MENU
    for observed in list(gmg_frame.observed_vgg_list[1:]):
        _fire_menu(gmg_frame, observed.id + 15000)
    assert _menu_labels(gmg_frame.m_obs_vgg_submenu) == ['vgg test data']


def test_deleting_an_observed_vgg_dataset_after_a_load_does_not_raise(gmg_module, gmg_frame):
    """
    Test the VGG delete item is bound to delete_obs_vgg, not delete_obs_grav

    replot_observed_vgg_data bound delete_obs_grav for a 15000 ID. That subtracts 11000, so deleting an
    observed VGG dataset after a load indexed observed_gravity_list[4000 + id] and raised IndexError.
    """
    observed = _attach_observed_data(gmg_module, gmg_frame, 'vgg')

    gmg_frame.save_model(None)
    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    _fire_menu(gmg_frame, gmg_frame.observed_vgg_list[0].id + 15000)

    assert gmg_frame.observed_vgg_list[0] is None
    assert _menu_labels(gmg_frame.m_obs_vgg_submenu) == []


@pytest.mark.parametrize('field', sorted(OBSERVED_PANELS))
def test_a_derivative_can_still_be_deleted_after_a_load(gmg_module, gmg_frame, field):
    """
    Test a horizontal derivative keeps a working delete item across a save and load

    The derivative handlers stored object.id as base + counter - the menu ID rather than the object ID -
    while every other handler stored the bare counter, so replot_* computed base + base + counter and
    the delete item stopped matching its dataset.
    """
    base = dict(gmg_module.Gmg.OBSERVED_DATA_MENU_BASES)[OBSERVED_PANELS[field][0]]

    _attach_observed_data(gmg_module, gmg_frame, field, data_type='derivative')

    gmg_frame.save_model(None)
    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    data_list = getattr(gmg_frame, OBSERVED_PANELS[field][0])
    assert data_list[0].id == 0  # THE BARE COUNTER, NOT base + IT

    _fire_menu(gmg_frame, data_list[0].id + base)
    assert data_list[0] is None


@pytest.mark.parametrize('field', sorted(OBSERVED_PANELS))
def test_a_model_saved_with_offset_ids_is_repaired_on_load(gmg_module, gmg_frame, model_path, field):
    """
    Test a .model file written when derivative IDs carried their menu base still deletes correctly

    Those files are already on disc, so load_model takes the offset back off.
    """
    list_name = OBSERVED_PANELS[field][0]
    base = dict(gmg_module.Gmg.OBSERVED_DATA_MENU_BASES)[list_name]

    _attach_observed_data(gmg_module, gmg_frame, field, data_type='derivative')
    gmg_frame.save_model(None)

    # IMITATE A FILE WRITTEN BEFORE THE IDs WERE NORMALISED
    with open(str(model_path), 'rb') as saved_file:
        saved = pickle.load(saved_file)
    for observed in saved[list_name]:
        observed.id += base
    with open(str(model_path), 'wb') as saved_file:
        pickle.dump(saved, saved_file, protocol=pickle.HIGHEST_PROTOCOL)

    gmg_frame.model_saved = True
    gmg_frame.load_model(None)

    data_list = getattr(gmg_frame, list_name)
    assert data_list[0].id == 0

    _fire_menu(gmg_frame, data_list[0].id + base)
    assert data_list[0] is None


def test_taking_a_topography_derivative_does_not_raise(gmg_module, gmg_frame, monkeypatch):
    """
    Test the topography derivative binds a handler that exists

    It bound self.delete_obs_topo, which is not defined anywhere, so the Bind raised AttributeError -
    and unlike its three siblings this handler has no try, so it never completed.
    """
    observed = _attach_observed_data(gmg_module, gmg_frame, 'topography')

    monkeypatch.setattr(gmg_module, 'HorizontalDerivative',
                        _derivative_dialog('topo derivative', observed.data))
    gmg_frame.take_topography_horizontal_derivative(None)

    assert _menu_labels(gmg_frame.m_topo_submenu) == ['topography test data', 'topo derivative']


def test_the_image_loaders_do_not_share_menu_ids(gmg_module, gmg_frame, tmp_path, monkeypatch):
    """
    Test deleting a grid leaves a loaded PNG alone

    Every image loader bound `counter + 1000`. wx binds menu IDs on the frame, so the later Bind won
    and one loader's delete item called another loader's handler.
    """
    _load_png(gmg_module, gmg_frame, monkeypatch, _write_png(tmp_path))
    _load_grid(gmg_module, gmg_frame, monkeypatch, _write_grid(tmp_path))

    assert len(gmg_frame.model_frame.images) == 2

    # DELETE THE GRID THROUGH ITS MENU ITEM, WHICH IS WHERE THE CLASH LIVED
    _fire_menu(gmg_frame, gmg_frame.grid_data_list[0].id + 5000)

    assert gmg_frame.grid_data_list[0] is None
    assert gmg_frame.png_data_list[0] is not None
    assert len(gmg_frame.model_frame.images) == 1


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


def test_dc_shift_defaults_off_and_raw(gmg_frame):
    """
    Test the DC shift is off by default and the plotted curve is the raw forward calculation
    """
    assert (gmg_frame.gravity_dc_mode, gmg_frame.vgg_dc_mode, gmg_frame.magnetic_dc_mode) == ('off',) * 3
    assert (gmg_frame.gravity_dc_value, gmg_frame.vgg_dc_value, gmg_frame.magnetic_dc_value) == (0.0,) * 3

    # NB: A LOOSE TOLERANCE, BECAUSE matplotlib STORES LINE DATA AS float32
    np.testing.assert_allclose(gmg_frame.pred_gravity_plot.get_ydata(), gmg_frame.predicted_gravity,
                               rtol=1.0e-6)
    assert 'DC' not in gmg_frame.gravity_frame.get_ylabel()


def test_dc_shift_auto_zeroes_the_mean_residual(gmg_frame, reported, tmp_path, monkeypatch):
    """
    Test auto mode fits the constant that zeroes the mean residual against the RMS input

    An observed dataset that is the calculated curve moved up by 5 must produce a DC of exactly +5, an RMS
    of zero, and a plotted (and exported) curve lying on the observed data. The raw predicted values must
    never change: the shift is applied at the points of consumption only.
    """
    raw = gmg_frame.predicted_gravity.copy()
    gmg_frame.obs_gravity_data_for_rms = np.column_stack((gmg_frame.xp * 0.001, raw + 5.0))
    gmg_frame.gravity_dc_mode = 'auto'
    gmg_frame.run_algorithms()

    assert abs(gmg_frame.gravity_dc_value - 5.0) < 1.0e-6
    assert gmg_frame.gravity_rms_value < 1.0e-4  # rms() ROUNDS TO 4 D.P.
    assert 'DC +5.00' in gmg_frame.gravity_frame.get_ylabel()
    np.testing.assert_array_equal(gmg_frame.predicted_gravity, raw)
    np.testing.assert_allclose(gmg_frame.pred_gravity_plot.get_ydata(), raw + 5.0, rtol=1.0e-6)

    # THE EXPORT IS WYSIWYG, SO IT MUST CONTAIN THE SHIFTED VALUES
    target = tmp_path / 'predicted_gravity.txt'
    monkeypatch.setattr(wx, 'FileDialog', _file_dialog_returning(target))
    gmg_frame.save_modelled_grav(None)
    exported = np.loadtxt(str(target))
    np.testing.assert_allclose(exported[:, 1], raw + 5.0, atol=1.0e-6)

    # SWITCHING BACK OFF MUST RESTORE THE RAW CURVE EXACTLY
    gmg_frame.gravity_dc_mode = 'off'
    gmg_frame.run_algorithms()
    assert gmg_frame.gravity_dc_value == 0.0
    np.testing.assert_allclose(gmg_frame.pred_gravity_plot.get_ydata(), raw, rtol=1.0e-6)


def test_dc_shift_manual_moves_curve_and_unreliable_markers(gmg_frame):
    """
    Test manual mode shifts the magnetic curve and the unreliable sample markers ride on it
    """
    raw = gmg_frame.predicted_nt.copy()
    gmg_frame.magnetic_dc_mode = 'manual'
    gmg_frame.magnetic_dc_manual = -12.5
    gmg_frame.run_algorithms()

    assert gmg_frame.magnetic_dc_value == -12.5
    assert 'DC -12.50' in gmg_frame.magnetic_frame.get_ylabel()
    np.testing.assert_allclose(gmg_frame.predicted_nt_plot.get_ydata(), raw - 12.5, rtol=1.0e-6)

    # RAISE THE LAYER THROUGH THE CALCULATION ELEVATION SO SOME SAMPLES BECOME UNRELIABLE, AND CHECK
    # THEIR MARKERS SIT ON THE SHIFTED CURVE (SEE test_interior_observation_points_are_marked)
    layer = gmg_frame.layer_list[gmg_frame.currently_active_layer_id]
    layer.y_nodes[:] = [-2.0, -2.0, 15.0, 15.0]
    gmg_frame.update_layer_data()
    gmg_frame.run_algorithms()

    inside = gmg_frame.unreliable_magnetic_observation_points()
    assert inside.any()
    np.testing.assert_allclose(gmg_frame.untrusted_nt_plot.get_ydata(),
                               (gmg_frame.predicted_nt + gmg_frame.magnetic_dc_value)[inside], rtol=1.0e-6)


def test_dc_shift_auto_with_nothing_to_fit_behaves_as_off(gmg_frame):
    """
    Test auto mode with no RMS input, or with the calculation switched off, applies no shift

    Fitting against an unset dataset is impossible, and fitting against the zero curve of a switched off
    calculation would produce a spurious mean(obs) shift.
    """
    gmg_frame.gravity_dc_mode = 'auto'
    gmg_frame.run_algorithms()
    assert gmg_frame.gravity_dc_value == 0.0

    raw = gmg_frame.predicted_gravity.copy()
    gmg_frame.obs_gravity_data_for_rms = np.column_stack((gmg_frame.xp * 0.001, raw + 5.0))
    gmg_frame.calc_grav_switch = False
    gmg_frame.run_algorithms()
    assert gmg_frame.gravity_dc_value == 0.0


def test_dc_shift_round_trips_and_old_files_default_off(gmg_frame, model_path):
    """
    Test the DC settings survive a save and load, and files saved before the feature load with it off
    """
    gmg_frame.gravity_dc_mode = 'manual'
    gmg_frame.gravity_dc_manual = 3.25
    gmg_frame.magnetic_dc_mode = 'auto'
    gmg_frame.run_algorithms()
    gmg_frame.save_model(None)

    gmg_frame.model_saved = True
    gmg_frame.load_model(None)
    assert (gmg_frame.gravity_dc_mode, gmg_frame.gravity_dc_manual) == ('manual', 3.25)
    assert gmg_frame.magnetic_dc_mode == 'auto'
    assert gmg_frame.gravity_dc_value == 3.25  # REDERIVED BY run_algorithms DURING THE LOAD

    # STRIP THE DC KEYS TO IMITATE A FILE WRITTEN BEFORE THE FEATURE EXISTED
    with open(str(model_path), 'rb') as saved_file:
        saved = pickle.load(saved_file)
    for key in ('gravity_dc_mode', 'gravity_dc_manual', 'vgg_dc_mode', 'vgg_dc_manual',
                'magnetic_dc_mode', 'magnetic_dc_manual'):
        del saved[key]
    with open(str(model_path), 'wb') as saved_file:
        pickle.dump(saved, saved_file, protocol=pickle.HIGHEST_PROTOCOL)

    gmg_frame.model_saved = True
    gmg_frame.load_model(None)
    assert (gmg_frame.gravity_dc_mode, gmg_frame.vgg_dc_mode, gmg_frame.magnetic_dc_mode) == ('off',) * 3
    assert (gmg_frame.gravity_dc_manual, gmg_frame.gravity_dc_value) == (0.0, 0.0)


def test_dc_shift_dialog_applies_and_warns_without_rms_input(gmg_frame, gmg_module, monkeypatch):
    """
    Test the DC Shift dialog handler stores the answer, and warns when auto has nothing to fit to
    """
    class _DCShift(_AcceptingDialog):
        mode, manual_value = 'auto', 0.0

        def ShowModal(self):
            return 1

    warnings = []
    monkeypatch.setattr(gmg_module, 'DCShiftDialog', _DCShift)
    monkeypatch.setattr(wx, 'MessageBox', lambda text, *args, **kwargs: warnings.append(text) or wx.YES)

    gmg_frame.set_gravity_dc_shift(None)

    assert gmg_frame.gravity_dc_mode == 'auto'
    assert gmg_frame.gravity_dc_value == 0.0  # NO RMS INPUT, SO NOTHING IS APPLIED
    assert any('RMS input' in text for text in warnings)

    # WITH AN RMS INPUT SET THERE IS NO WARNING AND THE SHIFT IS APPLIED
    warnings.clear()
    raw = gmg_frame.predicted_gravity.copy()
    gmg_frame.obs_gravity_data_for_rms = np.column_stack((gmg_frame.xp * 0.001, raw - 2.0))
    gmg_frame.set_gravity_dc_shift(None)
    assert not warnings
    assert abs(gmg_frame.gravity_dc_value - (-2.0)) < 1.0e-6
