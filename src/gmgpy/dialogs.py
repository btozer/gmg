"""
GMG EXTERNAL DIALOG BOXES GO HERE

EACH DIALOG BOX IS DEFINED AS A CLASS
"""

import wx
import matplotlib
matplotlib.use('WXAgg')
import numpy as np
from scipy import signal
from scipy import interpolate as ip
from scipy.ndimage import gaussian_filter1d

class NewModelDialog(wx.Dialog):
    """
    CREATE A NEW MODEL FRAME. RETURNS MODEL PARAMETERS AND CALCULATION INCREMENT

    **Paramaters**

    * User input : Wx.Dialog entries

    Returns:

        * x1 : float
        The start x value for the model frame
        * x2 : float
        The end x value for the model frame
        * z1 : float
        The start z value for the model frame
        * z2 : float
        The end z value for the model frame
        * xp_inc : float
        The increment at which predicted anomalies are calculated
    """

    def __init__(self, parent, id, title, m_x1=None, m_x2=None, m_z1=None, m_z2=None):
        """DIALOG BOX USED TO GATHER USER INPUT FOR CREATING A NEW MODEL"""
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.MAXIMIZE_BOX | wx.OK | wx.CANCEL
                                                          | wx.BORDER_RAISED)

        self.floating_panel = wx.Panel(self, -1)
        self.set_button = False

        self.main_box = wx.BoxSizer(wx.HORIZONTAL)

        self.title_text = wx.StaticText(self.floating_panel, -1, "Enter the region boundary for the new model\n"
                                                                 "and the x-increment for anomaly calculations:",
                                        style=wx.ALIGN_CENTRE_HORIZONTAL)
        self.line1 = (wx.StaticLine(self.floating_panel), 0, wx.ALL | wx.EXPAND, 5)
        self.x_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.x_Text = wx.StaticText(self.floating_panel, -1, "X1, X2 (km)      ")
        self.new_x1 = wx.TextCtrl(self.floating_panel, -1, "0", size=(80, -1))
        self.new_x2 = wx.TextCtrl(self.floating_panel, -1, "0", size=(80, -1))
        self.x_sizer.AddMany([self.x_Text, self.new_x1, self.new_x2])
       
        self.z_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.z_Text = wx.StaticText(self.floating_panel, -1, "Z1, Z2 (km)      ")
        self.new_z1 = wx.TextCtrl(self.floating_panel, -1, "0", size=(80, -1))
        self.new_z2 = wx.TextCtrl(self.floating_panel, -1, "0", size=(80, -1))
        self.z_sizer.AddMany([self.z_Text, self.new_z1, self.new_z2])
       
        self.xp1_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.xp1_text = wx.StaticText(self.floating_panel, -1, "Calculation       \nSpacing (km)  ")
        self.xp1_inc = wx.TextCtrl(self.floating_panel, -1, "0", size=(80, -1))
        self.xp1_sizer.AddMany([self.xp1_text, self.xp1_inc])

        self.line2 = (wx.StaticLine(self.floating_panel), 0, wx.ALL | wx.EXPAND, 5)

        self.create_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.spacing = wx.StaticText(self.floating_panel, -1, "                           ")
        self.b_create_button = wx.Button(self.floating_panel, -1, "Create New Model")
        self.Bind(wx.EVT_BUTTON, self.create_button, self.b_create_button)
        self.create_button_sizer.AddMany([self.spacing, self.b_create_button])

        self.line3 = (wx.StaticLine(self.floating_panel), 0, wx.ALL | wx.EXPAND, 5)

        self.footer_text = wx.StaticText(self.floating_panel, -1,
                                         "NB. these values can be modified during modelling,\n"
                                         "it may be beneficial to begin with a coarse spacing.")
       
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
        """WHEN THE "Create Model" BUTTON IS PRESSED"""
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

        # SET FILE PATH IN CLASS PANEL
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
        elif self.type == 'vgg':
            self.parent.open_obs_vgg()
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


class LayerNameDialog(wx.Dialog):
    """SET MAGNETIC FIELD PARAMETERS"""

    def __init__(self, parent, id, title, current_name):
        wx.Dialog.__init__(self, parent, id, 'Rename layer', style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                                   | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        # OBSERVATION ELEVATION
        self.set_name = wx.StaticText(input_panel, -1, "New layer name:")
        self.set_name_text = wx.TextCtrl(input_panel, -1, current_name, size=(75, -1))
        self.set_name_text.SetInsertionPoint(0)

        # SET BUTTON
        self.b_set_button = wx.Button(input_panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.set_button, self.b_set_button)

        sizer = wx.FlexGridSizer(cols=2, hgap=7, vgap=7)
        sizer.AddMany([self.set_name, self.set_name_text, self.b_set_button])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_button(self, event):
        self.name = str(self.set_name_text.GetValue())
        self.EndModal(1)


class MagDialog(wx.Dialog):
    """SET MAGNETIC FIELD PARAMETERS"""

    def __init__(self, parent, id, title, mag_observation_elv):
        wx.Dialog.__init__(self, parent, id, 'Set Magnetic Field', style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                                         | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        # OBSERVATION ELEVATION
        self.set_mag_observation_elv = wx.StaticText(input_panel, -1, "Observation\nElevation:")
        self.mag_observation_elv_text = wx.TextCtrl(input_panel, -1, str(mag_observation_elv * 0.001), size=(75, -1))
        self.mag_observation_elv_units = wx.StaticText(input_panel, -1, "km")
        self.mag_observation_elv_text.SetInsertionPoint(0)

        # CREATE BOX SIZER
        sizer = wx.FlexGridSizer(cols=3, hgap=7, vgap=7)
        self.b_set_button_mag = wx.Button(input_panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.set_button_mag_elv, self.b_set_button_mag)

        # ADD ITEMS TO SIZER
        sizer.AddMany([self.set_mag_observation_elv, self.mag_observation_elv_text, self.mag_observation_elv_units,
                       self.b_set_button_mag])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_button_mag_elv(self, event):
        self.mag_observation_elv = float(self.mag_observation_elv_text.GetValue())
        self.EndModal(1)


class GravDialog(wx.Dialog):
    """POPOUT BOX TO LET THE USER SET THE ELEVATION AT WHICH GRAVITY ANOMALIES ARE CALCULATED"""

    def __init__(self, parent, id, title, grav_observation_elv):
        wx.Dialog.__init__(self, parent, id, 'Set Gravity elevation', style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                                            | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)

        # OBSERVATION ELEVATION
        self.set_grav_observation_elv = wx.StaticText(input_panel, -1, "Observation\nElevation:")
        self.grav_observation_elv_text = wx.TextCtrl(input_panel, -1, str(grav_observation_elv * 0.001), size=(75, -1))
        self.grav_observation_elv_units = wx.StaticText(input_panel, -1, "km")
        self.grav_observation_elv_text.SetInsertionPoint(0)

        # CREATE BOX SIZER
        sizer = wx.FlexGridSizer(cols=3, hgap=7, vgap=7)
        self.b_set_button_grav_elv = wx.Button(input_panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.set_button_grav_elv, self.b_set_button_grav_elv)

        # ADD ITEMS TO SIZER
        sizer.AddMany([self.set_grav_observation_elv, self.grav_observation_elv_text, self.grav_observation_elv_units,
                       self.b_set_button_grav_elv])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_button_grav_elv(self, event):
        """ON BUTTON PRESS - SET THE NEW GRAVITY ELEVATION"""
        self.grav_observation_elv = float(self.grav_observation_elv_text.GetValue())
        self.EndModal(1)


class PinchDialog(wx.Dialog):
    """PINCH A LAYER TO ANOTHER LAYER"""

    def __init__(self, parent, id, title, layer_list, currently_active_layer_id):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)
        self.layer_list = layer_list
        self.currently_active_layer_id = currently_active_layer_id
        self.p_start = wx.StaticText(input_panel, -1, "Pinch from (km):")
        self.p_start_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.p_start_Text.SetInsertionPoint(0)
        self.p_end = wx.StaticText(input_panel, -1, "Pinch to (km):")
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

        current_x = self.layer_list[self.currently_active_layer_id].x_nodes
        current_y = self.layer_list[self.currently_active_layer_id].y_nodes
        above_x = self.layer_list[self.currently_active_layer_id - 1].x_nodes
        above_y = self.layer_list[self.currently_active_layer_id - 1].y_nodes

        # SET ID FOR LAYER GETTING PINCHED
        self.layer_getting_pinched = self.currently_active_layer_id - 1

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

        # CLOSE POPOUT BOX
        self.EndModal(1)

    def down_set_button(self, event):
        p_start = float(self.p_start_Text.GetValue())
        p_end = float(self.p_end_Text.GetValue())

        current_x = self.layer_list[self.currently_active_layer_id].x_nodes
        current_y = self.layer_list[self.currently_active_layer_id].y_nodes
        below_x = self.layer_list[self.currently_active_layer_id + 1].x_nodes
        below_y = self.layer_list[self.currently_active_layer_id + 1].y_nodes

        # SET ID FOR LAYER GETTING PINCHED
        self.layer_getting_pinched = self.currently_active_layer_id + 1

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

        # CLOSE POPOUT BOX
        self.EndModal(1)


class DepinchDialog(wx.Dialog):
    """DEPINCH A LAYER"""

    def __init__(self, parent, id, title, layer_list, currently_active_layer_id, total_layer_count):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)
        self.layer_list = layer_list
        self.currently_active_layer_id = currently_active_layer_id
        self.total_layer_count = total_layer_count
        self.p_start = wx.StaticText(input_panel, -1, "Depinch from (km):")
        self.p_start_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.p_start_Text.SetInsertionPoint(0)
        self.p_end = wx.StaticText(input_panel, -1, "Depinch to (km):")
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

        current_x = self.layer_list[self.currently_active_layer_id].x_nodes
        current_y = self.layer_list[self.currently_active_layer_id].y_nodes

        if self.currently_active_layer_id != 0:
            above_x = self.layer_list[self.currently_active_layer_id - 1].x_nodes
            above_y = self.layer_list[self.currently_active_layer_id - 1].y_nodes
        if self.currently_active_layer_id != self.total_layer_count:
            below_x = self.layer_list[self.currently_active_layer_id + 1].x_nodes
            below_y = self.layer_list[self.currently_active_layer_id + 1].y_nodes

        # REMOVE PINCHED NODES
        self.depinched_x = [tup for i, tup in enumerate(current_x) if current_x[i] < p_start or current_x[
            i] > p_end]  # PASS X VALUES IF THEY ARE OUTSIDE THE PINCH RANGE
        self.depinched_y = [tup for i, tup in enumerate(current_y) if current_x[i] < p_start or current_x[
            i] > p_end]  # PASS Y VALUES IF THEY ARE OUTSIDE THE PINCH RANGE

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

        # CLOSE POPOUT
        self.EndModal(1)


class GaussianFilterDialog(wx.Dialog):
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
        self.filter_window = wx.StaticText(input_panel, -1, "Sigma")
        self.filter_window_text = wx.TextCtrl(input_panel, -1, "2")

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
        self.Bind(wx.EVT_BUTTON, self.gaussian_pass, self.b_apply_filter)

        # DEFINE SIZER'
        sizer = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        sizer.AddMany([self.obs_combo_list_text, self.obs_combo_list, self.output_name, self.output_name_text,
                       self.output_color, self.output_color_text, self.filter_window, self.filter_window_text,
                       self.b_apply_filter])
        input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def gaussian_pass(self, event):
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
                self.filtered_output[:, 1] = gaussian_filter1d(self.filter_input[:, 1], self.filter_length)
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
        self.filter_window_text = wx.TextCtrl(input_panel, -1, "1")

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

        # GET DATA TO FILTER
        for i in range(len(self.observed_list)):
            if self.observed_list[i] is not None:
                if self.observed_list[i].name == str(self.obs_to_filter_name):
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
        self.deriv = np.zeros(shape=((len(input_interpolated) - 1), 2))  # INITIALISE OUTPUT ARRAY

        for i in range(1, len(input_interpolated) - 1):
            self.deriv[i, 0] = input_interpolated[i, 0]

            # CALC USING FINITE DIFFERENCE METHOD
            self.deriv[i, 1] = abs((float(input_interpolated[i + 1, 1]) - float(input_interpolated[i - 1, 1]))) \
                               / (2. * float(self.x_inc))

        # SET FIRST AND LAST VALUES (WHICH ARE NOT CALCULATED BY FD METHOD) EQUAL TO SECOND AND SECOND-TO-LAST VALUES
        self.deriv[0, 0] = self.deriv[1, 0] - self.x_inc
        self.deriv[0, 1] = self.deriv[1, 1]
        self.deriv[-1, 0] = self.deriv[-2, 0] + self.x_inc
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

    def __init__(self, parent, id, title, layer_list, currently_active_layer_id):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        input_panel = wx.Panel(self, -1)
        self.layer_list = layer_list
        self.currently_active_layer_id = currently_active_layer_id
        self.bulk_panel_x = wx.StaticText(input_panel, -1, "X shift (km):")
        self.bulk_panel_x_Text = wx.TextCtrl(input_panel, -1, "0", size=(100, -1))
        self.bulk_panel_x_Text.SetInsertionPoint(0)
        self.bulk_panel_y = wx.StaticText(input_panel, -1, "Y shift (km):")
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
        """WHEN BUTTON IS PRESSED TO EXECUTE A BULK SHIFT"""
        # GET THE USER INPUT
        self.x_shift_value = float(self.bulk_panel_x_Text.GetValue())
        self.y_shift_value = float(self.bulk_panel_y_Text.GetValue())

        # GET THE CURRENT VALUES
        current_x = self.layer_list[self.currently_active_layer_id].x_nodes
        current_y = self.layer_list[self.currently_active_layer_id].y_nodes

        # BULK SHIFT NODES, SET TO ZERO (a.k.k 0.01) IF NEW VALUE IS < 0.01
        self.new_x = [x + self.x_shift_value if x + self.x_shift_value > 0.01 else 0.01 for x in current_x]
        self.new_y = [y + self.y_shift_value if y + self.y_shift_value > 0.01 else 0.01 for y in current_y]

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

    def __init__(self, parent, id, title, kind):
        wx.Dialog.__init__(self, parent, id, title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
                                                          | wx.MAXIMIZE_BOX | wx.MAXIMIZE_BOX)
        self.input_panel = wx.Panel(self, -1)

        self.kind = str(kind)  ## IS THE LAYER USER DEFINED OR LOADED FROM FILE?

        # CREATE FIXED LAYER BUTTON
        self.b_fixed = wx.Button(self.input_panel, -1, "Fixed layer")
        self.Bind(wx.EVT_BUTTON, self.set_fixed, self.b_fixed)
        # CREATE FLOATING LAYER BUTTON
        self.b_floating = wx.Button(self.input_panel, -1, "Floating layer")
        self.Bind(wx.EVT_BUTTON, self.set_floating, self.b_floating)
        #  DEFINE SIZER
        sizer = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        sizer.AddMany([self.b_fixed, self.b_floating])
        self.input_panel.SetSizerAndFit(sizer)
        sizer.Fit(self)

    def set_fixed(self, event):
        """ APPEND NEW LAYER BELOW LATEST FIXED LAYER"""
        self.fixed = True
        if self.kind == str('new'):
            floating_dialogbox = self.SetNewThickness(self, -1, "Set new layer thickness")
            answer = floating_dialogbox.ShowModal()
            self.new_thickness = floating_dialogbox.new_thickness
        self.EndModal(1)

    def set_floating(self, event):
        """ APPEND NEW LAYER BELOW LATEST FIXED LAYER"""
        self.fixed = False
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