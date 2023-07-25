"""
CLASSES THAT CONTAIN GMG POP OUT FRAMES ARE PLACED HERE
"""

import wx
import matplotlib
matplotlib.use('WXAgg')
from wx.lib.agw import floatspin as fs
import wx.grid as gridlib
import numpy as np
import struct

class CaptureCoordinates(wx.Frame):
    """CAPTURE MOUSE CLICK COORDINATES AND WRITE TO DISK FILE. RETURNS ASCII TEXT FILE"""

    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, None, wx.ID_ANY, 'Capture coordinates', size=(350, 500))
        self.input_panel = wx.Panel(self)

        # SET INSTANCE OF GMG CLASS TO RECEIVE NEW ATTRIBUTES
        self.parent = parent

        # BIND PROGRAM EXIT BUTTON WITH EXIT FUNCTION
        self.Bind(wx.EVT_CLOSE, self.on_close_button)

        # CREATE LIST CONTROL
        self.table = wx.ListCtrl(self.input_panel, size=(260, 500), style=wx.LC_REPORT)
        self.table.InsertColumn(0, 'X')
        self.table.InsertColumn(1, 'Y')
        self.table.SetColumnWidth(0, 130)
        self.table.SetColumnWidth(1, 130)

        # CREATE SAVE BUTTON
        self.save_btn = wx.Button(self.input_panel, label="Save coordinates")
        self.save_btn.Bind(wx.EVT_BUTTON, self.on_save)

        # ADD FEATURES TO SIZER
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.table, 0, wx.ALL | wx.EXPAND, 5)
        sizer.Add(self.save_btn, 0, wx.ALL | wx.CENTER, 5)
        self.input_panel.SetSizer(sizer)
        sizer.Fit(self)

    def on_save(self, event):
        # CREATE OUTPUT FILE
        save_file_dialog = wx.FileDialog(self, "Save XY data", "", "", "xy files (*.xy)|*.xy",
                                         wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # USER CHANGED THEIR MIND

        # NOW WRITE OUT THE DATA
        # GET OUTPUT FILENAME
        output_stream = save_file_dialog.GetPath()

        # GET DATA
        x = []
        y = []
        for i in range(self.table.GetItemCount()):
            x.append(float(self.table.GetItem(itemIdx=i, col=0).GetText()))
            y.append(float(self.table.GetItem(itemIdx=i, col=1).GetText()))

        # OUTPUT DATA
        np.savetxt(output_stream, list(zip(x, y)), delimiter=' ', fmt='%0.6f %0.6f')

    def on_close_button(self, event):
        self.parent.capture = False
        self.parent.toolbar.ToggleTool(self.parent.t_capture_coordinates.GetId(), False)
        self.Destroy()


class PlotSettingsDialog(wx.Frame):
    """CREATE AN EXTERNAL FIGURE PLOT FROM THE MODEL. RETURNS DRAW PARAMETERS"""

    def __init__(self, parent, id, title, model_aspect, dcanvas_aspect):
        wx.Frame.__init__(self, None, wx.ID_ANY, 'Figure Construction Menu')
        input_panel = wx.Panel(self, -1)

        # CREATE INSTANCE OF MAIN FRAME CLASS TO RECEIVE NEW ATTRIBUTES
        self.parent = parent

        self.model_aspect = model_aspect
        self.grav_frame_aspect = dcanvas_aspect

        # CREATE MAIN WINDOW FRAME
        self.main_box = wx.BoxSizer(wx.HORIZONTAL)

        # MAKE SIZER
        sizer = wx.GridBagSizer(hgap=3, vgap=3)
        r = 0  # CURRENT ROW
        c = 0  # CURRENT COLUMN

        # --------------------------------------------------------------------------------------------------------------
        # 0 Line Separator
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # FIGURE OUTPUT TITLE
        r += 1
        c = 0
        self.figure_output_title = wx.StaticText(input_panel, -1, "Set figure output path:")
        sizer.Add(self.figure_output_title, pos=(r, c), span=(1, 4), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 FIGURE OUTPUT FILE NAME
        r += 1
        c = 0
        self.file_path_text = wx.TextCtrl(input_panel, -1, value="model_figure", size=(150, -1))
        sizer.Add(self.file_path_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1
        # 1 FILE TYPE SELECTION
        self.file_types = ['pdf', 'png', 'eps', 'ps']
        self.file_type_text = wx.ComboBox(input_panel, -1, value='pdf', choices=self.file_types, size=(75, -1),
                                          style=wx.CB_DROPDOWN)
        sizer.Add(self.file_type_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2 FILE BUTTON
        self.b_file_path = wx.Button(input_panel, -1, "File...")
        self.Bind(wx.EVT_BUTTON, self.file_path, self.b_file_path)
        sizer.Add(self.b_file_path, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 3 USE TIGHT LAYOUT SWITCH
        self.use_tight_layout_checkbox = wx.CheckBox(input_panel, -1, "Use tight\nlayout?")
        sizer.Add(self.use_tight_layout_checkbox, pos=(r, c), span=(1, 1), flag=wx.ALIGN_BOTTOM, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0-3 STATIC LINE SEPARATOR
        r += 1
        c = 0
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 FIGURE FONT SIZE TEXT LABEL
        r += 1
        c = 0
        self.set_fs = wx.StaticText(input_panel, -1, "Text size:")
        sizer.Add(self.set_fs, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1 FIGURE FONT SIZE
        self.fs_text = fs.FloatSpin(input_panel, -1, min_val=0.01, max_val=20.0, increment=0.1, value=8., size=(75, -1))
        self.fs_text.SetFormat("%f")
        self.fs_text.SetDigits(3)
        sizer.Add(self.fs_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)

        # 2 MODEL ASPECT RATIO
        c += 1
        self.set_aspect_ratio = wx.StaticText(input_panel, -1, "Model aspect\nratio:")
        sizer.Add(self.set_aspect_ratio, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 3 MODEL ASPECT RATIO TEXT
        self.aspect_ratio_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=2000.0, increment=0.1,
                                              value=self.model_aspect, size=(75, -1))
        self.aspect_ratio_text.SetFormat("%f")
        self.aspect_ratio_text.SetDigits(2)
        sizer.Add(self.aspect_ratio_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        r += 1
        c = 0

        # 0 POINT SIZE
        self.set_ps = wx.StaticText(input_panel, -1, "Observed point size:")
        sizer.Add(self.set_ps, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1 MARKER SIZE TEXT
        self.ps_text = fs.FloatSpin(input_panel, -1, min_val=0.01, max_val=20.0, increment=0.1, value=0.5,
                                    size=(75, -1))
        self.ps_text.SetFormat("%f")
        self.ps_text.SetDigits(2)
        sizer.Add(self.ps_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2 LINE WIDTH

        self.set_lw = wx.StaticText(input_panel, -1, "Calculated\nline width:")
        sizer.Add(self.set_lw, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 3 LINE WIDTH TEXT
        self.lw_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=1.0,
                                    size=(75, -1))
        self.lw_text.SetFormat("%f")
        self.lw_text.SetDigits(2)
        sizer.Add(self.lw_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 FONT TYPE
        r += 1
        c = 0
        self.fonts = ['Times New Roman', 'Times', 'Courier', 'Courier New', 'Helvetica', 'Sans', 'verdana', 'Arial']
        self.set_font_type = wx.StaticText(input_panel, -1, "Font type:")
        sizer.Add(self.set_font_type, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1 FONT TYPE TEXT
        self.font_type_text = wx.ComboBox(input_panel, -1, value='Times New Roman', choices=self.fonts, size=(75, -1),
                                          style=wx.CB_DROPDOWN)
        sizer.Add(self.font_type_text, pos=(r, c), span=(1, 2), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LINE SEPARATOR
        r += 1
        c = 0
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 OBSERVED TOPO AXIS Y-MIN & Y-MAX
        r += 1
        c = 0
        self.topo_frame_text = wx.StaticText(input_panel, -1, "Observed topography\ny-axis max/min:")
        sizer.Add(self.topo_frame_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1
        self.topo_min_text = fs.FloatSpin(input_panel, -1, min_val=-2000., max_val=2000.0, increment=1, value=-100,
                                          size=(75, -1))
        self.topo_min_text.SetFormat("%f")
        self.topo_min_text.SetDigits(2)
        sizer.Add(self.topo_min_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2
        self.topo_max_text = fs.FloatSpin(input_panel, -1, min_val=-2000., max_val=2000.0, increment=1, value=100,
                                          size=(75, -1))
        self.topo_max_text.SetFormat("%f")
        self.topo_max_text.SetDigits(2)
        sizer.Add(self.topo_max_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 OBSERVED GRAVITY AXIS Y-MIN & Y-MAX
        r += 1
        c = 0
        self.grav_frame_text = wx.StaticText(input_panel, -1, "Observed gravity\ny-axis max/min:")
        sizer.Add(self.grav_frame_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1
        self.grav_min_text = fs.FloatSpin(input_panel, -1, min_val=-2000., max_val=2000.0, increment=1, value=-100,
                                          size=(75, -1))
        self.grav_min_text.SetFormat("%f")
        self.grav_min_text.SetDigits(2)
        sizer.Add(self.grav_min_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2
        self.grav_max_text = fs.FloatSpin(input_panel, -1, min_val=-2000., max_val=2000.0, increment=1, value=100,
                                          size=(75, -1))
        self.grav_max_text.SetFormat("%f")
        self.grav_max_text.SetDigits(2)
        sizer.Add(self.grav_max_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 OBSERVED GRAVITY AXIS Y-MIN & Y-MAX
        r += 1
        c = 0
        self.mag_frame_text = wx.StaticText(input_panel, -1, "Observed magnetic\ny-axis max/min:")
        sizer.Add(self.mag_frame_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1
        self.mag_min_text = fs.FloatSpin(input_panel, -1, min_val=-2000., max_val=2000.0, increment=1, value=-100,
                                         size=(75, -1))
        self.mag_min_text.SetFormat("%f")
        self.mag_min_text.SetDigits(2)
        sizer.Add(self.mag_min_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2
        self.mag_max_text = fs.FloatSpin(input_panel, -1, min_val=-2000., max_val=2000.0, increment=1, value=100,
                                         size=(75, -1))
        self.mag_max_text.SetFormat("%f")
        self.mag_max_text.SetDigits(2)
        sizer.Add(self.mag_max_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LINE SEPARATOR
        r += 1
        c = 0
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LAYER POLYGONS
        r += 1
        c = 0
        self.draw_polygons_checkbox = wx.CheckBox(input_panel, -1, " Draw layer polygons?")
        self.draw_polygons_checkbox.SetValue(True)
        sizer.Add(self.draw_polygons_checkbox, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1. POLYGON ALPHA
        self.set_poly_alpha = wx.StaticText(input_panel, -1, "Polygon\ntransparency:")
        sizer.Add(self.set_poly_alpha, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2 POLYGON ALPHA TEXT
        self.poly_alpha_text = fs.FloatSpin(input_panel, -1, min_val=0.01, max_val=1.0, increment=0.01, value=0.5,
                                            size=(75, -1))
        self.poly_alpha_text.SetFormat("%f")
        self.poly_alpha_text.SetDigits(2)
        sizer.Add(self.poly_alpha_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)

        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LAYER LINES
        r += 1
        c = 0
        self.draw_fixed_layers_checkbox = wx.CheckBox(input_panel, -1, " Draw fixed layers?")
        self.draw_fixed_layers_checkbox.SetValue(True)
        sizer.Add(self.draw_fixed_layers_checkbox, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1 LAYER LINE WIDTH
        self.set_layer_lw = wx.StaticText(input_panel, -1, "Layer line\nwidth:")
        sizer.Add(self.set_layer_lw, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2 LAYER LINE WIDTH TEXT
        self.layer_lw_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=0.1,
                                          size=(75, -1))
        self.layer_lw_text.SetFormat("%f")
        self.layer_lw_text.SetDigits(3)
        sizer.Add(self.layer_lw_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 FLOATING LAYER LINES
        r += 1
        c = 0
        self.draw_floating_layer_lines_checkbox = wx.CheckBox(input_panel, -1, " Draw floating layers?")
        self.draw_floating_layer_lines_checkbox.SetValue(True)
        sizer.Add(self.draw_floating_layer_lines_checkbox, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT,
                  border=10)
        c += 1
        # 1 LAYER LINE TRANSPARENCY
        self.set_layer_line_alpha = wx.StaticText(input_panel, -1, "Layer line\ntransparency:")
        sizer.Add(self.set_layer_line_alpha, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2 LAYER LINE TRANSPARENCY TEXT
        self.layer_line_alpha_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.7,
                                                  size=(75, -1))
        self.layer_line_alpha_text.SetFormat("%f")
        self.layer_line_alpha_text.SetDigits(2)
        sizer.Add(self.layer_line_alpha_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LINE SEPARATOR
        r += 1
        c = 0
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 13.0 DRAW COLORBAR LINES
        r += 1
        c = 0
        self.draw_colorbar_checkbox = wx.CheckBox(input_panel, -1, "Draw colorbar?")
        self.draw_colorbar_checkbox.SetValue(True)
        sizer.Add(self.draw_colorbar_checkbox, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 14.0 COLORBAR XY
        r += 1
        c = 0
        self.colorbar_xy = wx.StaticText(input_panel, -1, "Colorbar XY Location:")
        sizer.Add(self.colorbar_xy, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 14.1 COLORBAR XY
        self.colorbar_x_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.7,
                                            size=(75, -1))
        self.colorbar_x_text.SetFormat("%f")
        self.colorbar_x_text.SetDigits(2)
        sizer.Add(self.colorbar_x_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 14.2 COLORBAR XY
        self.colorbar_y_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.55,
                                            size=(75, -1))
        self.colorbar_y_text.SetFormat("%f")
        self.colorbar_y_text.SetDigits(2)
        sizer.Add(self.colorbar_y_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 15.0 COLORBAR SIZE
        r += 1
        c = 0
        self.colorbar_size = wx.StaticText(input_panel, -1, "Colorbar X/Y text size:")
        sizer.Add(self.colorbar_size, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 15.1
        self.colorbar_size_x_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.15,
                                                 size=(75, -1))
        self.colorbar_size_x_text.SetFormat("%f")
        self.colorbar_size_x_text.SetDigits(2)
        sizer.Add(self.colorbar_size_x_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 15.2
        self.colorbar_size_y_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=1.0, increment=0.1, value=0.005,
                                                 size=(75, -1))
        self.colorbar_size_y_text.SetFormat("%f")
        self.colorbar_size_y_text.SetDigits(3)
        sizer.Add(self.colorbar_size_y_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LINE SEPARATOR
        r += 1
        c = 0
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 DRAW XY DATA
        r += 1
        c = 0
        self.draw_xy_checkbox = wx.CheckBox(input_panel, -1, "Draw XY data?")
        self.draw_xy_checkbox.SetValue(True)
        sizer.Add(self.draw_xy_checkbox, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)

        # --------------------------------------------------------------------------------------------------------------
        # 0 XY POINT SIZE
        r += 1
        c = 0
        self.set_xy_size = wx.StaticText(input_panel, -1, "XY point size:")
        sizer.Add(self.set_xy_size, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1 XY POINT SIZE TEXT
        self.xy_size_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=1.5,
                                         size=(75, -1))
        self.xy_size_text.SetFormat("%f")
        self.xy_size_text.SetDigits(2)
        sizer.Add(self.xy_size_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1
        # 2 XY COLOR
        self.colors = ['red', 'orange', 'yellow', 'green', 'blue', 'grey', 'white', 'black']
        self.set_xy_color = wx.StaticText(input_panel, -1, "XY color:")
        sizer.Add(self.set_xy_color, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 3 XY COLOR TEXT
        self.xy_color_text = wx.ComboBox(input_panel, -1, value='black', choices=self.colors, size=(75, -1),
                                         style=wx.CB_DROPDOWN)
        sizer.Add(self.xy_color_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LINE SEPARATOR
        r += 1
        c = 0
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 DRAW WELLS
        r += 1
        c = 0
        self.draw_wells_checkbox = wx.CheckBox(input_panel, -1, "Draw wells?")
        self.draw_wells_checkbox.SetValue(True)
        sizer.Add(self.draw_wells_checkbox, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 WELL FONT SIZE
        r += 1
        c = 0
        self.set_well_font_size = wx.StaticText(input_panel, -1, "Well font size:")
        sizer.Add(self.set_well_font_size, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1 WELL FONT SIZE TEXT
        self.well_font_size_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=1.5,
                                                size=(75, -1))
        self.well_font_size_text.SetFormat("%f")
        self.well_font_size_text.SetDigits(2)
        sizer.Add(self.well_font_size_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 2 WELL LINE WIDTH

        self.set_well_lw = wx.StaticText(input_panel, -1, "Well line width:")
        sizer.Add(self.set_well_lw, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 3 WELL LINE WIDTH TEXT
        self.well_lw_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=0.1,
                                         size=(75, -1))
        self.well_lw_text.SetFormat("%f")
        self.well_lw_text.SetDigits(3)
        sizer.Add(self.well_lw_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LINE SEPARATOR
        r += 1
        c = 0
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 DRAW FAULTS
        r += 1
        c = 0
        self.draw_faults_checkbox = wx.CheckBox(input_panel, -1, "Draw Faults?")
        self.draw_faults_checkbox.SetValue(True)
        sizer.Add(self.draw_faults_checkbox, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 FAULTS LINE WIDTH
        r += 1
        c = 0
        self.set_faults_lw = wx.StaticText(input_panel, -1, "Faults line width:")
        sizer.Add(self.set_faults_lw, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1 FAULTS LINE WIDTH TEXT
        self.faults_lw_text = fs.FloatSpin(input_panel, -1, min_val=0.0, max_val=20.0, increment=0.1, value=0.1,
                                           size=(75, -1))
        self.faults_lw_text.SetFormat("%f")
        self.faults_lw_text.SetDigits(3)
        sizer.Add(self.faults_lw_text, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 LINE SEPARATOR
        r += 1
        c = 0
        line = wx.StaticLine(input_panel)
        sizer.Add(line, pos=(r, c), span=(1, 4), flag=wx.EXPAND | wx.ALIGN_LEFT, border=10)
        c += 1
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # 0 DRAW BUTTON
        r += 1
        c = 0
        self.b_draw_button = wx.Button(input_panel, -1, "Draw figure")
        self.Bind(wx.EVT_BUTTON, self.draw_button, self.b_draw_button)
        sizer.Add(self.b_draw_button, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        c += 1

        # 1 EXIT BUTTON
        self.b_exit = wx.Button(input_panel, -1, "Exit")
        self.Bind(wx.EVT_BUTTON, self.exit, self.b_exit)
        sizer.Add(self.b_exit, pos=(r, c), span=(1, 1), flag=wx.ALIGN_LEFT, border=10)
        # --------------------------------------------------------------------------------------------------------------

        # NOW ADD THE SIZER TO THE PANEL
        self.main_box.Add(sizer, proportion=1, flag=wx.ALL | wx.EXPAND, border=10)
        input_panel.SetSizerAndFit(self.main_box)
        self.main_box.Fit(self)

    def draw_button(self, event):
        self.file_path = str(self.file_path_text.GetValue())
        self.file_type = str(self.file_type_text.GetValue())
        self.use_tight_layout = self.use_tight_layout_checkbox.GetValue()
        #
        self.fs = float(self.fs_text.GetValue())  # FONT SIZE
        self.aspect_ratio = float(self.aspect_ratio_text.GetValue())  # MODEL ASPECT RATIO
        self.ps = float(self.ps_text.GetValue())  # OBSERVED POINT SIZE?
        self.lw = float(self.lw_text.GetValue())  # CALCUALTED LINE WIDTH
        #
        self.font_type = str(self.font_type_text.GetValue())
        #
        self.topo_frame_min = float(self.topo_min_text.GetValue())
        self.topo_frame_max = float(self.topo_max_text.GetValue())
        self.grav_frame_min = float(self.grav_min_text.GetValue())
        self.grav_frame_max = float(self.grav_max_text.GetValue())
        self.mag_frame_min = float(self.mag_min_text.GetValue())
        self.mag_frame_max = float(self.mag_max_text.GetValue())
        #
        self.draw_polygons = self.draw_polygons_checkbox.GetValue()
        self.polygon_alpha = float(self.poly_alpha_text.GetValue())
        #
        self.draw_fixed_layers = self.draw_fixed_layers_checkbox.GetValue()
        self.layer_line_width = float(self.layer_lw_text.GetValue())
        #
        self.draw_floating_layers = self.draw_floating_layer_lines_checkbox.GetValue()
        self.layer_line_alpha = float(self.layer_line_alpha_text.GetValue())
        #
        self.draw_colorbar = self.draw_colorbar_checkbox.GetValue()
        self.colorbar_x = float(self.colorbar_x_text.GetValue())
        self.colorbar_y = float(self.colorbar_y_text.GetValue())
        self.colorbar_size_x = float(self.colorbar_size_x_text.GetValue())
        self.colorbar_size_y = float(self.colorbar_size_y_text.GetValue())
        #
        self.draw_xy_data = self.draw_xy_checkbox.GetValue()
        self.xy_size = self.xy_size_text.GetValue()
        self.xy_color = str(self.xy_color_text.GetValue())
        #
        self.draw_wells = self.draw_wells_checkbox.GetValue()
        self.well_fs = float(self.well_font_size_text.GetValue())
        self.well_line_width = float(self.well_lw_text.GetValue())
        #
        self.draw_faults = self.draw_faults_checkbox.GetValue()
        self.faults_lw = float(self.faults_lw_text.GetValue())

        # CALL DRAW FUNCTION
        self.parent.draw_model()

    def exit(self, event):
        self.Destroy()

    def file_path(self, event):
        self.save_file_dialog = wx.FileDialog(self, "Save As", "", "", "Figure files (*.*)|*.*",
                                              wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if self.save_file_dialog.ShowModal() == wx.ID_CANCEL:
            return  # THE USER CHANGED IDEA...
        self.chosen_path = self.save_file_dialog.GetPath()
        self.file_path_text.SetValue(str(self.chosen_path))
        self.save_file_dialog.Destroy()


class AttributeEditor(wx.Frame):
    """OPENS A TABLE FOR VIEWING AND EDITING LABEL ATTRIBUTES"""

    def __init__(self, parent, id, title, tree_items, layer_list):
        wx.Frame.__init__(self, None, wx.ID_ANY, 'Attribute editor', size=(1000, 500))
        self.input_panel = wx.Panel(self)

        # CREATE INSTANCE OF MAIN FRAME CLASS TO RECEIVE NEW ATTRIBUTES
        self.parent = parent

        # SET VARIABLES
        self.tree_items = tree_items
        self.layer_list = layer_list

        # DEFINE ATTRIBUTE GRID
        self.attr_grid = gridlib.Grid(self.input_panel, -1, size=(1000, 500))
        self.attr_grid.CreateGrid(len(self.tree_items) - 1, 9)
        self.attr_grid.SetColLabelValue(0, 'Layer Name')
        self.attr_grid.SetColLabelValue(1, 'Density')
        self.attr_grid.SetColLabelValue(2, 'Reference density')
        self.attr_grid.SetColLabelValue(3, 'Susceptibility')
        self.attr_grid.SetColLabelValue(4, 'Angle A')
        self.attr_grid.SetColLabelValue(5, 'Angle B')
        self.attr_grid.SetColLabelValue(6, 'Angle C')
        self.attr_grid.SetColLabelValue(7, 'Earth Field')
        self.attr_grid.SetColLabelValue(8, 'Layer color')

        # SET COLUMN FORMATS
        self.attr_grid.SetColFormatFloat(1, 3, 2)
        self.attr_grid.SetColFormatFloat(2, 3, 2)
        self.attr_grid.SetColFormatFloat(3, 9, 7)
        self.attr_grid.SetColFormatFloat(4, 3, 2)
        self.attr_grid.SetColFormatFloat(5, 3, 2)
        self.attr_grid.SetColFormatFloat(6, 3, 2)
        self.attr_grid.SetColFormatFloat(7, 3, 2)

        # CREATE SET BUTTON
        self.b_set_attr_button = wx.Button(self.input_panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.set_attr_button, self.b_set_attr_button)

        # POPULATE ATTRIBUTE TABLE
        self.length = len(self.tree_items)

        for i in range(self.length - 1):
            self.attr_grid.SetCellValue(i, 0, self.tree_items[i + 1])
            self.attr_grid.SetCellValue(i, 1, str(float(self.layer_list[i + 1].density) / 1000.))
            self.attr_grid.SetCellValue(i, 2, str(float(self.layer_list[i + 1].reference_density) / 1000.))
            self.attr_grid.SetCellValue(i, 3, str(self.layer_list[i + 1].susceptibility))
            self.attr_grid.SetCellValue(i, 4, str(self.layer_list[i + 1].angle_a))
            self.attr_grid.SetCellValue(i, 5, str(self.layer_list[i + 1].angle_b))
            self.attr_grid.SetCellValue(i, 6, str(self.layer_list[i + 1].angle_c))
            self.attr_grid.SetCellValue(i, 7, str(self.layer_list[i + 1].earth_field))
            self.attr_grid.SetCellValue(i, 8, str(self.layer_list[i + 1].color))

        # SET SIZER
        for col in range(8):
            self.attr_grid.SetColSize(col, 12)

        self.attr_grid.AutoSize()
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.b_set_attr_button)
        sizer.Add(self.attr_grid)
        self.input_panel.SetSizer(sizer)
        sizer.Fit(self)

        # ACTION BINDINGS
        self.Bind(wx.EVT_CHAR_HOOK, self.on_key)
        # self.attr_grid.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.open_colour_box)
        self.attr_grid.Bind(wx.grid.EVT_GRID_CELL_LEFT_DCLICK, self.open_colour_box)
        # self.attr_grid.Bind(wx.EVT_SIZE, self.on_size)

    def open_colour_box(self, event):
        """IF THE ROW SELECTED IS 8 (THE LAYER COLOR) THEN OPEN THE COLOR SELECTION WIDGET"""
        if event.GetCol() == 8:
            row = event.GetRow()
            self.on_color_dlg(event, row)
        else:
            pass

    def on_color_dlg(self, event, row):
        """SET COLOUR FOR LAYER"""
        dlg = wx.ColourDialog(self)

        # ENSURE THE FULL COLOUR DIALOG IS DISPLAYED, NOT THE ABBREVIATED VERSION
        dlg.GetColourData().SetChooseFull(True)

        if dlg.ShowModal() == wx.ID_OK:
            rgb = dlg.GetColourData().GetColour().Get()
            rgb = rgb[0:3]
            html = struct.pack('BBB', *rgb).encode('hex')
            self.attr_grid.SetCellValue(row, 8, '#' + str(html))
        dlg.Destroy()

    def on_key(self, event):
        if event.ControlDown() and event.GetKeyCode() == 67:
            self.selection()
            self.copy()  # CALL COPY METHOD
        # SKIP OTHER KEY EVENTS
        if event.GetKeyCode():
            event.Skip()
            return

    def selection(self):
        # SHOW CELL SELECTION
        # IF SELECTION IS CELL
        if self.attr_grid.GetSelectedCells():
            print(("Selected cells " + str(self.GetSelectedCells())))
        # IF SELECTION IS BLOCK
        if self.attr_grid.GetSelectionBlockTopLeft():
            print(("Selection block top left " + str(self.attr_grid.GetSelectionBlockTopLeft())))
        if self.attr_grid.GetSelectionBlockbottomRight():
            print(("Selection block bottom right " + str(self.attr_grid.GetSelectionBlockbottomRight())))
        # IF SELECTION IS COL
        if self.attr_grid.GetSelectedCols():
            print(("Selected cols " + str(self.attr_grid.GetSelectedCols())))
        # IF SELECTION IS ROW
        if self.attr_grid.GetSelectedRows():
            print(("Selected rows " + str(self.attr_grid.GetSelectedRows())))

    def currentcell(self):
        # SHOW CURSOR POSITION
        row = self.attr_grid.GetGridCursorRow()
        col = self.attr_grid.GetGridCursorCol()
        cell = (row, col)
        print(("Current cell " + str(cell)))

    def copy(self):
        # NUMBER OF ROWS AND COLS
        rows = self.attr_grid.GetSelectionBlockbottomRight()[0][0] - \
               self.attr_grid.GetSelectionBlockTopLeft()[0][0] + 1
        cols = self.attr_grid.GetSelectionBlockbottomRight()[0][1] - \
               self.attr_grid.GetSelectionBlockTopLeft()[0][1] + 1

        # DATA VARIABLE CONTAIN TEXT THAT MUST BE SET IN THE CLIPBOARD
        data = ''

        # FOR EACH CELL IN SELECTED range APPEND THE CELL VALUE IN THE DATA VARIABLE
        # TABS '\t' FOR COLS AND '\r' FOR ROWS
        for r in range(rows):
            for c in range(cols):
                data = data + str(
                    self.attr_grid.GetCellValue(self.attr_grid.GetSelectionBlockTopLeft()[0][0] + r,
                                                self.attr_grid.GetSelectionBlockTopLeft()[0][1] + c))
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

    def set_attr_button(self, event):
        """WHEN THE SET ATTRIBUTES BUTTON IS PRESSED"""
        # RECREATE ARRAYS (INCLUDE VALUES FOR "LAYER 0)"
        self.tree_items = ['Layer 1']

        # SET NEW TREE NAMES AND LAYER ATTRIBUTES FOR LAYERS 1 to i
        for i in range(self.length - 1):
            self.tree_items.append(str(self.attr_grid.GetCellValue(i, 0)))
            self.layer_list[i + 1].density = float(self.attr_grid.GetCellValue(i, 1)) * 1000.
            self.layer_list[i + 1].reference_density = float(self.attr_grid.GetCellValue(i, 2)) * 1000.
            self.layer_list[i + 1].susceptibility = float(self.attr_grid.GetCellValue(i, 3))
            self.layer_list[i + 1].angle_a = float(self.attr_grid.GetCellValue(i, 4))
            self.layer_list[i + 1].angle_b = float(self.attr_grid.GetCellValue(i, 5))
            self.layer_list[i + 1].angle_c = float(self.attr_grid.GetCellValue(i, 6))
            self.layer_list[i + 1].earth_field = float(self.attr_grid.GetCellValue(i, 7))
            self.layer_list[i + 1].color = str(self.attr_grid.GetCellValue(i, 8))

        # UPDATE MAIN FRAME
        self.parent.attribute_set(self.tree_items, self.layer_list)

        # UPDATE GMG MODEL
        self.parent.update_layer_data()
        self.parent.run_algorithms()
        self.parent.draw()
