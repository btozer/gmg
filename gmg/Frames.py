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
    def __init__(coordinate_list, parent, id, title):
        wx.Frame.__init__(coordinate_list, None, wx.ID_ANY, 'Capture coordinates', size=(350, 500))
        coordinate_list.input_panel = wx.Panel(coordinate_list)

        # SET INSTANCE OF GMG CLASS TO RECEIVE NEW ATTRIBUTES
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
        np.savetxt(output_stream, list(zip(x, y)), delimiter=' ', fmt='%0.6f %0.6f')

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
        self.grav_frame_aspect = dcanvas_aspect

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
            print(("Selected cells " + str(attribute_edit.GetSelectedCells())))
        # IF SELECTION IS BLOCK
        if attribute_edit.attr_grid.GetSelectionBlockTopLeft():
            print(("Selection block top left " + str(attribute_edit.attr_grid.GetSelectionBlockTopLeft())))
        if attribute_edit.attr_grid.GetSelectionBlockbottomRight():
            print(("Selection block bottom right " + str(attribute_edit.attr_grid.GetSelectionBlockbottomRight())))
        # IF SELECTION IS COL
        if attribute_edit.attr_grid.GetSelectedCols():
            print(("Selected cols " + str(attribute_edit.attr_grid.GetSelectedCols())))
        # IF SELECTION IS ROW
        if attribute_edit.attr_grid.GetSelectedRows():
            print(("Selected rows " + str(attribute_edit.attr_grid.GetSelectedRows())))

    def currentcell(attribute_edit):
        # SHOW CURSOR POSITION
        row = attribute_edit.attr_grid.GetGridCursorRow()
        col = attribute_edit.attr_grid.GetGridCursorCol()
        cell = (row, col)
        print(("Current cell " + str(cell)))

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
