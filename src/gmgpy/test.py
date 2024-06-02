import wx
import wx.lib.agw.foldpanelbar as fpb

class MyFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, title="FoldPanelBar Example")
        
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)

        # Create a FoldPanelBar
        fold_panel_bar = fpb.FoldPanelBar(panel, agwStyle=fpb.FPB_VERTICAL)

        # Create a fold panel and add it to the FoldPanelBar
        fold_panel_item1 = fold_panel_bar.AddFoldPanel("Panel 1")

        # Add content to the fold panel
        label1 = wx.StaticText(fold_panel_item1, label="Content of Panel 1")
        fold_panel_item1.AddWindow(label1)

        # Create another fold panel and add it to the FoldPanelBar
        fold_panel_item2 = fold_panel_bar.AddFoldPanel("Panel 2")

        # Add content to the fold panel
        label2 = wx.StaticText(fold_panel_item2, label="Content of Panel 2")
        fold_panel_item2.AddWindow(label2)

        # Add the FoldPanelBar to the sizer
        sizer.Add(fold_panel_bar, proportion=1, flag=wx.EXPAND)

        panel.SetSizerAndFit(sizer)
        self.Fit()

app = wx.App()
frame = MyFrame()
frame.Show()
app.MainLoop()
