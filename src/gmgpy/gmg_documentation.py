import os
import webbrowser
from sys import platform
import wx

def open_documentation(event):
    """OPEN DOCUMENTATION HTML"""
    doc_url = "https://btozer.github.io/gmg/html/gmg_documentation.html"
    if platform == "linux" or platform == "linux2" or platform == "win32":
        # LINUX & WINDOWS
        webbrowser.open_new(doc_url)
    else:
        # OS X
        client = webbrowser.get("open -a /Applications/Safari.app %s")
        client.open(doc_url)

def about_gmg(event):
    """ SHOW SOFTWARE INFORMATION"""
    about = [
        "gmgpy is an Open Source Graphical User Interface (GUI) designed principally for modelling 2D potential "
        "field (gravity and magnetic) profiles. The software also includes functions for loading XY data, "
        "seismic reflection SEGY data and exploration well horizons. The software therefore provides an "
        "integrated geological/geophysical interpretation package. It is anticipated that gmgpy will be "
        "useful for undergraduate teaching purposes. \n \n"
        "Data I/O is made as simple as possible using space delimited ASCII text files. \n \n"
        "The project was instigated after failing to find an adequate open source option (in which the source "
        "code can be viewed and modified by the user) for performing 2D geophysical modeling tasks.\n \n      "
        "Inspiration came from l and The Generic Mapping Tools (GMT). \n \n"
        "gmgpy was initially developed at the University of Oxford 2014-2017. \n \n"
        "B. Tozer"]

    dlg = wx.MessageDialog(None, about[0], "About", wx.OK | wx.ICON_INFORMATION)
    result = dlg.ShowModal()
    dlg.Destroy()

def legal(event):
    """ SHOW LICENCE"""
    licence = ["Copyright 2015-2026 Brook Tozer \n\nRedistribution and use in source and binary forms, with or "
               "without modification, are permitted provided that the following conditions are met: \n \n"
               "1. Redistributions of source code must retain the above copyright notice, this list of conditions "
               "and the following disclaimer. \n\n2. Redistributions in binary form must reproduce the above "
               "copyright notice, this list of conditions and the following disclaimer in the documentation and/or "
               "other materials provided with the distribution. \n\n3. Neither the name of the copyright holder "
               "nor the names of its contributors may be used to endorse or promote products  derived from this "
               "software without specific prior written permission. \n\nTHIS SOFTWARE IS PROVIDED BY THE "
               "COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT "
               "NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE "
               "DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, "
               "INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, "
               "PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS "
               "INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,"
               " OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, "
               "EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."]

    dlg = wx.MessageDialog(None, licence[0], "BSD-3-Clause Licence", wx.OK | wx.ICON_INFORMATION)
    result = dlg.ShowModal()
    dlg.Destroy()