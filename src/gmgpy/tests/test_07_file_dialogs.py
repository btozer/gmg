"""
Check every save dialog in gmg offers more than one file type filter.

This guards against a wxWidgets regression rather than against gmg's own logic. On macOS, wxWidgets 3.2.9
(wxPython 4.2.5) raises

    C++ assertion "nIndex < m_nCount" failed at ./include/wx/arrstr.h(227) in Item():
    wxArrayString: index out of bounds

from ``wx.FileDialog.ShowModal`` itself when a **save** dialog is given a wildcard carrying only one filter. With
a single filter macOS never builds the file type popup, so the filter index is never validly set, and that is the
index that goes out of bounds when the dialog is dismissed. The save then fails before any of gmg's own code
runs. See wxWidgets issue 26148, a duplicate of 26129, introduced by pull request 25854, and wxPython Phoenix
issue 2860. wxWidgets 3.2.8 and wxPython 4.2.4 are unaffected.

:mod:`~gmgpy.dialogs`.\\ :func:`show_save_dialog` guarantees a second filter for anything routed through it, but a
few dialogs build their own, so the source is checked directly here.

.. note:: The fault itself cannot be reproduced in a test. It needs a real dialog to be dismissed, and the native
          macOS panel runs its own modal loop that a script cannot drive. The application tests also call
          ``wx.DisableAsserts()`` and replace ``wx.FileDialog`` outright, so neither would ever see it. Checking
          the wildcards is what is left, and it is what stops the fault being reintroduced.
"""
import ast
import os
import re

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# THE MODULES THAT BUILD DIALOGS
SOURCE_FILES = ['gmg.py', 'frames.py', 'dialogs.py']


def source_of(filename):
    """Read one of the package's source files."""
    with open(os.path.join(PACKAGE_DIR, filename)) as source_file:
        return source_file.read()


def filter_count(wildcard):
    """
    How many file type filters a wildcard string offers.

    A wildcard pairs a description with a pattern, "Description (*.ext)|*.ext", and joins the pairs with pipes,
    so a single filter carries one pipe and n filters carry 2n - 1.
    """
    return (wildcard.count('|') + 1) // 2


def save_dialog_wildcards(filename):
    """
    Yield the wildcard of every wx.FileDialog in the given file that is constructed for saving.

    The call is parsed rather than pattern matched, so that a wildcard split over several source lines, or passed
    by keyword, is still found.

    :param filename: str : Name of the source file within the package

    :returns: iterator of (int, str) : Line number and wildcard of each save dialog
    """
    tree = ast.parse(source_of(filename))

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not (isinstance(node.func, ast.Attribute) and node.func.attr == 'FileDialog'):
            continue

        # THE STYLE IS THE SIXTH POSITIONAL ARGUMENT OR THE style KEYWORD
        style = None
        if len(node.args) >= 6:
            style = node.args[5]
        for keyword in node.keywords:
            if keyword.arg == 'style':
                style = keyword.value
        if style is None or 'FD_SAVE' not in ast.unparse(style):
            continue

        # THE WILDCARD IS THE FIFTH POSITIONAL ARGUMENT OR THE wildcard KEYWORD
        wildcard = None
        if len(node.args) >= 5:
            wildcard = node.args[4]
        for keyword in node.keywords:
            if keyword.arg == 'wildcard':
                wildcard = keyword.value

        # ONLY LITERAL WILDCARDS CAN BE CHECKED HERE; ONE BUILT AT RUNTIME IS COVERED BY show_save_dialog
        if wildcard is None:
            continue
        try:
            value = ast.literal_eval(wildcard)
        except ValueError:
            continue

        yield node.lineno, value


@pytest.mark.parametrize('filename', SOURCE_FILES)
def test_every_save_dialog_offers_more_than_one_filter(filename):
    """
    Test no save dialog is built with a single file type filter

    A single filter makes the dialog raise on macOS before gmg's own save code is reached.
    """
    offenders = [(line, wildcard) for line, wildcard in save_dialog_wildcards(filename)
                 if filter_count(wildcard) < 2]

    assert offenders == [], '\n'.join(
        '%s:%d has only one filter: %r' % (filename, line, wildcard) for line, wildcard in offenders)


def test_the_shared_helper_adds_a_second_filter_itself():
    """
    Test show_save_dialog guarantees a second filter, so a new caller cannot reintroduce the fault

    The helper cannot be called here, since it would open a dialog, so the guard is read from its source.
    """
    helper = re.search(r'def show_save_dialog\(.*?\n(?=\n\ndef |\Z)', source_of('dialogs.py'), re.DOTALL)
    assert helper is not None, 'show_save_dialog has been renamed or removed'

    body = helper.group(0)
    assert "wildcard.count('|') < 3" in body, 'the single filter guard has gone'
    assert 'All files' in body, 'the guard no longer appends a second filter'
    assert 'wxAssertionError' in body, 'the assertion is no longer caught'


def test_the_helper_is_used_rather_than_bare_dialogs():
    """
    Test the save routines go through the shared helper

    Anything building its own wx.FileDialog misses the assertion recovery, so the helper should be the norm.
    """
    used = sum(source_of(filename).count('show_save_dialog(') for filename in SOURCE_FILES)

    # ONE DEFINITION IN dialogs.py, PLUS THE CALL SITES
    assert used >= 7, 'expected the save routines to route through show_save_dialog, found %d uses' % used
