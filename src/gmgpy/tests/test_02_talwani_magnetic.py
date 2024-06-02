from gmgpy import talwani_and_heirtzler
from gmgpy.polygon import Polygon
import numpy as np

def get_circle_points(r, x, z, N):
    """
    :param r: float : Radius of circle (m)
    :param x: float : X position of circle center
    :param z: float : Z position of circle center
    :param N: int : Number of nodes defining the 'circle'

    :returns: 2D array (float): C1:X coordinate; C2:Z coordinate

    """
    angle = np.linspace(0., 2.*np.pi, N)
    output_points = np.zeros(shape=(len(angle), 2))

    for i in range(len(angle)):
        output_points[i, 0] = x + (r * np.cos(angle[i]))
        output_points[i, 1] = z + (r * np.sin(angle[i]))

    return output_points

def test_talwani_magnetic():
    """
    Test the bott talwani and heirtzler algorithm
    """
    # Define x and z points for calculated anomaly profile
    xp = np.linspace(0., 4000., num=4001)
    zp = np.zeros_like(xp)

    # Create test cylinder nodes (defining the test polygon) 
    input_points = get_circle_points(1000., 2000., 2000., 360)
    s = 0.001
    f = 50000.
    a = 0.
    b = 0.
    c = 0.

    test_polygon =  Polygon(input_points,{'susceptibility': s, 
                                          'angle_a': a,
                            'angle_b': b, 'angle_c': c, 'f': f})
    polygons = [test_polygon]

    # Run bott test
    anomaly = talwani_and_heirtzler.nt(xp, zp, polygons)

if __name__ == '__main__':
    test_talwani_magnetic()