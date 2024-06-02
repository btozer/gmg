from gmgpy import kim_and_wessel
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

def test_kim_vgg():
    """
    Test the kim and wessel vgg algorithm
    """
    # Define x and z points for calculated anomaly profile
    xp = np.linspace(0., 4000., num=4001)
    zp = np.zeros_like(xp)

    # Create test cylinder nodes (defining the test polygon) 
    input_points = get_circle_points(1000., 2000., 2000., 360)
    density_contrast = 250
    test_polygon =  Polygon(input_points, {'density': density_contrast})
    polygons = [test_polygon]

    # Run bott test
    kim_and_wessel.gz(xp, zp, polygons)

if __name__ == '__main__':
    test_kim_vgg()