from mypymath.geometry.figures import Point, Plane, PolygonRegular, PrismRegular
from mypymath.linalg import Vector, Matrix3


class RegularPrismFeomGeoPlanes:
    """
    Produces a convertion from
      a regular prism made from planes and stored in .geo
    to
      PrismRegular that may be analyzed with mypymath module.
    """
    def __init__(self, planes):
        top = planes[0]
        bottom = planes[1]
        sides = planes[2:]
        top_poly_vertices = []
        bottom_poly_vertices = []
        for i in range(len(sides)):
            # find intersection of side[i], side[i-1] and top and
            #      intersection of side[i],  side[i-1] and bottom
            side_one = sides[i]
            side_two = sides[i - 1]
            # with top:
            det_system_matrix = Matrix3([side_one.a, side_one.b, side_one.c,
                                         side_two.a, side_two.b, side_two.c,
                                         top.a, top.b, top.c]).det()
            det_pt_x_matrix = Matrix3([-side_one.d, side_one.b, side_one.c,
                                       -side_two.d, side_two.b, side_two.c,
                                       -top.d, top.b, top.c]).det()
            det_pt_y_matrix = Matrix3([side_one.a, -side_one.d, side_one.c,
                                       side_two.a, -side_two.d, side_two.c,
                                       top.a, -top.d, top.c]).det()
            det_pt_z_matrix = Matrix3([side_one.a, side_one.b, -side_one.d,
                                       side_two.a, side_two.b, -side_two.d,
                                       top.a, top.b, -top.d]).det()
            if det_system_matrix == 0:
                print('Error in RegularPrismFromGeoPlanes:',
                      'det is null')
                return None
            pt_top = Point(det_pt_x_matrix / det_system_matrix,
                           det_pt_y_matrix / det_system_matrix,
                           det_pt_z_matrix / det_system_matrix)
            # with bottom:
            det_system_matrix = Matrix3([side_one.a, side_one.b, side_one.c,
                                         side_two.a, side_two.b, side_two.c,
                                         bottom.a, bottom.b, bottom.c]).det()
            det_pt_x_matrix = Matrix3([-side_one.d, side_one.b, side_one.c,
                                       -side_two.d, side_two.b, side_two.c,
                                       -bottom.d, bottom.b, bottom.c]).det()
            det_pt_y_matrix = Matrix3([side_one.a, -side_one.d, side_one.c,
                                       side_two.a, -side_two.d, side_two.c,
                                       bottom.a, -bottom.d, bottom.c]).det()
            det_pt_z_matrix = Matrix3([side_one.a, side_one.b, -side_one.d,
                                       side_two.a, side_two.b, -side_two.d,
                                       bottom.a, bottom.b, -bottom.d]).det()
            if det_system_matrix == 0:
                print('Error in RegularPrismFromGeoPlanes:',
                      'det is null')
                return None
            pt_bottom = Point(det_pt_x_matrix / det_system_matrix,
                              det_pt_y_matrix / det_system_matrix,
                              det_pt_z_matrix / det_system_matrix)
            top_poly_vertices.append(pt_top)
            bottom_poly_vertices.append(pt_bottom)
        top_facet = PolygonRegular(top_poly_vertices)
        bottom_facet = PolygonRegular(bottom_poly_vertices)
        self.prism_regular = PrismRegular(top_facet, bottom_facet)

class GeoReader:
    def __init__(self, fname):
        self.fname = fname
        self.cell_xlen = None
        self.cell_ylen = None
        self.cell_zlen = None

    def read_polygonal_cylinders_raw(self):
        """
        Reads geofile as a file containing PolygonalCylinders
        in a cubic cell. A lot of hardcode that usues the format
        of CppPolygons GeoWriter for easy parsing.
        Returns a list of PolygonalCylinders in a form of
        a list of PrismRegular.

        """
        pcs = [] # Here PolygonalCylinders will be stored
                 # pc is a very good abbreviation for me
                 # although it does not look this way.
        with open(self.fname, 'r') as f:
            self.fillers = []
            self.shells = []
            filler = None
            shell = None
            current_strusture = None
            for line in f:
                 # heading; does not matter
                 if line.startswith('algebraic3d'):
                     continue
                 # cubic cell definition
                 if line.startswith('solid cell'):
                     orthobrick_str = line.split(' = ')[1]
                     pt_coords_str = orthobrick_str.split('(')[1]
                     pt_coords_str = pt_coords_str.split(')')[0]
                     pt_coords_begin = (pt_coords_str.split(';')[0]).split()
                     pt_coords_end = (pt_coords_str.split(';')[1]).split()
                     if len(pt_coords_begin) != 3 or len(pt_coords_end) != 3:
                         print('Error in reader.read_polygonal_cylinders_raw:',
                               'something wrong with cell dimensions')
                         return None
                     self.cell_xlen = (float(pt_coords_end[0].split(',')[0]) -
                                       float(pt_coords_begin[0].split(',')[0]))
                     self.cell_ylen = (float(pt_coords_end[1].split(',')[0]) -
                                       float(pt_coords_begin[2].split(',')[0]))
                     self.cell_zlen = (float(pt_coords_end[2]) -
                                       float(pt_coords_begin[2].split(';')[0]))
                 # Filler particle definition
                 if line.startswith('solid polygonalCylinder'):
                     filler = []
                 # Interface layer definition
                 if line.startswith('solid pc'):
                     shell = []
                 # Filler of shell definition started here
                 if (line.startswith(' plane') or
                     line.startswith(' and plane') or # well, this is not good...
                     line.startswith('and plane')):
                     brackets = line.split('plane')[1]
                     plane = self.plane_from_brackets(brackets)
                     if filler is not None:
                         # Reading filler now
                         filler.append(plane)
                     elif shell is not None:
                         # Reading shell now
                         shell.append(plane)
                     else:
                         # Particle but not shell and not filler?!
                         print('Error in GeoReader.read_raw:',
                               'tag = asdfadsgqw')
                         return None
                 if line.endswith('and cell;\n'):
                     # It is the last plane from planes defining a particle.
                     if filler is not None:
                         self.fillers.append(filler)
                         filler = None
                     elif shell is not None:
                         self.shells.append(shell)
                         shell = None

    def plane_from_brackets(self, brackets):
        """
        Makes a plane from a line called 'brackets':
        brackets == (pt.x, pt.y, pt.z; normal.x, normal.y, normal.z)
        Returns a Plane that may be used in mypymath.
        """
        pt_coords = brackets.split('(')[1].split(';')[0].split(', ')
        ptx = float(pt_coords[0])
        pty = float(pt_coords[1])
        ptz = float(pt_coords[2])
        pt = Point(ptx, pty, ptz)
        normal_coords = brackets.split('(')[1].split(';')[1].split(')')[0]
        normal_coords = normal_coords.split(', ')
        nx = float(normal_coords[0])
        ny = float(normal_coords[1])
        nz = float(normal_coords[2])
        normal = Vector(nx, ny, nz)
        return Plane(pt, normal)
