#!/usr/bin/env python3

from mypymath.geometry.functions import DistanceCalculator

from geo_reader import GeoReader, RegularPrismFeomGeoPlanes


EPSILON = 0.001


def test_reader_read_polygonal_cylinders_raw():
    fname = 'test.geo' # CppPolygons produced this file.
                       # Only this convention about file format
                       # is supported now.
    reader = GeoReader(fname)
    reader.read_polygonal_cylinders_raw()
    fillers = reader.fillers
    shells = reader.shells
    filler_prisms = []
    shell_prisms = []
    for filler in fillers:
        filler_prism = RegularPrismFeomGeoPlanes(filler).prism_regular
        filler_prisms.append(filler_prism)
    #print(shells)
    for shell in shells:
        shell_prism = RegularPrismFeomGeoPlanes(shell).prism_regular
        shell_prisms.append(shell_prism)
    dc = DistanceCalculator()
    distances = []
    for shell1 in shell_prisms:
        for shell2 in shell_prisms:
            if shell1 is shell2:
                continue
            else:
                distances.append(dc.distance(shell1, shell2))
    min_distance = min(distances)
    if min_distance > EPSILON:
        print('no intersection, min_distance is', min_distance)
    else:
        print('intersection exists')


test_reader_read_polygonal_cylinders_raw()
