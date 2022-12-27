from cover import (
    cover_shape,
    _cover_point,
    _cover_multipoint,
    _cover_linestring,
    _cover_multilinestring,
    _cover_polygon,
    _cover_multipolygon,
    cover_shape_exact,
    _cover_polygon_exact,
    _cover_multipolygon_exact,
    index_to_point,
    index_to_polygon,
    antimeridian_handler,
)

from shapely.geometry import (
    Point,
    MultiPoint,
    LineString,
    MultiLineString,
    Polygon,
    MultiPolygon,
)

import pandas as pd


antimeridian_index = "82f387fffffffff"  # known to cross the antimeridian

point = Point([(16, 36)])

mpoint = MultiPoint([(16, 36), (18, 38)])

line = LineString([(16, 36), (18, 38)])

mline = MultiLineString([[(16, 36), (18, 38)], [(16, 38), (18, 36)]])

polygon = Polygon([(16, 36), (16, 38), (18, 38), (18, 36)])

mpolygon = MultiPolygon(
    [
        Polygon(
            shell=((16, 36), (16, 38), (18, 38), (18, 36)),
            holes=[((16.5, 36.5), (16.5, 37.5), (17.5, 37.5), (17.5, 36.5))],
        ),
        Polygon([(19, 37), (19, 38), (20, 38), (20, 37)]),
    ]
)


def test_cover_shape():
    assert cover_shape(point, 4) == ["843f315ffffffff"]
    assert cover_shape(mpoint, 4) == ["843f203ffffffff", "843f315ffffffff"]
    assert cover_shape(line, 4) == [
        "843f227ffffffff",
        "843f311ffffffff",
        "843f207ffffffff",
        "843f22bffffffff",
        "843f315ffffffff",
        "843f221ffffffff",
        "843f203ffffffff",
        "843f31bffffffff",
    ]
    assert cover_shape(mline, 4) == [
        "843f227ffffffff",
        "843f311ffffffff",
        "843f207ffffffff",
        "843f22bffffffff",
        "843f041ffffffff",
        "843f22dffffffff",
        "843f049ffffffff",
        "843f315ffffffff",
        "843f043ffffffff",
        "843f221ffffffff",
        "843f203ffffffff",
        "843f31bffffffff",
        "843f267ffffffff",
    ]
    assert cover_shape(polygon, 4) == [
        "843f045ffffffff",
        "843f311ffffffff",
        "843f207ffffffff",
        "843f22bffffffff",
        "843f235ffffffff",
        "843f22dffffffff",
        "843f049ffffffff",
        "843f313ffffffff",
        "843f229ffffffff",
        "843f04dffffffff",
        "843f23dffffffff",
        "843f225ffffffff",
        "843f31bffffffff",
        "843f267ffffffff",
        "843f227ffffffff",
        "843f041ffffffff",
        "843f223ffffffff",
        "843f319ffffffff",
        "843f221ffffffff",
    ]
    assert cover_shape(mpolygon, 4) == [
        "843f045ffffffff",
        "843f311ffffffff",
        "843f207ffffffff",
        "843f22bffffffff",
        "843f217ffffffff",
        "843f235ffffffff",
        "843f22dffffffff",
        "843f049ffffffff",
        "843f2e9ffffffff",
        "843f313ffffffff",
        "843f229ffffffff",
        "843f23dffffffff",
        "843f04dffffffff",
        "843f225ffffffff",
        "843f2e5ffffffff",
        "843f267ffffffff",
        "843f2edffffffff",
        "843f041ffffffff",
        "843f2e1ffffffff",
        "843f319ffffffff",
    ]


def test_cover_shape_exact():
    assert cover_shape_exact(
        polygon, resolution=4, precision_resolution=5, inside_only=False
    ) == [
        "843f235ffffffff",
        "843f049ffffffff",
        "843f229ffffffff",
        "843f04dffffffff",
        "843f35bffffffff",
        "843f265ffffffff",
        "843f267ffffffff",
        "843f227ffffffff",
        "843f041ffffffff",
        "843f205ffffffff",
        "843f045ffffffff",
        "843f311ffffffff",
        "843f207ffffffff",
        "843f22bffffffff",
        "843f317ffffffff",
        "843f353ffffffff",
        "843f22dffffffff",
        "843f31dffffffff",
        "843f313ffffffff",
        "843f23dffffffff",
        "843f225ffffffff",
        "843f31bffffffff",
        "843f047ffffffff",
        "843f315ffffffff",
        "843f223ffffffff",
        "843f043ffffffff",
        "843f04bffffffff",
        "843f203ffffffff",
        "843f319ffffffff",
        "843f239ffffffff",
        "843f221ffffffff",
        "843f263ffffffff",
    ]
    assert cover_shape_exact(
        mpolygon, resolution=4, precision_resolution=5, inside_only=False
    ) == [
        "843f235ffffffff",
        "843f2e9ffffffff",
        "843f049ffffffff",
        "843f229ffffffff",
        "843f04dffffffff",
        "843f35bffffffff",
        "843f23bffffffff",
        "843f2c5ffffffff",
        "843f233ffffffff",
        "843f2e5ffffffff",
        "843f265ffffffff",
        "843f267ffffffff",
        "843f059ffffffff",
        "843f041ffffffff",
        "843f205ffffffff",
        "843f2e1ffffffff",
        "843f045ffffffff",
        "843f207ffffffff",
        "843f311ffffffff",
        "843f22bffffffff",
        "843f217ffffffff",
        "843f353ffffffff",
        "843f22dffffffff",
        "843f317ffffffff",
        "843f31dffffffff",
        "843f313ffffffff",
        "843f23dffffffff",
        "843f2ebffffffff",
        "843f225ffffffff",
        "843f31bffffffff",
        "843f2edffffffff",
        "843f2e7ffffffff",
        "843f047ffffffff",
        "843f315ffffffff",
        "843f043ffffffff",
        "843f223ffffffff",
        "843f04bffffffff",
        "843f203ffffffff",
        "843f319ffffffff",
        "843f239ffffffff",
        "843f221ffffffff",
        "843f263ffffffff",
    ]


def test_index_to_point():
    point = index_to_point("84538abffffffff")
    wkt = [p.wkt for p in point]
    assert wkt == []

    point2 = index_to_point(["84538abffffffff", "84539b5ffffffff"])
    wkt2 = [p.wkt for p in point2]
    assert wkt2 == []

    point3 = index_to_point({"84538abffffffff", "84539b5ffffffff"})
    wkt3 = [p.wkt for p in point3]
    assert wkt3 == []

    point4 = index_to_point(pd.Series(["84538abffffffff", "84539b5ffffffff"]))
    wkt4 = [p.wkt for p in point4]
    assert wkt4 == []


def test_index_to_polygon():
    polygon = index_to_polygon("82f387fffffffff")
    wkt = [p.wkt for p in polygon]
    assert wkt == []

    am_polygon = index_to_polygon("82f387fffffffff", antimeridian_cutting=True)
    am_wkt = [p.wkt for p in am_polygon]
    assert am_wkt == []
