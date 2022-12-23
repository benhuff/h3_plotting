from shapely.geometry import (
    Point,
    MultiPoint,
    LineString,
    MultiLineString,
    Polygon,
    MultiPolygon,
)
from shapely.geometry.base import BaseGeometry
import h3


def _cover_point(shape: Point, resolution: int, buffer: int = None):

    hexes = set()

    hexes.update([h3.geo_to_h3(shape.x, shape.y, resolution=resolution)])

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.hex_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    return hexes


def _cover_multipoint(shape: MultiPoint, resolution: int, buffer: int = None):

    hexes = set()

    for geom in shape.geoms:
        hexes.update([h3.geo_to_h3(geom.x, geom.y, resolution=6)])

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.hex_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    return hexes


def _cover_linestring(shape: LineString, resolution: int, buffer: int = None):

    hexes = set()

    endpoint_hexes = [h3.geo_to_h3(t[0], t[1], resolution) for t in list(shape.coords)]

    for i in range(len(endpoint_hexes) - 1):
        hexes.update(h3.h3_line(endpoint_hexes[i + 1], endpoint_hexes[i]))

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    return hexes


def _cover_multilinestring(shape: MultiLineString, resolution: int, buffer: int = None):

    hexes = set()

    for geom in shape.geoms:
        endpoint_hexes = [
            h3.geo_to_h3(t[0], t[1], resolution) for t in list(geom.coords)
        ]

        for i in range(len(endpoint_hexes) - 1):
            hexes.update(h3.h3_line(endpoint_hexes[i + 1], endpoint_hexes[i]))

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    return hexes


def _cover_polygon(shape: Polygon, resolution: int, buffer: int = None):

    exteriors = set()
    interiors = set()

    exteriors.update(
        h3.polyfill(Polygon(shape.exterior).__geo_interface__, res=resolution)
    )
    for interior in shape.interiors:
        interiors.update(
            h3.polyfill(Polygon(interior).__geo_interface__, res=resolution)
        )

    hexes = exteriors - interiors

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    return hexes


def _cover_multipolygon(shape: MultiPolygon, resolution: int, buffer: int = None):

    exteriors = set()
    interiors = set()

    for geom in shape.geoms:
        exteriors.update(h3.polyfill(Polygon(geom).__geo_interface__, res=resolution))
        for interior in geom.interiors:
            interiors.update(
                h3.polyfill(Polygon(interior).__geo_interface__, res=resolution)
            )

    hexes = exteriors - interiors

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    return hexes


def cover_shape(shape: BaseGeometry, resolution: int, buffer: int = None):
    try:
        if isinstance(shape, Point):
            return _cover_point(shape=shape, resolution=resolution, buffer=buffer)

        elif isinstance(shape, MultiPoint):
            return _cover_multipoint(shape=shape, resolution=resolution, buffer=buffer)

        elif isinstance(shape, LineString):
            return _cover_linestring(shape=shape, resolution=resolution, buffer=buffer)

        elif isinstance(shape, MultiLineString):
            return _cover_multilinestring(
                shape=shape, resolution=resolution, buffer=buffer
            )

        elif isinstance(shape, Polygon):
            return _cover_polygon(shape=shape, resolution=resolution, buffer=buffer)

        elif isinstance(shape, MultiPolygon):
            return _cover_multipolygon(
                shape=shape, resolution=resolution, buffer=buffer
            )

        else:
            raise Exception("Unknown Shapely geometry type!")

    except:
        raise Exception(
            """
            Unknown Shapely geometry type! Currently supported geometry types: 
            Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon.
            """
        )
