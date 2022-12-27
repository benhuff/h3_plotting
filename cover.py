from shapely.geometry import (
    Point,
    MultiPoint,
    LineString,
    MultiLineString,
    Polygon,
    MultiPolygon,
)
from shapely.ops import split, transform

import h3
import numpy as np

from typing import Union, Literal, get_args
from pandas import Series
import warnings

_index_types = Literal["int", "str"]


def _cover_point(shape: Point, resolution: int, buffer: int = None):
    """
    Cover any Shapely Point geometry with H3 cell indexes at a specified H3 resolution.
    Optionally, buffer the H3 cell indexes with any number of additional k-rings.

    Args:
        shape (Point): Shapely Point geometry.
        resolution (int): H3 index resolution.
        buffer (int, optional): Number of additional k-rings. Defaults to None.

    Returns:
        list: List of H3 indexes.
    """
    hexes = set()

    hexes.update([h3.geo_to_h3(shape.y, shape.x, resolution=resolution)])

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    if len(hexes) == 0:
        warnings.warn(
            f"There were no H3 indexes for {shape} at resolution {resolution}"
        )

    return list(hexes)


def _cover_multipoint(shape: MultiPoint, resolution: int, buffer: int = None):
    """
    Cover any Shapely MultiPoint geometry with H3 cell indexes at a specified H3 resolution.
    Optionally, buffer the H3 cell indexes with any number of additional k-rings.

    Args:
        shape (MultiPoint): Shapely MultiPoint geometry.
        resolution (int): H3 index resolution.
        buffer (int, optional): Number of additional k-rings. Defaults to None.

    Returns:
        list: List of H3 indexes.
    """
    hexes = set()

    for geom in shape.geoms:
        hexes.update([h3.geo_to_h3(geom.y, geom.x, resolution=resolution)])

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    if len(hexes) == 0:
        warnings.warn(
            f"There were no H3 indexes for {shape} at resolution {resolution}"
        )

    return list(hexes)


def _cover_linestring(shape: LineString, resolution: int, buffer: int = None):
    """
    Cover any Shapely LineString geometry with H3 cell indexes at a specified H3 resolution.
    Optionally, buffer the H3 cell indexes with any number of additional k-rings.

    Args:
        shape (LineString): Shapely LineString geometry.
        resolution (int): H3 index resolution.
        buffer (int, optional): Number of additional k-rings. Defaults to None.

    Returns:
        list: List of H3 indexes.
    """
    hexes = set()

    endpoint_hexes = [h3.geo_to_h3(t[1], t[0], resolution) for t in list(shape.coords)]

    for i in range(len(endpoint_hexes) - 1):
        hexes.update(h3.h3_line(endpoint_hexes[i], endpoint_hexes[i + 1]))

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    if len(hexes) == 0:
        warnings.warn(
            f"There were no H3 indexes for {shape} at resolution {resolution}"
        )

    return list(hexes)


def _cover_multilinestring(shape: MultiLineString, resolution: int, buffer: int = None):
    """
    Cover any Shapely MultiLineString geometry with H3 cell indexes at a specified H3 resolution.
    Optionally, buffer the H3 cell indexes with any number of additional k-rings.

    Args:
        shape (MultiLineString): Shapely MultiLineString geometry.
        resolution (int): H3 index resolution.
        buffer (int, optional): Number of additional k-rings. Defaults to None.

    Returns:
        list: List of H3 indexes.
    """
    hexes = set()

    for geom in shape.geoms:
        endpoint_hexes = [
            h3.geo_to_h3(t[1], t[0], resolution) for t in list(geom.coords)
        ]

        for i in range(len(endpoint_hexes) - 1):
            hexes.update(h3.h3_line(endpoint_hexes[i], endpoint_hexes[i + 1]))

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    if len(hexes) == 0:
        warnings.warn(
            f"There were no H3 indexes for {shape} at resolution {resolution}"
        )

    return list(hexes)


def _cover_polygon(shape: Polygon, resolution: int, buffer: int = None):
    """
    Cover any Shapely Polygon geometry with H3 cell indexes at a specified H3 resolution.
    Optionally, buffer the H3 cell indexes with any number of additional k-rings.

    Args:
        shape (Polygon): Shapely Polygon geometry.
        resolution (int): H3 index resolution.
        buffer (int, optional): Number of additional k-rings. Defaults to None.

    Returns:
        list: List of H3 indexes.
    """
    exteriors = set()
    interiors = set()

    exteriors.update(
        h3.polyfill(Polygon(shape.exterior).__geo_interface__, res=resolution, geo_json_conformant=True)
    )
    for interior in shape.interiors:
        interiors.update(
            h3.polyfill(Polygon(interior).__geo_interface__, res=resolution, geo_json_conformant=True)
        )

    hexes = exteriors - interiors

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    if len(hexes) == 0:
        warnings.warn(
            f"There were no H3 indexes for {shape} at resolution {resolution}"
        )

    return list(hexes)


def _cover_multipolygon(shape: MultiPolygon, resolution: int, buffer: int = None):
    """
    Cover any Shapely MultiPolygon geometry with H3 cell indexes at a specified H3 resolution.
    Optionally, buffer the H3 cell indexes with any number of additional k-rings.

    Args:
        shape (MultiPolygon): Shapely MultiPolygon geometry.
        resolution (int): H3 index resolution.
        buffer (int, optional): Number of additional k-rings. Defaults to None.

    Returns:
        list: List of H3 indexes.
    """
    exteriors = set()
    interiors = set()

    for geom in shape.geoms:
        exteriors.update(h3.polyfill(Polygon(geom).__geo_interface__, res=resolution, geo_json_conformant=True))
        for interior in geom.interiors:
            interiors.update(
                h3.polyfill(Polygon(interior).__geo_interface__, res=resolution, geo_json_conformant=True)
            )

    hexes = exteriors - interiors

    if buffer:
        buffer_hexes = set()

        [buffer_hexes.update(h3.k_ring(h, buffer)) for h in list(hexes)]
        hexes = hexes.union(buffer_hexes)

    if len(hexes) == 0:
        warnings.warn(
            f"There were no H3 indexes for {shape} at resolution {resolution}"
        )

    return list(hexes)


def cover_shape(
    shape: Union[Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon],
    resolution: int,
    buffer: int = None,
    index_type: _index_types = "str",
):
    """
    Cover any Shapely geometry with H3 cell indexes at a specified H3 resolution.
    Optionally, buffer the H3 cell indexes with any number of additional k-rings.

    Args:
        shape (Union[Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon]): Shapely geometry.
        resolution (int): H3 index resolution.
        buffer (int, optional): Number of additional k-rings. Defaults to None.
        index_type (str, optional): String equal to 'str' or 'int' denoting either H3 integer indexes or H3 string indexes. Defaults to 'int'.

    Raises:
        Exception: Unknown Shapely geometry type.

    Returns:
        list or str or int: List of H3 indexes of type `index_type`.

    Examples:
        One use case would be to find the H3 cell indexes that cover a Shapely MultiPolygon geometry:
        ```
        >>> print(point)
        POINT (18.99551192594893 34.93768898538498)
        >>> cover.cover_shape(point, 6, index_type='int') # get H3 integer index
        [604951465911386111]
        >>> cover.cover_shape(point, 6, index_type='str') # get H3 string index
        ['86538272fffffff']
        ```
    """
    options = get_args(_index_types)
    assert (
        index_type in options
    ), f'{index_type} is not a supported index type. "str" and "int" are supported.'

    if isinstance(shape, Point):
        indexes = _cover_point(shape=shape, resolution=resolution, buffer=buffer)
        if index_type == "int":
            indexes = [h3.string_to_h3(i) for i in indexes]
        return indexes

    elif isinstance(shape, MultiPoint):
        indexes = _cover_multipoint(shape=shape, resolution=resolution, buffer=buffer)
        if index_type == "int":
            indexes = [h3.string_to_h3(i) for i in indexes]
        return indexes

    elif isinstance(shape, LineString):
        indexes = _cover_linestring(shape=shape, resolution=resolution, buffer=buffer)
        if index_type == "int":
            indexes = [h3.string_to_h3(i) for i in indexes]
        return indexes

    elif isinstance(shape, MultiLineString):
        indexes = _cover_multilinestring(
            shape=shape, resolution=resolution, buffer=buffer
        )
        if index_type == "int":
            indexes = [h3.string_to_h3(i) for i in indexes]
        return indexes

    elif isinstance(shape, Polygon):
        indexes = _cover_polygon(shape=shape, resolution=resolution, buffer=buffer)
        if index_type == "int":
            indexes = [h3.string_to_h3(i) for i in indexes]
        return indexes

    elif isinstance(shape, MultiPolygon):
        indexes = _cover_multipolygon(shape=shape, resolution=resolution, buffer=buffer)
        if index_type == "int":
            indexes = [h3.string_to_h3(i) for i in indexes]
        return indexes

    else:
        raise Exception(
            f"""
            Unknown Shapely geometry type! Currently supported geometry types: 
            Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon.
            {type(shape)} was provided.
            """
        )


def _cover_polygon_exact(
    shape: Polygon,
    resolution: int,
    precision_resolution: int = None,
    inside_only: bool = True,
):
    """
    Cover any Shapely Polygon geometry with H3 cell indexes at a specified H3 resolution.
    A precision H3 resolution is used to determine whether or not to include boundary H3 indexes.
    Whether to include these boundary H3 indexes is determined by the inside_only argument.

    Args:
        shape (Polygon): Shapely Polygon geometry.
        resolution (int): H3 index resolution.
        precision_resolution (int): H3 index resolution to use for precision determination of whether or not to include boundary H3 indexes.
        inside_only (bool, optional): Whether to include boundary H3 indexes is determined by the inside_only argument. Defaults to True.
        buffer (int, optional): Number of additional k-rings. Defaults to None.

    Returns:
        list: List of H3 indexes.
    """
    exteriors = set()
    exterior_boundaries = set()
    interiors = set()
    interior_boundaries = set()

    exteriors.update(cover_shape(shape.representative_point(), resolution=resolution))
    exteriors.update(cover_shape(Polygon(shape.exterior), resolution=resolution))
    exterior_boundaries.update(
        cover_shape(shape.boundary, resolution=precision_resolution)
    )
    for interior in shape.interiors:
        interiors.update(cover_shape(Polygon(interior), resolution=resolution))
        interior_boundaries.update(
            cover_shape(interior, resolution=precision_resolution)
        )

    hexes = exteriors - interiors

    boundaries = exterior_boundaries.union(interior_boundaries)
    boundaries = [h3.h3_to_parent(i, resolution) for i in list(boundaries)]

    if inside_only:
        hexes = hexes - set(boundaries)
    else:
        hexes = hexes.union(set(boundaries))

    if len(hexes) == 0:
        warnings.warn(
            f"There were no H3 indexes for {shape} at resolution {resolution}"
        )

    return list(hexes)


def _cover_multipolygon_exact(
    shape: MultiPolygon,
    resolution: int,
    precision_resolution: int = None,
    inside_only: bool = True,
):
    """
    Cover any Shapely MultiPolygon geometry with H3 cell indexes at a specified H3 resolution.
    A precision H3 resolution is used to determine whether or not to include boundary H3 indexes.
    Whether to include these boundary H3 indexes is determined by the inside_only argument.

    Args:
        shape (MultiPolygon): Shapely MultiPolygon geometry.
        resolution (int): H3 index resolution.
        precision_resolution (int): H3 index resolution to use for precision determination of whether or not to include boundary H3 indexes.
        inside_only (bool, optional): Whether to include boundary H3 indexes is determined by the inside_only argument. Defaults to True.
        buffer (int, optional): Number of additional k-rings. Defaults to None.

    Returns:
        list: List of H3 indexes.
    """
    exteriors = set()
    exterior_boundaries = set()
    interiors = set()
    interior_boundaries = set()

    for geom in shape.geoms:
        exteriors.update(
            cover_shape(geom.representative_point(), resolution=resolution)
        )
        exteriors.update(cover_shape(Polygon(geom.exterior), resolution=resolution))
        exterior_boundaries.update(
            cover_shape(geom.boundary, resolution=precision_resolution)
        )
        for interior in geom.interiors:
            interiors.update(cover_shape(Polygon(interior), resolution=resolution))
            interior_boundaries.update(
                cover_shape(interior, resolution=precision_resolution)
            )

    hexes = exteriors - interiors

    boundaries = exterior_boundaries.union(interior_boundaries)
    boundaries = [h3.h3_to_parent(i, resolution) for i in list(boundaries)]

    if inside_only:
        hexes = hexes - set(boundaries)
    else:
        hexes = hexes.union(set(boundaries))

    if len(hexes) == 0:
        warnings.warn(
            f"There were no H3 indexes for {shape} at resolution {resolution}"
        )

    return list(hexes)


def cover_shape_exact(
    shape: Union[Polygon, MultiPolygon],
    resolution: int,
    precision_resolution: int,
    inside_only: bool = True,
    index_type: _index_types = "str",
):
    """
    Cover any Shapely Polygon or MultiPolygon geometry with H3 cell indexes at a specified H3 resolution.
    A precision H3 resolution is used to determine whether or not to include boundary H3 indexes.
    Whether to include these boundary H3 indexes is determined by the inside_only argument.

    Args:
        shape (Union[Polygon, MultiPolygon]): Shapely geometry.
        resolution (int): H3 index resolution.
        precision_resolution (int): H3 index resolution to use for precision determination of whether or not to include boundary H3 indexes.
        inside_only (bool, optional): Whether to include boundary H3 indexes is determined by the inside_only argument. Defaults to True.
        index_type (str, optional): String equal to 'str' or 'int' denoting either H3 integer indexes or H3 string indexes. Defaults to 'int'.

    Raises:
        Exception: Unknown Shapely geometry type.

    Returns:
        list or str or int: List of H3 indexes of type `index_type`.

    Examples:
        One use case would be to find the H3 cell indexes that cover a Shapely MultiPolygon geometry:
        ```
        >>> print(polygon)
        POLYGON ((15.06998710098651 37.51339776673665, 15.07016793307378 37.51265325152787, 15.0717265334375 37.5127830299427, 15.07128306427364 37.51366756628664, 15.06998710098651 37.51339776673665))
        >>> cover.cover_shape_exact(polygon, resolution=10, precision_resolution=14, inside_only=False)
        ['8a52a2ccc297fff', '8a52a2ccc667fff', '8a52a2ccc74ffff', '8a52a2ccc66ffff']
        >>> cover.cover_shape_exact(polygon, resolution=10, precision_resolution=14, inside_only=False, index_type='int')
        [622950495352881151, 622950495356878847, 622950495357829119, 622950495356911615]
        ```
    """
    options = get_args(_index_types)
    assert (
        index_type in options
    ), f'{index_type} is not a supported index type. "str" and "int" are supported.'
    assert (
        precision_resolution > resolution
    ), f"The precision_resolution must be greater than the resolution. {resolution} is not greater than {precision_resolution}."

    if isinstance(shape, Polygon):
        indexes = _cover_polygon_exact(
            shape=shape,
            resolution=resolution,
            precision_resolution=precision_resolution,
            inside_only=inside_only,
        )

        if index_type == "int":
            indexes = [h3.string_to_h3(i) for i in indexes]
        return indexes

    elif isinstance(shape, MultiPolygon):
        indexes = _cover_multipolygon_exact(
            shape=shape,
            resolution=resolution,
            precision_resolution=precision_resolution,
            inside_only=inside_only,
        )
        if index_type == "int":
            indexes = [h3.string_to_h3(i) for i in indexes]
        return indexes

    else:
        raise Exception(
            f"""
            Unknown Shapely geometry type! Currently supported geometry types: 
            Polygon, MultiPolygon.
            {type(shape)} was provided.
            """
        )


def antimeridian_handler(shape: Polygon):
    """
    Cut a Shapely Polygon that cross the antimeridian ((180, 90), (180, -90)) into a MultiPolygon for visualizing in WGS84.
    If the given Polygon does not cross the

    Args:
        shape (Polygon): Shapely Polygon.

    Returns:
        Polygon or MultiPolygon: Shapely Polygon or MultiPolygon.
    """
    x = np.array(shape.exterior.coords.xy[0])
    y = np.array(shape.exterior.coords.xy[1])

    idxs_neg = np.where(x < 0)
    x[idxs_neg] += 360

    extended = Polygon(np.array(list(zip(x, y))))
    am = LineString(((180, 90), (180, -90)))

    if extended.intersects(am):
        return MultiPolygon(split(extended, am))
    else:
        return shape


def index_to_polygon(
    index: Union[list, set, Series, str, int], antimeridian_cutting: bool = False
):
    """
    Convert one or more H3 cell indexes into Shapely Polygons representing the geometry of the H3 cell boundary.

    Args:
        index (Union[list, set, Series, str, int]): One or more H3 cell indexes.

    Raises:
        Exception: Unsupported index input type.

    Returns:
        list: list of Shapely Polygon geometries. CRS of returned geometries is EPSG:4326.

    Examples:
        One use case would be to create a GeoDataFrame from one or more H3 cell indexes:
        ```
        >>> index = {'885391a507fffff', '88539c02cdfffff', '88539ccb47fffff'}
        >>> gdf = gpd.GeoDataFrame(geometry=cover.index_to_point(index), crs='epsg:4326')
        ```

        Another use case could be to add a geometry column to a DataFrame that has a column with H3 cell indexes:
        ```
        >>> df.loc[:,'geometry'] = cover.index_to_point(df['h3_index'])
        ```
    """
    if isinstance(index, list):
        if isinstance(index[0], int):
            index = [h3.h3_to_string(i) for i in index]
        index = index

    elif isinstance(index, set):
        index = list(index)
        if isinstance(index[0], int):
            index = [h3.h3_to_string(i) for i in index]
        index = index

    elif isinstance(index, Series):
        index = index.to_list()
        if isinstance(index[0], int):
            index = [h3.h3_to_string(i) for i in index]
        index = index

    elif isinstance(index, str):
        index = [index]

    elif isinstance(index, int):
        index = [h3.h3_to_string(index)]

    else:
        raise Exception(
            f"""
            Index is of unknown type! Currently supported types:
            List, Set, Series, String, Integer.
            {type(index)} was provided.
            """
        )

    if antimeridian_cutting:
        polygons = [
            antimeridian_handler(Polygon(h3.h3_to_geo_boundary(i, geo_json=True))) for i in index
        ]

    else:
        polygons = [Polygon(h3.h3_to_geo_boundary(i, geo_json=True)) for i in index]

    return polygons


def flip(x, y):
    "Swap x, y order."
    return y, x


def index_to_point(index: Union[list, set, Series, str, int]):
    """
    Convert one or more H3 cell indexes into Shapely Points representing the centroid of the H3 cell boundary.

    Args:
        index (Union[list, set, Series, str, int]): One or more H3 cell indexes.

    Raises:
        Exception: Unsupported index input type.

    Returns:
        list: list of Shapely Point geometries. CRS of returned geometries is EPSG:4326.

    Examples:
        One use case would be to create a GeoDataFrame from one or more H3 cell indexes:
        ```
        >>> index = {'885391a507fffff', '88539c02cdfffff', '88539ccb47fffff'}
        >>> gdf = gpd.GeoDataFrame(geometry=cover.index_to_point(index), crs='epsg:4326')
        ```

        Another use case could be to add a geometry column to a DataFrame that has a column with H3 cell indexes:
        ```
        >>> df.loc[:,'geometry'] = cover.index_to_point(df['h3_index'])
        ```
    """
    if isinstance(index, list):
        if isinstance(index[0], int):
            index = [h3.h3_to_string(i) for i in index]
        index = index

    elif isinstance(index, set):
        index = list(index)
        if isinstance(index[0], int):
            index = [h3.h3_to_string(i) for i in index]
        index = index

    elif isinstance(index, Series):
        index = index.to_list()
        if isinstance(index[0], int):
            index = [h3.h3_to_string(i) for i in index]
        index = index

    elif isinstance(index, str):
        index = [index]

    elif isinstance(index, int):
        index = [h3.h3_to_string(index)]

    else:
        raise Exception(
            f"""
            Index is of unknown type! Currently supported types:
            List, Set, Series, String, Integer.
            {type(index)} was provided.
            """
        )

    return [transform(flip, Point(h3.h3_to_geo(i))) for i in index]
