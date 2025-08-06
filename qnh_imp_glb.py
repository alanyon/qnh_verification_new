"""
Script to load the IMPROVER data and find the lowest mean sea level
pressure (MSLP) in each FOQNH region for the 50th percentile.

Functions:
    main: Main function calling other functions.
    find_lowest_mslp_in_polygon: Finds the lowest MSLP in given polygon.
    get_cube: Gets the IMPROVER data for the specified cycle point.
    get_regions: Gets the FOQNH regions from the shapefile.
"""
import os
import sys
from datetime import datetime, timedelta

import iris
import numpy as np
import shapefile as shp
from shapely.geometry import Point, Polygon

# Define environment constants
DATA_DIR = os.environ['DATA_DIR']
CYCLE_POINT = os.environ['CYCLE_POINT']
MOO_DIR = os.environ['MOO_DIR']

# Other constants
EX_DIR = f'{DATA_DIR}/grid_global/{CYCLE_POINT}'
EXTENT = [-10, 3, 48.5, 63.5]
REGION_ORDER = ['SKERRY', 'PORTREE', 'RATTRAY', 'TYNE', 'BELFAST', 'HOLYHEAD',
                'BARNSLEY', 'HUMBER', 'SCILLIES', 'WESSEX', 'CHATHAM',
                'PORTLAND', 'YARMOUTH', 'COTSWOLD', 'SHETLAND', 'ORKNEY',
                'MARLIN', 'PETREL', 'SKUA', 'PUFFIN']

# Stop iris future warnings
iris.FUTURE.date_microseconds = True


def main():
    """
    Main function to load the IMPROVER data and find lowest MSLP in
    each FOQNH region.

    Args:
        None
    Returns:
        None
    """
    # Get FOQNH region coordinates
    regions = get_regions()

    # Need valid time string 5 hours after cycle point
    cycle_dt = datetime.strptime(CYCLE_POINT, '%Y%m%dT%H%MZ')
    vt_dts = [cycle_dt + timedelta(hours=5), cycle_dt + timedelta(hours=6)]
    vt_str = vt_dts[0].strftime('%Y%m%d%H%MZ')

    # Extract data from MASS
    cube = get_cube(vt_dts)

    # Exit if no cube found
    if cube is None:
        sys.exit()

    # Collect data for all low percentiles
    for perc in [5, 10, 15, 20, 30, 40, 50]:

        # Get cube just required percentile
        cube_perc = cube.extract(iris.Constraint(percentile=perc))

        # Loop through each region
        for coords in regions.values():

            # Create a polygon from the coordinates
            polygon = Polygon(coords)

            # Find the lowest MSLP in the polygon for that percentile
            lowest_mslp = find_lowest_mslp_in_polygon(cube_perc, polygon)

            # Add min MSLP to a .dat file
            dat_fname = f'{DATA_DIR}/min_mslps/imp_{perc}_QNH_{vt_str2}.dat'
            with open(dat_fname, 'a', encoding='utf-8') as file:
                file.write(f'{lowest_mslp}\n')

            # Also do the min - 3 for 50th percentile
            if perc == 50:
                dat_fname = (f'{DATA_DIR}/min_mslps/imp_{perc}_minus_3_QNH_'
                             f'{vt_str}.dat')
                with open(dat_fname, 'a', encoding='utf-8') as file:
                    file.write(f'{lowest_mslp - 3}\n')

    # Remove extracted data
    os.system(f'rm -rf {EX_DIR}')


def find_lowest_mslp_in_polygon(cube, polygon):
    """
    Finds the lowest mean sea level pressure (MSLP) in a given polygon
    within the cube.

    Args:
        cube (iris.cube.Cube): The input cube containing MSLP data.
        polygon (shapely.geometry.Polygon): The polygon to check.
    Returns:
        float: The lowest MSLP value found within the polygon.
    """
    # Get latitude and longitude coordinates and flatten them
    lat = cube.coord('latitude').points
    lon = cube.coord('longitude').points
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    lat_flat = lat_grid.flatten()
    lon_flat = lon_grid.flatten()

    # Get data and flatten it
    data_flat = cube.data.flatten()

    # Create a list of values that are within the polygon
    values_in_polygon = [
        data for data, x, y in zip(data_flat, lon_flat, lat_flat)
        if polygon.contains(Point(x, y))
    ]

    # Find minimum value in the polygon (as an integer)
    min_val = int(min(values_in_polygon))

    return min_val


def get_cube(vt_dts):
    """
    Gets the IMPROVER data for the specified cycle point and returns as
    a cube.

    Args:
        vt_dts (list): List of valid time datetimes.
    Returns:
        cube(iris.cube.Cube): Cube containing the IMPROVER data.
    """
    # Get valid time strs
    vt_strs = [dt.strftime('%Y%m%dT%H%MZ') for dt in vt_dts]

    # Create directory to put data in
    if not os.path.exists(EX_DIR):
        os.makedirs(EX_DIR)

    # Ensure directory exists
    for os_num in ['OS45.2', 'OS46']:
        m_path = f'{MOO_DIR}/{os_num}/engl_suite_{CYCLE_POINT}/grid.tar'
        if os.system(f'moo ls {m_path} > /dev/null 2>&1') == 0:

            # If file exists in MASS, extract it
            os.system(f'moo get {m_path} {EX_DIR}/')

            # Untar required file
            mslp_fname_1 = (f'grid/percentile_extract_{vt_strs[0]}-'
                            'PT0005H00M-pressure_at_mean_sea_level.nc')
            mslp_fname_2 = (f'grid/percentile_extract_{vt_strs[1]}-'
                            'PT0006H00M-pressure_at_mean_sea_level.nc')
            os.system(f'tar -xvf {EX_DIR}/grid.tar -C {EX_DIR} {mslp_fname_1}')
            os.system(f'tar -xvf {EX_DIR}/grid.tar -C {EX_DIR} {mslp_fname_2}')

            # Load the data using iris
            cube_1 = iris.load_cube(f'{EX_DIR}/{mslp_fname_1}')
            cube_2 = iris.load_cube(f'{EX_DIR}/{mslp_fname_2}')

            # Constrain cube using extents
            cube_1 = cube_1.intersection(longitude=(EXTENT[0], EXTENT[1]),
                                         latitude=(EXTENT[2], EXTENT[3]))
            cube_2 = cube_2.intersection(longitude=(EXTENT[0], EXTENT[1]),
                                         latitude=(EXTENT[2], EXTENT[3]))

            # Convert units to hPa
            cube_1.convert_units('hPa')
            cube_2.convert_units('hPa')

            # Get cube with minimum MSLPs from the two cubes
            cube = iris.cube.CubeList([cube_1, cube_2]).merge_cube()
            cube = cube.collapsed(['longitude', 'latitude'], iris.analysis.MIN)

            return cube

    # If no file found, return None
    print(f'No IMPROVER data found for cycle point {CYCLE_POINT}')
    return None


def get_regions():
    """
    Gets the FOQNH regions from the shapefile and returns as a
    dictionary.

    Args:
        None
    Returns:
        regions(dict): Dictionary of region names and their coordinates.
    """
    # Load the shapefile containing the FOQNH regions
    shapes = shp.Reader(f'{DATA_DIR}/MO_ForecastQNHRegions.shp')

    # Create a dictionary to hold the region names and their coordinates
    regions = {}
    for shape_rec in shapes.shapeRecords():
        name = shape_rec.record[0]
        shape = shape_rec.shape
        coords = shape.points
        regions[name] = coords

    # Check names are all there and in same order used in verification
    msg = 'Region names do not match expected order or are missing.'
    assert set(regions.keys()) == set(REGION_ORDER), msg

    return regions


if __name__ == "__main__":
    t1 = datetime.now()
    main()
    t2 = datetime.now()
    # print time in seconds
    time_taken = (t2 - t1).total_seconds()
    print(f'Finished in {time_taken} seconds')
