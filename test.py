import iris
import cartopy.feature as feature
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import shapefile as shp

import iris.plot as iplt

# Define data directory
DATA_DIR = '/data/users/andre.lanyon/QNH'
# Other constants needed
IMP_DIFFS = {
    'uk': {'grid_dir': 'grid_uk',
           'dt_str': '20250414T1500Z-B20250414T1230Z',
           'x_coord': 'projection_x_coordinate',
           'y_coord': 'projection_y_coordinate',
           'title': 'UK IMPROVER'},
    'global': {'grid_dir': 'grid_global',
               'dt_str': '20250414T1500Z-PT0003H00M',
               'x_coord': 'longitude',
               'y_coord': 'latitude',
               'title': 'Global IMPROVER'}
}
TARGET_CRS = ccrs.PlateCarree()
PLOT_EXTENT = [-15, 5, 48, 65]
# Stop iris future warnings
iris.FUTURE.date_microseconds = True


def main():

    # Get FOQNH region coordinates
    regions = get_regions()

    # Get 50th percentile cubes for UK and global IMPROVER data and find
    # min and max MSLPs to ensure same colour scale for plotting
    cube_50ths = []
    vmin, vmax = float('inf'), float('-inf')
    for imp_type, imp_data in IMP_DIFFS.items():
    
        # Load the cube for the given IMPROVER type
        dt_str = imp_data['dt_str']
        grid_dir = imp_data['grid_dir']
        fname = f'percentile_extract_{dt_str}-pressure_at_mean_sea_level.nc'
        cube = iris.load_cube(f'{DATA_DIR}/{grid_dir}/{fname}')

        # Convert units to hPa
        cube.convert_units('hPa')

        # Convert lat/lon coords
        cube = convert_coords(cube, imp_type, imp_data)

        # Get cube just at 50th percentile and plot extents
        cube_50th = cube.extract(iris.Constraint(percentile=10))

        lat_aux = cube_50th.coord('latitude_platecarree').points
        lon_aux = cube_50th.coord('longitude_platecarree').points

        # Create a mask for the desired extent
        mask = (
            (lon_aux >= PLOT_EXTENT[0]) & (lon_aux <= PLOT_EXTENT[1]) &
            (lat_aux >= PLOT_EXTENT[2]) & (lat_aux <= PLOT_EXTENT[3])
        )

        # Mask the data
        masked_data = np.ma.masked_where(~mask, cube_50th.data)

        # Create a new cube with the masked data
        cube_50th.data = masked_data

        cube_50ths.append(cube_50th)

        # Update global min/max
        data_min = cube_50th.data.min()
        data_max = cube_50th.data.max()
        vmin = min(vmin, data_min)
        vmax = max(vmax, data_max)
    
    # Create figure and axis with a Transverse Mercator projection
    map_crs = ccrs.TransverseMercator()
    fig, axs = plt.subplots(1, 2, figsize=(25, 20), 
                            subplot_kw={'projection': map_crs})
    
    # Plot each cube on same colour scale
    for ax, cube, imp_data in zip(axs, cube_50ths, IMP_DIFFS.values()):

        # Plot data on map using latitude and longitude coordinates
        mapable = iplt.pcolormesh(cube, cmap='jet', vmin=vmin, vmax=vmax, 
                                  axes=ax)

        # Draw coastlines and background stuff
        ax.add_feature(feature.COASTLINE,linewidth=0.8)
        ax.add_feature(feature.BORDERS, linestyle=':',linewidth=0.4)
        ax.set_extent([-10, 3, 48.5, 63.5], crs=TARGET_CRS)
        ax.gridlines(draw_labels=True, linewidth=0.4)

        # Add title
        ax.set_title(f'{imp_data["title"]}', fontsize=28, weight='bold')

        # Plot FOQNH areas as polygons
        for area_name, coords in regions.items():

            # Create a polygon from the coordinates
            polygon = Polygon(coords)

            # Define the centroid and exterior of the polygon
            centroid = polygon.centroid
            x, y = polygon.exterior.xy

            # Find the lowest MSLP in the polygon
            lowest_mslp = find_lowest_mslp_in_polygon(cube, polygon)

            # Plot the polygon and annotate with area name and lowest MSLP
            ax.plot(x, y, color='red', linewidth=1, transform=TARGET_CRS)
            txt = f'{area_name}\n{lowest_mslp:.1f} hPa'
            ax.text(centroid.x, centroid.y, txt, 
                fontsize=12, ha='center', va='center', color='white',
                fontweight='bold',
                transform=TARGET_CRS)

    # Add colorbar
    plt.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.02])
    cbar = plt.colorbar(mapable, cax=cbar_ax, orientation='horizontal', pad=0.05, 
                        aspect=50, shrink=0.8)
    cbar.set_label('Pressure at Mean Sea Level (hPa)', fontsize=22)
    cbar.ax.tick_params(labelsize=16)


    # Save and close plot
    fig.savefig('cube_50th_plot_10_perc.png', bbox_inches='tight')
    plt.close()


def convert_coords(cube, imp_type, imp_data):

    if imp_type == 'uk':

        coord_system = cube.coord('projection_x_coordinate').coord_system

        source_crs = ccrs.LambertAzimuthalEqualArea(
            central_latitude=coord_system.latitude_of_projection_origin,
            central_longitude=coord_system.longitude_of_projection_origin,
            false_easting=coord_system.false_easting,
            false_northing=coord_system.false_northing,
            globe=ccrs.Globe(
                semimajor_axis=coord_system.ellipsoid.semi_major_axis,
                semiminor_axis=coord_system.ellipsoid.semi_minor_axis
            )
        )

    else:

        source_crs = ccrs.Geodetic(globe=ccrs.Globe(semimajor_axis=6371229.0, 
                                                    semiminor_axis=6371229.0))

    x = cube.coord(imp_data['x_coord']).points
    y = cube.coord(imp_data['y_coord']).points

    x2d, y2d = np.meshgrid(x, y)

    lonlat = TARGET_CRS.transform_points(source_crs, x2d, y2d)
    lon = lonlat[..., 0]
    lat = lonlat[..., 1]

    lat_coord = iris.coords.AuxCoord(lat, long_name='latitude_platecarree', 
                                     units='degrees')
    lon_coord = iris.coords.AuxCoord(lon, long_name='longitude_platecarree', 
                                     units='degrees')

    cube.add_aux_coord(lat_coord, (1, 2))
    cube.add_aux_coord(lon_coord, (1, 2))

    return cube


def find_lowest_mslp_in_polygon(cube, polygon):

    # Get latitude and longitude coordinates and flatten them
    lat = cube.coord('latitude_platecarree').points
    lon = cube.coord('longitude_platecarree').points
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    # Get data and flatten it
    data_flat = cube.data.flatten()

    # Create a list of values that are within the polygon
    values_in_polygon = [
        data for data, x, y in zip(data_flat, lon_flat, lat_flat)
        if polygon.contains(Point(x, y))
    ]

    # Find minimum value in the polygon
    min_val = min(values_in_polygon)

    return min_val


def get_regions():

    shapes = shp.Reader('MO_ForecastQNHRegions.shp')

    regions = {}
    for shape_rec in shapes.shapeRecords():
        name = shape_rec.record[0]
        shape = shape_rec.shape
        coords = shape.points
        regions[name] = coords

    return regions



def within_extent(cell):
    lon = cell.coord('longitude_platecarree').points
    lat = cell.coord('latitude_platecarree').points
    return ((PLOT_EXTENT[0] <= lon <= PLOT_EXTENT[1]) 
            and (PLOT_EXTENT[2] <= lat <= PLOT_EXTENT[3]))



if __name__ == "__main__":
    main()
    print('Finished')