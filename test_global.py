import iris
import cartopy.feature as feature
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import shapefile as shp

import iris.plot as iplt

# Define data directory
DATA_DIR = '/data/users/andre.lanyon/tafs/grid'


def main():

    # Get FOQNH region coordinates
    regions = get_regions()

    # dt_str = '20250415T1900Z-B20250414T1230Z'
    dt_str = '20250415T1000Z-PT0022H00M'
    fname = f'percentile_extract_{dt_str}-pressure_at_mean_sea_level.nc'
    cube = iris.load_cube(f'{DATA_DIR}/{fname}')

    # Convert units to hPa
    cube.convert_units('hPa')

    # Define the source CRS as a geographic coordinate system with a 
    # spherical globe
    source_crs = ccrs.Geodetic(globe=ccrs.Globe(semimajor_axis=6371229.0, 
                                                semiminor_axis=6371229.0))


    target_crs = ccrs.PlateCarree()

    x = cube.coord('longitude').points
    y = cube.coord('latitude').points
    x2d, y2d = np.meshgrid(x, y)

    lonlat = target_crs.transform_points(source_crs, x2d, y2d)
    lon = lonlat[..., 0]
    lat = lonlat[..., 1]

    lat_coord = iris.coords.AuxCoord(lat, long_name='latitude_platecarree', 
                                     units='degrees')
    lon_coord = iris.coords.AuxCoord(lon, long_name='longitude_platecarree', 
                                     units='degrees')

    cube.add_aux_coord(lat_coord, (1, 2))
    cube.add_aux_coord(lon_coord, (1, 2))

    # Get cube just at 50th percentile
    cube_50th = cube.extract(iris.Constraint(percentile=50))

    # Plot data on map using latitude and longitude coordinates
    map_crs = ccrs.TransverseMercator()
    fig, ax = plt.subplots(figsize=(20, 20), 
                           subplot_kw={'projection': map_crs})

    iplt.pcolormesh(cube_50th, cmap='jet')

    # Add colorbar
    cbar = plt.colorbar(ax=ax, orientation='horizontal', pad=0.05, 
                        aspect=50, shrink=0.8)
    cbar.set_label('Pressure at Mean Sea Level (hPa)', fontsize=14)

    # Draw coastlines and background stuff
    ax.add_feature(feature.COASTLINE,linewidth=0.8)
    ax.add_feature(feature.BORDERS, linestyle=':',linewidth=0.4)
    ax.set_extent([-10, 3, 48.5, 63.5], crs=target_crs)
    ax.gridlines(draw_labels=True, linewidth=0.4)

    # Plot FOQNH areas as polygons
    for area_name, coords in regions.items():

        # Create a polygon from the coordinates
        polygon = Polygon(coords)

        # Define the centroid and exterior of the polygon
        centroid = polygon.centroid
        x, y = polygon.exterior.xy

        # Find the lowest MSLP in the polygon
        lowest_mslp = find_lowest_mslp_in_polygon(cube_50th, polygon)

        # Plot the polygon and annotate with area name and lowest MSLP
        ax.plot(x, y, color='red', linewidth=1, transform=target_crs)
        txt = f'{area_name}\n{lowest_mslp:.1f} hPa'
        ax.text(centroid.x, centroid.y, txt, 
            fontsize=12, ha='center', va='center', color='white',
            fontweight='bold',
            transform=ccrs.PlateCarree())

    # Save and close plot
    fig.savefig('cube_50th_plot_global.png', bbox_inches='tight')
    plt.close()



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


if __name__ == "__main__":
    main()
    print('Finished')