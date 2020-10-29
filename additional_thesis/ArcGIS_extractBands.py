# Imports
from arcpy import env

# Environment Variables
env.workspace = 'F:/ArcGIS/Additional_Thesis'

# Create dictionary with names of bands and corresponding layer number of the image stack
dic = {'B4': 1, 'B5': 2, 'B6': 3, 'B7': 4, 'B8a': 5}

# Make Raster Layers (temporary) from each band
for key, value in dic.items():
	arcpy.management.MakeRasterLayer('F:/ArcGIS/Additional_Thesis/L2A/Sentinel 2/Mosaic/S220mstackClip.tif', 'S2_L1C_'+key+'.tif', '#', '#', str(value))

# Convert each band to Raster object and save it as GeoTIFF (make it permanent)
for key in dic.keys():	
	outraster = arcpy.Raster('S2_L1C_'+key+'.tif')
	outraster.save('F:/ArcGIS/Additional_Thesis/L2A/Sentinel 2/Mosaic/S2_L1C_'+key+'.tif')