from arcpy import env

# Environment Variables
env.workspace = 'F:/ArcGIS/Additional_Thesis/Level 1C'

# List raster names into a list
rasters = arcpy.ListRasters()

# delete some rasters probably
rasters.pop(2)

for item in rasters:

	arcpy.Delete_management(item) # delete rasters
