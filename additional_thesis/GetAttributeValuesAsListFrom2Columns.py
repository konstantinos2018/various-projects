# Imports
import arcpy

fc = 'F:/ArcGIS/Additional_Thesis/Geometries/griddedpoints_l8_maize_30m.shp'
field1 = 'S2A_NDVI_A'
field2 = 'LC08_NDVI_'

# Do things on point shapefile
cursor = arcpy.SearchCursor(fc, fields= field1 + ';' + field2)

# Initialize
S2 = []
L8 = []

# open txt file to write the outputs
f = open('F:/ArcGIS/Additional_Thesis/Level 1C/AROP_iCOR_NDVI_30m/NDVI_iCOR_AROP_30m.txt', 'w')

# parse every row and write it to a list (for S2 and L8)
for i, row in enumerate(cursor):
    S2.append(row.getValue(field1))
    L8.append(row.getValue(field2))
	
	# Write outputs
    f.write('{0},{1}\n'.format(S2[i], L8[i]))

# Close txt file
f.close()