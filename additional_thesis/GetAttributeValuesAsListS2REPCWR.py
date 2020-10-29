# Imports
import arcpy

from matplotlib import pylab as plt

# set environment variables
arcpy.env.workspace = 'F:/SNAP/iCOR/Landsat 8/CWR'

fc = 'F:/ArcGIS/Additional_Thesis/Geometries/griddedpoints_l8_maize_30m.shp'
field = ['RT_LC08_Br']

title = ''
for i, f in enumerate(field):
    if i == 0:
        break
    title = title + '{0};'.format(f)

title = title + '{0}'.format(field[0])

# Convert to list
titleList = title.split(';')

# Do things on point shapefile
cursor = arcpy.SearchCursor(fc, fields=title)

# open txt file to write the outputs
txtfile = open('F:/ArcGIS/Additional_Thesis/Miscellaneous/LC08_BrightTempB10.txt', 'w')

# parse every row and write it to a list (for S2 and L8)
for row in cursor:
    for j, f in enumerate(titleList):
        if j != 0:
            # Write outputs
            txtfile.write('{0},'.format(row.getValue(f)))
        else:
            txtfile.write('{0}\n'.format(row.getValue(f)))


# Close txt file
txtfile.close()



# # Plot
# plt.figure()
# plt.scatter(S2, L8, marker='.')
# plt.axis('tight')
# plt.xlabel('Sentinel 2')
# plt.ylabel('Landsat 8')
# plt.title('NDVI iCOR')
# plt.show()

# Get attribute fields as lists
# S2 = layer.GetField('S2A_Mosaic')
# L8 = layer.GetField('LC08_NDVI_')