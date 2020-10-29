# Imports
import arcpy

from matplotlib import pylab as plt

# set environment variables
arcpy.env.workspace = 'F:/ArcGIS/Additional_Thesis/Silvia_MOSES_NDVI/Kostas_share_2/20_06_2017/S2/collocate/'
'comparisons txt'

fc = 'F:/ArcGIS/Additional_Thesis/Geometries/griddedpoints_l8_maize_30m.shp'
field1 = ['b2_Albedo_', 'b2_h_c_2_1', 'b2_MOSE_CW', 'b2_MOSE__1', 'b2_MOSES_K', 'b2_MOSES_1', 'b2_MOSES_L']
field2 = ['Albedo_171', 'h_c_2_171', 'MOSES_CWD_', 'MOSES_CWD1', 'MOSES_Kc_A', 'MOSES_Kc_N', 'MOSES_LAI_']

title = ''
for i, f in enumerate(zip(field1, field2)):
    if i == 6:
        break
    title = title + '{0};{1};'.format(f[0], f[1])

title = title + '{0};{1}'.format(field1[6], field2[6])

# Convert to list
titleList = title.split(';')

# Do things on point shapefile
cursor = arcpy.SearchCursor(fc, fields=title)

# open txt file to write the outputs
txtfile = open('F:/ArcGIS/Additional_Thesis/Silvia_MOSES_NDVI/Kostas_share_2/20_06_2017/S2/collocate/comparisons txt'
               '/S2L8_DOS1_Biophysical_collocate_30m_.txt', 'w')

# parse every row and write it to a list (for S2 and L8)
for i, row in enumerate(cursor):
    for j, f in enumerate(titleList):
        if j != 13:
            # Write outputs
            txtfile.write('{0},'.format(row.getValue(f)))
        else:
            print j
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