# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Landsat_8_Final.py
# Created on: 2020-10-30 17:06:03.00000
#   (generated by ArcGIS/ModelBuilder)
# Description: 
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy


# Local variables:
L8L2A_B5clip1_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\L8L2A_B5clip1.tif"
Output_Coordinate_System = "PROJCS['WGS_1984_UTM_Zone_32N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',9.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]"
Output_Cell_Size = "20 20"
L8L2A_B5clip1_ProjectRaster_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\ReprojectedResampled\\L8L2A_B5clip1_ProjectRaster.tif"
L8L2A_B4clip1_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\L8L2A_B4clip1.tif"
L8L2A_B4clip1_ProjectRaster_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\ReprojectedResampled\\L8L2A_B4clip1_ProjectRaster.tif"
L8ndviReproj_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\ReprojectedResampled\\NDVI\\L8ndviReproj.tif"
Maize_within_Eo_shp = "F:\\ArcGIS\\Additional_Thesis\\Geometries\\Maize_within_Eo.shp"
L8ndviReproj_Maize_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\ReprojectedResampled\\NDVI\\L8ndviReproj_Maize.tif"

# Set Geoprocessing environments
arcpy.env.scratchWorkspace = "F:\\ArcGIS\\Additional_Thesis\\L2A"
arcpy.env.snapRaster = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Sentinel 2\\Mosaic\\QPB04Clip2.tif"
arcpy.env.workspace = "F:\\ArcGIS\\Additional_Thesis\\L2A"

# Process: Project Raster
arcpy.ProjectRaster_management(L8L2A_B5clip1_tif, L8L2A_B5clip1_ProjectRaster_tif, Output_Coordinate_System, "NEAREST", Output_Cell_Size, "", "", "PROJCS['WGS_1984_UTM_Zone_33N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',15.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]")

# Process: Project Raster (2)
arcpy.ProjectRaster_management(L8L2A_B4clip1_tif, L8L2A_B4clip1_ProjectRaster_tif, Output_Coordinate_System, "NEAREST", Output_Cell_Size, "", "", "PROJCS['WGS_1984_UTM_Zone_33N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',15.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]")

# Process: Raster Calculator
arcpy.gp.RasterCalculator_sa("(Float(\"%L8L2A_B5clip1_ProjectRaster.tif%\") - Float(\"%L8L2A_B4clip1_ProjectRaster.tif%\"))/(Float(\"%L8L2A_B5clip1_ProjectRaster.tif%\") + Float(\"%L8L2A_B4clip1_ProjectRaster.tif%\"))", L8ndviReproj_tif)

# Process: Clip
arcpy.Clip_management(L8ndviReproj_tif, "731960 4868700 800660 4933666,77474976", L8ndviReproj_Maize_tif, Maize_within_Eo_shp, "-3.402823e+038", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

