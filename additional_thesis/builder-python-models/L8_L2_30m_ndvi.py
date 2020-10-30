# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# L8_L2_30m_ndvi.py
# Created on: 2020-10-30 17:05:17.00000
#   (generated by ArcGIS/ModelBuilder)
# Description: 
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy


# Local variables:
L8L2A_B5clip1_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\L8L2A_B5clip1.tif"
L8L2A_B5clip1_Project_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\Reprojected\\L8L2A_B5clip1_Project.tif"
L8L2A_B4clip1_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\L8L2A_B4clip1.tif"
L8L2A_B4clip1_Project_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\Reprojected\\L8L2A_B4clip1_Project.tif"
L8_ndvi_Project_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\Reprojected\\L8_ndvi_Project.tif"
et0pm_polygDiss_shp = "F:\\ArcGIS\\Additional_Thesis\\Geometries\\et0pm_polygDiss.shp"
L8_ndvi_Project_Clip_tif = "F:\\ArcGIS\\Additional_Thesis\\L2A\\Landsat 8\\Clipped\\Reprojected\\NDVI\\L8_ndvi_Project_Clip.tif"

# Process: Project Raster (2)
arcpy.ProjectRaster_management(L8L2A_B5clip1_tif, L8L2A_B5clip1_Project_tif, "PROJCS['WGS_1984_UTM_Zone_32N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',9.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]", "NEAREST", "30 30", "", "", "PROJCS['WGS_1984_UTM_Zone_33N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',15.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]")

# Process: Project Raster
arcpy.ProjectRaster_management(L8L2A_B4clip1_tif, L8L2A_B4clip1_Project_tif, "PROJCS['WGS_1984_UTM_Zone_32N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',9.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]", "NEAREST", "30 30", "", "", "PROJCS['WGS_1984_UTM_Zone_33N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',15.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]")

# Process: Raster Calculator
arcpy.gp.RasterCalculator_sa("(Float(\"%L8L2A_B5clip1_Project.tif%\") - Float(\"%L8L2A_B4clip1_Project.tif%\"))/(Float(\"%L8L2A_B5clip1_Project.tif%\") + Float(\"%L8L2A_B4clip1_Project.tif%\"))", L8_ndvi_Project_tif)

# Process: Clip
arcpy.Clip_management(L8_ndvi_Project_tif, "731316,168484007 4867605,18488819 800766,168484008 4934565,18488819", L8_ndvi_Project_Clip_tif, et0pm_polygDiss_shp, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

