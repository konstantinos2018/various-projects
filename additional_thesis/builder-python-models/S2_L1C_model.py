# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# S2_L1C_model.py
# Created on: 2020-10-30 17:08:18.00000
#   (generated by ArcGIS/ModelBuilder)
# Description: 
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy


# Local variables:
T32TQQ_20170620T100031_B08_jp2 = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\Sentinel 2\\S2A_MSIL1C_20170620T100031_N0205_R122_T32TQQ_20170620T100453.SAFE\\GRANULE\\L1C_T32TQQ_A010415_20170620T100453\\IMG_DATA\\T32TQQ_20170620T100031_B08.jp2"
rectangle_shp = "F:\\ArcGIS\\Additional_Thesis\\AreaOfInterest\\rectangle.shp"
s2l1cb08qq_clip_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2l1cb08qq_clip.tif"
T32TQQ_20170620T100031_B04_jp2 = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\Sentinel 2\\S2A_MSIL1C_20170620T100031_N0205_R122_T32TQQ_20170620T100453.SAFE\\GRANULE\\L1C_T32TQQ_A010415_20170620T100453\\IMG_DATA\\T32TQQ_20170620T100031_B04.jp2"
s2l1cb04qq_clip_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2l1cb04qq_clip.tif"
s2l1cndviqq_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2l1cndviqq.tif"
s2l1cndviqp__2_ = s2l1cndviqq_tif
T32TQP_20170620T100031_B08_jp2 = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\Sentinel 2\\S2A_MSIL1C_20170620T100031_N0205_R122_T32TQP_20170620T100453.SAFE\\GRANULE\\L1C_T32TQP_A010415_20170620T100453\\IMG_DATA\\T32TQP_20170620T100031_B08.jp2"
s2l1cb08qp_clip_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2l1cb08qp_clip.tif"
T32TQP_20170620T100031_B04_jp2 = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\Sentinel 2\\S2A_MSIL1C_20170620T100031_N0205_R122_T32TQP_20170620T100453.SAFE\\GRANULE\\L1C_T32TQP_A010415_20170620T100453\\IMG_DATA\\T32TQP_20170620T100031_B04.jp2"
s2l1cb04qp_clip_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2l1cb04qp_clip.tif"
s2l1cndviqp_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2l1cndviqp.tif"
LC08_L1TP_191029_20170620_20170630_01_T1_B5_TIF = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\Landsat 8\\LC08_L1TP_191029_20170620_20170630_01_T1_B5.TIF"
l8l1cb5_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\Landsat 8\\l8l1cb5.tif"
LC08_L1TP_191029_20170620_20170630_01_T1_B4_TIF = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\Landsat 8\\LC08_L1TP_191029_20170620_20170630_01_T1_B4.TIF"
l8l1cb4_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\Landsat 8\\l8l1cb4.tif"
l8l1cndvi_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\l8l1cndvi.tif"
rectangle_shp__2_ = "F:\\ArcGIS\\Additional_Thesis\\AreaOfInterest\\rectangle.shp"
l8l1cndviclip_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\l8l1cndviclip.tif"
l8l1cndviclip_tif__2_ = l8l1cndviclip_tif
s2lc1_ndvi20m_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2lc1_ndvi20m.tif"
s2lc1_ndvi20mclip_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2lc1_ndvi20mclip.tif"
Output_Cell_Size = "30 30"
Resampling_Technique = "NEAREST"
s2lc1_ndvi30m_tif = "F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2lc1_ndvi30m.tif"
Output_Link_File = ""

# Process: Clip (4)
arcpy.Clip_management(T32TQQ_20170620T100031_B08_jp2, "723426,097038214 4863831,6003664 806997,156311637 4937300,66346392", s2l1cb08qq_clip_tif, rectangle_shp, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# Process: Clip (3)
arcpy.Clip_management(T32TQQ_20170620T100031_B04_jp2, "723426,097038214 4863831,6003664 806997,156311637 4937300,66346392", s2l1cb04qq_clip_tif, rectangle_shp, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# Process: Raster Calculator (2)
arcpy.gp.RasterCalculator_sa("(Float(\"%s2l1cb08qq_clip.tif%\") - Float(\"%s2l1cb04qq_clip.tif%\")) / (Float(\"%s2l1cb08qq_clip.tif%\") + Float(\"%s2l1cb04qq_clip.tif%\"))", s2l1cndviqq_tif)

# Process: Clip (2)
arcpy.Clip_management(T32TQP_20170620T100031_B08_jp2, "723426,097038214 4863831,6003664 806997,156311637 4937300,66346392", s2l1cb08qp_clip_tif, rectangle_shp, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# Process: Clip
arcpy.Clip_management(T32TQP_20170620T100031_B04_jp2, "723426,097038214 4863831,6003664 806997,156311637 4937300,66346392", s2l1cb04qp_clip_tif, rectangle_shp, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# Process: Raster Calculator
arcpy.gp.RasterCalculator_sa("(Float(\"%s2l1cb08qp_clip.tif%\") - Float(\"%s2l1cb04qp_clip.tif%\")) / (Float(\"%s2l1cb08qp_clip.tif%\") + Float(\"%s2l1cb04qp_clip.tif%\"))", s2l1cndviqp_tif)

# Process: Mosaic
arcpy.Mosaic_management("'F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2l1cndviqq.tif';'F:\\ArcGIS\\Additional_Thesis\\Level 1C\\s2l1cndviqp.tif'", s2l1cndviqp_tif, "LAST", "FIRST", "", "", "NONE", "0", "NONE")

# Process: Raster Calculator (4)
arcpy.gp.RasterCalculator_sa("SetNull(\"%LC08_L1TP_191029_20170620_20170630_01_T1_B5.TIF%\" == 0, \"%LC08_L1TP_191029_20170620_20170630_01_T1_B5.TIF%\")", l8l1cb5_tif)

# Process: Raster Calculator (3)
arcpy.gp.RasterCalculator_sa("SetNull(\"%LC08_L1TP_191029_20170620_20170630_01_T1_B4.TIF%\" == 0, \"%LC08_L1TP_191029_20170620_20170630_01_T1_B4.TIF%\")", l8l1cb4_tif)

# Process: Raster Calculator (5)
arcpy.gp.RasterCalculator_sa("(Float(\"%l8l1cb5.tif%\") - Float(\"%l8l1cb4.tif%\")) / (Float(\"%l8l1cb5.tif%\") + Float(\"%l8l1cb4.tif%\"))", l8l1cndvi_tif)

# Process: Clip (5)
arcpy.Clip_management(l8l1cndvi_tif, "241514,694364513 4859043,41602539 330198,804643095 4938392,49105477", l8l1cndviclip_tif, rectangle_shp__2_, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# Process: Clip (6)
arcpy.Clip_management(s2lc1_ndvi20m_tif, "723426,097038214 4863831,6003664 806997,156311637 4937300,66346392", s2lc1_ndvi20mclip_tif, rectangle_shp, "", "NONE", "NO_MAINTAIN_EXTENT")

# Process: Resample
arcpy.Resample_management(s2lc1_ndvi20mclip_tif, s2lc1_ndvi30m_tif, Output_Cell_Size, Resampling_Technique)

# Process: Register Raster
arcpy.RegisterRaster_management(l8l1cndviclip_tif, "REGISTER", s2lc1_ndvi30m_tif, "", "POLYORDER1", Output_Link_File, "")
