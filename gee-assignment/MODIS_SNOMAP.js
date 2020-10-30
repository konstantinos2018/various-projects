var MOD = ee.ImageCollection('MODIS/006/MOD10A1');

MOD = MOD
  .filterBounds(Pyr)
  .filterDate('2010-10-01', '2011-04-01');
print(MOD)
var median = MOD.select('NDSI').median();
print('MOD', MOD);
print('median', median);
Map.centerObject(Pyr);  
Map.addLayer(median.clip(Pyr), {min:-100, max: 250}, 'Modis snow');
  
var MOD2 = ee.ImageCollection('MODIS/MYD13A1');
print(MOD2);