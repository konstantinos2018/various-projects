////////////// MASK CLOUDS //////////////////
// Import landsat images
var l8 = ee.ImageCollection('LANDSAT/LC8_L1T_TOA_FMASK');

// Filter imageCollection
var l8Masked = l8
  .select(['B4', 'B3', 'B2', 'B6'])
  .filterDate('2013-10-01', '2014-04-01')
  .filterBounds(Pyr)
  .filterMetadata('CLOUD_COVER', 'less_than', 5);
  
//Export.table.toDrive(l8Masked, 'Properties', 'Properties')
print(l8Masked);

// Apply count()
var countImage = l8Masked.count();
// Centerv view on region of interest
Map.centerObject(Pyr);
// Add the masked layer using RGB true color
Map.addLayer(l8Masked, {bands:['B4','B3','B2'], min:0, max:0.5}, 'True color')
// Add the masked layer using the count image
Map.addLayer(countImage, {}, 'Masked Clouds', false)
print(countImage)

///////////// NDSI //////////////////
// Create NDSI function
function addNdsi(image) {
  var ndsi = image.normalizedDifference(['B3', 'B6']).rename('ndsi'); // compute NDSI
  return image.addBands(ndsi); // ndsi band
}

var median = l8Masked.map(addNdsi).median();


Map.addLayer(median.select('ndsi').clip(Pyr), {palette: '000000, ffffff', min:-0.05, max:0.6}, 'NDSI', true);
print('median', median);


// Assign NDSIs to a different variable
var snow = l8Masked.map(addNdsi);
var snow2 = snow.select('ndsi');
print('snow2',snow2);
Map.addLayer(snow, {}, 'snow');

// Statistics for SUN_AZIMUTH
var sunStats = ee.Feature(l8Masked.aggregate_stats('SUN_ELEVATION'));
print('Sun azimuth statistics: ', sunStats);


