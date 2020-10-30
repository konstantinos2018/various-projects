// COMPUTE SNOMAP 2015-2016 winter
var SentMasked = Sent2
  //.select(L8_NAMES, STD_NAMES)
  .filterDate('2015-10-01', '2016-04-01')
  .filterBounds(Pyr)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 9.5)
  //.filterMetadata('CLOUD_COVERAGE_ASSESSMENT', 'less_than', 4)
print('SentMasked', SentMasked);

// Apply maskcloud function on Landsat imageCollection
//var SentMasked = Sent2.map(maskClouds);

// Apply count()
var countImage = SentMasked.count();

Map.addLayer(SentMasked, {bands: ['B4', 'B3', 'B2'], min:-200, max:5000}, 'Sentinel Masked', true);
print('countImage', countImage);

// Create NDSI function
function addNdsi(image) {
  var ndsi = image.normalizedDifference(['B3', 'B11']).rename('ndsi')
  return image.addBands(ndsi);
}

var median = SentMasked.map(addNdsi).median();

var rgb = median.select(['B4', 'B3', 'B2']);
Map.addLayer(SentMasked, {bands: ['B4', 'B3', 'B2'], min:-200, max:5000}, 'RGB 2014-2015', true);
Map.centerObject(Pyr);
Map.addLayer(median.select('ndsi').clip(Pyr), {palette: '000000, ffffff', min:-0.05, max:0.25}, 'NDSI', true);
print('median', median);

// Statistics for SUN_AZIMUTH
var sunStats = SentMasked.aggregate_stats('MEAN_SOLAR_ZENITH_ANGLE');
print('Sun zenith statistics: ', sunStats);
