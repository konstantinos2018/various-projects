///////////////////////////////////////////////////////////////
// Set Region of interest
///////////////////////////////////////////////////////////////

//roi = geometry
Map.centerObject(roi)

// Visualize ROI
Map.addLayer(roi,{color:'red'},'ROI');

// Set dates
var start = '2015-07-01';
var end = '2017-03-30';


// Apply second function that accounts for masking
var maskCloud = function(image){
     return image.updateMask(image.select('Cloud').neq(1)).addBands(image.metadata('system:time_start').rename('time'))
}

// Function to QA bit data
var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    // Return the QA-bit image data
    return image.select([0], [newName]).bitwiseAnd(pattern).rightShift(start);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pre-process LANDSAT
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load collection
var Sent2 = ee.ImageCollection('COPERNICUS/S2');

// Filter collection based on location (roi) and date (2015-2016 winter)
var SentMasked = Sent2
  .filterBounds(roi)
  .filterDate(start,end)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 9.5)
print(SentMasked.getInfo())


 // Create function to mask out cloudy pixels
 var maskClouds = function(image) {
  // Select the QA band containing the cloud mask values
  var QA = image.select('QA60');
   // Return an image masking out cloudy areas (fmask=4) and shadow areas (fmask=2)
   var cloud = QA.neq(4).and(QA.neq(2)).rename('Cloud')
  return image.addBands(cloud).updateMask(QA.neq(1));
};

// Apply cloud mask function to filtered collection
var Sentmask = SentMasked.map(maskClouds);
print('Sentmask', Sentmask);
// Create function to convert RGB images to snow cover
var snomap = function(image){
  // Calculate Normalized Difference Snow Index
  var NDSI = image.normalizedDifference(['B3','B11']);
  // Convert NDSI to snow based on 0.4 threshold
  var snow = NDSI.gte(0.4);
  //var land = NDSI.lt(0.4);
  return image.addBands(snow.rename('Snow'));
}

// Apply the snomap function to the masked collection
var SentSnow = SentMasked.map(snomap);

// interpolate images ------------------------------------------------------
// code assisted (greatly) by Noel Gorelick
// finds each image's nearest neighbors, mosaics them, and then interpolates
// based on the image's time

// select bands of interest
var interpBands = ee.List(['time','Snow']);
var bandi = SentSnow.map(maskCloud).select(interpBands);

var time = 'system:time_start'
var lagRange = 10


// Looks for all images up to 'lagRange' days away from the current image.
var maxDiff = ee.Filter([
  ee.Filter.maxDifference(lagRange * (1000*60*60*24), time, null, time)])

// Images after, sorted in descending order (so closest is last).
var f1 = ee.Filter.and(maxDiff, ee.Filter.lessThanOrEquals(time, null, time))
var c1 = ee.Join.saveAll('after', time, false).apply(bandi, bandi, f1)

// Images before, sorted in ascending order (so closest is last).
var f2 = ee.Filter.and(maxDiff, ee.Filter.greaterThanOrEquals(time, null, time))
var c2 = ee.Join.saveAll('before', time, true).apply(c1, c1, f2)

var Sentinterpolated = ee.ImageCollection(c2.map(function(img) {
  img = ee.Image(img);
  var before = ee.ImageCollection.fromImages(ee.List(img.get('before'))).mosaic()
  var after = ee.ImageCollection.fromImages(ee.List(img.get('after'))).mosaic()

  // Compute the ratio between the image times.
  var x1 = before.select('time').double()
  var x2 = after.select('time').double()
  var now = ee.Image.constant(img.date().millis()).double();
  var ratio = now.subtract(x1).divide(x2.subtract(x1))  // this is zero anywhere x1 = x2
  // Compute the interpolated image.
  var interp = after.subtract(before).multiply(ratio).add(before)
  return interp.set('system:time_start', img.get('system:time_start'));
}))

// Landsat area
// TIME SERIES
print(ui.Chart.image.series(SentSnow.select(['Snow']), roi, ee.Reducer.sum(),1000)
  // .setChartType('ScatterChart')
  .setChartType('LineChart')
  .setOptions({
      title: 'Sentinel-2 snow area',
      hAxis: {'title': 'Time', min : new Date(start), max: new Date(end)},
      vAxis: {'title': 'Area [x1000 km]'},
      pointSize: 3,
      series: {
            0: { lineWidth: 1 },
            1: { lineWidth: 3 },
            2: { lineWidth: 3 }},
      crosshair: { trigger: 'both' },
      colors: ['#e41a1c','#4daf4a','#377eb8']
})
);

// Landsat fraction
print(ui.Chart.image.series(SentSnow.map(maskCloud).select(['Snow']), roi, ee.Reducer.mean(),1000)
  // .setChartType('ScatterChart')
  .setChartType('LineChart')
  //.setChartType('ColumnChart')
  .setOptions({
      title: 'Sentinel-2 area percentage',
      hAxis: {'title': 'Time', min : new Date(start), max: new Date(end)},
      vAxis: {'title': 'Area fraction [%]'},
      pointSize: 3,
      series: {
            0: { lineWidth: 3 },
            1: { lineWidth: 3 }},
      crosshair: { trigger: 'both' },
      colors:['#4daf4a','#377eb8']
}))

// Landsat 
Map.addLayer(Sentinterpolated.mean(),{bands:'Snow', min:0, max:1, palette:['#084594','#f7fbff']},'L8 snow frequency: interpolated') 
