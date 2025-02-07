// Import datasets
var s2_2A = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED");
var Nations = ee.FeatureCollection("FAO/GAUL/2015/level0");

// Define the area of interest: Canary Islands (part of Spain)
var spain_area = Nations.filter(ee.Filter.eq('ADM0_NAME', 'Spain')).geometry();

// Define the region of interest as a polygon for the Canary Islands
var Ari_islands_triangle = ee.Geometry.Polygon([
  [[-16, 29.4], [-16, 27.5], [-18.5, 27.5], [-16, 29.4]]
]);

// Intersect with Spain to limit to the Canary Islands
var Ari_islands = Ari_islands_triangle.intersection(spain_area);

// Add the region to the map
Map.addLayer(Ari_islands, {color: 'black'}, 'Ariagona Islands');
Map.centerObject(Ari_islands, 8);

// Cloud masking function for Sentinel-2 images
function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask)
    .divide(10000)
    .copyProperties(image, ['system:time_start']);
}

// Prepare the Sentinel-2 ImageCollection with cloud masking and date filtering
var s2_data = s2_2A
  .filterBounds(Ari_islands)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40))
  .map(maskS2clouds)
  .map(function(image) {
    return image.clip(Ari_islands);
  })
  .filter(ee.Filter.date('2015-01-01', '2024-12-31')) // Use the full time range
  .sort('system:time_start');

// Add indices to the images
function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndmi = image.normalizedDifference(['B8', 'B11']).rename('NDMI');
  var evi = image.expression(
    '(2.5 * ((NIR - RED)) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('B8'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
    }).rename('EVI');
  var savi = image.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5', {
      'NIR': image.select('B8'),
      'RED': image.select('B4')
    }).rename('SAVI');
  var tcVeg = image.expression(
    '-0.2848 * BLUE - 0.2435 * GREEN - 0.5436 * RED + 0.7243 * NIR + 0.0840 * SWIR1 - 0.1800 * SWIR2', {
      'BLUE': image.select('B2'),
      'GREEN': image.select('B3'),
      'RED': image.select('B4'),
      'NIR': image.select('B8'),
      'SWIR1': image.select('B11'),
      'SWIR2': image.select('B12')
    }).rename('TC_Veg');
  var tcWet = image.expression(
    '0.1509 * BLUE + 0.1973 * GREEN + 0.3279 * RED + 0.3406 * NIR - 0.7112 * SWIR1 - 0.4572 * SWIR2', {
      'BLUE': image.select('B2'),
      'GREEN': image.select('B3'),
      'RED': image.select('B4'),
      'NIR': image.select('B8'),
      'SWIR1': image.select('B11'),
      'SWIR2': image.select('B12')
    }).rename('TC_Wet');
  return image.addBands([ndvi, ndmi, evi, savi, tcVeg, tcWet]);
}

// Add indices to the ImageCollection
var survey = s2_data.map(addIndices);

// Calculate the mean of each index
var ndviMean = survey.select('NDVI').mean().reproject({crs: 'EPSG:32628', scale: 30});
var ndmiMean = survey.select('NDMI').mean().reproject({crs: 'EPSG:32628', scale: 30});
var eviMean = survey.select('EVI').mean().reproject({crs: 'EPSG:32628', scale: 30});
var saviMean = survey.select('SAVI').mean().reproject({crs: 'EPSG:32628', scale: 30});
var tcVegMean = survey.select('TC_Veg').mean().reproject({crs: 'EPSG:32628', scale: 30});
var tcWetMean = survey.select('TC_Wet').mean().reproject({crs: 'EPSG:32628', scale: 30});

// Export each index as a separate image
Export.image.toDrive({
  image: ndviMean,
  description: 'NDVI_Mean_2015_2024',
  folder: 'Environmental_Variables',
  region: Ari_islands,
  scale: 10,
  crs: 'EPSG:32628',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: ndmiMean,
  description: 'NDMI_Mean_2015_2024',
  folder: 'Environmental_Variables',
  region: Ari_islands,
  scale: 10,
  crs: 'EPSG:32628',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: eviMean,
  description: 'EVI_Mean_2015_2024',
  folder: 'Environmental_Variables',
  region: Ari_islands,
  scale: 10,
  crs: 'EPSG:32628',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: saviMean,
  description: 'SAVI_Mean_2015_2024',
  folder: 'Environmental_Variables',
  region: Ari_islands,
  scale: 10,
  crs: 'EPSG:32628',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: tcVegMean,
  description: 'TC_Veg_Mean_2015_2024',
  folder: 'Environmental_Variables',
  region: Ari_islands,
  scale: 10,
  crs: 'EPSG:32628',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: tcWetMean,
  description: 'TC_Wet_Mean_2015_2024',
  folder: 'Environmental_Variables',
  region: Ari_islands,
  scale: 10,
  crs: 'EPSG:32628',
  maxPixels: 1e13
});

// Add mean layers to the map for visualization
Map.addLayer(ndviMean, {min: 0, max: 1, palette: ['blue', 'green', 'yellow', 'red']}, 'NDVI Mean');
Map.addLayer(ndmiMean, {min: 0, max: 1, palette: ['blue', 'green', 'yellow', 'red']}, 'NDMI Mean');
Map.addLayer(eviMean, {min: 0, max: 1, palette: ['blue', 'green', 'yellow', 'red']}, 'EVI Mean');
Map.addLayer(saviMean, {min: 0, max: 1, palette: ['blue', 'green', 'yellow', 'red']}, 'SAVI Mean');
Map.addLayer(tcVegMean, {min: 0, max: 1, palette: ['blue', 'green', 'yellow', 'red']}, 'TC_Veg Mean');
Map.addLayer(tcWetMean, {min: 0, max: 1, palette: ['blue', 'green', 'yellow', 'red']}, 'TC_Wet Mean');
