// =====================================================================
// Sardoba Reservoir, April-May 2020: Google Earth Engine scripts
// Author: Sunnat Abdukaxorov
//
// Each numbered section below is a SEPARATE Earth Engine script.
// Open code.earthengine.google.com, clear the editor, paste ONE section
// and press Run. Results appear in the Console tab on the right.
// Study area: 40.3465 N, 68.4719 E (box of 20 or 40 km, see each section).
// =====================================================================


// ---------------------------------------------------------------------
// 1. RADAR, SIMPLE THRESHOLD (20 km box, descending passes, -16 dB)
//    Gives the "simple method" numbers in the paper (2020 and 2019).
// ---------------------------------------------------------------------
var center = ee.Geometry.Point([68.4719, 40.3465]);
var aoi = center.buffer(10000).bounds();
Map.centerObject(aoi, 11);

function s1(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select('VV');
}

function areaKm2(img) {
  var sum = img.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: aoi, scale: 10, maxPixels: 1e10
  }).get('VV');
  return ee.Number(sum).divide(1e6);
}

var thr = -16;

function analyse(year) {
  var before = s1(year + '-04-01', year + '-05-01');
  var after  = s1(year + '-05-01', year + '-05-20');
  var b = before.median().focal_median(30, 'circle', 'meters').clip(aoi);
  var a = after.median().focal_median(30, 'circle', 'meters').clip(aoi);
  var flooded = a.lt(thr).and(b.lt(thr).not());
  var drained = b.lt(thr).and(a.lt(thr).not());
  print(year + ' images before', before.size());
  print(year + ' images after', after.size());
  print(year + ' new water km2', areaKm2(flooded));
  print(year + ' lost water km2', areaKm2(drained));
}
analyse(2020);
analyse(2019);


// ---------------------------------------------------------------------
// 2. RADAR, CHANGE DETECTION + PATCH FILTER (40 km box, orbit 151)
//    Gives the "filtered" numbers: new water and lost water, 2020 and 2019.
// ---------------------------------------------------------------------
var track = 151;   // relative orbit; 49 is the other track over the area

var center = ee.Geometry.Point([68.4719, 40.3465]);
var aoi = center.buffer(20000).bounds();
Map.centerObject(aoi, 10);

function s1(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .filter(ee.Filter.eq('relativeOrbitNumber_start', track))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select('VV');
}

function areaKm2(img) {
  var sum = img.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: aoi, scale: 20, maxPixels: 1e10
  }).get('VV');
  return ee.Number(sum).divide(1e6);
}

var thr = -16;     // water threshold (dB)
var drop = -3;     // minimum brightness change (dB)

function analyse(year) {
  var before = s1(year + '-04-01', year + '-05-01');
  var after  = s1(year + '-05-01', year + '-05-20');
  print(year + ' images before', before.size());
  print(year + ' images after', after.size());

  var b = before.median().focal_median(30, 'circle', 'meters').clip(aoi);
  var a = after.median().focal_median(30, 'circle', 'meters').clip(aoi);
  var diff = a.subtract(b);

  var wet = a.lt(thr).and(b.lt(thr).not()).and(diff.lt(drop)).selfMask();
  var wetClean = wet.updateMask(wet.connectedPixelCount(100, true).gte(30));
  var dry = b.lt(thr).and(a.lt(thr).not()).and(diff.gt(-drop)).selfMask();
  var dryClean = dry.updateMask(dry.connectedPixelCount(100, true).gte(30));

  print(year + ' new water (filtered) km2', areaKm2(wetClean));
  print(year + ' lost water (filtered) km2', areaKm2(dryClean));
}
analyse(2020);
analyse(2019);


// ---------------------------------------------------------------------
// 3. SENSITIVITY TO THE WATER THRESHOLDS (radar and optical)
//    Gives Tables 3a and 3b in the paper.
// ---------------------------------------------------------------------
var center = ee.Geometry.Point([68.4719, 40.3465]);
var aoi = center.buffer(20000).bounds();

function areaKm2(img, band) {
  var sum = img.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: aoi, scale: 20, maxPixels: 1e10
  }).get(band);
  return ee.Number(sum).divide(1e6);
}

// ---- Radar (orbit 151, descending) ----
function s1(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi).filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .filter(ee.Filter.eq('relativeOrbitNumber_start', 151))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select('VV');
}
function comp(year, start, end) {
  return s1(year + start, year + end).median()
    .focal_median(30, 'circle', 'meters').clip(aoi);
}
var b20 = comp('2020', '-04-01', '-05-01'), a20 = comp('2020', '-05-01', '-05-20');
var b19 = comp('2019', '-04-01', '-05-01'), a19 = comp('2019', '-05-01', '-05-20');

function radarLost(b, a, t) {
  var dry = b.lt(t).and(a.lt(t).not()).and(a.subtract(b).gt(3)).selfMask();
  var clean = dry.updateMask(dry.connectedPixelCount(100, true).gte(30));
  return areaKm2(clean, 'VV');
}
[-14, -15, -16, -17, -18].forEach(function(t) {
  print('Radar 2020, threshold ' + t + ' dB, lost km2', radarLost(b20, a20, t));
  print('Radar 2019, threshold ' + t + ' dB, lost km2', radarLost(b19, a19, t));
});

// ---- Optical (NDWI threshold) ----
function s2(start, end) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi).filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30)).median();
}
var nb = s2('2020-04-10', '2020-05-01').normalizedDifference(['B3', 'B8']);
var na = s2('2020-05-01', '2020-05-25').normalizedDifference(['B3', 'B8']);
[-0.1, 0, 0.1, 0.2].forEach(function(t) {
  print('Optical 2020, NDWI > ' + t + ', lost km2', areaKm2(nb.gt(t).and(na.gt(t).not()), 'nd'));
});


// ---------------------------------------------------------------------
// 4. OPTICAL (SENTINEL-2), LOST WATER AT NDWI > 0
//    Gives the 28.9 km2 optical result, the scene dates and the layers.
// ---------------------------------------------------------------------
var center = ee.Geometry.Point([68.4719, 40.3465]);
var aoi = center.buffer(20000).bounds();
Map.centerObject(aoi, 11);

function s2(start, end) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30));
}

var before = s2('2020-04-10', '2020-05-01');
var after  = s2('2020-05-01', '2020-05-25');

print('S2 before dates', before.aggregate_array('system:time_start').map(function(t){
  return ee.Date(t).format('YYYY-MM-dd');}).join(', '));
print('S2 after dates', after.aggregate_array('system:time_start').map(function(t){
  return ee.Date(t).format('YYYY-MM-dd');}).join(', '));

var rgb = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000};
Map.addLayer(before.median().clip(aoi), rgb, 'S2 before (true colour)');
Map.addLayer(after.median().clip(aoi), rgb, 'S2 after (true colour)');

var wb = before.median().normalizedDifference(['B3', 'B8']).gt(0);
var wa = after.median().normalizedDifference(['B3', 'B8']).gt(0);
Map.addLayer(wa.and(wb.not()).selfMask().clip(aoi), {palette: ['red']}, 'S2 new water', false);
Map.addLayer(wb.and(wa.not()).selfMask().clip(aoi), {palette: ['blue']}, 'S2 lost water', false);

function areaKm2(img) {
  var sum = img.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: aoi, scale: 20, maxPixels: 1e10
  }).get('nd');
  return ee.Number(sum).divide(1e6);
}
print('S2 lost water km2', areaKm2(wb.and(wa.not())));
print('S2 new water km2', areaKm2(wa.and(wb.not())));


// ---------------------------------------------------------------------
// 5. FIGURES: exports the four map figures to Google Drive
//    (folder "sardoba_figures"). Run, then start each task in the Tasks tab.
// ---------------------------------------------------------------------
var center = ee.Geometry.Point([68.4719, 40.3465]);
var aoi = center.buffer(20000).bounds();
Map.centerObject(aoi, 11);

function s2(start, end) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi).filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .median().clip(aoi);
}
var rgb = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000};
var s2Before2020 = s2('2020-04-10', '2020-05-01');
var s2After2020  = s2('2020-05-01', '2020-05-25');
var s2Before2019 = s2('2019-04-10', '2019-05-01');

var track = 151, thr = -16, drop = -3;
function s1(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi).filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .filter(ee.Filter.eq('relativeOrbitNumber_start', track))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select('VV');
}
function lostWater(year) {
  var b = s1(year + '-04-01', year + '-05-01').median()
    .focal_median(30, 'circle', 'meters').clip(aoi);
  var a = s1(year + '-05-01', year + '-05-20').median()
    .focal_median(30, 'circle', 'meters').clip(aoi);
  var dry = b.lt(thr).and(a.lt(thr).not()).and(a.subtract(b).gt(-drop)).selfMask();
  return dry.updateMask(dry.connectedPixelCount(100, true).gte(30));
}

var blue = {palette: ['0000FF']};
var fig1 = s2Before2020.visualize(rgb);
var fig2 = s2After2020.visualize(rgb);
var fig3 = s2Before2020.visualize(rgb).blend(lostWater(2020).visualize(blue));
var fig4 = s2Before2019.visualize(rgb).blend(lostWater(2019).visualize(blue));

Map.addLayer(fig1, {}, 'Fig 1: optical before 2020');
Map.addLayer(fig2, {}, 'Fig 2: optical after 2020', false);
Map.addLayer(fig3, {}, 'Fig 3: radar lost water 2020', false);
Map.addLayer(fig4, {}, 'Fig 4: radar lost water 2019 (control)', false);

function exportFig(img, name) {
  Export.image.toDrive({
    image: img, description: name, folder: 'sardoba_figures',
    region: aoi, scale: 20, maxPixels: 1e9, fileFormat: 'GeoTIFF'
  });
}
exportFig(fig1, 'fig1_optical_before_2020');
exportFig(fig2, 'fig2_optical_after_2020');
exportFig(fig3, 'fig3_radar_lost_water_2020');
exportFig(fig4, 'fig4_radar_lost_water_2019');
// Figure 5 (-14 dB vs -16 dB) is a screenshot of the map made with section 6.


// ---------------------------------------------------------------------
// 6. FIGURE 5 MAP (-14 dB in red, -16 dB in blue on top) AND SCENE DATES
// ---------------------------------------------------------------------
var center = ee.Geometry.Point([68.4719, 40.3465]);
var aoi = center.buffer(20000).bounds();
Map.centerObject(aoi, 11);

function s1(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi).filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .filter(ee.Filter.eq('relativeOrbitNumber_start', 151))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select('VV');
}
function comp(start, end) {
  return s1(start, end).median().focal_median(30, 'circle', 'meters').clip(aoi);
}
var b = comp('2020-04-01', '2020-05-01');
var a = comp('2020-05-01', '2020-05-20');
function lost(t) {
  var dry = b.lt(t).and(a.lt(t).not()).and(a.subtract(b).gt(3)).selfMask();
  return dry.updateMask(dry.connectedPixelCount(100, true).gte(30));
}
Map.addLayer(lost(-14), {palette: ['red']}, 'Lost water at -14 dB');
Map.addLayer(lost(-16), {palette: ['blue']}, 'Lost water at -16 dB');

function dates(start, end) {
  return s1(start, end).aggregate_array('system:time_start')
    .map(function(t) { return ee.Date(t).format('YYYY-MM-dd'); }).join(', ');
}
print('2020 before', dates('2020-04-01', '2020-05-01'));
print('2020 after',  dates('2020-05-01', '2020-05-20'));
print('2019 before', dates('2019-04-01', '2019-05-01'));
print('2019 after',  dates('2019-05-01', '2019-05-20'));
