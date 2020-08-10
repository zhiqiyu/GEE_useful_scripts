/**
* @file This is the original GEE code for the inland fishpond mapping in Singra 
*       Upazila of Bangladesh research.
* @author Zhiqi Yu <zyu2@gmu.edu>
* @license MIT
*/

// image assets imports
var s2 = ee.ImageCollection("COPERNICUS/S2"),
    bangladesh = ee.FeatureCollection("users/CSISS_LCLUC/bangladesh"),
    jrc = ee.ImageCollection("JRC/GSW1_1/YearlyHistory");


function getSentinelImagesFromTile(tile, bound, year, cloud){
  var pool = s2.filterMetadata('MGRS_TILE', 'equals', tile)
                  .filterBounds(bound)
                  .filterDate(String(year) + '-01-01', String(year+1) + '-01-01')
                  .filterMetadata('PROCESSING_BASELINE', 'not_ends_with', '2')
                  .filterMetadata('PRODUCT_ID', 'starts_with', 'S2A_OPER')
                  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', cloud);
  return pool;
}

function getIndexSeq(imgCollection, name, bands){
  // {ndvi: ['B8', 'B4'], mndwi: ['B3', 'B11']}
  var ndCollection = imgCollection.map(function(img){
    var nd = img.normalizedDifference(bands).rename(ee.String(img.get('system:index')).cat('_'+name));
    return nd;
  });
  return ndCollection.toBands();
}


function otsu(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });

  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
}


function awei_nsh(img, bands){
  return img.expression('(4*(g-swir1) - 0.25*nir - 2.75*swir2)/(g+nir+swir1+swir2)', {
    'g': img.select(bands['green']),
    'nir': img.select(bands['nir']),
    'swir1': img.select(bands['swir1']),
    'swir2': img.select(bands['swir2'])
  });
}


function add_spatial_attr(featureCollec, bound){
  // Add object-based spatial features to features
  var water_feature = featureCollec.map(function(f){
    var perimeter = f.perimeter(1);
    var area = f.area(1);
    var convex_hull = f.convexHull(1);
    var convex_area = convex_hull.area(1);
    var convex_perim = convex_hull.perimeter(1);
    var shaped = f.set('area', area)
            .set('perimeter', perimeter)
            .set('ipq', area.multiply(4*Math.PI).divide(perimeter.pow(2)))
            .set('soli', area.divide(convex_area))
            .set('pfd', perimeter.divide(4).log().multiply(2).divide(area.log()))
            .set('conv', convex_perim.divide(perimeter))
            .set('sqp', ee.Number(1).subtract(area.sqrt().multiply(4).divide(perimeter)))
    return shaped;
  });
  return water_feature;
}

function otsu_mask(image, mask, bins, scale, is_gt){
  var histogram = image.mask(mask).reduceRegion({
    reducer: ee.Reducer.histogram(bins, 0.001),
    geometry: image.geometry(), 
    scale: scale,
    maxPixels: 1e13
  });
  
  var masks = image.bandNames().map(function(name){
    var key = ee.String(name);
    var hist = ee.Dictionary(histogram.get(key));
    
    // compute threshold for the current band
    var thres = otsu(hist);
    
    // use threshold to segment the band
    var band = image.select(ee.String(name));
    if(is_gt){
      return band.gte(thres);
    }
    else{
      return band.lt(thres);
    }
  });
  
  var water_mask = ee.ImageCollection(masks).toBands().reduce(ee.Reducer.allNonZero());
  return water_mask;
}

//////////////////////////////////////////////////////////////////////

// ------------------------
// -- Data Preprocessing -- 
// ------------------------
// band names
var s2bands = {'green': 'B3', 'red': 'B4', 'nir': 'B8', 'swir1': 'B11', 'swir2': 'B12'};

// study area
var bound = ee.Feature(bangladesh.filterMetadata('ADM3_EN', 'equals', 'Singra').first());
var bound_geom = bound.geometry().geometries().get(2);
// var comilla = ee.Feature(bangladesh.filterMetadata('ADM2_EN', 'equals', 'Comilla').first());
// Map.addLayer(comilla)

var scale = 10;


// Sentinel2
var pool15 = getSentinelImagesFromTile('45RYH', bound.geometry(), 2016, 10);
print(pool15)

// index images
var ndvi15 = getIndexSeq(pool15, 'ndvi', [s2bands['nir'], s2bands['red']]).clip(bound);
var ndwi15 = getIndexSeq(pool15, 'ndwi', [s2bands['green'], s2bands['nir']]).clip(bound);
var mndwi15 = getIndexSeq(pool15, 'mndwi', [s2bands['green'], s2bands['swir1']]).clip(bound);
var awei15 = pool15.map(function(img){
  return awei_nsh(img, s2bands).rename(ee.String(img.get('GENERATION_TIME')).cat('_awei'))
}).toBands().clip(bound);

Export.image.toDrive({
  image: ndwi15,
  description: "ndwi15",
  folder: "EarthEngine",
  region: bound_geom,
  scale: scale,
  maxPixels: 1e13
});
Export.image.toDrive({
  image: mndwi15,
  description: "mndwi15",
  folder: "EarthEngine",
  region: bound_geom,
  scale: scale,
  maxPixels: 1e13
});
Export.image.toDrive({
  image: awei15,
  description: "awei15",
  folder: "EarthEngine",
  region: bound_geom,
  scale: scale,
  maxPixels: 1e13
});


// -------------------------------
// -- Pixel selection technique --
// -------------------------------
var naive_mask = mndwi15.gt(0).reduce(ee.Reducer.allNonZero());

Export.image.toDrive({
  image: naive_mask,
  description: "naive_mask",
  folder: "EarthEngine",
  region: bound_geom,
  scale: scale,
  maxPixels: 1e13
});

var pre_feature = naive_mask.selfMask().reduceToVectors({
  geometry: bound.geometry(),
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});

var pre_buffered = pre_feature.map(function(f){
  return f.buffer(5 * scale);  // was 5 * scale f.area(1).divide(Math.PI).sqrt().multiply(0.2)
});

var pre_mask = pre_buffered.reduceToImage({
  properties: ['count'],
  reducer: ee.Reducer.anyNonZero()
}).clip(bound);

Export.image.toDrive({
  image: pre_mask,
  description: "pre_mask",
  folder: "EarthEngine",
  region: bound_geom,
  scale: scale,
  maxPixels: 1e10
});

// ------------------------
// -- Spectral filtering --
// ------------------------
var ndwi_mask15 = otsu_mask(ndwi15, pre_mask, 500, scale, true);
var mndwi_mask15 = otsu_mask(mndwi15, pre_mask, 500, scale, true);
var awei_mask15 = otsu_mask(awei15, pre_mask, 500, scale, true);

Export.image.toDrive({
  image: ndwi15.mask(pre_mask),
  description: "ndwi_masked",
  folder: 'EarthEngine',
  region: bound_geom,
  scale: scale, 
  maxPixels: 1e13
});

Export.image.toDrive({
  image: mndwi15.mask(pre_mask),
  description: "mndwi_masked",
  folder: 'EarthEngine',
  region: bound_geom,
  scale: scale, 
  maxPixels: 1e13
});
Export.image.toDrive({
  image: awei15.mask(pre_mask),
  description: "awei_masked",
  folder: 'EarthEngine',
  region: bound_geom,
  scale: scale, 
  maxPixels: 1e13
});
              
var voted15 = ee.ImageCollection([ndwi_mask15, mndwi_mask15, awei_mask15])
              .toBands()
              .reduce(ee.Reducer.mode());

Export.image.toDrive({
  image: voted15,
  description: "voted15",
  folder: 'EarthEngine',
  region: bound_geom,
  scale: scale, 
  maxPixels: 1e13
});


/*******************
 * Spatial Filtering
 *******************/
var water_vector15 = voted15.selfMask().reduceToVectors({
  geometry: bound.geometry(),
  scale: scale,
  eightConnected: true,
  maxPixels: 1e13
});

// add spatial attributes to each feature
var water_shaped15 = add_spatial_attr(water_vector15, bound);
print('total water number: ', water_shaped15.size());


Export.table.toDrive({
  collection: water_shaped15,
  description: "water_shaped15",
  folder: "EarthEngine",
  fileFormat: 'SHP'
});


// logistic regresion
var pond_voted15 = water_shaped15.map(function(f) {
  var w0 = -10.136;
  var w1 = 4.459;
  var w2 = 2.666;
  var w3 = 2.035;
  var w4 = 4.277;
  var w5 = -3.027;
  var ipq = ee.Number(f.get('ipq')).multiply(w1)
  var soli = ee.Number(f.get('soli')).multiply(w2)
  var pfd = ee.Number(f.get('pfd')).multiply(w3)
  var conv = ee.Number(f.get('conv')).multiply(w4)
  var sqp = ee.Number(f.get('sqp')).multiply(w5)
  var weighted_sum = ipq.add(soli).add(pfd).add(conv).add(sqp).add(w0)
  var res = ee.Algorithms.If(weighted_sum.gte(0), 1, 0)
  return f.set('result', res);
}).filterMetadata('result', 'equals', 1);


// decision tree method 
// var pond_voted15 = water_shaped15.map(function(f){
//   var dec1 = ee.Number(f.get('conv')).gt(0.924);
//   var dec2 = ee.Number(f.get('conv')).gt(0.873);
//   var dec3 = ee.Number(f.get('ipq')).gt(0.503);
//   var dec4 = ee.Number(f.get('pfd')).gt(1.035);
  
//   var res = ee.Algorithms.If(dec1, 1, 
//               ee.Algorithms.If(dec2, ee.Algorithms.If(dec3, 1, 0), 0));
//   return f.set('result', res);
// }).filterMetadata('result', 'equals', 1);

// area threshold
var area_up_thres = 1000000;
var area_bot_thres = 200;

var pond15 = pond_voted15.filterMetadata('area', 'less_than', area_up_thres)
                                .filterMetadata('area', 'greater_than', area_bot_thres);

print('filtered water number: ', pond15.size());
print("total area of ponds: ", pond15.reduceColumns(ee.Reducer.sum(),
                                                    ['area']));

Export.table.toDrive({
  collection: pond15,
  description: "pond_result_5_dt",
  folder: "EarthEngine",
  fileFormat: 'SHP'
});

Export.table.toDrive({
  collection: pond15.map(function(f){return f.centroid(1)}),
  description: "pond_center_5_dt",
  folder: "EarthEngine",
  fileFormat: 'SHP'
});

// print(ui.Chart.feature.histogram(water_shaped15, 'ipq'))


// var jrc15 = jrc.filterMetadata('year', 'equals', 2016).clip(bangladesh.union().geometry()).toBands().gt(1).selfMask();
// jrc15 = jrc15.reduceToVectors({
//   geometry: bound.geometry(),
//   scale: scale,
//   eightConnected: true,
//   maxPixels: 1e13
// })

Map.setOptions('SATELLITE');
Map.centerObject(bound, 13);
Map.addLayer(ndvi15, {}, 'ndvi', false)
Map.addLayer(ndwi15, {}, 'ndwi', false)
Map.addLayer(mndwi15, {}, 'mdnwi', false)
// Map.addLayer(ndvi_mask00, {}, 'ndvi_mask', false);
// Map.addLayer(ndbi_mask00, {}, 'ndbi_mask', false);
// Map.addLayer(mndwi_mask00, {}, 'mndwi_mask', false);
// Map.addLayer(voted00, {}, 'voted00');
Map.addLayer(voted15, {}, 'voted15', false);
Map.addLayer(water_shaped15, {color: '#EEEEEE'}, 'water_feature15', false)
Map.addLayer(pond15, {color: '#EEEEEE'}, 'pond15')
// Map.addLayer(jrc15, {}, 'jrc', false);

