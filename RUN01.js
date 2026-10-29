// ================================================================
// RUN01_V3.5 — QA MODES + PREPROCESSING + SEASONAL PCA + DAILY LAND/WATER TIME SERIES
// Spain saltmarsh unsupervised classification pipeline
//
// PURPOSE
//   Run ONE estuary and ONE season at a time. The script:
//   1) loads Sentinel-2 L2A SR images;
//   2) links COPERNICUS/S2_CLOUD_PROBABILITY correctly;
//   3) applies scale factors and resharpening (B5, B11, B12 → 10 m);
//   4) applies cloud/shadow, tree (>2 m), glint and water masks to EVERY
//      individual image; DPTM, Urban Atlas and Coastal Zones exclusions are
//      already embedded in the filtered estuary geometry loaded as roi;
//   5) calculates individual-image indices for the phenological time series:
//        ndvi, savi, gndvi, ndmi, grvi, valid;
//   6) exports each valid DATE as an image child of a pre-created
//      ImageCollection Asset;
//   7) builds one seasonal median composite, Gaussian-smoothed ONLY for
//      PCA/NDVI fallback products;
//   8) exports seasonal COMPOSITE, PCA3, NDVI, NDVI_CLASS and QA;\n//   9) displays one representative highest-coverage daily index image for QA.
//   9) maps only the first valid temporal image for index QA and stores
//      product metadata used by RUN02_V3 preflight validation.
//
// IMPORTANT ASSET PREPARATION
//   A) Create this regular Asset folder before running:
//      projects/ee-cassie-mario/assets/saltmarsh_v3_seasons
//
//   B) Create this as an IMAGE COLLECTION asset before running: 
//      projects/ee-cassie-mario/assets/saltmarsh_v3_timeseries
//
//      The script exports one Image asset per valid Sentinel-2 DATE below
//      that ImageCollection. These Assets preserve system:time_start, est_id,
//      season_tag and date_key so RUN02_V3 can load/filter/sort them directly.
//
//   C) The filtered estuary Asset used as roi must already exclude DPTM,
//      Urban Atlas and Coastal Zones areas that cannot be saltmarsh.
//      The tree-canopy mask remains an individual-image mask in GEE.
//
// IMPORTANT TREE-CANOPY ASSUMPTION
//   The supplied canopy Asset is a vector FeatureCollection of polygons with
//   tree canopy >2 m. All features are painted as exclusion areas.
//
// Mario Luiz Mascagni — RUN01_V3.4 (QA_ONLY / QA_AND_EXPORT + Otsu-buffer diagnostics)
// ================================================================


// ================================================================
// 1. USER CONTROLS
// ================================================================
  
// ---------------------------------------------------------------
// 1.1. Study estuary
// ---------------------------------------------------------------
var EST_ID_QUERY = 'marismas_de_joyel'; // 'marismas_de_joyel' 'tina_menor' 'ria_de_oyambre'

// Reference boundary used only for visual QA.
var ORIGINAL_ESTUARIES_ASSET =
  'projects/evapedraza-491311/assets/mitiga_saltmarsh/original_estuaries_spain_dissolv_est_id_1';

// Effective V3 processing area, prepared in QGIS:
// original estuary minus DPTM / Urban Atlas / Coastal Zones exclusions.

//'projects/evapedraza-491311/assets/mitiga_saltmarsh/estuaries_spain_copernicus_filtered_1'
//'projects/evapedraza-491311/assets/mitiga_saltmarsh/original_estuaries_spain_dissolv_est_id_1'
var ESTUARIES_ASSET =
  'projects/evapedraza-491311/assets/mitiga_saltmarsh/original_estuaries_spain_dissolv_est_id_1';

var originalEstuaries = ee.FeatureCollection(ORIGINAL_ESTUARIES_ASSET);
var originalEstuary = ee.Feature(
  originalEstuaries.filter(ee.Filter.eq('EST_ID', EST_ID_QUERY)).first()
);

var estuaries = ee.FeatureCollection(ESTUARIES_ASSET);
var estuary = ee.Feature(
  estuaries.filter(ee.Filter.eq('EST_ID', EST_ID_QUERY)).first()
);

// Effective processing geometry used throughout RUN01_V3.
var roi = estuary.geometry();
var originalRoi = originalEstuary.geometry();

Map.centerObject(roi, 13);
Map.addLayer(originalRoi, {color: '808080'}, 'Original estuary boundary', false);
Map.addLayer(roi, {color: 'ffff00'}, 'Filtered classification area', false);


// ---------------------------------------------------------------
// 1.2. Temporal controls — run one season at a time
// ---------------------------------------------------------------
var RUN01_MODE = 'QA_AND_EXPORT'; // 'QA_ONLY' | 'QA_AND_EXPORT'

var CLOUDY_PIXEL_PERCENTAGE_MAX = 30;

var OTSU_WATER_INDEX = 'MNDWI'; // 'MNDWI' | 'NDWI' | 'AWEI'
var OTSU_LAND_BUFFER = 0.06

var YEAR = 2024;
var ACTIVE_SEASON_INDEX = 2;  // 0 = SPR | 1 = SUM | 2 = AUT

// End dates are EXCLUSIVE in filterDate(). No overlap between seasons.
// AUT ends on 1 December so the phenology series covers March–November.

var SEASONS = [
  {tag: "SPR", start: "2024-03-01", end: "2024-06-01"},
  //{tag: "SPR", start: "2025-04-01", end: "2025-05-15"},
  {tag: "SUM", start: "2024-06-01", end: "2024-09-01"},
  //{tag: "SUM", start: "2025-06-22", end: "2025-09-24"},
  {tag: "AUT", start: "2024-09-01", end: "2024-12-01"}
  //{tag: "AUT", start: "2025-09-15", end: "2025-11-15"}
];

var SEASON = SEASONS[ACTIVE_SEASON_INDEX];
var SEASON_TAG = SEASON.tag;
var START_DATE = SEASON.start;
var END_DATE = SEASON.end;


// ---------------------------------------------------------------
// 1.3. Execution mode, QA and exports
// ---------------------------------------------------------------
// QA_ONLY      = calculate all preprocessing, Otsu diagnostics, seasonal
//                products and QA layers, but create NO export Tasks.
// QA_AND_EXPORT = create seasonal and time-series Tasks only after the
//                server-side preflight returns qa_pass = 1.
//var RUN01_MODE = 'QA_ONLY'; // 'QA_ONLY' | 'QA_AND_EXPORT'

// These controls are applied only when RUN01_MODE = 'QA_AND_EXPORT'.
var DO_EXPORT_ASSET = true;
var DO_EXPORT_DRIVE = false;

// Conservative structural checks. They deliberately do not prescribe the
// Otsu threshold or buffer value; those should be judged from QA maps.
var QA_MIN_CLOUD_MATCHED_IMAGES = 1;
var QA_MIN_VALID_DAILY_DATES = 1;
var QA_MIN_SEASONAL_VALID_PIXELS = 1;
var QA_REQUIRE_COMPOSITE_OTSU_WITHIN_CLAMP = true;

// Visual diagnostics specifically for tuning OTSU_LAND_BUFFER. They show
// pixels that are retained as land only because the positive buffer expands
// the Otsu land threshold.
var QA_ADD_OTSU_BUFFER_LAYERS = true;
var QA_PRINT_OTSU_BUFFER_SUMMARY = true;

// Leave blank to preserve the Sentinel-2 native projection of the selected
// estuary/tile. This avoids hard-coding EPSG:32630 for all Spanish estuaries.
// Set a CRS string only if you deliberately need a fixed export CRS.
var EXPORT_CRS = '';
var EXPORT_SCALE = 10;
var EXPORT_MAX_PIXELS = 1e13;

var DRIVE_FOLDER = 'GEE_saltmarsh_exports';

// This must be a normal Asset folder.
var SEASON_EXPORT_FOLDER =
  //'projects/evapedraza-491311/assets/mitiga_saltmarsh/run01_mitiga_saltmarsh';
  'projects/evapedraza-491311/assets/Joyel2024_run01';
// This must be pre-created as an ImageCollection Asset, not a normal folder.
// It preserves the original six-band index contract used by RUN02.
var TIME_SERIES_COLLECTION_ASSET =
  //'projects/evapedraza-491311/assets/mitiga_saltmarsh/saltmarsh_new_run01_timeseries';
  'projects/evapedraza-491311/assets/Joyel2024_vegindex_timeseries_run01';
// This must also be pre-created as an ImageCollection Asset. It is separate
// from the index collection so RUN02_V3 retains its existing band contract.
var LAND_WATER_TIME_SERIES_COLLECTION_ASSET =
  //'projects/evapedraza-491311/assets/mitiga_saltmarsh/saltmarsh_new_run01_landwater_timeseries';
  'projects/evapedraza-491311/assets/Joyel2024_run01_landwater'
// Individual-date exports can produce several Tasks per season.
// -1 = export every available date. Set a positive value only while testing.
var MAX_TS_EXPORTS_PER_RUN = -1;
var MAX_LAND_WATER_TS_EXPORTS_PER_RUN = -1;

// Optional seasonal RGB QA exports. The full COMPOSITE contains B2...B12
// already, so these are only convenience products for visual inspection.
var EXPORT_RGB_QA = true;


// ---------------------------------------------------------------
// 1.4. Debug / map display
// ---------------------------------------------------------------
var DEBUG = true;
var DEBUG_LEVEL = 1; // 0=min | 1=essential | 2=moderate | 3=heavy

var DO_ADD_MAPS = true;

// Adds ONE representative daily time-series image to the map for visual QA.
// The selected image is the valid daily product with the HIGHEST spatial
// coverage inside static_allowed_mask, not merely the earliest date with one
// surviving pixel. It never creates extra export Tasks.
var DO_ADD_FIRST_TS_INDEX_QC = true;
var FIRST_TS_QC_PRINT_STATS = true;


// ---------------------------------------------------------------
// 1.5. Sentinel-2 cloud / shadow / smoothing settings
// ---------------------------------------------------------------
// Scene metadata filter: intentionally permissive because actual masking is
// applied pixel-by-pixel using S2_CLOUD_PROBABILITY.
// var CLOUDY_PIXEL_PERCENTAGE_MAX = 60;

var PIX_CLOUD_PROB_THRESHOLD = 60;
var CLOUD_ERODE_PIXELS = 1;
var CLOUD_DILATE_PIXELS = 3;

var NIR_DARK_THRESHOLD = 0.20;
var CLOUD_SHADOW_DISTANCE_PIXELS = 40; // ≈400 m at 10 m pixels
var CLOUD_SHADOW_BUFFER_METERS = 100;

// Gaussian blur is used ONLY on the seasonal spectral composite for PCA and
// NDVI fallback products. Individual temporal index images are not blurred.
var DO_SEASONAL_GAUSSIAN_BLUR = true;
var GAUSSIAN_RADIUS_PIXELS = 1; // 3 × 3 neighborhood
var GAUSSIAN_SIGMA_PIXELS = 0.5;


// ---------------------------------------------------------------
// 1.6. Static environmental / administrative masks
// ---------------------------------------------------------------
// DPTM and Urban Atlas / Coastal Zones code is retained below only as
// historical documentation. It is disabled because those exclusions are
// already incorporated into ESTUARIES_ASSET in QGIS.

/*
// ---------------------------
// 1.6.1. DPTM allowed-area mask
// ---------------------------
// REQUIRED: replace with the final polygon/raster asset representing the
// area INSIDE the official DPTM administrative limit.
var USE_DPTM_MASK = true;
var DPTM_ALLOWED_ASSET =
  'REPLACE_WITH_YOUR_DPTM_ALLOWED_AREA_ASSET';

// 'FEATURE_COLLECTION' or 'IMAGE'
var DPTM_SOURCE_TYPE = 'FEATURE_COLLECTION';

// FEATURE_COLLECTION mode:
// true = all features represent allowed DPTM area.
// false = filter by DPTM_CLASS_FIELD / DPTM_CLASS_VALUES before painting.
var DPTM_USE_ALL_FEATURES = true;
var DPTM_CLASS_FIELD = '';
var DPTM_CLASS_VALUES = [];

// IMAGE mode only:
var DPTM_IMAGE_BAND = '';       // blank = first band
var DPTM_IMAGE_THRESHOLD = 0;   // pixels > threshold are allowed
*/

// ---------------------------
// 1.6.2. Tree canopy > 2 m mask
// ---------------------------
var USE_TREE_MASK = true;
var TREE_ASSET =
  'projects/evapedraza-491311/assets/mitiga_saltmarsh/spain_estuaries_canopy_tree_gt2m_1';

// Uploaded as a vector Shapefile/FeatureCollection of polygons representing
// tree canopy higher than 2 m. Every feature is an exclusion polygon.
var TREE_SOURCE_TYPE = 'FEATURE_COLLECTION';
var TREE_USE_ALL_FEATURES = true;

// Not used in FEATURE_COLLECTION mode. Keep blank for clarity.
var TREE_CLASS_FIELD = '';
var TREE_CLASS_VALUES = [];
var TREE_IMAGE_BAND = '';
var TREE_IMAGE_THRESHOLD = 0;

/*
// ---------------------------
// 1.6.3. Urban Atlas + coastal land-cover exclusion mask
// ---------------------------
var USE_URBAN_MASK = true;
var URBAN_ASSET =
  'projects/ee-ihset-mario/assets/coastal_zones_2021_atlas_2018_simplified';

// Expected source type for a simplified Urban Atlas / coastal LULC polygon.
var URBAN_SOURCE_TYPE = 'FEATURE_COLLECTION';

// Set true only when every feature in the Asset is an exclusion polygon.
var URBAN_USE_ALL_FEATURES = true;

// Use these only when URBAN_USE_ALL_FEATURES = false.
var URBAN_CLASS_FIELD = '';
var URBAN_CLASS_VALUES = [];

// IMAGE mode only (unused with the current default).
var URBAN_IMAGE_BAND = '';
var URBAN_IMAGE_THRESHOLD = 0;
*/

// ---------------------------------------------------------------
// 1.7. Conservative B8-based sunglint guard
// ---------------------------------------------------------------
// This is a WATER-ONLY mask, not a broad B8 threshold over the entire scene.
// It removes very bright NIR water candidates while preserving vegetation:
//   glint candidate = B8 > threshold AND MNDWI > threshold.
// Tune only after visual QA over representative estuaries.
var DO_GLINT_WATER_MASK = true;
var GLINT_B8_MIN = 0.15;
var GLINT_MNDWI_MIN = 0.15;


// ---------------------------------------------------------------
// 1.8. Water mask: selectable Otsu index by individual image
// ---------------------------------------------------------------
// Select the index used ONLY to derive the Otsu land/water threshold:
//   'MNDWI' = (B3 - B11) / (B3 + B11), default and legacy behaviour
//   'NDWI'  = (B3 - B8)  / (B3 + B8), useful where SWIR1 is unreliable
//   'AWEI'  = AWEIsh, a non-normalised water-enhancement index
//
// All three indices increase towards water. Consequently the rule remains:
//     land = selected_water_index <= threshold + OTSU_LAND_BUFFER
//
// The selection does NOT alter the RUN02 six-band time-series contract. The
// land/water ImageCollection retains `mndwi` as a diagnostic band for
// backward compatibility; its `water_threshold` band is documented through
// image properties as the threshold of the selected Otsu index.
//var OTSU_WATER_INDEX = 'MNDWI'; // 'MNDWI' | 'NDWI' | 'AWEI'

var OTSU_WATER_INDEX_OPTIONS = {
  MNDWI: {
    band: 'mndwi',
    label: 'MNDWI (B3-B11)/(B3+B11)',
    fallbackThreshold: 0.05,
    clampMin: -0.40,
    clampMax: 0.40
  },
  NDWI: {
    band: 'ndwi',
    label: 'NDWI (B3-B8)/(B3+B8)',
    fallbackThreshold: 0.05,
    clampMin: -0.40,
    clampMax: 0.40
  },
  AWEI: {
    band: 'awei',
    label: 'AWEIsh',
    fallbackThreshold: 0.00,
    clampMin: -1.50,
    clampMax: 1.50
  }
};

OTSU_WATER_INDEX = String(OTSU_WATER_INDEX).toUpperCase();
if (!OTSU_WATER_INDEX_OPTIONS.hasOwnProperty(OTSU_WATER_INDEX)) {
  throw new Error(
    'Invalid OTSU_WATER_INDEX: ' + OTSU_WATER_INDEX +
    '. Use MNDWI, NDWI or AWEI.'
  );
}

var ACTIVE_OTSU_WATER_INDEX = OTSU_WATER_INDEX_OPTIONS[OTSU_WATER_INDEX];
var OTSU_WATER_INDEX_BAND = ACTIVE_OTSU_WATER_INDEX.band;
var OTSU_WATER_INDEX_LABEL = ACTIVE_OTSU_WATER_INDEX.label;

// These active values are derived from the selected index configuration.
// They retain the V3.4 values when OTSU_WATER_INDEX = 'MNDWI'.
var WATER_FALLBACK_THRESHOLD = ACTIVE_OTSU_WATER_INDEX.fallbackThreshold;
var OTSU_CLAMP_MIN = ACTIVE_OTSU_WATER_INDEX.clampMin;
var OTSU_CLAMP_MAX = ACTIVE_OTSU_WATER_INDEX.clampMax;

var USE_OTSU_WATER_THRESHOLD = true;

// Positive buffer makes land retention more permissive because land is
// defined as selected_water_index <= threshold + buffer. The magnitude is in
// the units of the selected index. Values tuned for MNDWI/NDWI should not be
// transferred automatically to AWEI; inspect the QA layers when testing AWEI.
// var OTSU_LAND_BUFFER = 0.15;

var OTSU_SCALE = 10;
var OTSU_MAX_BUCKETS = 256;
var OTSU_MIN_BUCKET_WIDTH = 0.001;
var OTSU_TILE_SCALE = 4;
var OTSU_MIN_VALID_PIXELS = 200;

// Optional QA comparison layers. The selected index is always mapped when
// QA_ADD_OTSU_BUFFER_LAYERS=true; enable this only to compare alternatives.
var QA_ADD_ALL_WATER_INDEX_LAYERS = false;

// Apply a second seasonal-composite water mask before PCA and NDVI_class.
var APPLY_COMPOSITE_LAND_MASK = true;


// ---------------------------------------------------------------
// 1.9. Indices and PCA
// ---------------------------------------------------------------
var SPECTRAL_BANDS = ['B2','B3','B4','B5','B8','B11','B12'];
var TIME_SERIES_INDEX_BANDS = ['ndvi','savi','gndvi','ndmi','grvi'];
var PCA_BANDS = ['B2','B3','B4','B5','B8','B11','B12'];

// ----------------------------------------------------------------
// 1.9.1. Export contract for RUN02_V3 preflight QA
// ----------------------------------------------------------------
// These values do not change the raster pixels. They make every exported
// product self-describing so RUN02 can reject missing, outdated or
// incompatible Assets before the Fourier fitting step.
var PROCESSING_VERSION = 'RUN01_V3.5';
var SCHEMA_VERSION = 'saltmarsh_v3.1';

// Keep this label as "original_estuary_validation" while comparing V3 with
// the legacy NDVI-proxy workflow. Change it only when ESTUARIES_ASSET is
// switched to the filtered QGIS classification polygons.
var ROI_VERSION = 'original_estuary_validation';

var TS_BAND_SCHEMA = 'ndvi,savi,gndvi,ndmi,grvi,valid';

// Daily land-water collection for RUN04/RUN06. These bands are kept in a
// separate ImageCollection so RUN02 continues to load only the six index bands.
var LAND_WATER_TS_BANDS = [
  'analysis_observed',
  'land_observed',
  'water_observed',
  'mndwi',
  'water_threshold',
  'otsu_used',
  'valid'
];

var LAND_WATER_TS_BAND_SCHEMA =
  'analysis_observed,land_observed,water_observed,mndwi,' +
  'water_threshold,otsu_used,valid';

var COMPOSITE_BAND_SCHEMA = 'B2,B3,B4,B5,B8,B11,B12';
var PCA_BAND_SCHEMA = 'PC1,PC2,PC3';
var QA_BAND_SCHEMA =
  'obs_count,static_allowed_mask,land_mask_composite,' +
  'valid_season,n_source_images,n_dates_exported,' +
  'composite_water_threshold';

var EXPORT_CRS_LABEL = EXPORT_CRS === '' ? 'native_sentinel2' : EXPORT_CRS;


// ---------------------------------------------------------------
// 1.10. Visualisation
// ---------------------------------------------------------------
var RGB_BANDS = ['B4','B3','B2'];
var VIS_RGB = {bands: RGB_BANDS, min: 0.02, max: 0.25};
var VIS_NDVI = {
  min: -1, max: 1,
  palette: ['313695','74add1','e0f3f8','ffffff','fee090','f46d43','a50026']
};
var VIS_MNDWI = {
  min: OTSU_CLAMP_MIN,
  max: OTSU_CLAMP_MAX,
  palette: ['8c510a','d8b365','f6e8c3','c7eae5','5ab4ac','01665e']
};
var VIS_BINARY = {min: 0, max: 1, palette: ['ff00ff']};
var VIS_PCA = {
  min: -0.2, max: 0.2,
  palette: ['313695','4575b4','74add1','abd9e9','e0f3f8','ffffff',
            'fee090','fdae61','f46d43','d73027','a50026']
};
var NDVI_CLASS_PALETTE = [
  '000000', // 0 = invalid / masked
  '2c7bb6', // 1 = 0.0–0.2
  'abd9e9', // 2 = 0.2–0.4
  'ffffbf', // 3 = 0.4–0.6
  'fdae61', // 4 = 0.6–0.8
  'd7191c'  // 5 = >=0.8
];
var VIS_NDVI_CLASS = {min: 0, max: 5, palette: NDVI_CLASS_PALETTE};


// ================================================================
// 2. HELPERS
// ================================================================

function dbg(level) {
  if (!DEBUG || DEBUG_LEVEL < level) return;
  print.apply(null, Array.prototype.slice.call(arguments, 1));
}

function assetPath(name) {
  return SEASON_EXPORT_FOLDER + '/' + name + '_' + EST_ID_QUERY + '_' +
    SEASON_TAG + YEAR;
}

// Metadata contract shared by all seasonal products.
function baseProductMetadata(productType, bandSchema) {
  return {
    est_id: EST_ID_QUERY,
    season_tag: SEASON_TAG,
    year: YEAR,
    season_start: START_DATE,
    season_end_exclusive: END_DATE,

    product_type: productType,
    band_schema: bandSchema,
    processing_version: PROCESSING_VERSION,
    schema_version: SCHEMA_VERSION,

    roi_version: ROI_VERSION,
    source_roi_asset: ESTUARIES_ASSET,
    source_original_roi_asset: ORIGINAL_ESTUARIES_ASSET,
    source_tree_asset: TREE_ASSET,
    static_mask_mode: 'roi_plus_tree_gt2m',

    s2_collection: 'COPERNICUS/S2_SR_HARMONIZED',
    cloud_probability_collection: 'COPERNICUS/S2_CLOUD_PROBABILITY',
    otsu_water_index: OTSU_WATER_INDEX,
    otsu_water_band: OTSU_WATER_INDEX_BAND,
    otsu_water_index_label: OTSU_WATER_INDEX_LABEL,
    otsu_fallback_threshold: WATER_FALLBACK_THRESHOLD,
    otsu_land_buffer: OTSU_LAND_BUFFER,
    otsu_clamp_min: OTSU_CLAMP_MIN,
    otsu_clamp_max: OTSU_CLAMP_MAX,
    export_scale_m: EXPORT_SCALE,
    export_crs: EXPORT_CRS_LABEL
  };
}

function makeExportArgs(image, description, assetId, regionGeom, pyramidingPolicy) {
  var args = {
    image: ee.Image(image),
    description: description,
    assetId: assetId,
    region: regionGeom,
    scale: EXPORT_SCALE,
    maxPixels: EXPORT_MAX_PIXELS,
    pyramidingPolicy: pyramidingPolicy || {'.default': 'mean'}
  };

  if (EXPORT_CRS !== '') args.crs = EXPORT_CRS;
  return args;
}

function makeDriveExportArgs(image, description, regionGeom) {
  var args = {
    image: ee.Image(image),
    description: description,
    folder: DRIVE_FOLDER,
    region: regionGeom,
    scale: EXPORT_SCALE,
    maxPixels: EXPORT_MAX_PIXELS,
    fileFormat: 'GeoTIFF'
  };

  if (EXPORT_CRS !== '') args.crs = EXPORT_CRS;
  return args;
}

function printProjection(label, image) {
  print(label + ' projection:', ee.Image(image).projection());
}

function getFirstBand(image, preferredBandName) {
  image = ee.Image(image);
  return preferredBandName === '' ? image.select([0]) : image.select(preferredBandName);
}

function featureCollectionMask(assetId, useAllFeatures, classField, classValues, outName) {
  var fc = ee.FeatureCollection(assetId);

  if (!useAllFeatures) {
    fc = fc.filter(ee.Filter.inList(classField, classValues));
  }

  return ee.Image(0).byte()
    .paint(fc, 1)
    .rename(outName)
    .clip(roi)
    .unmask(0)
    .eq(1);
}

function binaryMaskFromConfiguredAsset(config) {
  var mask;

  if (config.sourceType === 'FEATURE_COLLECTION') {
    mask = featureCollectionMask(
      config.asset,
      config.useAllFeatures,
      config.classField,
      config.classValues,
      config.outName
    );
  } else {
    mask = getFirstBand(ee.Image(config.asset), config.imageBand)
      .unmask(0)
      .gt(config.imageThreshold)
      .rename(config.outName)
      .clip(roi);
  }

  return mask.unmask(0).toByte();
}

function getStaticMasks() {
  // DPTM, Urban Atlas and Coastal Zones were already applied in QGIS.
  // The filtered estuary geometry (roi) is therefore the base allowed area.
  var classificationArea = ee.Image.constant(1)
    .rename('classification_area_mask')
    .clip(roi)
    .unmask(0)
    .toByte();

  // Keep the canopy tree mask in GEE because it contains small, discontinuous
  // internal exclusions that should not fragment the estuary polygon in QGIS.
  var treeExcluded = ee.Image.constant(0)
    .rename('tree_excluded')
    .clip(roi)
    .unmask(0)
    .toByte();

  if (USE_TREE_MASK) {
    treeExcluded = binaryMaskFromConfiguredAsset({
      asset: TREE_ASSET,
      sourceType: TREE_SOURCE_TYPE,
      useAllFeatures: TREE_USE_ALL_FEATURES,
      classField: TREE_CLASS_FIELD,
      classValues: TREE_CLASS_VALUES,
      imageBand: TREE_IMAGE_BAND,
      imageThreshold: TREE_IMAGE_THRESHOLD,
      outName: 'tree_excluded'
    });
  }

  // The filtered classification polygon is enforced on every image, and tree
  // pixels are removed on top of it.
  var staticAllowed = classificationArea.eq(1)
    .and(treeExcluded.eq(0))
    .rename('static_allowed_mask')
    .unmask(0)
    .toByte();

  return {
    classificationArea: classificationArea,
    treeExcluded: treeExcluded,
    staticAllowed: staticAllowed
  };
}

function commonValidMask(image, bandNames) {
  return ee.Image(image)
    .select(bandNames)
    .mask()
    .reduce(ee.Reducer.min())
    .rename('valid_mask');
}

// Assigned after getStaticMasks() is built. It is the common denominator for
// daily spatial coverage, independent of clouds/water in a particular date.
var staticAllowedPixelCount = ee.Number(0);

function addImageValidCount(image) {
  image = ee.Image(image);
  var valid = image.select('valid');

  var count = ee.Number(valid.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 20,
    maxPixels: EXPORT_MAX_PIXELS,
    tileScale: 4
  }).get('valid', 0));

  var validFractionStatic = ee.Number(ee.Algorithms.If(
    staticAllowedPixelCount.gt(0),
    count.divide(staticAllowedPixelCount),
    0
  ));

  return image.set({
    valid_pixel_count: count,
    valid_fraction_static: validFractionStatic
  });
}

function addLayerIf(condition, image, vis, name, shown) {
  if (condition) Map.addLayer(image, vis, name, shown);
}


// ================================================================
// 3. SENTINEL-2 PREPROCESSING FUNCTIONS
// ================================================================

// ---------------------------------------------------------------
// 3.1. Scale Sentinel-2 SR reflectance to 0–1
// ---------------------------------------------------------------
function applyScaleFactors(image) {
  var optical = image.select('B.*').multiply(0.0001);
  return image.addBands(optical, null, true);
}


// ---------------------------------------------------------------
// 3.2. Sharpen B5, B11 and B12 with B8
// ---------------------------------------------------------------
var geeSharp = require('users/aazuspan/geeSharp:geeSharp');

function sharpened(image) {
  image = ee.Image(image);

  var sharpB5 = geeSharp.sharpen(image.select(['B5']), image.select(['B8']))
    .rename('B5_sharp');
  var sharpB11 = geeSharp.sharpen(image.select(['B11']), image.select(['B8']))
    .rename('B11_sharp');
  var sharpB12 = geeSharp.sharpen(image.select(['B12']), image.select(['B8']))
    .rename('B12_sharp');

  return image.addBands([sharpB5, sharpB11, sharpB12]);
}

function selectSharpenedSpectralBands(image) {
  return ee.Image(image).select(
    ['B2','B3','B4','B5_sharp','B8','B11_sharp','B12_sharp'],
    ['B2','B3','B4','B5','B8','B11','B12']
  ).copyProperties(image, image.propertyNames());
}


// ---------------------------------------------------------------
// 3.3. Link and apply S2 cloud probability / cloud shadows
// ---------------------------------------------------------------
function linkCloudProbability(collection) {
  var cloudProbability = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
    .filterBounds(roi)
    .filterDate(START_DATE, END_DATE);

  var joined = ee.Join.saveFirst('cloud_probability_image').apply({
    primary: ee.ImageCollection(collection),
    secondary: cloudProbability,
    condition: ee.Filter.equals({
      leftField: 'system:index',
      rightField: 'system:index'
    })
  });

  // Images without a matching cloud-probability raster are discarded rather
  // than silently passing through unmasked.
  return ee.ImageCollection(joined)
    .filter(ee.Filter.notNull(['cloud_probability_image']));
}

function addCloudBands(image) {
  image = ee.Image(image);
  var cloudProbability = ee.Image(image.get('cloud_probability_image'))
    .select('probability')
    .rename('cloud_probability');

  var clouds = cloudProbability.gt(PIX_CLOUD_PROB_THRESHOLD).rename('clouds');

  return image.addBands([cloudProbability, clouds]);
}

function addShadowBands(image) {
  image = ee.Image(image);

  var darkPixels = image.select('B8')
    .lt(NIR_DARK_THRESHOLD)
    .rename('dark_pixels');

  var shadowAzimuth = ee.Number(90)
    .subtract(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

  var cloudProjection = image.select('clouds')
    .directionalDistanceTransform(
      shadowAzimuth,
      CLOUD_SHADOW_DISTANCE_PIXELS
    )
    .select('distance')
    .mask()
    .rename('cloud_transform');

  var shadows = cloudProjection.and(darkPixels).rename('shadows');

  return image.addBands([darkPixels, cloudProjection, shadows]);
}

function applyCloudShadowMask(image) {
  image = ee.Image(image);

  var cloudShadow = image.select('clouds')
    .or(image.select('shadows'))
    .focal_min(CLOUD_ERODE_PIXELS, 'circle', 'pixels')
    .focal_max(CLOUD_DILATE_PIXELS, 'circle', 'pixels')
    .focal_max(CLOUD_SHADOW_BUFFER_METERS / 10, 'circle', 'pixels')
    .rename('cloud_shadow_mask');

  return image
    .addBands(cloudShadow)
    .updateMask(cloudShadow.not());
}


// ---------------------------------------------------------------
// 3.4. Apply static masks on every individual image
// ---------------------------------------------------------------
function applyStaticMasks(image, staticAllowedMask) {
  return ee.Image(image).updateMask(ee.Image(staticAllowedMask).eq(1));
}


// ---------------------------------------------------------------
// 3.5. Conservative B8 water-only sunglint guard
// ---------------------------------------------------------------
function applyGlintWaterMask(image) {
  image = ee.Image(image);

  if (!DO_GLINT_WATER_MASK) {
    return image.addBands(ee.Image(0).rename('glint_flag').toByte());
  }

  var preGlintMndwi = image.normalizedDifference(['B3', 'B11'])
    .rename('pre_glint_mndwi');

  var glintFlag = image.select('B8').gt(GLINT_B8_MIN)
    .and(preGlintMndwi.gt(GLINT_MNDWI_MIN))
    .rename('glint_flag')
    .toByte();

  return image
    .addBands([preGlintMndwi, glintFlag])
    .updateMask(glintFlag.not());
}


// ---------------------------------------------------------------
// 3.6. Otsu helpers and individual water mask
// ---------------------------------------------------------------
function otsu(histogram) {
  histogram = ee.Dictionary(histogram);

  var counts = ee.Array(histogram.get('histogram', [1, 1]));
  var means = ee.Array(histogram.get('bucketMeans', [-0.1, 0.1]));
  var size = ee.Number(means.length().get([0]));
  var total = ee.Number(counts.reduce(ee.Reducer.sum(), [0]).get([0]));
  var sum = ee.Number(means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]));
  var globalMean = sum.divide(total);

  var indices = ee.List.sequence(1, size.subtract(1));

  var bss = indices.map(function(i) {
    i = ee.Number(i);

    var aCounts = counts.slice(0, 0, i);
    var aCount = ee.Number(aCounts.reduce(ee.Reducer.sum(), [0]).get([0]));
    var aMeans = means.slice(0, 0, i);
    var aMean = ee.Number(aMeans.multiply(aCounts)
      .reduce(ee.Reducer.sum(), [0]).get([0]))
      .divide(aCount);

    var bCount = total.subtract(aCount);
    var bMean = ee.Number(ee.Algorithms.If(
      bCount.gt(0),
      sum.subtract(aCount.multiply(aMean)).divide(bCount),
      globalMean
    ));

    return ee.Number(ee.Algorithms.If(
      aCount.gt(0).and(bCount.gt(0)),
      aCount.multiply(aMean.subtract(globalMean).pow(2))
        .add(bCount.multiply(bMean.subtract(globalMean).pow(2))),
      -1
    ));
  });

  var candidateMeans = means.slice(0, 1, size);
  return ee.Number(ee.Array(candidateMeans).sort(bss).get([-1]));
}

function addWaterFeatures(image) {
  image = ee.Image(image);

  // All indices are retained as diagnostics. OTSU_WATER_INDEX_BAND selects
  // which one controls the binary land/water decision below.
  var mndwi = image.normalizedDifference(['B3','B11']).rename('mndwi');
  var ndwi = image.normalizedDifference(['B3','B8']).rename('ndwi');

  // AWEIsh is non-normalised; use the dedicated fallback/clamp configuration
  // above rather than reusing MNDWI/NDWI numeric expectations.
  var awei = image.expression(
    '4 * (green - swir1) - (0.25 * nir + 2.75 * swir2)',
    {
      green: image.select('B3'),
      nir: image.select('B8'),
      swir1: image.select('B11'),
      swir2: image.select('B12')
    }
  ).rename('awei');

  return image.addBands([mndwi, ndwi, awei]);
}

function selectedOtsuWaterIndex(image) {
  return ee.Image(image)
    .select(OTSU_WATER_INDEX_BAND)
    .rename('water_index');
}

// Returns compact threshold diagnostics so `otsu_used` is recorded directly
// rather than inferred from a possible numerical equality with fallback.
function getOtsuThresholdInfo(waterIndex, regionGeom) {
  waterIndex = ee.Image(waterIndex).rename('water_index');

  var validPixels = ee.Number(waterIndex.mask().reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: regionGeom,
    scale: OTSU_SCALE,
    maxPixels: EXPORT_MAX_PIXELS,
    tileScale: OTSU_TILE_SCALE
  }).get('water_index', 0));

  var reduced = waterIndex.reduceRegion({
    reducer: ee.Reducer.histogram({
      maxBuckets: OTSU_MAX_BUCKETS,
      minBucketWidth: OTSU_MIN_BUCKET_WIDTH
    }),
    geometry: regionGeom,
    scale: OTSU_SCALE,
    bestEffort: true,
    maxPixels: EXPORT_MAX_PIXELS,
    tileScale: OTSU_TILE_SCALE
  });

  var hist = ee.Dictionary(reduced.get('water_index', ee.Dictionary({})));
  var hasHistogram = hist.contains('histogram');
  var hasEnoughValidPixels = validPixels.gte(OTSU_MIN_VALID_PIXELS);

  var canUseOtsu = ee.Number(ee.Algorithms.If(
    hasHistogram,
    hasEnoughValidPixels,
    0
  ));

  var otsuUsed = ee.Number(USE_OTSU_WATER_THRESHOLD ? 1 : 0)
    .multiply(canUseOtsu);

  var rawThreshold = ee.Number(ee.Algorithms.If(
    otsuUsed.eq(1),
    otsu(hist),
    WATER_FALLBACK_THRESHOLD
  ));

  var threshold = rawThreshold.max(OTSU_CLAMP_MIN).min(OTSU_CLAMP_MAX);

  return ee.Dictionary({
    raw_threshold: rawThreshold,
    threshold: threshold,
    otsu_used: otsuUsed,
    valid_pixel_count: validPixels,
    histogram_available: ee.Number(ee.Algorithms.If(hasHistogram, 1, 0))
  });
}

function getOtsuThreshold(waterIndex, regionGeom) {
  return ee.Number(
    getOtsuThresholdInfo(waterIndex, regionGeom).get('threshold')
  );
}

// ----------------------------------------------------------------
// 3.6.1 Preserve clear land/water observations BEFORE masking water
// ----------------------------------------------------------------
// analysis_observed is available after ROI/tree, cloud/shadow and glint
// masks, but before the Otsu land/water mask. It therefore distinguishes
// “observed as water” from “not observed because of cloud/shadow/glint”.
function addLandWaterObservationBands(image) {
  image = ee.Image(image);

  var waterIndex = selectedOtsuWaterIndex(image);

  var analysisObserved = commonValidMask(image, SPECTRAL_BANDS)
    .rename('analysis_observed')
    .unmask(0)
    .toByte()
    .clip(roi);

  var otsuInfo = getOtsuThresholdInfo(waterIndex, roi);
  var otsuThreshold = ee.Number(otsuInfo.get('threshold'));
  var otsuUsedFlag = ee.Number(otsuInfo.get('otsu_used'));
  var finalThreshold = otsuThreshold.add(OTSU_LAND_BUFFER);

  var landObserved = analysisObserved.eq(1)
    .and(waterIndex.lte(finalThreshold))
    .rename('land_observed')
    .unmask(0)
    .toByte();

  var waterObserved = analysisObserved.eq(1)
    .and(waterIndex.gt(finalThreshold))
    .rename('water_observed')
    .unmask(0)
    .toByte();

  // Contract note: `water_threshold` stores the final threshold of the
  // SELECTED Otsu index; see water_index_for_otsu image properties.
  var waterThreshold = ee.Image.constant(finalThreshold)
    .rename('water_threshold')
    .toFloat()
    .updateMask(analysisObserved.eq(1));

  var otsuUsed = ee.Image.constant(otsuUsedFlag)
    .rename('otsu_used')
    .toByte()
    .updateMask(analysisObserved.eq(1));

  return image.addBands([
    analysisObserved,
    landObserved,
    waterObserved,
    waterThreshold,
    otsuUsed
  ]).set({
    water_threshold: finalThreshold,
    otsu_threshold: otsuThreshold,
    otsu_threshold_raw: otsuInfo.get('raw_threshold'),
    otsu_used: otsuUsedFlag,
    water_index_for_otsu: OTSU_WATER_INDEX,
    water_index_band_for_otsu: OTSU_WATER_INDEX_BAND,
    water_index_label_for_otsu: OTSU_WATER_INDEX_LABEL,
    water_threshold_definition:
      'selected_otsu_index_threshold_plus_land_buffer'
  });
}

// ----------------------------------------------------------------
// 3.6.2 Use preserved land_observed to mask the index time series
// ----------------------------------------------------------------
function applyIndividualWaterMask(image) {
  image = ee.Image(image);

  var landMask = image.select('land_observed')
    .rename('land_mask')
    .toByte();

  return image
    .addBands(landMask, null, true)
    .updateMask(landMask.eq(1));
}


// ---------------------------------------------------------------
// 3.7. Individual time-series vegetation indices
// ---------------------------------------------------------------
function addTimeSeriesIndices(image) {
  image = ee.Image(image);

  var ndvi = image.normalizedDifference(['B8','B4']).rename('ndvi');

  var savi = image.expression(
    '1.5 * ((nir - red) / (nir + red + 0.5))',
    {nir: image.select('B8'), red: image.select('B4')}
  ).rename('savi');

  var gndvi = image.normalizedDifference(['B8','B3']).rename('gndvi');
  var ndmi = image.normalizedDifference(['B8','B11']).rename('ndmi');
  var grvi = image.normalizedDifference(['B3','B4']).rename('grvi');

  var indexImage = image.addBands([ndvi, savi, gndvi, ndmi, grvi]);
  var valid = commonValidMask(indexImage, TIME_SERIES_INDEX_BANDS)
    .rename('valid')
    .unmask(0)
    .toByte();

  return indexImage.addBands(valid);
}


// ================================================================
// 4. BUILD INDIVIDUAL-IMAGE COLLECTION
// ================================================================

var s2Base = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUDY_PIXEL_PERCENTAGE_MAX));

var staticMasks = getStaticMasks();

// Spatial denominator shared by all daily images. This lets the script
// distinguish a truly useful image from one with only a few surviving pixels.
staticAllowedPixelCount = ee.Number(
  staticMasks.staticAllowed.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: EXPORT_SCALE,
    maxPixels: EXPORT_MAX_PIXELS,
    tileScale: 4
  }).get('static_allowed_mask', 0)
);

print('Static allowed pixels (coverage denominator):', staticAllowedPixelCount);

function addRun01ImageMetadata(image) {
  image = ee.Image(image);

  return image.set({
    est_id: EST_ID_QUERY,
    season_tag: SEASON_TAG,
    year: YEAR,
    date_key: ee.Date(image.get('system:time_start')).format('YYYYMMdd'),
    source_system_index: image.get('system:index')
  });
}

// Build the date list from the cloud-probability-matched collection BEFORE
// per-pixel Otsu processing. This keeps the callback lightweight and prevents
// an index-export failure when land/water diagnostics are computationally heavy.
var s2WithCloudProbability = linkCloudProbability(s2Base)
  .map(addRun01ImageMetadata);

var landWaterDateKeysServer = ee.List(
  s2WithCloudProbability.aggregate_array('date_key')
).distinct().sort();

var s2PreWater = s2WithCloudProbability
  .map(applyScaleFactors)
  .map(sharpened)
  .map(selectSharpenedSpectralBands)
  .map(addCloudBands)
  .map(addShadowBands)
  .map(applyCloudShadowMask)
  .map(function(image) {
    return applyStaticMasks(image, staticMasks.staticAllowed);
  })
  .map(applyGlintWaterMask);

// This branch preserves the water state for RUN04/RUN06.
var s2LandWater = s2PreWater
  .map(addWaterFeatures)
  .map(addLandWaterObservationBands);

// The original RUN02 input branch remains semantically unchanged: only land
// observations contribute to vegetation-index time series.
var s2Processed = s2LandWater
  .map(applyIndividualWaterMask)
  .map(addTimeSeriesIndices)
  .map(addRun01ImageMetadata);

print('RUN01_V3.3 estuary:', EST_ID_QUERY);
print('RUN01_V3.3 season:', SEASON_TAG, START_DATE, END_DATE);
print('S2 scenes before cloud-probability join:', s2Base.size());
print('S2 scenes with cloud-probability match:', s2WithCloudProbability.size());
print('S2 scenes after full individual index processing:', s2Processed.size());

// Do not call getInfo() here. Earth Engine uses deferred execution, and the
// processed-scene count is already printed above as a server-side diagnostic.
// If no images remain, later outputs will reveal the empty collection.


// ================================================================
// 5. BUILD ONE TEMPORAL IMAGE PER ACQUISITION DATE
//    Handles same-day overlap from adjacent Sentinel-2 tiles.
// ================================================================

function buildDailyMedian(dateKey) {
  dateKey = ee.String(dateKey);

  var sameDate = s2Processed.filter(ee.Filter.eq('date_key', dateKey));
  var dailyIndices = sameDate.select(TIME_SERIES_INDEX_BANDS).median();
  var dailyMask = commonValidMask(dailyIndices, TIME_SERIES_INDEX_BANDS);

  dailyIndices = dailyIndices.updateMask(dailyMask);

  var valid = dailyMask.rename('valid').unmask(0).toByte();
  var dailyDate = ee.Date.parse('YYYYMMdd', dateKey);

  return dailyIndices.addBands(valid)
    .set({
      'system:time_start': dailyDate.millis(),
      date_key: dateKey,
      est_id: EST_ID_QUERY,
      season_tag: SEASON_TAG,
      year: YEAR,
      source_images_same_date: sameDate.size()
    });
}

var dateKeys = ee.List(s2Processed.aggregate_array('date_key')).distinct().sort();
var timeSeriesDaily = ee.ImageCollection.fromImages(
  dateKeys.map(buildDailyMedian)
).map(addImageValidCount)
 .filter(ee.Filter.gt('valid_pixel_count', 0))
 .sort('system:time_start');

print('Unique acquisition dates with valid pixels:', timeSeriesDaily.size());
print('Time-series dates:', timeSeriesDaily.aggregate_array('date_key'));

// Keep date keys server-side here. They are transferred to the browser ONCE,
// asynchronously, only at the end of the script to create export Tasks.
var timeSeriesDateKeysServer = timeSeriesDaily
  .sort('system:time_start')
  .aggregate_array('date_key');


// ----------------------------------------------------------------
// Daily land/water image per source date
// ----------------------------------------------------------------
// “Land wins” for overlapping same-day Sentinel-2 tiles. The `valid` band is
// intentionally set to the daily land mask, matching the eligibility rule of
// the vegetation-index collection without rebuilding that costly collection.
function buildDailyLandWater(dateKey) {
  dateKey = ee.String(dateKey);

  var sameDate = s2LandWater.filter(ee.Filter.eq('date_key', dateKey));
  var dailyDate = ee.Date.parse('YYYYMMdd', dateKey);

  var analysisObserved = sameDate.select('analysis_observed')
    .max()
    .unmask(0)
    .rename('analysis_observed')
    .toByte();

  var landObserved = sameDate.select('land_observed')
    .max()
    .unmask(0)
    .rename('land_observed')
    .toByte();

  var waterObserved = analysisObserved.eq(1)
    .and(landObserved.eq(0))
    .rename('water_observed')
    .toByte();

  var dailyMndwi = sameDate.select('mndwi')
    .median()
    .rename('mndwi')
    .updateMask(analysisObserved.eq(1))
    .toFloat();

  var dailyWaterThreshold = sameDate.select('water_threshold')
    .median()
    .rename('water_threshold')
    .updateMask(analysisObserved.eq(1))
    .toFloat();

  var dailyOtsuUsed = sameDate.select('otsu_used')
    .max()
    .rename('otsu_used')
    .updateMask(analysisObserved.eq(1))
    .toByte();

  // Equivalent eligibility definition for the index series: a clear,
  // static, non-glint observation classified as land on that day.
  var indexEligible = landObserved.rename('valid').toByte();

  return ee.Image.cat([
    analysisObserved,
    landObserved,
    waterObserved,
    dailyMndwi,
    dailyWaterThreshold,
    dailyOtsuUsed,
    indexEligible
  ]).set({
    'system:time_start': dailyDate.millis(),
    date_key: dateKey,
    est_id: EST_ID_QUERY,
    season_tag: SEASON_TAG,
    year: YEAR,
    source_images_same_date: sameDate.size(),
    aggregation_rule: 'land_wins_same_date',
    index_valid_definition: 'daily_land_observed',
    water_index_for_otsu: OTSU_WATER_INDEX,
    water_index_band_for_otsu: OTSU_WATER_INDEX_BAND,
    water_index_label_for_otsu: OTSU_WATER_INDEX_LABEL,
    water_threshold_definition:
      'selected_otsu_index_threshold_plus_land_buffer'
  });
}

var timeSeriesLandWaterDaily = ee.ImageCollection.fromImages(
  landWaterDateKeysServer.map(buildDailyLandWater)
).sort('system:time_start');


// ================================================================
// 6. SEASONAL COMPOSITE, PCA, NDVI AND QA
// ================================================================

// Raw seasonal spectral composite before the second composite land mask.
var seasonalCompositeRaw = s2PreWater
  .select(SPECTRAL_BANDS)
  .median()
  .clip(roi);

var exportRefProjection = seasonalCompositeRaw.select('B2').projection();

function gaussianSmoothSpectral(image) {
  var smoothed = ee.Image(image).select(SPECTRAL_BANDS).convolve(
    ee.Kernel.gaussian({
      radius: GAUSSIAN_RADIUS_PIXELS,
      sigma: GAUSSIAN_SIGMA_PIXELS,
      units: 'pixels',
      normalize: true
    })
  ).rename(SPECTRAL_BANDS);

  return smoothed.setDefaultProjection(exportRefProjection);
}

var seasonalCompositeSmooth = ee.Image(ee.Algorithms.If(
  DO_SEASONAL_GAUSSIAN_BLUR,
  gaussianSmoothSpectral(seasonalCompositeRaw),
  seasonalCompositeRaw.select(SPECTRAL_BANDS)
)).clip(roi).setDefaultProjection(exportRefProjection);

// Composite water mask is a final safeguard for PCA and fallback NDVI class.
// The selected index is controlled by OTSU_WATER_INDEX; MNDWI/NDWI/AWEI are
// all still available in compositeWaterFeatures for visual comparison.
var compositeWaterFeatures = addWaterFeatures(seasonalCompositeSmooth);
var compositeSelectedWaterIndex = selectedOtsuWaterIndex(compositeWaterFeatures)
  .rename('composite_selected_water_index');

var compositeOtsuInfo = getOtsuThresholdInfo(
  compositeSelectedWaterIndex, roi
);
var compositeOtsuThresholdRaw = ee.Number(
  compositeOtsuInfo.get('raw_threshold')
);
var compositeOtsuThreshold = ee.Number(
  compositeOtsuInfo.get('threshold')
);
var compositeOtsuUsed = ee.Number(compositeOtsuInfo.get('otsu_used'));
var compositeFinalThreshold = compositeOtsuThreshold.add(OTSU_LAND_BUFFER);

var compositeMndwi = compositeWaterFeatures.select('mndwi')
  .rename('composite_mndwi');
var compositeNdwi = compositeWaterFeatures.select('ndwi')
  .rename('composite_ndwi');
var compositeAwei = compositeWaterFeatures.select('awei')
  .rename('composite_awei');

// Otsu-only land mask and the added transition strip are kept only for QA.
// With a positive buffer, pixels in compositeBufferRetainedLand would be
// water under raw Otsu but are deliberately retained as land in RUN01.
var compositeLandMaskOtsuOnly = compositeSelectedWaterIndex
  .lte(compositeOtsuThreshold)
  .rename('land_mask_otsu_only')
  .toByte();

var compositeBufferRetainedLand = compositeSelectedWaterIndex
  .gt(compositeOtsuThreshold)
  .and(compositeSelectedWaterIndex.lte(compositeFinalThreshold))
  .rename('buffer_retained_land')
  .toByte();

var compositeLandMask = compositeSelectedWaterIndex
  .lte(compositeFinalThreshold)
  .rename('land_mask_composite')
  .toByte();

var compositeBufferRetainedLandCount = ee.Number(
  compositeBufferRetainedLand.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: EXPORT_SCALE,
    maxPixels: EXPORT_MAX_PIXELS,
    tileScale: 4
  }).get('buffer_retained_land', 0)
);

var seasonalWorkingImage = ee.Image(ee.Algorithms.If(
  APPLY_COMPOSITE_LAND_MASK,
  seasonalCompositeSmooth.updateMask(compositeLandMask.eq(1)),
  seasonalCompositeSmooth
)).updateMask(staticMasks.staticAllowed.eq(1))
  .clip(roi)
  .setDefaultProjection(exportRefProjection);

var seasonalNdvi = seasonalWorkingImage
  .normalizedDifference(['B8','B4'])
  .rename('ndvi')
  .setDefaultProjection(exportRefProjection);

function classifyNdvi(ndviImage) {
  ndviImage = ee.Image(ndviImage);

  return ee.Image(0)
    .where(ndviImage.gte(0.0).and(ndviImage.lt(0.2)), 1)
    .where(ndviImage.gte(0.2).and(ndviImage.lt(0.4)), 2)
    .where(ndviImage.gte(0.4).and(ndviImage.lt(0.6)), 3)
    .where(ndviImage.gte(0.6).and(ndviImage.lt(0.8)), 4)
    .where(ndviImage.gte(0.8), 5)
    .rename('ndvi_class')
    .toByte();
}

var seasonalNdviClass = classifyNdvi(seasonalNdvi)
  .unmask(0)
  .rename('ndvi_class')
  .toByte();


// ---------------------------------------------------------------
// 6.1. PCA on masked seasonal composite
// ---------------------------------------------------------------
var PCA_INFO = {eigenValues: null, variancePercent: null, cumulativePercent: null};

function addPCA(image, regionGeom) {
  image = ee.Image(image);

  var pcaImage = image.select(PCA_BANDS).toFloat();
  var validMask = commonValidMask(pcaImage, PCA_BANDS);
  pcaImage = pcaImage.updateMask(validMask);

  var meansDict = pcaImage.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: regionGeom,
    scale: EXPORT_SCALE,
    bestEffort: true,
    maxPixels: EXPORT_MAX_PIXELS,
    tileScale: 4
  });

  var means = ee.Image.constant(PCA_BANDS.map(function(bandName) {
    return ee.Number(ee.Algorithms.If(meansDict.get(bandName), meansDict.get(bandName), 0));
  })).rename(PCA_BANDS);

  var centered = pcaImage.subtract(means).updateMask(validMask);

  var covarianceDict = centered.toArray().reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: regionGeom,
    scale: EXPORT_SCALE,
    bestEffort: true,
    maxPixels: EXPORT_MAX_PIXELS,
    tileScale: 4
  });

  var covariance = ee.Array(covarianceDict.get('array'));
  var eigens = covariance.eigen();
  var eigenValues = eigens.slice(1, 0, 1);
  var eigenVectors = eigens.slice(1, 1);

  var eigenList = ee.Array(eigenValues).project([0]).toList();
  var eigenTotal = ee.Number(eigenList.reduce(ee.Reducer.sum()));
  var variancePercent = eigenList.map(function(v) {
    return ee.Number(v).divide(eigenTotal).multiply(100);
  });
  var cumulativePercent = ee.List(variancePercent).iterate(function(v, previous) {
    previous = ee.List(previous);
    var last = ee.Number(ee.Algorithms.If(previous.size().gt(0), previous.get(-1), 0));
    return previous.add(last.add(ee.Number(v)));
  }, ee.List([]));

  PCA_INFO.eigenValues = eigenList;
  PCA_INFO.variancePercent = variancePercent;
  PCA_INFO.cumulativePercent = cumulativePercent;

  var arrayImage = centered.toArray().toArray(1);
  var pcNames = ee.List.sequence(1, PCA_BANDS.length).map(function(i) {
    return ee.String('PC').cat(ee.Number(i).int());
  });

  var pcs = ee.Image(eigenVectors)
    .matrixMultiply(arrayImage)
    .arrayProject([0])
    .arrayFlatten([pcNames]);

  return image.addBands(pcs.toFloat(), null, true);
}

// Seasonal PCA validity diagnostic. Keep it server-side to avoid a synchronous
// getInfo() request; it will be shown in the Console by print().
var seasonalValidMask = commonValidMask(seasonalWorkingImage, PCA_BANDS);
var seasonalValidCount = ee.Number(seasonalValidMask.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: roi,
  scale: EXPORT_SCALE,
  maxPixels: EXPORT_MAX_PIXELS,
  tileScale: 4
}).get('valid_mask', 0));

var seasonalPca = addPCA(seasonalWorkingImage, roi)
  .select(['PC1','PC2','PC3'])
  .setDefaultProjection(exportRefProjection)
  .clip(roi);


// ---------------------------------------------------------------
// 6.2. Seasonal QA stack
// ---------------------------------------------------------------
var obsCount = timeSeriesDaily.select('valid')
  .sum()
  .rename('obs_count')
  .unmask(0)
  .toUint16();

var validSeason = obsCount.gt(0)
  .and(staticMasks.staticAllowed.eq(1))
  .and(compositeLandMask.eq(1))
  .rename('valid_season')
  .unmask(0)
  .toByte();

var nSourceImages = ee.Image.constant(s2Processed.size())
  .rename('n_source_images')
  .clip(roi)
  .toUint16();

var nDatesExported = ee.Image.constant(timeSeriesDaily.size())
  .rename('n_dates_exported')
  .clip(roi)
  .toUint16();

var compositeWaterThreshold = ee.Image.constant(compositeFinalThreshold)
  .rename('composite_water_threshold')
  .clip(roi)
  .toFloat();

var seasonalQa = obsCount
  .addBands(staticMasks.staticAllowed.rename('static_allowed_mask').unmask(0).toByte())
  .addBands(compositeLandMask.unmask(0).toByte())
  .addBands(validSeason)
  .addBands(nSourceImages)
  .addBands(nDatesExported)
  .addBands(compositeWaterThreshold)
  .set(baseProductMetadata('seasonal_qa', QA_BAND_SCHEMA))
  .set({
    n_source_images: s2Processed.size(),
    n_dates_exported: timeSeriesDaily.size(),
    composite_water_threshold: compositeFinalThreshold,
    composite_water_index_for_otsu: OTSU_WATER_INDEX,
    composite_water_index_band_for_otsu: OTSU_WATER_INDEX_BAND,
    composite_otsu_threshold_raw: compositeOtsuThresholdRaw,
    composite_otsu_threshold: compositeOtsuThreshold,
    composite_otsu_used: compositeOtsuUsed
  });


// ================================================================
// 7. MAPS AND DEBUG
// ================================================================

if (DO_ADD_MAPS) {
  /*
  Map.addLayer(staticMasks.classificationArea.selfMask(),
   {min: 0, max: 1, palette: ['00ffff']},
    'Filtered classification area (QGIS)', false);
  Map.addLayer(staticMasks.treeExcluded.selfMask(),
    {min: 0, max: 1, palette: ['006400']},
    'Tree > 2 m excluded', false);
  Map.addLayer(staticMasks.staticAllowed.selfMask(),
    {min: 0, max: 1, palette: ['ffff00']},
    'Static allowed mask: filtered area minus trees', false);
*/
  Map.addLayer(seasonalCompositeRaw.select(RGB_BANDS), VIS_RGB,
    'Seasonal RGB pre-water-mask', false);
  //Map.addLayer(seasonalWorkingImage.select(RGB_BANDS), VIS_RGB,
  //  'Seasonal RGB masked / smoothed', false);

  if (QA_ADD_OTSU_BUFFER_LAYERS) {
    var selectedWaterIndexVis = {
      min: OTSU_CLAMP_MIN,
      max: OTSU_CLAMP_MAX,
      palette: ['8c510a', 'f6e8c3', 'f7f7f7', 'c7eae5', '01665e']
    };

    Map.addLayer(compositeSelectedWaterIndex, selectedWaterIndexVis,
      'QA Otsu | selected index: ' + OTSU_WATER_INDEX, false);
    Map.addLayer(compositeLandMaskOtsuOnly.selfMask(),
      {min: 0, max: 1, palette: ['00bcd4']},
      'QA Otsu | composite land before buffer', false);
    Map.addLayer(compositeBufferRetainedLand.selfMask(),
      {min: 0, max: 1, palette: ['ffeb00']},
      'QA Otsu | land retained only by buffer', true);

    if (QA_ADD_ALL_WATER_INDEX_LAYERS) {
      /*
      Map.addLayer(compositeMndwi, {
        min: -0.30, max: 0.30,
        palette: ['8c510a', 'f6e8c3', 'f7f7f7', 'c7eae5', '01665e']
      }, 'QA Otsu comparison | MNDWI', false);
      Map.addLayer(compositeNdwi, {
        min: -0.30, max: 0.30,
        palette: ['8c510a', 'f6e8c3', 'f7f7f7', 'c7eae5', '01665e']
      }, 'QA Otsu comparison | NDWI', false);
      Map.addLayer(compositeAwei, {
        min: -1.00, max: 1.00,
        palette: ['8c510a', 'f6e8c3', 'f7f7f7', 'c7eae5', '01665e']
      }, 'QA Otsu comparison | AWEIsh', false);
    }
*/
  }

  Map.addLayer(compositeLandMask.selfMask(),
    {min: 0, max: 1, palette: ['00ffff']}, 'Composite land mask', true);
  Map.addLayer(seasonalNdvi, VIS_NDVI, 'Seasonal NDVI', false);
  Map.addLayer(seasonalNdviClass.updateMask(seasonalNdviClass.neq(0)),
    VIS_NDVI_CLASS, 'Seasonal NDVI class (fallback)', false);
  /*
  Map.addLayer(seasonalPca.select('PC1'), VIS_PCA, 'PCA PC1', false);
  Map.addLayer(seasonalPca.select('PC2'), VIS_PCA, 'PCA PC2', false);
  Map.addLayer(seasonalPca.select('PC3'), VIS_PCA, 'PCA PC3', false);
  Map.addLayer(obsCount, {min: 0, max: 15, palette: ['000000','2c7bb6','ffffbf','d7191c']},
    'Observation count (dates)', false);
  Map.addLayer(validSeason.selfMask(),
    {min: 0, max: 1, palette: ['00ff00']}, 'Valid seasonal domain', false);
*/
  }
}

// First-image temporal index QC is added asynchronously together with the
// export tasks, after the one permitted client-side transfer of date keys.
//
if (DEBUG_LEVEL >= 1) {
  //print('PCA eigenvalues:', PCA_INFO.eigenValues);
  //print('PCA variance (%) by component:', PCA_INFO.variancePercent);
  print('PCA cumulative variance (%):', PCA_INFO.cumulativePercent);
  print('Composite Otsu selected water index:', {
    index: OTSU_WATER_INDEX,
    band: OTSU_WATER_INDEX_BAND,
    label: OTSU_WATER_INDEX_LABEL
  });
  print('Composite Otsu threshold (raw before clamp):', compositeOtsuThresholdRaw);
  print('Composite Otsu threshold (clamped):', compositeOtsuThreshold);
  print('Composite Otsu + buffer threshold:', compositeFinalThreshold);
  print('Composite Otsu used:', compositeOtsuUsed);
  print('Seasonal valid pixels for PCA:', seasonalValidCount);

  if (QA_PRINT_OTSU_BUFFER_SUMMARY) {
    print('QA Otsu-buffer diagnostic:', {
      otsu_water_index: OTSU_WATER_INDEX,
      otsu_water_band: OTSU_WATER_INDEX_BAND,
      otsu_fallback_threshold: WATER_FALLBACK_THRESHOLD,
      otsu_clamp_min: OTSU_CLAMP_MIN,
      otsu_clamp_max: OTSU_CLAMP_MAX,
      otsu_land_buffer: OTSU_LAND_BUFFER,
      composite_otsu_threshold_raw: compositeOtsuThresholdRaw,
      composite_otsu_threshold: compositeOtsuThreshold,
      composite_final_land_threshold: compositeFinalThreshold,
      composite_otsu_used: compositeOtsuUsed,
      pixels_retained_only_by_buffer: compositeBufferRetainedLandCount
    });
  }
}

if (DEBUG_LEVEL >= 2) {
  printProjection('Seasonal composite B2', seasonalWorkingImage.select('B2'));
  printProjection('PCA3', seasonalPca.select('PC1'));
  print('Processed source-image threshold values:', s2Processed.aggregate_array('water_threshold'));
  print('Temporal valid-pixel counts:', timeSeriesDaily.aggregate_array('valid_pixel_count'));
  print('Temporal valid fractions of static area:',
    timeSeriesDaily.aggregate_array('valid_fraction_static'));
  print('RUN02 expected temporal band schema:', TS_BAND_SCHEMA);
  print('RUN02 expected seasonal schemas:', {
    composite: COMPOSITE_BAND_SCHEMA,
    pca3: PCA_BAND_SCHEMA,
    qa: QA_BAND_SCHEMA
  });
}


// ================================================================
// 8. EXPORTS
// ================================================================

var prefix = EST_ID_QUERY + '_' + SEASON_TAG;

var compositeExport = seasonalWorkingImage.select(SPECTRAL_BANDS)
  .toFloat()
  .set(baseProductMetadata('seasonal_composite', COMPOSITE_BAND_SCHEMA));

var pcaExport = seasonalPca.select(['PC1','PC2','PC3'])
  .toFloat()
  .set(baseProductMetadata('seasonal_pca3', PCA_BAND_SCHEMA));

var ndviExport = seasonalNdvi.rename('ndvi')
  .toFloat()
  .set(baseProductMetadata('seasonal_ndvi', 'ndvi'));

var ndviClassExport = seasonalNdviClass.rename('ndvi_class')
  .unmask(0)
  .toByte()
  .set(baseProductMetadata('seasonal_ndvi_class', 'ndvi_class'));

var rgbRawExport = seasonalCompositeRaw.select(RGB_BANDS)
  .toFloat()
  .set(baseProductMetadata('seasonal_rgb_raw_qa', 'B4,B3,B2'));

var rgbMaskedExport = seasonalWorkingImage.select(RGB_BANDS)
  .toFloat()
  .set(baseProductMetadata('seasonal_rgb_masked_qa', 'B4,B3,B2'));

function createSeasonalAssetExportTasks() {
  if (!DO_EXPORT_ASSET) return;

  // 8.1 Seasonal products
  /*
  Export.image.toAsset(makeExportArgs(
    compositeExport,
    'R01V3_COMP_' + prefix,
    assetPath('RUN01V3_COMPOSITE'),
    roi,
    {'.default': 'mean'}
  ));
  */

  Export.image.toAsset(makeExportArgs(
    pcaExport,
    'R01V3_PCA3_' + prefix,
    assetPath('RUN01V3_PCA3'),
    roi,
    {PC1: 'mean', PC2: 'mean', PC3: 'mean', '.default': 'mean'}
  ));

  Export.image.toAsset(makeExportArgs(
    ndviExport,
    'R01V3_NDVI_' + prefix,
    assetPath('RUN01V3_NDVI'),
    roi,
    {ndvi: 'mean', '.default': 'mean'}
  ));

  Export.image.toAsset(makeExportArgs(
    ndviClassExport,
    'R01V3_NDVIC_' + prefix,
    assetPath('RUN01V3_NDVI_CLASS'),
    roi,
    {ndvi_class: 'mode', '.default': 'mode'}
  ));

  Export.image.toAsset(makeExportArgs(
    seasonalQa,
    'R01V3_QA_' + prefix,
    assetPath('RUN01V3_QA'),
    roi,
    {
      obs_count: 'mean',
      static_allowed_mask: 'max',
      land_mask_composite: 'max',
      valid_season: 'max',
      n_source_images: 'mean',
      n_dates_exported: 'mean',
      composite_water_threshold: 'mean',
      '.default': 'mean'
    }
  ));

  if (EXPORT_RGB_QA) {
    Export.image.toAsset(makeExportArgs(
      rgbRawExport,
      'R01V3_RGBRAW_' + prefix,
      assetPath('RUN01V3_RGB_RAW'),
      roi,
      {'.default': 'mean'}
    ));

    Export.image.toAsset(makeExportArgs(
      rgbMaskedExport,
      'R01V3_RGBMASK_' + prefix,
      assetPath('RUN01V3_RGB_MASKED'),
      roi,
      {'.default': 'mean'}
    ));
  }

  print('RUN01 seasonal Asset export Tasks created for:', prefix);
  print('RUN01 index and land/water time-series Tasks are created by independent asynchronous callbacks.');
  print('RUN01 index time-series parent ImageCollection:', TIME_SERIES_COLLECTION_ASSET);
  print('RUN01 land/water time-series parent ImageCollection:',
    LAND_WATER_TIME_SERIES_COLLECTION_ASSET);
}

// ----------------------------------------------------------------
// 8.2 Individual DATE image exports to the pre-created ImageCollection
// ----------------------------------------------------------------
// Export.image.toAsset() needs client-side strings for task description and
// assetId. Therefore, date keys must cross from the server to the browser.
// Use ONE asynchronous evaluate() request instead of calling getInfo() once
// per image. This is lighter and avoids blocking the Code Editor UI.
function createTimeSeriesExportTasks(dateKeysClient) {
  dateKeysClient = dateKeysClient || [];

  if (dateKeysClient.length === 0) {
    print('WARNING: no valid time-series date keys were returned. No TS export Tasks created.');
    return;
  }

  var nToExport = MAX_TS_EXPORTS_PER_RUN < 0 ?
    dateKeysClient.length :
    Math.min(MAX_TS_EXPORTS_PER_RUN, dateKeysClient.length);

  for (var i = 0; i < nToExport; i++) {
    var dateKeyClient = String(dateKeysClient[i]);

    // Re-select the server-side image by a client-side date key. This preserves
    // all spectral/index processing lazily on the Earth Engine server.
    var tsImage = ee.Image(
      timeSeriesDaily
        .filter(ee.Filter.eq('date_key', dateKeyClient))
        .first()
    );

    var tsName = 'RUN01V3_TS_' + EST_ID_QUERY + '_' +
      SEASON_TAG + '_' + dateKeyClient;

    var tsExport = tsImage.select(TIME_SERIES_INDEX_BANDS.concat(['valid']))
      .set(baseProductMetadata('timeseries_daily_indices', TS_BAND_SCHEMA))
      .set({
        'system:time_start': tsImage.get('system:time_start'),
        date_key: dateKeyClient,
        source_images_same_date: tsImage.get('source_images_same_date'),
        valid_pixel_count: tsImage.get('valid_pixel_count'),
        valid_fraction_static: tsImage.get('valid_fraction_static')
      });

    Export.image.toAsset(makeExportArgs(
      tsExport,
      'R01V3_TS_' + SEASON_TAG + '_' + dateKeyClient,
      TIME_SERIES_COLLECTION_ASSET + '/' + tsName,
      roi,
      {
        ndvi: 'mean',
        savi: 'mean',
        gndvi: 'mean',
        ndmi: 'mean',
        grvi: 'mean',
        valid: 'max',
        '.default': 'mean'
      }
    ));
  }

  print('Time-series date keys received:', dateKeysClient);
  print('Time-series ImageCollection export Tasks created:', nToExport);
}


function createLandWaterTimeSeriesExportTasks(dateKeysClient) {
  dateKeysClient = dateKeysClient || [];

  if (dateKeysClient.length === 0) {
    print('WARNING: no land/water date keys were returned. No land/water Tasks created.');
    return;
  }

  var nToExport = MAX_LAND_WATER_TS_EXPORTS_PER_RUN < 0 ?
    dateKeysClient.length :
    Math.min(MAX_LAND_WATER_TS_EXPORTS_PER_RUN, dateKeysClient.length);

  for (var i = 0; i < nToExport; i++) {
    var dateKeyClient = String(dateKeysClient[i]);

    var lwImage = ee.Image(
      timeSeriesLandWaterDaily
        .filter(ee.Filter.eq('date_key', dateKeyClient))
        .first()
    );

    var lwName = 'RUN01V3_LW_' + EST_ID_QUERY + '_' +
      SEASON_TAG + '_' + dateKeyClient;

    var lwExport = lwImage.select(LAND_WATER_TS_BANDS)
      .set(baseProductMetadata(
        'timeseries_daily_landwater',
        LAND_WATER_TS_BAND_SCHEMA
      ))
      .set({
        'system:time_start': lwImage.get('system:time_start'),
        date_key: dateKeyClient,
        source_images_same_date: lwImage.get('source_images_same_date'),
        aggregation_rule: lwImage.get('aggregation_rule'),
        index_valid_definition: lwImage.get('index_valid_definition'),
        index_time_series_collection: TIME_SERIES_COLLECTION_ASSET
      });

    Export.image.toAsset(makeExportArgs(
      lwExport,
      'R01V3_LW_' + SEASON_TAG + '_' + dateKeyClient,
      LAND_WATER_TIME_SERIES_COLLECTION_ASSET + '/' + lwName,
      roi,
      {
        analysis_observed: 'max',
        land_observed: 'max',
        water_observed: 'max',
        mndwi: 'mean',
        water_threshold: 'mean',
        otsu_used: 'max',
        valid: 'max',
        '.default': 'mean'
      }
    ));
  }

  print('Land/water date keys received:', dateKeysClient);
  print('Land/water ImageCollection export Tasks created:', nToExport);
}

// ----------------------------------------------------------------
// 8.3 Representative temporal image: visual QA for exported index content
// ----------------------------------------------------------------
// IMPORTANT:
// The temporal collection retains ALL dates with at least one valid pixel.
// That is correct for pixel-wise Fourier fitting because a sparse date can
// still contribute a real observation at the pixels where it is valid.
//
// For map-based visual QA, however, we select the DAILY image with the
// highest valid_fraction_static. This avoids displaying the earliest date
// when it contains only a few isolated unmasked pixels.
var emptyQcTsImage = ee.Image.constant([0, 0, 0, 0, 0, 0])
  .rename(TIME_SERIES_INDEX_BANDS.concat(['valid']))
  .updateMask(ee.Image(0))
  .set({
    date_key: 'NO_VALID_DATE',
    valid_pixel_count: 0,
    valid_fraction_static: 0
  });

var qcTimeSeriesImage = ee.Image(ee.Algorithms.If(
  timeSeriesDaily.size().gt(0),
  timeSeriesDaily.sort('valid_fraction_static', false).first(),
  emptyQcTsImage
));

function addRepresentativeTimeSeriesIndexQc(qcImage) {
  if (!DO_ADD_FIRST_TS_INDEX_QC) return;

  qcImage = ee.Image(qcImage);

  var indexVis = {
    min: -1,
    max: 1,
    palette: [
      '313695', '74add1', 'e0f3f8',
      'ffffff',
      'fee090', 'f46d43', 'a50026'
    ]
  };

  Map.addLayer(
    qcImage.select('ndvi'),
    indexVis,
    'QC TS highest coverage | NDVI',
    false
  );

  Map.addLayer(
    qcImage.select('savi'),
    indexVis,
    'QC TS highest coverage | SAVI',
    false
  );

  Map.addLayer(
    qcImage.select('gndvi'),
    indexVis,
    'QC TS highest coverage | GNDVI',
    false
  );

  Map.addLayer(
    qcImage.select('ndmi'),
    indexVis,
    'QC TS highest coverage | NDMI',
    false
  );

  Map.addLayer(
    qcImage.select('grvi'),
    indexVis,
    'QC TS highest coverage | GRVI',
    false
  );

  Map.addLayer(
    qcImage.select('valid').selfMask(),
    {min: 0, max: 1, palette: ['00ff00']},
    'QC TS highest coverage | valid',
    false
  );

  print('QC representative time-series metadata:', qcImage.toDictionary([
    'system:time_start',
    'date_key',
    'est_id',
    'season_tag',
    'year',
    'source_images_same_date',
    'valid_pixel_count',
    'valid_fraction_static'
  ]));

  if (FIRST_TS_QC_PRINT_STATS) {
    var stats = qcImage.select(TIME_SERIES_INDEX_BANDS).reduceRegion({
      reducer: ee.Reducer.minMax().combine({
        reducer2: ee.Reducer.mean(),
        sharedInputs: true
      }),
      geometry: roi,
      scale: EXPORT_SCALE,
      maxPixels: EXPORT_MAX_PIXELS,
      tileScale: 4
    });

    print('QC representative time-series index statistics:', stats);
  }
}

// This QC runs directly from the server-side collection. It does NOT need
// getInfo() or evaluate(), and it does not affect the export collection.
addRepresentativeTimeSeriesIndexQc(qcTimeSeriesImage);

// ----------------------------------------------------------------
// 8.4 Export orchestration after RUN01 preflight
// ----------------------------------------------------------------
function createTimeSeriesAssetExportTasks() {
  if (!DO_EXPORT_ASSET) return;

  // Index tasks are intentionally independent from land/water tasks. Thus,
  // even if a future land/water diagnostic fails, RUN02 input exports remain
  // available and the existing pipeline is not blocked.
  timeSeriesDateKeysServer.evaluate(function(dateKeysClient, failure) {
    if (failure !== null && failure !== undefined) {
      print('ERROR retrieving INDEX date keys; no index TS Tasks were created:', failure);
      return;
    }

    if (!Array.isArray(dateKeysClient)) {
      print('ERROR: index date-key request returned an unexpected value. No index TS Tasks were created:',
        dateKeysClient);
      return;
    }

    try {
      createTimeSeriesExportTasks(dateKeysClient);
    } catch (callbackError) {
      print('ERROR while creating index TS Tasks:', String(callbackError));
    }
  });

  // This list is derived before expensive per-pixel Otsu work. It should be a
  // lightweight metadata request and cannot block the index-series callback.
  landWaterDateKeysServer.evaluate(function(dateKeysClient, failure) {
    if (failure !== null && failure !== undefined) {
      print('ERROR retrieving LAND/WATER date keys; no land/water Tasks were created:', failure);
      return;
    }

    if (!Array.isArray(dateKeysClient)) {
      print('ERROR: land/water date-key request returned an unexpected value. No land/water Tasks were created:',
        dateKeysClient);
      return;
    }

    try {
      createLandWaterTimeSeriesExportTasks(dateKeysClient);
    } catch (callbackError) {
      print('ERROR while creating land/water Tasks:', String(callbackError));
    }
  });
}

function createDriveExportTasks() {
  if (!DO_EXPORT_DRIVE) return;

  Export.image.toDrive(makeDriveExportArgs(
    compositeExport, 'RUN01V3_COMPOSITE_' + prefix, roi
  ));
  Export.image.toDrive(makeDriveExportArgs(
    pcaExport, 'RUN01V3_PCA3_' + prefix, roi
  ));
  Export.image.toDrive(makeDriveExportArgs(
    ndviExport, 'RUN01V3_NDVI_' + prefix, roi
  ));
  Export.image.toDrive(makeDriveExportArgs(
    ndviClassExport, 'RUN01V3_NDVI_CLASS_' + prefix, roi
  ));
  Export.image.toDrive(makeDriveExportArgs(
    seasonalQa, 'RUN01V3_QA_' + prefix, roi
  ));

  print('RUN01 Drive export Tasks created for:', prefix);
}

function createRun01ExportTasks() {
  if (RUN01_MODE !== 'QA_AND_EXPORT') {
    print('RUN01 is QA_ONLY. No RUN01 Tasks created.');
    return;
  }

  if (!DO_EXPORT_ASSET && !DO_EXPORT_DRIVE) {
    print('RUN01 QA passed, but DO_EXPORT_ASSET and DO_EXPORT_DRIVE are both false. No Tasks created.');
    return;
  }

  createSeasonalAssetExportTasks();
  createTimeSeriesAssetExportTasks();
  createDriveExportTasks();

  print('RUN01 export mode:', RUN01_MODE);
  print('RUN01 export switches:', {
    asset: DO_EXPORT_ASSET,
    drive: DO_EXPORT_DRIVE,
    rgb_qa: EXPORT_RGB_QA,
    max_index_ts_exports: MAX_TS_EXPORTS_PER_RUN,
    max_landwater_ts_exports: MAX_LAND_WATER_TS_EXPORTS_PER_RUN
  });
}

// Run01 QA checks structural readiness only. Otsu thresholds/buffer values
// are intentionally not treated as fixed scientific truth; evaluate them in
// the diagnostic layers before moving from QA_ONLY to QA_AND_EXPORT.
var seasonalValidFractionStatic = ee.Number(ee.Algorithms.If(
  staticAllowedPixelCount.gt(0),
  seasonalValidCount.divide(staticAllowedPixelCount),
  0
));

// Check the raw threshold before clamping. This flags scenes where Otsu is
// being strongly constrained by the configured range even though a clamped
// threshold remains available for safe execution.
var compositeOtsuWithinClamp = compositeOtsuThresholdRaw
  .gte(OTSU_CLAMP_MIN)
  .and(compositeOtsuThresholdRaw.lte(OTSU_CLAMP_MAX));

var run01QaConditions = staticAllowedPixelCount.gt(0)
  .and(s2WithCloudProbability.size().gte(QA_MIN_CLOUD_MATCHED_IMAGES))
  .and(timeSeriesDaily.size().gte(QA_MIN_VALID_DAILY_DATES))
  .and(seasonalValidCount.gte(QA_MIN_SEASONAL_VALID_PIXELS));

run01QaConditions = ee.Number(ee.Algorithms.If(
  QA_REQUIRE_COMPOSITE_OTSU_WITHIN_CLAMP,
  run01QaConditions.and(compositeOtsuWithinClamp),
  run01QaConditions
));

var run01Preflight = ee.Dictionary({
  run01_mode: RUN01_MODE,
  est_id: EST_ID_QUERY,
  year: YEAR,
  season_tag: SEASON_TAG,
  season_start: START_DATE,
  season_end_exclusive: END_DATE,

  static_allowed_pixels: staticAllowedPixelCount,
  n_source_images: s2Base.size(),
  n_cloud_probability_matched_images: s2WithCloudProbability.size(),
  n_processed_images: s2Processed.size(),
  n_valid_daily_dates: timeSeriesDaily.size(),
  n_landwater_source_dates: landWaterDateKeysServer.size(),
  seasonal_valid_pixels: seasonalValidCount,
  seasonal_valid_fraction_static: seasonalValidFractionStatic,

  use_otsu_water_threshold: USE_OTSU_WATER_THRESHOLD ? 1 : 0,
  otsu_water_index: OTSU_WATER_INDEX,
  otsu_water_band: OTSU_WATER_INDEX_BAND,
  otsu_water_index_label: OTSU_WATER_INDEX_LABEL,
  otsu_fallback_threshold: WATER_FALLBACK_THRESHOLD,
  otsu_clamp_min: OTSU_CLAMP_MIN,
  otsu_clamp_max: OTSU_CLAMP_MAX,
  composite_otsu_threshold_raw: compositeOtsuThresholdRaw,
  composite_otsu_threshold: compositeOtsuThreshold,
  composite_otsu_used: compositeOtsuUsed,
  otsu_land_buffer: OTSU_LAND_BUFFER,
  composite_final_land_threshold: compositeFinalThreshold,
  composite_otsu_within_clamp: compositeOtsuWithinClamp,
  pixels_retained_only_by_buffer: compositeBufferRetainedLandCount,

  qa_min_cloud_matched_images: QA_MIN_CLOUD_MATCHED_IMAGES,
  qa_min_valid_daily_dates: QA_MIN_VALID_DAILY_DATES,
  qa_min_seasonal_valid_pixels: QA_MIN_SEASONAL_VALID_PIXELS,
  qa_require_composite_otsu_within_clamp:
    QA_REQUIRE_COMPOSITE_OTSU_WITHIN_CLAMP ? 1 : 0,
  qa_pass: run01QaConditions
});

// Single lightweight callback. Only after qa_pass=1 can QA_AND_EXPORT create
// any tasks. QA_ONLY is safe for testing OTSU_WATER_INDEX and OTSU_LAND_BUFFER repeatedly.
run01Preflight.evaluate(function(preflightClient, failure) {
  if (failure !== null && failure !== undefined) {
    print('RUN01 preflight evaluation failed. No Tasks were created:', failure);
    return;
  }

  print('RUN01 preflight summary:', preflightClient);

  if (!preflightClient || Number(preflightClient.qa_pass) !== 1) {
    print('RUN01 QA FAILED. No Tasks were created.');
    print('Inspect scene/date counts, seasonal valid domain and selected Otsu-index diagnostics.');
    return;
  }

  print('RUN01 QA PASSED.');
  createRun01ExportTasks();
});

// ================================================================
// END RUN01_V3.5 — selectable Otsu water index + QA modes
// ================================================================

