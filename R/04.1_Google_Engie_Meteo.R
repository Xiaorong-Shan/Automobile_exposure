/*******************************
 * https://code.earthengine.google.com/
 * NLDAS-2 FORA0125 (hourly) annual mean wind (u/v)
 * Export: state x year -> 1 GeoTIFF with 2 bands (Wind_E, Wind_N)
 *******************************/

// ===== CONFIG =====
var startYear = 2011;
var endYear   = 2020;

var stateAbbrs = ["CA","TX","VA","MD","NJ","IL","NY","DC","MA"];

// OPTIONAL: limit for testing (set null to run all)
var onlyTheseStates = ["MA"];  // e.g., ["CA"] or null
var onlyTheseYears  = null;  // e.g., [2017] or null

var outFolder = "NLDAS_Annual_Wind_State";

// NLDAS-2 FORA0125 hourly
var ic = ee.ImageCollection("NASA/NLDAS/FORA0125_H002");

// IMPORTANT: actual wind bands in this collection
var windBandsIn = ["wind_u", "wind_v"];     // u (east-west), v (north-south)
var windBandsOut = ["Wind_E", "Wind_N"];    // your preferred names

// NLDAS native ~14km
var exportScale = 14000;

// TIGER state boundaries (includes DC), field: STUSPS
var states = ee.FeatureCollection("TIGER/2018/States");

// ===== HELPERS =====
function getStateRegion(abbr) {
  return states
    .filter(ee.Filter.eq("STUSPS", abbr))
    .geometry()
    .buffer(1000);
}

function annualMeanWindUV(year, regionGeom) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end   = start.advance(1, "year");

  // Annual mean of hourly u/v, then rename to Wind_E/Wind_N
  return ic
    .filterDate(start, end)
    .select(windBandsIn)
    .mean()
    .rename(windBandsOut)
    .clip(regionGeom)
    .set({year: year});
}

// ===== MAIN LOOP =====
var years = [];
for (var y = startYear; y <= endYear; y++) years.push(y);

stateAbbrs.forEach(function(abbr) {
  if (onlyTheseStates && onlyTheseStates.indexOf(abbr) === -1) return;

  var region = getStateRegion(abbr);

  years.forEach(function(yr) {
    if (onlyTheseYears && onlyTheseYears.indexOf(yr) === -1) return;

    var img = annualMeanWindUV(yr, region);

    var desc = "NLDAS_" + abbr + "_" + yr + "_Wind_EN_mean";

    Export.image.toDrive({
      image: img,                 // 2 bands in one file: Wind_E, Wind_N
      description: desc,
      folder: outFolder,
      fileNamePrefix: desc,
      region: region,
      scale: exportScale,
      maxPixels: 1e13
      // don't set crs; let EE keep native projection
    });
  });
});

print("Export tasks created. Run them in the Tasks tab.");
