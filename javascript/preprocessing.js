/**
 * @file This script contains functions to preprocess GEE image data
 * @author Zhiqi Yu <yzq1027@hotmail.com>
 * @license MIT 
 */


/**
 * This function combines MODIS Terra and Aqua data to generate one dataset. 
 * Originally written for mitigating missing value problem.
 * 
 * @param {ee.Image} terra - multi-band Terra image, number of bands can be 1 to n.
 * @param {ee.Image} aqua - multi-band Aqua image, number of bands should match terra bands.
 */
function combineModis(terra, aqua) {
    var bandwise_and = terra.and(aqua);

    /**
     * if band-wise-and returns 1, then two bands all have values, take the mean,
     * if not, then at least one band has missing value, add them up as the final value
     */
    var multiplicative = bandwise_and.add(1).pow(-1);   // from [0, 1] -> [1, 2] -> [1, 1/2] 
    var combined = terra.add(aqua).multiply(multiplicative);
    return combined;
}

