# projectionpasta
A python script for reprojecting raster maps

Created 2023 by Amadea de Silva and Nikolai Hersfeldt

As part of the [Worldbuilding Pasta](https://worldbuildingpasta.blogspot.com/) project

Hardly the most detailed or user-friendly map projection software out there, but comes with a couple key advantages:
- Both the input and output maps can use any of the supported projections and any arbitrary aspect (orientation of the map relative to the globe's coordinates)
- The color mode of the original map is always preserved (presuming it's one that Pillow supports)

It is available in 2 forms: First, a standalone .exe that can be run on Windows without any dependencies but, due to filesize restrictions, lacks interpolation;

Second, the raw python script, which depends on the following packages that should be available through pip:
- numpy
- Pillow
- scipy   (if using interpolation; won't be imported otherwise)

Either will run as a command-line program, and the latter can potentially be imported and used as a python function.

The following projections are currently supported:
- Equirectangular / Plate Caree
- Sinusoidal
- Mollweide
- Hammer
- Aitoff
- Winkel Tripel
- Kavrayskiy VII
- Wagner VI
- Ortelius Oval
- Nicolosi Globular
- Eckert IV
- Azimuthal Equidistant
- Orthographic
- Stereographic
- Lambert Azimuthal Equal-Area
- Mercator (truncated at 85.05 latitude to form a square)
- Gall Stereographic
- Miller Cylindrical

If there's a projection you really need for something you can raise an issue and I can look into it (Robinson would be nice but is a right pain to implement, Cylindrical equal-area I can add if someone really wants it, fancier stuff like conic and polyhedral projections are probably possible but would be a fair bit of work, and one of these days I might get around to adding the global version of Nicolosi)

The azimuthal projections (Nicolosi, Azimuthal Equidistant, Orthographic, Stereographic, Lambert) can be projected as single global maps (where possible), a single hemispheric map, or a bihemispheric map with maps of each hemisphere placed side-by-side.

Both input and output map can have any arbitrary aspect, which is defined by specifying the latitude and longitude of the map's center and then the clockwise rotation of the globe around that center point, relative to the default of having the north pole directly above the center. Projecting out to a given output aspect and then taking that map as an input with that same information for its aspect should always allow you to return to the original aspect.

By default, each pixel is copied from the input map to the output by a nearest-neighbor approach; for each pixel of the output map, the script finds the corresponding position on the input map, finds the nearest pixel, and directly copies it to the output map. As of v1.4 there are now also interpolation options from scipy:
- Nearest neighbor is in principle the same as default behavior, but probably a bit more optimized and reliable; directly copies pixels without any blurring or averaging, so best for preserving colors with hard edges.
- Linear is generally the fastest and most reliable interpolation, but can make the map pretty blurry over multiple steps of reprojection
- Cubic should be more precise but is very slow to run on large maps and can introduce odd pixel artifacts on complex maps, particularly with color.
- PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) is even slower but less prone to artifacting.
- Spline is reasonably fast but perhaps the most prone to introduce artifacts on color maps. For the moment, spline interpolation with a rectangular input map (equirectangular, mercator, Gall, Miller) is the only case where the interpolation treats the input map properly as a spherical surface rather than a flat map, and so best for avoiding seams at the map edges.

The output map is saved with the same color mode as the input map. Generally interpolation works best with greyscale and probably won't work sensibly at all if the input map uses a color palette.

Input maps should be cropped to the exact edges of the map. Where input and output projections have the same aspect ratio, the output will have the same dimensions as the input; where aspect ratios vary, the output's dimensions will always be increased from the input's, never decreased. Not that these dimensions are determined assuming the input has a correct aspect ratio, so e.g. if you put in an equirectangular that doesn't have a 2:1 aspect ratio, the output will similarly distorted dimensions.

You can also now choose to trim the output map dimensions; in this case you'll be prompted to specify the pixel indices for the edges of the output map, and then all further projection and interpolation steps will be applied only to that section, which can save time on large maps you only want a section of.

Output maps can also be cropped or uncropped (though the distinction only matters with non-rectangular projections): Cropping removes areas not usually shown, such that each point of the globe only appears once on the map, but taking a cropped map as input and projecting again may sometimes leave missing points at the edges; this can be avoided with uncropped maps, but for some output projections (those using the iterative method described below) inaccuracies and unusual noise artifacts will appear in the areas usually cropped (though these usually shouldn't impact future reprojections from that map).

Several of the supported projections (Aitoff, Winkel Tripel, Ortelius Oval, Nicolosi) are defined by functions that specify map position from given latitude, longitude coordinates, which can be used for input maps here (because the script works backwards from the output map to the input), but no closed-form solution exists for determining lat, lon from map position. For these cases, use as an output map is allowed using the methodology from ["An iterative approach for inverse transformation of map projections"](https://doi.org/10.1080/15230406.2016.1200492), Bildirici 2016. This determines lat, lon from map position with an iterative approach that is repeated until either the error in every pixel is less than 1 millionth of a radian (excluding pixels inside the usually cropped area which never seem to converge properly for some of these projections, hence the noise that appears there) or a maximum number of iterations has passed (20 by default), in which case a warning will appear.

Mollweide and Eckert IV have no closed-form solutions to determine map position from lat, lon, and when used as input have their own iterative methods with the same tolerances.
