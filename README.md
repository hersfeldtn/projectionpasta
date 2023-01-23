# projectionpasta
A python script for reprojecting raster maps
Created 2023 by Amadea de Silva and Nikolai Hersfeldt
As part of the [Worldbuilding Pasta](https://worldbuildingpasta.blogspot.com/) project

Hardly the most detailed or user-friendly map projection software out there, but comes with a couple key advantages:
- Both the input and output maps can use any of the supported projections and any arbitrary aspect (orientation of the map relative to the globe's coordinates)
- The color mode of the original map is always preserved (presuming it's one that Pillow supports)

The script depends on the following packages, which should be available through pip:
- numpy
- Pillow

It can then be run as a command-line program, or potentially imported and used as a python function.

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
- Azimuthal Equidistant
- Orthographic
- Stereographic
- Mercator (truncated at 85.05 latitude to form a square)
- Gall Stereographic
- Miller Cylindrical

I'm working on Lambert Azimuthal Equal-Area and Nicolosi Globular but haven't worked out the bugs yet, and I'll probably throw in Cylindrical Equal-Area at some point

The azimuthal projections (Azimuthal Equidistant, Orthographic, Stereographic) can be projected as single global maps (where possible), a single hemispheric map, or a bihemispheric map with maps of each hemisphere placed side-by-side.

Both input and output map can have any arbitrary aspect, which is defined by specifying the latitude and longitude of the map's center and then the clockwise rotation of the globe around that center point, relative to the default of having the north pole directly above the center. Projecting out to a given output aspect and then taking that map as an input with that same information for its aspect should always allow you to return to the original aspect.

Output maps can be cropped or uncropped (though the distinction only matters with non-rectangular projections): Cropping removes areas not usually shown, such that each point of the globe only appears once on the map, but taking a cropped map as input and projecting again may sometimes missing points at the edges; this can be avoided with uncropped maps, but for some output projections (those using the iterative method described below) inaccuracies and unusual noise artifacts will appear in the areas usually cropped (though these shouldn't impact future reprojections from that map).

Data is copied from the input map to the output by a nearest-neighbor approach; for each pixel of the output map, the script finds the corresponding position on the input map, finds the nearest pixel, and directly copies it to the output map, without any interpolation or averaging function. The output map is saved with the same color mode as the input map.

Where input and output projections have the same aspect ratio, the output will have the same dimensions as the input; where aspect ratios vary, the output's dimensions will be increased from the input's.

Several of the supported projections (Aitoff, Winkel Tripel, Ortelius Oval) are defined by functions that specify map position from given latitude, longitude coordinates, which can be used for input maps here (because the script works backwards from the output map to the input), but no closed-form solution exists for determining lat, lon from map position. For these cases, use as an output map is allowed using the methodology from ["An iterative approach for inverse transformation of map projections"](https://doi.org/10.1080/15230406.2016.1200492), Bildirici 2016. This determines lat, lon from map position with an iterative approach that is repeated until either the error in every pixel is less than 1 millionth of a radian (excluding pixels inside the usually cropped area which never seem to converge properly, hence the noise that appears there) or a maximum number of iterations has passed (20 by default), in which case a warning will appear.

Mollweide has no closed-form solution to determine map position from lat, lon, and when used as input has its own iterative method with the same tolerances.
