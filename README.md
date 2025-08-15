# projectionpasta
A python script for reprojecting raster maps

Created 2023 by Mads de Silva and Nikolai Hersfeldt

As part of the [Worldbuilding Pasta](https://worldbuildingpasta.blogspot.com/) project

Hardly the most detailed or user-friendly map projection software out there, but comes with a couple key advantages:
- Both the input and output maps can use any of the supported projections and any arbitrary aspect (orientation of the map relative to the globe's coordinates)
- The color mode of the original map is always preserved (presuming it's one that Pillow supports)

It is available in 2 forms: First, a standalone .exe that can be run on Windows without any dependencies but, due to filesize restrictions, lacks interpolation;

Second, the raw python script, which depends on the following packages that should be available through pip:
- numpy
- Pillow
- scipy   (if using interpolation; won't be imported otherwise)

Either will run as a command-line program, and the latter can potentially be imported and used as a python function. I've only tested these on Windows, but there's no particular reason to think the python version shouldn't work on other operating systems.

The following projections are currently supported:
- Equirectangular / Plate Caree
- Mercator
- Gall Stereographic
- Miller Cylindrical
- Cylindrical Equal Area
- Sinusoidal
- Mollweide
- Hammer
- Eckert IV
- Equal Earth
- Winkel Tripel
- Robinson
- Wagner VI
- Kavrayskiy VII
- Natural Earth
- Aitoff
- Azimuthal Equidistant
- Lambert Azimuthal Equal-Area
- Stereographic
- Orthographic
- Equidistant Conic
- Albers Equal-Area Conic
- Lambert Conformal Conic
- Nicolosi Globular
- Ortelius Oval
- Pseudostereographic
- Pseudoorthographic

If there's a projection you really need for something you can raise an issue and I can look into it (check the list of planned projections below).

## Basic Use

On running the script or .exe, a command-line window will open and you'll be shown a list of supported projections and prompted for information on the input and output maps. Filenames should include extensions like .png, and projections should be selected by their index code as shown in the list. You'll then be prompted for the map aspect, defined by the longitude and latitude of the map center and a third rotation from north; note that this changes the globe orientation of the map relative to the map center rather than the shape of the resulting image (probably best to just mess around with it to see what it does). These should all be in degrees, inputting 0 for all just projects the map without any change in aspect.

A few projections may prompt you with additional options:
- Azimuthal and conic projections can be projected as single global maps (where possible), a single hemisphere, or bihemispheric maps with maps of each hemisphere placed side-by-side
- Conic projections and Cylindrical Equal Area require reference latitudes of best shape.
- Mercator and Stereographic and Lambert Conformal Conic in global mode require maximum latitudes at which maps will be truncated.

The script should then run, informing you of the progress, and save the resulting map image when complete. The output configuration of one map (projection, aspect, etc.) can then be used as the input configuration for projecting it again.

Input maps should be cropped to the exact edges of the map, with no margins. Where input and output projections have the same aspect ratio, the output will have the same dimensions as the input; where aspect ratios vary, the output's dimensions will always be increased from the input's, never decreased. The output map is saved with the same color mode as the input map, and thus e.g. will not reduce the level of detail for 16-bit heightmaps. Non-visible areas of the output will be filled in with 0, which usually comes out as black, but may vary if using a color palette. By default, pixels are projected with a nearest-neighbor approach, which copies pixels directly from input map to output map without any blurring or interpolation.

## Advanced Use

Additional configuration options are available in the def_opts listed at the front of the script or the external projpasta_options.cfg file (which also works for the .exe version). The comments within explain all the options, but some highlights include:
- Options to control the output aspect ratio or image size
- Options to control truncation of the visible map area
- Cropping of visible area or output image size based on pixel indices or lon/lat coordinates
- Graticules (still semi-experimental, seems to work fine within the map area but has some edge issues)
- Alternate interpolation types (python script version only)
- Options to set the input and output configuration and skip the command-line input.

Functions in the python script can also be run directly. 3 main functions are intended for easy use:
- Proj_Image takes an image file and produces a projected image file; essentially what running from command line does, but accessible as a function, with these parameters:
  - file_in: input file path string
  - file_out: output file path string, 'projp_out.png' by default
  - proj_in: input projection name string or index int, equirectangular by default
  - proj_out: as above for output
  - aspect_in: input projection aspect, tuple or list of (center lon, center lat, 3rd rotation), (0,0,0) by default
  - aspect_out: as above for output
  - opts: a dictionary of any custom options (as seen in the def_opts list), any unspecified will be taken from def_opts or projpasta_options.cfg if present.
  - deg: treat aspects as degrees, True by default
  - get_index: function returns a projection index, which can be used to project other images with all the same configuration without requiring new calculations, False by default
  - index: an index to be used from a previous run of this function, as returned by get_index, None by default
- Proj_Array takes any 2- or 3-dimensional numpy array (with dimensions 0 and 1 being y and x, respectively) and returns a projected array, with these parameters:
  - proj_in, proj_out, aspect_in, aspect_out, opts, deg, get_index, and index as described above
  - data_in: the input numpy array
- Find_index takes any collection of x and y coordinates for an input map (where both range -1 to 1, left to right and top to bottom) and returns projected coordinates, with these parameters:
  - opts and deg as above
  - proj1, proj2, aspect1, and aspect2 like above, projecting from 1 to 2, though note that usual "backwards" projection actuall projects from output map to input map (to find the corresponding input position for each output pixel)
  - x1, y1: the x and y inputs, with can be single numbers, lists, or arrays, but must be ordered such that each pair of values describes one point on the map
  - get_lon: in addition to x and y outputs, return lon and lat coordinates for these points for, in order, map 1 treated as normal aspect, true coordinates based on aspect1, and map 2 treated as normal aspect

There's also an addition Globe_gif function to produce a rotating globe .gif image using the Orthographic projaction, with these parameters:
  - file_in, file_out, proj_in, opts, deg as above
  - init_aspect for the initial output aspect (treating the input as having aspect (0,0,0)), with only center longitude changed in subsequent frames to simulate globe rotation, (0,0,0) by default
  - frames: number of frames per total rotation, 36 by default
  - duration: length for complete rotation in milliseconds, 100 by default
  - loop: number of times to loop, 0 by default for infinite looping

## To Do

Outstanding issues:
- Winkel Tripel still has some small problem areas at the corners if untruncated
- Nicolosi Globular still lacks global mode, had some promising progress but nothing usable so far
- Experimental "forward" projection generally unreliable, needs work, should at least have some cropping implemented but I'm not sure I can do graticules
- Graticules can be marked along the edges of some projections in undesireable ways

 Future plans (no guarantees):
 - Map gores
 - tkinter GUI
 - output graticules only
 - anything that might speed up calculations a bit
 - New projections:
   - Peirce Quincuncial
   - 2-point equidistant
   - Chamberlin Trimetric
   - Polygonals maybe? Could be a real pain but no one will know I'm a real map nerd if I can't make a Cahill butterfly
   - HealPIX if I can be bothered
   - Maybe I should just have a general option for slicing maps together at a given latitude?

## Sources

- Most projection formulae taken from Wikipedia
- Some additional formulae (particularly conic inverses) from ["Map Projections - A Working Manual"](https://kartoweb.itc.nl/geometrics/Publications/Map%20Projections%20-%20A%20Working%20manual%20-%20by%20J.P.%20Snyder.pdf)  by John P. Snyder, 1987
- General iterative inverse projection method used for Equal Earth, Winkel Tripel, Natural Earth, Nicolosi Globular, and Ortelius Oval from Bildirici 2016, ["An iterative approach for inverse transformation of map projections"](https://www.tandfonline.com/doi/full/10.1080/15230406.2016.1200492)
- Multiquadric interpolation for Robinson from Ipbuker 2013, ["Numerical Evaluation of the Robinson Projection"](https://www.tandfonline.com/doi/abs/10.1559/1523040041649425)
