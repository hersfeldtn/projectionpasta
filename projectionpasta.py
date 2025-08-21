from PIL import Image
import numpy as np
import os
import configparser
import math as ma

ver_num = "2.0.1"

#Increases potential map size to 1 million x 1 million pixels, which is probably more than anyone's RAM could actually handle
Image.MAX_IMAGE_PIXELS = 1e12

#Dictionary of default options
def_opts = {

    'alt_config': None,                 #Alternative config file to load options from

    #Projection-specific options

    'azim_type_in': 'global_if',        #'global' : project azimuthal or conic maps as single global maps
                                        #   'hem' : project as single hemispheres
                                        #   'bihem' : project as paired hemispheres
                                        #   'global_if' : project as global if possible, bihem otherwise
    'azim_type_out': 'global_if',       #as above for output
    'azim_type': None,                  #azim type for both input and output; overwrites both above if not None
    'ref_in': None,                     #reference values for input proojection (e.g. maximum lat, reference lats)
    'ref_out': None,                    #as above for output
    'ref': None,                        #reference values for both input and output, overwrites both above if not None

    #Output options

    'force_ratio': None,                #Force output map width/height aspect ratio (overwriting per-projection ratios; not relative to input)
    'relative_ratio': False,            #scale output map by relative ratios of input and output projections; if false, output ratio determined purely by output projection
    'force_scale': None,                #Force output scale; if single value, will be taken as width and height will be scaled to appropriate ratio; if pair, will be taken as (width,height)
        
    'truncate_in': True,                #treat input map as truncated to a single world surface or hemisphere
    'truncate_out': True,               #truncate output map to a single world surface or hemisphere
    
    'crop_in': None,                    #crop input map to (left, right, top, bottom) pixel indices of given input map, limiting visible area on output
    'crop_out': None,                   #crop output map to (left, right, top, bottom) pixel indices of typical global output map
    'is_crop_in': None,                 #treat input map as having been cropped to these indices; (left, right, top, bottom, full map width, full map height), the last 2 referring to the size of the original input map
                                        #   (i.e. use to reproject map after using crop_out or crop_vis, with the first 4 indices being the same)
    
    'crop_coords_true': None,           #limit output map to within (min lon, max lon, min lat, max lat) based on given aspect information
    'crop_coords_in': None,             #limit output map to within (min lon, max lon, min lat, max lat) on input map, ignoring given aspect and treating as normal aspect (with 0,0 center)
    'crop_coords_out': None,            #limit output map to within (min lon, max lon, min lat, max lat) on output map, treating as normal aspect
                                        #   all crop_coords just hide areas and don't change output image size

    'crop_vis': False,                  #crop output to smallest box containing visible areas after all other cropping applied, and report indices used (for later use of is_crop_in)
                                        #   otherwise only crop_out actually crops map rather than just hiding areas

    'graticules': None,                 #Add graticules
                                        #   'true' : mark true lon/lat based on output aspect
                                        #   'in': mark lon/lat treating input map as normal aspect, ignoring the given input aspect
                                        #   'out': mark lon/lat on output map, treating it as normal aspect
    'grat_lon': 30,                     #list of longitudes to mark with graticules, or single interval between graticules (starting at antimeridian and proceding east), in degrees
    'grat_lat': 30,                     #as above for latitudes (with interval starting at south pole and proceding north)
    'grat_color':(100,100,100,50),      #graticule color in RGBA

    #Procedure configuration

    'interp_type': 'nearest',           #type of interpolation to apply
                                        #   'none' : no interpolation and scipy not required; backward projection only
                                        #   backward projection supports all scipy interpn methods ('nearest','linear','slinear','cubic''quintic','pchip','spline2fd')
                                        #   forward projection supports only scipy griddata methods ('nearest','linear','cubic')
    'proj_direction': 'backward',       #'backward': always project backward (interpolate on input map)
                                        #   'forward': always project forward (interpolate on output map)
                                        #   'avoid_iter': backward by default, but avoid iterative methods if possible
    'use_sym': True,                    #use projection symmetry to speed calculations; appropriate only for global maps
    'avoid_seam': True,                 #use extra steps to avoid seams at edge of original map
    'tolerance': 1e-6,                  #tolerated maximum proportional error on iterative methods
    'max_iter': 20,                     #maximum iterations for iterative methods

    #Internal use; should only be used if accessing internal functions directly
    
    'in': True,                         #trigger to indicate if particular function is being used for input or output map (for selecting azim type and ref)
    'global_in': True,                  #trigger for if input map should be treated as complete globe surface
    'skip_update': False,               #trigger to stop updating opts from file once it's already been done
    'skip_setup': False,                #trigger to skip setup, used only by configs

}



#Remapping functions for each projection:
#coords determines (lon,lat) coordinates for a given (x,y) position on the map
#pos determines (x,y) position for given (lon,lat) coordinates
#vis determines if an (x,y) position on the output should be shown on non-rectangular maps
#pres determines if an (x,y) position is present on the input
#typ incorporates extra information required for specific projection types
#rat is the map width/height ratio
#int determines type of interpolation
#x and y have ranges of (-1, 1)
#lon is (-pi,pi) or (-pi/2, pi/2) as appropriate
#lat is (-pi/2, pi/2)
#with (0,0) at the map center

#First, some generic functions used for multiple projections

def From_pol(rh,th):   #Converts from polar around equator to lon,lat
    rh *= ma.pi
    lat = np.arcsin(-np.sin(rh)*np.cos(th))
    lon = np.arctan2(np.sin(th)*np.sin(rh),np.cos(rh))
    return lon,lat
def To_pol(lon,lat):    #Converts from lon,lat to polar around equator
    rh = np.arccos(np.cos(lat)*np.cos(lon))/ma.pi
    th = np.arctan2(np.sin(lon)*np.cos(lat),np.sin(lat))
    return rh,th
def Get_ref(opts=def_opts, def_ref=None, deg=True):         #determine appropriate ref
    if opts['ref'] is not None:
        ref = opts['ref']
    elif opts['in']:
        ref = opts['ref_in']
    else:
        ref = opts['ref_out']
    if ref is None:
        ref = def_ref
    if deg:
        ref = np.radians(ref)
    return ref

def der(lon,lat,t,f,opts):   #Approximate partial derivatives by the secant method
    x1,y1 = f(lon,lat+t,opts)
    x2,y2 = f(lon,lat-t,opts)
    x3,y3 = f(lon+t,lat,opts)
    x4,y4 = f(lon-t,lat,opts)
    dxla = (x1-x2)/(2*t)
    dyla = (y1-y2)/(2*t)
    dxlo = (x3-x4)/(2*t)
    dylo = (y3-y4)/(2*t)
    return dxla,dxlo,dyla,dylo
def Inv_coords(x,y,fpos,fcoords,vis,opts=def_opts): #General iterative inverse function from Bildirici 2016: https://doi.org/10.1080/15230406.2016.1200492
    print("  (using iterative method, may take a bit)")
    t = opts['tolerance']
    imax = opts['max_iter']
    lon,lat = fcoords(x,y,opts)
    #lon = np.fmod(lon + ma.pi, 2*ma.pi) - ma.pi
    vis1 = vis(x,y,None,None,opts)
    i=1
    while True:
        x1,y1 = fpos(lon,lat,opts)
        dxla,dxlo,dyla,dylo = der(lon,lat,t,fpos,opts)
        div = dxla*dylo - dyla*dxlo
        divz=np.where(div!=0,True,False)
        div1=np.where(divz,div,1)
        difx = x1-x
        dify = y1-y
        dlon = np.where(divz,(dify*dxla - difx*dyla)/div1,0)
        dlat = np.where(divz,(difx*dylo - dify*dxlo)/div1,0)
        if np.amax(np.where(vis1,np.abs(dlon),0)) < t:      #only check within the region that'll typically be visible, there are some glitchy regions outside I haven't managed to eliminate
            if np.amax(np.where(vis1,np.abs(dlat),0)) < t:
                break
        lon = lon-dlon
        lat = lat-dlat
        i += 1
        if i>imax:
            print("  Reached maximum of "+str(imax)+" iterations without converging, outputting result")
            break
    return lon,lat



def Def_vis(x,y,lat=None,lon=None,opts=None):     #all map visible
    return True
def Coords_vis(x,y,lon,lat,opts=None):  #mark map edge by coordinates
    return np.where(np.logical_or(np.abs(lon) > ma.pi, np.abs(lat) > ma.pi/2), False, True)
def Edge_vis(x,y=None,lat=None,lon=None, f_coords=None,f_pos=None, opts=def_opts):
    if lat is None:
        lon, lat = f_coords(x,y,opts)
        return np.where(np.abs(lon) > ma.pi, False, True)
    xmax, y1 = f_pos(ma.pi,np.where(np.abs(lat)>ma.pi/2,ma.pi/2,lat),opts)
    return np.where(np.abs(x) > xmax, False, True)
def Circ_vis(x,y,lat=None,lon=None,opts=None):        #for circular or elliptical maps
    return np.where(np.sqrt(x**2 + y**2) > 1, False, True)
def Pill_vis(x,y,lat=None,lon=None,opts=None):        #2 semicircles flanking a square (at 2:1 aspect ratio)
    return np.where(np.abs(x)>1/2,Circ_vis(1-abs(x*2),y,opts),True)


#Lists and dictionaries of projection properties
sections = {}           #Section headers
proj_list = []          #projection names
info_list = {}          #extra info for each projection
coordsl = {}            #coords functions (x,y to lon,lat)
posl = {}               #pos functions  (lon,lat to x,y)
visl = {}               #vis functions  (portion of map visible)
ratl = {}               #aspect ratios
wrapl = {}              #wrap types (for map edge padding)
syml = {}               #symmetry types
aziml = []              #azimuthal projections
conicl = []             #conic projections
nonglobal = []          #projections that don't cover whole globe
iter_coords = []        #projections with iterative coords functions
iter_pos = []           #projections with iterative pos functions
optlists = {}           #dictionary of dictionaries of additional options for specific projections
optlists_global = {}    #option lists only shown for global version of azimuthal or conic projections

sections['Equirectangular'] = 'Cylindrical'


proj_list.append("Equirectangular")
info_list['Equirectangular'] = "2:1 rectangle (cylindrical); equidistant (vertical axis)"

def Equi_coords(x,y,opts=None):
    lon = ma.pi*x
    lat = ma.pi*y/2
    return lon,lat
def Equi_pos(lon,lat,opts=None):
    x = lon/ma.pi
    y = 2*lat/ma.pi
    return x,y
coordsl['Equirectangular'] = Equi_coords
posl['Equirectangular'] = Equi_pos
visl['Equirectangular'] = Def_vis
ratl['Equirectangular'] = 2
wrapl['Equirectangular'] = 'rect'
syml['Equirectangular'] = 'sym4lat'


proj_list.append('Mercator')
info_list['Mercator'] = 'variable rectangle (cylindrical); conformal; truncated poles'

def Merc_coords(x,y,opts=def_opts):
    ref = Get_ref(opts, Merc_def_ref)
    y1 = y * ma.log(ma.tan(ma.pi/4 + ref/2))/ma.pi
    lon = ma.pi*x
    lat = np.sign(y1)*(2*np.arctan(np.exp(abs(y1)*ma.pi))-ma.pi/2)
    return lon,lat
def Merc_pos(lon,lat,opts=def_opts):
    ref = Get_ref(opts, Merc_def_ref)
    x = lon/ma.pi
    y = np.sign(lat)*np.log(np.tan(ma.pi/4 + abs(lat)/2))/ma.pi
    y = y / ma.log(ma.tan(ma.pi/4 + ref/2))*ma.pi
    return x,y
def Merc_vis(x,y,lat=None,lon=None,opts=def_opts):
    y1 = (y * x) / x    #hack to make quickvis work because I can't be bothered making it more elegant
    return np.where(abs(y1) > 1, False, True)
def Merc_rat(opts=def_opts):
    ref = Get_ref(opts, Merc_def_ref)
    return ma.pi/ma.log(ma.tan(ma.pi/4 + ref/2))
Merc_def_ref = np.degrees(2*ma.atan(ma.exp(ma.pi))-ma.pi/2)
coordsl['Mercator'] = Merc_coords
posl['Mercator'] = Merc_pos
visl['Mercator'] = Merc_vis
ratl['Mercator'] = Merc_rat
wrapl['Mercator'] = 'rect'
syml['Mercator'] = 'sym4lat'
nonglobal.append('Mercator')
optlists['Mercator'] = {
    'Select maximum latitude to truncate map (0-90 degrees from equator)': 1,
    '85.05: forms a 1:1 square': Merc_def_ref
}


proj_list.append('Gall Stereographic')
info_list['Gall Stereographic'] = '1.301:1 rectangle (cylindrical)'

def Gallst_coords(x,y,opts=None):
    lon = x*ma.pi
    lat = 2*np.arctan(y)
    return lon,lat
def Gallst_pos(lon,lat,opts=None):
    x = lon/ma.pi
    y = np.tan(lat/2)
    return x,y
coordsl['Gall Stereographic'] = Gallst_coords
posl['Gall Stereographic'] = Gallst_pos
visl['Gall Stereographic'] = Def_vis
ratl['Gall Stereographic'] = (ma.pi/ma.sqrt(2))/(1+ma.sqrt(2)/2)
wrapl['Gall Stereographic'] = 'rect'
syml['Gall Stereographic'] = 'sym4lat'


proj_list.append('Miller Cylindrical')
info_list['Miller Cylindrical'] = '1.364:1 rectangle (cylindrical)'

def Mill_coords(x,y,opts=None):
    lon = x*ma.pi
    y1 = y*ma.asinh(ma.tan(ma.pi*2/5))
    lat = 5/4*np.arctan(np.sinh(y1))
    return lon,lat
def Mill_pos(lon,lat,opts=None):
    x = lon/ma.pi
    y1 = np.arcsinh(np.tan(lat*4/5))
    y = y1/ma.asinh(ma.tan(ma.pi*2/5))
    return x,y
coordsl['Miller Cylindrical'] = Mill_coords
posl['Miller Cylindrical'] = Mill_pos
visl['Miller Cylindrical'] = Def_vis
ratl['Miller Cylindrical'] = ma.pi/(5/4*ma.asinh(ma.tan(ma.pi*2/5)))
wrapl['Miller Cylindrical'] = 'rect'
syml['Miller Cylindrical'] = 'sym4lat'


proj_list.append('Cylindrical Equal-Area')
info_list['Cylindrical Equal-Area'] = 'variable rectangle (cylindrical); equal-area'

def Cyleq_coords(x,y,opts=None):
    lon = x * ma.pi
    lat = np.arcsin(y)
    return lon, lat
def Cyleq_pos(lon,lat,opts=None):
    x = lon / ma.pi
    y = np.sin(lat)
    return x, y
def Cyleq_rat(opts=def_opts):
    ref = Get_ref(opts, Cyleq_def_ref)
    return ma.pi * ma.cos(ref)**2
Cyleq_def_ref = 37.5
coordsl['Cylindrical Equal-Area'] = Cyleq_coords
posl['Cylindrical Equal-Area'] = Cyleq_pos
visl['Cylindrical Equal-Area'] = Def_vis
ratl['Cylindrical Equal-Area'] = Cyleq_rat
wrapl['Cylindrical Equal-Area'] = 'rect'
syml['Cylindrical Equal-Area'] = 'sym4lat'
optlists['Cylindrical Equal-Area'] = {
    'Select reference latitude of best shape (0-90 degrees from equator)': 1,
    '0: Lambert, aspect ratio 3.141:1': 0,
    '30: Behrmann, 2.356:1': 30,
    '37.071: Smyth/Craster, 2:1': ma.degrees(ma.acos(ma.sqrt(2/ma.pi))),
    '37.5: Hobo-Dyer, 1.977:1': 37.5,
    '45: Gall-Peters, 1.571:1': 45,
    '50: Balthasart, 1.298:1': 50,
    '55.654: Tobler, 1:1': ma.degrees(ma.acos(ma.sqrt(1/ma.pi)))
}


sections['Sinusoidal'] = 'Pseudocylindrical/azimuthal Equal-Area'


proj_list.append('Sinusoidal')
info_list['Sinusoidal'] = '2:1 sinusoid (pseudocylindrical); equal-area, equidistant (horizontal axis)'

def Sin_coords(x,y,opts=None):
    lat = ma.pi*y/2
    lon = ma.pi*x/np.cos(lat)
    return lon,lat
def Sin_pos(lon,lat,opts=None):
    x = np.cos(lat)*lon/ma.pi
    y = 2*lat/ma.pi
    return x,y
def Sin_vis(x,y,lat=None,lon=None,opts=None):
    #if lat is None:
    return np.where(abs(x) > np.cos(ma.pi*y/2), False, True)
    #return Coords_vis(x,y,lat,lon,opts)
coordsl['Sinusoidal'] = Sin_coords
posl['Sinusoidal'] = Sin_pos
visl['Sinusoidal'] = Sin_vis
ratl['Sinusoidal'] = 2
wrapl['Sinusoidal'] = 'xwrap'
syml['Sinusoidal'] = 'sym4lat'


proj_list.append('Mollweide')
info_list['Mollweide'] = '2:1 ellipse (pseudocylindrical); equal-area'

def Moll_coords(x,y,opts=None):
    th = np.arcsin(y)
    lat = np.arcsin((2*th + np.sin(2*th))/ma.pi)
    #lon = np.fmod(ma.pi*x/(np.cos(th)), ma.pi)
    lon = ma.pi*x/(np.cos(th))
    return lon,lat
def Moll_pos(lon,lat,opts=def_opts):
    print("  (using iterative method, may take a bit)")
    t = opts['tolerance']
    imax= opts['max_iter']
    i = 1
    th = lat
    err = 1
    while err > t:   #no closed-form solution, so iterate until maximum error is < tolerance
        th = th - (2*th + np.sin(2*th) - ma.pi*np.sin(lat)) / (2 + 2*np.cos(2*th))
        err = np.amax(np.abs((2*th+np.sin(2*th))/(ma.pi*np.sin(lat))-1))
        i+=1
        if i>imax:
            print("  Reached maximum of "+str(imax)+" iterations without converging, outputting result")
            print(f"   max remaining error: {err}")
            break
    x = lon*np.cos(th)/ma.pi
    y = np.sin(th)
    return x,y
coordsl['Mollweide'] = Moll_coords
posl['Mollweide'] = Moll_pos
visl['Mollweide'] = Circ_vis
ratl['Mollweide'] = 2
wrapl['Mollweide'] = 'xwrap'
syml['Mollweide'] = 'sym4lat'
iter_pos.append('Mollweide')


proj_list.append('Hammer')
info_list['Hammer'] = '2:1 ellipse (pseudoazimuthal); equal-area'

def Hammer_coords(x,y,opts=None):
    x1 = x*ma.sqrt(2)*2
    y1 = y*ma.sqrt(2)
    z = np.sqrt(1 - (x1/4)**2 - (y1/2)**2)
    lon = 2*np.arctan( z*x1 / (2*(2*z**2-1)) )
    lat = np.arcsin(z*y1)
    return lon,lat
def Hammer_pos(lon,lat,opts=None):
    x = np.cos(lat)*np.sin(lon/2) / np.sqrt(1+np.cos(lat)*np.cos(lon/2))
    y = np.sin(lat) / np.sqrt(1+np.cos(lat)*np.cos(lon/2))
    return x,y
coordsl['Hammer'] = Hammer_coords
posl['Hammer'] = Hammer_pos
visl['Hammer'] = Circ_vis
ratl['Hammer'] = 2
wrapl['Hammer'] = 'xwrap'
syml['Hammer'] = 'sym4'


proj_list.append('Eckert IV')
info_list['Eckert IV'] = '2:1 oval (pseudocylindrical); equal-area'

def Eckiv_coords(x,y,opts=None):
    r = 1 / (2 * ma.sqrt(ma.pi / (4 + ma.pi)))
    th = np.arcsin(y * np.sqrt(4+ma.pi) / (2*np.sqrt(ma.pi)*r))
    lat = np.arcsin((th + np.sin(th)*np.cos(th) + 2*np.sin(th)) / (2 + ma.pi/2))
    lon = 2 * x * np.sqrt(4*ma.pi + ma.pi**2) / (2 * r * (1 + np.cos(th)))
    return lon,lat
def Eckiv_pos(lon,lat,opts=def_opts):
    print("  (using iterative method, may take a bit)")
    t = opts['tolerance']
    imax= opts['max_iter']
    i = 1
    th = lat/2
    err = 1
    while err > t:   #no closed-form solution, so iterate until maximum error is < tolerance
        th = th - (th + np.sin(th)*np.cos(th) + 2*np.sin(th) - (2 + ma.pi/2) * np.sin(lat)) / (2 * np.cos(th) * (1 + np.cos(th)))
        err = np.amax(np.abs((th + np.sin(th)*np.cos(th) + 2*np.sin(th)) - ((2 + ma.pi/2) * np.sin(lat))))
        i+=1
        if i>imax:
            print("  Reached maximum of "+str(imax)+" iterations without converging, outputting result")
            print(f"   max remaining error: {err}")
            break
    r = 1
    x = lon * (1 + np.cos(th)) / (2 * ma.pi)#(1 / np.sqrt(4*ma.pi + ma.pi**2)) * lon * (1 + np.cos(th)) *r 
    y = np.sin(th)#2 * np.sqrt(ma.pi / (4 + ma.pi)) * np.sin(th) *r
    return x,y
coordsl['Eckert IV'] = Eckiv_coords
posl['Eckert IV'] = Eckiv_pos
visl['Eckert IV'] = Pill_vis
ratl['Eckert IV'] = 2
wrapl['Eckert IV'] = 'xwrap'
syml['Eckert IV'] = 'sym4lat'
iter_pos.append('Eckert IV')


proj_list.append('Equal Earth')
info_list['Equal Earth'] = '1.845:1 ovalish (pseudocylindrical); equal-area'

def Eqear_coords(x,y,opts=def_opts):
    return Inv_coords(x,y,Eqear_pos,Wag_coords,Def_vis,opts)
def Eqear_pos(lon,lat,opts=None):
    th = np.arcsin(np.sin(lat) * ma.sqrt(3)/2)
    x = lon * np.cos(th)
    x /= (9*Eqear_c[3]*th**8 + 7*Eqear_c[2]*th**6 + 3*Eqear_c[1]*th**2 + Eqear_c[0])
    y = Eqear_c[3]*th**9 + Eqear_c[2]*th**7 + Eqear_c[1]*th**3 + Eqear_c[0]*th
    x /= (ma.pi / Eqear_c[0])
    y /= Eqear_rat_denom
    return x,y
def Eqear_vis(x,y,lon=None,lat=None,opts=None):
    if lat is None:
        return Edge_vis(x,y,f_coords=Eqear_coords,opts=opts)
    return Coords_vis(x,y,lon,lat,opts)
Eqear_c = [1.340264, -0.081106, 0.000893, 0.003796]
coordsl['Equal Earth'] = Eqear_coords
posl['Equal Earth'] = Eqear_pos
visl['Equal Earth'] = Eqear_vis
Eqear_rat_num = (2*ma.sqrt(3)*ma.pi / 3*Eqear_c[0])
Eqear_rat_th = ma.asin(ma.sqrt(3)/2)
Eqear_rat_denom = Eqear_c[3]*Eqear_rat_th**9 + Eqear_c[2]*Eqear_rat_th**7 + Eqear_c[1]*Eqear_rat_th**3 + Eqear_c[0]*Eqear_rat_th
Eqear_rat = 0.5 * Eqear_rat_num / Eqear_rat_denom
ratl['Equal Earth'] = Eqear_rat
wrapl['Equal Earth'] = 'xwrap'
syml['Equal Earth'] = 'sym4lat'
iter_coords.append('Equal Earth')


sections['Winkel Tripel'] = 'Pseudocylindrical/azimuthal Compromise'


proj_list.append('Winkel Tripel')
info_list['Winkel Tripel'] = '1.637:1 ovalish (pseudoazimuthal)'

def Wink_coords(x,y,opts=def_opts):
    return Inv_coords(x,y,Wink_pos,Ait_coords,Wink_vis,opts)
def Wink_pos(lon,lat,opts=None):
    al = np.arccos(np.cos(lat)*np.cos(lon/2))/ma.pi
    x = (lon/ma.pi + np.cos(lat)*np.sin(lon/2) / np.sinc(al))/(1+ma.pi/2)
    y = (lat + np.sin(lat) / np.sinc(al))/ma.pi
    return x,y
def Wink_vis(x,y,lon=None,lat=None,opts=None):
    if lat is None:
        x1 = np.where(abs(x) > 1, 1, x)
        lat1 = np.arccos(2*(np.abs(x1)*(1+ma.pi/2)-1)/ma.pi)
        ymax = (lat1 + ma.pi/2 * np.sin(lat1))/ma.pi
        return np.where (np.abs(y) > ymax,False,True)
    return Coords_vis(x,y,lon,lat,opts)
coordsl['Winkel Tripel'] = Wink_coords
posl['Winkel Tripel'] = Wink_pos
visl['Winkel Tripel'] = Wink_vis
ratl['Winkel Tripel'] = (1+ma.pi/2)/(ma.pi/2)
wrapl['Winkel Tripel'] = 'xwrap'
syml['Winkel Tripel'] = 'sym4'
iter_coords.append('Winkel Tripel')


proj_list.append('Robinson')
info_list['Robinson'] = '1.972:1 ovalish (pseudocylindrical)'

def Rob_coords(x,y,opts=None):  #Multiquadric interpolation from Ipbuker 2013 https://doi.org/10.1559/1523040041649425
    A_st = np.zeros_like(y)                     #note: there are a some seeming errors and ambiguities in the paper
    B_st = abs(y) * 1.3523                      # where equations 22 and 23 says one should take the sum from j=1 to j=18 of a coefficient * abs(5*j - lat)
    ca = np.asarray(Rob_A) * 0.8487             # it should be from j=0 to j=18, 5 degrees should be converted to radians and lat given in radians,
    cb = np.asarray(Rob_B) * 1.3523             # and lat also given as abs value (nested within the abs) and then resulting y converted back to negative as appropriate
    for j in range(19):                         # and similarly for inverse projection, abs(y) should be given and resulting lat converted back
        A_st += Rob_cm[j] * abs(cb[j] - B_st)
    lon = (x * (0.8487 * ma.pi))/A_st
    lat = np.zeros_like(y)
    for j in range(19):
        lat += Rob_cn[j] * np.sqrt((ca[j] - A_st)**2 + (cb[j] - B_st)**2)
    lat *= np.sign(y)
    return lon,lat
def Rob_pos(lon,lat,opts=None):
    x = np.zeros_like(lat)
    y = np.zeros_like(lat)
    k = np.radians(5)
    alat = abs(lat)
    for j in range(19):
        x += Rob_cp[j] * abs(k*j - alat)
        y += Rob_cq[j] * abs(k*j - alat)
    x = lon * x / (0.8487 * ma.pi) #/ (ma.pi)
    y = np.sign(lat) * y / 1.3523 #/ (ma.pi)
    return x,y
def Rob_vis(x,y,lon=None,lat=None,opts=None):
    if lat is None:
        A_st = np.zeros_like(y)
        B_st = abs(y) * 1.3523
        cb = np.asarray(Rob_B) * 1.3523
        for j in range(19):
            A_st += Rob_cm[j] * abs(cb[j] - B_st)
        return np.where(abs(x) > A_st/0.8487, False, True)
    return Coords_vis(x,y,lon,lat,opts)
Rob_A  = [ 1.0000, 0.9986, 0.9954, 0.9900, 0.9822,
           0.9730, 0.9600, 0.9427, 0.9216, 0.8962,
           0.8679, 0.8350, 0.7986, 0.7597, 0.7186,
           0.6732, 0.6213, 0.5722, 0.5322]
Rob_B  = [ 0.0000, 0.0620, 0.1240, 0.1860, 0.2480,
           0.3100, 0.3720, 0.4340, 0.4958, 0.5571,
           0.6176, 0.6769, 0.7346, 0.7903, 0.8435,
           0.8936, 0.9394, 0.9761, 1.0000]
Rob_cp = [ 0.40711579454, -0.00875326537, -0.01069796348, -0.01167039606, -0.00680782592,
          -0.01847822803, -0.02090931959, -0.01847842619, -0.02090971277, -0.01410147990,
          -0.02236858853, -0.01701955610, -0.01215649454, -0.01069792545, -0.02090967766,
          -0.03160740722,  0.01361549135,  0.04425022432,  0.60843116534]
Rob_cq = [ 0.91083562255, -0.00000589975,  0.00000564852, -0.00000557909,  0.00000555879,
          -0.00000001291, -0.00000546138, -0.00154708482, -0.00387351841, -0.00619324913,
          -0.00930492848, -0.01239340212, -0.01549814705, -0.01937169560, -0.02401844414,
          -0.03331171624, -0.07051393824, -0.09917388904,  0.24527101656]
Rob_cm = [ 0.47371661113, -0.00911028522, -0.01113479305, -0.01214704697, -0.00708577740,
          -0.01923282436, -0.02176345915, -0.01957843209, -0.02288586729, -0.01676092031,
          -0.02731224791, -0.02386224240, -0.02119239013, -0.02327513775, -0.04193330922,
          -0.07123235442, -0.06423048161, -0.10536278437,  1.00598851957]
Rob_cn = [ 1.07729625255, -0.00012324928, -0.00032923415, -0.00056627609, -0.00045168290,
          -0.00141388769, -0.00211521349, -0.00083658786,  0.00073523299,  0.00349045186,
           0.00502041018,  0.00860101415,  0.01281238969,  0.01794606372,  0.02090220870,
           0.02831504310,  0.11177176318,  0.28108668066, -0.45126573496]
coordsl['Robinson'] = Rob_coords
posl['Robinson'] = Rob_pos
visl['Robinson'] = Rob_vis
ratl['Robinson'] = ma.pi * 0.8487/1.3523
wrapl['Robinson'] = 'xwrap'
syml['Robinson'] = 'sym4lat'


proj_list.append('Wagner VI')
info_list['Wagner VI'] = '2:1 ovalish (pseudocylindrical)'

def Wag_coords(x,y,opts=None):
    lat = y*ma.pi/2
    ph = np.arcsin(lat*ma.sqrt(3)/ma.pi)
    lon = x*ma.pi/np.cos(ph)
    return lon,lat
def Wag_pos(lon,lat,opts=None):
    y = lat*2/ma.pi
    x = lon/ma.pi*np.sqrt(1-3*(lat/ma.pi)**2)
    return x,y
def Wag_vis(x,y,lon=None,lat=None,opts=None):
    return np.where(abs(x) > np.sqrt(1-3*(y/2)**2),False,True)
coordsl['Wagner VI'] = Wag_coords
posl['Wagner VI'] = Wag_pos
visl['Wagner VI'] = Wag_vis
ratl['Wagner VI'] = 2
wrapl['Wagner VI'] = 'xwrap'
syml['Wagner VI'] = 'sym4lat'


proj_list.append('Kavrayskiy VII')
info_list['Kavrayskiy VII'] = '1.732:1 ovalish (pseudocylindrical); stretched version of Wagner VI'

coordsl['Kavrayskiy VII'] = Wag_coords
posl['Kavrayskiy VII'] = Wag_pos
visl['Kavrayskiy VII'] = Wag_vis
ratl['Kavrayskiy VII'] = ma.sqrt(3)
wrapl['Kavrayskiy VII'] = 'xwrap'
syml['Kavrayskiy VII'] = 'sym4lat'


proj_list.append('Natural Earth')
info_list['Natural Earth'] = '1.923:1 ovalish (pseudocylindrical)'

def Natear_coords(x,y,opts=def_opts):
    return Inv_coords(x,y,Natear_pos,Wag_coords,Def_vis,opts)
def Natear_pos(lon,lat,opts=None):
    lat2 = lat**2
    lat10 = lat**10
    x = lon * (Natear_ca[0] + Natear_ca[1]*lat2 + Natear_ca[2]*lat**4 + Natear_ca[3]*lat10 + Natear_ca[4]*lat**12)
    y = lat * (Natear_cb[0] + Natear_cb[1]*lat2 + Natear_cb[2]*lat**6 + Natear_cb[3]*lat**8 + Natear_cb[4]*lat10)
    x = x / (Natear_ca[0] * ma.pi)
    y = y / (Natear_rat_denom * ma.pi/2)
    return x,y
def Natear_vis(x,y,lon=None,lat=None,opts=None):
    if lat is None:
        return Edge_vis(x,y,f_coords=Natear_coords,opts=opts)
    return Coords_vis(x,y,lon,lat,opts)
Natear_ca = [0.870700, -0.131979, -0.013791, 0.003971, -0.001529]
Natear_cb = [1.007226,  0.015085, -0.044475, 0.028874, -0.005916]
coordsl['Natural Earth'] = Natear_coords
posl['Natural Earth'] = Natear_pos
visl['Natural Earth'] = Natear_vis
Natear_rat_denom = (Natear_cb[0] + Natear_cb[1]*(ma.pi/2)**2 + Natear_cb[2]*(ma.pi/2)**6 + Natear_cb[3]*(ma.pi/2)**8 + Natear_cb[4]*(ma.pi/2)**10)
Natear_rat = 2 * Natear_ca[0] / Natear_rat_denom
ratl['Natural Earth'] = Natear_rat
wrapl['Natural Earth'] = 'xwrap'
syml['Natural Earth'] = 'sym4lat'
iter_coords.append('Natural Earth')



proj_list.append('Aitoff')
info_list['Aitoff'] = '2:1 ellipse (pseudoazimuthal)'

#def Ait_guess(x,y,opts=None):  #Special routine for Ait_coords initial guess
#    return np.where(Circ_vis(x,y,opts),Hammer_coords(x,y,opts),Wag_coords(x,y,opts))
#def Ait_coords(x,y,opts=None):
#    return Inv_coords(x,y,Ait_pos,Ait_guess,Circ_vis,opts)
def Ait_coords(x,y,opts=None):
    rh = np.sqrt((x/2)**2+(y/2)**2)
    th = np.arctan2(x/2,-y/2)
    lon, lat = From_pol(rh,th)
    return lon*2, lat
def Ait_pos(lon,lat,opts=None):
    al = ma.pi/2*np.sinc(np.arccos(np.cos(lat)*np.cos(lon/2))/ma.pi)
    x = np.cos(lat)*np.sin(lon/2)/al
    y = np.sin(lat)/al
    return x,y
coordsl['Aitoff'] = Ait_coords
posl['Aitoff'] = Ait_pos
visl['Aitoff'] = Circ_vis
ratl['Aitoff'] = 2
wrapl['Aitoff'] = 'xwrap'
syml['Aitoff'] = 'sym4'


sections['Azimuthal Equidistant'] = 'Azimuthal'


def Get_azim_type(opts=def_opts, no_glob=False):     #determine appropriate azim type from opts
    if opts['azim_type'] is not None:
        hem = opts['azim_type']
    elif opts['in']:
        hem = opts['azim_type_in']
    else:
        hem = opts['azim_type_out']
    if hem == 'global_if' and no_glob:
        hem = 'bihem'
    elif hem == 'global' and no_glob:
        print(f'    Warning: global type selected for {"input" if opts["in"] else "output"} projection not possible; using bihemisphere instead')
        hem = 'bihem'
    return hem

def Onehem_coords(x,y,f,s,opts=def_opts):
    return f(x*s,y*s,opts)
def Bihem_coords(x,y,f,s,opts=def_opts):
    lon,lat=Onehem_coords(np.where(x>0,2*x-1,2*x+1),y,f,s,opts)
    lon = np.where(x>0,lon+ma.pi/2,lon-ma.pi/2)
    return lon,lat
def Azim_coords(x,y,f,s,opts=def_opts,no_glob=False):
    hem = Get_azim_type(opts,no_glob)
    if hem in ('global', 'global_if'):
        return f(x,y,opts)
    elif hem == 'hem':
        return Onehem_coords(x,y,f,s,opts)
    elif hem == 'bihem':
        return Bihem_coords(x,y,f,s,opts)
    else:
        raise Exception('Error: invalid azim type')

def Onehem_pos(lon,lat,f,s,opts=def_opts):
    x,y = f(lon,lat,opts)
    return x/s, y/s
def Bihem_pos(lon,lat,f,s,opts=def_opts):
    lon1 = np.where(lon>0,lon,lon+ma.pi)-ma.pi/2
    x,y = Onehem_pos(lon1,lat,f,s,opts)
    x /= 2
    return np.where(lon>0,x+0.5,x-0.5),y
def Azim_pos(lon,lat,f,s,opts=def_opts,no_glob=False):
    hem = Get_azim_type(opts,no_glob)
    if hem in ('global', 'global_if'):
        return f(lon,lat,opts)
    elif hem == 'hem':
        return Onehem_pos(lon,lat,f,s,opts)
    elif hem == 'bihem':
        return Bihem_pos(lon,lat,f,s,opts)
    else:
        raise Exception('Error: invalid azim type')
    
def Azim_vis(x,y,lon=None,lat=None,opts=def_opts,no_glob=False):
    hem = Get_azim_type(opts,no_glob)
    if hem in ('global', 'global_if','hem'):
        return Circ_vis(x,y,opts)
    elif hem == 'bihem':
        return Circ_vis(1-abs(x*2),y,opts)
    else:
        raise Exception('Error: invalid azim type')
def Azim_noglob_vis(x,y,lon=None,lat=None,opts=def_opts):
    return Azim_vis(x,y,lon,lat,opts,no_glob=True)

def Azim_rat(opts=def_opts, no_glob=False):
    hem = Get_azim_type(opts, no_glob)
    if hem in ('global', 'global_if', 'hem'):
        return 1
    elif hem == 'bihem':
        return 2
    else:
        raise Exception('Error: invalid azim type')
def Azim_noglob_rat(opts=def_opts):
    return Azim_rat(opts, no_glob=True)


proj_list.append('Azimuthal Equidistant')
info_list['Azimuthal Equidistant'] = '1:1 circle (azimuthal); equidistant (from center)'

def Azimeq_coords1(x,y,opts=None):
    rh = np.sqrt(x**2+y**2)
    th = np.arctan2(x,-y)
    return From_pol(rh,th)
def Azimeq_pos1(lon,lat,opts=None):
    rh,th = To_pol(lon,lat)
    x = rh*np.sin(th)
    y = rh*np.cos(th)
    return x,y
def Azimeq_coords(x,y,opts=def_opts):
    return Azim_coords(x,y,Azimeq_coords1,1/2,opts)
def Azimeq_pos(lon,lat,opts=def_opts):
    return Azim_pos(lon,lat,Azimeq_pos1,1/2,opts)
coordsl['Azimuthal Equidistant'] = Azimeq_coords
posl['Azimuthal Equidistant'] = Azimeq_pos
visl['Azimuthal Equidistant'] = Azim_vis
ratl['Azimuthal Equidistant'] = Azim_rat
wrapl['Azimuthal Equidistant'] = 'xwrap'
syml['Azimuthal Equidistant'] = 'sym4'
aziml.append('Azimuthal Equidistant')


proj_list.append('Lambert Azimuthal Equal-Area')
info_list['Lambert Azimuthal Equal-Area'] = '1:1 circle (azimuthal); equal-area, preserves 3d chord distances'

def Lamb_coords1(x,y,opts=None):
    rh = np.sqrt(x**2+y**2)
    rh = np.fmod(rh,1)
    th = np.arctan2(x,-y)
    rh1 = 2*np.arcsin(rh) / ma.pi
    return From_pol(rh1,th)
def Lamb_pos1(lon,lat,opts=None):
    rh,th = To_pol(lon,lat)
    rh1 = np.sin(rh*ma.pi/2)
    x = rh1*np.sin(th)
    y = rh1*np.cos(th)
    return x,y
def Lamb_coords(x,y,opts=def_opts):
    return Azim_coords(x,y,Lamb_coords1,ma.sqrt(2)/2,opts)
def Lamb_pos(lon,lat,opts=def_opts):
    return Azim_pos(lon,lat,Lamb_pos1,ma.sqrt(2)/2,opts)
coordsl['Lambert Azimuthal Equal-Area'] = Lamb_coords
posl['Lambert Azimuthal Equal-Area'] = Lamb_pos
visl['Lambert Azimuthal Equal-Area'] = Azim_vis
ratl['Lambert Azimuthal Equal-Area'] = Azim_rat
wrapl['Lambert Azimuthal Equal-Area'] = 'xwrap'
syml['Lambert Azimuthal Equal-Area'] = 'sym4'
aziml.append('Lambert Azimuthal Equal-Area')


proj_list.append('Stereographic')
info_list['Stereographic'] = '1:1 circle (azimuthal); conformal, preserves circles; non-global only'

def Stereo_coords1(x,y,opts=None):
    r = np.sqrt(x**2+y**2)
    th = np.arctan2(x,-y)
    rh = (np.arctan(2*r) - ma.pi/4)*2/ma.pi + 1/2
    return From_pol(rh,th)
def Stereo_pos1(lon,lat,opts=None):
    rh,th = To_pol(lon,lat)
    r = np.tan(ma.pi/4 + ma.pi*(rh-1/2)/2)/2
    x = r*np.sin(th)
    y = r*np.cos(th)
    return x,y
def Stereo_coords(x,y,opts=def_opts):
    s = 1/2
    hem = Get_azim_type(opts)
    if hem == 'global':
        ref = Get_ref(opts, Stereo_def_ref)
        s = ma.tan(ref/2)/2
        return Onehem_coords(x,y,Stereo_coords1,s,opts|{'azim_type':hem})
    return Azim_coords(x,y,Stereo_coords1,s,opts,no_glob=True)
def Stereo_pos(lon,lat,opts=def_opts):
    s = 1/2
    hem = Get_azim_type(opts)
    if hem == 'global':
        hem = 'hem'
        ref = Get_ref(opts, Stereo_def_ref)
        s = ma.tan(ref/2)/2
        return Onehem_pos(lon,lat,Stereo_pos1,s,opts|{'azim_type':hem})
    return Azim_pos(lon,lat,Stereo_pos1,s,opts,no_glob=True)
Stereo_def_ref = 90
coordsl['Stereographic'] = Stereo_coords
posl['Stereographic'] = Stereo_pos
visl['Stereographic'] = Azim_vis
ratl['Stereographic'] = Azim_rat
wrapl['Stereographic'] = 'xwrap'
syml['Stereographic'] = 'sym4'
aziml.append('Stereographic')
nonglobal.append('Stereographic')
optlists_global['Stereographic'] = {
    'Select maximum latitude range to truncate map (0-180 degrees from center)': 1
}


proj_list.append('Orthographic')
info_list['Orthographic'] = '1:1 circle (azimuthal); appearance of 3d globe; hemispheres only'

def Ortho_coords1(x,y,opts=None):
    rh = np.sqrt((2*x)**2+(2*y)**2)
    c = np.arcsin(np.fmod(rh,1))
    lat = np.arcsin(2*y)#*np.sin(c)/rh)
    #lon = np.arcsin(2*x / np.sqrt(1-(2*y)**2))
    lon = np.arctan(2*x*np.sin(c)/(rh*np.cos(c)))
    return lon,lat
def Ortho_pos1(lon,lat,opts=None):
    x = np.cos(lat)*np.sin(lon)/2
    x = np.where(np.abs(lon) > ma.pi/2, 1, x)  #ensures outside hemisphere doesn't index into map area
    y = np.sin(lat)/2
    return x,y
def Ortho_coords(x,y,opts=def_opts):
    return Azim_coords(x,y,Ortho_coords1,1/2,opts,no_glob=True)
def Ortho_pos(lon,lat,opts=def_opts):
    return Azim_pos(lon,lat,Ortho_pos1,1/2,opts,no_glob=True)
coordsl['Orthographic'] = Ortho_coords
posl['Orthographic'] = Ortho_pos
visl['Orthographic'] = Azim_noglob_vis
ratl['Orthographic'] = Azim_noglob_rat
wrapl['Orthographic'] = 'xwrap'
syml['Orthographic'] = 'sym4lat'
aziml.append('Orthographic')
nonglobal.append('Orthographic')


sections['Equidistant Conic'] = 'Conic'


proj_list.append('Equidistant Conic')
info_list['Equidistant Conic'] = 'partial circular arc (conic); equidistant (from top edge)'

def Coniceq_coords(x,y,opts=def_opts):
    n, G, max_x, max_y, hem = Coniceq_consts(opts)
    x1 = x * max_x
    if hem == 'hem':
        y1 = (y + 1) * max_y/2
    elif hem == 'bihem':
        y1 = np.abs(y * max_y)
    else:
        y1 = y * (max_y/2 + ma.pi/4) + (max_y/2 - ma.pi/4)
    rho = np.sign(n) * np.sqrt(x1**2 + (G-y1)**2)
    lat = G - rho
    th = np.arctan2(x1, (G - y1),)
    lon = th / n
    if hem == 'bihem':
        lat *= np.sign(y)
    return lon, lat
def Coniceq_pos(lon,lat,opts=def_opts):
    n, G, max_x, max_y, hem = Coniceq_consts(opts)
    rho = G - (np.abs(lat) if hem == 'bihem' else lat)
    th = n*lon
    x = rho * np.sin(th)
    y = G - rho * np.cos(th)
    x /= max_x
    if hem == 'hem':
        y = (y / max_y) * 2 - 1
    elif hem == 'bihem':
        y = y / max_y
        y *= np.sign(lat)
    else:
        y = (y - (max_y/2 - ma.pi/4)) / (max_y/2 + ma.pi/4)
    return x, y
def Coniceq_vis(x,y,lon=None,lat=None,opts=def_opts,get_far=False):
    if lat is None:
        n, G, max_x, max_y, hem = Coniceq_consts(opts)
        x1 = x * max_x
        if hem == 'hem':
            y1 = (y + 1) * max_y/2
        elif hem == 'bihem':
            y1 =np.abs(y * max_y)
        else:
            y1 = y * (max_y/2 + ma.pi/4) + (max_y/2 - ma.pi/4)
        rho_s = x1**2 + (G-y1)**2
        close = rho_s < (G - ma.pi/2)**2
        if hem in ('hem','bihem'):
            far = rho_s > G**2
        else:
            far = rho_s > (G + ma.pi/2)**2
        sides = np.arctan2(np.abs(x1) , (G - y1)) > ma.pi*n
        vis = np.logical_not(close | far | sides)
        if get_far:
            return vis,far
        return vis
    hem = Get_azim_type(opts)
    lat1 = lat.copy()
    if hem in ('hem','bihem'):
        if hem == 'bihem':
            lat1 *= np.sign(y)  #flip sign in bottom map half
        lat1 = lat1*2 - ma.pi/2   #shift lat so that whole southern hemisphere is treated as off-map
    far = lat1 < -ma.pi/2
    vis = Coords_vis(x,y,lon,lat1,opts)
    if get_far:
        return vis,far
    return vis
def Coniceq_rat(opts=def_opts):
    n, G, max_x, max_y, hem = Coniceq_consts(opts)
    if hem == 'hem':
        return max_x*2 / (max_y)
    if hem == 'bihem':
        return max_x / max_y
    return max_x*2 / (max_y + ma.pi/2)
def Coniceq_consts(opts=def_opts):      #common constants used in all cases
    ref = Get_ref(opts,Coniceq_def_ref)
    if ref[0] == ref[1]:
        n = ma.sin(ref[0])
    else:
        n = (ma.cos(ref[0]) - ma.cos(ref[1])) / (ref[1] - ref[0])
    G = ma.cos(ref[0])/n + ref[0]
    hem = Get_azim_type(opts)
    if n > 1/2:     #determines if map curves more than a half circle
        if hem in ('hem','bihem'):
            max_x = G
            max_y = G - G * ma.cos(n * ma.pi)
        else:
            max_x = G + ma.pi/2     #rightmost edge of outer perimeter
            max_y = G - (G + ma.pi/2) * ma.cos(n * ma.pi)   #right end of outer perimeter
    else:
        if hem in ('hem','bihem'):
            max_x = G * ma.sin(n * ma.pi)
        else:
            max_x = (G + ma.pi/2) * ma.sin(n * ma.pi)       #right end of outer perimeter
        max_y = G - (G - ma.pi/2) * ma.cos(n * ma.pi)   #right end of inner perimeter
    return n, G, max_x, max_y, hem
Coniceq_def_ref = (15,45)
coordsl['Equidistant Conic'] = Coniceq_coords
posl['Equidistant Conic'] = Coniceq_pos
visl['Equidistant Conic'] = Coniceq_vis
ratl['Equidistant Conic'] = Coniceq_rat
wrapl['Equidistant Conic'] = 'conic_xwrap'
syml['Equidistant Conic'] = 'symx'
conicl.append('Equidistant Conic')
optlists['Equidistant Conic'] = {
    'Select reference latitudes of best shape (-90-90 degrees from equator, mean must be positive)': 2
}


proj_list.append('Albers Equal-Area Conic')
info_list['Albers Equal-Area Conic'] = 'partial circular arc (conic); equal-area'

def Conicalb_coords(x,y,opts=def_opts):
    n, C, rho0, max_x, max_y, min_y, hem = Conicalb_consts(opts)
    x1 = x * max_x
    if hem == 'bihem':
        y1 = np.abs(y) * max_y
    else:
        y1 = y * (max_y - min_y)/2 + (max_y + min_y)/2
    rho = np.sqrt(x1**2 + (rho0 - y1)**2)
    th = np.arctan2(x1 , (rho0 - y1))
    lon = th / n
    sinla = (C - (rho**2 * n**2)) / (2*n)
    sinla1 = np.maximum(sinla,-1)
    sinla1 = np.minimum(sinla1,1)       #avoid arcsin errors
    lat = np.arcsin(sinla1)
    lat = np.where(sinla > 1, 0.01+ma.pi/2, lat)
    lat = np.where(sinla < -1, -0.01-ma.pi/2, lat)   #make sure visible area clips properly
    if hem == 'bihem':
        lat *= np.sign(y)
    return lon, lat
def Conicalb_pos(lon,lat,opts=def_opts):
    n, C, rho0, max_x, max_y, min_y, hem = Conicalb_consts(opts)
    rho = np.sqrt(C - 2 * n * np.sin(np.abs(lat) if hem == 'bihem' else lat)) / n
    th = n*lon
    x = rho * np.sin(th)
    y = rho0 - rho * np.cos(th)
    x /= max_x
    if hem == 'bihem':
        y = np.sign(lat) * y / max_y
    else:
        y = (y - (max_y + min_y)/2) / ((max_y - min_y)/2)
    return x, y
def Conicalb_vis(x,y,lon=None,lat=None,opts=def_opts,get_far=False):
    if lat is None:
        n, C, rho0, max_x, max_y, min_y, hem = Conicalb_consts(opts)
        x1 = x * max_x
        if hem == 'bihem':
            y1 = np.abs(y) * max_y
        else:
            y1 = y * (max_y - min_y)/2 + (max_y + min_y)/2
        rho_s = x1**2 + (rho0 - y1)**2
        close = rho_s < (C - 2*n) / n**2
        if hem in ('hem','bihem'):
            far = rho_s > C / n**2
        else:
            far = rho_s > (C + 2*n) / n**2
        sides = np.arctan2(np.abs(x1) , (rho0 - y1)) > n*ma.pi
        vis = np.logical_not(close | far | sides)
        if get_far:
            return vis,far
        return vis
    hem = Get_azim_type(opts)
    lat1 = lat.copy()
    if hem in ('hem','bihem'):
        if hem == 'bihem':
            lat1 *= np.sign(y)  #flip sign in bottom map half
        lat1 = lat1*2 - ma.pi/2   #shift lat so that whole southern hemisphere is treated as off-map
    far = lat1 < -ma.pi/2
    vis = Coords_vis(x,y,lon,lat1,opts)
    if get_far:
        return vis,far
    return vis
def Conicalb_rat(opts=def_opts):
    n, C, rho0, max_x, max_y, min_y, hem = Conicalb_consts(opts)
    if hem == 'bihem':
        return max_x / max_y
    return max_x*2 / (max_y - min_y)
def Conicalb_consts(opts=def_opts):      #common constants used in all cases
    ref = Get_ref(opts,Coniceq_def_ref)
    if ref[0] == ref[1]:
        n = ma.sin(ref[0])
    else:
        n = (ma.sin(ref[0]) + ma.sin(ref[1])) / 2
    C = ma.cos(ref[0])**2 + 2 * n * ma.sin(ref[0])
    rho0 = np.sqrt(C)/n
    hem = Get_azim_type(opts)
    if hem in ('hem','bihem'):
        rhomax = rho0
        min_y = 0
    else:
        rhomax = ma.sqrt(C + 2*n) / n
        min_y = rho0 - rhomax
    if n > 1/2:     #determines if map curves more than a half circle
        max_x = rhomax     #rightmost edge of outer perimeter
        max_y = rho0 - rhomax * ma.cos(n*ma.pi)   #right end of outer perimeter
    else:
        max_x = rhomax * ma.sin(n*ma.pi)       #right end of outer perimeter
        max_y = rho0 - (ma.sqrt(C - 2*n) / n) * ma.cos(n*ma.pi)   #right end of inner perimeter
    
    return n, C, rho0, max_x, max_y, min_y, hem
Conicalb_def_ref = (15,45)
coordsl['Albers Equal-Area Conic'] = Conicalb_coords
posl['Albers Equal-Area Conic'] = Conicalb_pos
visl['Albers Equal-Area Conic'] = Conicalb_vis
ratl['Albers Equal-Area Conic'] = Conicalb_rat
wrapl['Albers Equal-Area Conic'] = 'conic_xwrap'
syml['Albers Equal-Area Conic'] = 'symx'
conicl.append('Albers Equal-Area Conic')
optlists['Albers Equal-Area Conic'] = {
    'Select reference latitudes of best shape (-90-90 degrees from equator, mean must be positive)': 2
}


proj_list.append('Lambert Conformal Conic')
info_list['Lambert Conformal Conic'] = 'partial circular arc (conic); conformal'


def Coniclamb_coords(x,y,opts=def_opts):
    n, F, rhomax, max_x, max_y, min_y, hem = Coniclamb_consts(opts)
    x1 = x * max_x
    if hem == 'bihem':
        y1 = np.abs(y) * max_y
    else:
        y1 = y * (max_y - min_y)/2 + (max_y + min_y)/2
    rho = np.sign(n) * np.sqrt(x1**2 + (F - y1)**2)
    th = np.arctan2(x1 , (F - y1))
    lon = th / n
    lat = 2 * np.arctan((F/rho)**(1/n)) - ma.pi / 2
    if hem == 'bihem':
        lat *= np.sign(y)
    return lon, lat
def Coniclamb_pos(lon,lat,opts=def_opts):
    n, F, rhomax, max_x, max_y, min_y, hem = Coniclamb_consts(opts)
    rho = F / np.tan(ma.pi/4 + (np.abs(lat) if hem == 'bihem' else lat)/2)**n
    th = n * lon
    x = rho * np.sin(th)
    y = F - rho * np.cos(th)
    x /= max_x
    if hem == 'bihem':
        y = np.sign(lat) * y / max_y
    else:
        y = (y - (max_y + min_y)/2) / ((max_y - min_y)/2)
    return x, y
def Coniclamb_vis(x,y,lon=None,lat=None,opts=def_opts,get_far=False):
    if lat is None:
        n, F, rhomax, max_x, max_y, min_y, hem = Coniclamb_consts(opts)
        x1 = x * max_x
        if hem == 'bihem':
            y1 = np.abs(y) * max_y
        else:
            y1 = y * (max_y - min_y)/2 + (max_y + min_y)/2
        rho_s = x1**2 + (F - y1)**2
        far = rho_s > rhomax**2
        sides = np.arctan2(np.abs(x1) , (F - y1)) > n*ma.pi
        vis = np.logical_not(far | sides)
        if get_far:
            return vis,far
        return vis
    hem = Get_azim_type(opts)
    lat1 = lat.copy()
    if hem in ('hem','bihem'):
        if hem == 'bihem':
            lat1 *= np.sign(y)  #flip sign in bottom map half
        lat1 = lat1*2 - ma.pi/2   #shift lat so that whole southern hemisphere is treated as off-map
    else:
        ref = Get_ref(opts,Coniclamb_def_ref)
        lat1 = np.where(lat > ref[0], lat, -2)
    far = lat1 < -ma.pi/2
    vis = Coords_vis(x,y,lon,lat1,opts)
    if get_far:
        return vis,far
    return vis
def Coniclamb_rat(opts=def_opts):
    n, F, rhomax, max_x, max_y, min_y, hem = Coniclamb_consts(opts)
    if hem == 'bihem':
        return max_x / max_y
    return max_x*2 / (max_y - min_y)
def Coniclamb_consts(opts=def_opts):
    ref = Get_ref(opts,Coniclamb_def_ref)
    if len(ref) < 3:
        ref = (0,ref[0],ref[1])
    n = (ma.log(ma.cos(ref[1])/ma.cos(ref[2])) /
         ma.log(ma.tan(ma.pi/4 + ref[2]/2) / ma.tan(ma.pi/4 + ref[1]/2)))
    F = ma.cos(ref[1]) * (ma.tan(ma.pi/4 + ref[1]/2) ** n) / n
    hem = Get_azim_type(opts)
    if hem in ('hem','bihem'):
        rhomax = F
        min_y = 0
    else:
        rhomax = F / (ma.tan(ma.pi/4 + (ref[0])/2) ** n)
        min_y = F - rhomax
    if n > 1/2:     #determines if map curves more than a half circle
        max_x = rhomax     #rightmost edge of outer perimeter
        max_y = F - rhomax * ma.cos(n*ma.pi)   #right end of outer perimeter
    else:
        max_x = rhomax * ma.sin(n*ma.pi)       #right end of outer perimeter
        max_y = F    #top pole
    
    return n, F, rhomax, max_x, max_y, min_y, hem
    
Coniclamb_def_ref = (0,15,45)
coordsl['Lambert Conformal Conic'] = Coniclamb_coords
posl['Lambert Conformal Conic'] = Coniclamb_pos
visl['Lambert Conformal Conic'] = Coniclamb_vis
ratl['Lambert Conformal Conic'] = Coniclamb_rat
wrapl['Lambert Conformal Conic'] = 'conic_xwrap'
syml['Lambert Conformal Conic'] = 'symx'
conicl.append('Lambert Conformal Conic')
nonglobal.append('Lambert Conformal Conic')
optlists['Lambert Conformal Conic'] = {
    'Select reference latitudes of best shape (-90-90 degrees from equator, mean must be positive)': 2
}
optlists_global['Lambert Conformal Conic'] = {
    'Select minimum latitude to truncate map (-90-90 degrees from equator)': 1
}


sections['Nicolosi Globular'] = 'Miscellaneous / Special-Use'


proj_list.append('Nicolosi Globular')
info_list['Nicolosi Globular'] = '1:1 circle (polyconic); hemispheres only'

#Bit messy right now, seems to be working fine for hem and bihem but having trouble getting global working


# '1.354:1 clamshell for when I get global working

def Nic_coords1(x,y,opts=def_opts):     #wip better Nicolosi inverse that should for global maps, but seems to have some sign issues
    print("  (using iterative method, may take a bit)")
    t = opts['tolerance']
    imax = opts['max_iter']#*100
    hem = Get_azim_type(opts)
    if hem == 'global':
        x1 = x * (1/2 + ma.sqrt(2))
        y1 = y * ma.sqrt(2)
    x2 = x1**2
    x2y2 = x2 + y1**2
    lon = ((x2y2 - 1 + np.sqrt((1-x2y2)**2 + 4 * x2)) / (2 * x1)) * ma.pi/2
    err = 1
    i = 1
    lat = y1 * ma.pi/2
    while err > t:
        #print(np.max(lat))
        sinla = np.sin(lat)
        cosla = np.sin(lat)
        lat2 = lat**2
        l1 = x2y2 * (ma.pi * sinla - 2 * lat) * ma.pi + 4 * lat2 * (y1*(2*ma.pi) - sinla) + 2 * ma.pi * lat - (ma.pi**2) * y1
        l2 = x2y2 * (ma.pi * cosla - 2) * ma.pi + 8 * lat * (y1 - sinla) - 4 * lat2 * cosla + 2 * ma.pi
        adjust = l1/l2
        lat = (lat - adjust)#*np.sign(y))
        #lat = lat * np.sign(lat) * np.sign(y)
        #lat *= np.sign(lat) * np.sign(y)
        err = np.max(np.abs(adjust))
        i+=1
        if i>imax:
            print("  Reached maximum of "+str(imax)+" iterations without converging, outputting result")
            print(f"   max remaining error: {err}")
            break
    #lat = y * ma.pi/2
    #lat *= np.sign(y)
    return lon, lat

def Nic_coords2(x,y,opts=def_opts):
    return Inv_coords(x,y,Nic_pos,Azimeq_coords,Nic_vis,opts|{'azim_type':Get_azim_type(opts,no_glob=True)})
def Nic_pos1(lon,lat,opts=def_opts):
    lat0 = lat == 0
    lon0 = lon == 0
    latmax = np.abs(lat) == ma.pi/2
    lonmax = np.abs(lon) == ma.pi/2
    #lat_lim = lat0 & lon0 & latmax & lonmax
    b = ma.pi/(2*lon) - 2*lon/ma.pi
    c = 2 * lat / ma.pi
    sinla = np.sin(lat)
    d = (1-c**2) / (sinla - c)
    b2 = b**2
    d2 = d**2
    b2d2 = 1 + b2/d2
    d2b2 = 1 + d2/b2
    M = (b*sinla/d - b/2) / b2d2
    N = (d2*sinla/b2 + d/2) / d2b2
    x = M + np.sign(lon) * np.sqrt(M**2 + np.cos(lat)**2 / b2d2)
    y = N + np.sign(lon) * np.sign(-lat * b) * np.sqrt(N**2 - ((d2b2-1)*sinla**2 + d*sinla - 1) / d2b2)    #'just do nicolosi' I thought, 'it'll be fun and it seems pretty straightforward, right?'
    x[lon0 & latmax] = 0
    x[lat0] = lon[lat0] / ma.pi
    x[lonmax] = lon[lonmax] * np.cos(lat[lonmax]) / ma.pi
    y[lon0] = lat[lon0] * 2/ma.pi
    y[latmax] = np.sign(lat[latmax])
    y[lat0] = 0
    y[lonmax] = sinla[lonmax]
    if Get_azim_type(opts) == 'global':
        x = x / (1/2 + ma.sqrt(2))
        y = y / (ma.sqrt(2))
    else:
        x = x / 2
        y = y / 2
    return x, y

def Nic_pos2(lon,lat,opts=def_opts):    #old nicolosi method that won't work for global maps but I'm keeping it for reference
    R = 2 / ma.pi
    lat0 = np.where(lat == 0, True, False)
    lon0 = np.where(lon == 0, True, False)
    latmax = np.where(np.abs(lat) == ma.pi/2, True, False)
    lonmax = np.where(np.abs(lon) == ma.pi/2, True, False)
    lat1 = np.where(lat0 | latmax, lat+1e-8, lat)
    lon1 = np.where(lon0 | lonmax, lon+1e-8, lon)
    sinla = np.sin(lat1)
    b = ma.pi/(2*lon1) - 2*lon1/ma.pi
    c = 2*lat1/ma.pi
    d = (1-c**2) / (sinla - c)
    b2 = b**2
    d2 = d**2
    b2d2 = 1 + b2/d2
    d2b2 = 1 + d2/b2
    M = (b*sinla/d - b/2) / b2d2
    N = (d2*sinla/b2 + d/2) / d2b2
    x1 = np.sqrt(M**2 + np.cos(lat1)**2 / b2d2)
    x = np.where(lon>0, M+x1, M-x1) * R * ma.pi / 2
    y1 = np.sqrt(N**2 - (d2/b2 * sinla**2 + d*sinla - 1) / d2b2)
    x = np.where(lon0 | latmax,0,x)
    x = np.where(lat0 | lonmax,np.cos(lat)*lon*R,x)
    y = np.where(lon0 | lat0,R*lat,y)
    y = np.where(lonmax,np.sin(lat)*R*ma.pi/2,y)
    y = np.where(latmax,R*lat,y)
    x = np.where(lon > ma.pi/2, ma.pi - x, np.where(lon < -ma.pi/2, -ma.pi - x, x))
    return x,y
def Nic_coords(x,y,opts=def_opts):
    hem = Get_azim_type(opts,no_glob=True)
    if hem in ('hem'):
        return Azim_coords(x,y,Nic_coords2,1/2,opts,no_glob=True)
    else:
        x1 = np.abs(x) * 2 - 1
        lon,lat = Onehem_coords(x1,y,Nic_coords2,1/2,opts|{'azim_type':'hem'})
        lon *= np.sign(x)
        return lon,lat
        x1,x2 = np.array_split(x,2,1)
        y1,y2 = np.array_split(y,2,1)
        lon1,lat1 = Onehem_coords(2*x1+1,y1,Nic_coords1,1,opts|{'azim_type':'hem'})
        lon1 = lon1 - ma.pi/2
        lon2,lat2 = Onehem_coords(2*x2-1,y2,Nic_coords1,1,opts|{'azim_type':'hem'})
        lon2 = lon2 + ma.pi/2
        lon = np.concatenate((lon1,lon2),1)
        lat = np.concatenate((lat1,lat2),1)
        #lon = np.where(x>0,lon+ma.pi/2,lon-ma.pi/2)
        return lon,lat
def Nic_pos(lon,lat,opts=def_opts):
    return Azim_pos(lon,lat,Nic_pos1,1/2,opts,no_glob=True)
def Nic_vis(x,y,lon=None,lat=None,opts=def_opts):
    hem = Get_azim_type(opts,no_glob=True)
    x1 = np.copy(x)
    y1 = np.copy(y)
    if hem == 'global':
        if True:#lat is None:
            x1 = ((1/2 + ma.sqrt(2)) * np.abs(x1) - 1/2) / ma.sqrt(2)
        else:
            return Coords_vis(x,y,lon,lat,opts)
    return Azim_vis(x1,y1,lon,lat,opts,no_glob=True)
def Nic_rat(opts=def_opts):
    hem = Get_azim_type(opts,no_glob=True)
    if hem == 'global':
        return (1/2 + ma.sqrt(2))/ma.sqrt(2)
    return Azim_rat(opts,no_glob=True)

coordsl['Nicolosi Globular'] = Nic_coords2
posl['Nicolosi Globular'] = Nic_pos
visl['Nicolosi Globular'] = Nic_vis#Azim_noglob_vis
ratl['Nicolosi Globular'] = Nic_rat#Azim_noglob_rat
wrapl['Nicolosi Globular'] = 'xwrap'
syml['Nicolosi Globular'] = 'sym4'
aziml.append('Nicolosi Globular')
nonglobal.append('Nicolosi Globular')
iter_coords.append('Nicolosi Globular')


proj_list.append('Ortelius Oval')
info_list['Ortelius Oval'] = '2:1 oval (pseudocylindrical)'

def Ort_coords1(x,y,opts=def_opts):
    lon,lat = Inv_coords(x,y,Ort_pos1,Hammer_coords,Circ_vis,opts)
    lat = y*ma.pi/2
    return lon,lat
def Ort_coords2(x,y,opts=def_opts):
    lon1,lat1 = Inv_coords(x,y,Ort_pos2,Hammer_coords,Circ_vis,opts)
    lat = y*ma.pi/2
    lon2 = np.abs(x)*ma.pi + ma.pi/2 - np.sqrt(ma.pi**2/4 - lat**2)
    lon2 = np.where(x>0, lon2, -lon2)
    lon = np.where(np.sqrt((x*2)**2+y**2)>1,lon2,lon1)
    return lon,lat
def Ort_pos1(lon,lat,opts=None):
    y = lat*2/ma.pi
    abslon = np.abs(lon)
    F = (ma.pi**2/(4*abslon)+abslon)/2
    x = abslon - F + np.sqrt(F**2 - (y*ma.pi/2)**2)
    x = 2 * np.sign(lon) * x/ma.pi 
    return x,y
def Ort_pos2(lon,lat,opts=None):
    y = lat*2/ma.pi
    abslon = np.abs(lon)
    F = (ma.pi**2/(4*abslon)+abslon)/2
    x = np.where(abslon>ma.pi/2,
                 np.sqrt(ma.pi**2/4 - lat**2) + abslon - ma.pi/2,
                 abslon - F + np.sqrt(F**2 - (y*ma.pi/2)**2))
    x = np.where(lon>0, x/ma.pi, -x/ma.pi)
    return x,y
def Ort_coords(x,y,opts=def_opts):
    if Get_azim_type(opts) in ('hem', 'bihem'):
        return Azim_coords(x,y,Ort_coords1,1,opts)
    return Ort_coords2(x,y,opts)
def Ort_pos(lon,lat,opts=def_opts):
    if Get_azim_type(opts) in ('hem', 'bihem'):
        return Azim_pos(lon,lat,Ort_pos1,1,opts)
    return Ort_pos2(lon,lat,opts)
def Ort_vis(x,y,lat=None,lon=None,opts=def_opts):
    if Get_azim_type(opts) in ('bihem','hem'):
        return Azim_vis(x,y,lat,lon,opts)
    return Pill_vis(x,y,opts)
def Ort_rat(opts=def_opts):
    if Get_azim_type(opts) in ('hem'):
        return 1
    return 2
coordsl['Ortelius Oval'] = Ort_coords
posl['Ortelius Oval'] = Ort_pos
visl['Ortelius Oval'] = Ort_vis
ratl['Ortelius Oval'] = Ort_rat
wrapl['Ortelius Oval'] = 'xwrap'
syml['Ortelius Oval'] = 'sym4lat'
aziml.append('Ortelius Oval')
iter_coords.append('Ortelius Oval')


proj_list.append('Pseudostereographic')
info_list['Pseudostereographic'] = '2:1 ellipse (pseudoazimuthal)'

def Psstereo_coords(x,y,opts=def_opts):
    lon, lat = Onehem_coords(x,y,Stereo_coords1,1/2,opts)
    return lon*2, lat
def Psstereo_pos(lon,lat,opts=def_opts):
    return Onehem_pos(lon/2,lat,Stereo_pos1,1/2,opts)
coordsl['Pseudostereographic'] = Psstereo_coords
posl['Pseudostereographic'] = Psstereo_pos
visl['Pseudostereographic'] = Circ_vis
ratl['Pseudostereographic'] = 2
wrapl['Pseudostereographic'] = 'xwrap'
syml['Pseudostereographic'] = 'sym4'


proj_list.append('Pseudoorthographic')
info_list['Pseudoorthographic'] = '2:1 ellipse (pseudoazimuthal)'

def Psortho_coords(x,y,opts=def_opts):
    lon, lat = Onehem_coords(x,y,Ortho_coords1,1/2,opts)
    return lon*2, lat
def Psortho_pos(lon,lat,opts=def_opts):
    return Onehem_pos(lon/2,lat,Ortho_pos1,1/2,opts)
coordsl['Pseudoorthographic'] = Psortho_coords
posl['Pseudoorthographic'] = Psortho_pos
visl['Pseudoorthographic'] = Circ_vis
ratl['Pseudoorthographic'] = 2
wrapl['Pseudoorthographic'] = 'xwrap'
syml['Pseudoorthographic'] = 'sym4'


# Experimental stuff I might add into the official list later:

#proj_list.append('Pseudonicolosi')
info_list['Pseudonicolosi'] = '2:1 ellipse (polyconic?)'

def Psnic_coords(x,y,opts=def_opts):
    lon, lat = Onehem_coords(x,y,Nic_coords2,1,opts|{'azim_type':'hem'})
    return lon*2, lat
def Psnic_pos(lon,lat,opts=def_opts):
    return Onehem_pos(lon/2,lat,Nic_pos1,1/2,opts)
coordsl['Pseudonicolosi'] = Psnic_coords
posl['Pseudonicolosi'] = Psnic_pos
visl['Pseudonicolosi'] = Circ_vis
ratl['Pseudonicolosi'] = 2
wrapl['Pseudonicolosi'] = 'xwrap'
syml['Pseudonicolosi'] = 'sym4'


#proj_list.append('Sheared Azimuthal Bihemisphere')
info_list['Sheared Azimuthal Bihemisphere'] = ''

def Azimeqsh_coords(x,y,opts=def_opts):
    x1 = x - 0.5 - np.sign(x) * np.sqrt(1 - y**2)/2
    lon, lat = Bihem_coords(x1,y,Azimeq_coords1,1/2,opts)
    lat = np.where(np.abs(x1) < 0.5, y * ma.pi/2, lat)
    return lon, lat

def Azimeqsh_pos(lon,lat,opts=def_opts):
    x,y = Bihem_pos(lon,lat,Azimeq_pos1,1/2,opts)
    x1 = x + 0.5 + np.sign(x) * np.sqrt(1 - y**2)/2
    y = np.where(np.abs(lon) < ma.pi/2, lat * 2/ma.pi, y)
    #x2,y2 = Hammer_pos(lon,lat,opts)
    #x = (x1 + x2)/2
    #y = (y + y2)/2
    return x1,y
coordsl['Sheared Azimuthal Bihemisphere'] = Azimeqsh_coords
posl['Sheared Azimuthal Bihemisphere'] = Azimeqsh_pos
visl['Sheared Azimuthal Bihemisphere'] = Circ_vis
ratl['Sheared Azimuthal Bihemisphere'] = 2
wrapl['Sheared Azimuthal Bihemisphere'] = 'xwrap'
syml['Sheared Azimuthal Bihemisphere'] = 'sym4'

#proj_list.append('Sheared Equal-Area Bihemisphere')
info_list['Sheared Equal-Area Bihemisphere'] = ''

def Lambsh_coords(x,y,opts=def_opts):
    x1 = x - 0.5 - np.sign(x) * np.sqrt(1 - y**2)/2
    lon, lat = Bihem_coords(x1,y,Lamb_coords1,ma.sqrt(2)/2,opts)
    lon1, lat1 = Onehem_coords(0,y,Lamb_coords1,ma.sqrt(2)/2,opts)
    lat = np.where(np.abs(x1) < 0.5, lat1, lat)
    return lon, lat

def Lambsh_pos(lon,lat,opts=def_opts):
    x,y = Bihem_pos(lon,lat,Lamb_pos1,ma.sqrt(2)/2,opts)
    x1 = x + 0.5 + np.sign(x) * np.sqrt(1 - y**2)/2
    x2, y1 = Onehem_pos(0,lat,Lamb_pos1,ma.sqrt(2)/2,opts)
    y = np.where(np.abs(lon) < ma.pi/2, y1, y)
    #x2,y2 = Hammer_pos(lon,lat,opts)
    #x = (x1 + x2)/2
    #y = (y + y2)/2
    return x1,y
coordsl['Sheared Equal-Area Bihemisphere'] = Lambsh_coords
posl['Sheared Equal-Area Bihemisphere'] = Lambsh_pos
visl['Sheared Equal-Area Bihemisphere'] = Circ_vis
ratl['Sheared Equal-Area Bihemisphere'] = 2
wrapl['Sheared Equal-Area Bihemisphere'] = 'xwrap'
syml['Sheared Equal-Area Bihemisphere'] = 'sym4'

#proj_list.append('Sheared Stereographic Bihemisphere')
info_list['Sheared Stereographic Bihemisphere'] = ''

def Stereosh_coords(x,y,opts=def_opts):
    x1 = x - 0.5 - np.sign(x) * np.sqrt(1 - y**2)/2
    lon, lat = Bihem_coords(x1,y,Stereo_coords1,1/2,opts)
    lon1, lat1 = Onehem_coords(0,y,Stereo_coords1,1/2,opts)
    lat = np.where(np.abs(x1) < 0.5, lat1, lat)
    return lon, lat

def Stereosh_pos(lon,lat,opts=def_opts):
    x,y = Bihem_pos(lon,lat,Stereo_pos1,1/2,opts)
    x1 = x + 0.5 + np.sign(x) * np.sqrt(1 - y**2)/2
    x2, y1 = Onehem_pos(0,lat,Stereo_pos1,1/2,opts)
    y = np.where(np.abs(lon) < ma.pi/2, y1, y)
    #x2,y2 = Hammer_pos(lon,lat,opts)
    #x = (x1 + x2)/2
    #y = (y + y2)/2
    return x1,y
coordsl['Sheared Stereographic Bihemisphere'] = Stereosh_coords
posl['Sheared Stereographic Bihemisphere'] = Stereosh_pos
visl['Sheared Stereographic Bihemisphere'] = Circ_vis
ratl['Sheared Stereographic Bihemisphere'] = 2
wrapl['Sheared Stereographic Bihemisphere'] = 'xwrap'
syml['Sheared Stereographic Bihemisphere'] = 'sym4'

#proj_list.append('Sheared Orthographic Bihemisphere')
info_list['Sheared Orthographic Bihemisphere'] = ''

def Orthosh_coords(x,y,opts=def_opts):
    x1 = x - 0.5 - np.sign(x) * np.sqrt(1 - y**2)/2
    lon, lat = Bihem_coords(x1,y,Ortho_coords1,1/2,opts)
    lon1, lat1 = Onehem_coords(0,y,Ortho_coords1,1/2,opts)
    #lat = np.where(np.abs(x1) < 0.5, lat1, lat)
    return lon, lat

def Orthosh_pos(lon,lat,opts=def_opts):
    x,y = Bihem_pos(lon,lat,Ortho_pos1,1/2,opts)
    x1 = x + 0.5 + np.sign(x) * np.sqrt(1 - y**2)/2
    x2, y1 = Onehem_pos(0,lat,Ortho_pos1,1/2,opts)
    #y = np.where(np.abs(lon) < ma.pi/2, y1, y)
    #x2,y2 = Hammer_pos(lon,lat,opts)
    #x = (x1 + x2)/2
    #y = (y + y2)/2
    return x1,y
coordsl['Sheared Orthographic Bihemisphere'] = Orthosh_coords
posl['Sheared Orthographic Bihemisphere'] = Orthosh_pos
visl['Sheared Orthographic Bihemisphere'] = Circ_vis
ratl['Sheared Orthographic Bihemisphere'] = 2
wrapl['Sheared Orthographic Bihemisphere'] = 'xwrap'
syml['Sheared Orthographic Bihemisphere'] = 'sym4'



for p in proj_list:
    if p not in info_list:
        info_list[p] = ''
    if p not in visl:
        visl[p] = Def_vis
    if p not in ratl:
        ratl[p] = 2
    if p not in wrapl:
        wrapl[p] = 'none'
    if p not in syml:
        syml[p] = 'none'


    

#Rotates from one aspect to another; written almost exclusively by Mads de Silva

def to_cart(lat, lon):
    return np.array([np.cos(lon)*np.cos(lat), np.sin(lon)*np.cos(lat), np.sin(lat)])

def from_cart(point):
    lat = np.arcsin(point[2,:])
    lon = np.arctan2(point[1,:], point[0,:])
    return lat, lon

def construct_matrix(aspect):

    new_lon = aspect[0]
    new_lat = aspect[1]
    n_rot = aspect[2]

    #set up matrix -> constructed s.t. y-axis clockwise rotation -> z-axis anticlockwise rotation
    rotation = np.empty((3,3))

    #set matrix values
    rotation[0,0] = np.cos(new_lon)*np.cos(new_lat)
    rotation[0,1] = -np.sin(new_lon)*np.cos(n_rot) + np.sin(n_rot)*np.sin(new_lat)*np.cos(new_lon)
    rotation[0,2] = -np.sin(new_lat)*np.cos(new_lon)*np.cos(n_rot) - np.sin(n_rot)*np.sin(new_lon)

    rotation[1,0] = np.sin(new_lon)*np.cos(new_lat)
    rotation[1,1] = np.cos(new_lon)*np.cos(n_rot) + np.sin(n_rot)*np.sin(new_lat)*np.sin(new_lon)
    rotation[1,2] = -np.sin(new_lat)*np.sin(new_lon)*np.cos(n_rot) + np.sin(n_rot)*np.cos(new_lon)

    rotation[2,0] = np.sin(new_lat)
    rotation[2,1] = -np.sin(n_rot)*np.cos(new_lat)
    rotation[2,2] = np.cos(new_lat)*np.cos(n_rot)
    return rotation


def find_point(tar_lat, tar_lon, inverse):
    
    tar_point = to_cart(tar_lat, tar_lon)
    
    #apply inverse rotation matrix to target point to find in new basis
    tar_point_nb = np.matmul(inverse, tar_point)
    
    return from_cart(tar_point_nb)

def rev_find_point(nb_tar_lat, nb_tar_lon, rotation):

    tar_point = to_cart(nb_tar_lat, nb_tar_lon)

    #apply rotation matrix to target point to find in new basis
    tar_point_ob = np.matmul(rotation, tar_point)
    
    return from_cart(tar_point_ob)

def Rotate_from(lon, lat, aspect):

    if aspect[1] != 0 or aspect[2] != 0:    #Use rotation matrix only where necessary
        
        rotation = construct_matrix(aspect)
        latr, lonr = rev_find_point(lat.flatten(), lon.flatten(), rotation)
        lon = np.reshape(lonr, lon.shape)
        lat = np.reshape(latr, lat.shape)
    elif aspect[0] != 0:  #If only rotated by longitude, a simple frameshift can be used
        lon += aspect[0]
        if aspect[0] > 0:
            lon = np.where(lon > ma.pi, lon-2*ma.pi, lon)
        else:
            lon = np.where(lon < -ma.pi, lon+2*ma.pi, lon)

    return lon, lat

def Rotate_to(lon, lat, aspect):
    if aspect[1] != 0 or aspect[2] != 0:
        inverse = np.transpose(construct_matrix(aspect))    #set inverse matrix -> rotation matrix so inverse is its transpose
        latr, lonr = find_point(lat.flatten(), lon.flatten(), inverse)
        lon = np.reshape(lonr, lon.shape)
        lat = np.reshape(latr, lat.shape)
    elif aspect[0] != 0:
        lon -= aspect[0]
        if aspect[0] < 0:
            lon = np.where(lon > ma.pi, lon-2*ma.pi, lon)
        else:
            lon = np.where(lon < -ma.pi, lon+2*ma.pi, lon)
    
    return lon, lat

#Slices global map into smaller pieces for quicker projection
#proj_prev: previous projection applied to same set, to ensure compatible symmetry handling
def sym_slice(x,y, proj, proj_prev=None):
    lenx = 0
    leny = 0
    symproj = syml[proj]
    if proj_prev is not None:
        if symproj == 'sym4lat' and syml[proj_prev] != 'sym4lat':
            symproj = 'sym4'
        if symproj == 'sym4' and syml[proj_prev] not in ('sym4', 'sym4lat'):
            symproj = 'sym2'
        if symproj == 'sym2' and syml[proj_prev] not in ('sym2', 'sym4', 'sym4lat'):
            symproj = 'none'
    if symproj in ('symx', 'sym4', 'sym4lat'):
        lenx = x.shape[1]
        halfx = ma.ceil(lenx/2)   #preserve middle column for odd resolutions
        x = x[:,:halfx]     #slice in half on x axis
        y = y[:,:halfx]
        if symproj in ('sym4', 'sym4lat'):
            leny = y.shape[0]
            halfy = ma.ceil(leny/2)
            x = x[:halfy,:]     #slice in half on y axis
            if symproj == 'sym4lat':
                y = y[:halfy,0]          #extract single slice
                y = np.expand_dims(y,1)     #keep 2-dimensional for proper broadcasting
            else:
                y = y[:halfy,:]
    elif symproj != 'none':
        print(f'    Warning: invalid symmetry type "{symproj}" found for {proj}')
    return x, y, lenx, leny

#rejoin sliced map to original resolution
#proj_prev: previous projection applied to same set
def sym_join(x,y,proj,lenx,leny,proj_prev=None):

    
    symproj = syml[proj]
    if proj_prev is not None:
        if symproj == 'sym4lat' and syml[proj_prev] != 'sym4lat':
            symproj = 'sym4'
        if symproj == 'sym4' and syml[proj_prev] not in ('sym4', 'sym4lat'):
            symproj = 'sym2'
        if symproj == 'sym2' and syml[proj_prev] not in ('sym2', 'sym4', 'sym4lat'):
            symproj = 'none'

    if symproj in ('symx','sym4','sym4lat'):
        if symproj in ('sym4','sym4lat'):
            if symproj == 'sym4lat':
                halfx = ma.ceil(lenx/2)
                halfy = ma.ceil(leny/2)
                y = np.broadcast_to(y, (halfy,halfx))
            x2 = np.flip(x,0)
            y2 = np.flip(y,0)
            if leny % 2:
                x2 = x2[1:,:]
                y2 = y2[1:,:]
            x = np.concatenate((x,x2), 0)
            y = np.concatenate((y,-y2), 0)
        x2 = np.flip(x,1)
        y2 = np.flip(y,1)
        if lenx % 2:
            x2 = x2[:,1:]
            y2 = y2[:,1:]
        x = np.concatenate((x,-x2), 1)
        y = np.concatenate((y,y2), 1)
    return x,y


#Quicker vis function taking advantage of symmetry
def Quickvis(proj,x,y,lon=None,lat=None,opts=def_opts,get_far=False):
    visf = visl[proj]
    if visf == Def_vis:
        if get_far:
            return True, False
        return True
    if opts['use_sym']:
        x, y, lenx, leny = sym_slice(x,y, proj, proj_prev=None)
        if lat is not None:
            lon, lat, lenx, leny = sym_slice(lon,lat, proj, proj_prev=None)
    if get_far and proj in conicl:
        vis, far = visf(x,y,lon,lat,opts,get_far=True)
    else:
        vis = visf(x,y,lon,lat,opts)
        far = np.zeros((3,3))   #create dummy far array just so I don't have to keep checking if one is needed
    if opts['use_sym']:
        symproj = syml[proj]
        if symproj in ('symx','sym4','sym4lat'):
            if symproj in ('sym4','sym4lat'):
                vis2 = np.flip(vis,0)
                #far2 = np.flip(vis,0)      #shouldn't need these because all conics are symx, but leaving them in case that changes
                if leny % 2:
                    vis2 = vis2[1:,:]
                    #far2 = far2[1:,:]
                vis = np.concatenate((vis,vis2), 0)
                #far = np.concatenate((far,far2), 0)
            vis2 = np.flip(vis,1)
            far2 = np.flip(far,1)
            if lenx % 2:
                vis2 = vis2[:,1:]
                far2 = far2[:,1:]
            vis = np.concatenate((vis,vis2), 1)
            far = np.concatenate((far,far2), 1)
    if get_far:
        return vis, far
    return vis


#Finds the corresponding index in proj2 in aspect2 for every point in proj1 in aspect1
#Note: default backwards projection treats proj1 as output projection and proj2 as input projection
#   deg: interpret input aspects as degrees, otherwise as radians
#   get_coords: return lon and lat for both projections
#   use_sym: use symmetry to speed projection of global maps

def Find_index(x1, y1, proj1=0, proj2=0, aspect1=(0,0,0), aspect2=(0,0,0), opts=def_opts, deg=True, get_lon=False):
    proj1, proj2, opts, aspect1, aspect2 = Prep_pars(proj1, proj2, opts, aspect1, aspect2, deg)
    if opts['proj_direction'] == 'forward':
        print(" Working forwards from input to output projection")
        opts['in'] = True
    else:
        print(" Working backwards from output to input projection")
        opts['in'] = False
    coords = coordsl[proj1]
    pos = posl[proj2]
    print(f"  Determining lat/lon from {proj1} map...")
    if opts['use_sym']:
        x2, y2, lenx, leny = sym_slice(x1,y1,proj1)
    else:
        x2 = x1
        y2 = y1
    lon1, lat1 = coords(x2,y2,opts)
    if opts['use_sym']:
        lon1, lat1 = sym_join(lon1,lat1,proj1,lenx,leny)
    lon2 = np.copy(lon1)
    lat2 = np.copy(lat1)
    if np.max(np.abs(aspect1)) > 0:
        print(f"  Rotating from {'input' if opts['in'] else 'output'} orientation of {np.degrees(aspect1)}...")
        lon2, lat2 = Rotate_from(lon2, lat2, aspect1)
    opts['in'] = not opts['in']
    lon3 = np.copy(lon2)
    lat3 = np.copy(lat2)
    if np.max(np.abs(aspect2)) > 0:
        print(f"  Rotating to {'input' if opts['in'] else 'output'} orientation of {np.degrees(aspect2)}...")
        lon3, lat3 = Rotate_to(lon3, lat3, aspect2)
    print(f"  Determining corresponding position on {proj2} map...")
    lon4 = np.copy(lon3)
    lat4 = np.copy(lat3)
    lon4 = np.remainder(lon4 + ma.pi, 2*ma.pi) - ma.pi      #longitude loops
    if opts['use_sym'] and np.max(np.abs(aspect1)) == 0 and np.max(np.abs(aspect2)) == 0:
        lon4, lat4, lenx, leny = sym_slice(lon4,lat4,proj2,proj1)
    x2, y2 = pos(lon4,lat4,opts)
    if opts['use_sym'] and np.max(np.abs(aspect1)) == 0 and np.max(np.abs(aspect2)) == 0:
        x2, y2 = sym_join(x2,y2,proj2,lenx,leny,proj1)
    if get_lon:
        return x2, y2, lon1, lat1, lon2, lat2, lon3, lat3
    return x2, y2   


#Control functions

#Decide how an option value string should be interpreted and returns the appropriate data type
def interp_opt(opt_val):
    opt_type = 0
    if not isinstance(opt_val,str):
        return opt_val
    for c in opt_val:
        if c.isnumeric() or c=='-':     #treat as int if all numbers or minus sign
            pass
        elif c=='.':        #if also has a decimal, treat as float
            opt_type = 1
        else:               #if has other characters, treat as string
            opt_type = 2
            break
    if opt_type == 2:   #strings have a few special cases:
        if opt_val in ("none","None","NONE"):  #interpret as none type
            opt_key = None
        elif opt_val in ("true","True","TRUE"):   #interpret as boolean true
            opt_key = True
        elif opt_val in ("false","False","FALSE"):    #interpret as boolean false
            opt_key = False
        elif "(" in opt_val and "," in opt_val:        #attempt to interpret as tuple
            try:
                opt_key = tuple([int(v) for v in opt_val[1:-1].split(',')])
            except:
                opt_key = str(opt_val)      #interpret as string if tuple failed
        else:
            opt_key = str(opt_val)
    elif opt_type == 1:
        opt_key = float(opt_val)
    else:
        opt_key = int(opt_val)
    return opt_key

#load options from config file
# cfg_file: config file to load options from
def Load_options(cfg_file):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(cfg_file)
    in_opts = {}
    for section in cfg.sections():
        for k, v in cfg.items(section):
            in_opts.update({k: interp_opt(v)})
    return in_opts

#get options list, with combination of input and default options
def Get_opts(opts, def_opts=def_opts):
    if opts is None:
        opts = def_opts
    else:
        try:
            if opts['skip_update']:
                return opts
        except:
            pass
    opts_out = def_opts
    config = 'projpasta_options.cfg'
    if def_opts['alt_config'] is not None:
        config = def_opts['alt_config']
    if def_opts['alt_config'] is not None:
        config = opts['alt_config']
    if os.path.exists(config):
        print(f'Loading options from {config}...')
        opts_out.update(Load_options(config))
    opts_out.update(opts)
    opts['skip_update'] = True

    return opts_out

#get appropriate projection, allowing for string input
def Get_proj(proj):
    if isinstance(proj, int):
        return proj_list[proj]
    return proj

#Common function to process input parameters
def Prep_pars(proj_in, proj_out, opts, aspect_in=None, aspect_out=None, deg=False):
    opts = Get_opts(opts)
    proj_in = Get_proj(proj_in)
    proj_out = Get_proj(proj_out)
    if deg:
        aspect_in = np.radians(aspect_in)
        aspect_out = np.radians(aspect_out)
        return proj_in, proj_out, opts, aspect_in, aspect_out
    return proj_in, proj_out, opts

#pad edges of visible area to avoid seams on projection and interpolation
#   ptype: type of padding
#       rect: simply copy over left and right sides of array
#       xwrap: copy over left and right sides of non-rectangular visible area
#       reproj: reproject pixels from visible area
#   proj_in: projection of array
#   data_in: data array to pad
#   x_in: data x coords
#   y_in: data y coords
#   padding: pixels of padding to add
#   vis_f: vis function to use to determine visible area
#   opts: options for use with vis function
def pad(ptype, proj_in, vis_in, data_in, x_in, y_in, pd=2, opts=def_opts, far=None):
                
    x_in = np.concatenate((x_in[-pd:] - 2,       #add padding pixels to either side
                            x_in,
                            x_in[:pd] + 2))

    data_in1 = data_in[:,-pd:] if len(data_in.shape) == 2 else data_in[:,-pd:,:]
    if ptype == 'rect':
        data_in2 = data_in[:,:pd] if len(data_in.shape) == 2 else data_in[:,:pd,:]   #simple wrap-around for rectangular
    else:
        data_in1 = np.zeros_like(data_in1)
        data_in2 = data_in1
    data_in = np.concatenate((data_in1, data_in, data_in2), 1)
    if ptype in ('xwrap', 'conic_xwrap'):           #for more complex x-wrapping, mark out 2-pixel margin around input map
        
        vis1 = np.zeros_like(vis_in[:,:pd])
        vis_in = np.concatenate((vis1,vis_in,vis1),1)

        vis_n = np.logical_not(vis_in)

        visl = np.logical_and(np.roll(vis_in,-1,1), vis_n)
        visr = np.logical_and(np.roll(vis_in,1,1), vis_n)
        visu = np.logical_and(np.roll(vis_in,-1,0), vis_n)
        visd = np.logical_and(np.roll(vis_in,1,0), vis_n)

        if ptype == 'conic_xwrap' and far is not None: #special procedure for conic hemispheres
            hem = Get_azim_type(opts|{'in':True})
            if hem in ('hem','bihem') or proj_in in nonglobal:
                far1 = np.zeros_like(far[:,:pd])
                far = np.concatenate((far1,far,far1),1)
                if hem == 'bihem':
                    data_f = np.flip(data_in,0)
                    data_in[visr] = np.roll(data_f,1,1)[visr]   #wrap outside edges between each hemisphere
                    data_in[visl] = np.roll(data_f,-1,1)[visl]
                    data_in[visd] = np.roll(data_f,1,0)[visd]
                    data_in[visu] = np.roll(data_f,-1,0)[visu]

                farn = np.logical_not(far)
                visl = np.logical_and(visl,farn)    #make sure regular xwrapping doesn't affect outside edges
                visr = np.logical_and(visr,farn)
                visu = np.logical_and(visu,farn)
                visd = np.logical_and(visd,farn)

        data_f = np.flip(data_in,1)

        data_in[visd] = np.roll(data_f,1,0)[visd]
        data_in[visu] = np.roll(data_f,-1,0)[visu]
        data_in[visr] = np.roll(data_f,1,1)[visr]
        data_in[visl] = np.roll(data_f,-1,1)[visl]


    
    if ptype not in ('none', 'rect', 'xwrap', 'conic_xwrap'):
        print(f'    Warning: invalid wrap type "{ptype}" found for {proj_in}')

    #im = Image.fromarray(data_in)       #debug function
    #im.save('pad_test.png')
    
    return data_in, x_in, y_in


#find appropriate locations for continuous graticules
def Graticules(lon,lat,g_lon,g_lat,width=1):
    lon1 = np.remainder(lon + ma.pi, 2 * ma.pi) - ma.pi
    lat1 = np.remainder(lat + ma.pi/2, ma.pi) - ma.pi/2
    grats = np.zeros(lon.shape,dtype=bool)
    for g_l, l in zip((g_lon,g_lat),(lon1,lat1)):
        if g_l is not None:
            for g in g_l:
                g = ma.radians(g)
                if g == -ma.pi:       #ensure proper treatment of antimeridian
                    al = np.where(l > 0, l - 2*ma.pi, l)
                    diff = al - g
                else:
                    diff = l-g
                adiff = np.abs(diff)
                sdiff = np.sign(diff)
                for n in range(2):
                    over = sdiff + np.roll(sdiff,1,n)
                    over = np.where(np.abs(diff - np.roll(diff,1,n)) > ma.pi/2, 2, over)  #avoid issues with longitude wrapping
                    grats1 = over == 0
                    grats2 = grats1 & (adiff > np.roll(adiff,1,n))
                    grats1 = grats1 & np.logical_not(grats2)
                    n_grats = (sdiff == 0) | grats1 | np.roll(grats2,-1,n)
                    grats = grats | n_grats
    return grats




#Reproject data array
#   index: precomputed index(/ces), to avoid repeating work; generated if not provided
#   deg: interpret input aspects as degrees
#   get_index: return index(/ces) along with projected array
def Proj_Array(data_in, proj_in=0, proj_out=0, aspect_in=(0,0,0), aspect_out=(0,0,0), opts=def_opts,
               deg=True, index=None, get_index=False):
    proj_in, proj_out, opts = Prep_pars(proj_in, proj_out, opts)

    
    if opts['interp_type'] == 'none':
        opts['proj_direction'] = 'backward'
    elif opts['proj_direction'] == 'avoid_iter':
        if proj_out in iter_coords:
            if proj_in in iter_pos:
                opts['proj_direction'] = 'backward'
            else:
                opts['proj_direction'] = 'forward'
        else:
            opts['proj_direction'] = 'backward'
    
    if proj_in in aziml or conicl:
        hem = Get_azim_type(opts|{'in':True})
        if hem == 'hem' or (hem=='global' and proj_in in nonglobal):
            opts['global_in'] = False
            if proj_in in aziml:
                opts['avoid_seam'] = False

    if proj_in in nonglobal:
        opts['global_in'] = False



    shape_in = np.asarray(data_in.shape)
    if opts['is_crop_in'] is not None:
        crop = opts['is_crop_in']
        if len(crop) != 6:
            print('    Warning: input for is_crop_in has incorrect number of values, skipping uncrop')
            opts['is_crop_in'] = None
        else:
            opts['global_in'] = False
            if len(shape_in) > 2:
                data_full = np.zeros((crop[5],crop[4],shape_in[2]),dtype=data_in.dtype)
            else:
                data_full = np.zeros((crop[5],crop[4]),dtype=data_in.dtype)
            data_full[crop[2]:crop[3],crop[0]:crop[1]] = data_in[:]
            data_in = data_full
            shape_in[0] = crop[5]
            shape_in[1] = crop[4]

    x_in = np.linspace(-1,1,shape_in[1],False) + 1/shape_in[1]
    y_in = np.linspace(1,-1,shape_in[0],False) - 1/shape_in[0]

    if opts['force_ratio'] is not None:
        ratio = opts['force_ratio']
    else:
        ratio = ratl[proj_out]
        if callable(ratio):
            ratio = ratio(opts|{'in':False})
        if opts['relative_ratio']:
            rat_in = ratl[proj_in]
            if callable(rat_in):
                rat_in = rat_in(opts|{'in':True})
            ratio *= (shape_in[1]/shape_in[0]) / rat_in
    
    
    shape_out = np.copy(shape_in)
    if opts['force_scale'] is not None:
        shape_out1 = opts['force_scale']
        if isinstance(shape_out1, int):
            shape_out1 = [shape_out1, round(shape_out1/ratio)]
        shape_out[0] = shape_out1[1]
        shape_out[1] = shape_out1[0]
        
    else:
        rat_rel = ratio / (shape_in[1] / shape_in[0])
        if rat_rel > 1:
            shape_out[1] = round(shape_out[1] * rat_rel)
        elif rat_rel < 1:
            shape_out[0] = round(shape_out[0] / rat_rel)


    data_out = np.zeros((shape_out), data_in.dtype)
    
    x_out = np.linspace(-1,1,shape_out[1],False) + 1/shape_out[1]
    y_out = np.linspace(1,-1,shape_out[0],False) - 1/shape_out[0]

    if opts['crop_out'] is not None:
        crop = opts['crop_out']
        if len(crop) != 4:
            print('    Warning: input for crop_out has incorrect number of values, skipping crop')
            opts['crop_out'] = None
        else:
            opts['use_sym'] = False
            print('   Cropping output area;')
            print(f'    For future use of is_crop_in use ({crop[0]},{crop[1]},{crop[2]},{crop[3]},{shape_out[1]},{shape_out[0]})')
            data_out = data_out[crop[2]:crop[3],crop[0]:crop[1]]
            x_out = x_out[crop[0]:crop[1]]
            y_out = y_out[crop[2]:crop[3]]
    
    x_out, y_out = np.meshgrid(x_out, y_out)




    if opts['proj_direction'] == 'backward':

        if index is None:
            index = {}
        try:
            x_ind = index['x_ind']
            y_ind = index['y_ind']
        except:
            x_ind, y_ind, lon_out, lat_out, lon_tru, lat_tru, lon_in, lat_in = Find_index(x_out, y_out, proj_out, proj_in, aspect_out, aspect_in, opts, deg, get_lon=True)
            index['x_ind'] = x_ind
            index['y_ind'] = y_ind
        
        if opts['avoid_seam'] and opts['truncate_in']:
            print('  Padding input map to avoid seams...')
            ptype = wrapl[proj_in]
            try:
                vis_in = index['vis_in']
                far = index['far']
            except:
                x_in1, y_in1 = np.meshgrid(x_in,y_in)
                vis_in, far = Quickvis(proj_in,x_in1,y_in1,None,None,opts=opts|{'in':True},get_far=True)
                index['vis_in'] = vis_in
                index['far'] = far
            data_in, x_in, y_in = pad(ptype, proj_in, vis_in, data_in, x_in, y_in, opts=opts, far=far)
            if opts['is_crop_in'] is not None:      #allow for padding of partial input maps but erase anything outside input area
                crop = opts['is_crop_in']
                data_in1 = np.zeros_like(data_in)
                data_in1[crop[2]:crop[3],crop[0]:crop[1]] = data_in[crop[2]:crop[3],crop[0]:crop[1]]
                data_in = data_in1

        if opts['interp_type'] == 'none':
            print('  Copying data from nearest input pixels...')
            x_ind1 = np.clip(np.rint((x_ind+1)/2 * (data_in.shape[1]) - 0.5 + round((data_in.shape[1] - shape_in[1])/2)),0,data_in.shape[1]-1)
            y_ind1 = np.clip(np.rint((1-y_ind)/2 * (data_in.shape[0]) - 0.5 + round((data_in.shape[0] - shape_in[0])/2)),0,data_in.shape[0]-1)
            data_out = data_in[y_ind1.astype(np.int64),x_ind1.astype(np.int64)]
        
        else:
            print('  Interpolating data from input map...')
            from scipy.interpolate import interpn
            points = (y_in,x_in)
            xi = np.stack((y_ind, x_ind), -1)
            if len(data_out.shape) > 2:
                for n in range(data_out.shape[2]):
                    data_out[:,:,n] = interpn(points, data_in[:,:,n].astype(float), xi, method=opts['interp_type'], bounds_error=False, fill_value=None).astype(data_in.dtype)
            else:
                data_out[:] = interpn(points, data_in.astype(float), xi, method=opts['interp_type'], bounds_error=False, fill_value=None).astype(data_in.dtype)
        
        print(' Projection complete')

        #add graticules
        if opts['graticules'] is not None:
            print('  Adding graticules...')
            try:
                grat = index['grat']
            except:
                if opts['graticules'] == 'in':
                    g_lo = lon_in
                    g_la = lat_in
                elif opts['graticules'] == 'out':
                    g_lo = lon_out
                    g_la = lat_out
                elif opts['graticules'] == 'true':
                    g_lo = lon_tru
                    g_la = lat_tru
                else:
                    print('    Warning: invalid entry for graticules, marking none')
                    g_lo = lon_out
                    g_la = lat_out
                    opts['grat_lon'] == None    #graticules function still runs, but will produce not graticules
                    opts['grat_lat'] == None
                gr_lon = opts['grat_lon']       #constructing list of graticules to add
                gr_lat = opts['grat_lat']
                for n in range(2):
                    if n == 0:
                        gr = gr_lon
                    else:
                        gr = gr_lat
                    if gr == None:
                        continue
                    try:
                        gr = [ma.radians(g) for g in gr]
                    except:
                        if n == 0:
                            max = 180
                        else:
                            max = 90
                        gr = np.arange(-max, max, gr)
                        gr = gr.tolist()
                        if opts['truncate_out'] or wrapl[proj_out] == 'rect':   #various procedures to avoid graticule on edge of map
                            to_remove = []
                            hem = Get_azim_type(opts|{'in':True},no_glob = True if (proj_out in nonglobal) else False)
                            if proj_out == 'Lambert Conformal Conic' and (opts['azim_type_out'] == 'global' or opts['azim_type'] == 'global'):
                                hem = 'global'
                            if proj_out == 'Stereographic' and (opts['azim_type_out'] == 'global' or opts['azim_type'] == 'global'):
                                pass
                            elif proj_out in aziml and hem == 'global' and proj_out not in ('Nicolosi Globular','Ortelius Oval'):
                                if (n==0 and abs(aspect_out[1]) == 90 and
                                    (opts['graticules'] == 'true' or
                                    (opts['graticules'] == 'in' and abs(aspect_in[1]) == 90))):
                                    to_remove.append(-max)
                            elif (opts['graticules'] == 'out' or 
                                (np.max(np.abs(aspect_out)) == 0 and
                                (opts['graticules'] == 'true' or
                                (opts['graticules'] == 'in' and np.max(np.abs(aspect_in)) == 0)))):
                                to_remove.append(-max)
                                if ((n == 0 and proj_out in aziml and hem == 'bihem') or
                                    (n == 1 and proj_out in conicl and hem in ('hem','bihem'))):
                                    to_remove.append(0)
                                elif n == 0 and proj_out in aziml and hem == 'hem':
                                    to_remove.append(-90)
                                    to_remove.append(90)
                            elif np.max(np.abs(aspect_out[1:2])) == 0:
                                if opts['graticules'] == 'true':
                                    if n == 0:
                                        to_remove.append((aspect_out[0]) % 360 - 180)
                                        if proj_out in aziml:
                                            if hem == 'bihem':
                                                to_remove.append((aspect_out[0] + 180) % 360 - 180)
                                            elif hem == 'hem':
                                                to_remove.append((aspect_out[0] + 90) % 360 - 180)
                                                to_remove.append((aspect_out[0] + 270) % 360 - 180)
                                    else:
                                        to_remove.append(-max)
                                        if proj_out in conicl and hem in ('hem','bihem'):
                                            to_remove.append(0)
                                elif opts['graticules'] == 'in' and np.max(np.abs(aspect_in[1:2])) == 0:
                                    if n == 0:
                                        diff = aspect_out[0] - aspect_in[0]
                                        if diff > 0:
                                            to_remove.append((diff) % 360 - 180)
                                        if proj_out in aziml:
                                            if hem == 'bihem':
                                                to_remove.append((diff + 180) % 360 - 180)
                                            elif hem == 'hem':
                                                to_remove.append((diff + 90) % 360 - 180)
                                                to_remove.append((diff + 270) % 360 - 180)
                                    else:
                                        to_remove.append(-max)
                                        if proj_out in conicl and hem in ('hem','bihem'):
                                            to_remove.append(0)
                            for r in to_remove:
                                if r in gr:
                                    gr.remove(r)
                    if n == 0:
                        gr_lon = gr
                    else:
                        gr_lat = gr

                grat = Graticules(g_lo,g_la,gr_lon,gr_lat)
                index['grat'] = grat
        
            g_col = opts['grat_color']
            g_alpha = g_col[3] / 255
            g_col = np.asarray(g_col[0:3]).astype(float)
            if len(data_out.shape) > 2:
                for n in range(data_out.shape[2]):
                    data_out[grat] = (data_out[grat].astype(float) * (1-g_alpha) + g_col * g_alpha).astype(data_out.dtype)
            else:
                data_out[grat] = (data_out[grat].astype(float) * (1-g_alpha) + g_col[0] * g_alpha).astype(data_out.dtype)

        # choose areas to black out
        if ((opts['truncate_in'] and not opts['global_in']) |
            (opts['truncate_out'] and wrapl[proj_out] != 'rect') |
            (opts['crop_in'] is not None) |
            (opts['crop_out'] is not None) |
            (opts['crop_coords_true'] is not None) |
            (opts['crop_coords_in'] is not None) |
            (opts['crop_coords_out'] is not None)):
            print('  Trimming visible map area...')

        try:
            vis_all = index['vis_all']
        except:
            vis_all = np.ones_like(x_ind).astype(bool)
            vis_all = np.where((np.abs(x_ind) > 1) | (np.abs(y_ind) > 1), False, vis_all)   #areas outside input image
            if opts['truncate_in'] and not opts['global_in']:
                vis_all = np.where(visl[proj_in](x_ind, y_ind, None, None, opts=opts|{'in':True}), vis_all, False) #truncate input map
                if opts['is_crop_in']:
                    if opts['crop_in'] is None:
                        opts['crop_in'] = opts['is_crop_in'][:4]    #use crop_in method to limit to area of input
                    else:
                        crop = opts['crop_in']
                        ncrop = opts['is_crop_in']
                        for n in range(4):
                            if n in (0,2):
                                crop[n] = max(crop[n],ncrop[n])
                            else:
                                crop[n] = min(crop[n],ncrop[n])
                        opts['crop_in'] = crop

            if opts['truncate_out']:
                vis_all = np.where(Quickvis(proj_out, x_out, y_out, lon_out, lat_out, opts|{'in':False}), vis_all, False)   #truncate output map
            
            if opts['crop_in'] is not None:
                crop = opts['crop_in']
                if len(crop) != 4:
                    print('    Warning: input for crop_in has incorrect number of values, skipping crop')
                else:
                    lenx = len(x_in)
                    leny = len(y_in)
                    if crop[0] <= 0:
                        minx = -2
                    else:
                        minx = x_in[crop[0]+1] + 1/lenx
                    if crop[1] >= lenx:
                        maxx = 2
                    else:
                        maxx = x_in[crop[1]+1] + 1/lenx
                    if crop[2] <= 0:
                        maxy = 2
                    else:
                        maxy = y_in[crop[2]] + 1/leny
                    if crop[3] >= leny:
                        miny = -2
                    else:
                        miny = y_in[crop[3]] + 1/leny
                    vis_all = np.where((x_ind < minx) |
                                    (x_ind > maxx) |
                                    (y_ind < miny) |
                                    (y_ind > maxy),
                                    False, vis_all)
            for cropt in ('crop_coords_true','crop_coords_in','crop_coords_out'):
                crop = opts[cropt]
                if crop is None:
                    continue
                if len(crop) != 4:
                    print(f'    Warning: input for {cropt} has incorrect number of values, skipping crop')
                    continue
                if cropt == 'crop_coords_true':
                    lon = lon_tru
                    lat = lat_tru
                elif cropt == 'crop_coords_in':
                    lon = lon_in
                    lat = lat_in
                else:
                    lon = lon_out
                    lat = lat_out
                if crop[0] <= -180:
                    minlon = -1000
                else:
                    minlon = ma.radians(crop[0])
                if crop[1] >= 180:
                    maxlon = 1000
                else:
                    maxlon = ma.radians(crop[1])
                if crop[2] <= -90:
                    minlat = -1000
                else:
                    minlat = ma.radians(crop[2])
                if crop[3] >= 90:
                    maxlat = 1000
                else:
                    maxlat = ma.radians(crop[3])
                if minlon > maxlon:
                    vis_all = np.where((lon <= maxlon) | (lon >= minlon), vis_all, False)
                else:
                    vis_all = np.where((lon < minlon) | (lon > maxlon), False, vis_all)
                vis_all = np.where((lat < minlat) | (lat > maxlat), False, vis_all)
            index['vis_all'] = vis_all
        if len(data_out.shape) > 2:
            vis_all = np.expand_dims(vis_all, -1)
        
        data_out = np.where(vis_all, data_out, 0)

        if opts['crop_vis']:
            vis_r = np.nonzero(vis_all)
            y_min = np.min(vis_r[0])
            y_max = np.max(vis_r[0]+1)
            x_min = np.min(vis_r[1])
            x_max = np.max(vis_r[1]+1)
            if opts['crop_out'] is not None:
                crop = opts['crop_out']
                y_minr = y_min + crop[2]
                y_maxr = y_max + crop[2]
                x_minr = x_min + crop[0]
                x_maxr = x_max + crop[0]
            else:
                y_minr = y_min
                y_maxr = y_max
                x_minr = x_min
                x_maxr = x_max
            print(f'   Cropping to visible area')
            print(f'    For future use of is_crop_in use ({x_minr},{x_maxr},{y_minr},{y_maxr},{shape_out[1]},{shape_out[0]})')
            data_out = data_out[y_min:y_max,x_min:x_max]





            



    else:
        lon_in=None
        lat_out=None
        x_in, y_in = np.meshgrid(x_in, y_in)
        if index is None:
            index = {}
        try:
            x_ind = index['x_ind']
            y_ind = index['y_ind']
        except:
            x_ind, y_ind, lon_in, lat_in, lon_tru, lat_tru, lon_out, lat_out = Find_index(x_in, y_in, proj_in, proj_out, aspect_in, aspect_out, opts, deg, get_lon=True)
            index['x_ind'] = x_ind
            index['y_ind'] = y_ind
        try:
            vis_in = index['vis_in']
        except:
            vis_in = Quickvis(proj_in, x_in, y_in, lon_in, lat_in, opts|{'in':True})
            index['vis_in'] = vis_in

        x_ind = x_ind[vis_in]
        y_ind = y_ind[vis_in]
            
        x_ind = np.ravel(x_ind)
        y_ind = np.ravel(y_ind)



        if opts['avoid_seam']:
            try:
                extras = index['extras']
                x_ext = index['x_ext']
                y_ext = index['y_ext']
            except:
                extras = np.abs(lon_in) > ma.radians(170)
                lonext = lon_out[extras]
                lonext = np.where(lonext > 0, lonext - 2*ma.pi, lonext + 2*ma.pi)
                latext = lat_in[extras]
                x_ext, y_ext = posl[proj_out](lonext, latext, opts|{'in':False})
                x_ext = np.ravel(x_ext)
                y_ext = np.ravel(y_ext)
                index['extras'] = extras
                index['x_ext'] = x_ext
                index['y_ext'] = y_ext

            x_ind = np.concatenate((x_ind, x_ext))
            y_ind = np.concatenate((y_ind, y_ext))

        
        if len(data_in.shape) > 2:
            data_in1 = []
            for n in range(data_in.shape[2]):
                data_n = data_in[:,:,n]
                data_n = data_n[vis_in]
                data_in1.append(np.ravel(data_n))
        else:
            data_in1 = np.ravel(data_in[vis_in])
        
        if opts['avoid_seam']:
            if len(data_in.shape) > 2:
                for n in range(data_in.shape[2]):
                    datan = data_in[:,:,n]
                    data_in1[n] = np.concatenate((data_in1[n], datan[extras]))
            else:
                data_in1 = np.concatenate((data_in1, data_in[extras]))
        
        
        from scipy.interpolate import griddata

        points = (y_ind, x_ind)
        xi = np.stack((y_out, x_out), -1)
        if len(data_out.shape) > 2:
            for n in range(data_out.shape[2]):
                data_out[:,:,n] = griddata(points, data_in1[n].astype(float), xi, method=opts['interp_type']).astype(data_in.dtype)
        else:
            data_out[:] = griddata(points, data_in1.astype(float), xi, method=opts['interp_type']).astype(data_in.dtype)

        lon_out = None
        lat_out = None



    
        if opts['truncate_out']:
            try:
                vis_out = index['vis_out']
            except:
                vis_out = Quickvis(proj_out, x_out, y_out, lon_out, lat_out, opts|{'in':False})
                index['vis_out'] = vis_out
            if len(data_out.shape) > 2:
                vis_out = np.expand_dims(vis_out, -1)
            data_out = np.where(vis_out, data_out, 0)
    
        


    #g = Image.fromarray(grat)#.astype(np.uint8))
    #g.save('grat.png')

    print(' Operation complete')
    
    if get_index:
        return data_out, index
    return data_out




def Proj_Image(file_in, file_out='projp_out.png', proj_in=0, proj_out=0, aspect_in=(0,0,0), aspect_out=(0,0,0), opts=def_opts,
               deg=True, index=None, get_index=False):
    proj_in, proj_out, opts = Prep_pars(proj_in, proj_out, opts)


    map_in = Image.open(file_in)
    if map_in.mode == 'P' and opts['graticules'] != True:   #convert RGB to allow overlaying graticules
        map_in = map_in.convert('RGB')
    data_in = np.asarray(map_in)

    data_out = Proj_Array(data_in, proj_in, proj_out, aspect_in, aspect_out, opts,
                          deg=deg, index=index, get_index=get_index)
    print('  Saving image...')

    map_out = Image.fromarray(data_out, mode=map_in.mode)
    if map_in.mode == "P":
        map_out.putpalette(map_in.getpalette())
    map_out.save(file_out)
    print(' Map saved to '+file_out)

    if get_index:
        return index
    return




def Globe_gif(file_in, file_out='globe.gif', proj_in=0, init_aspect=(0,0,0), frames=36, duration=100, loop=0, prograde=True, opts=def_opts, deg=True):
    proj_in, proj_out, opts = Prep_pars(proj_in, "Orthographic", opts)
    map_in = Image.open(file_in)
    if map_in.mode == 'P' and opts['graticules'] != True:   #convert RGB to allow overlaying graticules
        map_in = map_in.convert('RGB')
    data_in = np.asarray(map_in)

    lon_step = 360/frames * -1 if prograde else 1
    ims = []
    aspect_out = list(init_aspect)
    opts['azim_type_out'] = 'hem'
    for f in range(frames):
        print(f'Frame {f}')
        dat = Proj_Array(data_in, proj_in, proj_out, (0,0,0), aspect_out, opts, deg)
        aspect_out[0] = aspect_out[0] + lon_step
        im = Image.fromarray(dat, mode=map_in.mode)
        if map_in.mode == "P":
            im.putpalette(map_in.getpalette())
        ims.append(im)
    ims[0].save(file_out, save_all=True, append_images=ims[1:],duration=duration,loop=loop)

    return




#Legacy functions to support pre-2.0 input functions

def Main(file_in, file_out, proj_in=0, proj_out=0, lon_in=0, lat_in=0, rot_in=0, lon_out=0, lat_out=0, rot_out=0,
         tol=1e-6, imax=20, hem_in=0, hem_out=0, trim=0, trunc=False, interp=0, aviter = False):
    MainDeg(file_in, file_out, proj_in, proj_out, ma.degrees(lon_in), ma.degrees(lat_in), ma.degrees(rot_in), ma.degrees(lon_out), ma.degrees(lat_out), ma.degrees(rot_out),
         tol, imax, hem_in, hem_out, trim, trunc, interp, aviter)
    return

#Runs the main function but directly takes angle inputs in degrees      
def MainDeg(file_in, file_out, proj_in=0, proj_out=0, lon_in=0, lat_in=0, rot_in=0, lon_out=0, lat_out=0, rot_out=0,
         tol=1e-6, imax=20, hem_in=0, hem_out=0, trim=0, trunc=False, interp=0, aviter = False):
    aspect_in = [lon_in, lat_in, rot_in]
    aspect_out = [lon_out, lat_out, rot_out]
    if hem_in == 0:
        hem_in = 'global_if'
    elif hem_in == 1:
        hem_in = 'hem'
    elif hem_in == 2:
        hem_in = 'bihem'
    if hem_out == 0:
        hem_out = 'global_if'
    elif hem_out == 1:
        hem_out = 'hem'
    elif hem_out == 2:
        hem_out = 'bihem'
    opts = {
        'tolerance': tol,
        'max_iter': imax,
        'azim_type_in': hem_in,
        'azim_type_out': hem_out,
        'truncate': trunc,
        'proj_direction': 'avoid_iter' if aviter else 'backward'
    }
    Proj_Image(file_in,file_out,proj_in,proj_out,aspect_in,aspect_out,opts,deg=True)
    return

#Input    


def Print_list(proj_list=proj_list,info_list=info_list,sections=sections):
    print('Projection Options and Codes:')
    
    maxlen = np.max([len(p) for p in proj_list])
    for p in range(len(proj_list)):
        pint = str(p)
        proj = proj_list[p]
        if proj in sections:
            endl = ''
            for i in range(80-len(sections[proj])):
                endl += '-'
            print(f'  ----  {sections[proj]}  {endl}')
        if len(pint) < 2:
            pint = ' ' + pint
        space = '  '
        for i in range(maxlen-len(proj)):
            space += ' '
        info = info_list[proj]
        if False:            #colors for windows command line, only used for .exe version
            os.system('color')
            if 'equidistant' in info:
                info = info.replace('equidistant','\033[36mequidistant\033[0m')
            if 'equal-area' in info:
                info = info.replace('equal-area','\033[32mequal-area\033[0m')
            if 'conformal' in info:
                info = info.replace('conformal','\033[33mconformal\033[0m')
        print(f'{pint}: {proj}{space}{info}')
       # print(f'     {info_list[p]}')
    return

#Standard input routine
# prompt: text prompt on screen
# 
def Inprompt(prompt,f,list=None):
    inp = input(prompt)
    while True:
        try:
            var = f(inp)
            if list is not None:
                var = list[var]
            return var
        except:
            inp = input(" Invalid input, try again: ")


if __name__ == "__main__":
    print(f'''Projection Pasta version {ver_num}
For Reprojection of maps between arbitrary aspects
Made 2023 by Mads de Silva and Nikolai Hersfeldt
''')
    opts = Get_opts(None)
    if opts['skip_setup']:
        print("Skipping setup per 'skip_setup' option")
        file_in = opts['file_in']
        file_out = opts['file_out']
        proj_in = opts['proj_in']
        proj_out = opts['proj_out']
        aspect_in = opts['aspect_in']
        aspect_out = opts['aspect_out']

    else:
        print('')
        Print_list()
        for i in range(2):
            print("Input Image" if i==0 else "Output Image")
            while True:
                file = input(" Filename: ")
                if i==1 or os.path.exists(file):
                    break
                print("  No file found at "+str(file))
            proj = Inprompt(" Projection: ",int,proj_list)
            if proj in aziml:
                print(" Projection subtype: ")
                if proj in nonglobal:
                    if proj in optlists_global:
                        print("  0: large map of custom size")
                else:
                    print("  0: single global map")
                print("  1: single hemisphere")
                print("  2: bihemisphere; 2 maps of opposite hemispheres")
                hem = Inprompt("  Choose subtype: ", int, ['global','hem','bihem'])
                opts['azim_type_in' if i==0 else 'azim_type_out'] = hem
            prev_ref = None
            for ol in (optlists_global, optlists):
                if proj in ol:
                    if ol == optlists_global and opts['azim_type_in' if i==0 else 'azim_type_out'] != 'global':
                        continue
                    optlist = ol[proj]
                    optkeys = list(optlist.keys())
                    print(' ' + optkeys[0])
                    if len(optkeys) == 1:
                        optchoose = 0
                    else:
                        print('  0: Enter custom value')
                        for k in range(1,len(optkeys)):
                            print(f'  {k}: {optkeys[k]}')
                        optchoose = Inprompt("  Choose option: ",int,[n for n in range(len(optkeys))])
                    if optchoose == 0:
                        numref = optlist[optkeys[0]]
                        if numref == 1:
                            cref = Inprompt("  Input custom value (in degrees): ",float)
                        else:
                            cref = []
                            for n in range(numref):
                                cref.append(Inprompt(f"  Input custom value {n+1} (in degrees): ",float))
                    else:
                        cref = optlist[optkeys[optchoose]]
                    if prev_ref is not None:
                        ref = [prev_ref]
                        for n in cref:
                            ref.append(n)
                    else:
                        ref = cref
                    prev_ref = ref
                    opts['ref_in' if i==0 else 'ref_out'] = ref
            lon = Inprompt(" Center longitude (-180 to 180): ",float)
            lat = Inprompt(" Center latitude (-90 to 90): ",float)
            rot = Inprompt(" Clockwise rotation from north (0 to 360): ",float)
            print('')

            if i==0:
                file_in = file
                proj_in = proj
                aspect_in = (lon,lat,rot)
            else:
                file_out = file
                proj_out = proj
                aspect_out = (lon,lat,rot)

    print("""
Working...""")

    Proj_Image(file_in,file_out,proj_in,proj_out,aspect_in,aspect_out,opts,deg=True)

    z = input("Press enter to close")


