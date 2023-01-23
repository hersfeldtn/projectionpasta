from PIL import Image
import numpy as np
import os
import math as ma

#Remapping functions for each projection:
#coords determines (lon,lat) coordinates for a given (x,y) position on the map
#pos determines (x,y) position for given (lon,lat) coordinates
#vis determines if an (x,y) position on the output should be shown on non-rectangular maps
#pres determines if an (x,y) position is present on the input
#typ incorporates extra information required for specific projection types
#rat is the map width/height ratio
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

def der(lon,lat,t,f):   #Approximate partial derivatives by the secant method
    x1,y1 = f(lon,lat+t,t)
    x2,y2 = f(lon,lat-t,t)
    x3,y3 = f(lon+t,lat,t)
    x4,y4 = f(lon-t,lat,t)
    dxla = (x1-x2)/(2*t)
    dyla = (y1-y2)/(2*t)
    dxlo = (x3-x4)/(2*t)
    dylo = (y3-y4)/(2*t)
    return dxla,dxlo,dyla,dylo
def Inv_coords(x,y,t,imax,f,f1,vis,ex): #General iterative inverse function from Bildirici 2016: https://doi.org/10.1080/15230406.2016.1200492
    print("  (using iterative method, may take a bit)")
    lon,lat = f1(x,y,ex)
    vis1 = vis(x,y,ex)
    i=1
    while True:
        x1,y1 = f(lon,lat,t)
        dxla,dxlo,dyla,dylo = der(lon,lat,t,f)
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

def Onehem_coords(x,y,ex,f,s):    #Set of general functions for maps that can be mapped as hemispheres
    return f(x*s,y*s,ex)
def Bihem_coords(x,y,ex,f,s):
    lon,lat=Onehem_coords(np.where(x>0,2*x-1,2*x+1),y,ex,f,s)
    lon = np.where(x>0,lon+ma.pi/2,lon-ma.pi/2)
    return lon,lat
def Hem_coords(x,y,hem,ex,f,s):
    if hem == 0:
        return f(x,y,ex)
    elif hem == 1:
        return Onehem_coords(x,y,ex,f,s)
    elif hem == 2:
        return Bihem_coords(x,y,ex,f,s)
    else:
        raise Exception("Invalid output map subtype selection (must be 0, 1, or 2)")

def Onehem_pos(lon,lat,ex,f,s):
    x,y = f(lon,lat,ex)
    return x/s, y/s
def Bihem_pos(lon,lat,ex,f,s):
    lon1 = np.where(lon>0,lon,lon+ma.pi)-ma.pi/2
    x,y = Onehem_pos(lon1,lat,ex,f,s)
    x /= 2
    return np.where(lon>0,x+0.5,x-0.5),y
def Hem_pos(lon,lat,hem,ex,f,s):
    if hem == 0:
        return f(lon,lat,ex)
    elif hem == 1:
        return Onehem_pos(lon,lat,ex,f,s)
    elif hem == 2:
        return Bihem_pos(lon,lat,ex,f,s)
    else:
        raise Exception("Invalid input map subtype selection (must be 0, 1, or 2)")    

def Def_vis(x,y,ex):
    return True
def Edge_vis(x,lat,f):
    xmax, y1 = f(ma.pi,np.where(np.abs(lat)>ma.pi/2,ma.pi/2,lat),ex)
    return np.where(np.abs(x) > xmax, False, True)
def Circ_vis(x,y,ex):
    return np.where(np.sqrt(x**2 + y**2) > 1, False, True)
def Hem_vis(x,y,hem):
    if hem == 0 or hem == 1:
        return Circ_vis(x,y,0)
    elif hem == 2:
        return Circ_vis(1-abs(x*2),y,0)

def Def_pres(x,y,ex):
    return True
def Hem_pres(x,y,hem):
    if hem == 0 or hem == 2:
        return True
    elif hem == 1:
        return np.where(np.sqrt(x**2 + y**2) > 1, False, True)

def Def_typ(tol,imax,hem):
    return 0
def Hemfull_typ(tol,imax,hem):
    return hem
def Hemonly_typ(tol,imax,hem):
    if hem == 0:
        hem = 1
    return hem
def Iter_typ(tol,imax,hem):
    return (tol,imax)
def HemfullIter_typ(tol,imax,hem):
    return (tol,imax,hem)
    

#Then, functions by projection
def Equi_coords(x,y,ex):
    lon = ma.pi*x
    lat = ma.pi*y/2
    return lon,lat
def Equi_pos(lon,lat,ex):
    x = lon/ma.pi
    y = 2*lat/ma.pi
    return x,y
Equi_vis = Def_vis
Equi_pres = Def_pres
Equi_typ = Def_typ
Equi_rat = 2

def Sin_coords(x,y,ex):
    lat = ma.pi*y/2
    lon = ma.pi*x/np.cos(lat)
    return lon,lat
def Sin_pos(lon,lat,ex):
    x = np.cos(lat)*lon/ma.pi
    y = 2*lat/ma.pi
    return x,y
def Sin_vis(x,y,lat):
    return np.where(abs(x) > np.cos(ma.pi*y/2), False, True)
Sin_pres = Def_pres
Sin_typ = Def_typ
Sin_rat = 2

def Moll_coords(x,y,ex):
    th = np.arcsin(y)
    lat = np.arcsin((2*th + np.sin(2*th))/ma.pi)
    lon = np.fmod(ma.pi*x/(np.cos(th)), ma.pi)
    return lon,lat
def Moll_pos(lon,lat,ex):
    print("  (using iterative method, may take a bit)")
    t = ex[0]
    imax= ex[1]
    i=1
    th = lat
    while np.amax(np.abs(2*th+np.sin(2*th)/ma.pi*np.sin(lat)-1)) > t:   #no closed-form solution, so iterate until maximum error is < tolerance
        th = th - (2*th + np.sin(2*th) - ma.pi*np.sin(lat)) / (2 + 2*np.cos(2*th))
        i+=1
        if i>imax:
            print("  Reached maximum of "+str(imax)+" iterations without converging, outputting result")
            break
    x = lon*np.cos(th)/ma.pi
    y = np.sin(th)
    return x,y
Moll_vis = Circ_vis
Moll_pres = Def_pres
Moll_typ = Iter_typ
Moll_rat = 2

def Hammer_coords(x,y,ex):
    x1 = x*ma.sqrt(2)*2
    y1 = y*ma.sqrt(2)
    z = np.sqrt(1 - (x1/4)**2 - (y1/2)**2)
    lon = 2*np.arctan( z*x1 / (2*(2*z**2-1)) )
    lat = np.arcsin(z*y1)
    return lon,lat
def Hammer_pos(lon,lat,ex):
    x = np.cos(lat)*np.sin(lon/2) / np.sqrt(1+np.cos(lat)*np.cos(lon/2))
    y = np.sin(lat) / np.sqrt(1+np.cos(lat)*np.cos(lon/2))
    return x,y
Hammer_vis = Circ_vis
Hammer_pres = Def_pres
Hammer_typ = Def_typ
Hammer_rat = 2

def Ait_guess(x,y,ex):  #Special routine for Ait_coords initial guess
    return np.where(Circ_vis(x,y,ex),Hammer_coords(x,y,ex),Wag_coords(x,y,ex))
def Ait_coords(x,y,ex):
    return Inv_coords(x,y,ex[0],ex[1],Ait_pos,Ait_guess,Ait_vis,0)
def Ait_pos(lon,lat,ex):
    al = ma.pi/2*np.sinc(np.arccos(np.cos(lat)*np.cos(lon/2))/ma.pi)
    x = np.cos(lat)*np.sin(lon/2)/al
    y = np.sin(lat)/al
    return x,y
Ait_vis = Circ_vis
Ait_pres = Def_pres
Ait_typ = Iter_typ
Ait_rat = 2

def Wink_coords(x,y,ex):
    return Inv_coords(x,y,ex[0],ex[1],Wink_pos,Wag_coords,Wink_vis,0)
def Wink_pos(lon,lat,ex):
    al = np.arccos(np.cos(lat)*np.cos(lon/2))/ma.pi
    x = (lon/ma.pi + np.cos(lat)*np.sin(lon/2) / np.sinc(al))/(1+ma.pi/2)
    y = (lat + np.sin(lat) / np.sinc(al))/ma.pi
    return x,y
def Wink_vis(x,y,ex):
    lat1 = np.arccos(2*(np.abs(x)*(1+ma.pi/2)-1)/ma.pi)
    ymax = (lat1 + ma.pi/2 * np.sin(lat1))/ma.pi
    return np.where (np.abs(y) > ymax,False,True)
Wink_pres = Def_pres
Wink_typ = Iter_typ
Wink_rat = (1+ma.pi/2)/(ma.pi/2)

def Wag_coords(x,y,ex):
    lat = y*ma.pi/2
    ph = np.arcsin(lat*ma.sqrt(3)/ma.pi)
    lon = x*ma.pi/np.cos(ph)
    return lon,lat
def Wag_pos(lon,lat,ex):
    y = lat*2/ma.pi
    x = lon/ma.pi*np.sqrt(1-3*(lat/ma.pi)**2)
    return x,y
def Wag_vis(x,y,ex):
    return np.where(abs(x) > np.sqrt(1-3*(y/2)**2),False,True)
Wag_pres = Def_pres
Wag_typ = Def_typ
Wag_rat = 2

Kav_coords = Wag_coords
Kav_pos = Wag_pos
Kav_vis = Wag_vis
Kav_pres = Wag_pres
Kav_typ = Wag_typ
Kav_rat = ma.sqrt(3)

def Ort_coords(x,y,ex):
    lon1,lat1 = Inv_coords(x,y,ex[0],ex[1],Ort_pos,Hammer_coords,Circ_vis,0)
    lat = y*ma.pi/2
    lon2 = np.abs(x)*ma.pi + ma.pi/2 - np.sqrt(ma.pi**2/4 - lat**2)
    lon2 = np.where(x>0, lon2, -lon2)
    lon = np.where(np.sqrt((x*2)**2+y**2)>1,lon2,lon1)
    return lon,lat
def Ort_pos(lon,lat,ex):
    y = lat*2/ma.pi
    abslon = np.abs(lon)
    F = (ma.pi**2/(4*abslon)+abslon)/2
    x = np.where(abslon>ma.pi/2,
                 np.sqrt(ma.pi**2/4 - lat**2) + abslon - ma.pi/2,
                 abslon - F + np.sqrt(F**2 - (y*ma.pi/2)**2))
    x = np.where(lon>0, x/ma.pi, -x/ma.pi)
    return x,y
def Ort_vis(x,y,ex):
    return np.where(np.abs(x)>1/2,Hem_vis(x,y,2),True)
Ort_pres = Def_pres
Ort_typ = Iter_typ
Ort_rat = 2

def Nic_coords1(x,y,ex):
    return Inv_coords(x,y,ex[0],ex[1],Nic_pos,Azim_coords,Nic_vis,ex[2])
def Nic_pos1(lon,lat,ex):
    lat0 = np.where(lat == 0, True, False)
    lon0 = np.where(lon == 0, True, False)
    latmax = np.where(np.abs(lat) == ma.pi/2, True, False)
    lonmax = np.where(np.abs(lon) == ma.pi/2, True, False)
    lat1 = np.where(lat0 | latmax, lat+1e-4, lat)
    lon1 = np.where(lon0 | latmax, lon+1e-4, lon)
    sinla = np.sin(lat1)
    b = ma.pi/(2*lon1) - 2*lon1/ma.pi
    c = 2*lat1/ma.pi
    d = (1-c**2) / (sinla - c)
    M = (b*sinla/d - b/2) / (1 + b**2/d**2)
    N = (d**2*sinla/b**2 + d/2) / (1 + d**2/b**2)
    x1 = np.sqrt(M**2 + np.cos(lat1)**2 / (1 + b**2/d**2))
    x = np.where(lon>0, M+x1, M-x1)
    y1 = np.sqrt(N**2 - (d**2/b**2 * sinla**2 + d*sinla - 1) / (1 + d**2/b**2))
    y = np.where(lat>0, N+y1, N-y1)
    x = np.where(lon0 | latmax,0,x)
    x = np.where(lat0 | lonmax,np.cos(lat)*lon,x)/2
    y = np.where(lon0 | lat0,lat,y)
    y = np.where(lonmax,np.sin(lat)/2,y)
    y = np.where(latmax,lat,y)/2
    return x,y
def Nic_coords(x,y,ex):
    return Hem_coords(x,y,ex[2],ex,Nic_coords1,1/2)
def Nic_pos(x,y,ex):
    return Hem_pos(x,y,ex[2],ex,Nic_pos1,1/2)
Nic_vis = Hem_vis
def Nic_pres(x,y,ex):
    return Hem_pres(x,y,ex[2])
Nic_typ = HemfullIter_typ
Nic_rat = 1

def Azim_coords1(x,y,ex):
    rh = np.sqrt(x**2+y**2)
    th = np.arctan2(x,-y)
    return From_pol(rh,th)
def Azim_pos1(lon,lat,ex):
    rh,th = To_pol(lon,lat)
    x = rh*np.sin(th)
    y = rh*np.cos(th)
    return x,y
def Azim_coords(x,y,hem):
    return Hem_coords(x,y,hem,0,Azim_coords1,1/2)
def Azim_pos(x,y,hem):
    return Hem_pos(x,y,hem,0,Azim_pos1,1/2)
Azim_vis = Hem_vis
Azim_pres = Hem_pres
Azim_typ = Hemfull_typ
Azim_rat = 1

def Ortho_coords1(x,y,ex):
    rh = np.sqrt((2*x)**2+(2*y)**2)
    c = np.arcsin(np.fmod(rh,1))
    lat = np.arcsin(2*y*np.sin(c)/rh)
    lon = np.arctan(2*x*np.sin(c)/(rh*np.cos(c)))
    return lon,lat
def Ortho_pos1(lon,lat,ex):
    x = np.cos(lat)*np.sin(lon)/2
    y = np.sin(lat)/2
    return x,y
def Ortho_coords(x,y,hem):
    return Hem_coords(x,y,hem,0,Ortho_coords1,1/2)
def Ortho_pos(x,y,hem):
    return Hem_pos(x,y,hem,0,Ortho_pos1,1/2)
Ortho_vis = Hem_vis
Ortho_pres = Hem_pres
Ortho_typ = Hemonly_typ
Ortho_rat = 1

def Stereo_coords1(x,y,ex):
    r = np.sqrt(x**2+y**2)
    th = np.arctan2(x,-y)
    rh = (np.arctan(2*r) - ma.pi/4)*2/ma.pi + 1/2
    return From_pol(rh,th)
def Stereo_pos1(lon,lat,ex):
    rh,th = To_pol(lon,lat)
    r = np.tan(ma.pi/4 + ma.pi*(rh-1/2)/2)/2
    x = r*np.sin(th)
    y = r*np.cos(th)
    return x,y
def Stereo_coords(x,y,hem):
    return Hem_coords(x,y,hem,0,Stereo_coords1,1/2)
def Stereo_pos(x,y,hem):
    return Hem_pos(x,y,hem,0,Stereo_pos1,1/2)
Stereo_vis = Hem_vis
Stereo_pres = Hem_pres
Stereo_typ = Hemonly_typ
Stereo_rat = 1

def Lamb_coords1(x,y,ex):
    r = np.sqrt(x**2+y**2)
    th = np.arctan2(x,-y)
    rh = np.arcsin(np.fmod(r,1))*2/ma.pi
    return From_pol(rh,th)
def Lamb_pos1(lon,lat,ex):
    rh,th = To_pol(lon,lat)
    r = np.sin(rh*ma.pi/2)
    x = r*np.sin(th)
    y = r*np.cos(th)
    return x,y
def Lamb_coords(x,y,hem):
    return Hem_coords(x,y,hem,0,Lamb_coords1,ma.sqrt(2)/2)
def Lamb_pos(x,y,hem):
    return Hem_pos(x,y,hem,0,Lamb_pos1,ma.sqrt(2)/2)
Lamb_vis = Hem_vis
Lamb_pres = Hem_pres
Lamb_typ = Hemfull_typ
Lamb_rat = 1
    

def Merc_coords(x,y,ex):
    lon = ma.pi*x
    lat = y/abs(y)*(2*np.arctan(np.exp(abs(y)*ma.pi))-ma.pi/2)
    return lon,lat
def Merc_pos(lon,lat,ex):
    x = lon/ma.pi
    y = lat/abs(lat)*np.log(np.tan(ma.pi/4 + abs(lat)/2))/ma.pi
    return x,y
Merc_vis = Def_vis
def Merc_pres(x,y,ex):
    return np.where(abs(y) > 1, False, True)
Merc_typ = Def_typ
Merc_rat = 1

def Gallst_coords(x,y,ex):
    lon = x*ma.pi
    lat = 2*np.arctan(y)
    return lon,lat
def Gallst_pos(lon,lat,ex):
    x = lon/ma.pi
    y = np.tan(lat/2)
    return x,y
Gallst_vis = Def_vis
Gallst_pres = Def_pres
Gallst_typ = Def_typ
Gallst_rat = (ma.pi/ma.sqrt(2))/(1+ma.sqrt(2)/2)

def Mill_coords(x,y,ex):
    lon = x*ma.pi
    y1 = y*ma.asinh(ma.tan(ma.pi*2/5))
    lat = 5/4*np.arctan(np.sinh(y1))
    return lon,lat
def Mill_pos(lon,lat,ex):
    x = lon/ma.pi
    y1 = np.arcsinh(np.tan(lat*4/5))
    y = y1/ma.asinh(ma.tan(ma.pi*2/5))
    return x,y
Mill_vis = Def_vis
Mill_pres = Def_pres
Mill_typ = Def_typ
Mill_rat = ma.pi/(5/4*ma.asinh(ma.tan(ma.pi*2/5)))

#Lists referring to the functions in same order for easy reference
#Nicolosi and Lambert still buggy
coordsl = [
    Equi_coords,
    Sin_coords,
    Moll_coords,
    Hammer_coords,
    Ait_coords,
    Wink_coords,
    Kav_coords,
    Wag_coords,
    Ort_coords,
    #Nic_coords,
    Azim_coords,
    Ortho_coords,
    Stereo_coords,
    #Lamb_coords,
    Merc_coords,
    Gallst_coords,
    Mill_coords
    ]

posl = [
    Equi_pos,
    Sin_pos,
    Moll_pos,
    Hammer_pos,
    Ait_pos,
    Wink_pos,
    Kav_pos,
    Wag_pos,
    Ort_pos,
    #Nic_pos,
    Azim_pos,
    Ortho_pos,
    Stereo_pos,
    #Lamb_pos,
    Merc_pos,
    Gallst_pos,
    Mill_pos
    ]

visl = [
    Equi_vis,
    Sin_vis,
    Moll_vis,
    Hammer_vis,
    Ait_vis,
    Wink_vis,
    Kav_vis,
    Wag_vis,
    Ort_vis,
    #Nic_vis,
    Azim_vis,
    Ortho_vis,
    Stereo_vis,
    #Lamb_vis,
    Merc_vis,
    Gallst_vis,
    Mill_vis
    ]

presl = [
    Equi_pres,
    Sin_pres,
    Moll_pres,
    Hammer_pres,
    Ait_pres,
    Wink_pres,
    Kav_pres,
    Wag_pres,
    Ort_pres,
    #Nic_pres,
    Azim_pres,
    Ortho_pres,
    Stereo_pres,
    #Lamb_pres,
    Merc_pres,
    Gallst_pres,
    Mill_pres
    ]

ratl = [
    Equi_rat,
    Sin_rat,
    Moll_rat,
    Hammer_rat,
    Ait_rat,
    Wink_rat,
    Kav_rat,
    Wag_rat,
    Ort_rat,
    #Nic_rat,
    Azim_rat,
    Ortho_rat,
    Stereo_rat,
    #Lamb_rat,
    Merc_rat,
    Gallst_rat,
    Mill_rat
    ]

typl = [
    Equi_typ,
    Sin_typ,
    Moll_typ,
    Hammer_typ,
    Ait_typ,
    Wink_typ,
    Kav_typ,
    Wag_typ,
    Ort_typ,
    #Nic_typ,
    Azim_typ,
    Ortho_typ,
    Stereo_typ,
    #Lamb_typ,
    Merc_typ,
    Gallst_typ,
    Mill_typ
    ]

#Rotates from one aspect to another; written almost exclusively by Amadea de Silva

def to_cart(lat, lon):
    return np.array([np.cos(lon)*np.cos(lat), np.sin(lon)*np.cos(lat), np.sin(lat)])

def from_cart(point):
    lat = np.arcsin(point[2])
    lon = np.arctan2(point[1], point[0])
    return lat, lon

def construct_matrix(new_lat, new_lon, n_rot):
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

def Rotate(lon, lat, lon_in, lat_in, rot_in, lon_out, lat_out, rot_out):

    if lat_out != 0 or rot_out != 0:    #Use rotation matrix only where necessary
        print(" Rotating from output orientation...")
        rotation = construct_matrix(lat_out, lon_out, rot_out)
        with np.nditer([lon, lat], op_flags=[['readwrite'],['readwrite']]) as co:
            for lo, la in co:
                la[...], lo[...] = rev_find_point(la, lo, rotation) 
    elif lon_out != 0:  #If only rotated by longitude, a simple frameshift can be used
        print(" Rotating from output orientation...")
        lon -= lon_out
        if lon_out < 0:
            lon = np.where(lon > ma.pi, lon-2*ma.pi, lon)
        else:
            lon = np.where(lon < -ma.pi, lon+2*ma.pi, lon)

    if lat_in != 0 or rot_in != 0:
        print(" Rotating to input orientation...")
        inverse = np.transpose(construct_matrix(lat_in, lon_in, rot_in))    #set inverse matrix -> rotation matrix so inverse is its transpose
        with np.nditer([lon, lat], op_flags=[['readwrite'],['readwrite']]) as co:
            for lo, la in co:
                la[...], lo[...] = find_point(la, lo, inverse)
    elif lon_in != 0:
        print(" Rotating to input orientation...")
        lon += lon_in
        if lon_in < 0:
            lon = np.where(lon > ma.pi, lon-2*ma.pi, lon)
        else:
            lon = np.where(lon < -ma.pi, lon+2*ma.pi, lon)
    
    return lon, lat
    

#Finds the corresponding pixel index in input for every index in the output

def Find_index(x1, y1, coords, pos, lon_in, lat_in, rot_in, ex_in, lon_out, lat_out, rot_out, ex_out):
    print(" Determining lat/lon of output pixels...")
    lon1, lat1 = coords(x1,y1,ex_out)
    lon2, lat2 = Rotate(lon1, lat1, lon_in, lat_in, rot_in, lon_out, lat_out, rot_out)
    print(" Determining corresponding input pixels...")
    x2, y2 = pos(lon2,lat2,ex_in)
    return y2, x2

#Main function

def Main(file_in, file_out, proj_in=0, proj_out=0, lon_in=0, lat_in=0, rot_in=0, lon_out=0, lat_out=0, rot_out=0,
         tol=1e-6, imax=10, hem_in=0, hem_out=0, trunc=1):
    coords = coordsl[proj_out]
    pos = posl[proj_in]
    if trunc == 1:
        vis = visl[proj_out]
    else:
        vis = Def_vis
    pres = presl[proj_in]
    ex_in = typl[proj_in](tol,imax,hem_in)
    ex_out = typl[proj_out](tol,imax,hem_out)
    rat = ratl[proj_in]/ratl[proj_out]
    if hem_in == 2:
        rat *= 2
    if hem_out == 2:
        rat /= 2
    map_in = Image.open(file_in)
    mapw = map_in.width
    maph = map_in.height
    mapw_out = mapw
    maph_out = maph
    if rat > 1:
        maph_out = int(maph*rat)    #increase rather than decrease resolution as necessary to reach new aspect ratio
    elif rat < 1:
        mapw_out = int(mapw/rat)
    tol = min(tol, 0.1/mapw, 0.1/maph)

    data_in = np.asarray(map_in) 
    
    #if data_in.ndim > 2:
    #    blank = np.zeros(data_in[0,0].shape)    #Create blank color appropriate to image color mode
    #else:
    #    blank = 0
    x_out = np.linspace(-1+1/mapw_out,1-1/mapw_out,mapw_out)    #All maps are treated internally as squares
    y_out = np.linspace(1-1/maph_out,-1+1/maph_out,maph_out)
    x_out,y_out = np.meshgrid(x_out,y_out)
    y_in,x_in = Find_index(x_out, y_out, coords, pos, lon_in, lat_in, rot_in, ex_in, lon_out, lat_out, rot_out, ex_out)
    print(" Remapping pixels...")
    x_in2 = np.where(np.abs(x_in)>1, (np.abs(x_in-1) % 2 - 1)* np.where(x_in>0,1,-1), x_in)  #x wraps around
    y_in2 = np.where(np.abs(y_in)>1, np.abs(np.abs(y_in-1)%4-2)-1, y_in)    #y inverts
    #print(np.amax(x_in))
    #print(np.amax(y_in))
    vis1 = vis(x_out,y_out,ex_out)
    pres1 = pres(x_in,y_in,ex_in)
    if data_in.ndim > 2:    #various procedures to allow the visible and present masks to be cast out to the appropriate array shape
        vis2 = np.expand_dims(vis1,-1)
        vis3 = vis2
        pres2 = np.expand_dims(pres1,-1)
        pres3 = pres2
        if data_in.shape[2] > 1:
            for i in range(data_in.shape[2]-1):
                vis3 = np.concatenate((vis3,vis2),-1)
                pres3 = np.concatenate((pres3,pres2),-1)
    else:
        vis3=vis1
        pres3=pres1
            
    data_out = data_in[
            (np.rint((1-y_in2) * maph/2 - 0.5).astype(np.int64),
            np.rint((x_in2+1) * mapw/2 - 0.5).astype(np.int64))]
    blank = np.zeros_like(data_out)
    data_out = np.where(
        vis3 & pres3,
        data_out,
        blank).astype(data_in.dtype)
    print(" Outputting image...")
    map_out = Image.fromarray(data_out, mode=map_in.mode)
    if map_in.mode == "P":
        map_out.putpalette(map_in.getpalette())
    map_out.save(file_out)
    print("Finished; saved to "+file_out)
            


#Input    

def Inprompt(prompt,f):
    inp = input(prompt)
    while True:
        try:
            var = f(inp)
            return var
        except:
            inp = input(" Invalid input, try again: ")

if __name__ == "__main__":
    hem_in=0
    hem_out=0
    print('''
Projection Pasta
For Reprojection of maps between arbitrary aspects
Made 2022 by Amadea de Silva and Nikolai Hersfeldt

Projection Options and Codes (with profile):
  0: Equirectangular/Plate Caree (2:1 rectangle)
  1: Sinusoidal (2:1 sinusoid)
  2: Mollweide (2:1 ellipse)
  3: Hammer (2:1 ellipse)
  4: Aitoff (2:1 ellipse)
  5: Winkel Tripel (1.637:1 ovalish)
  6: Kavrayskiy VII (1.732:1 ovalish)
  7: Wagner VI (2:1 ovalish)
  8: Ortelius Oval (2:1 oval)
  9: Azimuthal Equidistant (1:1 circle)
 10: Orthographic (1:1 hemisphere)
 11: Stereographic (1:1 hemisphere)
 12: Mercator truncated to square (1:1 square)
 13: Gall Stereographic (1.301:1 rectangle)
 14: Miller Cylindrical (1.364:1 rectangle)

''')
    print("Input Image")
    while True:
        file_in = input(" Filename: ")
        if os.path.exists(file_in):
            break
        print("  No file found at "+str(file_in))
    proj_in = Inprompt(" Projection: ",int)
    if typl[proj_in] == Hemfull_typ or typl[proj_in] == HemfullIter_typ:
        hem_in = Inprompt("  0 for global, 1 for hemisphere, 2 for bihemisphere: ",int)
    elif typl[proj_in] == Hemonly_typ:
        hem_in = Inprompt("  1 for hemisphere, 2 for bihemisphere: ",int)
    lon_in = ma.radians(Inprompt(" Center longitude (-180 to 180): ",float))
    lat_in = ma.radians(Inprompt(" Center latitude (-90 to 90): ",float))
    rot_in = ma.radians(Inprompt(" Clockwise rotation from north (0 to 360): ",float))

    print("""
Output Image""")
    file_out = input(" Filename: ")
    proj_out = Inprompt(" Projection: ",int)
    if typl[proj_out] == Hemfull_typ or typl[proj_out] == HemfullIter_typ:
        hem_out = Inprompt("  0 for global, 1 for hemisphere, 2 for bihemisphere: ",int)
    elif typl[proj_out] == Hemonly_typ:
        hem_out = Inprompt("  1 for hemisphere, 2 for bihemisphere: ",int)
    lon_out = ma.radians(Inprompt(" Center longitude (-180 to 180): ",float))
    lat_out = ma.radians(Inprompt(" Center latitude (-90 to 90): ",float))
    rot_out = ma.radians(Inprompt(" Clockwise rotation from north (0 to 360): ",float))
    trunc = input("""
Crop map edges to single globe surface?
 Cropped maps may have missing pixels on edges when reprojected a second time
 Uncropped maps may have odd noise on edges in some projections (outside of the usually cropped area)
 y/n: """)
    if 'y' in trunc or 'Y' in trunc or '1' in trunc:
        trunc = 1
    else:
        trunc = 0
    print("""
Working...""")

    Main(file_in, file_out, proj_in, proj_out, lon_in, lat_in, rot_in, lon_out, lat_out, rot_out,hem_in=hem_in,hem_out=hem_out,trunc=trunc)



