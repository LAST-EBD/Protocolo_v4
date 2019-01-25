
import os, shutil, re, time, subprocess, pandas, rasterio, sys, urllib, fiona, sqlite3, math, pymongo
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal, gdalconst
from datetime import datetime, date

class Product(object):
    
    
    '''Esta clase genera los productos deinundacion, turbidez del agua y ndvi de las escenas normalizadas'''
    
        
    def __init__(self, ruta_nor):
        
        self.ruta_escena = ruta_nor
        self.escena = os.path.split(self.ruta_escena)[1]
        self.raiz = os.path.split(os.path.split(self.ruta_escena)[0])[0]
        print(self.raiz)
        self.nor = os.path.join(self.raiz, os.path.join('nor', self.escena))
        print(self.nor)
        self.ori = os.path.join(self.raiz, os.path.join('ori', self.escena))
        self.data = os.path.join(self.raiz, 'data')
        self.temp = os.path.join(self.raiz, 'temp')
        self.productos = os.path.join(self.raiz, 'pro')
        self.pro_esc = os.path.join(self.productos, self.escena)
        os.makedirs(self.pro_esc, exist_ok=True)
        
        if 'l8oli' in self.ruta_escena:
            self.sat = 'L8'
        elif 'l7etm' in self.escena:
            self.sat = 'L7'
        elif 'l5tm' in self.ruta_escena:
            self.sat =  'L5'
        elif 'l4tm' in self.ruta_escena:
            self.sat =  'L5'
        else:
            print('No identifico el satelite')
        
        print(self.sat)
        
        if self.sat == 'L8':

            for i in os.listdir(self.nor):
                if re.search('tif$', i):
                    
                    banda = i[-6:-4].lower()
                                        
                    if banda == 'b2':
                        self.blue = os.path.join(self.nor, i)
                    elif banda == 'b3':
                        self.green = os.path.join(self.nor, i)
                    elif banda == 'b4':
                        self.red = os.path.join(self.nor, i)
                    elif banda == 'b5':
                        self.nir = os.path.join(self.nor, i)
                    elif banda == 'b6':
                        self.swir1 = os.path.join(self.nor, i)
                    elif banda == 'b7':
                        self.swir2 = os.path.join(self.nor, i)
                    elif banda == 'k4':
                        self.fmask = os.path.join(self.nor, i)
                    
        else:

            for i in os.listdir(self.nor):
                if re.search('tif$', i):
                    
                    banda = i[-6:-4].lower()
                                        
                    if banda == 'b1':
                        self.blue = os.path.join(self.nor, i)
                    elif banda == 'b2':
                        self.green = os.path.join(self.nor, i)
                    elif banda == 'b3':
                        self.red = os.path.join(self.nor, i)
                    elif banda == 'b4':
                        self.nir = os.path.join(self.nor, i)
                    elif banda == 'b5':
                        self.swir1 = os.path.join(self.nor, i)
                    elif banda == 'b7':
                        self.swir2 = os.path.join(self.nor, i)
                    elif banda == 'k4':
                        self.fmask = os.path.join(self.nor, i)
        
        #Insertamos la cobertura de nubes en la BD
        connection = pymongo.MongoClient("mongodb://localhost")
        db=connection.teledeteccion
        landsat = db.landsat
        
        
        try:
        
            landsat.update_one({'_id':self.escena}, {'$set':{'Productos': []}},  upsert=True)
            
        except Exception as e:
            print("Unexpected error:", type(e), e)
            
        print('escena importada para productos correctamente')
        
        
        
    def ndvi(self):

        outfile = os.path.join(self.productos, self.escena + '_ndvi_.tif')
        print(outfile)
        
        with rasterio.open(self.nir) as nir:
            NIR = nir.read()
            
        with rasterio.open(self.red) as red:
            RED = red.read()

        num = NIR.astype(float)-RED.astype(float)
        den = NIR+RED
        ndvi = np.true_divide(num, den)
                
        profile = nir.meta
        profile.update(nodata=-9999)
        profile.update(dtype=rasterio.float32)

        with rasterio.open(outfile, 'w', **profile) as dst:
            dst.write(ndvi.astype(rasterio.float32))
            
        #Insertamos la cobertura de nubes en la BD
        connection = pymongo.MongoClient("mongodb://localhost")
        db=connection.teledeteccion
        landsat = db.landsat
        
        
        try:
        
            landsat.update_one({'_id':self.escena}, {'$set':{'Productos': ['NDVI']}},  upsert=True)
            
        except Exception as e:
            print("Unexpected error:", type(e), e)
            
        print('NDVI Generado')
        
        
        
    def flood(self):
        
        waterMask = os.path.join(self.data, 'water_mask_turb.tif')
        outfile = os.path.join(self.productos, self.escena + '_flood.tif')
        print(outfile)
        
        with rasterio.open(waterMask) as wmask:
            WMASK = wmask.read()
                        
        with rasterio.open(self.fmask) as fmask:
            FMASK = fmask.read()
            
        with rasterio.open(self.swir1) as swir1:
            SWIR1 = swir1.read()
            

        flood = np.where(((FMASK != 2) & (FMASK != 4)) & ((SWIR1 != 0) & (SWIR1 <= 1200)) & (WMASK > 0), 1, 0)
        
        
        profile = swir1.meta
        profile.update(nodata=0)
        profile.update(dtype=rasterio.ubyte)

        with rasterio.open(outfile, 'w', **profile) as dst:
            dst.write(flood.astype(rasterio.ubyte))
            
        #Insertamos la cobertura de nubes en la BD
        connection = pymongo.MongoClient("mongodb://localhost")
        db=connection.teledeteccion
        landsat = db.landsat
        
        
        try:
        
            landsat.update_one({'_id':self.escena}, {'$set':{'Productos': ['Flood']}},  upsert=True)
            
        except Exception as e:
            print("Unexpected error:", type(e), e)
            
        print('Fllod Mask Generada')
        
        
        
    def turbidity(self, flood):
        
        waterMask = os.path.join(self.data, 'water_mask_turb.tif')
        outfile = os.path.join(self.productos, self.escena + '_turbidity.tif')
        print(outfile)
        
        with rasterio.open(flood) as flood:
            FLOOD = flood.read()
        
        with rasterio.open(waterMask) as wmask:
            WMASK = wmask.read()
            
        with rasterio.open(self.blue) as blue:
            BLUE = blue.read()
            BLUE = np.where(BLUE == 0, 1, BLUE)
            BLUE = np.true_divide(BLUE, 10000)
                        
        with rasterio.open(self.green) as green:
            GREEN = green.read()
            GREEN = np.where(GREEN == 0, 1, GREEN)
            GREEN = np.true_divide(GREEN, 10000)
            GREEN_R = np.where((GREEN<0.1), 0.1, GREEN)
            GREEN_RECLASS = np.where((GREEN_R>=0.4), 0.4, GREEN_R)

        with rasterio.open(self.red) as red:
            RED = red.read()
            RED = np.where(RED == 0, 1, RED)
            RED = np.true_divide(RED, 10000)
            RED_RECLASS = np.where((RED>=0.2), 0.2, RED)
            
        with rasterio.open(self.nir) as nir:
            NIR = nir.read()
            NIR = np.where(NIR == 0, 1, NIR)
            NIR = np.true_divide(NIR, 10000)
            NIR_RECLASS = np.where((NIR>0.5), 0.5, NIR)
            
        with rasterio.open(self.swir1) as swir1:
            SWIR1 = swir1.read()
            SWIR1 = np.where(SWIR1 == 0, 1, SWIR1)
            SWIR1 = np.true_divide(SWIR1, 10000)
            SWIR_RECLASS = np.where((SWIR1>=0.09), 0.9, SWIR1)
        
        
        #Turbidez para la el rio
        rio = (-4.3 + (85.22 * GREEN_RECLASS) - (455.9 * np.power(GREEN_RECLASS,2)) \
            + (594.58 * np.power(GREEN_RECLASS,3)) + (32.3 * RED) - (15.36 * NIR_RECLASS)  \
            + (21 * np.power(NIR_RECLASS,2))) - 0.01        
        #RIO = np.power(math.e, rio)
        
        #Turbidez para la marisma        
        marisma = (4.1263574 + (18.8113118 * RED_RECLASS) - (32.2615219 * SWIR_RECLASS) \
        - 0.0114108989999999 * np.true_divide(BLUE, NIR)) - 0.01
        #MARISMA = np.power(math.e, marisma)
        
        
        TURBIDEZ = np.where(((FLOOD == 1) & (WMASK == 1)), marisma, 
                             np.where(((FLOOD == 1) & (WMASK == 2)), rio, 0))
        
        profile = swir1.meta
        profile.update(nodata=0)
        profile.update(dtype=rasterio.float32)
                             
        with rasterio.open(outfile, 'w', **profile) as dst:
            dst.write(TURBIDEZ.astype(rasterio.float32))
            
        #Insertamos la cobertura de nubes en la BD
        connection = pymongo.MongoClient("mongodb://localhost")
        db=connection.teledeteccion
        landsat = db.landsat
        
        
        try:
        
            landsat.update_one({'_id':self.escena}, {'$set':{'Productos': ['Turbidity']}},  upsert=True)
            
        except Exception as e:
            print("Unexpected error:", type(e), e)
            
        print('Turbidity Mask Generada')
