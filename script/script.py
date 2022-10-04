#-----------------------------------------------------------------------------------------

#Libreria y funciones para manipulación y visualización de archivos .TIF
import rasterio
from rasterio import plot
from rasterio.mask import mask
from rasterio.plot import show

#-----------------------------------------------------------------------------------------

#Libreria para manipulación de archivos .shp
import fiona

#---------------------------------------------------------------------------------------

#Libreria para trabajar con datos geograficos
import gdal

#---------------------------------------------------------------------------------------

#Funcion para el reescalamiento de imagenes
from skimage.transform import rescale

#---------------------------------------------------------------------------------------

#Libreria para manipulacion y visualizacion de datos
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------------------

#Libreria de manejo de sistema para el manejo de paths
import os

#---------------------------------------------------------------------------------------

#Funcion trigonometrica
from math import sin, cos

#---------------------------------------------------------------------------------------

#Funciones creadas en el archivo Funciones.py
from Funciones import corte_landsat, corte_sentinel, radiancia_reflectanciaTOA
from Funciones import index_varios, clorofila_index, mask_agua, fai
from Funciones import inter_agua_landsat, inter_agua_sentinel
from Funciones import porcentaje_posicion_alga, distancia, dispersion_alga, metrica
from Funciones import percent_algae_in_landsat, percent_algae_in_sentinel

#---------------------------------------------------------------------------------------

#Libreria para ignorar warnings
import warnings
warnings.filterwarnings('ignore', category = FutureWarning)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#---------------------------------------- SCRIPT ----------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

#Path de la imagen satelital descargada
imagePath = 'Imagen_Estudio/Imagen_In/'

#Se ingresa el path del archivo .shp que contiene el area de interes
#este archivo debe estar georeferenciado 
shp_clip = 'Imagen_Estudio/Zona_Interes/zona_interes.shp'

#Path de salida de las nuevas bandas centradas en nuestro area de interes
cortePath = 'Imagen_Estudio/Cortes/'

#Aca vemos el contenido de la carpeta, para ver si es imagen sentinel o landsat
bandList = [band for band in os.listdir(imagePath) if 'SR_B' in band or '_B0' in band or '_B1' in band]
elemento = bandList[0]



#----------------------------------------------------------------------------------------
#--------------------------------------- Landsat ---------------------------------------
#----------------------------------------------------------------------------------------

if elemento[-3:] == 'TIF':
	#Se obtiene el archivo con los metadatos de la imagen satelital
	path_dict = [band for band in os.listdir(imagePath) if 'MTL.txt' in band][0]
	
	print('Realizando corte en zona de interes...')
	#Usamos la funcion de corte y creamos lista con los paths de los cortes
	cortesPath = corte_landsat(shp_clip, imagePath, cortePath, bandList)
	cortesPath.sort()
	cortesPath = cortesPath[1:]

	#Creamos un diccionario con los metadatos de la imagen
	# keys = metadatos ; values = valores de los metadatos
	dict_meta = {}

	with open(imagePath + path_dict, 'r') as _:
		for line in _:
		    line = line.strip()
		    if line != 'END':
		        key, value = line.split('=')
		        dict_meta[key] = value
		        
	#A las bandas, las cuales tenemos su path, les sacamos los valores reflectanciaTOA
	lista_bands_f = radiancia_reflectanciaTOA(cortesPath, cortePath, dict_meta)
	
	#Se abre una banda auxiliar para obtener metadatos al crear imagenes de los indices
	banda_aux = rasterio.open(cortePath + cortesPath[0])
	
	
	print('Calculando indices...')
	
	
#--------------------------------------- Indices ---------------------------------------

	#--------------------------------- NDVI ---------------------------------
	#Bandas espectrales necesarias
	nir_copy = lista_bands_f[4].copy()
	red_copy = lista_bands_f[3].copy()
	
	#Calculo del indice NDVI
	ndvi = index_varios(nir_copy, red_copy)

	#Se forma y se guarda la imagen con el calculo del indice NDVI
	ndviImage = rasterio.open('Imagen_Estudio/Index_Out/ndvi_landsat.tiff', 'w', driver = 'Gtiff',
				               width = banda_aux.width, height = banda_aux.height, count = 1,
				               crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	ndviImage.write(ndvi, 1)
	ndviImage.close()

	#------------------------------- NDWI Gao -------------------------------
	#Bandas espectrales necesarias
	nir_copy_g = lista_bands_f[4].copy()
	swir_copy_g = lista_bands_f[5].copy()
	
	#Calculo del indice NDWI Gao
	ndwi_g = index_varios(nir_copy_g, swir_copy_g)

	#Se forma, guarda y cierra la imagen con el calculo del NDWI Gao
	ndwiImage_g = rasterio.open('Imagen_Estudio/Index_Out/ndwi_g_landsat.tiff', 'w', driver = 'Gtiff',
				                width = banda_aux.width, height = banda_aux.height, count = 1,
				                crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	ndwiImage_g.write(ndwi_g, 1)
	ndwiImage_g.close()

	#---------------------------- NDWI McFeeters ----------------------------
	#Bandas espectrales necesarias
	green_copy_m = lista_bands_f[2].copy()
	nir_copy_m = lista_bands_f[4].copy()
	
	#Calculo del indice NDWI McFeeters
	ndwi_m = index_varios(green_copy_m, nir_copy_m)

	#Se forma, guarda y cierra la imagen con el calculo del NDWI McFeeters
	ndwiImage_m = rasterio.open('Imagen_Estudio/Index_Out/ndwi_m_landsat.tiff', 'w', driver = 'Gtiff',
				                width = banda_aux.width, height = banda_aux.height, count = 1,
				                crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	ndwiImage_m.write(ndwi_m, 1)
	ndwiImage_m.close()

	#------------------------------- NDWI Xu -------------------------------
	#Bandas espectrales necesarias
	green_copy_x = lista_bands_f[2].copy()
	swir_copy_x = lista_bands_f[5].copy()
	
	#Calculo del indice NDWI Xu
	ndwi_x = index_varios(green_copy_x, swir_copy_x)

	#Se forma, guarda y cierra la imagen con el calculo del NDWI Xu
	ndwiImage_x = rasterio.open('Imagen_Estudio/Index_Out/ndwi_x_landsat.tiff', 'w', driver = 'Gtiff',
				                width = banda_aux.width, height = banda_aux.height, count = 1,
				                crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	ndwiImage_x.write(ndwi_x, 1)
	ndwiImage_x.close()	

	#--------------------------------- FAI ---------------------------------
	# factor = (lambdaNIR - lambdaRED) / (lambdaSWIR - lambdaRED)
	factor = (865 - 655)/(1610 - 655)

	#Bandas espectrales necesarias
	RED_F = lista_bands_f[3].copy()
	NIR_F = lista_bands_f[4].copy()
	SWIR_F = lista_bands_f[5].copy()

	#Calculo del indice FAI, normalizamos y se deja el indice entre -1 y 1
	fai = fai(RED_F, NIR_F, SWIR_F, factor)
	fai = (fai / np.linalg.norm(fai)) * 100

	#Aca creamos la imagen con el indice fai, se guardara en la carpeta output
	faiImage = rasterio.open('Imagen_Estudio/Index_Out/fai_landsat.tiff', 'w', driver = 'Gtiff',
			                 width = banda_aux.width, height = banda_aux.height, count = 1,
			                 crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	faiImage.write(fai, 1)
	faiImage.close()
	
	#Se calcula la metrica 2 para el indice FAI
	metric2_fai = percent_algae_in_landsat('fai', fai, 0.7)	
	
	#--------------------------------- SEI ---------------------------------
	#Bandas espectrales necesarias
	NIR_S = lista_bands_f[4].copy()
	SWIR_S = lista_bands_f[5].copy()
	
	#Calculo del indice SEI
	sei = index_varios(NIR_S, SWIR_S)
	
	seiImage = rasterio.open('Imagen_Estudio/Index_Out/sei_landsat.tiff', 'w', driver = 'Gtiff',
			                 width = banda_aux.width, height = banda_aux.height, count = 1,
			                 crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	seiImage.write(sei, 1)
	seiImage.close()
	
	#Se calcula la metrica 2 para el indice SEI
	metric2_sei = percent_algae_in_landsat('sei', sei, 0.7)	
	
	#--------------------------------- GCI ---------------------------------
	#Bandas espectrales necesarias
	NIR_G = lista_bands_f[4].copy()
	GREEN_G = lista_bands_f[2].copy()
	
	#Calculo del indice SEI
	gci = clorofila_index(NIR_G, GREEN_G)
	
	gciImage = rasterio.open('Imagen_Estudio/Index_Out/gci_landsat.tiff', 'w', driver = 'Gtiff',
			                 width = banda_aux.width, height = banda_aux.height, count = 1,
			                 crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	gciImage.write(gci, 1)
	gciImage.close()
	
	#--------------------------------- RCI ---------------------------------
	#Bandas espectrales necesarias
	NIR_R = lista_bands_f[4].copy()
	RED_R = lista_bands_f[3].copy()
	
	#Calculo del indice SEI
	rci = clorofila_index(NIR_R, RED_R)
	
	rciImage = rasterio.open('Imagen_Estudio/Index_Out/rci_landsat.tiff', 'w', driver = 'Gtiff',
			                 width = banda_aux.width, height = banda_aux.height, count = 1,
			                 crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	rciImage.write(rci, 1)
	rciImage.close()	
	
	
	print('Indices calculados.\n')
	print('Realizando intersecciones...')
	
	
#----------------------------------- Intersecciones -----------------------------------

	#Se forma una mascara de zonas de agua
	agua = mask_agua(ndwi_g, ndwi_m, ndwi_x, ndvi, 'landsat')
	
	#---------------------------- FAI en el agua ----------------------------
	
	#Se realiza una interseccion con zonas de agua y el indice fai
	fai_inter, mask_fai = inter_agua_landsat('fai', fai, agua)

	#Aca creamos la imagen con la interseccion del indice fai
	faiImage = rasterio.open('Imagen_Estudio/Index_Out/fai_inter_landsat.tiff', 'w', driver = 'Gtiff',
				               width = banda_aux.width, height = banda_aux.height, count = 1,
				               crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	faiImage.write(fai_inter, 1)
	faiImage.close()
	
	#---------------------------- SEI en el agua ----------------------------
	
	#Se realiza una interseccion con zonas de agua y el indice sei	
	sei_inter, mask_sei = inter_agua_landsat('sei', sei, agua)

	#Aca creamos la imagen con la interseccion del indice sei
	seiImage = rasterio.open('Imagen_Estudio/Index_Out/sei_inter_landsat.tiff', 'w', driver = 'Gtiff',
				               width = banda_aux.width, height = banda_aux.height, count = 1,
				               crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	seiImage.write(sei_inter, 1)
	seiImage.close()	

	#---------------------------- GCI en el agua ----------------------------
	#Se realiza una interseccion con zonas de agua y el indice gci	
	gci_inter, mask_gci = inter_agua_landsat('gci', gci, agua)

	#Aca creamos la imagen con la interseccion del indice sei
	gciImage = rasterio.open('Imagen_Estudio/Index_Out/gci_inter_landsat.tiff', 'w', driver = 'Gtiff',
				               width = banda_aux.width, height = banda_aux.height, count = 1,
				               crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	gciImage.write(gci_inter, 1)
	gciImage.close()

	#---------------------------- RCI en el agua ----------------------------
	#Se realiza una interseccion con zonas de agua y el indice gci	
	rci_inter, mask_rci = inter_agua_landsat('rci', rci, agua)

	#Aca creamos la imagen con la interseccion del indice sei
	rciImage = rasterio.open('Imagen_Estudio/Index_Out/rci_inter_landsat.tiff', 'w', driver = 'Gtiff',
				               width = banda_aux.width, height = banda_aux.height, count = 1,
				               crs = banda_aux.crs, transform = banda_aux.transform, dtype = 'float64')
	rciImage.write(rci_inter, 1)
	rciImage.close()


	print('Intersecciones realizadas.')



#----------------------------------------------------------------------------------------
#--------------------------------------- Sentinel ---------------------------------------
#----------------------------------------------------------------------------------------

else:
	print('Realizando corte en zona de interes...')
	
	#Usamos la funcion de corte y creamos lista con los paths de los cortes
	#junto a datos necesarios para la formacion de imagenes
	cortesPath, outImage, outTransform = corte_sentinel(shp_clip, imagePath, cortePath, bandList)
	cortesPath.sort()

	#Se abren las bandas y se agregan a una lista
	list_bands = [ (rasterio.open(cortePath + band).read(1).astype('float64')) / 10000 for band in cortesPath]

	#Reescalamos las bandas para tener resolucion de 10 metros
	banda3 = list_bands[2]
	banda11 = list_bands[10]
	
	#Factor de reescalamiento
	escala = banda3.shape[0] / banda11.shape[0]
	
	#Reescalamiento
	banda11 = rescale(banda11, escala, anti_aliasing = False)
	
	#Formamos lista final de bandas para trabajar
	lista_bands_f = list_bands[0:10] + [banda11, list_bands[11]]

	
	print('Calculando indices...')

	
#--------------------------------------- Indices ---------------------------------------

	#--------------------------------- NDVI ---------------------------------
	#Bandas espectrales necesarias
	nir_copy = lista_bands_f[7].copy()
	red_copy = lista_bands_f[3].copy()
	
	#Calculo del indice NDVI
	ndvi = index_varios(nir_copy, red_copy)

	#Se forma, guarda y cierra la imagen con el calculo del NDVI
	ndviImage = rasterio.open('Imagen_Estudio/Index_Out/ndvi_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	ndviImage.write(ndvi, 1)
	ndviImage.close()	

	#------------------------------- NDWI Gao -------------------------------
	#Bandas espectrales necesarias
	nir_copy = lista_bands_f[7].copy()
	swir_copy = lista_bands_f[10].copy()
	
	#Calculo del indice NDWI Gao
	ndwi_g = index_varios(nir_copy, swir_copy)
	
	#Guardaremos la imagen de NDWI Gao
	ndwiImage_g = rasterio.open('Imagen_Estudio/Index_Out/ndwi_g_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	ndwiImage_g.write(ndwi_g, 1)
	ndwiImage_g.close()
	
	#---------------------------- NDWI McFeeters ----------------------------
	#Bandas espectrales necesarias
	green_copy = lista_bands_f[2].copy()
	nir_copy = lista_bands_f[7].copy()
	
	#Calculo del indice NDWI McFeeters
	ndwi_m = index_varios(green_copy, nir_copy)
	
	#Guardaremos la imagen de NDWI McFeeters
	ndwiImage_m = rasterio.open('Imagen_Estudio/Index_Out/ndwi_m_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	ndwiImage_m.write(ndwi_m, 1)
	ndwiImage_m.close()
	
	#------------------------------- NDWI Xu -------------------------------
	#Bandas espectrales necesarias
	green_copy = lista_bands_f[2].copy()
	swir_copy = lista_bands_f[10].copy()
	
	#Calculo del indice NDWI Xu
	ndwi_x = index_varios(green_copy, swir_copy)
	
	#Guardaremos la imagen de NDWI MXu
	ndwiImage_x = rasterio.open('Imagen_Estudio/Index_Out/ndwi_x_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	ndwiImage_x.write(ndwi_x, 1)
	ndwiImage_x.close()
	
	#--------------------------------- FAI ---------------------------------
	#Factor = (lambdaNIR - lambdaRED) / (lambdaSWIR - lambdaRED)
	factor = (833 - 665) / (1610.4 - 665)
	
	#Bandas espectrales necesarias
	RED = lista_bands_f[3].copy()
	NIR = lista_bands_f[7].copy()
	SWIR = lista_bands_f[10].copy()
	
	#Calculo del indice FAI
	fai = fai(RED, NIR, SWIR, factor)
	fai = (fai / np.linalg.norm(fai)) * 100
	
	#Guardaremos el indice fai
	faiImage = rasterio.open('Imagen_Estudio/Index_Out/fai_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	faiImage.write(fai, 1)
	faiImage.close()
	
	#Se calcula la metrica 2 para el indice FAI
	metric2_fai = percent_algae_in_sentinel('fai', fai, 0.7)

	#--------------------------------- SEI ---------------------------------
	#Bandas espectrales necesarias
	NIR = lista_bands_f[7].copy()
	SWIR = lista_bands_f[10].copy()
	
	#Calculo del indice SEI
	sei = index_varios(NIR, SWIR)
	
	seiImage = rasterio.open('Imagen_Estudio/Index_Out/sei_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	seiImage.write(sei, 1)
	seiImage.close()
	
	#Se calcula la metrica 2 para el indice SEI
	metric2_sei = percent_algae_in_sentinel('sei', sei, 0.7)	
	
	#--------------------------------- GCI ---------------------------------
	#Bandas espectrales necesarias
	NIR_G = lista_bands_f[7].copy()
	GREEN_G = lista_bands_f[2].copy()
	
	#Calculo del indice SEI
	gci = clorofila_index(NIR_G, GREEN_G)
	gci = (gci / np.linalg.norm(gci)) * 10
	
	gciImage = rasterio.open('Imagen_Estudio/Index_Out/gci_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	gciImage.write(gci, 1)
	gciImage.close()
	
	#--------------------------------- RCI ---------------------------------
	#Bandas espectrales necesarias
	NIR_R = lista_bands_f[7].copy()
	RED_R = lista_bands_f[3].copy()
	
	#Calculo del indice SEI
	rci = clorofila_index(NIR_R, RED_R)
	rci = (rci / np.linalg.norm(rci)) * 10
	
	rciImage = rasterio.open('Imagen_Estudio/Index_Out/rci_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	rciImage.write(rci, 1)
	rciImage.close()	

	
	print('Indices calculados.\n')
	print('Realizando intersecciones...')

	
#----------------------------------- Intersecciones -----------------------------------


	#Se forma una mascara de zonas de agua
	agua = mask_agua(ndwi_g, ndwi_m, ndwi_x, ndvi, 'sentinel')
	
	#---------------------------- FAI en el agua ----------------------------
	
	#Se realiza una interseccion con zonas de agua y el indice fai
	fai_inter, mask_fai = inter_agua_sentinel('fai', fai, agua)

	#Aca creamos la imagen con la interseccion del indice fai
	faiImage = rasterio.open('Imagen_Estudio/Index_Out/fai_inter_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	faiImage.write(fai_inter, 1)
	faiImage.close()
	
	#---------------------------- SEI en el agua ----------------------------
	
	#Se realiza una interseccion con zonas de agua y el indice sei	
	sei_inter, mask_sei = inter_agua_sentinel('sei', sei, agua)

	#Aca creamos la imagen con la interseccion del indice sei
	seiImage = rasterio.open('Imagen_Estudio/Index_Out/sei_inter_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	seiImage.write(sei_inter, 1)
	seiImage.close()	
	
	#---------------------------- GCI en el agua ----------------------------
	#Se realiza una interseccion con zonas de agua y el indice gci	
	gci_inter, mask_gci = inter_agua_sentinel('gci', gci, agua)

	#Aca creamos la imagen con la interseccion del indice sei
	gciImage = rasterio.open('Imagen_Estudio/Index_Out/gci_inter_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	gciImage.write(gci_inter, 1)
	gciImage.close()

	#---------------------------- RCI en el agua ----------------------------
	#Se realiza una interseccion con zonas de agua y el indice gci	
	rci_inter, mask_rci = inter_agua_sentinel('rci', rci, agua)

	#Aca creamos la imagen con la interseccion del indice sei
	rciImage = rasterio.open('Imagen_Estudio/Index_Out/rci_inter_senti.tiff', 'w', driver = 'Gtiff',
								width = outImage.shape[2], height = outImage.shape[1], count = 1,
								transform = outTransform, dtype = 'float64')
	rciImage.write(rci_inter, 1)
	rciImage.close()	


	print('Intersecciones realizadas.')



#----------------------------------------------------------------------------------------
#--------------------------------------- Metricas ---------------------------------------
#----------------------------------------------------------------------------------------

print('Calculando Metricas:\n')

#Calculo de metrica 1
p1 = 0.7
metric_fai = metrica(mask_fai, agua, p1)
metric_sei = metrica(mask_sei, agua, p1)

#Se promedia sei y fai en las metricas
metric_one = (metric_fai + metric_sei) / 2

print('\nEl valor del indice de dispersion es:', round(metric_one, 4))

#Si supera el umbral se tira alerta por parte de la metrica 1
if metric_one >= 0.3:
	metric_one = 1
	
#Se promedia sei y fai en las metricas
metric_two = (metric2_fai + metric2_sei) / 2

print('\nEl porcentaje de pixeles en la parte superior del rango es:', round(metric_two, 4))

#Si supera el umbral se tira alerta por parte de la metrica 1
if metric_two >= 0.5:
	metric_two = 1
	
#Alerta
if metric_one == 1 and metric_two == 1:
	print('\n\nEmitir alerta, ambas metricas pasaron el umbral!')
else:
	print('\n\nNo se emite alerta.')


#----------------------------------------------------------------------------------------
#--------------------------------------- Graficos ---------------------------------------
#----------------------------------------------------------------------------------------

#------------------------------ Grafico FAI ------------------------------
fig_f, ax_f = plt.subplots(1, 1, figsize = (8, 8))
plt.title('FAI en el agua')
plot.show(fai_inter, ax = ax_f)
plt.yticks([])
plt.xticks([])
fig_f.tight_layout()

#------------------------------ Grafico SEI ------------------------------
fig_s, ax_s = plt.subplots(1, 1, figsize = (8, 8))
plt.title('SEI en el agua')
plot.show(sei_inter, ax = ax_s)
plt.yticks([])
plt.xticks([])
fig_s.tight_layout()

#------------------------------ Grafico GCI ------------------------------
fig_g, ax_g = plt.subplots(1, 1, figsize = (8, 8))
plt.title('GCI en el agua')
plot.show(gci_inter, ax = ax_g)
plt.yticks([])
plt.xticks([])
fig_g.tight_layout()

#------------------------------ Grafico GCI ------------------------------
fig_r, ax_r = plt.subplots(1, 1, figsize = (8, 8))
plt.title('RCI en el agua')
plot.show(rci_inter, ax = ax_r)
plt.yticks([])
plt.xticks([])
fig_r.tight_layout()

plt.show()





















