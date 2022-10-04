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

#Libreria para ignorar warnings
import warnings
warnings.filterwarnings('ignore', category = FutureWarning)

#---------------------------------------------------------------------------------------
#-------------------------------------- Funciones --------------------------------------
#---------------------------------------------------------------------------------------

def corte_landsat(shp_path, image_path, output_path, bands_path):
	'''
	Argumentos:
	
		shp_path = path del archivo .shp del area de interes

		image_path = path de las imagenes

		output_path = path de carpeta donde quedaran las nuevas bandas

		bands_path = arreglo de los nombres de las bandas (imagen)

	Return:
	
		outputs_path = arreglo con los path de las nuevas imagenes cortadas

	Esta función corta las bandas dadas por imagenes landsat, centrandose 
	en un area de interes dado
	'''

	#Lista donde se guardaran las bandas cortadas
	outputs_path = []

	#Se recorre las bandas segun su path
	for band in bands_path:

		#Se realiza corte en base a .shp previo y se guarda como .tiff en carpeta
		options = gdal.WarpOptions(cutlineDSName = shp_path, cropToCutline = True)
		outBand = gdal.Warp(srcDSOrSrcDSTab = image_path + band,
							destNameOrDestDS = output_path + band[-6:-4] + '_corte' + band[-4:],
							options = options)
		outBand = None

		#Agregamos paths de bandas cortadas en arreglo
		outputs_path.append(band[-6:-4] + '_corte' + band[-4:])

	print('Corte en zona de interes realizado.')

	return outputs_path

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def corte_sentinel(shp_path, imagePath, outputPath, bands_path):
	'''
	Argumentos:
	
		shp_path = path del archivo .shp del area de interes

		imagePath = path de las imagenes

		outputPath = path de carpeta donde quedaran las nuevas bandas

		bands_path = arreglo de los nombres de las bandas (imagen)

	Return: 
	
		outputs_path = arreglo con los path de las nuevas imagenes cortadas

		outImage = Matriz de pixeles mascara

		outTransform = Metadatos de georeferenciacion

	Esta función corta las bandas dadas por imagenes landsat, centrandose 
	en un area de interes dado
	'''
	#Lista donde se guardaran las bandas cortadas
	outputs_path = []

	#Obtenemos el area de interes con sus respectivos datos de georeferenciacion
	aoiFile =fiona.open(shp_path)
	aoiGem = [aoiFile[0]['geometry']]

	#Corte de las bandas en la zona de interes
	for band in bands_path:
		# Abrimos las imagenes
		rasterPath = os.path.join(imagePath,band)
		rasterBand = rasterio.open(rasterPath)

		# Cortamos las imagenes en el area de interes
		outImage, outTransform = mask(rasterBand, aoiGem, crop=True)

		# Cargamos información meta
		outMeta = rasterBand.meta
		outMeta.update({"driver": 'JP2OpenJPEG', "height": outImage.shape[1],
						"width": outImage.shape[2], "transform": outTransform})

		#Se guarda la imagen
		outPath = outputPath + band[-7:-4] + '_corte' + band[-4:]
		outRaster = rasterio.open(outPath, "w", **outMeta)
		outRaster.write(outImage)
		outRaster.close()

		#Agregamos paths de bandas cortadas en arreglo
		outputs_path.append(band[-7:-4] + '_corte' + band[-4:])

	print('Corte en zona de interes realizado.')

	return outputs_path, outImage, outTransform

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def radiancia_reflectanciaTOA(bands_path, outputPath, dictionary):
	'''
	Argumentos:
	
		bands_path = arreglo con los path a los cuales se les hara una 
		 			 conversion de radiancia
		 
		outputPath = path de carpeta donde quedaran las nuevas bandas

		dictionary = diccionario con los metadatos de las bandas

	Return: 
	
		bands_ref = arreglo de matrices con los valores de los pixeles 
					despues de la conversion

	Esta funcion realiza una conversion de radiancia, para mejorar ciertos 
	factores que pueden causar errores en la visualizacion de la imagen 
	satelital para terminar realizando una conversion a reflectancia en el 
	techo de la atmosfera.
	'''
	#Lista donde se guardaran matrices de bandas
	bands_rad = []
	
	#Valores ESUN para landsat
	ESUN = [0, 2067, 1893, 1603, 972.6, 245.0, 79.72]
	
	#Obtenemos metadatos, angulo solar de elevacion y la distancia
	phi = float(dictionary['SUN_ELEVATION '])
	distance = float(dictionary['EARTH_SUN_DISTANCE '])

	#Se recorre las bandas segun su path
	for band in bands_path:
		
		#Se obtienen metadatos de la imagen para la radiancia
		rad_mul = float(dictionary['RADIANCE_MULT_BAND_'+ band[1] +' '])
		rad_sum = float(dictionary['RADIANCE_ADD_BAND_'+ band[1] +' '])

		#Se transforma la banda de bytes a float64 para el calculo
		banda = rasterio.open(outputPath + band).read(1).astype('float64')

		#Se forma matriz de mismas dimensiones para entregar
		new_band = banda.copy()

		#Se realiza el calculo de la radiancia y agregamos al arreglo
		new_band = new_band * rad_mul + rad_sum
		bands_rad.append(new_band)
	
	#PAra recorrer la lista considerando indica la primera banda	
	len_matrices = len(bands_rad)
	bands_ref = [bands_rad[0]]
	
	#Factor de TOA de valores que no cambian segun banda
	factor = (np.pi * (distance ** 2)) / cos((90 - phi) * np.pi / 180)
	
	for i in range(1, len_matrices):
		
		#Se crea matriz auxiliar
		new_band = bands_rad[i].copy()
		
		#Calculo de reflectancia TOA
		new_band = (factor * new_band) / ESUN[i]
		
		bands_ref.append(new_band)
		
	print('Radiancia y Reflectancia TOA realizadas.')

	return bands_ref

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def index_varios(band_1, band_2):
	'''
	Argumentos:
	
		band_1 = matriz de pixeles de la banda que lleva el factor 
				 positivo en la ecuacion
		band_2 = matriz de pixeles de la banda que lleva el factor 
				 positivo y negativo en la ecuacion

	Return: 
	
		INDEX = matriz con el calculo del indice para cada pixel

	Esta funcion realiza el calculo de los indices que requieran 2 bandas
	y la formula que sigan es la suma de estas como denominador y la resta
	en el numerador como los siguientes ejemplos: NDVI, GCI, SEI, NDWI, etc.
	'''
	#Matriz del indice a entregar
	INDEX = np.empty_like(band_1)

	rows, cols = INDEX.shape

	#Se recorre pixel a pixel
	for k in range(0, rows):
		for j in range(0, cols):

			#Valor del denominador en formula de indice
			denominador = band_1[k][j] + band_2[k][j]

			#Vemos si el denominador se indefine
			if denominador != 0:

				#Se realiza el calculo del indice ndwi
				INDEX[k][j] = (band_1[k][j] - band_2[k][j]) / denominador

			#Si el indice se indefine le damos el siguiente valor
			else:
				INDEX[k][j] = -0.8

	return INDEX

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def fai(band_1, band_2, band_3, factor):
	'''
	Argumentos:
	
		band_1 = matriz de pixeles de la banda RED
		 
		band_2 = matriz de pixeles de la banda NIR
		 
		band_3 = matriz de pixeles de la banda SWIR

	Return: 
	
		fai = matriz con el calculo del fai para cada pixel

	Esta funcion realiza el calculo de ndwi, y se puede ajustar al metodo
	de GAO o al de McFeeters
	'''
	#Copiamos una banda
	fai = band_2.copy()

	#Calculamos el indice FAI
	fai = fai - (band_1 + (band_3 - band_1) * factor)

	return fai

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def clorofila_index(band_1, band_2):
	'''
	Argumentos:
		band_1 = matriz de la banda que ira en el numerador
		
		band_2 = matriz de la banda que ira en el denominador
		
	Return:
		banda = matriz con el calculo del indice de clorofila para 
				cada pixel
		
	Esta funcion calcula los indices de clorofila verde y rojo
	'''
	#Dimensionar y obtenemos banda a entregar
	rows, cols = band_1.shape
	banda = np.zeros((rows, cols))
	
	#Para reemplazar en valores indefinidos
	den_mean = np.mean(band_1) / np.mean(band_2)
	
	#Recorremos matrices
	for i in range(0, rows):
		for j in range(0, cols):
		
			#Descartamos valores indefinidos y le damos el siguiente valor
			if band_2[i][j] == 0:
				banda[i][j] = den_mean - 1
				
			#En otro caso se emplea la formula del indice
			else:
				banda[i][j] = (band_1[i][j] / band_2[i][j]) - 1
	
	return banda

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def inter_agua_landsat(index, matrix_index, matrix_agua):
	'''
	Argumentos:
	
		index = string con el nombre del indice para establecer rangos

		matrix_index = matriz (banda) del indice al que se quiere ver 
		la interseccion

		matrix_agua = matriz binaria con pixeles de agua realzados

	Return:
	
		matriz = matriz de interseccion, pixeles que cumplen los rangos
				 tienen mayor visibilidad
		
		mask = mascara formada para el calculo de una metrica

	Esta funcion realiza la interseccion de matrices, en especifico, se
	asegura de formar dos matriz donde una sea para la visualización y 
	nos muestre un realce de los valores del indice requerido cuando cumple 
	los rangos establecidos para estar en zonas de agua, mientras que la 2da 
	es una matriz mascara para el calculo de alguna metrica.
	'''
	#Se crea una matriz dimensionada y una matriz mascara de valores nulos
	rows, cols = matrix_agua.shape
	matriz = matrix_index.copy()
	mask = np.zeros((rows, cols))
	
	#Se establecen rangos para los indices que

	#Indice FAI
	if index.lower() == 'fai':
		cota_min, cota_max = -0.15, -0.03
		suma = np.max(matrix_index) * 2

	#Indice SEI
	elif index.lower() == 'sei':
		cota_min, cota_max = 0.08, 0.2
		suma = np.max(matrix_index) * 2
	
	#Indice GCI
	elif index.lower() == 'gci':
		cota_min, cota_max = -0.5, -0.1
		suma = np.max(matrix_index) * 2
		
	#Indice RCI
	elif index.lower() == 'rci':
		cota_min, cota_max = -0.5, -0.1
		suma = np.max(matrix_index) * 2
	
	#No hay registro del indice pedido
	else:
		'indice no reconocido'
	
	#Se recorren las matrices de las bandas
	for i in range(0, rows):
		for j in range(0, cols):
			
			#Se verifican los pixeles dentro de los rangos dados para cada indice
			if matrix_agua[i][j] == 1  and (cota_min <= matrix_index[i][j] and matrix_index[i][j] <= cota_max):
			
				#Cambiamos los valores en la matriz para visualizacion y la mascara
				matriz[i][j] = matriz[i][j] + suma
				mask[i][j] = 1
			
	return matriz, mask
	
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def inter_agua_sentinel(index, matrix_index, matrix_agua):
	'''
	Argumentos:
	
		index = string con el nombre del indice para establecer rangos

		matrix_index = matriz (banda) del indice al que se quiere ver 
		la interseccion

		matrix_agua = matriz binaria con pixeles de agua realzados

	Return:
	
		matriz = matriz de interseccion, pixeles que cumplen los rangos
				 tienen mayor visibilidad
		
		mask = mascara formada para el calculo de una metrica

	Esta funcion realiza la interseccion de matrices, en especifico, se
	asegura de formar dos matriz donde una sea para la visualización y 
	nos muestre un realce de los valores del indice requerido cuando cumple 
	los rangos establecidos para estar en zonas de agua, mientras que la 2da 
	es una matriz mascara para el calculo de alguna metrica.
	'''
	#Se crea una matriz dimensionada y una matriz mascara de valores nulos
	rows, cols = matrix_agua.shape
	matriz = matrix_index.copy()
	mask = np.zeros((rows, cols))
	
	#Se establecen rangos para los indices que

	#Indice FAI
	if index.lower() == 'fai':
		cota_min, cota_max = -0.26, -0.1
		suma = np.max(matrix_index) * 1.25

	#Indice SEI
	elif index.lower() == 'sei':
		cota_min, cota_max = 0.08, 0.24
		suma = np.max(matrix_index) * 2
	
	#Indice GCI
	elif index.lower() == 'gci':
		cota_min, cota_max = -0.07, -0.008
		suma = np.max(matrix_index)
		
	#Indice RCI
	elif index.lower() == 'rci':
		cota_min, cota_max = -0.07, -0.005
		suma = np.max(matrix_index)
	
	#No hay registro del indice pedido
	else:
		'indice no reconocido'
	
	#Se recorren las matrices de las bandas
	for i in range(0, rows):
		for j in range(0, cols):
			
			#Se verifican los pixeles dentro de los rangos dados para cada indice
			if matrix_agua[i][j] == 1  and (cota_min <= matrix_index[i][j] and matrix_index[i][j] <= cota_max):
			
				#Cambiamos los valores en la matriz para visualizacion y la mascara
				matriz[i][j] = matriz[i][j] + suma
				mask[i][j] = 1
			
	return matriz, mask
	
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def mask_agua(matrix_gao, matrix_mf, matrix_xu, matrix_ndvi, sensor):
	'''
	Argumentos:
	
		matrix_gao = Matriz del indice NDWI Gao
		
		matrix_mf = Matriz del indice NDWI McFeeters
		
		matrix_xu = Matriz del indice NDWI Xu
		
		sensor = Satelite del cual se obtuvo la imagen
		
	Return:

		MASK = matriz binaria con pixeles en zonas de agua de
			   valor 1 y 0 para los demas
			   
	Esta funcion entrega una mascara para asegurar que se encuentra
	en zonas de agua, a traves de rangos dados para los indices
	NDVI y NDWI en 3 metodos, Gao, McFeeters y Xu
	'''
	#Se ve que sensor es
	if sensor == 'landsat':
		cota_ndwi, cota_ndvi = 0, 0.05
	
	else:
		cota_ndwi, cota_ndvi = -0.08, 0.03
	
	#Dimensiones de las matrices
	rows, cols = matrix_gao.shape
	
	#Mascara que se entrega de los pixeles de agua
	MASK = np.zeros((rows, cols))

	#Se recorre la matrix
	for i in range(0, rows):
		for j in range(0, cols):

			#Si los indices para el agua cumplen rangos
			if (matrix_gao[i][j] >= cota_ndwi and matrix_mf[i][j] >= cota_ndwi and matrix_xu[i][j] >= cota_ndwi) or matrix_ndvi[i][j] < cota_ndvi:

				#Mascara toma valores 1 si cumple rangos
				MASK[i][j] = 1

	return MASK
	
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def porcentaje_posicion_alga(alga_index, agua_index):
	'''
	Argumentos:
		
		alga_index = Matriz binaria (mask) de indice de algas
		
		agua_index = Matriz binaria (mask) de indice de agua
		
	Return:
		
		porcentaje = Porcentaje de pixeles de alga en zonas de agua
		
	Esta función entrega un porcentaje de la cantidad de pixeles de 
	algas en zonas de agua, es decir, saca un porcentaje del total de
	area ocupada por algas de un total de area ocupada por agua
	'''
	#Arreglo de de pixeles de alga
	list_alga = []
	
	#Se cuentan los pixeles de agua
	count_agua = np.count_nonzero(agua_index == 1)
	
	#Para que no se indefina (este es un caso extremo)
	if count_agua == 0:
		count_agua = 1
	
	#Dimensiones de matriz de indice de alga		
	rows, cols = alga_index.shape
	
	#Recorremos matriz
	for i in range(0, rows):
		for j in range(0, cols):
		
			#Condicion de presencia de alga
			if alga_index[i][j] == 1:
				list_alga.append((i, j))
	
	porcentaje = float(len(list_alga) / count_agua)
	
	return porcentaje, list_alga

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def distancia(z_1, z_2):
	'''
	Argumentos:
	
		z_1 = Vector 2D con posiciones en el plano cartesiano
		
		z_2 = Vector 2D con posiciones en el plano cartesiano
		
	Return:
	
		distance = Distancia Euclidiana
		
	Esta funcion calcula la distancia euclidiana de 2 vectores
	dados en el plano cartesiano
	'''
	#Calculo de la distancia Euclidiana
	distance = ((z_1[0] - z_2[0]) ** 2 + (z_1[1] - z_2[1]) ** 2) ** 0.5
	
	return distance

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def dispersion_alga(alga_list):
	'''
	Argumentos:
	
		alga_list = lista con vectores de posicion de algas
		
	Return:
	
		dispersion = Dispersion 
		
	Esta funcion calcula la dispersion que hay entre los
	pixeles de algas
	'''
	#Suma de distancias y distancia maxima entre pixeles de alga
	suma, d_max = 0, 0
	
	#Largo de la lista
	N = len(alga_list)
	
	#Se recorre la lista de posiciones de algas
	for vector_1 in alga_list:
		for vector_2 in alga_list:
			
			#Se calcula la distancia entre los vectores
			d = distancia(vector_1, vector_2)
			
			#Se suma al total de las distancias
			suma = suma + d
			
			#Se guarda la distancia maxima
			if d > d_max:
				d_max = d

	#Se calcula la dispersion
	dispersion = suma / (d_max * N * (N - 1))
	
	return dispersion

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def metrica(alga_index, agua_index, peso_1):
	'''
	Argumentos:
	
		alga_index = Matriz binaria (mask) de indice de algas
		
		agua_index = Matriz binaria (mask) de indice de agua
		
		peso_1 = Importancia que se le da al porcentaje y dispersion
				 0 <= peso_1 <= 1
	
	Return:
		
		valor = valor de la metrica propuesta
		
	Esta funcion calcula la metrica propuesta para la alerta, la cual
	considera el porcentaje de algas respecto al agua y su dispersion
	'''
	#Se calcula el porcentaje y una lista de posiciones de alga
	porcentaje, list_alga = porcentaje_posicion_alga(alga_index, agua_index)
	
	if porcentaje > 0.7:
		return 'Indice contaminado'
	
	#Se calcula la dispersion de la lista de posiciones de algas
	dispersion = dispersion_alga(list_alga)

	#Calculo de la metrica
	valor_metrica = (peso_1 * porcentaje) + ((1 - peso_1) * (1 - dispersion))

	return valor_metrica

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def percent_algae_in_landsat(alga_index, alga_matrix, percent):
	'''
	Argumentos:
	
		alga_index = String para reconocer indice de alga
	
		alga_matrix = Matriz de indice de alga
		
		percent = Porcentaje del rango que sera la nueva cota
				  minima
		
	Return:
	
		porcentaje = Porcentaje de algas sobre rango
		
	Esta función forma un rango modificado segun sea el indice
	y forma un porcentaje de las algas dentro del rango modificado
	sobre todos los pixeles de algas dentro del rango original
	para imagenes landsat
	'''
	#Indice FAI
	if alga_index.lower() == 'fai':
		cota_min, cota_max = -0.15, -0.03

	#Indice SEI
	elif alga_index.lower() == 'sei':
		cota_min, cota_max = 0.08, 0.2
		
	#Contadores para formar el porcentaje
	count_original, count_modificado = 0, 0
	
	#Se forma una nueva cota minima
	cota_mod = cota_max * percent + cota_min * (1 - percent)
	
	#Dimensiones de la matriz
	rows, cols = alga_matrix.shape
	
	#Se recorre la matriz
	for i in range(0, rows):
		for j in range(0, cols):
			
			#Se verifica que este dentro del rango original
			if alga_matrix[i][j] >= cota_min and alga_matrix[i][j] <= cota_max:
				
				#Contamos para el rango original
				count_original = count_original + 1
				
				#Se verifica para el rango modificado
				if alga_matrix[i][j] >= cota_mod:
					count_modificado = count_modificado + 1
	
	#Se calcula el porcentaje
	porcentaje = count_modificado / count_original
	
	return porcentaje

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

def percent_algae_in_sentinel(alga_index, alga_matrix, percent):
	'''
	Argumentos:
	
		alga_index = String para reconocer indice de alga
	
		alga_matrix = Matriz de indice de alga
		
		percent = Porcentaje del rango que sera la nueva cota
				  minima
		
	Return:
	
		porcentaje = Porcentaje de algas sobre rango
		
	Esta función forma un rango modificado segun sea el indice
	y forma un porcentaje de las algas dentro del rango modificado
	sobre todos los pixeles de algas dentro del rango original
	para imagenes sentinel
	'''
	#Indice FAI
	if alga_index.lower() == 'fai':
		cota_min, cota_max = -0.26, -0.1

	#Indice SEI
	elif alga_index.lower() == 'sei':
		cota_min, cota_max = 0.08, 0.24
		
	#Contadores para formar el porcentaje
	count_original, count_modificado = 0, 0
	
	#Se forma una nueva cota minima
	cota_mod = cota_max * percent + cota_min * (1 - percent)
	
	#Dimensiones de la matriz
	rows, cols = alga_matrix.shape
	
	#Se recorre la matriz
	for i in range(0, rows):
		for j in range(0, cols):
			
			#Se verifica que este dentro del rango original
			if alga_matrix[i][j] >= cota_min and alga_matrix[i][j] <= cota_max:
				
				#Contamos para el rango original
				count_original = count_original + 1
				
				#Se verifica para el rango modificado
				if alga_matrix[i][j] >= cota_mod:
					count_modificado = count_modificado + 1
	
	#Se calcula el porcentaje
	porcentaje = count_modificado / count_original
	
	return porcentaje

