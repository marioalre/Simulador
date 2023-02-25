import urllib.request
import requests
from datetime import datetime

def descargar_archivo_mas_reciente():
    url = "https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field"
    # Abrir la URL y leer el contenido de la página web
    respuesta = urllib.request.urlopen(url)
    contenido = respuesta.read().decode()

    # Buscar todos los enlaces a archivos de texto
    enlaces_txt = []
    indice_inicio = 0
    while True:
        indice_inicio = contenido.find("href=\"", indice_inicio)
        if indice_inicio == -1:
            break
        indice_inicio += 6
        indice_fin = contenido.find(".txt", indice_inicio)
        if indice_fin == -1:
            break
        enlace = url + "/" + contenido[indice_inicio:indice_fin+4]
        enlaces_txt.append(enlace)

    # Analizar los nombres de archivo para obtener la fecha de cada archivo
    fechas = []
    for enlace in enlaces_txt:
        nombre_archivo = enlace.split("/")[-1]
        fecha_str = nombre_archivo.split("_")[1]
        fecha = datetime.strptime(fecha_str, "%Y.%m")
        fechas.append(fecha)

    # Encontrar la fecha más reciente y descargar el archivo correspondiente
    indice_archivo_mas_reciente = fechas.index(max(fechas))
    enlace_archivo_mas_reciente = enlaces_txt[indice_archivo_mas_reciente]
    archivo = enlace_archivo_mas_reciente.split("/")[-1]
    respuesta = requests.get(enlace_archivo_mas_reciente)
    with open(archivo, "wb") as f:
        f.write(respuesta.content)
    print(f"Se ha descargado el archivo {archivo} desde {url}.")

# Citando la fuente en el comentario de la función
# Esta función busca todos los archivos de texto 

descargar_archivo_mas_reciente()