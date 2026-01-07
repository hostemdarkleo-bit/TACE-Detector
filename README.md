Instructions in English

This repository contains the essential files for replicating the results of the manuscript titled:

"Semi-Automatic model based on topologic isophotes in Sobolev spaces to identify objects in astrophysical images: Base of an AI algorithm"

The functions in both Matlab and Python are in this repository, as well as the link to obtain the images, where we only need to specify or clarify the filters and bands used. All functions must be loaded to obtain results, since the TACE (Topological Adaptation in Contour Evolution) model was programmed modularly, and the functions are correlated, as can be seen in the TACE evolution expression, Eq. 13 of the manuscript. Matlab version R2024b was used (we also recommend downloading all available toolboxes from your Matlab installation, as some functions may not install automatically and could cause errors in script execution). There is no specific order for executing the functions; they only need to be loaded into the main function or the corresponding folder. The script “DetTACE.m” is where everything is executed, and the results are displayed.

DetTACE.m
IntFilImg.m
LogFiltermIm.m
create_mask.m
estimate_radius.m
SMosaics.m
TACE_Max_NObj.m
TACE_SegmentObj.m

Link to the images (raw or unprocessed data) of “El Gordo (The Fat)” One)”:
https://hla.stsci.edu/hlaview.html#Inventory|filterText%3D%24filterTypes%3D|query_string=El%20gordo&posfilename=&poslocalname=&posfilecount=&listdelimiter=whitespace&listformat=degrees&RA=15.718750&Dec=-49.249440&Radius=0.200000&inst-control=all&inst=ACS&inst=ACSGrism&inst=WFC3&inst=WFPC2&inst=NICMOS&inst=NICGRISM&inst=COS&inst=WFPC2-PC&inst=STIS&inst=FOS&inst=GHRS&imagetype=best&prop_id=&spectral_elt=&proprietary=both&preview=1&output_size=256&cutout_size=12.8|ra=&dec=&sr=&level=&image=&inst=ACS%2CACSGrism%2CWFC3%2CWFPC2%2CNICMOS%2CNICGRISM%2CCOS%2CWFPC2-PC%2CSTIS%2CFOS%2CGHRS&ds=

The link for the M82 image is corresponding to the file h_m82_i_s05_drz_sci.fits obtained from the work of Mutchler et al. [21]:
https://iopscience.iop.org/article/10.1086/511160

Graphical User Interface (GUI) for Object Detection in Astrophysical Images from Different Instruments (Hubble or James Webb Space Telescope).

This repository contains the Python algorithms that comprise a Python-based graphical user interface (GUI) designed for the detection and segmentation of objects such as galaxies or stars (as bright objects) in astronomical images. It highlights the large size of these unprocessed images (commonly in shades of gray) compared to everyday images. Our proposal allows loading FITS files, visualizing their structure, dividing them into user-defined blocks (mosaics), and applying adaptive segmentation methods to isolate the bright regions of interest (ROIs) that constitute our TACE model.


Our numerical approach in the Python GUI is designed for high-resolution images, approximately 13000 x 13000 pixels in size. These images are divided into tiles (equal-sized sections of the image to be processed) for efficient and interactive processing. By allowing the user to segment the image into manageable blocks, the tool significantly reduces processing time and memory usage.

Features: The image is first read or loaded in FITS format using astropy.io.fits. Automatic dimension detection: It reads the image dimensions upon loading and reports the size, shape, and total number of pixels. User-defined block division: The user can choose any block size into which the image to be processed will be divided equally between the width and height of the image. Central cropping extraction: The size of each image in the tile is the same, but we note that for a size of 400x400, we have a maximum numerical prediction of the approach in Python only. Maximum Detection Algorithm (Lumen_max): Identifies local maxima within the clipping, representing bright astronomical sources. Segmentation Modes: Adaptive Radius: Automatically adjusts the object radius.

Static Radius: Uses a mask radius that is calculated variationally from the isophotes set, Eq. 13 of the manuscript. Interactive Visualization: The graphical user interface displays the positions of the maxima and their corresponding binary masks. Included Backend Algorithms: Custom segmentation routines (lumen_segment_nobj, lumen_segment_nobj1) included in the repository.

Requirements: To run the algorithms in Python (which are in this repository with the *.py extension), you need:
Python 3.8+
Tkinter (bundled with Python)
NumPy
Matplotlib
Astropy
As


Instrucciones en Español

Este repositorio considera aquellos archivos indispensables para replicar los resultados del manuscrito titulado:
"Semi-Automatic model based on topologic isophotes in Sobolev spaces to identify objects in astrophysical images: Base of an AI algorithm”.

Las funciones tanto en Matlab como en Python se encuentran en este repositorio y de igual manera la liga para la obtención de las imágenes donde solo especificamos o claramos los filtros y bandas empleadas. Todas las funciones deben ser cargadas para la obtención de resultados ya que el modelo TACE (Topological Adaptation in Contour Evolution) se programó de manera modular en las funciones y las funciones están correlacionadas, como se puede entender de la expresión de evolución de TACE, Eq. 13 del manuscrito. Se empleó de Matlab la versión R2024b (adeicionalmente recomendamos que se descarguen todos los toolbox disponibles en la instalación de matlab ya que en ocasiones hay funciones que no instalan y puede ocasionar error en la ejecución de scripts de manera general) donde no existe un orden para ejecutar las funciones solo que estén cargadas en el main o en la carpeta, el script “DetTACE.m” es donde se ejecuta todo y se observan los resultados.
DetTACE.m
IntFilImg.m
LogFiltermIm.m
create_mask.m
estimate_radius.m
SMosaics.m
TACE_Max_NObj.m
TACE_SegmentObj.m

Link donde se encuentran las imágenes (datos crudos o sin tratamiento) de “El Gordo (The Fat One)”:
https://hla.stsci.edu/hlaview.html#Inventory%7CfilterText%3D%24filterTypes%3D%7Cquery_string=El%20gordo&posfilename=&poslocalname=&posfilecount=&listdelimiter=whitespace&listformat=degrees&RA=15.718750&Dec=-49.249440&Radius=0.200000&inst-control=all&inst=ACS&inst=ACSGrism&inst=WFC3&inst=WFPC2&inst=NICMOS&inst=NICGRISM&inst=COS&inst=WFPC2-PC&inst=STIS&inst=FOS&inst=GHRS&imagetype=best&prop_id=&spectral_elt=&proprietary=both&preview=1&output_size=256&cutout_size=12.8%7Cra=&dec=&sr=&level=&image=&inst=ACS%2CACSGrism%2CWFC3%2CWFPC2%2CNICMOS%2CNICGRISM%2CCOS%2CWFPC2-PC%2CSTIS%2CFOS%2CGHRS&ds=

El link para la imagen de M82 es correspondiente al archivo h_m82_i_s05_drz_sci.fits obtenida del trabajo de Mutchler et al. [21]:
https://iopscience.iop.org/article/10.1086/511160

Interfaz gráfica de usuario (GUI) para la detección de objetos en imágenes astrofísicas por diferentes instrumentos (Hubble o James Web).

Este repositorio contiene los algoritmos en Python que conforman una interfaz gráfica de usuario (GUi) basada en Python diseñada para la detección y segmentación de objetos como galaxias o estrellas (como objetos brillantes) en imágenes astronómicas; resaltando el gran tamaño de estas imágenes sin tratamiento (comúnmente en tonos de gris) respecto de imágenes cotidianas. Nuestra propuesta permite cargar archivos FITS, visualizar su estructura, dividirlos en bloques definidos por el usuario (mosaicos) y aplicar métodos de segmentación adaptativos para aislar las regiones luminosas de interés o ROÍ (Regiones De Interés—por su sigla en inglés) que constituyen nuestro modelo TACE.

Nuestra propuesta numérica en la GUI de Python está diseñada para imágenes de alta resolución, con un tamaño aproximado de 13000X13000 pixeles; las cuales se dividieron en mosaicos (divisiones de igual tamaño de la imagen a procesar) para procesarse de forma eficiente e interactiva. Al permitir al usuario segmentar la imagen en bloques manejables, la herramienta reduce significativamente el tiempo de procesamiento y el uso de memoria.

Características: Primero se lee o se carga la imagen en formato FITS mediante el uso de astropy.io.fits. Detección automática de dimensiones: Lee las dimensiones de la imagen al cargarla e informa sobre el tamaño, la forma y el número total de píxeles. División de bloques definida por el usuario: El usuario puede elegir cualquier tamaño de bloque en que se dividirá la imagen a procesar de forma equitativamente entre el ancho y la altura de la imagen. Extracción del Recorte Central: El tamaño de cada su imagen en el mosaico es del mismo tamaño pero resaltamos que para un tamaño de 400X400 tenemos una predicción numérica máxima de la propuesta solo en Python. Algoritmo de Detección Máxima (Lumen_max): Identifica los máximos locales dentro del recorte, representando fuentes astronómicas brillantes. Modos de Segmentación: Radio Adaptativo: Ajusta automáticamente el radio del objeto.

Radio Estático: Utiliza un radio de máscara que se calcula variacionalmente del conjunto de isofotas, Eq. 13 del manuscrito. Visualización Interactiva: La interfaz gráfica de usuario muestra las posiciones de los máximos y sus correspondientes máscaras binarias. Algoritmos de Backend Incluidos: Rutinas de segmentación personalizadas (lumen_segment_nobj, lumen_segment_nobj1) incluidas en el repositorio.

Requisitos: Para la ejecución de los algoritmos en Python (mismos que están en este repositorio con extensión *.py) necesitamos:
Python 3.8+
Tkinter (bundled with Python)
NumPy
Matplotlib
Astropy
Astropy Cutout2D utilities
tkinter.ttk widgets

