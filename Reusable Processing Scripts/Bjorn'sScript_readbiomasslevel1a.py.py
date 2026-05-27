import rasterio
import numpy as np
import os
from lxml import etree
from PIL import Image


nom = input("Nom acquisition ? ")
nom1 = nom.lower()
nom2 = nom1.split("_01")[0]
# Chemin du fichier TIFF et du XML d'annotation
TIFF_PATH_AMP = nom+"/measurement/"+nom2+"_i_abs.tiff"
TIFF_PATH_PHA = nom+"/measurement/"+nom2+"_i_phase.tiff"
ANNOT_XML = nom+"/annotation/"+nom2+"_annot.xml"

# Extraire la liste des polarisations depuis le XML
with open(ANNOT_XML, 'rb') as f:
    tree = etree.parse(f)
    pol_nodes = tree.xpath('//polarisationList/polarisation')
    polarisations = [p.text.strip().lower() for p in pol_nodes]

# Ouvrir le fichier TIFF amplitude
with rasterio.open(TIFF_PATH_AMP) as srcamp:
    print(f"Taille de l'image amplitude : {srcamp.width} x {srcamp.height}")
    print(f"Nombre de polarisations amplitude : {srcamp.count}")
    nsamples = srcamp.width
    nlines = srcamp.height
    matricesamp = srcamp.read()
    if len(polarisations) != srcamp.count:
        raise ValueError("Le nombre de polarisations dans le XML ne correspond pas au nombre de bandes du TIFF.")

# Ouvrir le fichier TIFF phase
with rasterio.open(TIFF_PATH_PHA) as srcpha:
    print(f"Taille de l'image phase : {srcpha.width} x {srcpha.height}")
    print(f"Nombre de polarisations phase : {srcpha.count}")
    matricespha = srcpha.read()
    if len(polarisations) != srcpha.count:
        raise ValueError("Le nombre de polarisations dans le XML ne correspond pas au nombre de bandes du TIFF.")

repert = "output"
os.makedirs(repert, exist_ok=True)

# Ecriture fichier config
fichier = open(repert+"/config.txt","w")
fichier.write("Nrow"+"\n")
fichier.write(str(nlines)+"\n")
fichier.write("---------"+"\n")
fichier.write("Ncol"+"\n")
fichier.write(str(nsamples)+"\n")
fichier.write("---------"+"\n")
fichier.write("PolarCase"+"\n")
fichier.write("monostatic"+"\n")
fichier.write("---------"+"\n")
fichier.write("PolarType"+"\n")
fichier.write("full"+"\n")
fichier.close()

imagecplx = np.zeros((nlines,nsamples),dtype=np.complex64)
# Lecture amplitude et phase polarisations et ecriture format ENVI
for i in range(len(polarisations)):
    if (i == 0):
        nomfich = "s11"
    if (i == 1):
        nomfich = "s12"
    if (i == 2):
        nomfich = "s21"
    if (i == 3):
        nomfich = "s22"
    # Ecriture fichier HDR
    fichier = open(repert+"/"+nomfich+".bin.hdr","w")
    fichier.write("ENVI"+"\n")
    fichier.write("description = {"+"\n")
    fichier.write("PolSARpro File Imported to ENVI}"+"\n")
    fichier.write("samples = "+str(nsamples)+"\n")
    fichier.write("lines = "+str(nlines)+"\n")
    fichier.write("bands = 1"+"\n")
    fichier.write("header offset = 0"+"\n")
    fichier.write("file type = ENVI Standard"+"\n")
    fichier.write("data type = 6"+"\n")
    fichier.write("interleave = bsq"+"\n")
    fichier.write("sensor type = Unknown"+"\n")
    fichier.write("byte order = 0"+"\n")
    fichier.write("band names = {"+"\n")
    fichier.write(nomfich+".bin }"+"\n")
    fichier.close()
    matamp = matricesamp[i]
    matpha = matricespha[i]
    matamp_log = 20.0 * np.log10(np.clip(matamp, 1e-4, 5.0))
    mu = np.mean(matamp_log)
    sigma = np.std(matamp_log)
    print("-----------------------------------------------------------------------")
    print("Amplitude polarisation",polarisations[i])
    print("  lin : min =",np.min(matamp)," max =",np.max(matamp)," moy =",np.mean(matamp),"sdv =",np.std(matamp))
    print("  dB  : min =",np.min(matamp_log)," max =",np.max(matamp_log)," moy =",mu,"sdv =",sigma)
    print("Phase polarisation",polarisations[i])
    print("        min =",np.min(matpha)," max =",np.max(matpha)," moy =",np.mean(matpha),"sdv =",np.std(matpha)) 
    # Ecriture fichier BIN (float 2x32bits real,imag LE)
    fichier = open(repert+"/"+nomfich+".bin","wb")
    print("Image RAW :",repert+"/"+nomfich+".bin")
    imagecplx = matamp * np.cos(matpha) + matamp * np.sin(matpha) * 1.0j
    imagecplx.tofile(fichier)
    fichier.close()
    # Normalisation sur 0-255
    nbs = 6.0
    vmin = mu - nbs * sigma
    vmax = mu + nbs * sigma
    matamp_norm = np.clip((matamp_log - vmin) / (vmax - vmin), 0.0, 1.0)
    matamp_uint8 = (matamp_norm * 255.0).astype(np.uint8)
    # Ecriture fichier PNG
    img = Image.fromarray(matamp_uint8, mode='L')
    img.save(repert+"/"+nomfich+".png")
    print("Image PNG :",repert+"/"+nomfich+".png")
