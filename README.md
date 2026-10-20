# AACclustMac
# indstalar estos paquetes:
install.packages(c("Rcpp", "cluster", "dbscan", "amap", "devtools", "roxygen2"))

### en el terminal Mac: 
brew install libomp
#### subir el paquete de mi pagina web :  AACclustMac

# Ubicarte en la carpeta PARENT de AACclustMac:
setwd("/Users/tarikfaouzi/Documents/Tarik_clustering")

# Regenerar documentación:
devtools::document("AACclustMac")

# Hacer la revision del paquete:
devtools::check("AACclustMac")

# Instalar el paquete:
devtools::install("AACclustMac")
