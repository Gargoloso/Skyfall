# Skyfall

Guanajuato Release 0.1 
Última actualización: 26/06/2018

Este repositorio contiene los scripts desarrollados y probados sobre los cubos de CALIFA para el análisis de líneas de emisión.
Los diferentes scripts deben correrse en un órden determinado para crear adecuadamente los archivos de salida que otros scripts
usarán como entrada.

Suponiendo que extrajimos las líneas de emisión con la última versión de los scripts de Marcel, el primer paso es desenrojecer 
las líneas. Podemos usar tanto la ley de extinción de Cardelli et al. (1989) como la de Calzetti et al. (2000) para galaxias 
starburst. Se pueden usar ambas a la vez ya que los resultados se guardan en archivos separados. Los scripts de cada ley indican
cual es en su nombre. Debe correrse obligatoriamente el script Base que desenrojece las líneas H-alfa, H-beta, [N II] 6548,6584
y [O III] 5007. Las demás líneas ([O I] 6300, [O III] 4363 y [S II] 6716,6731) se corren en scripts separados y son opcionales.

Una vez corrido el script Base, podemos hacer un mapa de la extinción en la línea de H-alfa A(H-alfa) en magnitudes usando los
scripts Cardelli_AHa_Maps.r y Calzetti_AHa_Maps.r.

Diagramas BPT.

Por ahora solo está disponible el script para los BPT-NII que se puede correr con el output de cualquiera de los dos scritp Base 
de corrección de las líneas de extinción.

Por favor, consultar la guía paso a paso para utilizar adecuadamente estos scripts.
