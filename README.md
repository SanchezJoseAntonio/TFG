# Scripts utilizados en el TFG: Análisis computacional de datos de secuenciación de iCLIP
## Jerarquía circAnno
El programa circAnno (https://github.com/sysu-software/circAnno) es una herramienta que permite la anotación de ARN circular con anotaciones conocidas. En este TFG, se utiliza este programa para anotar datos de PAR-CLIP. A cada tipo de anotación (ARN transferente, 3' UTR, 5' UTR ... ) se le asigna una prioridad distinta, por lo que es necesario que se le proporcionen los archivos a este programa de manera jerarquizada. Con este fin se utiliza este _script_. 
## Anotación de lecturas utilizando código en Python (Totalmente solapadas)
Una forma previa de llevar a cabo la anotación fue utilizando código en Python, aquí proporcionado. En este caso, solo se consideraron las lecturas totalmente solapadas con la anotación como anotadas, y el resto como intergénicas.
## Anotación de lecturas utilizando código en Python (Totalmente y parcialmente solapadas)
Otra forma previa de llevar a cabo la anotación fue utilizando código en Python, aquí proporcionado. En este caso, solo se consideraron las lecturas parcialmente solapadas se consideraron parcialmente intergénicas (% de lectura fuera de la anotación) y parcialmente anotadas (% de lectura dentro de la anotación).
## Archivo de distribución de lecturas
Este archivo contiene el código que permite obtener la distribución de las lecturas en las anotaciones según el porcentaje de la lectura necesario para que se considere anotada y no como intergénico.  
**Todos los _scripts_ escritos en Python requieren que las coordenadas genómicas de los archivos que se pasen como input estén ordenados numéricamente**
