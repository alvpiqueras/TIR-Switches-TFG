# TIR-Switch-TFG

Repositorio asociado al Trabajo de Fin de Grado (TFG):  
**Conmutadores fotónicos basados en reflexión total interna frustrada (FTIR) y rejillas sub-longitud de onda (SWG)**.

Este proyecto incluye el código fuente desarrollado en Python con el framework [Tidy3D](https://www.flexcompute.com/tidy3d/) para diseñar y simular dos dispositivos fotónicos:

- `CGC_swg_batch.py`: simula un CGC modificado con una región SWG, explorando el impacto del ciclo de trabajo `f` en la transmisión óptica.
- `Cross_switch_swg_batch.py`: simula un conmutador en cruz (tipo splitter), donde se analiza la influencia de variaciones del índice de refracción en el rendimiento del dispositivo.

------------------------------
## Estructura del repositorio ##

```

TIR-Switch-TFG/
├── CGC_swg_batch.ipynb # Notebook del diseño A (CGC)
├── Cross_switch_swg_batch.ipynb # Notebook del diseño B (Splitter)
├── CGC_swg_batch.py # Versión exportada a .py del CGC
├── Cross_switch_swg_batch.py # Versión exportada a .py del Splitter
├── README.md # Este archivo

```


-----------------------------
## Requisitos ##

Es recomendable usar un entorno virtual:

```bash
python -m venv tfg_env
source tfg_env/bin/activate   # en Windows: .\tfg_env\Scripts\activate
pip install tidy3d matplotlib pandas numpy
```

Si estás en JupyterLab, también puedes instalar las dependencias allí.
Se sugiere usar Python 3.8 o superior.

-------------------------------
## ¿Qué hace cada archivo? ##

CGC_swg_batch.py
Simula un acoplador CGC con una región SWG transversal. Se analiza cómo la variación del ciclo de trabajo 
f afecta la división del flujo óptico entre los puertos superior, inferior y reflejado. Se exploran valores de 
f entre 0.5 y 0.8.

Cross_switch_swg_batch.py
Simula un conmutador en cruz con SWG inclinada, analizando el efecto de cambiar el índice de refracción efectivo (mediante el parámetro 
"delta"). Se exploran valores de delta entre 0.005 y 0.01, imitando una posible modulación por dispersión de portadores.

