# TIR-Switch-TFG

Repositorio asociado al Trabajo de Fin de Grado (TFG):  
**Conmutadores fotónicos basados en reflexión total interna frustrada (FTIR) y rejillas sub-longitud de onda (SWG)**.

Este proyecto incluye el código fuente desarrollado en Python con el framework [Tidy3D](https://www.flexcompute.com/tidy3d/) para diseñar, simular y analizar dos dispositivos fotónicos:

- `CGC_swg_batch.py` / `CGC_swg.ipynb`:  
  Simula un CGC modificado con una región SWG, explorando el impacto del ciclo de trabajo (`duty`) en la transmisión óptica.  
  - **Batch script** (`.py`): barrido de `duty_range` y post-procesado.  
  - **Notebook** (`.ipynb`): versión interactiva con plots embebidos.

- `Cross_switch_swg_batch.py` / `Cross_switch_swg.ipynb`:  
  Simula un conmutador en cruz con SWG inclinada, analizando cómo varía la transmisión al cambiar el índice efectivo (`delta`).  
  - **Batch script** (`.py`): barrido de `delta` y post-procesado.  
  - **Notebook** (`.ipynb`): versión interactiva.

- `CGC_swg.py` y `Cross_switch_swg.py`:  
  Scripts auxiliares donde se extraen y calculan las **pérdidas de inserción**, **pérdidas de exceso** y **splitting ratios**, a partir de los flujos normalizados de Tidy3D.

---

## Estructura del repositorio

```plain
TIR-Switch-TFG/
├── CGC_swg_batch.ipynb           # Notebook diseño A (CGC + SWG)
├── CGC_swg_batch.py              # Script batch diseño A
├── CGC_swg.py                    # Cálculo de pérdidas y ratios diseño A
├── Cross_switch_swg_batch.ipynb  # Notebook diseño B (splitting + SWG)
├── Cross_switch_swg_batch.py     # Script batch diseño B
├── Cross_switch_swg.py           # Cálculo de pérdidas y ratios diseño B
├── README.md                     # Este archivo
```

Se recomienda crear un entorno virtual:

```
python -m venv tfg_env
source tfg_env/bin/activate    # Windows: .\tfg_env\Scripts\activate
pip install tidy3d matplotlib pandas numpy
```

* Python 3.8 o superior
* JupyterLab (opcional, para trabajar con los notebooks)
