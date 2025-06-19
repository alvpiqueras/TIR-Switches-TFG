#!/usr/bin/env python
# coding: utf-8

# ## Imports y definición de parámetros

# In[ ]:


import tidy3d as td
from tidy3d import web

import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass
import pandas as pd


# In[ ]:


@dataclass
class Params:

    wvlength: float = 1.55 # Longitud de onda de trabajo

    # Guías
    h: float = 2 # Separación entre I/O
    w: float = 0.55 # Anchura de las guías
    l: float = 5 # Largo de las guías
    L: float = l*0.2 # Longitud de acoplo
    W: float = w*0.7 # Ancho de acoplo
    a: float = 0.22 # Altura de las guías
    s: float = 0.2 # Separación entre estructuras
    angle: float = np.arctan((h+w/2)/l) # Ángulo mitad longitudinal
    angle_deg: float = angle*(360/(2*np.pi)) # Ángulo mitad longitudinal en grados

    # SWG
    unit_len: float = 0.05 # Longitud unitaria
    thickness: float = 0.22 # Altura de las guías
    duty: float = wvlength/(2*2.51610124) # Ciclo de trabajo

param = Params()


# In[ ]:


print(param.angle_deg)


# In[ ]:


print(param.duty)


# ## Generar las guías

# In[ ]:


from shapely.geometry import Polygon

def generate_guides(param):

    # Medio
    silicio = td.Medium(permittivity=12.25, name="Silicio")

    # Punto hacia donde convergen ambas ramas superiores
    center = np.array([0, param.W])
    half_w = param.w / 2

    def facet_endpoints(pin):

        u = center - pin
        u /= np.linalg.norm(u)
        # n es perpendicular a u
        n = np.array([-u[1], u[0]])
        return pin + n*half_w, pin - n*half_w

    # — FACET IZQUIERDO —
    pin_L = np.array([-(param.l + param.L/2), param.h + half_w])
    pL1, pL2 = facet_endpoints(pin_L)
    # Ordenamos por y -> (inferior, superior)
    low_L, high_L = sorted([pL1, pL2], key=lambda p: p[1])
    A, B = low_L, high_L

    # — FACET DERECHO —
    pin_R = np.array([ (param.l + param.L/2), param.h + half_w])
    pR1, pR2 = facet_endpoints(pin_R)
    # Ordenamos por y descendente -> (superior, inferior)
    high_R, low_R = sorted([pR1, pR2], key=lambda p: p[1], reverse=True)
    sup_R, inf_R = high_R, low_R

    # Vértices centrales de la X
    C = (0, param.W)
    F = ( param.L/2, 0)
    G = (-param.L/2, 0)

    # Montamos el perímetro en sentido horario:
    sup_vertices = [
        tuple(A), tuple(B), # facet izquierdo (inf→sup)
             C,             # centro
        tuple(sup_R),       # facet derecho (sup)
        tuple(inf_R),       # facet derecho (inf)
             F, G           # base
    ]

    # Si por algún redondeo aún se detectara autodiscreto, lo limpiamos:
    poly = Polygon(sup_vertices)
    if not poly.is_valid:
        poly = poly.buffer(0)
        sup_vertices = list(poly.exterior.coords)[:-1]

    sup = td.Structure(
        geometry=td.PolySlab(
            vertices = sup_vertices,
            axis = 2,
            slab_bounds = (-param.a/2, param.a/2),
        ),
        medium = silicio,
        name   = "sup_wv"
    )

    # Ramas inferiores (espejo en y y corrido por el gap s)
    inf_vertices = [(x, -y - param.s) for x, y in sup_vertices]
    inf = td.Structure(
        geometry=td.PolySlab(
            vertices = inf_vertices,
            axis = 2,
            slab_bounds = (-param.a/2, param.a/2),
        ),
        medium = silicio,
        name   = "inf_wv"
    )

    return [sup, inf]

guides = generate_guides(param)


# ## Generar la rejilla sub‑longitud de onda (SWG)

# In[ ]:


def create_region(param, position, size, dutyX, dutyY):

    # param: acceso a los parámetros
    # position: tupla que indica dónde centrar la rejilla
    # size: tamaño de la rejilla
    # dutyX, dutyY: ciclos de trabajo en X e Y (porcentaje de la celda unidad vacío). Si es una tupla se hace una variación progresiva

    # En esta lista guardamos el conjunto de las estructuras
    geom_list = [] 

    # Número de elementos en cada dirección
    nx, ny = tuple(map(lambda x: int(x/param.unit_len), size))

    # Tamaños iniciales y finales
    # En X:
    if isinstance(dutyX, tuple):
        L0, L1 = tuple(map(lambda x: param.unit_len * (1-x), dutyX))
    else:
        L0, L1 = tuple(map(lambda x: param.unit_len * (1-x), (dutyX, dutyX)))

    # En Y:    
    if isinstance(dutyY, tuple):
        H0, H1 = tuple(map(lambda x: param.unit_len * (1-x), dutyY))
    else:
        H0, H1 = tuple(map(lambda x: param.unit_len * (1-x), (dutyY, dutyY)))

    # Iteramos en cada dirección
    for j in range(ny):
        for i in range(nx):
            xpos = position[0] + (i+0.5-0.5*nx) * param.unit_len # Partiendo de la posición inicial, se centra en la celda unidad y 
            ypos = position[1] + (j+0.5-0.5*ny) * param.unit_len # se empieza a rellenar de izquierda a derecha

            L = L0 + (i/nx) * (L1-L0) # Largo (X) y ancho (Y) de los cuadrados. Si no hay variación de ciclo de trabajo,
            H = H0 + (i/nx) * (H1-H0) # entonces L1=L0; H1=H0 y los bloques son constantes.

            geom = td.Box(center=(xpos, ypos, 0), size=(L, H, param.thickness))
            geom_list.append(geom)

    return geom_list


def generate_swg(param, dutyX, dutyY):

    # Material
    silicio = td.Medium(permittivity=(3.5)**2, name="Silicio")

    # Lista de rejillas. En este caso solo haremos una
    geometry_group = []

    # Parámetros que pasar a la función "create_region" (position, size, dutyX, dutyY)
    region_params = [((0, -param.s/2), (param.L, param.s), dutyX, dutyY)]

    for pos, size, dX, dY in region_params:
        geometry_group += create_region(param, pos, size, dX, dY)

    swg = td.Structure(
        geometry=td.GeometryGroup(geometries=geometry_group),
        medium=silicio,
        name="subwavelength_grating"
    )
    return swg


# ## Monitores y fuentes

# In[ ]:


lambda0 = param.wvlength # Longittud de onda de trabajo
freq0   = td.C_0/lambda0 # Frecuencia de trabajo
fwidth  = freq0/10 # Ancho de pulso gaussiano
run_time = 40/fwidth # Tiempo de simulación en base a la duración del pulso

# Pequeño linspace de frecuencias
Nfreqs = 5
freqs  = np.linspace(freq0 - fwidth, freq0 + fwidth, Nfreqs)

# Extremos de las guías izquierda y derecha
pin_L = np.array([-(param.l + param.L/2),  param.h + param.w/2])
pin_R = np.array([ (param.l + param.L/2),  param.h + param.w/2])

source = [
    td.ModeSource(
        center= (pin_L[0] + 0.15, pin_L[1] - 0.05, 0),
        size=   (0, 2*param.w, 3*param.a),
        source_time=td.GaussianPulse(freq0=freq0, fwidth=fwidth),
        direction="+",
        mode_index=0,
        mode_spec=td.ModeSpec(num_modes=1, angle_theta=-param.angle),
        name="source"
    )
]

fields = td.FieldMonitor(
    fields=["Ex","Ey","Ez"],
    center=(0, 0, 0),
    size=(2.5*param.l, 3.2*(param.h+param.w/2), 0),
    freqs=[freq0],
    name="fields"
)

flux_sup = td.FluxMonitor(
    center=(pin_R[0] - 0.15, pin_R[1] - 0.1, 0),
    size=(0, 2*param.w, 3*param.a),
    freqs=list(freqs),
    name="flux_sup"
)

flux_inf = td.FluxMonitor(
    center=(pin_R[0] - 0.15, -pin_R[1] - param.s + 0.1, 0),
    size=(0, 2*param.w, 3*param.a),
    freqs=list(freqs),
    name="flux_inf"
)

flux_back = td.FluxMonitor(
    center=(pin_L[0] + 0.15, -pin_L[1] - param.s + 0.1, 0),
    size=(0, 2*param.w, 3*param.a),
    freqs=list(freqs),
    name="flux_back"
)

monitors = [fields, flux_back, flux_sup, flux_inf]


# ## Lote de simulaciones

# In[ ]:


sio2 = td.material_library["SiO2"]["Horiba"]

N_duty = 15
duty_range = np.linspace(0.5, 0.8, N_duty)  

sims = {}  

for i in range(N_duty):

    swg = generate_swg(param, dutyX= duty_range[i], dutyY= duty_range[i])

    sim = td.Simulation(
        structures = guides + [swg],
        size=(3*param.l, 6*(param.h+param.w/2), 10*param.a),
        medium = sio2,
        boundary_spec = td.BoundarySpec(
            x=td.Boundary.pml(),
            y=td.Boundary.pml(),
            z=td.Boundary.pml()
        ),
        sources = source,
        monitors = monitors,
        run_time = run_time,
        grid_spec = td.GridSpec.auto(min_steps_per_wvl=10)
    )

    sims[f"duty_{duty_range[i]:.5f}"] = sim


# In[ ]:


print(list(sims))


# In[ ]:


batch = web.Batch(simulations=sims)
batch.estimate_cost()


# In[ ]:


batch_results = batch.run(path_dir="data_duty_sweep")


# ## Visualización de resultados

# In[ ]:


flujos_sup = {}
flujos_inf = {}
flujos_back = {}

for task_name, sim_data in batch_results.items():

    flujo_sup = sim_data["flux_sup"].flux.sel(f=freq0, method="nearest").values
    flujo_inf = sim_data["flux_inf"].flux.sel(f=freq0, method="nearest").values
    flujo_back = sim_data["flux_back"].flux.sel(f=freq0, method="nearest").values

    flujos_sup [task_name] = flujo_sup
    flujos_inf [task_name] = flujo_inf
    flujos_back [task_name] = flujo_back       


# In[ ]:


# Crear lista de registros
records = []

for task_name in flujos_sup:
    duty = float(task_name.split("_")[1])
    records.append({
        "duty": duty,
        "flux_sup": flujos_sup[task_name].item(),
        "flux_inf": flujos_inf[task_name].item(),
        "flux_back": flujos_back[task_name].item()
    })

# Crear DataFrame y ordenar
df = pd.DataFrame(records).sort_values("duty")

# Obtener máximos
row_max_sup  = df.loc[df["flux_sup"].idxmax()]
row_max_inf  = df.loc[df["flux_inf"].idxmax()]
row_max_back = df.loc[df["flux_back"].idxmax()]

# Imprimir resultados
print(f"Máximo flujo superior: {row_max_sup['flux_sup']:.5f} en f = {row_max_sup['duty']:.5f}")
print(f"Máximo flujo inferior: {row_max_inf['flux_inf']:.5f} en f = {row_max_inf['duty']:.5f}")
print(f"Máximo flujo back:     {row_max_back['flux_back']:.5f} en f = {row_max_back['duty']:.5f}")

# Graficar
plt.figure(figsize=(8, 5))
plt.plot(df["duty"], df["flux_sup"], label="Flujo Superior")
plt.plot(df["duty"], df["flux_inf"], label="Flujo Inferior")
plt.plot(df["duty"], df["flux_back"], label="Flujo Back")

# Marcar máximos
plt.scatter(row_max_sup["duty"], row_max_sup["flux_sup"], color="blue", marker="o", label="Máx. Superior")
plt.scatter(row_max_inf["duty"], row_max_inf["flux_inf"], color="orange", marker="o", label="Máx. Inferior")
plt.scatter(row_max_back["duty"], row_max_back["flux_back"], color="green", marker="o", label="Máx. Back")

# Etiquetas y estilo
plt.xlabel("f")
plt.ylabel("Flujo (1550 nm)")
plt.title("Flujos vs f")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

