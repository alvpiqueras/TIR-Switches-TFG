#!/usr/bin/env python
# coding: utf-8

# ## Imports y definición de parámetros

# In[ ]:


import tidy3d as td
from tidy3d import web 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass


# In[ ]:


@dataclass

class Params:

    wvlength: float = 1.310 # Longitud de onda de trabajo (O band)

    # Splitter

    W: float = 0.55 # Ancho de las guías 
    R: float = W/2 # Parámetro auxiliar para los vértices
    L: float = 2 # Longitud de las guías
    H: float = 0.22 # Altura de las guías
    G: float = 0.1/(np.sqrt(2)) # Hueco entre dispositivos -> Normalizado para hacer el desplazamiento en x e y
    d: float = 0.05 # Distancia de separación para los vértices de la curva

    angle: float = np.pi/4  # Ángulo de reflexión
    angle_deg: float = (360/(2*np.pi)) * angle # Ángulo de reflexión en grados

    offset_x: float = 0.125      # Desplazamientos horizontales y verticales para tener
    offset_y: float = offset_x   # en cuenta el efecto Goos-Hänchen

    symmetric: bool = True  # Parámetro que controla si el conmutador es simétrico o no

    sup_shift: float = 0.1890 # Deslazamiento diagonal de la estructura superior

    # Capa de Nitruro de Silicio

    E: float = 2  # Espesor de la capa

    # SWG

    s: float = G # Ancho de la rejilla
    thickness: float = H # Altura de la rejilla
    unit_len: float = 0.05 # Longitud unitaria
    duty: float = wvlength/(2*2.49995398) # Ciclo de trabajo

param = Params()


# ## Estructuras

# In[ ]:


# Medios y frecuencia de trabajo

freq0 = td.C_0/param.wvlength
silicio = td.medium_from_nk(n=3.5, k=0, freq=freq0,name="Silicio")
sio2 = td.material_library["SiO2"]["Horiba"]
SiN = td.medium_from_nk(n=2.02, k=0, freq=freq0,name="Nitruro de Silicio")


# In[ ]:


def generate_splitter(param):

    # 1) Factores de desplazamiento 
    f1 = 3.5; f2 = 2.5
    dx1 = f1 * param.offset_x;   dy1 = f1 * param.offset_y
    dx2 = f2 * param.offset_x;   dy2 = f2 * param.offset_y
    dd  = 6.5 * param.d

    # 2) Vértices “base” de la parte inferior
    xA, yA = -param.L-param.R,  param.R
    xB, yB = -param.R,          param.R
    xC, yC =  param.R,         -param.R
    xD, yD =  param.R,         -param.L-param.R
    xE, yE = -param.R,         -param.L-param.R
    xF, yF = -param.R,         -param.R-2*param.d
    xG, yG = -param.R-2*param.d,-param.R
    xH, yH = -param.L-param.R, -param.R

    verts_inf = [
        (xA,       yA),
        (xB,       yB),
        (xC+dx1,   yC-dy1),
        (xD+dx1,   yD),
        (xE+dx2,   yE),
        (xF+dx2,   yF-dd),
        (xG-param.d, yG),
        (xH,        yH),
    ]

    inf = td.Structure(
        geometry=td.PolySlab(
            vertices    = verts_inf,
            slab_bounds = (-param.H/2, param.H/2),
            axis        = 2,
        ),
        medium = silicio,
        name   = "inf"
    )

    # 3) Desplazamiento en diagonal para la parte superior
    dx_sup =  param.sup_shift * np.cos(param.angle)
    dy_sup = -param.sup_shift * np.sin(param.angle)

    # 4) Vértices de la parte superior
    if param.symmetric:
        verts_sup = [
            (
                -x + param.G + dx_sup,
                -y + param.G + dy_sup
            )
            for x, y in verts_inf
        ]
    else:                                                               # Esta es la versión original no-simétrica
        verts_sup = [
            (-xA+param.G+param.offset_x, -yA+param.G-param.offset_y),
            (-xB+param.G+param.offset_x, -yB+param.G-param.offset_y),
            (-xC+param.G, -yC+param.G),                                 # No se mueven estos vértices para
            (-xD+param.G, -yD+param.G-param.offset_y),                  # mantener la proporción
            (-xE+param.G+param.offset_x, -yE+param.G-param.offset_y),
            (-xF+param.G+param.offset_x, -yF+param.G-param.offset_y+param.d),
            (-xG+param.G+param.offset_x+param.d, -yG+param.G-param.offset_y),
            (-xH+param.G+param.offset_x, -yH+param.G-param.offset_y)
        ]

    sup = td.Structure(
        geometry=td.PolySlab(
            vertices    = verts_sup,
            slab_bounds = (-param.H/2, param.H/2),
            axis        = 2,
        ),
        medium = silicio,
        name   = "sup"
    )

    return [inf, sup]

splitter = generate_splitter(param)


# In[ ]:


# Capa de nitruro de silicio

sin = td.Structure(
    geometry=td.Box(center=(0,0,(param.E+param.H)/2), size=(2.5*(param.L+param.R), 2.5*(param.L+param.R), param.E)),
    medium=SiN,
    name = "Nitruro de Silicio"
    )


# In[ ]:


def rotated_square(center, size, angle, thickness):

    cx, cy = center
    half = size / 2
    c, s = np.cos(angle), np.sin(angle)

    # Coordenadas del cuadrado no rotado
    base_vertices = [
        (-half, -half),
        ( half, -half),
        ( half,  half),
        (-half,  half)
    ]

    # Rotar cada vértice
    rotated_vertices = [
        (cx + c*x - s*y, cy + s*x + c*y) for x, y in base_vertices
    ]

    return td.PolySlab(
        vertices=rotated_vertices,
        axis=2,  
        slab_bounds=(-thickness/2, thickness/2)
    )


# In[ ]:


def create_region_rotated(param, position, size, dutyX, dutyY, angle):

    # param: acceso a los parámetros
    # position: tupla que indica dónde centrar la rejilla
    # size: tamaño de la rejilla
    # dutyX, dutyY: ciclos de trabajo en X e Y (porcentaje de la celda unidad vacío). Si es una tupla se hace una variación progresiva
    #angle: cuánto se rota la estructura

    geom_list = []

    nx, ny = tuple(map(lambda x: int(x/param.unit_len), size))

    # Initial and final element sizes
    if isinstance(dutyX, tuple):
        L0, L1 = tuple(map(lambda x: param.unit_len * (1-x), dutyX))
    else:
        L0, L1 = tuple(map(lambda x: param.unit_len * (1-x), (dutyX, dutyX)))

    if isinstance(dutyY, tuple):
        H0, H1 = tuple(map(lambda x: param.unit_len * (1-x), dutyY))
    else:
        H0, H1 = tuple(map(lambda x: param.unit_len * (1-x), (dutyY, dutyY)))

    c, s = np.cos(angle), np.sin(angle)

    for j in range(ny):
        for i in range(nx):

            # Posición sin rotar
            x0 = position[0] + (i+0.5-0.5*nx) * param.unit_len
            y0 = position[1] + (j+0.5-0.5*ny) * param.unit_len

            # Rotación 2D alrededor del origen 
            xr =  c*x0 - s*y0
            yr =  s*x0 + c*y0

            L = L0 + (i/nx)*(L1-L0)
            H = H0 + (i/nx)*(H1-H0)

            geom = rotated_square(
                    center=(xr, yr),
                    size=param.unit_len * (1 - dutyX), 
                    angle=angle,
                    thickness=param.thickness
                )
            geom_list.append(geom)

    return geom_list


# In[ ]:


def generate_swg_rotated(param, angle, delta):

    silicio_delta = td.medium_from_nk(n=(3.5-delta), k=0, freq=freq0,name="Silicio variable")

    geometry_group = []

    # Definición de los parámetros para la rejilla
    region_params = [((param.G-0.05, 3*param.G/4), ((np.sqrt(2))*param.W+6*param.unit_len, 2*param.G), param.duty, param.duty)]

    # Creación de los elementos individuales de la rejilla, rotándolos en el proceso
    # Quedan almacenados en "geometry_group"
    for pos, size, dX, dY in region_params:
        geometry_group += create_region_rotated(param, pos, size, dX, dY, angle)

    # Creación de la rejilla en base a los elementos individuales    
    swg = td.Structure(
        geometry=td.GeometryGroup(geometries=geometry_group),
        medium=silicio_delta,
        name="swg_rotated"
    )
    return swg


# ## Monitores y fuente

# In[ ]:


fwidth = freq0 / 10
run_time = 40 / fwidth

Nfreqs = 5
freqs = np.linspace(freq0-fwidth, freq0+fwidth, Nfreqs)

source = [
        td.ModeSource(
        center = (-param.L-param.R+0.05, 0, 0),
        size = (0, 2*param.W, 1.5*param.H),
        source_time = td.GaussianPulse(freq0=freq0,fwidth=fwidth),
        direction = "+",
        mode_index = 0,
        mode_spec = td.ModeSpec(num_modes=1),
        name = "source"
        )
    ]

field_monitor = td.FieldMonitor(
    fields=["Ex", "Ey", "Ez"],
    center = (0, 0, 0),
    size = (2*(param.L+param.R+param.G+param.offset_x),2*(param.L+param.R+param.G+param.offset_y), 0),
    freqs=[freq0],
    name = "fields"
)

flux_T = td.FluxMonitor(
    center=(param.L+param.R+param.G+param.offset_x-0.1, param.G-param.offset_y, 0),
    size=(0, param.W, 2*param.H),
    freqs=list(freqs),
    name = "flux_T"
    )

flux_C = td.FluxMonitor(
    center=(0, -param.L-param.R+0.1, 0),
    size=(param.W+2*param.offset_x, 0, 2*param.H),
    freqs=list(freqs),
    name = "flux_C"
    )

monitors = [field_monitor, flux_T, flux_C]


# ## Simulación

# In[ ]:


N_delta = 15                                    
delta_range = np.linspace(0.005, 0.01, N_delta) # Rango del parámetro y cantidad de pasos

sims = {}  # Para almacenar las simulaciones

# Bucle para la generación de simulaciones en base a un parámetro

for i in range(N_delta):

    swg = generate_swg_rotated(param, angle=param.angle , delta=delta_range[i])

    sim = td.Simulation(
        structures = splitter + [swg] + [sin],
        size=(4*(param.L+param.R), 4*(param.L+param.R), 2*param.E),
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

    sims[f"delta_{delta_range[i]:.5f}"] = sim


# In[ ]:


print(list(sims))


# In[ ]:


batch = web.Batch(simulations=sims)
batch.estimate_cost()


# In[ ]:


batch_results = batch.run(path_dir="splitter_delta_sweep")


# ## Visualización de datos

# In[ ]:


flujos_T = {}
flujos_C = {}

for task_name, sim_data in batch_results.items():

    flujo_T = sim_data["flux_T"].flux.sel(f=freq0, method="nearest").values
    flujo_C = sim_data["flux_C"].flux.sel(f=freq0, method="nearest").values

    flujos_T [task_name] = flujo_T
    flujos_C [task_name] = flujo_C


# In[ ]:


# Crear lista de registros
records = []

for task_name in flujos_T:
    delta = float(task_name.split("_")[1])
    records.append({
        "delta": delta,
        "flux_T": flujos_T[task_name].item(),
        "flux_C": flujos_C[task_name].item()
    })

# Crear DataFrame y ordenar
df = pd.DataFrame(records).sort_values("delta")

# Obtener máximos
row_max_T = df.loc[df["flux_T"].idxmax()]
row_max_C = df.loc[df["flux_C"].idxmax()]

# Imprimir resultados
print(f"Máximo flujo paralelo: {row_max_T['flux_T']:.5f} en delta = {row_max_T['delta']:.5f}")
print(f"Máximo flujo cruzado:  {row_max_C['flux_C']:.5f} en delta = {row_max_C['delta']:.5f}")

# Graficar
plt.figure(figsize=(8, 5))
plt.plot(df["delta"], df["flux_T"], label="Flujo Paralelo")
plt.plot(df["delta"], df["flux_C"], label="Flujo Cruzado")

# Marcar máximos
plt.scatter(row_max_T["delta"], row_max_T["flux_T"], color="blue", marker="o", label="Máx. Paralelo")
plt.scatter(row_max_C["delta"], row_max_C["flux_C"], color="orange", marker="o", label="Máx. Cruzado")

# Etiquetas y estilo
plt.xlabel(r"$\delta$")
plt.ylabel("Flujo (1310 nm)")
plt.title(r"Flujo vs $\delta$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

