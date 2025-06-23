#!/usr/bin/env python
# coding: utf-8

# ## Imports y definición de parámetros

import tidy3d as td
from tidy3d import web 

import matplotlib.pyplot as plt
import numpy as np
import math
from dataclasses import dataclass


@dataclass

class Params:

    wvlength: float = 1.55 # Longitud de onda de trabajo

    # Splitter

    wvlength: float = 1.310 # Longitud de onda de trabajo (O band)

    W: float = 0.55 # Ancho de las guías 
    R: float = W/2 # Parámetro auxiliar para los vértices
    L: float = 2 # Longitud de las guías
    H: float = 0.22 # Altura de las guías
    G: float = 0.1/(np.sqrt(2)) # Hueco entre guías -> Normalizado para hacer el desplazamiento en x e y
    d: float = 0.05 # Distancia de separación para los vértices de la curva

    angle: float = np.pi/4  # Ángulo de reflexión
    angle_deg: float = (360/(2*np.pi)) * angle # Ángulo de reflexión en grados

    offset_x: float = 0.125       # Desplazamientos horizontales y verticales para tener
    offset_y: float = offset_x   # en cuenta el efecto Goos-Hänchen

    symmetric: bool = True  # Parámetro que controla si el conmutador es simétrico o no

    sup_shift: float = 0.1890  # Deslazamiento diagonal de la estructura superior

    # Capa de Nitruro de Silicio

    E: float = 2  # Espesor de la capa

    # SWG

    s: float = G # Ancho de la rejilla
    thickness: float = H # Altura de la rejilla
    unit_len: float = 0.05 # Longitud unitaria de los
    duty: float = wvlength/(2*2.49995398) # Ciclo de trabajo

param = Params()


# ## Estructuras

# Medios y frecuencia de trabajo

freq0 = td.C_0/param.wvlength
silicio = td.medium_from_nk(n=3.5, k=0, freq=freq0,name="Silicio")
sio2 = td.material_library["SiO2"]["Horiba"]
SiN = td.medium_from_nk(n=2.02, k=0, freq=freq0,name="Nitruro de Silicio")


def generate_splitter(param):

    f1 = 3.5; f2 = 2.5
    dx1 = f1 * param.offset_x;   dy1 = f1 * param.offset_y
    dx2 = f2 * param.offset_x;   dy2 = f2 * param.offset_y
    dd  = 6.5 * param.d

    # Vértices “base” de la parte inferior
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

    # Desplazamiento en diagonal para la parte superior
    dx_sup =  param.sup_shift * np.cos(param.angle)
    dy_sup = -param.sup_shift * np.sin(param.angle)

    # Vértices de la parte superior
    if param.symmetric:
        verts_sup = [
            (
                -x + param.G + dx_sup,
                -y + param.G + dy_sup
            )
            for x, y in verts_inf
        ]
    else:                     # Esta es la versión original no-simétrica
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


# Capa de nitruro de silicio

sin = td.Structure(
    geometry=td.Box(center=(0,0,(param.E+param.H)/2), size=(2.5*(param.L+param.R), 2.5*(param.L+param.R), param.E)),
    medium=SiN,
    name = "Nitruro de Silicio"
    )


def rotated_square(center, size, angle, thickness):
    """
    Devuelve un PolySlab cuadrado de tamaño `size`,
    rotado por `angle` [rad] y centrado en `center`.
    """
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
        axis=2,  # z-direction
        slab_bounds=(-thickness/2, thickness/2)
    )


def create_region_rotated(param, position, size, dutyX, dutyY, angle):

    geom_list = []

    nx, ny = tuple(map(lambda x: int(x/param.unit_len), size))

    # Tamaños iniciales y finales de cada elemento 
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


def generate_swg_rotated(param, angle):

    geometry_group = []

    region_params = [((param.G-0.05+param.offset_x, 3*param.G/4), ((np.sqrt(2))*param.W+6*param.unit_len, 2*param.G), param.duty, param.duty)]

    for pos, size, dX, dY in region_params:
        geometry_group += create_region_rotated(param, pos, size, dX, dY, angle)

    swg = td.Structure(
        geometry=td.GeometryGroup(geometries=geometry_group),
        medium=SiN,
        name="swg_rotated"
    )
    return swg

swg = generate_swg_rotated(param, angle=-param.angle)


# ## Monitores y fuente

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
    center=(param.L+param.R+param.G+-0.1, param.G, 0),
    size=(0, 2*param.W, 2*param.H),
    freqs=list(freqs),
    name = "flux_T"
)

flux_C = td.FluxMonitor(
    center=(3*param.offset_x, -param.L-param.R+0.1, 0),
    size=(2*param.W, 0, 2*param.H),
    freqs=list(freqs),
    name = "flux_C"
)

monitors = [field_monitor, flux_T, flux_C]


# ## Simulación

sim = td.Simulation(
    structures=splitter + [swg] + [sin],
    center=(0,0,1),
    size = (4*(param.L+param.R), 4*(param.L+param.R), 2*param.E), 
    medium=sio2,
    boundary_spec = td.BoundarySpec(
        x=td.Boundary.pml(),
        y=td.Boundary.pml(),
        z=td.Boundary.pml()
    ),
    sources=source,
    monitors=monitors,
    run_time=run_time,
    grid_spec=td.GridSpec.auto(min_steps_per_wvl=10)
)


fig, (ax1, ax2) = plt.subplots(1 , 2, figsize = (20,8))

sim.plot(x=-param.L-param.R+0.05, ax=ax1)
sim.plot(z=0, ax=ax2)
plt.show()


#sim.plot_3d()


task_name = f"cross_switch_def"


task_id = web.upload(sim, task_name=task_name, folder_name="TFG_PhotonicSwitch")


web.start(task_id)
web.monitor(task_id)

sim_data = web.load(task_id, path="data/sim.hdf5")


# ## Visualización de datos

fields = sim_data["fields"]

Ex = fields.Ex.isel(f=0, z=0)
Ey = fields.Ey.isel(f=0, z=0)
Ez = fields.Ez.isel(f=0, z=0)

E_abs = (abs(Ex)**2 + abs(Ey)**2 + abs(Ez)**2)**0.5

plt.figure(figsize=(6,5))
E_abs.plot(x="x", y="y", cmap="magma")
plt.title("Magnitud total del campo eléctrico |E|")
plt.tight_layout()
plt.show()


# ## Cálculo de las pérdidas y ratio de división

# Extraer flujos (P_C puede ser negativo si el sentido es contrario)
P_T = sim_data['flux_T'].flux.sel(f=freq0, method="nearest").values  # Through
P_C = sim_data['flux_C'].flux.sel(f=freq0, method="nearest").values  # Cross

# Tomar valor absoluto para tener potencias positivas
P_T = np.abs(P_T)
P_C = np.abs(P_C)

# Potencia total de salida
P_out_tot = P_T + P_C

# Pérdidas de inserción totales
IL_tot = -10 * np.log10(np.maximum(P_out_tot, 1e-12))

# Pérdidas de inserción por canal
IL_T = -10 * np.log10(np.maximum(P_T, 1e-12))
IL_C = -10 * np.log10(np.maximum(P_C, 1e-12))

# Pérdidas de exceso por canal
IL_ideal = 10 * np.log10(2)
EL_T = IL_T - IL_ideal
EL_C = IL_C - IL_ideal

# 4) Splitting ratio (fracción de la salida total)
splitting_ratio_T = P_T / P_out_tot
splitting_ratio_C = P_C / P_out_tot

# Mostrar resultados
print(f"P_out_tot (%):         {P_out_tot*100:.4f}%")
print(f"IL_total (dB):         {IL_tot:.4f}")
print(f"IL_Through (dB):       {IL_T:.4f}")
print(f"IL_Cross (dB):         {IL_C:.4f}")
print(f"EL_Through (dB):       {EL_T:.4f}")
print(f"EL_Cross (dB):         {EL_C:.4f}")
print(f"Splitting Through:     {splitting_ratio_T:.6f}")
print(f"Splitting Cross:       {splitting_ratio_C:.6f}")

