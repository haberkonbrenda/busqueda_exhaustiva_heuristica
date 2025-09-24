# Búsqueda Exhaustiva / Heurística / Greedy
Implementación de búsqueda exhaustiva y heurística para posicionamiento de pieza

Script para comparar **búsqueda exhaustiva** vs **búsqueda heurística (A\*)** en el problema de relocalizar **A** desde **B** sobre la horizontal **H**, con **rotación θ opcional**. 

---

## 1) Requisitos

- Python **3.8+**
- (Opcional, para gráficos) `matplotlib`

Instalación rápida:

```bash
python3 -m pip install matplotlib
# En Linux, puede requerir Tk para abrir la ventana de gráficos:
# Ubuntu/Debian:
sudo apt-get install python3-tk
```

---

## 2) ¿Qué implementa?

### Modelo (espacio de estados)
- **Estado:** `(x, θ)`  
  - `x`: posición sobre la horizontal **H** (si no se activa rotación, solo `x`).  
  - `θ`: giro (grados) alrededor de Z, **opcional**.
- **Operadores:** `x ± ΔH` y, si hay rotación, `θ ± Δθ`.
- **Meta:** `x = x_A` (y `θ = θ_A` si corresponde).
- **Costo:** 1 por movimiento (intentos/unidad de avance).

### Métodos
- **Exhaustiva radial** (BFS por niveles)  
  - 1D: explora `+k·ΔH, −k·ΔH` alrededor de B.  
  - Con rotación: capas **Manhattan** alrededor de `(x_B, θ_B)` (forma “diamante”).  
  - ✅ Simple/robusta · ❌ Puede hacer muchos intentos.
- **A\*** (`f = g + h`) con **h = distancia Manhattan**  
  - 1D: `h = |x−x_A|/ΔH`  
  - 2D: `h = |x−x_A|/ΔH + |θ−θ_A|/Δθ`  
  - ✅ Con esta h admisible/consistente, A\* **minimiza** el nº de pasos/intento.
- **Greedy Best-First** (opcional)  
  - Ordena solo por **h** (ignora `g`).  
  - ✅ Suele ir directo · ❌ No garantiza óptimo.

### Métricas impresas
- **Path / orden de intentos:** secuencia de estados visitados hasta la meta.  
- **N:** nº de intentos.  
- **d:** distancia mínima (Manhattan) en pasos de B a A.  
- **P = d/N:** **factor de dispersión** (1.0 es ideal). 

---

## 3) Uso

```bash
python tp2_busqueda.py --method {exhaustive|astar|greedy} [opciones]
python tp2_busqueda.py --compare [opciones]
```

- `--method`: ejecuta **un** método.  
- `--compare`: corre **los tres** métodos en **la misma instancia** (misma B→A).

### Opciones comunes

| Opción         | Tipo | Default | Descripción |
|---|---:|---:|---|
| `--seed`       | int  | —       | Semilla para reproducir la instancia aleatoria (elige A). |
| `--max-offset` | int  | 20      | Máximo corrimiento \|A−B\| en **unidades de paso** (sobre x). |
| `--step`       | int  | 1       | Tamaño de paso **ΔH** (avance en x). Debe ser > 0. |
| `--verbose`    | flag | false   | Muestra cada intento/expansión. |
| `--plot`       | flag | false   | Grafica la trayectoria (requiere matplotlib). |

### Opciones de rotación (opcionales)

| Opción               | Tipo | Default | Descripción |
|---|---:|---:|---|
| `--enable-rotation`  | flag | false   | Activa dimensión **θ**. |
| `--theta-step`       | int  | 1       | Paso angular **Δθ** (grados). Debe ser > 0. |
| `--theta-maxdeg`     | int  | 10      | Al muestrear la meta, límite de \|θ_A\| en grados. |

> Si **no** usás `--enable-rotation`, el problema se resuelve en 1D (solo `x`).

---

---

## 4) Ejemplos

### A) Una corrida por método 

```bash
# Exhaustiva sin rotación
python tp2_busqueda.py --method exhaustive --seed 42 --step 1 --max-offset 20 --verbose --plot

# A* sin rotación (misma instancia que arriba si repetís seed y parámetros)
python tp2_busqueda.py --method astar --seed 42 --step 1 --max-offset 20 --verbose --plot

# Greedy (opcional) sin rotación
python tp2_busqueda.py --method greedy --seed 42 --step 1 --max-offset 20 --verbose --plot
```

### B) Comparar métodos en la misma instancia

```bash
python tp2_busqueda.py --compare --seed 123 --step 1 --max-offset 20 --verbose --plot
```

### C) Con rotación θ (2D)

```bash
# Exhaustiva con rotación
python tp2_busqueda.py --method exhaustive --enable-rotation --theta-step 1 --theta-maxdeg 10 --seed 7 --verbose --plot

# A* con rotación
python tp2_busqueda.py --method astar --enable-rotation --theta-step 1 --theta-maxdeg 10 --seed 7 --verbose --plot

# Comparar los tres en la misma instancia con rotación
python tp2_busqueda.py --compare --enable-rotation --theta-step 1 --theta-maxdeg 10 --seed 7 --verbose --plot
```

---

## 5) Interpretación de la salida

- **Path / orden de intentos:** p.ej. `[(0,0), (1,0), (2,0), …, (7,0)]`.  
- **N:** nº de intentos (tamaño del path).  
- **d:** pasos mínimos (si B=0 y A=7 con ΔH=1 ⇒ d=7; si además θ_B=0 y θ_A=−3 con Δθ=1 ⇒ d=10).  
- **P = d/N:** eficiencia de direccionamiento (Exhaustiva→ P<1; A\*→ tiende a P=1).

--- 

## 6) Troubleshooting

- **“matplotlib no disponible; no se graficará”**  
  Instalar para Python 3:
  ```bash
  python3 -m pip install matplotlib
  # En Linux (backend de ventanas):
  sudo apt-get install python3-tk
  ```

- **Reproducibilidad:**  
  Usar `--seed` para que cada método use la **misma instancia** B→A y puedas comparar **N, d, P** de forma justa.

---
