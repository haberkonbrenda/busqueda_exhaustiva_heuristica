#!/usr/bin/env python3

"""
TP2 - Búsqueda exhaustiva y heuristica

--> VER README <-- 

Descripción:
- Simula el problema B -> A sobre una línea (horizontal H) con paso fijo ΔH y, opcionalmente, un ángulo θ (giro en Z).
- Implementa:
    * Exhaustiva radial 1D o 2D (capas de distancia Manhattan alrededor de B)
    * A* (f = g + h) con h(n) = distancia Manhattan a la meta (admisible/consistente)
    * Greedy Best-First (opcional) con h(n) = distancia Manhattan
- Permite comparar métodos en la misma instancia (misma semilla/objetivo).
- Mide: N (intentos), d (distancia mínima en pasos), P = d/N (factor de dispersión).
- Opcional: gráficos de trayectoria; si hay rotación, se grafican x y θ en figuras separadas.

Uso (ejemplos):
    python tp2_busqueda.py --method exhaustive --seed 42 --max-offset 20 --step 1 --verbose
    python tp2_busqueda.py --method astar --seed 42 --max-offset 20 --step 1 --verbose --plot
    python tp2_busqueda.py --compare --seed 123

    # Con rotación (θ):
    python tp2_busqueda.py --compare --enable-rotation --theta-step 1 --theta-maxdeg 10 --seed 2024 --verbose --plot

    
"""
from __future__ import annotations
import argparse #elegir metodo, semilla, pasos
import math
import random 
import heapq #cola prioridad
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional, Set

try:
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False

#Nodo para cola de prioridad A* -> se ordena por f 
@dataclass(order=True)
class PQNode:
    f: int
    x: int = field(compare=False)
    theta: int = field(compare=False)
    g: int = field(compare=False)
    parent: Optional['PQNode'] = field(compare=False, default=None)

#Resultado de corrida, metodo B/A, path y métricas
@dataclass
class RunResult:
    method: str
    x_B: int
    x_A: int
    step: int
    path: List[Tuple[int,int]]  # sequence of (x, theta)
    expansions: int
    theta_enabled: bool
    theta_B: int
    theta_A: int
    theta_step: int

    @property
    def d(self) -> int:
        # Distancia mínima en pasos (Manhattan) entre B y A.
        dx = abs(self.x_A - self.x_B) // self.step
        if not self.theta_enabled:
            return dx
        dtheta = abs(self.theta_A - self.theta_B) // self.theta_step
        return dx + dtheta

    @property
    def N(self) -> int:
        #Número de intentos / nodos visitados (path incluye la meta).
        return len(self.path)

    @property
    def P(self) -> float:
        #Factor de dispersión: d / N.
        return (self.d / self.N) if self.N > 0 else float('nan')


class AlignmentProblem:
    """
    Modelo B -> A:
      - Estados: (x, theta) con theta opcional
      - Operadores: x ± ΔH; si hay rotación, además θ ± Δθ
      - Meta: x == x_A y (si aplica) θ == θ_A
    """
    def __init__(self, x_B: int, x_A: int, step: int = 1,
                 theta_enabled: bool=False, theta_B: int=0, theta_A: int=0, theta_step: int=1):
        assert step > 0, "step (ΔH) debe ser positivo"
        assert theta_step > 0, "theta_step debe ser positivo"
        # Normalizar a múltiplos de 'step' y 'theta_step'
        def snap(v: int, s: int) -> int:
            return int(round(v / s)) * s

        self.x_B = snap(x_B, step)
        self.x_A = snap(x_A, step)
        self.step = step

        self.theta_enabled = theta_enabled
        self.theta_B = snap(theta_B, theta_step)
        self.theta_A = snap(theta_A, theta_step)
        self.theta_step = theta_step

    def is_goal(self, x: int, theta: int) -> bool:
        if not self.theta_enabled:
            return x == self.x_A
        return x == self.x_A and theta == self.theta_A

    def heuristic(self, x: int, theta: int) -> int:
        # Manhattan discreta en la grilla (admisible/consistente con coste unitario por movimiento).
        hx = abs(self.x_A - x) // self.step
        if not self.theta_enabled:
            return hx
        ht = abs(self.theta_A - theta) // self.theta_step
        return hx + ht

    def neighbors(self, x: int, theta: int) -> List[Tuple[int,int]]:
        s = self.step
        nbrs = [(x + s, theta), (x - s, theta)]
        if self.theta_enabled:
            t = self.theta_step
            nbrs.extend([(x, theta + t), (x, theta - t)])
        return nbrs

    def start_state(self) -> Tuple[int,int]:
        return (self.x_B, self.theta_B if self.theta_enabled else 0)

    def goal_state(self) -> Tuple[int,int]:
        return (self.x_A, self.theta_A if self.theta_enabled else 0)


def exhaustive_radial(problem: AlignmentProblem, verbose: bool=False) -> RunResult:
    """
    Exhaustiva radial:
      - 1D: +k, -k desde B
      - 2D (con θ): capas de distancia Manhattan k alrededor de (x_B, θ_B) en orden "diamante"
    """
    x_B, theta_B = problem.start_state()
    x_A, theta_A = problem.goal_state()
    s = problem.step
    ts = problem.theta_step
    visited: Set[Tuple[int,int]] = set()
    order: List[Tuple[int,int]] = []

    def verify(x: int, th: int) -> bool:
        order.append((x, th))
        if verbose:
            msg = f"Verifico (x={x}"
            if problem.theta_enabled:
                msg += f", θ={th}) → "
            else:
                msg += ") → "
            msg += ("META" if problem.is_goal(x, th) else "no meta")
            print(f"[Exhaustiva] {msg}")
        return problem.is_goal(x, th)

    # k = 0 (origen)
    visited.add((x_B, theta_B))
    if verify(x_B, theta_B):
        return RunResult("exhaustive", problem.x_B, problem.x_A, s, order, len(order),
                         problem.theta_enabled, problem.theta_B, problem.theta_A, ts)

    k = 1
    while True:
        if not problem.theta_enabled:
            # 1D: derecha / izquierda
            cand = [(x_B + k*s, theta_B), (x_B - k*s, theta_B)]
        else:
            # 2D: capa L1 = k alrededor de (x_B, theta_B)
            cand = []
            for dx in range(-k, k+1):
                rem = k - abs(dx)
                # rem >= 0; th limites pueden ser +rem y -rem (evita duplicacion si rem==0)
                if rem == 0:
                    cand.append((x_B + dx*s, theta_B))
                else:
                    cand.append((x_B + dx*s, theta_B + rem*ts))
                    cand.append((x_B + dx*s, theta_B - rem*ts))

        # Verificar candidatos de la capa k
        for (xx, th) in cand:
            if (xx, th) in visited:
                continue
            visited.add((xx, th))
            if verify(xx, th):
                return RunResult("exhaustive", problem.x_B, problem.x_A, s, order, len(order),
                                 problem.theta_enabled, problem.theta_B, problem.theta_A, ts)
        k += 1


def astar(problem: AlignmentProblem, verbose: bool=False) -> RunResult:
    #A*: f = g + h, con h = distancia Manhattan a la meta.
    
    x_B, theta_B = problem.start_state()
    x_A, theta_A = problem.goal_state()
    s = problem.step
    ts = problem.theta_step

    if verbose:
        if problem.theta_enabled:
            print(f"[A*] Inicio en B=(x={x_B}, θ={theta_B}), meta A=(x={x_A}, θ={theta_A}), pasos ΔH={s}, Δθ={ts}")
        else:
            print(f"[A*] Inicio en B=x={x_B}, meta A=x={x_A}, paso ΔH={s}")

    openpq: List[PQNode] = []
    start = PQNode(f=problem.heuristic(x_B, theta_B), x=x_B, theta=theta_B, g=0, parent=None)
    heapq.heappush(openpq, start)
    g_cost: Dict[Tuple[int,int], int] = {(x_B, theta_B): 0}
    closed: Set[Tuple[int,int]] = set()
    visit_order: List[Tuple[int,int]] = []

    while openpq:
        node = heapq.heappop(openpq)
        x, th = node.x, node.theta
        if (x, th) in closed:
            continue
        closed.add((x, th))
        visit_order.append((x, th))
        if verbose:
            print(f"[A*] Expando (x={x}, θ={th}) (g={node.g}, h={problem.heuristic(x, th)}, f={node.f})")
        if problem.is_goal(x, th):
            return RunResult("astar", problem.x_B, problem.x_A, s, visit_order, len(visit_order),
                             problem.theta_enabled, problem.theta_B, problem.theta_A, ts)

        for (nx, nth) in problem.neighbors(x, th):
            tentative_g = g_cost[(x, th)] + 1  # coste unitario por movimiento
            if (nx, nth) in g_cost and tentative_g >= g_cost[(nx, nth)]:
                continue
            g_cost[(nx, nth)] = tentative_g
            f = tentative_g + problem.heuristic(nx, nth)
            heapq.heappush(openpq, PQNode(f=f, x=nx, theta=nth, g=tentative_g, parent=None))

    raise RuntimeError("A* no encontró solución (no debería pasar con movimientos libres)")


def greedy(problem: AlignmentProblem, verbose: bool=False) -> RunResult:
    
    #Greedy Best-First: ordena por h(x, θ), ignora g.
    
    x_B, theta_B = problem.start_state()
    x_A, theta_A = problem.goal_state()
    s = problem.step
    ts = problem.theta_step

    if verbose:
        if problem.theta_enabled:
            print(f"[Greedy] Inicio en B=(x={x_B}, θ={theta_B}), meta A=(x={x_A}, θ={theta_A}), pasos ΔH={s}, Δθ={ts}")
        else:
            print(f"[Greedy] Inicio en B=x={x_B}, meta A=x={x_A}, paso ΔH={s}")

    frontier: List[Tuple[int, Tuple[int,int]]] = []  # (h, (x,theta))
    heapq.heappush(frontier, (problem.heuristic(x_B, theta_B), (x_B, theta_B)))
    visited: Set[Tuple[int,int]] = set()
    order: List[Tuple[int,int]] = []

    while frontier:
        h, (x, th) = heapq.heappop(frontier)
        if (x, th) in visited:
            continue
        visited.add((x, th))
        order.append((x, th))
        if verbose:
            print(f"[Greedy] Verifico (x={x}, θ={th}) (h={h}) → {'META' if problem.is_goal(x, th) else 'no meta'}")
        if problem.is_goal(x, th):
            return RunResult("greedy", problem.x_B, problem.x_A, s, order, len(order),
                             problem.theta_enabled, problem.theta_B, problem.theta_A, ts)
        for (nx, nth) in problem.neighbors(x, th):
            if (nx, nth) not in visited:
                heapq.heappush(frontier, (problem.heuristic(nx, nth), (nx, nth)))

    raise RuntimeError("Greedy no encontró solución (no debería pasar con movimientos libres)")


def print_summary(res: RunResult):
    print("\n=== Resumen ===")
    print(f"Método: {res.method}")
    if res.theta_enabled:
        print(f"B = (x={res.x_B}, θ={res.theta_B}) | A = (x={res.x_A}, θ={res.theta_A}) | ΔH = {res.step} | Δθ = {res.theta_step}")
    else:
        print(f"B = x={res.x_B} | A = x={res.x_A} | ΔH = {res.step}")
    print(f"Camino / orden de intentos: {res.path}")
    print(f"N (intentos) = {res.N}")
    print(f"d (distancia mínima en pasos) = {res.d}")
    print(f"P = d/N = {res.P:.3f}")


def plot_path(res: RunResult, title_suffix: Optional[str]=None):
    if not HAVE_MPL:
        print("(matplotlib no disponible; no se graficará)")
        return
    import matplotlib.pyplot as plt
    xs = [p[0] for p in res.path]
    idx = list(range(len(xs)))
    plt.figure()
    ttl = f"Trayectoria x ({res.method})"
    if title_suffix:
        ttl += " - " + title_suffix
    plt.plot(idx, xs, marker='o')
    plt.xlabel("Índice de intentos")
    plt.ylabel("Posición x")
    plt.title(ttl)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    if res.theta_enabled:
        thetas = [p[1] for p in res.path]
        plt.figure()
        ttl = f"Trayectoria θ ({res.method})"
        if title_suffix:
            ttl += " - " + title_suffix
        plt.plot(idx, thetas, marker='o')
        plt.xlabel("Índice de intentos")
        plt.ylabel("Ángulo θ (grados)")
        plt.title(ttl)
        plt.grid(True)
        plt.tight_layout()
        plt.show()


def run_instance(method: str, seed: Optional[int], max_offset: int, step: int,
                 verbose: bool, do_plot: bool,
                 theta_enabled: bool, theta_step: int, theta_maxdeg: int) -> RunResult:
    rng = random.Random(seed)
    x_B = 0
    theta_B = 0

    # Sample random A (x and optional theta)
    def sample_xA():
        # A no coincide con B para que haya movimiento; elegimos en [-max_offset, max_offset] pero ≠ 0 y múltiplo de step
        while True:
            k = rng.randint(-max_offset, max_offset)
            if k != 0 and k % step == 0:
                return k

    def sample_thetaA():
        if not theta_enabled:
            return 0
        # Elegimos en [-theta_maxdeg, theta_maxdeg], múltiplo de theta_step; puede salir 0
        k = rng.randint(-theta_maxdeg, theta_maxdeg)
        # Snap to multiple of theta_step
        k = int(round(k / theta_step)) * theta_step
        return k

    x_A = sample_xA()
    theta_A = sample_thetaA()

    prob = AlignmentProblem(x_B=x_B, x_A=x_A, step=step,
                            theta_enabled=theta_enabled,
                            theta_B=theta_B, theta_A=theta_A, theta_step=theta_step)

    if method == "exhaustive":
        res = exhaustive_radial(prob, verbose=verbose)
    elif method == "astar":
        res = astar(prob, verbose=verbose)
    elif method == "greedy":
        res = greedy(prob, verbose=verbose)
    else:
        raise ValueError("Método desconocido")

    print_summary(res)
    if do_plot:
        suffix = f"B→A (seed={seed})"
        plot_path(res, title_suffix=suffix)
    return res


def run_compare(seed: Optional[int], max_offset: int, step: int, verbose: bool, do_plot: bool,
                theta_enabled: bool, theta_step: int, theta_maxdeg: int):
    rng = random.Random(seed)
    x_B = 0
    theta_B = 0

    def sample_xA():
        while True:
            k = rng.randint(-max_offset, max_offset)
            if k != 0 and k % step == 0:
                return k

    def sample_thetaA():
        if not theta_enabled:
            return 0
        k = rng.randint(-theta_maxdeg, theta_maxdeg)
        k = int(round(k / theta_step)) * theta_step
        return k

    x_A = sample_xA()
    theta_A = sample_thetaA()

    prob = AlignmentProblem(x_B=x_B, x_A=x_A, step=step,
                            theta_enabled=theta_enabled,
                            theta_B=theta_B, theta_A=theta_A, theta_step=theta_step)

    if verbose:
        if theta_enabled:
            print(f"[COMPARE] Instancia fija: B=(x={prob.x_B}, θ={prob.theta_B})  A=(x={prob.x_A}, θ={prob.theta_A})  ΔH={prob.step}  Δθ={prob.theta_step}")
        else:
            print(f"[COMPARE] Instancia fija: B=x={prob.x_B}  A=x={prob.x_A}  ΔH={prob.step}")

    res_exh = exhaustive_radial(prob, verbose=verbose)
    res_ast = astar(prob, verbose=verbose)
    res_grd = greedy(prob, verbose=verbose)

    print_summary(res_exh)
    print_summary(res_ast)
    print_summary(res_grd)

    if do_plot and HAVE_MPL:
        import matplotlib.pyplot as plt
        # Gráfico de X
        plt.figure()
        for res in [res_exh, res_ast, res_grd]:
            idx = list(range(len(res.path)))
            xs = [p[0] for p in res.path]
            plt.plot(idx, xs, marker='o', label=res.method)
        plt.xlabel("Índice de intentos")
        plt.ylabel("Posición x")
        ttl = f"Comparación de trayectorias en x (seed={seed})"
        plt.title(ttl)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

        # Si hay rotación, gráfico de θ (figura separada, sin subplots)
        if theta_enabled:
            plt.figure()
            for res in [res_exh, res_ast, res_grd]:
                idx = list(range(len(res.path)))
                ths = [p[1] for p in res.path]
                plt.plot(idx, ths, marker='o', label=res.method)
            plt.xlabel("Índice de intentos")
            plt.ylabel("Ángulo θ (grados)")
            ttl = f"Comparación de trayectorias en θ (seed={seed})"
            plt.title(ttl)
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()
    elif do_plot and not HAVE_MPL:
        print("(matplotlib no disponible; no se graficará)")


def main():
    parser = argparse.ArgumentParser(description="TP2 - Búsqueda exhaustiva y heurística")
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--method", choices=["exhaustive", "astar", "greedy"], help="Método a ejecutar en una instancia aleatoria")
    g.add_argument("--compare", action="store_true", help="Comparar métodos en la misma instancia aleatoria")

    parser.add_argument("--seed", type=int, default=None, help="Semilla aleatoria (para reproducibilidad)")
    parser.add_argument("--max-offset", type=int, default=20, help="Máximo corrimiento |A-B| en unidades de paso (x)")
    parser.add_argument("--step", type=int, default=1, help="Tamaño de paso ΔH (entero positivo)")
    parser.add_argument("--verbose", action="store_true", help="Mostrar pasos durante la búsqueda")
    parser.add_argument("--plot", action="store_true", help="Graficar trayectoria(s) (si matplotlib está disponible)")

    # Rotación (θ) opcional
    parser.add_argument("--enable-rotation", action="store_true", help="Habilitar dimensión de rotación θ")
    parser.add_argument("--theta-step", type=int, default=1, help="Paso angular Δθ (grados, entero positivo)")
    parser.add_argument("--theta-maxdeg", type=int, default=10, help="Máximo |θ_A| al samplear la meta (grados)")

    args = parser.parse_args()

    if args.method:
        run_instance(args.method, args.seed, args.max_offset, args.step, args.verbose, args.plot,
                     args.enable_rotation, args.theta_step, args.theta_maxdeg)
    else:
        run_compare(args.seed, args.max_offset, args.step, args.verbose, args.plot,
                    args.enable_rotation, args.theta_step, args.theta_maxdeg)


if __name__ == "__main__":
    main()
