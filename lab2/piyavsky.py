#!/usr/bin/env python3
import heapq
import math
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

GRID = 200001
STAGES = (4, 10, 30)

ENV = {"sin": np.sin, "cos": np.cos, "tan": np.tan, "atan": np.arctan,
       "sinh": np.sinh, "cosh": np.cosh, "tanh": np.tanh, "exp": np.exp,
       "log": np.log, "sqrt": np.sqrt, "abs": np.abs, "pi": math.pi, "e": math.e}


def parse_function(text):
    expr = text.split("=", 1)[1] if "=" in text else text
    code = compile(expr.strip().replace("^", "**"), "<f>", "eval")
    unknown = [n for n in code.co_names if n not in ENV and n != "x"]
    if unknown:
        raise ValueError("неизвестные имена в формуле: " + ", ".join(unknown))
    env = {"__builtins__": {}, **ENV}

    def f(x):
        with np.errstate(all="ignore"):
            return eval(code, env, {"x": x})
    return f


def scalar(f, t):
    y = float(f(np.float64(t)))
    if not math.isfinite(y):
        raise ValueError(f"функция не определена в точке x = {t:g}")
    return y


def estimate_lipschitz(f, a, b):
    x = np.linspace(a, b, GRID)
    y = np.broadcast_to(np.asarray(f(x), dtype=float), x.shape)
    if not np.all(np.isfinite(y)):
        raise ValueError("функция не определена или не ограничена на отрезке")
    slope = float(np.max(np.abs(np.diff(y) / np.diff(x))))
    return slope if slope > 0 else 1.0


def push(heap, L, x1, y1, x2, y2):
    R = (y1 + y2) / 2 - L * (x2 - x1) / 2
    v = (x1 + x2) / 2 + (y1 - y2) / (2 * L)
    heapq.heappush(heap, (R, v, x1, y1, x2, y2))


def solve(f, a, b, eps, L):
    pts = [(a, f(a)), (b, f(b))]
    heap = []
    push(heap, L, *pts[0], *pts[1])
    x_best, f_best = min(pts, key=lambda p: p[1])
    while f_best - heap[0][0] > eps:
        _, v, x1, y1, x2, y2 = heapq.heappop(heap)
        y = f(v)
        pts.append((v, y))
        if y < f_best:
            x_best, f_best = v, y
        push(heap, L, x1, y1, v, y)
        push(heap, L, v, y, x2, y2)
    return pts, x_best, f_best, heap[0][0]


def broken_line(pts, L):
    p = sorted(pts)
    X, Y = [p[0][0]], [p[0][1]]
    for (x1, y1), (x2, y2) in zip(p, p[1:]):
        X += [(x1 + x2) / 2 + (y1 - y2) / (2 * L), x2]
        Y += [(y1 + y2) / 2 - L * (x2 - x1) / 2, y2]
    return X, Y


def draw(ax, xs, ys, pts, L, ylim):
    X, Y = broken_line(pts, L)
    ax.plot(xs, ys, color="black", lw=1.3, label="$f(x)$")
    ax.plot(X, Y, color="tab:orange", lw=0.9, label="ломаная $F_n(x)$")
    ax.scatter(*zip(*pts), s=9, color="tab:blue", zorder=3, label="испытания")
    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(*ylim)
    ax.grid(alpha=0.3)


def plot(f, a, b, pts, L, x_best, f_best, text, summary):
    xs = np.linspace(a, b, 4000)
    ys = np.broadcast_to(np.asarray(f(xs), dtype=float), xs.shape)
    pad = 0.1 * (ys.max() - ys.min() or 1.0)
    ylim = (ys.min() - pad, ys.max() + pad)
    stages = [n for n in STAGES if n < len(pts)]
    fig = plt.figure(figsize=(8, 5.2 if stages else 3.6))
    gs = fig.add_gridspec(2 if stages else 1, max(len(stages), 1))
    top = fig.add_subplot(gs[0, :])
    draw(top, xs, ys, pts, L, ylim)
    top.plot([x_best], [f_best], marker="*", ms=16, color="red", ls="none",
             zorder=4, label="найденный минимум")
    top.axvline(x_best, color="red", ls="--", lw=0.8)
    top.set_title(summary, fontsize=10)
    top.set_xlabel("$x$")
    top.legend(loc="upper center", ncol=4, fontsize=9)
    for i, n in enumerate(stages):
        ax = fig.add_subplot(gs[1, i])
        draw(ax, xs, ys, pts[:n], L, ylim)
        ax.set_title(f"ломаная после $n = {n}$ испытаний", fontsize=10)
        ax.set_xlabel("$x$")
    fig.suptitle(text, fontsize=12)
    fig.tight_layout()
    plt.show()


def to_float(s):
    return float(s.replace(",", "."))


def read_input():
    text = input("Функция f(x): ")
    a = to_float(input("Левый конец отрезка a: "))
    b = to_float(input("Правый конец отрезка b: "))
    eps = to_float(input("Точность eps: "))
    return text, a, b, eps


def main():
    text, a, b, eps = read_input()
    try:
        if not a < b or eps <= 0:
            raise ValueError("нужно a < b и eps > 0")
        f = parse_function(text)
        start = time.perf_counter()
        L = estimate_lipschitz(f, a, b)
        pts, x_best, f_best, low = solve(lambda t: scalar(f, t), a, b, eps, L)
        elapsed = time.perf_counter() - start
    except (ValueError, TypeError, SyntaxError) as err:
        sys.exit(f"Ошибка: {err}")

    n = len(pts)
    print(f"Функция:  {text}")
    print(f"Отрезок:  [{a:g}, {b:g}],  eps = {eps:g}")
    print(f"Константа Липшица (оценка):  L = {L:.6f}")
    print(f"Аргумент минимума:           x* = {x_best:.6f}")
    print(f"Минимальное значение:        f(x*) = {f_best:.6f}")
    print(f"Нижняя оценка минимума:      {low:.6f}  (f(x*) - оценка = {f_best - low:.6f})")
    print(f"Итераций: {n - 2}  (вычислений функции N = {n}, включая концы отрезка)")
    print(f"Время счёта: {elapsed:.4f} с")

    summary = (f"$x^*\\approx{x_best:.6f}$,  $f(x^*)\\approx{f_best:.6f}$,  "
               f"итераций {n - 2} (N = {n}),  время {elapsed:.3f} с,  $L={L:.3f}$")
    plot(f, a, b, pts, L, x_best, f_best, text, summary)


if __name__ == "__main__":
    main()
