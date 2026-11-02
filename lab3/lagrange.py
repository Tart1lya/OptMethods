#!/usr/bin/env python3
import sympy as sp


def main():
    x, y, lam = sp.symbols("x y lam", real=True)
    f = x*y
    h = x**2 + y**2 - 8
    L = f + lam*h
    unknowns = (x, y, lam)
    equations = [sp.diff(L, v) for v in unknowns]

    print("Вариант 2")
    print(f"f = {f}")
    print(f"h = {h} = 0")
    print(f"grad h = {[sp.diff(h, v) for v in (x, y)]}")
    print(f"L = {L}")
    print("Система:")
    for v, equation in zip(unknowns, equations):
        print(f"dL/d{v} = {equation} = 0")

    solutions = sp.solve(equations, unknowns, dict=True)
    print(f"Число действительных решений: {len(solutions)}")
    for i, solution in enumerate(solutions, 1):
        coordinates = ", ".join(f"{v} = {solution[v]}" for v in unknowns)
        residual = sp.simplify(h.subs(solution))
        value = sp.simplify(f.subs(solution))
        print(f"{i}) {coordinates}; h = {residual}; f = {value}")


if __name__ == "__main__":
    main()
