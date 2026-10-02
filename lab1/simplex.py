#!/usr/bin/env python3
import sys
from decimal import Decimal
from fractions import Fraction as Fr

C = [2, 3, 1, 4]
A = [[1, 1, 1, 1],
     [2, 1, -1, 1],
     [0, 1, 2, 1]]
B = [10, 8, 5]
SIGNS = ["<=", "=", ">="]


def num(x, tex=False):
    if x.denominator == 1:
        return str(x.numerator)
    d = x.denominator
    for p in (2, 5):
        while d % p == 0:
            d //= p
    if d == 1:
        s = format(Decimal(x.numerator) / Decimal(x.denominator), "f")
        return s.replace(".", "{,}") if tex else s
    return f"{x.numerator}/{x.denominator}"


def canonical():
    n = len(C)
    rows = [[Fr(v) for v in row] for row in A]
    names = [f"x{j + 1}" for j in range(n)]
    for i, s in enumerate(SIGNS):
        if s == "=":
            continue
        names.append(f"x{len(names) + 1}")
        for row in rows:
            row.append(Fr(0))
        rows[i][-1] = Fr(1) if s == "<=" else Fr(-1)
    assert all(v >= 0 for v in B)
    c = [Fr(v) for v in C] + [Fr(0)] * (len(names) - n)
    return rows, [Fr(v) for v in B], c, names


def unit_basis(rows):
    m = len(rows)
    basis = [None] * m
    for j in range(len(rows[0])):
        col = [rows[i][j] for i in range(m)]
        if sorted(col) == [0] * (m - 1) + [1]:
            i = col.index(1)
            if basis[i] is None:
                basis[i] = j
    return basis


def simplex(rows, rhs, basis, c, names, active, art, tables, title):
    while True:
        D = {j: sum(c[basis[i]] * rows[i][j] for i in range(len(rows))) - c[j]
             for j in active}
        D0 = sum(c[basis[i]] * rhs[i] for i in range(len(rows)))
        t = dict(title=title, cols=sorted(active), basis=list(basis),
                 cB=[c[j] for j in basis], rows=[r[:] for r in rows],
                 rhs=rhs[:], D=D, D0=D0, k=None, r=None, ratios=None)
        tables.append(t)
        k = min(sorted(active), key=lambda j: D[j])
        if D[k] >= 0:
            return rows, rhs, basis, D0
        ratios = {i: rhs[i] / rows[i][k] for i in range(len(rows)) if rows[i][k] > 0}
        if not ratios:
            raise SystemExit("Целевая функция не ограничена сверху - решений нет")
        r = min(ratios, key=lambda i: (ratios[i], i))
        t.update(k=k, r=r, ratios=ratios)
        p = rows[r][k]
        rows[r] = [v / p for v in rows[r]]
        rhs[r] /= p
        for i in range(len(rows)):
            if i != r and rows[i][k] != 0:
                f = rows[i][k]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
                rhs[i] -= f * rhs[r]
        if basis[r] in art:
            active.discard(basis[r])
        basis[r] = k


def text_table(t, names, num_no):
    cols = t["cols"]
    head = ["Базис", "c_B", "b"] + [names[j] for j in cols]
    lines = [[names[t["basis"][i]], num(t["cB"][i]), num(t["rhs"][i])]
             + [num(t["rows"][i][j]) for j in cols] for i in range(len(t["rows"]))]
    lines.append(["Δ", "", num(t["D0"])] + [num(t["D"][j]) for j in cols])
    w = [max(len(r[i]) for r in [head] + lines) for i in range(len(head))]
    out = [f"Таблица {num_no}: {t['title']}"]
    for r in [head] + lines:
        out.append("  ".join(s.rjust(w[i]) for i, s in enumerate(r)))
    if t["k"] is None:
        out.append("Все Δ_j >= 0 - план оптимален.")
    else:
        rat = ", ".join(f"{names[t['basis'][i]]}: {num(v)}" for i, v in t["ratios"].items())
        out.append(f"Ведущий столбец {names[t['k']]} (Δ = {num(t['D'][t['k']])}); "
                   f"отношения b/a: {rat}; ведущая строка - {names[t['basis'][t['r']]]}; "
                   f"ведущий элемент {num(t['rows'][t['r']][t['k']])}")
    return "\n".join(out)


def latex_table(t, names, num_no):
    cols = t["cols"]
    sub = lambda s: s[0] + "_{" + s[1:] + "}"
    ncol = 3 + len(cols)
    out = [f"% Таблица {num_no}: {t['title']}",
           r"\begin{tabular}{|c|c|r|" + "r|" * len(cols) + "}", r"\hline",
           "Базис & $c_{\\text{Б}}$ & $b$ & " + " & ".join(f"${sub(names[j])}$" for j in cols)
           + r" \\ \hline"]
    for i in range(len(t["rows"])):
        cells = []
        for j in cols:
            s = num(t["rows"][i][j], tex=True)
            if t["k"] == j and t["r"] == i:
                s = r"\boxed{" + s + "}"
            cells.append(f"${s}$")
        out.append(f"${sub(names[t['basis'][i]])}$ & ${num(t['cB'][i], tex=True)}$ & "
                   f"${num(t['rhs'][i], tex=True)}$ & " + " & ".join(cells) + r" \\")
    out.append(r"\hline")
    out.append(r"$\Delta$ & & $" + num(t["D0"], tex=True) + "$ & "
               + " & ".join("$" + num(t["D"][j], tex=True) + "$" for j in cols)
               + r" \\ \hline")
    out.append(r"\end{tabular}")
    return "\n".join(out)


def main():
    latex = "--latex" in sys.argv
    rows, rhs, c, names = canonical()
    m, n_can = len(rows), len(names)

    basis = unit_basis(rows)
    art = set()
    for i in range(m):
        if basis[i] is None:
            for row in rows:
                row.append(Fr(0))
            rows[i][-1] = Fr(1)
            names.append(f"y{len(art) + 1}")
            basis[i] = len(names) - 1
            art.add(basis[i])
    c_aux = [Fr(0)] * n_can + [Fr(-1)] * len(art)

    tables = []

    r1, b1, bas1, F = simplex([r[:] for r in rows], rhs[:], basis[:], c_aux, names,
                              set(range(len(names))), art, tables, "вспомогательная задача")
    n_aux = len(tables)
    if F < 0 or any(j in art for j in bas1):
        raise SystemExit("max F < 0 (или y остались в базисе): допустимая область пуста - решений нет")

    rows_m = [r[:n_can] for r in r1]
    x_names = names[:n_can]
    _, rhs_m, bas_m, Z = simplex(rows_m, b1, bas1, c, x_names, set(range(n_can)), set(),
                                 tables, "основная задача")

    for no, t in enumerate(tables, 1):
        nm = names if no <= n_aux else x_names
        print((latex_table if latex else text_table)(t, nm, no))
        print()
    if latex:
        return

    x = [Fr(0)] * n_can
    for i, j in enumerate(bas_m):
        x[j] = rhs_m[i]
    print(f"Ответ:  x* = ({', '.join(num(v) for v in x[:len(C)])}),  Z(x*) = {num(Z)}")
    print("Дополнительные переменные:", ", ".join(f"{names[j]} = {num(x[j])}" for j in range(len(C), n_can)))

    for row, s, bi in zip(A, SIGNS, B):
        lhs = sum(a * v for a, v in zip(row, x))
        ok = {"<=": lhs <= bi, "=": lhs == bi, ">=": lhs >= bi}[s]
        print(f"  {num(lhs)} {s} {bi}: {'верно' if ok else 'НАРУШЕНО'}")
    assert sum(ci * v for ci, v in zip(C, x)) == Z

    try:
        from scipy.optimize import linprog
    except ImportError:
        return
    res = linprog([-v for v in C],
                  A_ub=[A[0], [-v for v in A[2]]], b_ub=[B[0], -B[2]],
                  A_eq=[A[1]], b_eq=[B[1]], bounds=[(0, None)] * len(C))
    print(f"scipy.linprog: x = ({', '.join(f'{v:g}' for v in res.x)}),  Z = {-res.fun:g}")


if __name__ == "__main__":
    main()
