import dearpygui.dearpygui as dpg
import math
import random


def mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    result = [[0]*m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    return result

def transpose(A):
    return [list(row) for row in zip(*A)]

def dot(v1, v2):
    return sum(a*b for a,b in zip(v1,v2))

def norm(v):
    return math.sqrt(dot(v,v))


def eigenvalues_2x2(A):
    a,b = A[0]
    c,d = A[1]
    trace = a + d
    det = a*d - b*c
    disc = trace**2 - 4*det
    if disc < 0:
        raise ValueError("Complex eigenvalues not supported")
    root = math.sqrt(disc)
    return [(trace+root)/2, (trace-root)/2]

def eigenvector_2x2(A, lam):
    a,b = A[0]
    c,d = A[1]
    m1 = a - lam
    m2 = d - lam
    if abs(b) > 1e-8:
        return [1, -m1/b]
    elif abs(c) > 1e-8:
        return [-m2/c, 1]
    return [1,0]


def qr_decomposition(A):
    n = len(A)
    Q = [[0]*n for _ in range(n)]
    R = [[0]*n for _ in range(n)]
    Acols = transpose(A)

    for i in range(n):
        v = Acols[i][:]
        for j in range(i):
            qj = [Q[row][j] for row in range(n)]
            R[j][i] = dot(qj, Acols[i])
            v = [v[k] - R[j][i]*qj[k] for k in range(n)]
        R[i][i] = norm(v)
        if R[i][i] == 0:
            continue
        q = [x/R[i][i] for x in v]
        for row in range(n):
            Q[row][i] = q[row]
    return Q, R

def eigenvalues_qr(A, iterations=60):
    Ak = [row[:] for row in A]
    for _ in range(iterations):
        Q,R = qr_decomposition(Ak)
        Ak = mat_mul(R,Q)
    return [round(Ak[i][i],6) for i in range(len(Ak))]


def power_iteration(A, steps=100):
    n = len(A)
    b = [random.random() for _ in range(n)]
    for _ in range(steps):
        b_new = [sum(A[i][j]*b[j] for j in range(n)) for i in range(n)]
        mag = norm(b_new)
        if mag == 0:
            break
        b = [x/mag for x in b_new]
    return [round(x,6) for x in b]


def compute_eigen(matrix):
    n = len(matrix)

    if n == 2:
        vals = eigenvalues_2x2(matrix)
        vecs = [eigenvector_2x2(matrix,v) for v in vals]
        return vals, vecs

    elif n == 3:
        vals = eigenvalues_qr(matrix)
        vec = power_iteration(matrix)
        return vals, [vec]

    else:
        vals = eigenvalues_qr(matrix)
        vec = power_iteration(matrix)
        return vals, [vec]


def calculate_callback():
    try:
        raw = dpg.get_value("matrix_input").strip()
        rows = raw.split("\n")
        matrix = [list(map(float, r.split())) for r in rows]

        n = len(matrix)
        for row in matrix:
            if len(row) != n:
                raise ValueError("Matrix must be square")

        vals, vecs = compute_eigen(matrix)

        result = f"Eigenvalues:\n{vals}\n\nEigenvectors:\n"
        for v in vecs:
            result += str(v) + "\n"

        dpg.set_value("output_text", result)

    except Exception as e:
        dpg.set_value("output_text", f"Error:\n{e}")



dpg.create_context()

with dpg.window(label="Eigen Calculator (Pure Math)", width=600, height=500):
    dpg.add_text("Enter square matrix (space separated, newline per row)")
    dpg.add_input_text(
        tag="matrix_input",
        multiline=True,
        width=550,
        height=150,
        default_value="1 2\n5 4"
    )

    dpg.add_button(label="Calculate", callback=calculate_callback)
    dpg.add_spacer(height=10)

    dpg.add_input_text(
        tag="output_text",
        multiline=True,
        width=550,
        height=220,
        readonly=True
    )

dpg.create_viewport(title="Eigen App", width=640, height=520)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()
