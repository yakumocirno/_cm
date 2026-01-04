# -*- coding: utf-8 -*-
"""
線性代數實作大全（不靠 numpy.linalg 分解）
1) 遞回(拉普拉斯展開)算 det
2) LU 分解（含部分樞紐 pivot）算 det
3) 驗證：LU / 特徵值分解 / SVD 分解 乘回原矩陣
4) 用特徵值分解做出 SVD（透過 A^T A 的特徵分解）
5) PCA（用 SVD 做主成分分析）

允許使用 numpy 做矩陣乘法/基本運算，但分解與核心算法自己寫。
"""

import math
import numpy as np

EPS = 1e-10

# ============================================================
# 小工具
# ============================================================

def is_close_matrix(A, B, tol=1e-8):
    return np.linalg.norm(A - B, ord='fro') <= tol

def normalize(v):
    n = np.linalg.norm(v)
    if n < EPS:
        return v
    return v / n

def eye(n):
    return np.eye(n, dtype=float)

# ============================================================
# 1) 遞回方式計算行列式：Laplace 展開（O(n!)，只適合小矩陣）
# ============================================================

def det_recursive(A):
    A = np.array(A, dtype=float)
    n = A.shape[0]
    assert A.shape[0] == A.shape[1], "det needs square matrix"

    if n == 1:
        return A[0, 0]
    if n == 2:
        return A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0]

    # 沿第一列展開
    total = 0.0
    for j in range(n):
        if abs(A[0, j]) < EPS:
            continue
        M = np.delete(np.delete(A, 0, axis=0), j, axis=1)  # minor
        cofactor = ((-1) ** (0 + j)) * A[0, j] * det_recursive(M)
        total += cofactor
    return total

# ============================================================
# 2) LU 分解（部分樞紐 pivoting）：PA = LU
#    用 LU 算 det：det(A) = det(P)^(-1) * det(L)*det(U)
#    L 對角 = 1 => det(L)=1；det(U)=U 對角乘積
#    det(P) = (-1)^(交換次數)
# ============================================================

def lu_decompose_partial_pivot(A):
    A = np.array(A, dtype=float)
    n = A.shape[0]
    assert A.shape[0] == A.shape[1]

    U = A.copy()
    L = eye(n)
    P = eye(n)
    swap_count = 0

    for k in range(n):
        # 找 pivot
        pivot = k + np.argmax(np.abs(U[k:, k]))
        if abs(U[pivot, k]) < EPS:
            # 奇異矩陣
            continue

        if pivot != k:
            # 交換 U 的列
            U[[k, pivot], :] = U[[pivot, k], :]
            # 交換 P 的列
            P[[k, pivot], :] = P[[pivot, k], :]
            # L 的前 k 欄也要交換（已經算過的部分）
            if k > 0:
                L[[k, pivot], :k] = L[[pivot, k], :k]
            swap_count += 1

        # 消去
        for i in range(k+1, n):
            if abs(U[k, k]) < EPS:
                continue
            factor = U[i, k] / U[k, k]
            L[i, k] = factor
            U[i, k:] -= factor * U[k, k:]

    return P, L, U, swap_count

def det_via_lu(A):
    P, L, U, swaps = lu_decompose_partial_pivot(A)
    detU = float(np.prod(np.diag(U)))
    detP = -1.0 if (swaps % 2 == 1) else 1.0
    # PA = LU  => det(P)*det(A) = det(L)*det(U)=det(U)
    # det(A) = det(U)/det(P)
    return detU / detP

# ============================================================
# 3) 特徵值分解：對稱矩陣用 Jacobi 旋轉法  A = Q Λ Q^T
#    （為了 PCA / A^T A 的特徵分解做 SVD）
# ============================================================

def jacobi_eigen_symmetric(A, max_iter=10000, tol=1e-12):
    A = np.array(A, dtype=float)
    n = A.shape[0]
    assert n == A.shape[1], "square only"
    # 必須近似對稱
    if np.linalg.norm(A - A.T, ord='fro') > 1e-8:
        raise ValueError("Jacobi needs symmetric matrix")

    Q = eye(n)
    S = A.copy()

    def offdiag_norm(M):
        return math.sqrt(np.sum((M - np.diag(np.diag(M)))**2))

    for _ in range(max_iter):
        if offdiag_norm(S) < tol:
            break

        # 找最大非對角元素
        i, j = 0, 1
        maxv = 0.0
        for r in range(n):
            for c in range(r+1, n):
                v = abs(S[r, c])
                if v > maxv:
                    maxv = v
                    i, j = r, c

        if maxv < tol:
            break

        # Jacobi 旋轉：消掉 S[i,j]
        if abs(S[i, i] - S[j, j]) < EPS:
            theta = math.pi / 4
        else:
            theta = 0.5 * math.atan2(2*S[i, j], (S[j, j] - S[i, i]))

        c = math.cos(theta)
        s = math.sin(theta)

        # 更新 S = G^T S G, Q = Q G
        G = eye(n)
        G[i, i] = c
        G[j, j] = c
        G[i, j] = s
        G[j, i] = -s

        S = G.T @ S @ G
        Q = Q @ G

    eigenvalues = np.diag(S).copy()
    eigenvectors = Q.copy()
    return eigenvalues, eigenvectors

def sort_eigs_desc(vals, vecs):
    idx = np.argsort(vals)[::-1]
    return vals[idx], vecs[:, idx]

# ============================================================
# 4) 自己做 SVD（用特徵值分解）：A^T A = V Σ^2 V^T
#    Σ = sqrt(特徵值)，U = A V Σ^{-1}
# ============================================================

def svd_from_eig(A, tol=1e-10):
    A = np.array(A, dtype=float)
    m, n = A.shape

    AtA = A.T @ A  # n x n 對稱
    vals, V = jacobi_eigen_symmetric(AtA, max_iter=20000, tol=1e-14)
    vals, V = sort_eigs_desc(vals, V)

    # 奇異值
    sig = np.sqrt(np.clip(vals, 0.0, None))

    # 建 Σ（m x n 常見；這裡回傳 full 版本：U(mxm), Σ(mxn), V(nxn)
    U = eye(m)
    Sigma = np.zeros((m, n), dtype=float)

    # 計算前 r 個非零奇異值的 U
    r = min(m, n)
    Ucols = []
    for k in range(r):
        if sig[k] > tol:
            uk = A @ V[:, k] / sig[k]
            uk = normalize(uk)
            Ucols.append(uk)
            Sigma[k, k] = sig[k]
        else:
            Sigma[k, k] = 0.0

    # 把已算出的 uk 填入 U 的前幾欄
    if len(Ucols) > 0:
        U[:, :len(Ucols)] = np.column_stack(Ucols)

    # 若需要補齊 U 的其餘正交基底（這裡用簡單 Gram-Schmidt）
    # 讓 U 近似正交
    for k in range(len(Ucols), m):
        v = np.random.randn(m)
        for j in range(k):
            v -= (U[:, j] @ v) * U[:, j]
        U[:, k] = normalize(v)

    return U, Sigma, V

# ============================================================
# 5) 驗證：LU / 特徵值分解 / SVD 分解可以乘回原矩陣
# ============================================================

def verify_lu(A):
    P, L, U, _ = lu_decompose_partial_pivot(A)
    left = P @ A
    right = L @ U
    return is_close_matrix(left, right, tol=1e-8), P, L, U

def verify_eig_symmetric(A):
    # 只驗證對稱矩陣
    vals, Q = jacobi_eigen_symmetric(A)
    Lam = np.diag(vals)
    recon = Q @ Lam @ Q.T
    return is_close_matrix(A, recon, tol=1e-7), vals, Q

def verify_svd(A):
    U, S, V = svd_from_eig(A)
    recon = U @ S @ V.T
    return is_close_matrix(A, recon, tol=1e-6), U, S, V

# ============================================================
# 6) PCA：用 SVD 做主成分分析
#    X: (samples x features)
#    steps:
#      - 中心化 Xc = X - mean
#      - SVD: Xc = U Σ V^T
#      - 主成分方向：V 的欄向量
#      - 變異量比例：Σ^2 / sum(Σ^2)
#      - 投影到前 k 維：Z = Xc @ V_k
# ============================================================

def pca_via_svd(X, k=2):
    X = np.array(X, dtype=float)
    n, d = X.shape
    mu = X.mean(axis=0, keepdims=True)
    Xc = X - mu

    U, S, V = svd_from_eig(Xc)  # 用我們自己的 SVD
    # S 是 m x n 形狀，但對角線才有值
    sing = np.diag(S)[:min(n, d)]

    # explained variance
    var = sing**2
    total = var.sum() if var.sum() > EPS else 1.0
    ratio = var / total

    V_k = V[:, :k]
    Z = Xc @ V_k  # 投影後座標
    return {
        "mean": mu.ravel(),
        "components": V_k,          # (d x k)
        "explained_variance_ratio": ratio[:k],
        "projected": Z              # (n x k)
    }

# ============================================================
# Demo
# ============================================================

if __name__ == "__main__":

    # --- 測試矩陣（一般矩陣，用來測 LU / SVD）---
    A = np.array([
        [2,  1,  3],
        [4,  1,  6],
        [1, -2,  0],
    ], dtype=float)

    print("A=\n", A)

    # 1) det recursive
    print("\n[1] det_recursive(A) =", det_recursive(A))

    # 2) det via LU
    print("[2] det_via_lu(A)     =", det_via_lu(A))

    # 3) verify LU
    ok_lu, P, L, U = verify_lu(A)
    print("\n[3] Verify LU:  P@A == L@U ? ", ok_lu)
    print("P=\n", P)
    print("L=\n", L)
    print("U=\n", U)

    # 4) verify SVD (從特徵值做 SVD)
    ok_svd, U2, S2, V2 = verify_svd(A)
    print("\n[4] Verify SVD (from eig):  A == U*S*V^T ? ", ok_svd)
    print("U=\n", U2)
    print("S=\n", S2)
    print("V=\n", V2)

    # --- 對稱矩陣（用來測特徵值分解）---
    B = np.array([
        [4, 1, 1],
        [1, 3, 0],
        [1, 0, 2],
    ], dtype=float)

    ok_eig, vals, Q = verify_eig_symmetric(B)
    print("\n[5] Verify Eigen (symmetric Jacobi): B == Q*Λ*Q^T ? ", ok_eig)
    print("eigenvalues =", vals)
    print("Q=\n", Q)

    # 6) PCA demo
    # 10 筆 3 維資料（samples x features）
    X = np.array([
        [2.5, 2.4, 0.5],
        [0.5, 0.7, 1.2],
        [2.2, 2.9, 0.3],
        [1.9, 2.2, 0.7],
        [3.1, 3.0, 0.2],
        [2.3, 2.7, 0.4],
        [2.0, 1.6, 0.9],
        [1.0, 1.1, 1.5],
        [1.5, 1.6, 1.0],
        [1.1, 0.9, 1.8],
    ], dtype=float)

    res = pca_via_svd(X, k=2)
    print("\n[6] PCA via SVD")
    print("mean =", res["mean"])
    print("components (d x k)=\n", res["components"])
    print("explained_variance_ratio =", res["explained_variance_ratio"])
    print("projected (n x k)=\n", res["projected"])
