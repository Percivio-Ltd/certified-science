#!/usr/bin/env python3
"""Standalone second-stack replay of the accepted Q4 primitive evidence.

This program deliberately imports no physics_proof_harness module and no prior
Q4 verifier.  It uses only the Python standard library.  The printed R/S
vectors are transcribed below.  Exact algebra uses a generic quotient-basis
implementation of Q(i,sqrt(2),sqrt(3),sqrt(5)); positivity is certified by a
different path from the maintained verifier:

1. enclose every algebraic matrix entry by rational dyadic intervals;
2. form a rational-complex midpoint matrix C and an exact row-sum error rho;
3. prove C-rho*I positive definite by rational-complex LDL.

For a Hermitian exact matrix M, ||M-C||_2 <= rho, so this proves M positive
definite without algebraic-number sign decisions.  The script also rebuilds
the normalized Choi matrices and weighted joint output from the printed
vectors, proves trace preservation exactly, and checks the claimed rational
joint matrix and characteristic polynomial exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from fractions import Fraction as Q
from pathlib import Path
from typing import Any, Iterable


CHECKPOINT_SHA256 = "6af36ce14e3cea5ed8749598574bf6acd1c6eed723008c0f6dcd343a1bb57a60"
FLOOR_CANDIDATE_SHA256 = "980b1af4df2f4984d8b8fce7822a67073cbd3bb3b45e8426c97c08a51a9de419"
JOINT_CLAIM_SHA256 = "a64ad02561a509422756859120ca2be04224ce310266aee352114982588ebf9d"
CAP_CANDIDATE_SHA256 = "5a8a5ec0d91cf6dcdd0ee9805ce25508e36d6ae601341d6904e1b2d07c10e9e0"
FIELD_DEGREES = (2, 3, 5, -1)  # sqrt(2), sqrt(3), sqrt(5), i
ZERO16 = (Q(0),) * 16


class VerificationError(Exception):
    """A fail-closed verification error."""


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_fraction(value: Any, label: str) -> Q:
    if not isinstance(value, str) or not value:
        raise VerificationError(f"{label} is not a rational string")
    try:
        result = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise VerificationError(f"{label} is not rational") from exc
    if str(result) != value:
        raise VerificationError(f"{label} is not canonical")
    return result


class K:
    """Generic 16-term basis for Q(sqrt(2),sqrt(3),sqrt(5),i)."""

    __slots__ = ("c",)

    def __init__(self, coefficients: tuple[Q, ...] = ZERO16):
        if len(coefficients) != 16:
            raise ValueError("field element must have 16 coefficients")
        self.c = coefficients

    @staticmethod
    def rational(value: Q | int) -> "K":
        out = list(ZERO16)
        out[0] = Q(value)
        return K(tuple(out))

    @staticmethod
    def gaussian(real: Q | int, imaginary: Q | int) -> "K":
        out = list(ZERO16)
        out[0] = Q(real)
        out[8] = Q(imaginary)
        return K(tuple(out))

    @staticmethod
    def omega(conjugated: bool = False) -> "K":
        out = list(ZERO16)
        out[0] = Q(-1, 2)
        out[10] = Q(-1 if conjugated else 1, 2)  # i*sqrt(3)/2
        return K(tuple(out))

    def __add__(self, other: "K") -> "K":
        return K(tuple(a + b for a, b in zip(self.c, other.c)))

    def __sub__(self, other: "K") -> "K":
        return K(tuple(a - b for a, b in zip(self.c, other.c)))

    def __neg__(self) -> "K":
        return K(tuple(-a for a in self.c))

    def __mul__(self, other: "K") -> "K":
        out = [Q(0)] * 16
        for left_mask, left in enumerate(self.c):
            if left == 0:
                continue
            for right_mask, right in enumerate(other.c):
                if right == 0:
                    continue
                overlap = left_mask & right_mask
                factor = 1
                for bit, square in enumerate(FIELD_DEGREES):
                    if overlap & (1 << bit):
                        factor *= square
                out[left_mask ^ right_mask] += left * right * factor
        return K(tuple(out))

    def scale(self, value: Q | int) -> "K":
        value = Q(value)
        return K(tuple(value * a for a in self.c))

    def conjugate(self) -> "K":
        return K(tuple(-a if mask & 8 else a for mask, a in enumerate(self.c)))

    def is_zero(self) -> bool:
        return all(value == 0 for value in self.c)

    def is_real(self) -> bool:
        return all(value == 0 for mask, value in enumerate(self.c) if mask & 8)

    def is_rational(self) -> bool:
        return all(value == 0 for value in self.c[1:])

    def as_rational(self) -> Q:
        if not self.is_rational():
            raise VerificationError("expected a rational field element")
        return self.c[0]

    def __eq__(self, other: object) -> bool:
        return isinstance(other, K) and self.c == other.c


K0 = K()
K1 = K.rational(1)


def sqrt_rational(value: Q) -> K:
    """Return the exact positive square root when it lies in the frozen field."""

    if value <= 0:
        raise VerificationError("square root input must be positive")
    combined = value.numerator * value.denominator
    square = 1
    residual = combined
    mask = 0
    for bit, prime in enumerate((2, 3, 5)):
        while residual % (prime * prime) == 0:
            residual //= prime * prime
            square *= prime
        if residual % prime == 0:
            residual //= prime
            mask |= 1 << bit
    root = math.isqrt(residual)
    if root * root != residual:
        raise VerificationError(f"sqrt({value}) is outside the frozen field")
    square *= root
    out = list(ZERO16)
    out[mask] = Q(square, value.denominator)
    result = K(tuple(out))
    if result * result != K.rational(value):
        raise VerificationError("internal exact square-root check failed")
    return result


def inverse_sqrt_rational(value: Q) -> K:
    return sqrt_rational(value).scale(1 / value)


def ksum(values: Iterable[K]) -> K:
    result = K0
    for value in values:
        result = result + value
    return result


def zero_matrix(size: int) -> list[list[K]]:
    return [[K0 for _ in range(size)] for _ in range(size)]


def identity_matrix(size: int) -> list[list[K]]:
    return [[K1 if row == column else K0 for column in range(size)] for row in range(size)]


def matrix_multiply(left: list[list[K]], right: list[list[K]]) -> list[list[K]]:
    rows, shared, columns = len(left), len(right), len(right[0])
    if any(len(row) != shared for row in left):
        raise VerificationError("matrix dimensions differ")
    return [
        [ksum(left[row][k] * right[k][column] for k in range(shared)) for column in range(columns)]
        for row in range(rows)
    ]


R_MATRICES = [
    [[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0]],
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
    [[1, 0, 0, 0], [0, "w", 0, 0], [0, 0, "w2", 0]],
    [[0, 1, 0, 0], [0, 0, "w2", 0], [0, 0, 0, "w"]],
    [[0, 0, 1, 0], [0, 0, 0, -1], [0, 0, 0, 0]],
    [[0, 0, 0, 1], [0, 0, 0, 0], [-1, 0, 0, 0]],
]
S_MATRICES = [
    [[0, 0, 0, 0], [1, 0, 0, 0], [0, -1, 0, 0]],
    [[1, 0, 0, 0], [0, "w2", 0, 0], [0, 0, "w", 0]],
    [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
    [[0, 1, 0, 0], [0, 0, "w", 0], [0, 0, 0, "w2"]],
    [[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]],
    [[0, 0, 0, 1], [0, 0, 0, 0], [1, 0, 0, 0]],
]


def printed_entry(value: Any, conjugated: bool) -> K:
    if value == "w":
        return K.omega(conjugated)
    if value == "w2":
        return K.omega(not conjugated)
    return K.rational(Q(value))


def printed_vectors(matrices: list[list[list[Any]]], conjugated: bool) -> list[list[K]]:
    vectors: list[list[K]] = []
    for matrix in matrices:
        vector = [K0] * 12
        for b in range(3):
            for a in range(4):
                vector[3 * a + b] = printed_entry(matrix[b][a], conjugated)
        vectors.append(vector)
    return vectors


def inner(left: list[K], right: list[K]) -> K:
    return ksum(left[index].conjugate() * right[index] for index in range(len(left)))


def projector(vectors: list[list[K]]) -> tuple[list[list[K]], list[Q]]:
    norms: list[Q] = []
    for row, vector in enumerate(vectors):
        for column in range(row):
            if not inner(vectors[column], vector).is_zero():
                raise VerificationError("printed basis is not orthogonal")
        norm = inner(vector, vector)
        if not norm.is_rational() or norm.as_rational() <= 0:
            raise VerificationError("printed basis norm is not positive rational")
        norms.append(norm.as_rational())
    result = zero_matrix(12)
    for vector, norm in zip(vectors, norms):
        for row in range(12):
            if vector[row].is_zero():
                continue
            for column in range(12):
                if not vector[column].is_zero():
                    result[row][column] = result[row][column] + (
                        vector[row] * vector[column].conjugate()
                    ).scale(1 / norm)
    return result, norms


def partial_trace_b(matrix: list[list[K]]) -> list[list[K]]:
    return [
        [ksum(matrix[3 * a + b][3 * ap + b] for b in range(3)) for ap in range(4)]
        for a in range(4)
    ]


def build_choi(projection: list[list[K]]) -> tuple[list[list[K]], list[Q]]:
    marginal = partial_trace_b(projection)
    weights: list[Q] = []
    for a in range(4):
        for ap in range(4):
            value = marginal[a][ap]
            if a != ap and not value.is_zero():
                raise VerificationError("input marginal is not diagonal")
            if a == ap:
                if not value.is_rational() or value.as_rational() <= 0:
                    raise VerificationError("input marginal diagonal is invalid")
                weights.append(value.as_rational())
    result = zero_matrix(12)
    for a in range(4):
        for b in range(3):
            row = 3 * a + b
            for ap in range(4):
                factor = inverse_sqrt_rational(weights[a] * weights[ap])
                for bp in range(3):
                    column = 3 * ap + bp
                    result[row][column] = projection[row][column] * factor
    traced = partial_trace_b(result)
    for a in range(4):
        for ap in range(4):
            expected = K1 if a == ap else K0
            if traced[a][ap] != expected:
                raise VerificationError("normalized Choi matrix is not exactly trace preserving")
    return result, weights


def partial_transpose_b(matrix: list[list[K]]) -> list[list[K]]:
    result = zero_matrix(12)
    for a in range(4):
        for b in range(3):
            for ap in range(4):
                for bp in range(3):
                    result[3 * a + b][3 * ap + bp] = matrix[3 * a + bp][3 * ap + b]
    return result


def parse_floor_q(candidate: dict[str, Any]) -> dict[str, list[list[K]]]:
    if candidate.get("schema_version") != "renyi-block-positive-candidate-v1":
        raise VerificationError("floor candidate schema differs")
    result: dict[str, list[list[K]]] = {}
    for witness in candidate.get("witnesses", []):
        name = witness.get("name")
        rows = witness.get("Q")
        if name not in ("A", "B") or not isinstance(rows, list) or len(rows) != 12:
            raise VerificationError("floor witness identity or shape differs")
        matrix: list[list[K]] = []
        for row in rows:
            if not isinstance(row, list) or len(row) != 12:
                raise VerificationError("floor witness row shape differs")
            parsed_row = []
            for entry in row:
                try:
                    (nr, dr), (ni, di) = entry
                    parsed_row.append(K.gaussian(Q(int(nr), int(dr)), Q(int(ni), int(di))))
                except Exception as exc:
                    raise VerificationError("floor witness entry differs") from exc
            matrix.append(parsed_row)
        result[name] = matrix
    if set(result) != {"A", "B"}:
        raise VerificationError("floor witness set differs")
    return result


def parse_cap_q(candidate: dict[str, Any]) -> tuple[Q, dict[str, list[list[K]]]]:
    if candidate.get("schema_version") != "q4-output-cap-candidate-v1":
        raise VerificationError("cap candidate schema differs")
    cap = canonical_fraction(candidate.get("cap"), "cap")
    result: dict[str, list[list[K]]] = {}
    for witness in candidate.get("witnesses", []):
        name = witness.get("name")
        if name not in ("A", "B"):
            raise VerificationError("cap witness identity differs")
        matrix = zero_matrix(12)
        seen: set[tuple[int, int]] = set()
        for index, entry in enumerate(witness.get("entries", [])):
            if not isinstance(entry, dict):
                raise VerificationError("cap entry is not an object")
            row, column = entry.get("row"), entry.get("column")
            if not isinstance(row, int) or not isinstance(column, int) or not (0 <= row <= column < 12):
                raise VerificationError("cap entry index differs")
            if (row, column) in seen:
                raise VerificationError("duplicate cap entry")
            seen.add((row, column))
            real = canonical_fraction(entry.get("real"), f"cap real {name}:{index}")
            imaginary = canonical_fraction(entry.get("imaginary"), f"cap imaginary {name}:{index}")
            if row == column and imaginary != 0:
                raise VerificationError("cap diagonal is not real")
            matrix[row][column] = K.gaussian(real, imaginary)
            matrix[column][row] = K.gaussian(real, -imaginary)
        result[name] = matrix
    if set(result) != {"A", "B"}:
        raise VerificationError("cap witness set differs")
    return cap, result


def assert_hermitian(matrix: list[list[K]], label: str) -> None:
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise VerificationError(f"{label} is not square")
    for row in range(size):
        for column in range(row, size):
            if matrix[row][column] != matrix[column][row].conjugate():
                raise VerificationError(f"{label} is not Hermitian")


def sqrt_interval(integer: int, bits: int) -> tuple[Q, Q]:
    denominator = 1 << bits
    target = integer << (2 * bits)
    lower_integer = math.isqrt(target)
    lower = Q(lower_integer, denominator)
    if lower_integer * lower_integer == target:
        return lower, lower
    return lower, Q(lower_integer + 1, denominator)


def scale_interval(interval: tuple[Q, Q], scalar: Q) -> tuple[Q, Q]:
    low, high = interval
    if scalar >= 0:
        return scalar * low, scalar * high
    return scalar * high, scalar * low


def add_interval(left: tuple[Q, Q], right: tuple[Q, Q]) -> tuple[Q, Q]:
    return left[0] + right[0], left[1] + right[1]


def field_rectangle(value: K, bits: int) -> tuple[tuple[Q, Q], tuple[Q, Q]]:
    real = (Q(0), Q(0))
    imaginary = (Q(0), Q(0))
    for mask, coefficient in enumerate(value.c):
        if coefficient == 0:
            continue
        radicand = 1
        for bit, prime in enumerate((2, 3, 5)):
            if mask & (1 << bit):
                radicand *= prime
        term = scale_interval(sqrt_interval(radicand, bits), coefficient)
        if mask & 8:
            imaginary = add_interval(imaginary, term)
        else:
            real = add_interval(real, term)
    return real, imaginary


RC = tuple[Q, Q]


def rc_add(left: RC, right: RC) -> RC:
    return left[0] + right[0], left[1] + right[1]


def rc_sub(left: RC, right: RC) -> RC:
    return left[0] - right[0], left[1] - right[1]


def rc_mul(left: RC, right: RC) -> RC:
    return left[0] * right[0] - left[1] * right[1], left[0] * right[1] + left[1] * right[0]


def rc_conjugate(value: RC) -> RC:
    return value[0], -value[1]


def rc_div_real(value: RC, divisor: Q) -> RC:
    if divisor == 0:
        raise VerificationError("zero rational LDL divisor")
    return value[0] / divisor, value[1] / divisor


def rational_complex_ldl(matrix: list[list[RC]]) -> tuple[bool, list[Q]]:
    size = len(matrix)
    lower = [[(Q(0), Q(0)) for _ in range(size)] for _ in range(size)]
    pivots: list[Q] = []
    for column in range(size):
        pivot = matrix[column][column]
        for k in range(column):
            product = rc_mul(rc_mul(lower[column][k], rc_conjugate(lower[column][k])), (pivots[k], Q(0)))
            pivot = rc_sub(pivot, product)
        if pivot[1] != 0 or pivot[0] <= 0:
            return False, pivots + [pivot[0]]
        pivots.append(pivot[0])
        for row in range(column + 1, size):
            value = matrix[row][column]
            for k in range(column):
                product = rc_mul(rc_mul(lower[row][k], rc_conjugate(lower[column][k])), (pivots[k], Q(0)))
                value = rc_sub(value, product)
            lower[row][column] = rc_div_real(value, pivots[column])
    return True, pivots


def prove_positive_definite(matrix: list[list[K]], label: str, bits: int) -> dict[str, Any]:
    assert_hermitian(matrix, label)
    size = len(matrix)
    center = [[(Q(0), Q(0)) for _ in range(size)] for _ in range(size)]
    errors = [[Q(0) for _ in range(size)] for _ in range(size)]
    for row in range(size):
        for column in range(row, size):
            real, imaginary = field_rectangle(matrix[row][column], bits)
            real_mid = (real[0] + real[1]) / 2
            imag_mid = (imaginary[0] + imaginary[1]) / 2
            if row == column:
                if not matrix[row][column].is_real():
                    raise VerificationError(f"{label} diagonal is not exactly real")
                imag_mid = Q(0)
            center[row][column] = (real_mid, imag_mid)
            center[column][row] = (real_mid, -imag_mid)
            error = (real[1] - real[0]) / 2 + (imaginary[1] - imaginary[0]) / 2
            errors[row][column] = errors[column][row] = error
    rho = max(sum(row) for row in errors)
    shifted = [row[:] for row in center]
    for index in range(size):
        shifted[index][index] = rc_sub(shifted[index][index], (rho, Q(0)))
    valid, pivots = rational_complex_ldl(shifted)
    if not valid:
        raise VerificationError(f"{label} interval-LDL certificate failed")
    return {
        "bits": bits,
        "entry_error_operator_bound": str(rho),
        "minimum_rational_ldl_pivot": str(min(pivots)),
        "method": "dyadic-entry-enclosure + row-sum-Weyl + rational-complex-LDL",
    }


def shifted(matrix: list[list[K]], value: Q) -> list[list[K]]:
    result = [row[:] for row in matrix]
    for index in range(len(result)):
        result[index][index] = result[index][index] + K.rational(value)
    return result


def matrix_subtract(left: list[list[K]], right: list[list[K]]) -> list[list[K]]:
    return [[left[row][column] - right[row][column] for column in range(len(left))] for row in range(len(left))]


def exact_joint_output(
    first: list[list[K]],
    second: list[list[K]],
    first_weights: list[Q],
    second_weights: list[Q],
) -> list[list[K]]:
    products = [first_weights[index] * second_weights[index] for index in range(4)]
    normalizer = sum(products)
    amplitudes = [[sqrt_rational(products[j] * products[k]) for k in range(4)] for j in range(4)]
    result = [[K0 for _ in range(9)] for _ in range(9)]
    for b1 in range(3):
        for b2 in range(3):
            row = 3 * b1 + b2
            for c1 in range(3):
                for c2 in range(3):
                    column = 3 * c1 + c2
                    value = K0
                    for j in range(4):
                        for k in range(4):
                            value = value + first[3 * j + b1][3 * k + c1] * second[3 * j + b2][3 * k + c2] * amplitudes[j][k]
                    result[row][column] = value.scale(1 / normalizer)
    return result


def rational_charpoly(matrix: list[list[Q]]) -> list[Q]:
    size = len(matrix)
    identity = [[Q(1) if row == column else Q(0) for column in range(size)] for row in range(size)]
    current = identity
    coefficients: list[Q] = []
    for k in range(1, size + 1):
        product = [
            [sum(matrix[row][t] * current[t][column] for t in range(size)) for column in range(size)]
            for row in range(size)
        ]
        coefficient = sum(product[index][index] for index in range(size)) / k
        coefficients.append(coefficient)
        current = [
            [product[row][column] - (coefficient if row == column else 0) for column in range(size)]
            for row in range(size)
        ]
    result = [Q(0)] * (size + 1)
    result[size] = Q(1)
    for k, coefficient in enumerate(coefficients, 1):
        result[size - k] = -coefficient
    return result


def polynomial_multiply(left: list[Q], right: list[Q]) -> list[Q]:
    result = [Q(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


SPECTRUM = [(Q(0), 1), (Q(5, 46), 2), (Q(18, 161), 2), (Q(51, 322), 2), (Q(39, 322), 2)]


def spectrum_polynomial() -> list[Q]:
    result = [Q(1)]
    for eigenvalue, multiplicity in SPECTRUM:
        for _ in range(multiplicity):
            result = polynomial_multiply(result, [-eigenvalue, Q(1)])
    return result


def root_bounds(value: Q, degree: int, bits: int) -> tuple[Q, Q]:
    denominator = 1 << bits
    target = value.numerator * (denominator**degree)
    low, high = 0, denominator
    while low < high:
        middle = (low + high + 1) // 2
        if (middle**degree) * value.denominator <= target:
            low = middle
        else:
            high = middle - 1
    lower = Q(low, denominator)
    if (low**degree) * value.denominator == target:
        return lower, lower
    return lower, Q(low + 1, denominator)


def endpoint_gap(bits: int = 224) -> dict[str, Any]:
    envelope = (Q(2, 3), Q(99097, 300000), Q(301, 100000))
    a_low = a_high = Q(0)
    for value in envelope:
        lower, upper = root_bounds(value, 22, bits)
        a_low += lower
        a_high += upper
    b_low = b_high = Q(0)
    for value, multiplicity in SPECTRUM:
        if value == 0:
            continue
        lower, upper = root_bounds(value, 22, bits)
        b_low += multiplicity * lower
        b_high += multiplicity * upper
    gap_low = a_low * a_low - b_high
    gap_high = a_high * a_high - b_low
    if gap_low <= 0:
        raise VerificationError("independent endpoint enclosure is not strict")
    return {
        "p": "1/22",
        "bits": bits,
        "trace_power_gap_lower": str(gap_low),
        "trace_power_gap_upper": str(gap_high),
        "strict": True,
    }


def read_frozen_inputs(
    floor_path: Path, joint_path: Path, cap_path: Path
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    floor_bytes = floor_path.read_bytes()
    joint_bytes = joint_path.read_bytes()
    if sha256(floor_bytes) != FLOOR_CANDIDATE_SHA256:
        raise VerificationError("floor candidate SHA-256 differs")
    if sha256(joint_bytes) != JOINT_CLAIM_SHA256:
        raise VerificationError("joint claim SHA-256 differs")
    cap_bytes = cap_path.read_bytes()
    if sha256(cap_bytes) != CAP_CANDIDATE_SHA256:
        raise VerificationError("cap candidate SHA-256 differs")
    return json.loads(floor_bytes), json.loads(joint_bytes), json.loads(cap_bytes)


def verify(
    floor_path: Path,
    joint_path: Path,
    cap_path: Path,
    bits: int,
    control: str = "none",
) -> dict[str, Any]:
    floor_candidate, joint_claim, cap_candidate = read_frozen_inputs(
        floor_path, joint_path, cap_path
    )

    r_vectors = printed_vectors(R_MATRICES, False)
    s_vectors = printed_vectors(S_MATRICES, False)
    sbar_vectors = printed_vectors(S_MATRICES, True)
    r_projection, r_norms = projector(r_vectors)
    s_projection, s_norms = projector(s_vectors)
    sbar_projection, sbar_norms = projector(sbar_vectors)

    orthogonal_product = matrix_multiply(r_projection, s_projection)
    if any(not entry.is_zero() for row in orthogonal_product for entry in row):
        raise VerificationError("printed R and S supports are not orthogonal")

    first_choi, first_weights = build_choi(r_projection)
    second_choi, second_weights = build_choi(sbar_projection)
    assert_hermitian(first_choi, "first normalized Choi matrix")
    assert_hermitian(second_choi, "second normalized Choi matrix")

    floor_q = parse_floor_q(floor_candidate)
    if control == "floor-decomposition":
        floor_q["A"][0][0] = floor_q["A"][0][0] + K.rational(1)
    floor_scale = Q(17, 16)
    floor_reports: dict[str, Any] = {}
    for name, choi in (("A", first_choi), ("B", second_choi)):
        q_matrix = [[entry.scale(floor_scale) for entry in row] for row in floor_q[name]]
        assert_hermitian(q_matrix, f"floor Q {name}")
        p_matrix = matrix_subtract(choi, partial_transpose_b(q_matrix))
        floor_reports[name] = {
            "P_minus_3_over_1000_I": prove_positive_definite(shifted(p_matrix, -Q(3, 1000)), f"floor P {name}", bits),
            "Q_minus_1_over_100000_I": prove_positive_definite(shifted(q_matrix, -Q(1, 100000)), f"floor Q {name}", bits),
        }

    cap, cap_q = parse_cap_q(cap_candidate)
    if control == "cap-decomposition":
        cap_q["A"][0][0] = cap_q["A"][0][0] + K.rational(1)
    if cap != Q(2, 3):
        raise VerificationError("cap is not 2/3")
    cap_reports: dict[str, Any] = {}
    for name, choi in (("A", first_choi), ("B", second_choi)):
        q_matrix = cap_q[name]
        assert_hermitian(q_matrix, f"cap Q {name}")
        p_matrix = matrix_subtract(shifted(zero_matrix(12), cap), choi)
        p_matrix = matrix_subtract(p_matrix, partial_transpose_b(q_matrix))
        cap_reports[name] = {
            "P_minus_1_over_500_I": prove_positive_definite(shifted(p_matrix, -Q(1, 500)), f"cap P {name}", bits),
            "Q_minus_1_over_500_I": prove_positive_definite(shifted(q_matrix, -Q(1, 500)), f"cap Q {name}", bits),
        }

    sigma = exact_joint_output(first_choi, second_choi, first_weights, second_weights)
    assert_hermitian(sigma, "weighted joint output")
    rational_sigma: list[list[Q]] = []
    for row in sigma:
        rational_sigma.append([entry.as_rational() for entry in row])
    claimed_matrix = [[Q(value) for value in row] for row in joint_claim.get("matrix", [])]
    if control == "joint-matrix":
        claimed_matrix[0][0] += Q(1, 1000)
    if rational_sigma != claimed_matrix:
        raise VerificationError("reconstructed joint matrix differs from the frozen claim")
    if sum(rational_sigma[index][index] for index in range(9)) != 1:
        raise VerificationError("joint output trace differs from one")
    if rational_charpoly(rational_sigma) != spectrum_polynomial():
        raise VerificationError("joint-output characteristic polynomial differs")

    return {
        "schema_version": "q4-independent-primitive-replay-v1",
        "valid": True,
        "control": control,
        "independence_boundary": {
            "imports_harness_code": False,
            "imports_prior_q4_verifier": False,
            "algebra_representation": "generic quotient-basis bitmasks",
            "positivity_path": "rational dyadic perturbation bound and rational-complex LDL",
        },
        "inputs": {
            "checkpoint_sha256": CHECKPOINT_SHA256,
            "floor_candidate_sha256": FLOOR_CANDIDATE_SHA256,
            "joint_claim_sha256": JOINT_CLAIM_SHA256,
            "cap_candidate_sha256": CAP_CANDIDATE_SHA256,
        },
        "printed_construction": {
            "R_norms": [str(value) for value in r_norms],
            "S_norms": [str(value) for value in s_norms],
            "Sbar_norms": [str(value) for value in sbar_norms],
            "R_S_supports_orthogonal": True,
            "normalized_choi_trace_preserving_exactly": True,
            "first_input_marginal": [str(value) for value in first_weights],
            "second_input_marginal": [str(value) for value in second_weights],
        },
        "floor": {
            "witness_scale": "17/16",
            "combined_floor": "301/100000",
            "channels": floor_reports,
        },
        "cap": {"value": "2/3", "channels": cap_reports},
        "joint_output": {
            "matrix_matches_frozen_claim_exactly": True,
            "trace": "1",
            "characteristic_polynomial_matches_spectrum": True,
            "spectrum": [
                {"value": str(value), "multiplicity": multiplicity}
                for value, multiplicity in SPECTRUM
            ],
        },
        "endpoint_control": endpoint_gap(),
        "claim_scope": (
            "Independently replays the shared floor/cap/joint-spectrum primitive layer "
            "and the strict p=1/22 endpoint. The all-p interval still uses the two "
            "separately maintained interval checkers."
        ),
    }


def self_test() -> None:
    omega = K.omega()
    if omega * omega * omega != K1 or omega.conjugate() != K.omega(True):
        raise VerificationError("field omega self-test failed")
    for value in (Q(2), Q(5, 2), Q(400, 81), Q(81, 16), Q(1)):
        root = sqrt_rational(value)
        if root * root != K.rational(value):
            raise VerificationError("field square-root self-test failed")
    positive = [[(Q(2), Q(0)), (Q(1), Q(1))], [(Q(1), Q(-1)), (Q(3), Q(0))]]
    negative = [[(Q(1), Q(0)), (Q(2), Q(0))], [(Q(2), Q(0)), (Q(1), Q(0))]]
    if not rational_complex_ldl(positive)[0] or rational_complex_ldl(negative)[0]:
        raise VerificationError("rational-complex LDL self-test failed")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--floor-candidate", type=Path, required=True)
    parser.add_argument("--joint-claim", type=Path, required=True)
    parser.add_argument("--cap-candidate", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument(
        "--control",
        choices=("none", "floor-decomposition", "cap-decomposition", "joint-matrix"),
        default="none",
        help="internal adverse control; every non-none value must be rejected",
    )
    arguments = parser.parse_args(argv)
    report: dict[str, Any]
    try:
        if arguments.bits < 80:
            raise VerificationError("at least 80 dyadic bits are required")
        self_test()
        report = verify(
            arguments.floor_candidate,
            arguments.joint_claim,
            arguments.cap_candidate,
            arguments.bits,
            arguments.control,
        )
        status = 0
    except (OSError, KeyError, ValueError, VerificationError) as exc:
        report = {
            "schema_version": "q4-independent-primitive-replay-v1",
            "valid": False,
            "control": arguments.control,
            "error": str(exc),
        }
        status = 1
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return status


if __name__ == "__main__":
    raise SystemExit(main())
