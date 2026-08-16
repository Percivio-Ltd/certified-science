"""Certified symmetric-rank-two checks for depolarizing-channel capacity.

The checker implements the Schur--Weyl block formula of
arXiv:2605.09138v2 without the paper code's interpolation basis.  A roots of
unity coefficient-extraction argument reduces every block entry to a finite
sum of rational numbers times square roots of positive rationals.  Square
roots are enclosed by exact rational bisection; the remaining spectral and
entropy enclosure uses exact dyadic replay, a polar-factor bound,
Gershgorin component counting, and rational logarithm bounds.

The participant payload is a little-endian ``complex128`` array with shape
``(2, n + 1)``.  Its rows are unnormalised Dicke-basis vectors.  The checker
interprets every component as the exact dyadic represented by its bytes and
normalises the resulting rank-two density by its exact trace.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import stat
from typing import Any, Mapping, Sequence

import numpy as np

from physics_proof_harness import depolarizing_capacity as base


SCHEMA = "physics-proof-harness/depolarizing-capacity-rank-two-v1"
RETURN_SCHEMA = f"{SCHEMA}/whole-return-v1"
TASK_ID = "depolarizing-capacity-threshold-improvement-rank-two-v1"
RUN_ID = "depolarizing-capacity-rank-two-n24-45-2026-07-23-v1"
CONDITION_ID = "exact-p0064912-full-symmetric-rank-two-n24-45-v1"
HANDOFF_ID = "owner-to-pro-depolarizing-capacity-rank-two-v1"
CHECKPOINT_ID = "checkpoint-0-prelaunch-rank-two-v1"
FAMILY_ID = "full-symmetric-rank-two-n24-45-v1"
TARGET_P = "0.064912"
TARGET_P_EXACT = Fraction(4057, 62500)
CONTROL_P = "0.063761"
CONTROL_P_EXACT = Fraction(63761, 1_000_000)
BASELINE_REPORTED_P = "0.06376"
ALLOWED_N = tuple(range(24, 46))
SQRT_BITS = 224
ACCUMULATION_DECIMAL_DIGITS = 50
LOG_GRID_DECIMAL_DIGITS = 70
LOG_TERMS = 96
MAX_PAYLOAD_BYTES = 2 * (45 + 1) * 16
MAX_RETURN_FILE_BYTES = 8_000_000
MAX_RETURN_TOTAL_BYTES = 12_000_000
EXPECTED_RETURN_FILES = (
    "RETURN.json",
    "coefficients.c128",
    "CHECKER-RESULT.json",
)
NONCLAIMS = (
    "no optimized threshold or family completeness",
    "no full quantum-capacity or zero-capacity theorem",
    "no superactivation claim",
    "no novelty, priority, or literature-completeness claim",
    "a no-hit is not evidence against a witness",
)


class RankTwoCapacityError(ValueError):
    """Raised when a candidate, certificate, or lineage binding fails closed."""


Gaussian = tuple[Fraction, Fraction]
ZERO: Gaussian = (Fraction(0), Fraction(0))


@dataclass(frozen=True)
class RealInterval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if self.lower > self.upper:
            raise RankTwoCapacityError("interval-direction")


@dataclass(frozen=True)
class ComplexBox:
    real: RealInterval
    imag: RealInterval


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_json_bytes(document: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("ascii")
        + b"\n"
    )


def _fraction_text(value: Fraction) -> str:
    return "0" if value == 0 else f"{value.numerator}/{value.denominator}"


def _directed_decimal(value: Fraction, *, upper: bool, digits: int = 50) -> str:
    scale = 10**digits
    if upper:
        scaled = -((-value.numerator * scale) // value.denominator)
    else:
        scaled = (value.numerator * scale) // value.denominator
    sign = "-" if scaled < 0 else ""
    absolute = abs(scaled)
    whole, fractional = divmod(absolute, scale)
    return f"{sign}{whole}.{fractional:0{digits}d}"


def _scaled_directed_integer(
    value: Fraction, *, upper: bool, digits: int = ACCUMULATION_DECIMAL_DIGITS
) -> int:
    scale = 10**digits
    if upper:
        return -((-value.numerator * scale) // value.denominator)
    return (value.numerator * scale) // value.denominator


def _gaussian_add(left: Gaussian, right: Gaussian) -> Gaussian:
    return left[0] + right[0], left[1] + right[1]


def _gaussian_conjugate(value: Gaussian) -> Gaussian:
    return value[0], -value[1]


def _gaussian_scale(value: Gaussian, factor: Fraction) -> Gaussian:
    return value[0] * factor, value[1] * factor


def _float_gaussian(value: complex) -> Gaussian:
    return Fraction.from_float(float(value.real)), Fraction.from_float(float(value.imag))


def _sqrt_bounds(value: Fraction, *, bits: int = SQRT_BITS) -> tuple[Fraction, Fraction]:
    if value < 0:
        raise RankTwoCapacityError("negative-square-root")
    if value == 0:
        return Fraction(0), Fraction(0)
    scale = 1 << bits
    target = value.numerator * scale * scale
    low = math.isqrt(target // value.denominator)
    while (low + 1) * (low + 1) * value.denominator <= target:
        low += 1
    while low * low * value.denominator > target:
        low -= 1
    lower = Fraction(low, scale)
    if low * low * value.denominator == target:
        return lower, lower
    return lower, Fraction(low + 1, scale)


def _abs_upper(value: Gaussian) -> Fraction:
    return _sqrt_bounds(value[0] * value[0] + value[1] * value[1])[1]


def _regular_file(path: Path, *, maximum_bytes: int) -> bytes:
    try:
        information = path.lstat()
    except FileNotFoundError as exc:
        raise RankTwoCapacityError("return-file-missing") from exc
    if not stat.S_ISREG(information.st_mode):
        raise RankTwoCapacityError("return-file-not-regular")
    if information.st_size > maximum_bytes:
        raise RankTwoCapacityError("return-file-oversized")
    value = path.read_bytes()
    if len(value) != information.st_size:
        raise RankTwoCapacityError("return-file-size-race")
    return value


def _load_json_bytes(value: bytes) -> dict[str, Any]:
    try:
        document = json.loads(value)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RankTwoCapacityError("json-invalid") from exc
    if not isinstance(document, dict):
        raise RankTwoCapacityError("json-not-object")
    if canonical_json_bytes(document) != value:
        raise RankTwoCapacityError("json-not-canonical")
    return document


def load_coefficients(payload: bytes, *, n: int) -> np.ndarray:
    if n not in ALLOWED_N:
        raise RankTwoCapacityError("block-length-out-of-family")
    expected = 2 * (n + 1) * 16
    if len(payload) != expected:
        raise RankTwoCapacityError("coefficient-payload-size")
    state = np.frombuffer(payload, dtype="<c16").astype(np.complex128, copy=True)
    state = state.reshape(2, n + 1)
    if not np.isfinite(state.real).all() or not np.isfinite(state.imag).all():
        raise RankTwoCapacityError("coefficient-nonfinite")
    return state


def _exact_state(state: np.ndarray) -> tuple[list[list[Gaussian]], Fraction]:
    exact = [
        [_float_gaussian(complex(state[row, column])) for column in range(state.shape[1])]
        for row in range(2)
    ]
    norm = sum(
        value[0] * value[0] + value[1] * value[1]
        for row in exact
        for value in row
    )
    if norm <= 0:
        raise RankTwoCapacityError("candidate-zero")
    gram00 = sum(value[0] * value[0] + value[1] * value[1] for value in exact[0])
    gram11 = sum(value[0] * value[0] + value[1] * value[1] for value in exact[1])
    gram01 = ZERO
    for left, right in zip(exact[0], exact[1], strict=True):
        product = base._gaussian_multiply(left, _gaussian_conjugate(right))
        gram01 = _gaussian_add(gram01, product)
    determinant = gram00 * gram11 - gram01[0] * gram01[0] - gram01[1] * gram01[1]
    if determinant <= 0:
        raise RankTwoCapacityError("candidate-not-rank-two")
    return exact, norm


@lru_cache(maxsize=None)
def _irrep_rational_coefficient(
    *,
    n: int,
    first_row: int,
    final_row: int,
    final_column: int,
    k: int,
    l: int,
    p: Fraction,
) -> Fraction:
    """Coefficient before the single positive square-root basis factor.

    With ``x=z_i`` and ``y=z_j^{-1}``,

    ``D_p(|0+x1><0+y1|)`` has entries
    ``1-2p+2pxy, (1-4p)y, (1-4p)x, 2p+(1-2p)xy`` and determinant
    ``2p(1-2p)(1+xy)^2``.  The roots-of-unity double Fourier contraction
    selects exactly the coefficient of ``x^k y^l``.
    """

    degree = 2 * first_row - n
    determinant_power = n - first_row
    raw_row = degree - final_column
    raw_column = degree - final_row
    if k - l != raw_row - raw_column:
        return Fraction(0)
    a0 = 1 - 2 * p
    a1 = 2 * p
    d0 = 2 * p
    d1 = 1 - 2 * p
    off = 1 - 4 * p
    determinant_scalar = 2 * p * (1 - 2 * p)
    total = Fraction(0)
    lower_s = max(0, raw_row + raw_column - degree)
    upper_s = min(raw_row, raw_column)
    for summation in range(lower_s, upper_s + 1):
        b_power = raw_column - summation
        c_power = raw_row - summation
        d_power = degree - raw_row - raw_column + summation
        prefix = (
            math.comb(raw_column, summation)
            * math.comb(degree - raw_column, raw_row - summation)
            * off ** (b_power + c_power)
        )
        for a_xy in range(summation + 1):
            a_term = (
                math.comb(summation, a_xy)
                * a0 ** (summation - a_xy)
                * a1**a_xy
            )
            for d_xy in range(d_power + 1):
                d_term = (
                    math.comb(d_power, d_xy)
                    * d0 ** (d_power - d_xy)
                    * d1**d_xy
                )
                determinant_xy = k - (a_xy + d_xy + c_power)
                if not 0 <= determinant_xy <= 2 * determinant_power:
                    continue
                if l != a_xy + d_xy + determinant_xy + b_power:
                    continue
                total += (
                    prefix
                    * a_term
                    * d_term
                    * math.comb(2 * determinant_power, determinant_xy)
                    * determinant_scalar**determinant_power
                )
    return total


@lru_cache(maxsize=None)
def _entry_kernel(
    n: int,
    first_row: int,
    final_row: int,
    final_column: int,
    p: Fraction,
) -> tuple[tuple[int, int, RealInterval], ...]:
    degree = 2 * first_row - n
    raw_row = degree - final_column
    raw_column = degree - final_row
    factorial_ratio = Fraction(
        math.factorial(raw_row) * math.factorial(degree - raw_row),
        math.factorial(raw_column) * math.factorial(degree - raw_column),
    )
    terms: list[tuple[int, int, RealInterval]] = []
    for k in range(n + 1):
        l = k - (raw_row - raw_column)
        if not 0 <= l <= n:
            continue
        rational = _irrep_rational_coefficient(
            n=n,
            first_row=first_row,
            final_row=final_row,
            final_column=final_column,
            k=k,
            l=l,
            p=p,
        )
        if rational == 0:
            continue
        square = factorial_ratio / (math.comb(n, k) * math.comb(n, l))
        sqrt_lower, sqrt_upper = _sqrt_bounds(square)
        if rational >= 0:
            factor = RealInterval(rational * sqrt_lower, rational * sqrt_upper)
        else:
            factor = RealInterval(rational * sqrt_upper, rational * sqrt_lower)
        terms.append((k, l, factor))
    return tuple(terms)


def _multiply_gaussian_real_interval(
    value: Gaussian, factor: RealInterval
) -> ComplexBox:
    real_products = (value[0] * factor.lower, value[0] * factor.upper)
    imag_products = (value[1] * factor.lower, value[1] * factor.upper)
    return ComplexBox(
        RealInterval(min(real_products), max(real_products)),
        RealInterval(min(imag_products), max(imag_products)),
    )


def _entry_box(
    exact_state: Sequence[Sequence[Gaussian]],
    norm: Fraction,
    *,
    n: int,
    first_row: int,
    reference_row: int,
    reference_column: int,
    final_row: int,
    final_column: int,
    p: Fraction,
) -> ComplexBox:
    re_lower = Fraction(0)
    re_upper = Fraction(0)
    im_lower = Fraction(0)
    im_upper = Fraction(0)
    for k, l, factor in _entry_kernel(
        n, first_row, final_row, final_column, p
    ):
        left = exact_state[reference_row][k]
        right = _gaussian_conjugate(exact_state[reference_column][l])
        density = base._gaussian_multiply(left, right)
        if density == ZERO:
            continue
        density = _gaussian_scale(density, Fraction(1, 1) / norm)
        term = _multiply_gaussian_real_interval(density, factor)
        re_lower += term.real.lower
        re_upper += term.real.upper
        im_lower += term.imag.lower
        im_upper += term.imag.upper
    return ComplexBox(
        RealInterval(re_lower, re_upper),
        RealInterval(im_lower, im_upper),
    )


def _box_center_radius(box: ComplexBox) -> tuple[Gaussian, Fraction]:
    center = (
        (box.real.lower + box.real.upper) / 2,
        (box.imag.lower + box.imag.upper) / 2,
    )
    re_radius = (box.real.upper - box.real.lower) / 2
    im_radius = (box.imag.upper - box.imag.lower) / 2
    return center, _sqrt_bounds(re_radius * re_radius + im_radius * im_radius)[1]


def _symmetrize_center_radius(
    center: list[list[Gaussian]], radius: list[list[Fraction]]
) -> tuple[list[list[Gaussian]], list[list[Fraction]]]:
    dimension = len(center)
    output_center = [[ZERO for _ in range(dimension)] for _ in range(dimension)]
    output_radius = [
        [Fraction(0) for _ in range(dimension)] for _ in range(dimension)
    ]
    for row in range(dimension):
        for column in range(row, dimension):
            combined = _gaussian_scale(
                _gaussian_add(center[row][column], _gaussian_conjugate(center[column][row])),
                Fraction(1, 2),
            )
            combined_radius = (radius[row][column] + radius[column][row]) / 2
            if row == column:
                combined_radius += abs(combined[1])
                combined = (combined[0], Fraction(0))
            output_center[row][column] = combined
            output_center[column][row] = _gaussian_conjugate(combined)
            output_radius[row][column] = combined_radius
            output_radius[column][row] = combined_radius
    return output_center, output_radius


def _reference_block(
    exact_state: Sequence[Sequence[Gaussian]],
    norm: Fraction,
    *,
    n: int,
    first_row: int,
    p: Fraction,
) -> tuple[list[list[Gaussian]], list[list[Fraction]]]:
    degree = 2 * first_row - n
    dimension = 2 * (degree + 1)
    center = [[ZERO for _ in range(dimension)] for _ in range(dimension)]
    radius = [[Fraction(0) for _ in range(dimension)] for _ in range(dimension)]
    for ref_row in range(2):
        for ref_column in range(2):
            for row in range(degree + 1):
                for column in range(degree + 1):
                    box = _entry_box(
                        exact_state,
                        norm,
                        n=n,
                        first_row=first_row,
                        reference_row=ref_row,
                        reference_column=ref_column,
                        final_row=row,
                        final_column=column,
                        p=p,
                    )
                    value, error = _box_center_radius(box)
                    joint_row = ref_row * (degree + 1) + row
                    joint_column = ref_column * (degree + 1) + column
                    center[joint_row][joint_column] = value
                    radius[joint_row][joint_column] = error
    return _symmetrize_center_radius(center, radius)


def _partial_trace_reference(
    joint_center: list[list[Gaussian]],
    joint_radius: list[list[Fraction]],
) -> tuple[list[list[Gaussian]], list[list[Fraction]]]:
    joint_dimension = len(joint_center)
    if joint_dimension % 2:
        raise RankTwoCapacityError("joint-dimension")
    dimension = joint_dimension // 2
    center = [[ZERO for _ in range(dimension)] for _ in range(dimension)]
    radius = [[Fraction(0) for _ in range(dimension)] for _ in range(dimension)]
    for row in range(dimension):
        for column in range(dimension):
            for reference in range(2):
                joint_row = reference * dimension + row
                joint_column = reference * dimension + column
                center[row][column] = _gaussian_add(
                    center[row][column], joint_center[joint_row][joint_column]
                )
                radius[row][column] += joint_radius[joint_row][joint_column]
    return _symmetrize_center_radius(center, radius)


def _frobenius_radius(radius: Sequence[Sequence[Fraction]]) -> Fraction:
    square = sum(value * value for row in radius for value in row)
    return _sqrt_bounds(square)[1]


def _matrix_add_error(
    center: Sequence[Sequence[Gaussian]],
    interval_error: Fraction,
) -> tuple[list[tuple[Fraction, Fraction, int]], dict[str, str]]:
    """Enclose the eigenvalues of a Hermitian matrix near an exact center."""

    dimension = len(center)
    numerical = np.empty((dimension, dimension), dtype=np.complex128)
    exact_numerical: list[list[Gaussian]] = []
    conversion_errors: list[list[Gaussian]] = []
    for row in range(dimension):
        numerical_row: list[Gaussian] = []
        error_row: list[Gaussian] = []
        for column in range(dimension):
            rounded = complex(float(center[row][column][0]), float(center[row][column][1]))
            numerical[row, column] = rounded
            exact_rounded = _float_gaussian(rounded)
            numerical_row.append(exact_rounded)
            error_row.append(
                (
                    center[row][column][0] - exact_rounded[0],
                    center[row][column][1] - exact_rounded[1],
                )
            )
        exact_numerical.append(numerical_row)
        conversion_errors.append(error_row)
    numerical = (numerical + numerical.conj().T) / 2
    exact_numerical = [
        [_float_gaussian(numerical[row, column]) for column in range(dimension)]
        for row in range(dimension)
    ]
    conversion_errors = [
        [
            (
                center[row][column][0] - exact_numerical[row][column][0],
                center[row][column][1] - exact_numerical[row][column][1],
            )
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]
    matrix_error = base._frobenius_upper(conversion_errors) + interval_error
    _, vectors = np.linalg.eigh(numerical)
    vector_exact = [
        [_float_gaussian(vectors[row, column]) for column in range(dimension)]
        for row in range(dimension)
    ]
    gram = base._matrix_multiply_exact(base._adjoint_exact(vector_exact), vector_exact)
    for index in range(dimension):
        gram[index][index] = _gaussian_add(
            gram[index][index], (Fraction(-1), Fraction(0))
        )
    orthogonality = base._frobenius_upper(gram)
    if orthogonality >= 1:
        raise RankTwoCapacityError("eigenbasis-not-invertible")
    sqrt_minus = _sqrt_bounds(1 - orthogonality)
    sqrt_plus = _sqrt_bounds(1 + orthogonality)
    polar_distance = max(1 - sqrt_minus[0], sqrt_plus[1] - 1)
    vector_norm = sqrt_plus[1]
    center_norm = base._frobenius_upper(exact_numerical)
    transform_error = (
        matrix_error
        + polar_distance * center_norm * (1 + vector_norm)
    )
    transformed = base._matrix_multiply_exact(
        base._matrix_multiply_exact(base._adjoint_exact(vector_exact), exact_numerical),
        vector_exact,
    )
    discs: list[tuple[Fraction, Fraction]] = []
    maximum_off_diagonal = Fraction(0)
    for row in range(dimension):
        diagonal = transformed[row][row]
        if diagonal[1] != 0:
            raise RankTwoCapacityError("transform-complex-diagonal")
        off_diagonal = sum(
            _abs_upper(transformed[row][column])
            for column in range(dimension)
            if column != row
        )
        maximum_off_diagonal = max(maximum_off_diagonal, off_diagonal)
        disc_radius = off_diagonal + dimension * transform_error
        discs.append((diagonal[0] - disc_radius, diagonal[0] + disc_radius))
    discs.sort()
    components: list[tuple[Fraction, Fraction, int]] = []
    for lower, upper in discs:
        if not components or lower > components[-1][1]:
            components.append((lower, upper, 1))
        else:
            old_lower, old_upper, count = components[-1]
            components[-1] = (old_lower, max(old_upper, upper), count + 1)
    return components, {
        "algebraic_entry_interval_frobenius": _fraction_text(interval_error),
        "center_to_binary64_frobenius": _fraction_text(
            base._frobenius_upper(conversion_errors)
        ),
        "total_matrix_frobenius": _fraction_text(matrix_error),
        "basis_orthogonality_frobenius": _fraction_text(orthogonality),
        "polar_unitary_distance": _fraction_text(polar_distance),
        "unitary_transform_spectral_error": _fraction_text(transform_error),
        "maximum_exact_off_diagonal_row_sum": _fraction_text(maximum_off_diagonal),
    }


def _entropy_from_components(
    components: Sequence[tuple[Fraction, Fraction, int]]
) -> tuple[tuple[Fraction, Fraction], list[dict[str, Any]]]:
    lower = Fraction(0)
    upper = Fraction(0)
    documents: list[dict[str, Any]] = []
    for component_lower, component_upper, count in components:
        term_lower, term_upper = _entropy_range(component_lower, component_upper)
        lower += count * term_lower
        upper += count * term_upper
        documents.append(
            {
                "lower": _fraction_text(component_lower),
                "upper": _fraction_text(component_upper),
                "multiplicity": count,
                "entropy_term_lower": _directed_decimal(term_lower, upper=False),
                "entropy_term_upper": _directed_decimal(term_upper, upper=True),
            }
        )
    return (lower, upper), documents


def _atanh_log_bounds(
    value: Fraction,
) -> tuple[Fraction, Fraction]:
    """Directed log bounds on ``[1,2]`` without giant rational GCDs."""

    if not 1 <= value <= 2:
        raise RankTwoCapacityError("log-reduced-domain")
    parameter = (value - 1) / (value + 1)
    if parameter == 0:
        return Fraction(0), Fraction(0)
    square = parameter * parameter
    scale = 10**LOG_GRID_DECIMAL_DIGITS
    lower_scaled = 0
    upper_scaled = 0
    term = parameter
    for index in range(LOG_TERMS):
        summand = 2 * term / (2 * index + 1)
        lower_scaled += (
            summand.numerator * scale // summand.denominator
        )
        upper_scaled += -(
            (-summand.numerator * scale) // summand.denominator
        )
        term *= square
    tail = 2 * term / ((2 * LOG_TERMS + 1) * (1 - square))
    upper_scaled += -((-tail.numerator * scale) // tail.denominator)
    return Fraction(lower_scaled, scale), Fraction(upper_scaled, scale)


@lru_cache(maxsize=None)
def _log_two_bounds() -> tuple[Fraction, Fraction]:
    return _atanh_log_bounds(Fraction(2))


def _log_bounds(value: Fraction) -> tuple[Fraction, Fraction]:
    if not 0 < value <= 1:
        raise RankTwoCapacityError("log-domain")
    if value == 1:
        return Fraction(0), Fraction(0)
    reduced = value
    powers = 0
    while reduced < 1:
        reduced *= 2
        powers += 1
    reduced_log = _atanh_log_bounds(reduced)
    log_two = _log_two_bounds()
    return (
        reduced_log[0] - powers * log_two[1],
        reduced_log[1] - powers * log_two[0],
    )


def _entropy_term_bounds(value: Fraction) -> tuple[Fraction, Fraction]:
    if value == 0:
        return Fraction(0), Fraction(0)
    logarithm = _log_bounds(value)
    numerator_lower = -value * logarithm[1]
    numerator_upper = -value * logarithm[0]
    log_two = _log_two_bounds()
    return numerator_lower / log_two[1], numerator_upper / log_two[0]


def _entropy_range(
    lower: Fraction, upper: Fraction
) -> tuple[Fraction, Fraction]:
    lower = max(Fraction(0), lower)
    upper = min(Fraction(1), upper)
    if lower > upper:
        raise RankTwoCapacityError("spectrum-outside-probability-domain")
    at_lower = _entropy_term_bounds(lower)
    at_upper = _entropy_term_bounds(upper)
    minimum = min(at_lower[0], at_upper[0])
    maximum = max(at_lower[1], at_upper[1])
    if (
        lower <= Fraction(367_880, 1_000_000)
        and upper >= Fraction(367_879, 1_000_000)
    ):
        maximum = max(maximum, Fraction(1))
    return minimum, maximum


def _specht_dimension(n: int, first_row: int) -> int:
    return (
        math.comb(n + 1, first_row + 1)
        * (2 * first_row - n + 1)
        // (n + 1)
    )


def certified_coherent_information(
    state: np.ndarray,
    *,
    p_per_pauli: str,
    candidate_id: str,
    coefficient_payload_sha256: str,
) -> dict[str, Any]:
    n = state.shape[1] - 1
    if n not in ALLOWED_N:
        raise RankTwoCapacityError("block-length-out-of-family")
    try:
        p = Fraction(Decimal(p_per_pauli))
    except (ArithmeticError, ValueError) as exc:
        raise RankTwoCapacityError("noise-parameter-invalid") from exc
    if not 0 <= p <= Fraction(1, 3):
        raise RankTwoCapacityError("noise-parameter-range")
    exact_state, norm = _exact_state(state)
    output_entropy_lower_scaled = 0
    output_entropy_upper_scaled = 0
    joint_entropy_lower_scaled = 0
    joint_entropy_upper_scaled = 0
    trace_center = Fraction(0)
    trace_radius = Fraction(0)
    blocks: list[dict[str, Any]] = []
    for first_row in range((n + 1) // 2, n + 1):
        degree = 2 * first_row - n
        joint_center, joint_radius = _reference_block(
            exact_state,
            norm,
            n=n,
            first_row=first_row,
            p=p,
        )
        output_center, output_radius = _partial_trace_reference(
            joint_center, joint_radius
        )
        output_interval_error = _frobenius_radius(output_radius)
        joint_interval_error = _frobenius_radius(joint_radius)
        output_components, output_budgets = _matrix_add_error(
            output_center, output_interval_error
        )
        joint_components, joint_budgets = _matrix_add_error(
            joint_center, joint_interval_error
        )
        output_entropy, output_component_documents = _entropy_from_components(
            output_components
        )
        joint_entropy, joint_component_documents = _entropy_from_components(
            joint_components
        )
        multiplicity = _specht_dimension(n, first_row)
        output_entropy_lower_scaled += _scaled_directed_integer(
            multiplicity * output_entropy[0], upper=False
        )
        output_entropy_upper_scaled += _scaled_directed_integer(
            multiplicity * output_entropy[1], upper=True
        )
        joint_entropy_lower_scaled += _scaled_directed_integer(
            multiplicity * joint_entropy[0], upper=False
        )
        joint_entropy_upper_scaled += _scaled_directed_integer(
            multiplicity * joint_entropy[1], upper=True
        )
        block_trace = sum(output_center[index][index][0] for index in range(degree + 1))
        block_trace_radius = sum(
            output_radius[index][index] for index in range(degree + 1)
        )
        trace_center += multiplicity * block_trace
        trace_radius += multiplicity * block_trace_radius
        blocks.append(
            {
                "first_row": first_row,
                "degree": degree,
                "specht_dimension": str(multiplicity),
                "trace_center": _fraction_text(block_trace),
                "trace_radius": _fraction_text(block_trace_radius),
                "output_spectrum_components": output_component_documents,
                "joint_spectrum_components": joint_component_documents,
                "output_error_budgets": output_budgets,
                "joint_error_budgets": joint_budgets,
            }
        )
    if not trace_center - trace_radius <= 1 <= trace_center + trace_radius:
        raise RankTwoCapacityError("global-trace-not-enclosed")
    accumulation_scale = 10**ACCUMULATION_DECIMAL_DIGITS
    output_entropy_lower = Fraction(
        output_entropy_lower_scaled, accumulation_scale
    )
    output_entropy_upper = Fraction(
        output_entropy_upper_scaled, accumulation_scale
    )
    joint_entropy_lower = Fraction(
        joint_entropy_lower_scaled, accumulation_scale
    )
    joint_entropy_upper = Fraction(
        joint_entropy_upper_scaled, accumulation_scale
    )
    coherent_lower = output_entropy_lower - joint_entropy_upper
    coherent_upper = output_entropy_upper - joint_entropy_lower
    classification = (
        "certified-positive"
        if coherent_lower > 0
        else "certified-negative"
        if coherent_upper < 0
        else "certified-straddles-zero"
    )
    return {
        "schema": f"{SCHEMA}/evaluation-v1",
        "task_id": TASK_ID,
        "family_id": FAMILY_ID,
        "candidate_id": candidate_id,
        "coefficient_payload_sha256": coefficient_payload_sha256,
        "n": n,
        "p_per_pauli": p_per_pauli,
        "p_exact": _fraction_text(p),
        "conventions": {
            "channel": "(1-3p)rho+p(XrhoX+YrhoY+ZrhoZ)",
            "p": "each-nonidentity-Pauli-probability",
            "input": "rank-two-density-from-two-Dicke-rows-normalized-by-exact-trace",
            "entropy": "von-Neumann-log-base-2-bits",
        },
        "classification": classification,
        "coherent_information_bits": {
            "lower": _directed_decimal(coherent_lower, upper=False),
            "upper": _directed_decimal(coherent_upper, upper=True),
        },
        "output_entropy_bits": {
            "lower": _directed_decimal(output_entropy_lower, upper=False),
            "upper": _directed_decimal(output_entropy_upper, upper=True),
        },
        "joint_entropy_bits": {
            "lower": _directed_decimal(joint_entropy_lower, upper=False),
            "upper": _directed_decimal(joint_entropy_upper, upper=True),
        },
        "global_trace": {
            "center": _fraction_text(trace_center),
            "radius": _fraction_text(trace_radius),
        },
        "blocks": blocks,
        "resources": {
            "irrep_blocks": len(blocks),
            "maximum_output_dimension": n + 1,
            "maximum_joint_dimension": 2 * (n + 1),
            "square_root_bits": SQRT_BITS,
            "directed_accumulation_decimal_digits": (
                ACCUMULATION_DECIMAL_DIGITS
            ),
        },
        "claim_boundary": (
            "candidate-specific positive coherent information only; no optimized "
            "threshold, family completeness, full capacity, zero-capacity, "
            "superactivation, novelty, or priority claim"
        ),
    }


def frozen_bindings(target_condition_sha256: str) -> dict[str, str]:
    return {
        "task_id": TASK_ID,
        "run_id": RUN_ID,
        "condition_id": CONDITION_ID,
        "handoff_id": HANDOFF_ID,
        "checkpoint_id": CHECKPOINT_ID,
        "family_id": FAMILY_ID,
        "target_condition_sha256": target_condition_sha256,
    }


def validate_whole_return(
    return_directory: Path,
    *,
    target_condition_sha256: str,
) -> dict[str, Any]:
    if not return_directory.is_dir() or return_directory.is_symlink():
        raise RankTwoCapacityError("return-directory-invalid")
    names = tuple(sorted(path.name for path in return_directory.iterdir()))
    if names != tuple(sorted(EXPECTED_RETURN_FILES)):
        raise RankTwoCapacityError("return-inventory")
    raw_files: dict[str, bytes] = {}
    total = 0
    for name in EXPECTED_RETURN_FILES:
        value = _regular_file(
            return_directory / name, maximum_bytes=MAX_RETURN_FILE_BYTES
        )
        raw_files[name] = value
        total += len(value)
    if total > MAX_RETURN_TOTAL_BYTES:
        raise RankTwoCapacityError("return-total-oversized")
    returned = _load_json_bytes(raw_files["RETURN.json"])
    expected_bindings = frozen_bindings(target_condition_sha256)
    if returned.get("schema") != RETURN_SCHEMA:
        raise RankTwoCapacityError("return-schema")
    if returned.get("bindings") != expected_bindings:
        raise RankTwoCapacityError("return-bindings")
    if returned.get("nonclaims") != list(NONCLAIMS):
        raise RankTwoCapacityError("return-nonclaims")
    if returned.get("target_p_per_pauli") != TARGET_P:
        raise RankTwoCapacityError("return-target-p")
    n = returned.get("n")
    if not isinstance(n, int) or n not in ALLOWED_N:
        raise RankTwoCapacityError("return-n")
    inventory = returned.get("files")
    if not isinstance(inventory, dict) or set(inventory) != {
        "coefficients.c128",
        "CHECKER-RESULT.json",
    }:
        raise RankTwoCapacityError("return-file-bindings")
    for name in inventory:
        if inventory[name] != {
            "sha256": sha256_bytes(raw_files[name]),
            "size_bytes": len(raw_files[name]),
        }:
            raise RankTwoCapacityError("return-file-binding-mismatch")
    state = load_coefficients(raw_files["coefficients.c128"], n=n)
    recomputed = certified_coherent_information(
        state,
        p_per_pauli=TARGET_P,
        candidate_id=str(returned.get("candidate_id")),
        coefficient_payload_sha256=sha256_bytes(raw_files["coefficients.c128"]),
    )
    if canonical_json_bytes(recomputed) != raw_files["CHECKER-RESULT.json"]:
        raise RankTwoCapacityError("returned-checker-result-mismatch")
    expected_verdict = (
        "VALID_CERTIFICATE"
        if recomputed["classification"] == "certified-positive"
        else "NO_VALID_CERTIFICATE"
    )
    if returned.get("verdict") != expected_verdict:
        raise RankTwoCapacityError("return-verdict")
    return {
        "schema": f"{SCHEMA}/validated-return-v1",
        "bindings": expected_bindings,
        "candidate_id": returned.get("candidate_id"),
        "n": n,
        "classification": recomputed["classification"],
        "verdict": expected_verdict,
        "coherent_information_bits": recomputed["coherent_information_bits"],
        "coefficient_payload_sha256": sha256_bytes(raw_files["coefficients.c128"]),
        "checker_result_sha256": sha256_bytes(raw_files["CHECKER-RESULT.json"]),
        "target_authorized": expected_verdict == "VALID_CERTIFICATE",
        "terminal_semantics": (
            "candidate-specific positive coherent information at exact "
            "p=4057/62500; no broader threshold-optimality or capacity claim"
        ),
    }
