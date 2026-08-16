"""Directed coherent-information checks for the qubit depolarizing channel.

The implementation evaluates rank-two inputs supported on the symmetric
subspace without constructing a 2**n dimensional matrix.  It uses the
Schur--Weyl blocks from arXiv:2605.09138v2 and a second tensor-power basis as a
roundoff adverse control.  The returned intervals are directed binary64
enclosures under the numerical contract below; they are not exact-arithmetic
proofs of an optimized threshold.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import io
import json
import math
from pathlib import Path
import pickletools
from typing import Any, Mapping
import zipfile

import numpy as np


SCHEMA = "physics-proof-harness/depolarizing-capacity/v1"
CHANNEL_ID = "qubit-pauli-depolarizing-per-pauli-v1"
CHANNEL_FORMULA = "(1-3*p)*rho+p*(X*rho*X+Y*rho*Y+Z*rho*Z)"
NOISE_PARAMETER = "p-is-each-nonidentity-Pauli-probability;total-error=3*p"
ENTROPY_CONVENTION = "von-Neumann-entropy-log-base-2-bits"
INPUT_CONVENTION = "two-by-(n+1)-complex-Dicke-amplitudes;rho=sum_i|v_i><v_i|"
TOLERANCE_FLOOR = 1.0e-12
DISCREPANCY_SAFETY_FACTOR = 8.0
CERTIFIED_DIRECT_MAX_N = 5
CERTIFIED_ROOT_BITS = 192
CERTIFIED_LOG_TERMS = 96
CERTIFIED_SCHEMA = f"{SCHEMA}/certified-direct-v1"
WHOLE_RETURN_SCHEMA = f"{SCHEMA}/whole-return-v1"
CERTIFIED_TASK_ID = "depolarizing-capacity-threshold-improvement-v1"
CERTIFIED_FAMILY_ID = "exact-sector-rank2-n5-v1"
CERTIFIED_P = "0.064912"
CERTIFIED_RUN_ID = "depolarizing-capacity-n5-certified-run-2026-07-23-v1"
CERTIFIED_CONDITION_ID = "exact-p0064912-sector-rank2-n5-v1"
CERTIFIED_HANDOFF_ID = "owner-to-pro-depolarizing-capacity-n5-v1"
CERTIFIED_CHECKPOINT_ID = "checkpoint-0-prelaunch-certified-v1"
RETURN_NONCLAIMS = [
    "no n=45 threshold claim",
    "no optimized threshold",
    "no full quantum capacity or zero-capacity premise",
    "no superactivation, novelty, or literature-completeness claim",
]


class CapacityCheckError(ValueError):
    """Raised when a source, convention, candidate, or interval fails closed."""


@dataclass(frozen=True)
class DirectedValue:
    lower: float
    value: float
    upper: float

    def document(self) -> dict[str, float]:
        return {"lower": self.lower, "value": self.value, "upper": self.upper}


Gaussian = tuple[Fraction, Fraction]
ZERO_GAUSSIAN: Gaussian = (Fraction(0), Fraction(0))


def _fraction_text(value: Fraction) -> str:
    return "0" if value == 0 else f"{value.numerator}/{value.denominator}"


def _directed_decimal(value: Fraction, *, upper: bool, digits: int = 50) -> str:
    scale = 10**digits
    scaled = (
        -((-value.numerator * scale) // value.denominator)
        if upper
        else (value.numerator * scale) // value.denominator
    )
    sign = "-" if scaled < 0 else ""
    absolute = abs(scaled)
    whole, fractional = divmod(absolute, scale)
    return f"{sign}{whole}.{fractional:0{digits}d}"


def _gaussian_add(left: Gaussian, right: Gaussian) -> Gaussian:
    return left[0] + right[0], left[1] + right[1]


def _gaussian_multiply(left: Gaussian, right: Gaussian) -> Gaussian:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def _gaussian_conjugate(value: Gaussian) -> Gaussian:
    return value[0], -value[1]


def _gaussian_scale(value: Gaussian, factor: Fraction) -> Gaussian:
    return value[0] * factor, value[1] * factor


def _float_gaussian(value: complex) -> Gaussian:
    return Fraction.from_float(float(value.real)), Fraction.from_float(float(value.imag))


def _sqrt_bounds_fraction(
    value: Fraction, *, bits: int = CERTIFIED_ROOT_BITS
) -> tuple[Fraction, Fraction]:
    if value < 0:
        raise CapacityCheckError("certified-square-root-negative")
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


def _gaussian_abs_upper(value: Gaussian) -> Fraction:
    return _sqrt_bounds_fraction(value[0] * value[0] + value[1] * value[1])[1]


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


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


def canonical_json_sha256(document: Mapping[str, Any]) -> str:
    return sha256_bytes(canonical_json_bytes(document))


def validate_exact_bindings(
    actual: Mapping[str, str], expected: Mapping[str, str]
) -> None:
    if dict(actual) != dict(expected):
        raise CapacityCheckError("task-candidate-binding-mismatch")


def load_torch_state_base64(
    path: Path, *, expected_sha256: str, n: int
) -> np.ndarray:
    """Load the one-tensor author artifact without executing its pickle.

    PyTorch's ZIP format stores the complex128 payload as a raw storage.  The
    pickle is inspected opcode-by-opcode only to reject an unexpected schema;
    it is never unpickled.
    """

    encoded = "".join(path.read_text(encoding="ascii").split())
    try:
        archive = base64.b64decode(encoded, validate=True)
    except ValueError as exc:
        raise CapacityCheckError("source-state-base64-invalid") from exc
    if sha256_bytes(archive) != expected_sha256:
        raise CapacityCheckError("source-state-sha256-mismatch")

    with zipfile.ZipFile(io.BytesIO(archive)) as source:
        names = source.namelist()
        pickle_names = [name for name in names if name.endswith("/data.pkl")]
        storage_names = [name for name in names if name.endswith("/data/0")]
        byteorder_names = [name for name in names if name.endswith("/byteorder")]
        version_names = [name for name in names if name.endswith("/version")]
        if not all(
            len(group) == 1
            for group in (
                pickle_names,
                storage_names,
                byteorder_names,
                version_names,
            )
        ):
            raise CapacityCheckError("source-state-zip-schema-mismatch")
        program = source.read(pickle_names[0])
        allowed_globals = {
            "torch._utils _rebuild_tensor_v2",
            "torch ComplexDoubleStorage",
            "collections OrderedDict",
        }
        seen_globals: set[str] = set()
        for opcode, argument, _ in pickletools.genops(program):
            if opcode.name in {"EXT1", "EXT2", "EXT4", "STACK_GLOBAL"}:
                raise CapacityCheckError("source-state-pickle-dynamic-global")
            if opcode.name == "GLOBAL":
                seen_globals.add(str(argument))
        if seen_globals != allowed_globals:
            raise CapacityCheckError("source-state-pickle-global-mismatch")
        if b"cpu" not in program or b"storage" not in program:
            raise CapacityCheckError("source-state-pickle-storage-mismatch")
        byteorder = source.read(byteorder_names[0]).decode("ascii")
        if byteorder not in {"little", "big"}:
            raise CapacityCheckError("source-state-byteorder-invalid")
        if source.read(version_names[0]).strip() != b"3":
            raise CapacityCheckError("source-state-version-invalid")
        storage = source.read(storage_names[0])

    expected_values = 2 * (n + 1)
    if len(storage) != expected_values * 16:
        raise CapacityCheckError("source-state-dimension-mismatch")
    dtype = "<c16" if byteorder == "little" else ">c16"
    state = np.frombuffer(storage, dtype=dtype).astype(np.complex128, copy=True)
    return state.reshape(2, n + 1)


def load_complex_storage_base64(
    path: Path, *, expected_sha256: str, n: int
) -> np.ndarray:
    """Load an exact extracted little-endian complex128 storage payload."""

    encoded = "".join(path.read_text(encoding="ascii").split())
    try:
        storage = base64.b64decode(encoded, validate=True)
    except ValueError as exc:
        raise CapacityCheckError("source-storage-base64-invalid") from exc
    if sha256_bytes(storage) != expected_sha256:
        raise CapacityCheckError("source-storage-sha256-mismatch")
    if len(storage) != 2 * (n + 1) * 16:
        raise CapacityCheckError("source-storage-dimension-mismatch")
    return (
        np.frombuffer(storage, dtype="<c16")
        .astype(np.complex128, copy=True)
        .reshape(2, n + 1)
    )


def validate_state(state: np.ndarray, *, n: int, atol: float = 5.0e-13) -> None:
    value = np.asarray(state, dtype=np.complex128)
    if value.shape != (2, n + 1):
        raise CapacityCheckError("channel-input-dimension-mismatch")
    if not np.isfinite(value.real).all() or not np.isfinite(value.imag).all():
        raise CapacityCheckError("candidate-nonfinite")
    trace = float(np.sum(np.abs(value) ** 2))
    if abs(trace - 1.0) > atol:
        raise CapacityCheckError("candidate-not-normalized")
    gram_eigenvalues = np.linalg.eigvalsh(value @ value.conj().T)
    if gram_eigenvalues[0] <= 1.0e-12:
        raise CapacityCheckError("candidate-not-rank-two")


def validate_conventions(conventions: Mapping[str, Any]) -> None:
    expected = {
        "channel_id": CHANNEL_ID,
        "channel_formula": CHANNEL_FORMULA,
        "noise_parameter": NOISE_PARAMETER,
        "entropy": ENTROPY_CONVENTION,
        "input": INPUT_CONVENTION,
    }
    if dict(conventions) != expected:
        raise CapacityCheckError("convention-binding-mismatch")


def conventions_document() -> dict[str, str]:
    return {
        "channel_id": CHANNEL_ID,
        "channel_formula": CHANNEL_FORMULA,
        "noise_parameter": NOISE_PARAMETER,
        "entropy": ENTROPY_CONVENTION,
        "input": INPUT_CONVENTION,
    }


def _fibonacci_basis(count: int, phase_offset: float) -> np.ndarray:
    indices = np.arange(count, dtype=np.float64)
    golden = (1.0 + math.sqrt(5.0)) / 2.0
    theta = np.arccos(1.0 - 2.0 * (indices + 0.5) / count)
    azimuth = 2.0 * math.pi * indices / golden + phase_offset
    return np.column_stack(
        [
            np.cos(theta / 2.0),
            np.exp(1j * azimuth) * np.sin(theta / 2.0),
        ]
    )


def _tensor_power_dicke_basis(single_qubits: np.ndarray, n: int) -> np.ndarray:
    result = np.empty((n + 1, n + 1), dtype=np.complex128)
    for row in range(n + 1):
        for weight in range(n + 1):
            result[row, weight] = (
                math.sqrt(math.comb(n, weight))
                * single_qubits[row, 0] ** (n - weight)
                * single_qubits[row, 1] ** weight
            )
    return result


def _dual_coefficients(basis: np.ndarray, density: np.ndarray) -> np.ndarray:
    left_solved = np.linalg.solve(basis.T, density)
    return np.linalg.solve(basis.conj().T, left_solved.T).T


def _symmetric_power_batch(matrices: np.ndarray, degree: int) -> np.ndarray:
    """Evaluate the degree-``degree`` polynomial GL(2) irrep in a fixed basis."""

    batch = matrices.shape[0]
    powers: list[np.ndarray] = []
    for entry in (
        matrices[:, 0, 0],
        matrices[:, 0, 1],
        matrices[:, 1, 0],
        matrices[:, 1, 1],
    ):
        table = np.empty((degree + 1, batch), dtype=np.complex128)
        table[0] = 1.0
        for exponent in range(1, degree + 1):
            table[exponent] = table[exponent - 1] * entry
        powers.append(table)
    a_power, b_power, c_power, d_power = powers

    result = np.empty(
        (batch, degree + 1, degree + 1), dtype=np.complex128
    )
    factorials = [math.factorial(index) for index in range(degree + 1)]
    for row in range(degree + 1):
        for column in range(degree + 1):
            value = np.zeros(batch, dtype=np.complex128)
            for summation in range(
                max(0, row + column - degree), min(row, column) + 1
            ):
                value += (
                    math.comb(column, summation)
                    * math.comb(degree - column, row - summation)
                    * a_power[summation]
                    * b_power[column - summation]
                    * c_power[row - summation]
                    * d_power[degree - row - column + summation]
                )
            scale = math.sqrt(
                factorials[row]
                * factorials[degree - row]
                / (factorials[column] * factorials[degree - column])
            )
            result[:, row, column] = scale * value

    # The paper's displayed formula indexes basis weights in descending order
    # and uses the transpose matrix convention.  Convert to ascending Dicke
    # weight so reference-system restacking preserves positivity.
    return result[:, ::-1, ::-1].transpose(0, 2, 1)


def _specht_dimension(n: int, first_row: int) -> int:
    return (
        math.comb(n + 1, first_row + 1)
        * (2 * first_row - n + 1)
        // (n + 1)
    )


def _channel_outer_products(
    single_qubits: np.ndarray, p_per_pauli: float
) -> np.ndarray:
    outer = np.einsum(
        "ai,bj->abij", single_qubits, single_qubits.conj()
    )
    traces = np.einsum("abii->ab", outer)
    output = (1.0 - 4.0 * p_per_pauli) * outer
    output[:, :, 0, 0] += 2.0 * p_per_pauli * traces
    output[:, :, 1, 1] += 2.0 * p_per_pauli * traces
    return output


def _blocks_for_basis(
    state: np.ndarray, p_per_pauli: float, phase_offset: float
) -> tuple[list[dict[str, Any]], float, float]:
    n = state.shape[1] - 1
    count = n + 1
    single_qubits = _fibonacci_basis(count, phase_offset)
    basis = _tensor_power_dicke_basis(single_qubits, n)
    condition_number = float(np.linalg.cond(basis))

    purification = np.einsum("ik,jl->ijkl", state, state.conj())
    dual = np.empty_like(purification)
    for left in range(2):
        for right in range(2):
            dual[left, right] = _dual_coefficients(
                basis, purification[left, right]
            )

    channel_outputs = _channel_outer_products(single_qubits, p_per_pauli)
    flattened_outputs = channel_outputs.reshape(count * count, 2, 2)
    determinants = (
        channel_outputs[:, :, 0, 0] * channel_outputs[:, :, 1, 1]
        - channel_outputs[:, :, 1, 0] * channel_outputs[:, :, 0, 1]
    )

    blocks: list[dict[str, Any]] = []
    total_weight = 0.0
    maximum_negative_mass = 0.0
    for first_row in range((n + 1) // 2, n + 1):
        degree = 2 * first_row - n
        irrep = _symmetric_power_batch(flattened_outputs, degree)
        irrep *= determinants.reshape(-1, 1, 1) ** (n - first_row)
        irrep = irrep.reshape(
            count, count, degree + 1, degree + 1
        )
        reference_blocks = np.einsum(
            "ijab,abxy->ijxy", dual, irrep, optimize=True
        )
        reference_blocks = (
            reference_blocks
            + reference_blocks.conj().transpose(1, 0, 3, 2)
        ) / 2.0
        output = reference_blocks[0, 0] + reference_blocks[1, 1]
        joint = reference_blocks.transpose(0, 2, 1, 3).reshape(
            2 * (degree + 1), 2 * (degree + 1)
        )
        trace = float(output.trace().real)
        if trace <= 0.0:
            raise CapacityCheckError("nonpositive-irrep-trace")
        output_normalized = (output + output.conj().T) / (2.0 * trace)
        joint_normalized = (joint + joint.conj().T) / (2.0 * trace)
        output_eigenvalues = np.linalg.eigvalsh(output_normalized).real
        joint_eigenvalues = np.linalg.eigvalsh(joint_normalized).real
        negative_mass = max(
            float(-np.minimum(output_eigenvalues, 0.0).sum()),
            float(-np.minimum(joint_eigenvalues, 0.0).sum()),
        )
        maximum_negative_mass = max(maximum_negative_mass, negative_mass)
        output_eigenvalues = np.maximum(output_eigenvalues, 0.0)
        joint_eigenvalues = np.maximum(joint_eigenvalues, 0.0)
        output_eigenvalues /= output_eigenvalues.sum()
        joint_eigenvalues /= joint_eigenvalues.sum()
        dimension = _specht_dimension(n, first_row)
        weight = dimension * trace
        total_weight += weight
        blocks.append(
            {
                "first_row": first_row,
                "degree": degree,
                "specht_dimension": dimension,
                "trace": trace,
                "weight": weight,
                "output_eigenvalues": output_eigenvalues,
                "joint_eigenvalues": joint_eigenvalues,
            }
        )
    return blocks, condition_number, maximum_negative_mass


def _entropy(values: np.ndarray) -> float:
    positive = values[values > 0.0]
    return float(-np.sum(positive * np.log2(positive)))


def _entropy_interval(values: np.ndarray, radius: float) -> DirectedValue:
    lower = 0.0
    upper = 0.0
    critical = 1.0 / math.e

    def term(value: float) -> float:
        if value <= 0.0 or value >= 1.0:
            return 0.0
        return -value * math.log2(value)

    for value in values:
        left = max(0.0, float(value) - radius)
        right = min(1.0, float(value) + radius)
        candidates = [term(left), term(right)]
        lower += min(candidates)
        if left <= critical <= right:
            candidates.append(term(critical))
        upper += max(candidates)
    central = _entropy(values)
    if not lower <= central <= upper:
        raise CapacityCheckError("entropy-interval-direction-failure")
    return DirectedValue(lower, central, upper)


def _product_interval(
    weight: DirectedValue, value: DirectedValue
) -> DirectedValue:
    products = (
        weight.lower * value.lower,
        weight.lower * value.upper,
        weight.upper * value.lower,
        weight.upper * value.upper,
    )
    central = weight.value * value.value
    return DirectedValue(min(products), central, max(products))


def evaluate_coherent_information(
    state: np.ndarray,
    *,
    p_per_pauli: float | str,
    candidate_id: str,
    source_state_sha256: str,
) -> dict[str, Any]:
    """Return a directed coherent-information enclosure and full block trace."""

    n = state.shape[1] - 1
    validate_state(state, n=n)
    p_decimal = (
        p_per_pauli
        if isinstance(p_per_pauli, str)
        else format(p_per_pauli, ".17g")
    )
    try:
        p_value = float(Decimal(p_decimal))
    except (ValueError, ArithmeticError) as exc:
        raise CapacityCheckError("noise-parameter-invalid-decimal") from exc
    if not math.isfinite(p_value) or not 0.0 <= p_value <= 1.0 / 3.0:
        raise CapacityCheckError("noise-parameter-out-of-range")

    primary, condition_primary, negative_primary = _blocks_for_basis(
        state, p_value, 0.0
    )
    secondary, condition_secondary, negative_secondary = _blocks_for_basis(
        state, p_value, 0.371
    )
    if len(primary) != len(secondary):
        raise CapacityCheckError("dual-basis-block-count-mismatch")

    contributions: list[DirectedValue] = []
    block_documents: list[dict[str, Any]] = []
    maximum_spectral_discrepancy = 0.0
    maximum_entropy_radius = 0.0
    for first, second in zip(primary, secondary, strict=True):
        if first["first_row"] != second["first_row"]:
            raise CapacityCheckError("dual-basis-block-identity-mismatch")
        output_discrepancy = float(
            np.max(
                np.abs(
                    first["output_eigenvalues"]
                    - second["output_eigenvalues"]
                )
            )
        )
        joint_discrepancy = float(
            np.max(
                np.abs(
                    first["joint_eigenvalues"]
                    - second["joint_eigenvalues"]
                )
            )
        )
        spectral_discrepancy = max(output_discrepancy, joint_discrepancy)
        maximum_spectral_discrepancy = max(
            maximum_spectral_discrepancy, spectral_discrepancy
        )
        radius = max(
            TOLERANCE_FLOOR,
            DISCREPANCY_SAFETY_FACTOR * spectral_discrepancy,
            DISCREPANCY_SAFETY_FACTOR * negative_primary,
            DISCREPANCY_SAFETY_FACTOR * negative_secondary,
        )
        maximum_entropy_radius = max(maximum_entropy_radius, radius)
        output_entropy = _entropy_interval(
            first["output_eigenvalues"], radius
        )
        joint_entropy = _entropy_interval(
            first["joint_eigenvalues"], radius
        )
        entropy_difference = DirectedValue(
            output_entropy.lower - joint_entropy.upper,
            output_entropy.value - joint_entropy.value,
            output_entropy.upper - joint_entropy.lower,
        )
        weight_discrepancy = abs(first["weight"] - second["weight"])
        weight_radius = max(
            2.0e-14, DISCREPANCY_SAFETY_FACTOR * weight_discrepancy
        )
        weight = DirectedValue(
            max(0.0, first["weight"] - weight_radius),
            first["weight"],
            first["weight"] + weight_radius,
        )
        contribution = _product_interval(weight, entropy_difference)
        contributions.append(contribution)
        block_documents.append(
            {
                "first_row": first["first_row"],
                "degree": first["degree"],
                "specht_dimension": str(first["specht_dimension"]),
                "weight": weight.document(),
                "spectrum_radius": radius,
                "output_spectrum": first["output_eigenvalues"].tolist(),
                "joint_spectrum": first["joint_eigenvalues"].tolist(),
                "output_entropy_bits": output_entropy.document(),
                "joint_entropy_bits": joint_entropy.document(),
                "contribution_bits": contribution.document(),
            }
        )

    total = DirectedValue(
        math.fsum(item.lower for item in contributions),
        math.fsum(item.value for item in contributions),
        math.fsum(item.upper for item in contributions),
    )
    total_weight = math.fsum(item["weight"] for item in primary)
    if abs(total_weight - 1.0) > 5.0e-12:
        raise CapacityCheckError("irrep-weights-not-normalized")
    if not total.lower <= total.value <= total.upper:
        raise CapacityCheckError("coherent-information-interval-direction-failure")
    classification = (
        "directed-positive"
        if total.lower > 0.0
        else "directed-negative"
        if total.upper < 0.0
        else "interval-straddles-zero"
    )
    return {
        "schema": f"{SCHEMA}/evaluation",
        "candidate_id": candidate_id,
        "source_state_sha256": source_state_sha256,
        "n": n,
        "channel_uses": n,
        "input_dimension": 2**n,
        "input_representation_dimension": n + 1,
        "code_rank": 2,
        "p_per_pauli": p_decimal,
        "total_nonidentity_probability": str(Decimal(p_decimal) * 3),
        "conventions": conventions_document(),
        "coherent_information_bits": total.document(),
        "classification": classification,
        "trace": {
            "irrep_weight_sum": total_weight,
            "blocks": block_documents,
        },
        "numerics": {
            "semantics": (
                "directed-binary64-dual-basis-enclosure;not-exact-arithmetic"
            ),
            "primary_basis_condition_number": condition_primary,
            "secondary_basis_condition_number": condition_secondary,
            "maximum_spectral_discrepancy": maximum_spectral_discrepancy,
            "maximum_spectrum_radius": maximum_entropy_radius,
            "maximum_negative_spectral_mass": max(
                negative_primary, negative_secondary
            ),
            "discrepancy_safety_factor": DISCREPANCY_SAFETY_FACTOR,
            "spectrum_radius_floor": TOLERANCE_FLOOR,
        },
        "cost": {
            "tensor_power_bases": 2,
            "tensor_basis_size": n + 1,
            "channel_outer_products_per_basis": (n + 1) ** 2,
            "irrep_count": len(primary),
            "sum_output_block_dimensions": sum(
                item["degree"] + 1 for item in primary
            ),
            "sum_joint_block_dimensions": 2
            * sum(item["degree"] + 1 for item in primary),
            "maximum_output_block_dimension": n + 1,
            "maximum_joint_block_dimension": 2 * (n + 1),
        },
        "claim_boundary": (
            "candidate-specific coherent information only; no optimized "
            "threshold, full quantum capacity, novelty, or zero-capacity claim"
        ),
    }


def compare_decimal_strict(candidate: str, baseline: str) -> bool:
    """Strict, non-rounded threshold comparison for bound target documents."""

    return Decimal(candidate) > Decimal(baseline)


def validate_evaluation(document: Mapping[str, Any]) -> None:
    if document.get("schema") != f"{SCHEMA}/evaluation":
        raise CapacityCheckError("evaluation-schema-mismatch")
    validate_conventions(document.get("conventions", {}))
    interval = document.get("coherent_information_bits")
    if not isinstance(interval, Mapping):
        raise CapacityCheckError("evaluation-interval-missing")
    lower = float(interval["lower"])
    value = float(interval["value"])
    upper = float(interval["upper"])
    if not lower <= value <= upper:
        raise CapacityCheckError("evaluation-interval-direction-failure")
    expected = (
        "directed-positive"
        if lower > 0.0
        else "directed-negative"
        if upper < 0.0
        else "interval-straddles-zero"
    )
    if document.get("classification") != expected:
        raise CapacityCheckError("evaluation-classification-mismatch")


def validate_candidate_family(
    state: np.ndarray,
    *,
    family_id: str,
    candidate_id: str,
    source_state_sha256: str | None,
    allowed_author_seeds: Mapping[str, str],
    zero_atol: float = 1.0e-14,
) -> None:
    n = state.shape[1] - 1
    validate_state(state, n=n)
    if family_id == "dense-rank2-symmetric-n45-v1":
        if n != 45:
            raise CapacityCheckError("dense-family-wrong-block-length")
        if candidate_id.startswith("author-"):
            expected = allowed_author_seeds.get(candidate_id)
            if expected is None or source_state_sha256 != expected:
                raise CapacityCheckError("fabricated-author-family-membership")
        return
    if family_id == "sparse-two-dicke-weights-n45-v1":
        if n != 45:
            raise CapacityCheckError("sparse-family-wrong-block-length")
        active = np.flatnonzero(np.max(np.abs(state), axis=0) > zero_atol)
        if len(active) != 2:
            raise CapacityCheckError("sparse-family-support-mismatch")
        return
    raise CapacityCheckError("unknown-family-id")


def repetition_state(n: int) -> np.ndarray:
    state = np.zeros((2, n + 1), dtype=np.complex128)
    state[0, 0] = 1.0 / math.sqrt(2.0)
    state[1, n] = 1.0 / math.sqrt(2.0)
    return state


def certified_sector_repetition_state(n: int) -> np.ndarray:
    """Return exact-dyadic per-string sector amplitudes for a repetition code."""

    state = np.zeros((2, n + 1), dtype=np.complex128)
    state[0, 0] = 1.0
    state[1, n] = 1.0
    return state


def _certified_sector_state(
    state: np.ndarray,
) -> tuple[list[list[Gaussian]], Fraction]:
    value = np.asarray(state, dtype=np.complex128)
    if value.ndim != 2 or value.shape[0] != 2:
        raise CapacityCheckError("certified-sector-state-shape")
    n = value.shape[1] - 1
    if not 1 <= n <= CERTIFIED_DIRECT_MAX_N:
        raise CapacityCheckError("certified-direct-block-length-out-of-range")
    if not np.isfinite(value.real).all() or not np.isfinite(value.imag).all():
        raise CapacityCheckError("certified-sector-state-nonfinite")
    exact = [[_float_gaussian(value[row, weight]) for weight in range(n + 1)] for row in range(2)]
    normalization = Fraction(0)
    gram = [[ZERO_GAUSSIAN for _ in range(2)] for _ in range(2)]
    for left in range(2):
        for right in range(2):
            entry = ZERO_GAUSSIAN
            for weight in range(n + 1):
                product = _gaussian_multiply(
                    exact[left][weight],
                    _gaussian_conjugate(exact[right][weight]),
                )
                entry = _gaussian_add(
                    entry,
                    _gaussian_scale(product, Fraction(math.comb(n, weight))),
                )
            gram[left][right] = entry
    if gram[0][0][1] != 0 or gram[1][1][1] != 0:
        raise CapacityCheckError("certified-sector-gram-not-hermitian")
    normalization = gram[0][0][0] + gram[1][1][0]
    if normalization <= 0:
        raise CapacityCheckError("certified-sector-zero-state")
    determinant = _gaussian_add(
        _gaussian_multiply(gram[0][0], gram[1][1]),
        _gaussian_scale(_gaussian_multiply(gram[0][1], gram[1][0]), Fraction(-1)),
    )
    if determinant[1] != 0 or determinant[0] <= 0:
        raise CapacityCheckError("certified-sector-state-not-rank-two")
    return exact, normalization


def _apply_exact_depolarizing(
    operator: list[list[Gaussian]],
    *,
    n: int,
    p: Fraction,
) -> list[list[Gaussian]]:
    dimension = 1 << n
    if len(operator) != dimension or any(len(row) != dimension for row in operator):
        raise CapacityCheckError("certified-channel-operator-shape")
    same = 1 - 2 * p
    flip = 2 * p
    off_diagonal = 1 - 4 * p
    output = [[ZERO_GAUSSIAN for _ in range(dimension)] for _ in range(dimension)]
    for left in range(dimension):
        for right in range(dimension):
            coefficient = operator[left][right]
            if coefficient == ZERO_GAUSSIAN:
                continue
            choices: list[tuple[int, int, Fraction]] = [(0, 0, Fraction(1))]
            for bit in range(n):
                left_bit = (left >> bit) & 1
                right_bit = (right >> bit) & 1
                next_choices: list[tuple[int, int, Fraction]] = []
                if left_bit != right_bit:
                    for out_left, out_right, factor in choices:
                        next_choices.append(
                            (
                                out_left | (left_bit << bit),
                                out_right | (right_bit << bit),
                                factor * off_diagonal,
                            )
                        )
                else:
                    for out_left, out_right, factor in choices:
                        next_choices.append(
                            (
                                out_left | (left_bit << bit),
                                out_right | (right_bit << bit),
                                factor * same,
                            )
                        )
                        changed = 1 - left_bit
                        next_choices.append(
                            (
                                out_left | (changed << bit),
                                out_right | (changed << bit),
                                factor * flip,
                            )
                        )
                choices = next_choices
            for out_left, out_right, factor in choices:
                output[out_left][out_right] = _gaussian_add(
                    output[out_left][out_right],
                    _gaussian_scale(coefficient, factor),
                )
    return output


def _exact_output_and_joint(
    state: np.ndarray,
    p: Fraction,
) -> tuple[list[list[Gaussian]], list[list[Gaussian]], Fraction]:
    exact, normalization = _certified_sector_state(state)
    n = state.shape[1] - 1
    dimension = 1 << n
    blocks: list[list[list[list[Gaussian]]]] = []
    for reference_left in range(2):
        row_blocks: list[list[list[Gaussian]]] = []
        for reference_right in range(2):
            operator = [[ZERO_GAUSSIAN for _ in range(dimension)] for _ in range(dimension)]
            for left in range(dimension):
                left_value = exact[reference_left][left.bit_count()]
                for right in range(dimension):
                    right_value = exact[reference_right][right.bit_count()]
                    operator[left][right] = _gaussian_scale(
                        _gaussian_multiply(left_value, _gaussian_conjugate(right_value)),
                        Fraction(1, 1) / normalization,
                    )
            row_blocks.append(_apply_exact_depolarizing(operator, n=n, p=p))
        blocks.append(row_blocks)
    output = [
        [
            _gaussian_add(blocks[0][0][left][right], blocks[1][1][left][right])
            for right in range(dimension)
        ]
        for left in range(dimension)
    ]
    joint_dimension = 2 * dimension
    joint = [[ZERO_GAUSSIAN for _ in range(joint_dimension)] for _ in range(joint_dimension)]
    for reference_left in range(2):
        for reference_right in range(2):
            for left in range(dimension):
                for right in range(dimension):
                    joint[reference_left * dimension + left][reference_right * dimension + right] = (
                        blocks[reference_left][reference_right][left][right]
                    )
    return output, joint, normalization


def _verify_exact_density(matrix: list[list[Gaussian]]) -> None:
    dimension = len(matrix)
    if not dimension or any(len(row) != dimension for row in matrix):
        raise CapacityCheckError("certified-density-shape")
    trace = Fraction(0)
    for row in range(dimension):
        for column in range(dimension):
            if matrix[row][column] != _gaussian_conjugate(matrix[column][row]):
                raise CapacityCheckError("certified-density-not-hermitian")
        if matrix[row][row][1] != 0:
            raise CapacityCheckError("certified-density-complex-diagonal")
        trace += matrix[row][row][0]
    if trace != 1:
        raise CapacityCheckError("certified-density-trace-not-one")


def _matrix_multiply_exact(
    left: list[list[Gaussian]],
    right: list[list[Gaussian]],
) -> list[list[Gaussian]]:
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    if any(len(row) != middle for row in left) or any(len(row) != columns for row in right):
        raise CapacityCheckError("certified-matrix-product-shape")
    output = [[ZERO_GAUSSIAN for _ in range(columns)] for _ in range(rows)]
    for row in range(rows):
        for column in range(columns):
            value = ZERO_GAUSSIAN
            for index in range(middle):
                value = _gaussian_add(
                    value,
                    _gaussian_multiply(left[row][index], right[index][column]),
                )
            output[row][column] = value
    return output


def _adjoint_exact(matrix: list[list[Gaussian]]) -> list[list[Gaussian]]:
    return [
        [_gaussian_conjugate(matrix[column][row]) for column in range(len(matrix))]
        for row in range(len(matrix[0]))
    ]


def _frobenius_upper(matrix: list[list[Gaussian]]) -> Fraction:
    square = sum(
        value[0] * value[0] + value[1] * value[1]
        for row in matrix
        for value in row
    )
    return _sqrt_bounds_fraction(square)[1]


def _atanh_log_bounds(value: Fraction, terms: int) -> tuple[Fraction, Fraction]:
    if not 1 <= value <= 2:
        raise CapacityCheckError("certified-reduced-log-domain")
    parameter = (value - 1) / (value + 1)
    square = parameter * parameter
    partial = Fraction(0)
    term = parameter
    for index in range(terms):
        partial += term / (2 * index + 1)
        term *= square
    lower = 2 * partial
    if parameter == 0:
        return lower, lower
    tail = 2 * term / ((2 * terms + 1) * (1 - square))
    return lower, lower + tail


def _log_bounds_fraction(value: Fraction) -> tuple[Fraction, Fraction]:
    if not 0 < value <= 1:
        raise CapacityCheckError("certified-log-domain")
    if value == 1:
        return Fraction(0), Fraction(0)
    reduced = value
    powers = 0
    while reduced < 1:
        reduced *= 2
        powers += 1
    reduced_log = _atanh_log_bounds(reduced, CERTIFIED_LOG_TERMS)
    log_two = _atanh_log_bounds(Fraction(2), CERTIFIED_LOG_TERMS)
    return (
        reduced_log[0] - powers * log_two[1],
        reduced_log[1] - powers * log_two[0],
    )


def _entropy_term_bounds(value: Fraction) -> tuple[Fraction, Fraction]:
    if value == 0:
        return Fraction(0), Fraction(0)
    logarithm = _log_bounds_fraction(value)
    numerator_lower = -value * logarithm[1]
    numerator_upper = -value * logarithm[0]
    log_two = _atanh_log_bounds(Fraction(2), CERTIFIED_LOG_TERMS)
    return numerator_lower / log_two[1], numerator_upper / log_two[0]


def _entropy_range(lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    lower = max(Fraction(0), lower)
    upper = min(Fraction(1), upper)
    if lower > upper:
        raise CapacityCheckError("certified-spectrum-outside-probability-domain")
    at_lower = _entropy_term_bounds(lower)
    at_upper = _entropy_term_bounds(upper)
    minimum = min(at_lower[0], at_upper[0])
    maximum = max(at_lower[1], at_upper[1])
    # 1/e is bracketed by these exact decimals.  A crossing interval uses the
    # harmless universal one-bit upper bound rather than a transcendental
    # stationary-point approximation.
    e_inverse_lower = Fraction(367_879, 1_000_000)
    e_inverse_upper = Fraction(367_880, 1_000_000)
    if lower <= e_inverse_upper and upper >= e_inverse_lower:
        maximum = max(maximum, Fraction(1))
    return minimum, maximum


def _certified_spectrum_components(
    matrix: list[list[Gaussian]],
) -> tuple[list[tuple[Fraction, Fraction, int]], dict[str, str]]:
    _verify_exact_density(matrix)
    dimension = len(matrix)
    center = np.empty((dimension, dimension), dtype=np.complex128)
    center_exact: list[list[Gaussian]] = []
    entry_error: list[list[Gaussian]] = []
    for row in range(dimension):
        center_row: list[Gaussian] = []
        error_row: list[Gaussian] = []
        for column in range(dimension):
            value = matrix[row][column]
            rounded = complex(float(value[0]), float(value[1]))
            center[row, column] = rounded
            exact_rounded = _float_gaussian(rounded)
            center_row.append(exact_rounded)
            error_row.append((value[0] - exact_rounded[0], value[1] - exact_rounded[1]))
        center_exact.append(center_row)
        entry_error.append(error_row)
    center = (center + center.conj().T) / 2
    center_exact = [[_float_gaussian(center[row, column]) for column in range(dimension)] for row in range(dimension)]
    entry_error = [
        [
            (
                matrix[row][column][0] - center_exact[row][column][0],
                matrix[row][column][1] - center_exact[row][column][1],
            )
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]
    matrix_error = _frobenius_upper(entry_error)
    _, vectors = np.linalg.eigh(center)
    vector_exact = [[_float_gaussian(vectors[row, column]) for column in range(dimension)] for row in range(dimension)]
    gram = _matrix_multiply_exact(_adjoint_exact(vector_exact), vector_exact)
    for index in range(dimension):
        gram[index][index] = _gaussian_add(gram[index][index], (Fraction(-1), Fraction(0)))
    orthogonality_defect = _frobenius_upper(gram)
    if orthogonality_defect >= 1:
        raise CapacityCheckError("certified-eigenbasis-not-invertible")
    sqrt_one_minus = _sqrt_bounds_fraction(1 - orthogonality_defect)
    sqrt_one_plus = _sqrt_bounds_fraction(1 + orthogonality_defect)
    polar_distance = max(1 - sqrt_one_minus[0], sqrt_one_plus[1] - 1)
    vector_norm = sqrt_one_plus[1]
    center_norm = _frobenius_upper(center_exact)
    transform_error = (
        matrix_error
        + polar_distance * center_norm * (1 + vector_norm)
    )
    transformed = _matrix_multiply_exact(
        _matrix_multiply_exact(_adjoint_exact(vector_exact), center_exact),
        vector_exact,
    )
    discs: list[tuple[Fraction, Fraction]] = []
    maximum_off_diagonal = Fraction(0)
    for row in range(dimension):
        if transformed[row][row][1] != 0:
            raise CapacityCheckError("certified-transform-complex-diagonal")
        off_diagonal = sum(
            _gaussian_abs_upper(transformed[row][column])
            for column in range(dimension)
            if column != row
        )
        maximum_off_diagonal = max(maximum_off_diagonal, off_diagonal)
        radius = off_diagonal + dimension * transform_error
        center_value = transformed[row][row][0]
        discs.append((center_value - radius, center_value + radius))
    discs.sort()
    components: list[tuple[Fraction, Fraction, int]] = []
    for lower, upper in discs:
        if not components or lower > components[-1][1]:
            components.append((lower, upper, 1))
        else:
            old_lower, old_upper, count = components[-1]
            components[-1] = (old_lower, max(old_upper, upper), count + 1)
    if sum(count for _, _, count in components) != dimension:
        raise CapacityCheckError("certified-spectrum-multiplicity")
    return components, {
        "exact_matrix_to_binary64_frobenius": _fraction_text(matrix_error),
        "approximate_basis_orthogonality_frobenius": _fraction_text(orthogonality_defect),
        "polar_unitary_distance": _fraction_text(polar_distance),
        "unitary_transform_spectral_error": _fraction_text(transform_error),
        "maximum_exact_off_diagonal_row_sum": _fraction_text(maximum_off_diagonal),
    }


def _certified_entropy(
    matrix: list[list[Gaussian]],
) -> tuple[tuple[Fraction, Fraction], dict[str, Any]]:
    components, budgets = _certified_spectrum_components(matrix)
    lower = Fraction(0)
    upper = Fraction(0)
    component_documents: list[dict[str, Any]] = []
    for component_lower, component_upper, count in components:
        term_lower, term_upper = _entropy_range(component_lower, component_upper)
        lower += count * term_lower
        upper += count * term_upper
        component_documents.append(
            {
                "lower": _fraction_text(component_lower),
                "upper": _fraction_text(component_upper),
                "multiplicity": count,
                "entropy_term_lower": _directed_decimal(term_lower, upper=False),
                "entropy_term_upper": _directed_decimal(term_upper, upper=True),
            }
        )
    return (lower, upper), {"components": component_documents, "error_budgets": budgets}


def certified_coherent_information(
    state: np.ndarray,
    *,
    p_per_pauli: str,
    candidate_id: str,
    coefficient_payload_sha256: str,
    family_id: str = "exact-sector-rank2-n5-v1",
) -> dict[str, Any]:
    """Rigorous direct-state enclosure for the frozen exact n<=5 family.

    Each coefficient is interpreted as the exact binary rational represented
    by its complex128 bytes.  It is the common computational-basis amplitude
    for strings in that weight sector.  The checker normalizes by the exact
    binomially weighted squared norm, applies the exact rational tensor channel,
    and encloses Hermitian spectra through an exact polar/Gershgorin argument.
    """

    try:
        p = Fraction(Decimal(p_per_pauli))
    except (ArithmeticError, ValueError) as exc:
        raise CapacityCheckError("certified-noise-parameter-invalid") from exc
    if not 0 <= p <= Fraction(1, 3):
        raise CapacityCheckError("certified-noise-parameter-range")
    n = state.shape[1] - 1
    if family_id != "exact-sector-rank2-n5-v1" or n != 5:
        raise CapacityCheckError("certified-family-binding")
    output, joint, normalization = _exact_output_and_joint(state, p)
    output_entropy, output_trace = _certified_entropy(output)
    joint_entropy, joint_trace = _certified_entropy(joint)
    coherent_lower = output_entropy[0] - joint_entropy[1]
    coherent_upper = output_entropy[1] - joint_entropy[0]
    classification = (
        "certified-positive"
        if coherent_lower > 0
        else "certified-negative"
        if coherent_upper < 0
        else "certified-straddles-zero"
    )
    return {
        "schema": CERTIFIED_SCHEMA,
        "candidate_id": candidate_id,
        "coefficient_payload_sha256": coefficient_payload_sha256,
        "family_id": family_id,
        "n": n,
        "p_per_pauli": p_per_pauli,
        "p_exact": _fraction_text(p),
        "coefficient_semantics": (
            "exact complex128 dyadics; common computational-basis amplitude "
            "within each Hamming-weight sector; exact checker normalization"
        ),
        "raw_binomial_normalization": _fraction_text(normalization),
        "coherent_information_bits": {
            "lower": _directed_decimal(coherent_lower, upper=False),
            "upper": _directed_decimal(coherent_upper, upper=True),
        },
        "classification": classification,
        "output_entropy_bits": {
            "lower": _directed_decimal(output_entropy[0], upper=False),
            "upper": _directed_decimal(output_entropy[1], upper=True),
        },
        "joint_entropy_bits": {
            "lower": _directed_decimal(joint_entropy[0], upper=False),
            "upper": _directed_decimal(joint_entropy[1], upper=True),
        },
        "output_spectrum_certificate": output_trace,
        "joint_spectrum_certificate": joint_trace,
        "arithmetic": {
            "input_coefficients": "exact-binary-rational-from-complex128-bytes",
            "noise_parameter": "exact-decimal-rational",
            "channel_and_accumulation": "exact-gaussian-rational",
            "hermitian_eigenvalues": "exact-polar-transform-plus-Gershgorin-components",
            "entropy": "exact-rational-atanh-log-enclosures",
            "block_weights": "not-applicable-direct-full-state-formulation",
            "root_bits": CERTIFIED_ROOT_BITS,
            "log_terms": CERTIFIED_LOG_TERMS,
        },
        "claim_boundary": (
            "candidate-specific coherent information at exact p in the frozen "
            "n=5 sector-symmetric family only; no n=45 or optimized-threshold claim"
        ),
    }


def validate_certified_evaluation(document: Mapping[str, Any]) -> None:
    if document.get("schema") != CERTIFIED_SCHEMA:
        raise CapacityCheckError("certified-evaluation-schema")
    interval = document.get("coherent_information_bits")
    if not isinstance(interval, Mapping):
        raise CapacityCheckError("certified-evaluation-interval")
    lower = Fraction(str(interval["lower"]))
    upper = Fraction(str(interval["upper"]))
    if lower > upper:
        raise CapacityCheckError("certified-evaluation-direction")
    expected = (
        "certified-positive"
        if lower > 0
        else "certified-negative"
        if upper < 0
        else "certified-straddles-zero"
    )
    if document.get("classification") != expected:
        raise CapacityCheckError("certified-evaluation-classification")


def _strict_json(path: Path) -> dict[str, Any]:
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in items:
            if key in value:
                raise CapacityCheckError("whole-return-duplicate-json-key")
            value[key] = item
        return value

    try:
        document = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=pairs,
            parse_constant=lambda value: (_ for _ in ()).throw(
                CapacityCheckError(f"whole-return-nonfinite-{value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise CapacityCheckError("whole-return-json-invalid") from exc
    if not isinstance(document, dict):
        raise CapacityCheckError("whole-return-json-not-object")
    return document


def certified_return_bindings(target_condition_sha256: str) -> dict[str, str]:
    return {
        "task_id": CERTIFIED_TASK_ID,
        "target_condition_sha256": target_condition_sha256,
        "family_id": CERTIFIED_FAMILY_ID,
        "run_id": CERTIFIED_RUN_ID,
        "condition_id": CERTIFIED_CONDITION_ID,
        "handoff_id": CERTIFIED_HANDOFF_ID,
        "checkpoint_id": CERTIFIED_CHECKPOINT_ID,
        "p_per_pauli": CERTIFIED_P,
    }


def validate_whole_return(
    return_directory: Path,
    *,
    target_condition_sha256: str,
) -> dict[str, Any]:
    """Recompute and bind an exact three-file participant return."""

    root = return_directory.resolve()
    expected_files = {"RETURN.json", "coefficients.c128", "CHECKER-RESULT.json"}
    actual_files = {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file()
    }
    if actual_files != expected_files:
        raise CapacityCheckError("whole-return-file-inventory")
    return_path = root / "RETURN.json"
    return_document = _strict_json(return_path)
    if return_path.read_bytes() != canonical_json_bytes(return_document):
        raise CapacityCheckError("return-json-noncanonical")
    expected_keys = {
        "schema",
        "task_id",
        "target_condition_sha256",
        "family_id",
        "run_id",
        "condition_id",
        "handoff_id",
        "checkpoint_id",
        "candidate_id",
        "coefficient_payload_path",
        "coefficient_payload_sha256",
        "p_per_pauli",
        "checker_result_path",
        "checker_result_sha256",
        "cost_telemetry",
        "verdict",
        "nonclaims",
    }
    if set(return_document) != expected_keys:
        raise CapacityCheckError("whole-return-schema-keys")
    if return_document["schema"] != WHOLE_RETURN_SCHEMA:
        raise CapacityCheckError("whole-return-schema")
    expected_bindings = certified_return_bindings(target_condition_sha256)
    actual_bindings = {key: return_document[key] for key in expected_bindings}
    validate_exact_bindings(actual_bindings, expected_bindings)
    if (
        return_document["coefficient_payload_path"] != "coefficients.c128"
        or return_document["checker_result_path"] != "CHECKER-RESULT.json"
    ):
        raise CapacityCheckError("whole-return-path-binding")
    coefficient_bytes = (root / "coefficients.c128").read_bytes()
    if len(coefficient_bytes) != 2 * 6 * 16:
        raise CapacityCheckError("whole-return-coefficient-shape")
    coefficient_sha = sha256_bytes(coefficient_bytes)
    if return_document["coefficient_payload_sha256"] != coefficient_sha:
        raise CapacityCheckError("whole-return-coefficient-hash")
    state = (
        np.frombuffer(coefficient_bytes, dtype="<c16")
        .astype(np.complex128, copy=True)
        .reshape(2, 6)
    )
    candidate_id = return_document["candidate_id"]
    if not isinstance(candidate_id, str) or not candidate_id:
        raise CapacityCheckError("whole-return-candidate-id")
    recomputed = certified_coherent_information(
        state,
        p_per_pauli=CERTIFIED_P,
        candidate_id=candidate_id,
        coefficient_payload_sha256=coefficient_sha,
        family_id=CERTIFIED_FAMILY_ID,
    )
    result_path = root / "CHECKER-RESULT.json"
    result_bytes = result_path.read_bytes()
    if return_document["checker_result_sha256"] != sha256_bytes(result_bytes):
        raise CapacityCheckError("whole-return-result-hash")
    supplied = _strict_json(result_path)
    if result_bytes != (
        json.dumps(
            supplied,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("ascii")
        + b"\n"
    ):
        raise CapacityCheckError("whole-return-result-not-canonical")
    if supplied != recomputed:
        raise CapacityCheckError("whole-return-result-recompute")
    validate_certified_evaluation(supplied)
    expected_verdict = (
        "CERTIFIED_POSITIVE_WITNESS"
        if recomputed["classification"] == "certified-positive"
        else "NO_CERTIFIED_WITNESS"
    )
    if return_document["verdict"] != expected_verdict:
        raise CapacityCheckError("whole-return-verdict")
    if return_document["nonclaims"] != RETURN_NONCLAIMS:
        raise CapacityCheckError("whole-return-nonclaims")
    if not isinstance(return_document["cost_telemetry"], Mapping):
        raise CapacityCheckError("whole-return-cost-telemetry")
    return {
        "schema": f"{WHOLE_RETURN_SCHEMA}/validation",
        "status": "VALID",
        "bindings": expected_bindings,
        "candidate_id": candidate_id,
        "coefficient_payload_sha256": coefficient_sha,
        "checker_result_sha256": sha256_bytes(result_bytes),
        "classification": recomputed["classification"],
        "verdict": expected_verdict,
        "recomputed": True,
        "exact_file_inventory": sorted(expected_files),
        "nonclaims": RETURN_NONCLAIMS,
    }
