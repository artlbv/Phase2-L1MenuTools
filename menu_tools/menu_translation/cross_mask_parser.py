"""Parser for the `cross_masks` expressions of MenuTools rate-table menus.

Only a restricted grammar, corresponding to cuts the P2GT emulator can apply,
is accepted. Each mask is a conjunction (`&`) of comparisons with a constant:

    legI.deltaR(legJ)            -> PairCut "dR"
    (legI+legJ).mass             -> PairCut "mass"
    abs(legI.eta-legJ.eta)       -> PairCut "dEta"
    abs(legI.z0-legJ.z0)         -> PairCut "dz"
    legI.charge*legJ.charge      -> PairCut "chargeProduct"
    abs(legI.eta)                -> LegCut  "absEta"

Anything else raises a TranslationError naming the expression.
"""

import ast
import re
from dataclasses import dataclass
from typing import Union

from menu_tools.menu_translation.seed_model import TranslationError


@dataclass(frozen=True)
class PairCut:
    i: int
    j: int
    quantity: str
    op: str  # "<" or ">"
    value: float


@dataclass(frozen=True)
class LegCut:
    i: int
    quantity: str
    op: str
    value: float


Cut = Union[PairCut, LegCut]

_OPS = {ast.Lt: "<", ast.LtE: "<", ast.Gt: ">", ast.GtE: ">"}
_ABS_DIFF_QUANTITIES = {"eta": "dEta", "z0": "dz"}


def _leg_index(node: ast.expr) -> int:
    if isinstance(node, ast.Name):
        match = re.fullmatch(r"leg(\d+)", node.id)
        if match:
            return int(match.group(1))
    raise TranslationError(f"expected legN, got '{ast.unparse(node)}'")


def _leg_attribute(node: ast.expr) -> tuple[int, str]:
    if isinstance(node, ast.Attribute):
        return _leg_index(node.value), node.attr
    raise TranslationError(f"expected legN.<quantity>, got '{ast.unparse(node)}'")


def _constant(node: ast.expr) -> float:
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_constant(node.operand)
    raise TranslationError(f"expected a number, got '{ast.unparse(node)}'")


def _quantity(node: ast.expr) -> tuple[str, tuple[int, ...]]:
    """(quantity, leg indices) of the left-hand side of a comparison"""
    # legI.deltaR(legJ)
    if (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "deltaR"
        and len(node.args) == 1
    ):
        return "dR", (_leg_index(node.func.value), _leg_index(node.args[0]))
    # (legI+legJ).mass
    if (
        isinstance(node, ast.Attribute)
        and node.attr == "mass"
        and isinstance(node.value, ast.BinOp)
        and isinstance(node.value.op, ast.Add)
    ):
        return "mass", (_leg_index(node.value.left), _leg_index(node.value.right))
    # legI.charge*legJ.charge
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        (i, qi), (j, qj) = _leg_attribute(node.left), _leg_attribute(node.right)
        if qi == qj == "charge":
            return "chargeProduct", (i, j)
    # abs(...)
    if (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "abs"
        and len(node.args) == 1
    ):
        arg = node.args[0]
        if isinstance(arg, ast.BinOp) and isinstance(arg.op, ast.Sub):
            (i, qi), (j, qj) = _leg_attribute(arg.left), _leg_attribute(arg.right)
            if qi == qj and qi in _ABS_DIFF_QUANTITIES:
                return _ABS_DIFF_QUANTITIES[qi], (i, j)
        elif isinstance(arg, ast.Attribute):
            i, q = _leg_attribute(arg)
            if q == "eta":
                return "absEta", (i,)
    raise TranslationError(f"unsupported quantity '{ast.unparse(node)}'")


def _conjunction_terms(node: ast.expr) -> list[ast.expr]:
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitAnd):
        return _conjunction_terms(node.left) + _conjunction_terms(node.right)
    if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And):
        return [t for value in node.values for t in _conjunction_terms(value)]
    return [node]


def parse_cross_mask(expression: str) -> tuple[list[Cut], list[str]]:
    """Cuts of one `cross_masks` entry, plus warnings."""
    try:
        tree = ast.parse(expression.strip(), mode="eval")
    except SyntaxError as e:
        raise TranslationError(f"cannot parse '{expression}': {e}")

    cuts: list[Cut] = []
    warnings: list[str] = []
    terms = _conjunction_terms(tree.body)
    while terms:
        term = terms.pop(0)
        if not (isinstance(term, ast.Compare) and len(term.ops) == 1):
            raise TranslationError(f"unsupported expression '{ast.unparse(term)}'")
        right = term.comparators[0]
        # `x < 1 & (y > 0)` is parsed by python as `x < (1 & (y > 0))`
        if isinstance(right, ast.BinOp) and isinstance(right.op, ast.BitAnd):
            warnings.append(
                f"'{expression}': '&' binds tighter than '<'/'>' in python, "
                f"interpreted as '{ast.unparse(term.left)} "
                f"{_OPS[type(term.ops[0])]} {ast.unparse(right.left)}' "
                f"& '{ast.unparse(right.right)}' - please add parentheses"
            )
            terms = _conjunction_terms(right.right) + terms
            right = right.left
        op_type = type(term.ops[0])
        if op_type not in _OPS:
            raise TranslationError(f"unsupported comparison in '{ast.unparse(term)}'")
        op = _OPS[op_type]
        quantity, legs = _quantity(term.left)
        value = _constant(right)
        if len(legs) == 1:
            cuts.append(LegCut(legs[0], quantity, op, value))
        else:
            cuts.append(PairCut(legs[0], legs[1], quantity, op, value))
    return cuts, warnings
