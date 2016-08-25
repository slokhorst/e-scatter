import numbers
import collections
from . import (units)


class Predicate:
    def __init__(self, f):
        self.f = f

    def __call__(self, v):
        return self.f(v)

    def __and__(self, other):
        return Predicate(lambda v: self(v) and other(v))

    def __or__(self, other):
        return Predicate(lambda v: self(v) or other(v))


@Predicate
def is_integer(v):
    return isinstance(v, int)


@Predicate
def is_string(v):
    return isinstance(v, str)


@Predicate
def is_number(v):
    return isinstance(v, numbers.Number)


@Predicate
def is_none(v):
    return v is None


def has_units_like(u):
    @Predicate
    def _has_units(v):
        return isinstance(v, units.Quantity) and \
            v.dimensionality == u.dimensionality

    return _has_units


is_energy = has_units_like(units.J)
is_length = has_units_like(units.m)
is_volume = has_units_like(units.m**3)


def in_range(a, b):
    @Predicate
    def _in_range(v):
        return v >= a and v < b

    return _in_range


def is_(a):
    @Predicate
    def _is_(v):
        return v == a

    return _is_


@Predicate
def is_iterable(v):
    return isinstance(v, collections.Iterable)


def is_list_of(p):
    @Predicate
    def _is_list_of(v):
        if not isinstance(v, list):
            return False

        for e in v:
            if not p(e):
                return False

        return True

    return _is_list_of