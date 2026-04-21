from collections.abc import Iterable
from itertools import islice

from more_itertools import padded


def up_to_n[T, S](_iter: Iterable[T], n: int, fillvalue: S) -> Iterable[T | S]:
    return islice(padded(_iter, n=n, fillvalue=fillvalue), n)
