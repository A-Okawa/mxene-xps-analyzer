import pandas as pd
import pytest

from mxene_xps.io import normalize_peak_table, parse_two_column


def test_header_and_noheader_and_delimiters():
    body = "\n".join(f"{450 + i * 0.1:.2f}{d}{100 + i}" for i, d in zip(range(30), [",", "\t", ";", " "] * 10))
    s1 = parse_two_column("BE,intensity\n" + body)
    s2 = parse_two_column(body)
    assert s1.had_header and not s2.had_header and len(s1.be) == len(s2.be) == 30


def test_descending_sorted_and_nan_dropped():
    lines = [f"{480 - i * 0.5},{i}" for i in range(25)] + ["nan,1", "abc,def", "481,"]
    s = parse_two_column("\n".join(lines))
    assert s.was_descending and (s.be[1:] > s.be[:-1]).all() and s.n_dropped == 3


def test_too_few_points():
    with pytest.raises(ValueError, match="at least 20"):
        parse_two_column("1,2\n3,4\n")


def test_peak_table_required_columns():
    with pytest.raises(ValueError, match="missing required column"):
        normalize_peak_table(pd.DataFrame({"sample": ["a"], "component": ["x"], "area": [1]}))


def test_element_inferred():
    t = normalize_peak_table(pd.DataFrame({"sample": ["a"], "core_level": ["Nb 3d"], "component": ["x"],
                                           "area": [1.0], "role": ["framework_metal"]})).df
    assert t["element"].item() == "Nb"
