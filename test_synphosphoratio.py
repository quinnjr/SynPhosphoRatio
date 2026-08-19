from __future__ import annotations

import math
import sys
import types
from pathlib import Path

import pytest

PLUMA = Path(__file__).resolve().parent.parent / "PluMA"
sys.path.insert(0, str(PLUMA))

stub = types.ModuleType("PyPluMA")
stub._prefix = ""
stub.prefix = lambda: stub._prefix
sys.modules["PyPluMA"] = stub

from SynPhosphoRatioPlugin import SynPhosphoRatioPlugin, _fmt  # noqa: E402


def build_inputs(tmp_path: Path, total_rows: list[str], phospho_rows: list[str],
                  extra_params: str = "") -> Path:
    (tmp_path / "total.csv").write_text(
        "\n".join(['"","value"'] + total_rows) + "\n")
    (tmp_path / "phospho.csv").write_text(
        "\n".join(['"","value"'] + phospho_rows) + "\n")
    params = tmp_path / "params.txt"
    params.write_text(
        "total\ttotal.csv\nphospho\tphospho.csv\n" + extra_params
    )
    return params


def run(params: Path, tmp_path: Path) -> list[str]:
    stub._prefix = str(tmp_path)
    plugin = SynPhosphoRatioPlugin()
    plugin.input(str(params))
    plugin.run()
    out = tmp_path / "out.tsv"
    plugin.output(str(out))
    return out.read_text().splitlines()


def test_default_log_inputs_delogs_before_dividing(tmp_path: Path) -> None:
    # total log2 value 3.0 -> linear 2**3-1 = 7.0
    # phospho log2 value 2.0 -> linear 2**2-1 = 3.0
    # ratio = 3.0 / 7.0
    params = build_inputs(
        tmp_path,
        ['PD_001,3.0'],
        ['PD_001,2.0'],
    )
    lines = run(params, tmp_path)
    assert lines[0] == "sample\tps129_ratio\ttotal\tphospho"
    sample, ratio, total, phospho = lines[1].split("\t")
    assert sample == "PD_001"
    assert float(ratio) == pytest.approx(3.0 / 7.0)
    assert float(total) == pytest.approx(7.0)
    assert float(phospho) == pytest.approx(3.0)


def test_log_inputs_false_uses_raw_values(tmp_path: Path) -> None:
    params = build_inputs(
        tmp_path,
        ['PD_001,7.0'],
        ['PD_001,3.0'],
        extra_params="log_inputs\tfalse\n",
    )
    lines = run(params, tmp_path)
    sample, ratio, total, phospho = lines[1].split("\t")
    assert float(ratio) == pytest.approx(3.0 / 7.0)
    assert float(total) == pytest.approx(7.0)
    assert float(phospho) == pytest.approx(3.0)


def test_zero_guard_used_when_total_near_zero(tmp_path: Path) -> None:
    # log_inputs=false so total stays 0.0 (below default zero_guard 1e-6)
    params = build_inputs(
        tmp_path,
        ['PD_001,0.0'],
        ['PD_001,5.0'],
        extra_params="log_inputs\tfalse\n",
    )
    lines = run(params, tmp_path)
    sample, ratio, total, phospho = lines[1].split("\t")
    assert float(ratio) == pytest.approx(5.0 / 1e-6)
    assert float(total) == pytest.approx(0.0)


def test_custom_zero_guard_overrides_default(tmp_path: Path) -> None:
    params = build_inputs(
        tmp_path,
        ['PD_001,0.0'],
        ['PD_001,5.0'],
        extra_params="log_inputs\tfalse\nzero_guard\t0.5\n",
    )
    lines = run(params, tmp_path)
    _, ratio, _, _ = lines[1].split("\t")
    assert float(ratio) == pytest.approx(5.0 / 0.5)


def test_only_shared_samples_are_output(tmp_path: Path) -> None:
    params = build_inputs(
        tmp_path,
        ['PD_001,1.0', 'ONLY_TOTAL,2.0'],
        ['PD_001,1.0', 'ONLY_PHOSPHO,2.0'],
        extra_params="log_inputs\tfalse\n",
    )
    lines = run(params, tmp_path)
    samples = [line.split("\t")[0] for line in lines[1:]]
    assert samples == ["PD_001"]


def test_no_shared_samples_raises_value_error(tmp_path: Path) -> None:
    params = build_inputs(
        tmp_path,
        ['PD_001,1.0'],
        ['CTRL_001,1.0'],
    )
    with pytest.raises(ValueError, match="no shared samples"):
        run(params, tmp_path)


def test_quoted_sample_ids_are_stripped(tmp_path: Path) -> None:
    params = build_inputs(
        tmp_path,
        ['"PD_001",1.0'],
        ['"PD_001",1.0'],
        extra_params="log_inputs\tfalse\n",
    )
    lines = run(params, tmp_path)
    sample = lines[1].split("\t")[0]
    assert sample == "PD_001"


def test_fmt_formats_non_finite_values_as_na() -> None:
    assert _fmt(float("nan")) == "NA"
    assert _fmt(float("inf")) == "NA"
    assert _fmt(float("-inf")) == "NA"
    assert _fmt(1.5) == "1.5"


def test_fmt_matches_output_column_formatting(tmp_path: Path) -> None:
    params = build_inputs(
        tmp_path,
        ['PD_001,1.0'],
        ['PD_001,0.5'],
        extra_params="log_inputs\tfalse\n",
    )
    lines = run(params, tmp_path)
    _, ratio, total, phospho = lines[1].split("\t")
    assert ratio == _fmt(0.5 / 1.0)
    assert total == _fmt(1.0)
    assert phospho == _fmt(0.5)
