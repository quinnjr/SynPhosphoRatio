import math
import os

import PyIO
import PyPluMA


class SynPhosphoRatioPlugin:

    def input(self, filename):
        self.parameters = PyIO.readParameters(filename)
        prefix = PyPluMA.prefix()
        self.total_path = os.path.join(prefix, self.parameters["total"])
        self.phospho_path = os.path.join(prefix, self.parameters["phospho"])
        self.zero_guard = float(self.parameters.get("zero_guard", "1e-6"))
        self.log_inputs = self.parameters.get("log_inputs", "true").strip().lower() == "true"

    def run(self):
        total = _read_single_feature_csv(self.total_path)
        phospho = _read_single_feature_csv(self.phospho_path)
        shared = [s for s in total if s in phospho]
        if not shared:
            raise ValueError(
                "SynPhosphoRatio: no shared samples between "
                + self.total_path + " and " + self.phospho_path
            )
        self.rows = []
        for sample in shared:
            t = total[sample]
            p = phospho[sample]
            if self.log_inputs:
                t_lin = max(2.0 ** t - 1.0, 0.0)
                p_lin = max(2.0 ** p - 1.0, 0.0)
            else:
                t_lin, p_lin = t, p
            denom = t_lin if t_lin > self.zero_guard else self.zero_guard
            ratio = p_lin / denom
            self.rows.append((sample, ratio, t_lin, p_lin))

    def output(self, filename):
        with open(filename, "w") as out:
            out.write("sample\tps129_ratio\ttotal\tphospho\n")
            for sample, ratio, t, p in self.rows:
                out.write(
                    sample + "\t" + _fmt(ratio) + "\t" + _fmt(t) + "\t" + _fmt(p) + "\n"
                )


def _read_single_feature_csv(path):
    values = {}
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            parts = [p.strip() for p in line.rstrip("\n").split(",")]
            if len(parts) < 2:
                continue
            sample = parts[0].strip('"')
            try:
                values[sample] = float(parts[1])
            except ValueError:
                continue
    return values


def _fmt(v):
    if math.isnan(v) or math.isinf(v):
        return "NA"
    return repr(v)
