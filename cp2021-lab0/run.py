from dataclasses import dataclass
import itertools as it
import os
import shlex
import subprocess

import pandas as pd

# =============================================================================
# GRID CLASS
# =============================================================================


@dataclass
class Grid:
    params_grid: dict

    def __iter__(self):
        items = sorted(self.params_grid.items())
        if not items:
            yield {}
        else:
            keys, values = zip(*items)
            for v in it.product(*values):
                params = dict(zip(keys, v))
                yield params


# =============================================================================
# CONF
# =============================================================================

COMPILE_TPL = "{CC} -O{O} -fopenmp  -Wall -Wextra -o _out/matmul_{CC}_{O} matmul.c"

RUN_TPL = "_out/matmul_{CC}_{O}"


PARAM_GRID = {"CC": ["gcc", "clang"], "O": [0, 1, 2, 3]}


# =============================================================================
# RUN
# =============================================================================


def parse_out(out, params):
    rows = []
    for line in out.splitlines():
        implementation, gflops = line.replace("GFLOPS", "").split(" run:", 1)
        row = {"Implementation": implementation, "GFLOPS": float(gflops)}
        row.update(params)
        rows.append(row)
    return rows


def run():
    grid = Grid(PARAM_GRID)
    rows = []
    for params in grid:
        compile_cmd = COMPILE_TPL.format(**params)
        print("[COMPILING]", compile_cmd)
        args = shlex.split(compile_cmd)
        subprocess.run(args)

        run_cmd = RUN_TPL.format(**params)
        print("[RUNNING]", run_cmd)
        cmd = subprocess.run(run_cmd, capture_output=True)
        out = cmd.stdout.decode()

        srows = parse_out(out, params)
        rows.extend(srows)
        print("-" * 80)

    print("[WRITING] resume.csv")
    df = pd.DataFrame(rows)
    df.to_csv("resume.csv", index=False)


def load_df():
    if not os.path.exists("resume.csv"):
        run()
    return pd.read_csv("resume.csv")


if __name__ == "__main__":
    run()
