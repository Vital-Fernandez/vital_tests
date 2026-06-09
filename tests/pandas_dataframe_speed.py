import pandas as pd
import numpy as np
import timeit

# ── Config ────────────────────────────────────────────────────────────────────
N_ROWS   = 100_000
N_REPEAT = 5       # number of timing runs per method (best-of is reported)
N_ITER   = 3       # timeit iterations per run

# ── Setup ─────────────────────────────────────────────────────────────────────
rng = np.random.default_rng(42)
df  = pd.DataFrame({'A': rng.random(N_ROWS), 'B': rng.random(N_ROWS)})

methods = {
    "Vectorization (pandas)": lambda: df['A'] + df['B'],

    "NumPy (.to_numpy())": lambda: (
        df['A'].to_numpy() + df['B'].to_numpy()
    ),

    "itertuples()": lambda: [
        row.A + row.B for row in df.itertuples()
    ],

    "apply(axis=1)": lambda: df.apply(
        lambda row: row['A'] + row['B'], axis=1
    ),

    "iterrows()": lambda: [
        row['A'] + row['B'] for _, row in df.iterrows()
    ],

    # .loc loop kept but isolated — it's dramatically slower
    "Index loop (.loc)": lambda: [
        df.loc[i, 'A'] + df.loc[i, 'B'] for i in df.index
    ],
}

# ── Benchmark ─────────────────────────────────────────────────────────────────
results = {}
for name, fn in methods.items():
    # timeit.repeat runs fn N_ITER times, repeated N_REPEAT times
    times   = timeit.repeat(fn, number=N_ITER, repeat=N_REPEAT)
    best    = min(times) / N_ITER          # best single-run time
    results[name] = best

# ── Report ────────────────────────────────────────────────────────────────────
sorted_results = sorted(results.items(), key=lambda x: x[1])
fastest = sorted_results[0][1]

print(f"\n{'Method':<25} | {'Time (s)':>10} | {'Rel. speed':>12}")
print("-" * 55)

slow_methods = {}
threshold = fastest * 500   # separate methods >500× slower

for method, t in sorted_results:
    if t > threshold:
        slow_methods[method] = t
        continue
    rel = t / fastest
    bar = "█" * max(1, int(rel * 3))
    print(f"{method:<25} | {t:>10.5f} | {rel:>10.1f}×  {bar}")

if slow_methods:
    print("\n── Slow methods (shown separately to preserve scale) ──")
    for method, t in slow_methods.items():
        rel = t / fastest
        print(f"  {method:<23} | {t:>10.5f} | {rel:>10.1f}×")

print(f"\n(best of {N_REPEAT} runs × {N_ITER} iterations, {N_ROWS:,} rows)\n")