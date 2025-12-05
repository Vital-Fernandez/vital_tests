import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import GaussianModel, LinearModel
from lmfit import Model, fit_report

def make_array_model(cont_array, prefix='cont_', allow_scale=False):
    """
    Build a Model that returns the supplied continuum array.
    If allow_scale=True, include a 'scale' parameter (free by default).
    If allow_scale=False, include 'scale' but freeze it at 1.0 so it's constant.
    """
    cont_array = np.asarray(cont_array)

    def array_component(x, scale=1.0):
        # x is ignored for values, but we check shape for safety
        if np.shape(x) != np.shape(cont_array):
            raise ValueError("continuum array and x must have the same shape")
        return scale * cont_array

    m = Model(array_component, prefix=prefix)
    params = m.make_params(scale=1.0)
    if not allow_scale:
        params[f'{prefix}scale'].set(vary=False)  # freeze at 1.0 -> truly constant
    return m, params

def fit_gaussian_with_continuum(
    x, y,
    mode="linear",                 # "linear" | "array-fixed" | "array-scaled"
    continuum_array=None
):
    """
    Fit a Gaussian plus a continuum, where the continuum can be:
      - 'linear'        : LinearModel (slope + intercept are fitted)
      - 'array-fixed'   : supplied array added as a fixed (constant) component
      - 'array-scaled'  : supplied array with a fitted scale factor

    Returns
    -------
    result : lmfit.ModelResult
    components : dict of component evaluations
    """
    x = np.asarray(x)
    y = np.asarray(y)

    gmod = GaussianModel(prefix='g_')
    params = gmod.guess(y, x=x)

    if mode == "linear":
        lmod = LinearModel(prefix='lin_')
        model = gmod + lmod
        params.update(lmod.make_params(lin_intercept=np.median(y), lin_slope=0.0))

    elif mode in ("array-fixed", "array-scaled"):
        if continuum_array is None:
            raise ValueError("Provide continuum_array for array-fixed/array-scaled modes.")
        allow_scale = (mode == "array-scaled")
        amod, aparams = make_array_model(continuum_array, prefix='cont_', allow_scale=allow_scale)
        model = gmod + amod
        params.update(aparams)

    else:
        raise ValueError("mode must be 'linear', 'array-fixed', or 'array-scaled'.")

    result = model.fit(y, params, x=x)
    components = result.eval_components(x=x)

    # For consistency, construct a total best_fit too (lmfit already sets this)
    # but in array modes you may want to access the continuum component explicitly:
    #   components['cont_'] is the fixed/scaled array contribution
    return result, components

# --- Example usage ---
if __name__ == "__main__":
    rng = np.random.default_rng(1)
    x = np.linspace(-5, 5, 401)
    continuum_true = 0.25*x + 2.0
    gauss_true = (7.5/(0.6*np.sqrt(2*np.pi))) * np.exp(-(x-0.7)**2/(2*0.6**2))
    y_true = continuum_true + gauss_true
    y = y_true + rng.normal(0, 0.15, size=x.size)

    # A) Linear continuum fit
    res_lin, comps_lin = fit_gaussian_with_continuum(x, y, mode="linear")

    # B) Fixed array continuum (added as a constant component)
    res_fix, comps_fix = fit_gaussian_with_continuum(
        x, y, mode="array-fixed", continuum_array=continuum_true
    )

    # C) Scaled array continuum (same shape but with free scale)
    res_scl, comps_scl = fit_gaussian_with_continuum(
        x, y, mode="array-scaled", continuum_array=continuum_true * 0.9
    )

    print("== Linear continuum ==\n", fit_report(res_lin))
    print("== Array (fixed) ==\n", fit_report(res_fix))
    print("== Array (scaled) ==\n", fit_report(res_scl))

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
    for ax, res, comps, title in [
        (axes[0], res_lin, comps_lin, "LinearModel"),
        (axes[1], res_fix, comps_fix, "Array (fixed)"),
        (axes[2], res_scl, comps_scl, "Array (scaled)"),
    ]:
        ax.plot(x, y, '.', ms=3, label='data')
        ax.plot(x, res.best_fit, '-', lw=2, label='best-fit')
        ax.plot(x, comps['g_'], '--', lw=2, label='Gaussian')
        # component key for the array model is the prefix 'cont_'
        if 'lin_' in comps:
            ax.plot(x, comps['lin_'], ':', lw=2, label='Linear continuum')
        if 'cont_' in comps:
            ax.plot(x, comps['cont_'], ':', lw=2, label='Array continuum')
        ax.set_title(title)
        ax.legend()

    plt.tight_layout()
    plt.show()
