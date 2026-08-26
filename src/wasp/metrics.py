"""Legacy Xarray verification metrics for significant wave height fields."""

import numpy as np
import xarray as xr

def compute_mad(sim, obs):
    """Return mean absolute difference along the ``time`` dimension."""
    return np.abs(sim - obs).mean(dim="time")

def compute_madp(sim, obs):
    """Return mean absolute separation between empirical percentile curves."""
    percentiles = np.arange(101)
    sim_p = np.nanpercentile(sim, percentiles, axis=0)
    obs_p = np.nanpercentile(obs, percentiles, axis=0)
    return np.nanmean(np.abs(sim_p - obs_p), axis=0)

def metrics(data):
    """Compute the legacy SWH metric dataset.

    ``data`` must provide Xarray variables ``SWH_mod`` and ``SWH_sat`` with a
    ``time`` dimension. The returned dataset contains bias, observation count,
    RMSE, means, normalized bias/RMSE, MAD, MADP, and MADC.
    """
    bias = data['SWH_mod'] - data['SWH_sat']
    obs = data['SWH_sat']
    obs_mean = obs.mean(dim='time')

    # Criação do Dataset de saída
    result = xr.Dataset()

    # Bias
    bias_da = bias.mean(dim='time')
    bias_da.name = "bias"
    result['bias'] = bias_da

    # Número de observações
    nobs_da = bias.count(dim='time')
    nobs_da.name = "nobs"
    result['nobs'] = nobs_da

    # RMSE
    rmse_da = np.sqrt((bias**2).mean(dim='time'))
    rmse_da.name = "rmse"
    result['rmse'] = rmse_da

    # Médias
    model_hs_da = data['SWH_mod'].mean(dim='time')
    model_hs_da.name = "model_hs"
    result['model_hs'] = model_hs_da

    sat_hs_da = obs_mean
    sat_hs_da.name = "sat_hs"
    result['sat_hs'] = sat_hs_da

    # Normalizações em %
    nbias_da = 100 * (bias_da / obs_mean)
    nbias_da.name = "nbias"
    result['nbias'] = nbias_da

    nrmse_da = 100 * (rmse_da / obs_mean)
    nrmse_da.name = "nrmse"
    result['nrmse'] = nrmse_da

    # MAD
    mad_da = compute_mad(data['SWH_mod'], data['SWH_sat'])
    if not isinstance(mad_da, xr.DataArray):
        mad_da = xr.DataArray(mad_da, coords=bias_da.coords, dims=bias_da.dims, name="mad")
    mad_da.name = "mad"
    result['mad'] = mad_da

    # MADP
    madp_arr = compute_madp(data['SWH_mod'], data['SWH_sat'])
    madp_da = xr.DataArray(
        madp_arr,
        coords={k: v for k, v in data['SWH_mod'].coords.items() if k != 'time'},
        dims=[d for d in data['SWH_mod'].dims if d != 'time'],
        name="madp"
    )
    result['madp'] = madp_da

    # MADC
    madc_da = mad_da + madp_da
    madc_da.name = "madc"
    result['madc'] = madc_da

    return result

