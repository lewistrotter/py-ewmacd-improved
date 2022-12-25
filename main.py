import os
import time
import multiprocessing

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from numpy.lib.stride_tricks import sliding_window_view as slide

from typing import Union
from functools import partial


def _harmonic_matrix(
        timeseries_0_to_2pi: np.ndarray,  # todo what is this? is it time to radians?
        num_harmonics_sin: int,
        num_harmonics_cos: int
) -> np.ndarray:
    """
    Creates a harmonic (gram) matrix for time series values using n harmonics
    (e.g., 2 sin, 2 cos). Matrix where t = time, K = num harmonics (e.g., cos, sin = 2) and
    M = size of vector:
        X = | 1   | sin(t1) | cos(t2) | ... | sin(Kt1) | cos(Kt1) |
            | ... | ...     | ...     | ... | ...      | ...      |
            | 1   | sin(tM) | cos(tM) | ... | sin(KtM) | cos(KtM) |

    :param timeseries_0_to_2pi:
    :param num_harmonics_sin:
    :param num_harmonics_cos:
    :return:
    """

    # create array (a column) of constants (1)
    con = np.full(len(timeseries_0_to_2pi), 1)

    # create sin column(s), one for each harmonic, todo storing coefficients
    sin = np.tile([[i + 1] for i in np.arange(num_harmonics_sin)], len(timeseries_0_to_2pi))
    sin = np.sin(sin * timeseries_0_to_2pi).T

    # create cos column(s), one for each harmonic, todo storing coefficients
    cos = np.tile([[i + 1] for i in np.arange(num_harmonics_cos)], len(timeseries_0_to_2pi))
    cos = np.cos(cos * timeseries_0_to_2pi).T

    # combine all columns together
    X = np.column_stack([con, sin, cos])  # todo this should be 1, sin, cos, sin, cos etc
    return X


def _harmonic_regression(
        responses: np.ndarray,
        doys: np.ndarray,
        num_harmonics_sin: int,
        num_harmonics_cos: int,
        anomaly_threshold: float = 1.5,
        values_cleaned: bool = True
) -> tuple:
    """

    :param responses:
    :param doys:
    :param num_harmonics_sin:
    :param num_harmonics_cos:
    :param anomaly_threshold:
    :param values_cleaned:
    :return:
    """

    # todo do we ever not use cleaned values?
    if values_cleaned is False:
        missing_idx = np.flatnonzero(np.isnan(responses))

        if len(missing_idx) > 0:
            responses = np.delete(responses, missing_idx)  # can we just subset instead of flatnonzero?
            doys = np.delete(doys, missing_idx)

    # init regression parameters (assumed things are clean up to here!)
    beta = np.full((1 + num_harmonics_sin + num_harmonics_cos), np.nan)
    r2 = None
    rmse = None

    # scale doys by converting to doys range 0 to 2pi (radians)
    timeseries_0_to_2pi = doys * 2 * np.pi / 365

    # generate harmonic matrix using n harmonics
    X = _harmonic_matrix(
        timeseries_0_to_2pi=timeseries_0_to_2pi,
        num_harmonics_sin=num_harmonics_sin,
        num_harmonics_cos=num_harmonics_cos
    )

    # if design matrix is sufficient rank and non-singular...
    if len(responses) > (1 + num_harmonics_sin + num_harmonics_cos):

        # todo add comment
        det = np.linalg.det(X.T.dot(X))
        if np.abs(det) >= 0.001:

            # apply initial harmonic regression algorithm
            a = X.T.dot(X)                       # calc coeff matrix
            v = X.T.dot(responses)               # calc ordinate values (via responses)
            fits = X.dot(np.linalg.solve(a, v))  # perform least squares method

            # do shewhart x-bar anomaly (cloud) filtering
            residuals = responses - fits                               # calc residuals (real - fits)
            stdv = residuals.std(ddof=1)                               # calc stdv (ddof 1 to match R)
            screen = (np.abs(residuals) > (anomaly_threshold * stdv))  # get residual outliers beyond threshold
            keeps = np.flatnonzero(screen == 0)                        # isolate non-outlier idxs

            # with outlier residuals removed, recompute harmonic coeffs
            if len(keeps) > (1 + num_harmonics_sin + num_harmonics_cos):

                # subset only clean time series and harmonic matrix
                keeps_X = X[keeps]
                keeps_responses = responses[keeps]

                # re-compute refined harmonic coefficients via least squares
                a = keeps_X.T.dot(keeps_X)          # calc coeff matrix again
                v = keeps_X.T.dot(keeps_responses)  # calc ordinate values (via responses)
                beta = np.linalg.solve(a, v)        # least squares
                fits = np.dot(keeps_X, beta)        # perform least squares method

                # compute r-squared (ssr: sum squared regression, sst: total sum squares)
                ssr = np.square(keeps_responses - fits).sum()
                sst = np.square(keeps_responses - keeps_responses.sum() / len(keeps_responses)).sum()
                r2 = 1 - (ssr / sst)

                # compute rmse
                rmse = np.sum(np.square(keeps_responses - fits))

    return beta, r2, rmse


def _optimise_harmonic_regression(
        decimal_years_cleaned: np.ndarray,
        doys_cleaned: np.ndarray,
        arr_cleaned: np.ndarray,
        train_fit_min_quality: float,
        min_train_length: int,
        max_train_length: int,
        num_harmonics_sin: int,
        num_harmonics_cos: int,
        xbar_limit_1: float
) -> np.ndarray:
    """"""

    # todo handle out of bounds for when < 0 veg is thresholded improve efficiancy
    # todo might wanna move this out to parent and do this there?
    #if len(arr_cleaned) <= min_train_length:
        #raise ValueError('Not enough clean training data.')

    # get minimum index where date >= user-specified minimum train length and > one full year of data
    rule_one = decimal_years_cleaned >= decimal_years_cleaned[min_train_length - 1]  # length is not an index, so - 1
    rule_two = decimal_years_cleaned - decimal_years_cleaned[0] > 1.0
    min_history_idx = np.min(np.flatnonzero(rule_one & rule_two))

    # todo finish this
    #if np.isinf(min_history_bound):
        #min_history_bound = 1

    # get maximum training history bound index, -1 to convert to index
    max_history_idx = np.min([len(arr_cleaned) - min_history_idx, max_train_length])

    # todo finish this
    #if np.isinf(np.max(history_bound_candidates)):
        #history_bound_candidates = len(decimal_years_cleaned)

    # build array of training windows (i.e. moving window)
    arr_wins = slide(arr_cleaned, min_train_length)[0:max_history_idx]
    doy_wins = slide(doys_cleaned, min_train_length)[0:max_history_idx]

    history_idx, fit_quality = 0, 0.0
    for history_idx, arrs in enumerate(zip(arr_wins, doy_wins), start=min_history_idx):

        # extract values and doys arrays
        vals, doys = arrs

        # remove nans if exist
        test_responses = vals[np.flatnonzero(~np.isnan(vals))]

        # perform harmonic regression
        result = _harmonic_regression(
            responses=test_responses,
            doys=doys,
            num_harmonics_sin=1,  # todo why 1 hardcode?
            num_harmonics_cos=1,  # todo why 1 hardcode?
            anomaly_threshold=xbar_limit_1,
            values_cleaned=True   # todo do we need this?
        )

        # extract r-squared from result (tuple is beta, r-squared, rmse)
        fit_quality = result[1]

        # set to 0 if None
        if fit_quality is None:
            fit_quality = 0.0

        # break out if we got adequate fit
        if fit_quality >= train_fit_min_quality:
            break

    return history_idx, min_history_idx


def _ewma_chart(
        arr: np.ndarray,
        stdv_history: float,
        lambda_value: float,
        lambda_bounds: float,
        rounding: bool
) -> np.ndarray:
    """"""

    # todo comment
    ewma = np.full(len(arr), np.nan)

    # init ewma outputs with the first present residual
    ewma[0] = arr[0]

    # append new ewma values for all remaining values
    for i in range(1, len(arr)):
        ewma[i] = ewma[i - 1] * (1 - lambda_value) + lambda_value * arr[i]

    # ewma upper control limit:  this is the threshold which dictates when the chart signals a disturbance # todo clean this up
    # decreasing lambda_value makes chart bars "plateau" later, increasing makes them plateau sooner
    ucl = stdv_history * lambda_bounds * np.sqrt(lambda_value / (2 - lambda_value) * (1 - np.power(1 - lambda_value, 2 * np.arange(1, len(arr) + 1))))

    # todo comment
    if rounding is True:
        # integer value for ewma output relative to control limit (rounded towards 0).
        # a value of +/-1 represents the weakest disturbance signal
        # set int16 to remove -0s
        output = (np.sign(ewma) * np.floor(np.abs(ewma / ucl))).astype('int16')  # todo prolly dont need int16 here
    else:
        # ewma outputs in terms of residual scales
        output = np.round(ewma, 0)

    return output


def _persistence_filter(
        arr: np.ndarray,
        persistence: float
) -> np.ndarray:
    """Use persistence as the threshold for keeping only values for which
    disturbance is sustained."""

    # todo comment
    # todo change var names, bit obtuse
    tmp4 = np.full(len(arr), 0)

    # ensure sufficient data for tmp2 todo what is tmp1, tmp2, tmp3, tmp4 etc
    if persistence > 1 and len(arr) > persistence:

        # get disturbance direction
        tmp_sign = np.sign(arr)

        # get dates of which direction changes  # todo modified this a lot, check it works
        diffs = np.flatnonzero(np.diff(tmp_sign))
        shift_points = np.concatenate([[0], diffs, [len(tmp_sign) - 1]])

        # counting the consecutive dates in which directions are sustained todo comment
        tmp3 = np.full(len(tmp_sign), 0)

        # todo comment
        for i in range(len(tmp_sign)):
            tmp3_lo, tmp3_hi = 0, 0

            # todo comment
            while (i + 1) - tmp3_lo > 0:  # todo added + 1
                if tmp_sign[i] - tmp_sign[i - tmp3_lo] == 0:
                    tmp3_lo += 1
                else:
                    break

            # todo comment
            while tmp3_hi + (i + 1) <= len(tmp_sign):  # todo added + 1
                if tmp_sign[i + tmp3_hi] - tmp_sign[i] == 0:
                    tmp3_hi += 1
                else:
                    break

            # todo comment
            tmp3[i] = tmp3_lo + tmp3_hi - 1

        # todo comment
        tmp4 = np.full(len(tmp3), 0)
        tmp3[0:int(persistence)] = persistence   # todo set persistence to int when defining first time, set inputs all to int
        arr[0:int(persistence)] = 0.0            # todo set persistence to int when defining first time, set inputs all to int

        # if sustained dates are long enough, keep; otherwise set to previous sustained state
        for i in range(int(persistence), len(tmp3)):  # todo set persistence to int when defining first time, set inputs all to int
            if tmp3[i] < persistence and tmp3[0:i].max() >= persistence:
                idx = np.flatnonzero(tmp3[0:i + 1] >= persistence).max()  # todo no longer making matrix, just using two arrs. check fine
                tmp4[i] = arr[0:i + 1][idx]  # todo good check here as well
            else:
                tmp4[i] = arr[i]

        return tmp4


def _backfill(
        non_missing: np.ndarray,
        non_missing_idx: np.ndarray,
        with_missing: np.ndarray  # todo can these var names be clearer?
) -> np.ndarray:
    """backfill missing data"""

    # todo comment
    with_missing[non_missing_idx] = non_missing

    # if the first date of arr was missing/filtered,
    # then assign the EWMA output as 0 (no disturbance).
    if np.isnan(with_missing[0]):
        with_missing[0] = 0

    # if we have ewma information for the first date,
    # then for each missing/filtered date in the record,
    # fill with the last known ewma value
    for step in range(1, len(with_missing)):  # todo exclude first element
        if np.isnan(with_missing[step]):
            with_missing[step] = with_missing[step - 1]

    return with_missing


def _clean_data(
        years: np.ndarray,
        doys: np.ndarray,
        arr: np.ndarray,
        history_idx: np.ndarray,
        precedents: np.ndarray,
        num_harmonics_sin: int,
        num_harmonics_cos: int,
        lambda_value: float,           # was lambda
        lambda_bounds: float,          # was lambdaSigs
        xbar_limit_1: float,
        xbar_limit_2: int,
        #lambda_value: float,   # was lambda
        #lambda_bounds: float,  # was lambdaSigs
        rounding: bool,
        persistence: float
) -> np.ndarray:
    """"""

    # todo ?
    num_doys = len(doys)
    out_values = np.full(len(doys), np.nan)
    residual_out_values = np.full(len(doys), np.nan)
    beta = np.full((1 + num_harmonics_sin + num_harmonics_cos), np.nan)  # todo this needs to be a column

    # extract training idxs, values and testing values
    train_idxs = np.arange(0, history_idx + 1)
    train_values = arr[train_idxs]
    test_values = np.delete(arr, train_idxs)

    if len(train_values) > 0:

        # perform harmonic regression
        result = _harmonic_regression(
            responses=train_values[(history_idx - precedents):history_idx + 1],  # todo added last + 1, check this
            doys=doys[train_idxs][(history_idx - precedents):history_idx + 1],   # todo added last + 1, check this
            num_harmonics_sin=num_harmonics_sin,
            num_harmonics_cos=num_harmonics_cos,
            anomaly_threshold=xbar_limit_1
        )

        # extract beta from harmonic regression result
        beta = result[0]

        if not np.isnan(beta).all():

            # scale doys by converting to doys range 0 to 2pi (radians)
            timeseries_0_to_2pi = doys * 2 * np.pi / 365

            # generate harmonic matrix using n harmonics
            X_all = _harmonic_matrix(
                timeseries_0_to_2pi=timeseries_0_to_2pi,
                num_harmonics_sin=num_harmonics_sin,
                num_harmonics_cos=num_harmonics_cos
            )

            # extract residuals for all present data based on training coefficients
            residuals = arr - X_all.dot(beta)
            residual_out_values = residuals.copy()

            # extract only training residuals out
            train_residuals = residuals[train_idxs]

            # set up testing residuals
            test_residuals = np.array([], dtype='float64')
            if len(residuals) > len(train_residuals):
                test_residuals = np.delete(residuals, train_idxs)

            # set first estimate of historical standard deviation
            train_stdv = train_residuals.std(ddof=1)

            # set up training residual indexes
            residuals_idxs = np.arange(0, len(residuals))
            residual_train_idxs = residuals_idxs[train_idxs]

            # set up testing residual indexes
            residual_test_idxs = np.array([], dtype='float64')
            if len(residuals_idxs) > len(residual_train_idxs):
                residual_test_idxs = np.delete(residuals_idxs, train_idxs)

            # modify standard deviation estimates based on anomalous readings in training data
            # ucl = upper control limit
            # note that we don't want to filter out the changes in the testing data, so xBarLimit2 is much larger!
            ucl = np.concatenate([np.full(len(residual_train_idxs), xbar_limit_1),
                                  np.full(len(residual_test_idxs), xbar_limit_2)]) * train_stdv

            # keep only indexes, residuals for which we have some vegetation and aren't anomalously far from 0 in the residuals
            clean_idxs = residuals_idxs[np.abs(residuals) < ucl]
            clean_residuals = residuals[clean_idxs]

            # update the training stdv estimate; this is the all-important modifier for the EWMA control limits
            clean_train_stdv = (train_residuals[np.abs(train_residuals) < ucl[residual_train_idxs]]).std(ddof=1)

            # if np.isnan(clean_train_stdv):
            # todo this check

            # perform ewma
            chart_output = _ewma_chart(
                arr=clean_residuals,
                stdv_history=clean_train_stdv,
                lambda_value=lambda_value,
                lambda_bounds=lambda_bounds,
                rounding=rounding
            )

            # keeping only values for which a disturbance is sustained, using persistence as the threshold
            # todo this needs a quick check, changed some things from r
            persistence_output = _persistence_filter(
                arr=chart_output,
                persistence=persistence
            )

            # imputing for missing values screened out as anomalous at the control limit stage
            empty_arr = np.full(len(residuals), np.nan)
            out_arr = _backfill(
                non_missing=persistence_output,
                non_missing_idx=clean_idxs,
                with_missing=empty_arr
            )

    return out_arr, residual_out_values, beta


def _ewmacd_core_func(
        arr:                   np.ndarray,
        dates:                 np.ndarray,
        train_start_year:      Union[int, None],
        #train_end_year:       Union[int, None],   # todo disabled for now
        test_end_year:         Union[int, None],
        min_train_length:      Union[int, None],
        max_train_length:      Union[int, None],
        train_fit_min_quality: float,
        num_harmonics_sin:     int,
        num_harmonics_cos:     Union[int, None],
        xbar_limit_1:          float,               # todo what this do
        xbar_limit_2:          int,                 # todo what this do
        lambda_value:          float,               # was lambda
        lambda_bounds:         float,               # was lambdaSigs
        rounding:              bool,
        persistence_per_year:  float,
        threshold:             float,               # todo veg threshold
        reverse_order:         bool,
) -> np.ndarray:
    """"""

    #print(f'Working on arr: {arr}')
    #print(f'Working with dates: {dates}')

    # todo do the sort here

    # # check if values array is valid (type, shape)
    # if not isinstance(arr, np.ndarray):
    #     raise TypeError(f"Values array must be type numpy ndarray.")
    # elif arr is None or len(arr) == 0:
    #     raise ValueError(f"Values array is empty.")
    # elif len(arr.shape) != 1:
    #     raise ValueError(f"Values array must be 1-dimensional.")

    # check if dates is adequate
    # if not isinstance(dates, np.ndarray):
    #     raise TypeError(f"Dates array must be type numpy ndarray.")
    # elif dates is None or len(dates) == 0:
    #     raise ValueError(f"Input dates array is empty.")
    # elif len(dates.shape) != 1:
    #     raise ValueError(f"Input dates array must be 1-dimensional.")

    # check if array and dates are same size
    if len(arr) != len(dates):
        raise ValueError(f"Input date and value arrays must be same size.")

    # todo check type for datetime or datetime64
    # get years, doys from dates array
    dates = pd.DatetimeIndex(dates)
    years = dates.year.to_numpy()
    doys = dates.dayofyear.to_numpy()

    # build nan array for analysis period (was subset earlier)
    #arr_empty = np.full(len(dates), np.nan)

    # init break point tracker todo what does it do?
    breaks_tracker = np.arange(0, len(arr))
    breaks_start = np.array([], dtype='int64')  # todo we could not use this with python
    breaks_end = np.array([], dtype='int64')    # todo we could not use this with python

    # init first beta todo what does beta do?
    beta_first = np.full((1 + num_harmonics_sin + num_harmonics_cos), np.nan)

    # apply reverse-toggling, if requested todo why use this? disabled for now
    #if reverse_order is True:
        #arr = np.flip(arr)

    # todo ensure dates, years, doys, arr have no nans (must all be numerics)
    # r script did this via as.numeric...

    # convert dates to decimal dates for ordering
    decimal_years = years + doys / 365

    # no longer sorting vals to time, done in outer func
    # arr = arr[np.argsort(decimal_years)] # etc...

    # no longer doing training end here, do in outer func...
    #

    # no longer doing trim of data to train start and test end, we did in outer func
    # was trims, doys, years, yearsforannualoutput, arr, breaktracker

    # init arrays to store missing values and values under fitting threshold a priori todo what this mean? is this diff to arr_empty?
    arr_nans = np.full(len(arr), np.nan)
    arr_nans_and_residuals = np.full(len(arr), np.nan)  # todo is this ever used?

    # init clean parameters
    clean_idxs = np.flatnonzero(~np.isnan(arr) & (arr > threshold))
    arr_clean = arr[clean_idxs]
    years_clean = years[clean_idxs]
    doys_clean = doys[clean_idxs]
    decimal_years_clean = decimal_years[clean_idxs]
    breaks_tracker_clean = breaks_tracker[clean_idxs]

    # leave if no valid data
    if len(arr_clean) == 0:
        #raise NotImplemented('No implemented yet')  # todo implement this
        return arr_nans, arr_nans, beta_first, np.array([np.nan]), np.array([np.nan])  # todo temp solution
        #output = # todo
        #return

    # check if enough clean data  # todo added this, might want in
    if len(arr_clean) <= min_train_length:
        return arr_nans, arr_nans, beta_first, np.array([np.nan]), np.array([np.nan])
        #raise ValueError('Not enough clean training data.')

    # init vegetation persistence
    persistence = len(arr_clean) / len(np.unique(years_clean))
    persistence = np.ceil(persistence * persistence_per_year)

    # todo we no longer use static or dynamic here, selected in main func
    # todo done the min train length, max train length out of func

    # update clean decimal years  todo do we need this now that we dont trim earlier?
    decimal_years_clean = years_clean + doys_clean / 365

    # hreg = harmonic regression. this method optimises hreg before applying
    # # todo historybound comes out as last element if fit not reached, set to nan instead?
    optimal_results = _optimise_harmonic_regression(
        decimal_years_cleaned=decimal_years_clean,
        doys_cleaned=doys_clean,
        arr_cleaned=arr_clean,
        train_fit_min_quality=train_fit_min_quality,
        min_train_length=min_train_length,
        max_train_length=max_train_length,
        num_harmonics_sin=1,  # todo why 1?
        num_harmonics_cos=1,  # todo why 1?
        xbar_limit_1=xbar_limit_1
    )

    # optimal_result is a tuple of history_idx, fit_quality
    history_idx, train_precedents = optimal_results

    # create break point start and end
    # todo, do we really need to init breaks start above?
    breaks_start_idx = np.append(breaks_start, breaks_tracker_clean[0])
    breaks_end_idx = np.append(breaks_end, breaks_tracker_clean[history_idx])

    # if nothing returned, return empty array
    if np.isnan(history_idx):
        return arr_nans

    # clean pixel
    clean_results = _clean_data(
        years=years_clean,
        doys=doys_clean,
        arr=arr_clean,
        history_idx=history_idx,
        precedents=train_precedents,
        num_harmonics_sin=num_harmonics_sin,
        num_harmonics_cos=num_harmonics_cos,
        xbar_limit_1=1.5,     # todo accept this from parent func inputs
        xbar_limit_2=20,      # todo accept this from parent func inputs
        lambda_value=0.3,     # was lambda # todo accept this from parent func inputs
        lambda_bounds=3.0,    # was lambdaSigs # todo accept this from parent func inputs
        rounding=True,        # todo accept this from parent func inputs
        persistence=persistence,
    )

    # unpack elements from result (clean data (keeps), clean residuals, beta)
    run_keeps = clean_results[0]
    run_keeps_residuals = clean_results[1]
    first_beta = clean_results[2]  # todo already have var called beta_first, do we wanna override for readabilty?

    # post-processing of non-missing signals filtered by persistence
    final_chart_arr = _backfill(
        non_missing=run_keeps,
        non_missing_idx=clean_idxs,
        with_missing=arr_nans.copy()  # todo added copy. may want to do a copy in back fill instead
    )

    final_residuals_arr = _backfill(
        non_missing=run_keeps_residuals,
        non_missing_idx=clean_idxs,
        with_missing=arr_nans_and_residuals.copy()  # todo added copy. may want to do a copy in back fill instead
    )

    # todo simple output?
    # todo

    return final_chart_arr, final_residuals_arr, first_beta, breaks_start_idx, breaks_end_idx


def ewmacd(
        ds: xr.Dataset,  # update to support data array
        var: str,
        train_start_year: Union[int, None] = None,
        #train_end_year: Union[int, None] = None,    # todo disabled for now
        test_end_year: Union[int, None] = None,
        min_train_length: Union[int, None] = None,
        max_train_length: Union[int, None] = None,
        train_fit_min_quality: float = 0.8,
        num_harmonics_sin: int = 2,
        num_harmonics_cos: Union[int, None] = 2,
        xbar_limit_1: float = 1.5,
        xbar_limit_2: int = 20,
        lambda_value: float = 0.3,
        lambda_bounds: float = 3.0,
        rounding: bool = True,
        persistence_per_year: float = 1.0,
        threshold: float = 0.0,
        reverse_order: bool = False,
        summarise: str = 'raw'
) -> xr.Dataset:
    """
    Perform EWMACD change detection on a Xarry Dataset.

    :param xr.Dataset values: Dataset of 2D time-series data with x, y, time dimensions.
                          Dataset must contain x, y and time dimensions.
    :param str var: Dataset variable name to apply change detection. Variable
                    name must be provided and must exist in dataset.
    :param str method: Change detection method. Supports 'static' and 'dynamic'
                       methods. Default method is set to 'static'.
    :param int train_start_year: Year of the start of training period. This year
                                 is included in training period.
    :param int test_end_year: Year of testing period's end. The testing period
                              is the period that will be compared to the training
                              period. Setting this to None will use all available
                              years in dataset.
    :param int min_train_length: Minimum number of dates to use when training the
                                 harmonic regression model. Use this to fine tune
                                 the minimum number of dates the harmonic regression
                                 model will use when optimising. Set to None to
                                 let the optimiser solve this itself.
    :param int max_train_length: Maximum number of dates to use when training the
                                 harmonic regression model. Use this to fine tune
                                 the maximum number of dates the harmonic regression
                                 model will use when optimising. Set to None to
                                 let the optimiser solve this itself.
    :param float train_fit_min_quality: Minimum harmonic regression training fit.
                                        The harmonic regression optimiser will
                                        iterate through dates within the training
                                        period until it reaches this quality fit.
                                        Once reached, the harmonic model is used.
    :param int num_harmonics_sin: Number of harmonic sine functions to use
                                  during harmonic regression.
    :param int num_harmonics_cos: Number of harmonic cosine functions to use
                                  during harmonic regression. Set this to
                                  None to use the same number as sine harmonics.
    :param float xbar_limit_1: Threshold for initial training xBar limit.
    :param int xbar_limit_2: Threshold for running xBar limit.
    :param float lambda_value: The lambda tuning parameter weighting new years
                               versus the running average.
    :param float lambda_bounds: EWMA control bounds, in units of standard
                                deviation.
    :param rounding: Perform rounding for EWMA.
    :param persistence_per_year: Fraction of observations per year required in
                                 order for a change to be registered on the
                                 chart.
    :param threshold: Threshold for vegetation. Values below this are considered
                      non-vegetation. Threshold value must be in units of index
                      captured in the specified var parameter.

    :return: xr.Dataset
    """

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # check and prepare numerous input parameters
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # check if xr dataset is valid (type, shape)
    if not isinstance(ds, xr.Dataset):
        raise TypeError('Input dataset must be an Xarray Dataset type.')
    elif 'time' not in ds:
        raise ValueError('Input dataset does not have time coordinates/dimension.')
    elif ('x' not in ds or 'y' not in ds) and ('lon' not in ds or 'lat' not in ds):
        raise ValueError('Input dataset does not have x, y, lon, or lat coordinates/dimensions.')

    # check if requested var is valid
    if var not in ds:
        raise ValueError(f"Variable '{var}' not in input dataset.")

    # prepare training start year (first year if none or not exist)
    if train_start_year not in ds['time.year']:
        train_start_year = np.min(ds['time.year'].values)
        print(f'Inadequate training start year provided, using {train_start_year}.')

    # prepare training end year (first year + 3 if none or not exist)  # todo disabled for now
    #if train_end_year not in ds['time.year']:
        #train_end_year = train_start_year + 3
        #print(f'Inadequate training end year provided, using {train_end_year}.')

    # prepare testing end year (last year + 1 if none or not exist)
    if test_end_year not in ds['time.year']:
        test_end_year = np.max(ds['time.year'].values) + 1
        print(f'Inadequate testing end year provided, using {test_end_year}.')

    # prepare min, max training length
    if min_train_length is None:
        min_train_length = (1 + num_harmonics_sin + num_harmonics_cos) * 3
        print(f'No minimum training length provided, using {min_train_length}.')

    # prepare max training length
    if max_train_length is None:
        max_train_length = min_train_length * 2
        print(f'No maximum training length provided, using {max_train_length}.')

    # check min, max train length valid
    if min_train_length >= max_train_length:
        raise ValueError('Training length minimum cannot be <= maximum length.')

    # check if training fit quality is valid
    if train_fit_min_quality <= 0.0 or train_fit_min_quality > 1.0:
        raise ValueError('Training fit quality must be between 0.0 and 1.0.')

    # todo check if sin, cos can mismatch
    # prepare number of sine and cosine harmonics
    if num_harmonics_sin <= 0 or num_harmonics_cos <= 0:
        raise ValueError('Number of sine and cosine harmonics must be > 0.')

    # todo check if cos none, set to sin num

    # prepare initial training xbar limit 1 and 2
    if xbar_limit_1 <= 0:
        raise ValueError('Xbar limit 1 must > 0.0.')
    elif isinstance(xbar_limit_2, float) or xbar_limit_2 <= 0.0:
        raise ValueError('Xbar limit 2 must be an integer >= 1.')

    # prepare lambda value and bounds
    if lambda_value <= 0.0 or lambda_bounds <= 0.0:
        raise ValueError('Lambda value and bounds must be > 0.0.')

    # prepare ewma rounding
    if rounding not in [True, False]:
        raise ValueError('Rounding must be either True or False.')

    # check persistence is valid
    if persistence_per_year <= 0.0:
        raise ValueError('Persistence must be > 0.0.')

    # check threshold is valid
    if threshold is None:
        raise ValueError('Threshold not provided.')

    # prepare reverse order
    if reverse_order not in [True, False]:
        raise ValueError('Reverse order must be either True or False.')

    # prepare summarise is valid
    summaries = ['raw', 'mean', 'median', 'extreme', 'sign_mean', 'sign_median']
    if summarise not in summaries:
        raise ValueError('Summary method not supported.')

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # prepare xr dataset
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # create a copy of the xr dataset and ensure sort by time
    ds_tmp = ds[[var]].copy(deep=True).sortby('time')

    # subset xr dataset to train start and test end
    ds_tmp = ds_tmp.sel(time=slice(str(train_start_year),
                                   str(test_end_year)))

    # extract 1-dimensional numpy array of dates
    dates = ds_tmp['time'].values

    # stack dataset depending on x, y or lon, lat
    dims = ['y', 'x'] if 'x' in ds and 'y' in ds else ['lon', 'lat']
    ds_tmp = ds_tmp.stack(z=dims)
    das = ds_tmp[var].transpose().values

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # perform ewmacd per-pixel via multiprocessing
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # prepare parameter kwargs for convenience
    kwargs = {
        'dates':                 dates,
        'train_start_year':      train_start_year,
        #'train_end_year':       train_end_year,  # todo disabled
        'test_end_year':         test_end_year,
        'min_train_length':      min_train_length,
        'max_train_length':      max_train_length,
        'train_fit_min_quality': train_fit_min_quality,
        'num_harmonics_sin':     num_harmonics_sin,
        'num_harmonics_cos':     num_harmonics_cos,
        'xbar_limit_1':          xbar_limit_1,
        'xbar_limit_2':          xbar_limit_2,
        'lambda_value':          lambda_value,
        'lambda_bounds':         lambda_bounds,
        'rounding':              rounding,
        'persistence_per_year':  persistence_per_year,
        'threshold':             threshold,
        'reverse_order':         reverse_order
    }

    # start timer
    start = time.time()

    # map and process each ewmacd onto each pixel
    # results is tuple of chart, residuals, coeffs (beta), break start, break end values
    cpus = multiprocessing.cpu_count() - 1
    with multiprocessing.Pool(processes=cpus) as pool:
        results = pool.map(partial(_ewmacd_core_func, **kwargs), das)

    # notify user of execution time
    print(time.time() - start)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # convert results and add to xarray dataset
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # add chart values from results to new stacked variable
    ds_tmp['charts'] = xr.full_like(ds_tmp[var], np.nan)
    ds_tmp['charts'].data = np.array([e[0] for e in results]).T

    # add residual values from results to new stacked variable
    ds_tmp['residuals'] = xr.full_like(ds_tmp[var], np.nan)
    ds_tmp['residuals'].data = np.array([e[1] for e in results]).T

    # add coefficients (beta) values from results to new stacked variable  todo implement this
    #ds_tmp['coefficients'] = xr.full_like(ds_tmp[var], np.nan)
    #ds_tmp['coefficients'].data = np.array([e[2] for e in results]).T

    # unstack back to original shape
    ds_tmp = ds_tmp.transpose().unstack()

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # summarise dataset
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if summarise == 'raw':
        return ds_tmp
    elif summarise == 'mean':
        ds_tmp = ds_tmp.resample(time='1YS').mean()
        ds_tmp['charts'] = ds_tmp['charts'].round()  # todo round all or just charts?
        return ds_tmp
    elif summarise == 'median':
        ds_tmp = ds_tmp.resample(time='1YS').median()
        ds_tmp['charts'] = ds_tmp['charts'].round()  # todo round all or just charts?
        return ds_tmp
    elif summarise == 'extreme':
        ds_min = ds_tmp.resample(time='1YS').min()
        ds_max = ds_tmp.resample(time='1YS').max()
        ds_tmp = xr.where(abs(ds_min) > ds_max, ds_min, ds_max)
        ds_tmp['charts'] = ds_tmp['charts'].round()
        return ds_tmp
    elif summarise == 'sign_mean':
        raise NotImplemented('Not yet implemented.')
    elif summarise == 'sign_median':
        raise NotImplemented('Not yet implemented.')

    # todo do monthly output
    # monthly mean per year
    #yds = ds_tmp.time.dt.strftime('%Y-%m')
    #yds = np.array([np.datetime64(e + '-' + '01') for e in yds.values])
    #ds_tmp = ds_tmp.assign_coords({'time': yds}).groupby('time').mean()

    # monthly median per year
    #yds = ds_tmp.time.dt.strftime('%Y-%m')
    #yds = np.array([np.datetime64(e + '-' + '01') for e in yds.values])
    #ds_tmp = ds_tmp.assign_coords({'time': yds}).groupby('time').median()



    return None



# todo cube
def convert_to_cube():
    ...
    # import datetime
    #
    # arr = arr.rename({'band': 'time'})
    #
    # dts = []
    # for y, d in zip(years, doys):
    #     dt = datetime.datetime(int(y), 1, 1) + datetime.timedelta(int(d) - 1)
    #     dt = dt.strftime('%Y-%m-%d')
    #     dts.append(np.datetime64(dt))
    #
    # ds = arr.to_dataset(name='veg_idx')
    # ds['time'] = dts
    #
    # ds.to_netcdf(r'C:\Users\Lewis\PycharmProjects\py-ewmacd\tests\data\angel_island.nc')

def tifs_to_cube(folder):
    """

    :param folder:
    :return:
    """

    # todo: clean this up
    # load tifs
    das = []
    files = os.listdir(folder)
    for file in files:
        if file.endswith('.tif'):
            dt = file.split('.')[0]
            dt = pd.to_datetime(dt, format='%Y-%m-%d').to_datetime64()

            da = xr.open_dataset(os.path.join(folder, file))
            da = da.to_array().squeeze(drop=True)

            da = da.assign_coords({'time': dt})
            da = da.expand_dims('time')

            das.append(da)

    ds = xr.concat(das, 'time').to_dataset(name='veg_idx')
    ds = ds.sortby('time')

    return ds


if __name__ == '__main__':

    # todo: 99% comparbable to r
    # todo: some issues with missing pixels (set to nan here instead of 0 in r)
    # todo: also, a few small missing pixels early on - model sensitivty differences?


    # load dataset with x, y and time dims
    ds = xr.open_dataset(r'./tests/data/angel_island.nc')

    # load folder of tifs as cube
    #ds = tifs_to_cube(folder=r'./tests/data/bindi')

    #ds = ds.isel(x=slice(0, 250), y=slice(0, 250))

    ds = ewmacd(
        ds=ds,
        var='veg_idx',
        train_fit_min_quality=0.95,
        #persistence_per_year=1.0,
        #num_harmonics_sin=3,
        summarise='raw'
    )

    ds.to_netcdf('py_out.nc')

    #fig = plt.figure(figsize=[10, 8])
    #ds['charts'].isel(time=30).plot(cmap='RdYlBu', robust=False)
    #plt.show()

    # ds_r = xr.open_dataset('r_out.nc')
    # ds_p = xr.open_dataset('py_out.nc')
    #
    # ds_r = ds_r.rename({'band': 'time'})
    # ds_r = ds_r.rename({'band_data': 'charts'})
    #
    # i = 13
    # fig = plt.figure(figsize=[10, 8])
    # ds_r['charts'].isel(time=i).plot(cmap='RdYlBu', robust=True)
    # plt.show()
    #
    # fig = plt.figure(figsize=[10, 8])
    # ds_p['charts'].isel(time=i).plot(cmap='RdYlBu', robust=True)
    # plt.show()
    #
    # fig = plt.figure(figsize=[10, 8])
    # ds_d = ds_r['charts'].isel(time=i) - ds_p['charts'].isel(time=i)
    # ds_d.plot(cmap='RdYlBu', robust=True)
    # plt.show()