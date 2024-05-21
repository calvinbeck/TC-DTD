#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis
"""


import plotfunction
import functions


def run(df_sel, step, steps, data, r_squared, data_period,
        r_squared_period, depth, depth_vec, levels, dates_list, Delta_t,
        site, sel_average, save_figure):
    """
    Perform main analysis for selected sampling interval or averaging time.
    Calculate gradients and then perform linear regression for each day.
    """

    if sel_average:
        resample_string = str(step*Delta_t) + 'min'
        df = df_sel.resample(resample_string).mean()
    else:
        df = df_sel.iloc[::step]

    # Calculate second spacial derivative:
    df_d2Tdx2 = functions.calc_d2Tdx2(df, depth, depth_vec)
    # Calculate first temporal derivative:
    df_dTdt = functions.calc_dTdt(df, depth)

    # perform linear regression for complete period at once 
    for level in range(1, levels+1):
        param, intercept, r2 = functions.linres(df_dTdt, df_d2Tdx2,
                                                depth[level])
        data_period[level-1, step-1] = param
        r_squared_period[level-1, step-1] = r2

    if step == 1:
        plotfunction.gradients(df_d2Tdx2, df_dTdt, depth_vec, site, step,
                               save_figure)

    days = 0
    # perform linear regression (day by day)
    for idx in dates_list:
        df_dTdt_day = df_dTdt[df_dTdt.index.date == idx]
        df_d2Tdx2_day = df_d2Tdx2[df_d2Tdx2.index.date == idx]
        for level in range(1, levels+1):
            param, intercept, r2 = functions.linres(df_dTdt_day, df_d2Tdx2_day,
                                                    depth[level])
            data[level-1, days, step-1] = param
            r_squared[level-1, days, step-1] = r2
        # display mean layer temperature
        if step == 1:
            plotfunction.mean_layer_temperature(df, depth_vec, idx, days,
                                                dates_list, site, save_figure)
        days = days + 1

    return data, r_squared, data_period, r_squared_period
