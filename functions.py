#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions
"""

from datetime import timedelta, datetime
import numpy as np
import statsmodels.api as sm
import pandas as pd
from scipy.sparse import diags
from tqdm import tqdm
import sys
import os

import equations


def read_csv_file(path):
    df = pd.read_csv(path, delimiter=",")
    df['datetime'] = pd.to_datetime(df['datetime'])  # convert sring to dt
    index_list = df['datetime']
    df = df.set_index('datetime')
    df = df.loc[:, df.columns != 'T_air']  # neglect column with air themp.
    return df, index_list


def select_time_period(df, index_list, sel_period=False, from_date=None,
                       to_date=None):
    df['datetime'] = df.index
    if sel_period:
        from_date = datetime.strptime(from_date, '%Y-%m-%d')
        to_date = datetime.strptime(to_date, '%Y-%m-%d')
        df = df[df['datetime'].between(from_date, to_date)]
        if df.empty:
            sys.exit("Sorry, but there is no data in the selected period!")
    index_list = df.index
    day_count = df.datetime.resample('1D').count()
    if day_count[0] != np.max(day_count):
        df = df.drop(df[df.index.date == day_count.index[0].date()].index)
        print('First day neglected due to incomplete daily data.')
    if day_count[-1] != np.max(day_count):
        df = df.drop(df[df.index.date == day_count.index[-1].date()].index)
        print('Last day neglected due to incomplete daily data.')
    df = df.loc[:, df.columns != 'datetime']
    return df, index_list


def linres(df_dTdt, df_d2Tdx2, select_depth):
    Y = df_dTdt[select_depth]
    X = df_d2Tdx2[select_depth]
    X = sm.add_constant(X)
    result = sm.OLS(Y, X, missing='drop').fit()

    if len(result.params) == 1:
        slope = result.params[0]
        intercept = np.nan
    else:
        intercept = result.params[0]
        slope = result.params[1]
    return slope, intercept, result.rsquared


def sampling_interval(index_list, max_minutes):
    Delta_t = (index_list[1]-index_list[0]).total_seconds()/60
    n_max = int((max_minutes/Delta_t)+1)
    steps = np.arange(1, n_max, 1, dtype='int')
    steps_list = np.arange(Delta_t, max_minutes+1, Delta_t, dtype='int')
    return Delta_t, steps, steps_list


def depth_from_string(depth):  # EDIT check if columns are named correctly
    depth_vec = np.zeros(len(depth))
    for n in range(len(depth)):  # loop over different depths
        depth_vec[n] = float(depth[n][2:6])  # select depth from string
        if len(depth[n]) >= 8:
            depth_vec[n] = float(depth[n][2:7])
    return depth_vec

def list_duplicates(seq):  # EDIT
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in seq if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)

def calc_d2Tdx2(df, depth, depth_vec):
    df_d2Tdx2 = df.copy()
    for height in range(len(depth)):
        if height == 0 or height+1 == len(depth):
            df_d2Tdx2[depth[height]] = np.NaN  # boundary condition
        else:
            Delta_x = np.array([depth_vec[height-1], depth_vec[height],
                                depth_vec[height+1]])
            df_d2Tdx2[depth[height]] = equations.d2Tdx2(df[depth[height-1]],
                                                        df[depth[height]],
                                                        df[depth[height+1]],
                                                        Delta_x)
    return df_d2Tdx2


def calc_dTdt(df, depth):
    df = df.transpose()
    df_dTdt = df.copy()
    timestamp = list(df.columns)
    Delta_t = (timestamp[1]-timestamp[0]).total_seconds()
    for time in range(len(timestamp)):
        if time == 0 or time+1 == len(timestamp):  # boundary condition
            df_dTdt[timestamp[time]] = np.NaN
        else:
            df_dTdt[timestamp[time]] = equations.dTdt(df[timestamp[time-1]],
                                                      df[timestamp[time+1]],
                                                      Delta_t)
    df = df.transpose()
    df_dTdt = df_dTdt.transpose()
    return df_dTdt

def dates(df):
    sdate = df.index.date[0]   # start date
    edate = df.index.date[-1]   # end date
    delta = edate - sdate       # as timedelta
    dates_list = []

    for i in range(delta.days + 1):
        dates_list.append(sdate + timedelta(days=i))

    return dates_list


def export_csv_file(Tout, dt, mtype, depths, layer_thickness, days,
                    neglect_days, sel_model):
    """export as csv file"""
    depths[0] = 0
    df = pd.DataFrame(data=np.transpose(Tout))
    new_names = []

    for name in range(len(depths)):
        new_names.append('T_' + str("{:.3f}".format(depths[name])) + 'm')
    df.columns = new_names
    df = df.drop('T_0.000m', 1)  # EDIT

    # create datetime timeseries for dataframe:
    start_date = datetime(2021, 1, 1)
    end_date = start_date + timedelta(days=days) + timedelta(seconds=dt)
    delta_time = timedelta(seconds=dt)
    time_arr = np.arange(start_date, end_date, delta_time).astype(datetime)
    df.insert(0, 'datetime', time_arr)
    exclude_period = start_date + timedelta(days=neglect_days)
    df = df[~df['datetime'].between(start_date, exclude_period)]

    model_depth = '_depth_' + str(round(layer_thickness * 100)) + 'cm'
    if not os.path.exists('saved_model_data/'):
        os.makedirs('saved_model_data/')
    df.to_csv('saved_model_data/ModelData_' + sel_model + model_depth + '.csv',
              index=False)


def crank_nicolson(T, s, nx, nt, time, dt, top_mean_T, osci_amplitude,
                   sel_model):
    """Use crank nicolson scheme to calculate heat diffusion"""
    Tout = np.zeros((nx, len(time)))
    T0 = T[0]
    T1 = T[-1]
    A = diags([-0.5*s, 1+s, -0.5*s], [-1, 0, 1], shape=(nx-2, nx-2)).toarray()
    B1 = diags([0.5*s, 1-s, 0.5*s], [-1, 0, 1], shape=(nx-2, nx-2)).toarray()
    for n in tqdm(range(1, len(time)), position=0, leave=True):
        T0, T = temperature_function(T, n, time, nx, dt, top_mean_T,
                                     osci_amplitude, sel_model)
        Tn = T
        Tn[0] = T0
        B = np.dot(Tn[1:-1], B1)
        B[0] = B[0]+0.5*s*(T0 + T0)
        B[-1] = B[-1]+0.5*s*(T1 + T1)
        T[1:-1] = np.linalg.solve(A, B)
        Tout[:, n] = T
    print(s)
    return Tout


def temperature_function(T, n, time, nx, dt, top_mean_T, osci_amplitude,
                         sel_model):
    """select different intia functions"""
    if sel_model == "skewed_sine":
        T0 = top_mean_T - osci_amplitude*np.cos(2*np.pi*n*dt/(24*60*60) -0.5*np.cos(2*np.pi*n*dt/(24*60*60))) #squed sine curve
    elif sel_model == "pure_sine":
        T0 = top_mean_T - osci_amplitude*np.cos(2*np.pi*n*dt/(24*60*60))
    elif sel_model == "example_data_1":
        real_data = np.array(pd.read_csv('example_temperature_data/test_data_1.csv'))
    elif sel_model == "example_data_2":
        real_data = np.array(pd.read_csv('example_temperature_data/test_data_2.csv'))
    elif sel_model == "example_data_3":
        real_data = np.array(pd.read_csv('example_temperature_data/test_data_3.csv'))
    else:
        T0 = 15 + (0.25*np.cos(3*np.pi*n/len(time)+np.pi/3)+1.25)*10*np.sin(2*np.pi*n/(24*60)) + 4*n/len(time)*np.sin(7*np.pi*n/len(time)+np.pi/7)
    if sel_model[:7] == "example":
        T0 = real_data[n]
        if n == 1:
            T = np.linspace(np.mean(real_data), 0, num=nx)
    return T0, T
