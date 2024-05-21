#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotfunctions
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os


def temperature_timeseries(dates_list, df, depth, site, save_figure,
                           from_date=None, to_date=None):
    plt.figure(figsize=(12, 6))
    cmap = cm.inferno_r
    plt.plot(df.index, df[depth[0]]*0, color="Black", linestyle='--',
             label='_nolegend_')
    for d in range(len(depth)):
        if len(depth) <= 8:
            plt.plot(df.index, df[depth[d]], color=cmap((d+1)/(len(depth)+1)),
                     label=depth[d][2:])
        else:
            label_depths = np.round(np.linspace(0, len(depth), num=6), 0)
            if d in label_depths:
                plt.plot(df.index, df[depth[d]],
                         color=cmap((d+1)/(len(depth)+1)),
                         label=depth[d][2:])
                plt.plot(df.index, df[depth[d]],
                         color=cmap((d+1)/(len(depth)+1)),
                         alpha=0, label='⋯')
            else:
                plt.plot(df.index, df[depth[d]],
                         color=cmap((d+1)/(len(depth)+1)),
                         label='_nolegend_')
    plt.xlim(dates_list[0], dates_list[-1])
    plt.ylabel("Temperature (°C)")
    if len(depth) <= 8:
        plt.legend(loc='lower center', handletextpad=1, handlelength=2,
                   shadow=True, ncol=10)
    else:
        plt.legend(loc='lower center', handletextpad=0, handlelength=2,
                   shadow=True, ncol=10)
    plt.grid(linestyle='--')
    plt.title('Temperature cycle in debris layer (' + site + ')')
    if not os.path.exists('saved_figures/' + site + '_figures'):
        os.makedirs('saved_figures/' + site + '_figures')
    if from_date is None or to_date is None:
        plt.savefig('saved_figures/' + site + '_figures/' + site
                    + '_temperatur_timeseries.pdf')
    else:
        plt.savefig('saved_figures/' + site + '_figures/' + site
                    + '_temperatur_timeseries_' + from_date + '-' + to_date
                    + '.pdf')
    plt.show()


def kappa_values_timeseries(dates, data, r_squared, depth, site, sel_average,
                            save_figure):
    cmap = cm.inferno_r
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=False)
    for n in range(2):
        if n == 1:
            data = r_squared
        for level in range(len(data[:, 0])):
            if len(depth) <= 8:
                axs[n].plot(dates, data[level], label=depth[level+1][2:],
                            marker='x', color=cmap((level+2)/(len(depth)+1)))
            else:
                label_depths = np.round(np.linspace(0, len(depth), num=6), 0)
                if level in label_depths:
                    axs[n].plot(dates, data[level], label=depth[level+1][2:],
                                marker='x',
                                color=cmap((level+2)/(len(depth)+1)))
                    axs[n].plot(dates, data[level], label='⋯', marker='x',
                                color=cmap((level+2)/(len(depth)+1)),
                                alpha=0)
                else:
                    axs[n].plot(dates, data[level], marker='x',
                                color=cmap((level+2)/(len(depth)+1)),
                                label='_nolegend_')
        if n == 0:
            axs[n].set_title("Timeseries for measurement site: " + site)
            axs[n].set_ylabel(r'Thermal diffusivity $\kappa$ (m²/s)')
            data = data[~ np.isnan(data)]
            max_val = np.max(data)
            axs[n].set_ylim(0, 1.15*max_val)
            axs[n].text(dates[int(len(dates)/2)], 1.10*max_val,
                        r'Thermal diffusivity by day',
                        verticalalignment='top',
                        horizontalalignment='center',
                        fontsize=14, fontweight='bold',)
        elif n == 1:
            axs[n].set_ylim(0, 1.15)
            axs[n].set_ylabel('Coefficient of determination R²')
            axs[n].text(dates[int(len(dates)/2)], 1.1,
                        ' R² to corresponding linear regression',
                        verticalalignment='top',
                        horizontalalignment='center',
                        fontsize=14, fontweight='bold')
    axs[1].legend(loc='lower center', shadow=True, ncol=10)
    for ax in axs:
        ax.set_xlim(dates[0], dates[-1])
        ax.label_outer()
        ax.grid(linestyle='--')
    if not os.path.exists('saved_figures/' + site + '_figures'):
        os.makedirs('saved_figures/' + site + '_figures')
    if sel_average:
        plt.savefig('saved_figures/' + site + '_figures/' + site
                    + '_kappa_timeseries_avg.pdf')
    else:
        plt.savefig('saved_figures/' + site + '_figures/' + site
                    + '_kappa_timeseries.pdf')
    plt.show()


def resampling_time(data, data_period, r_squared, r_squared_period, steps_list,
                    site, depth, sel_average, save_figure):
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=False)
    for n in range(2):
        if n == 1:
            data = r_squared
            data_period = r_squared_period
        avg = np.average(data[depth, :, :], axis=0)
        std = np.std(data[depth, :, :], axis=0)
        for day in range(len(data[depth, :, 0])):
            if day == 0 and n == 1:
                axs[n].plot(steps_list, data[depth, day, :], color='black',
                            alpha=0.2, label='Single day')
            else:
                axs[n].plot(steps_list, data[depth, day, :], color='black',
                            alpha=0.2)
        if n == 1:
            axs[n].plot(steps_list, avg, color='#003361', label='Mean')
            axs[n].plot(steps_list, avg + std, color='#f39200',
                        label='Mean ± σ')
            axs[n].plot(steps_list, avg - std, color='#f39200')
            axs[n].plot(steps_list, data_period[depth, :], color='red',
                        label='Complete period')
        else:
            axs[n].plot(steps_list, avg, color='#003361')
            axs[n].plot(steps_list, avg + std, color='#f39200',)
            axs[n].plot(steps_list, avg - std, color='#f39200')
            axs[n].plot(steps_list, data_period[depth, :], color='red')
        if n == 0:
            axs[n].set_title("Timeseries for measurement site: " + site)
            axs[n].set_ylabel(r'Thermal diffusivity $\kappa$ (m²/s)')
            if len(steps_list) >= 200:
                num = np.where(steps_list >= 200)
                num = num[0][0]
            else:
                num = len(steps_list)
            max_val = np.percentile(data[depth, :, :num], 90)
            axs[n].set_ylim(0, 1.3*max_val)
            if sel_average:
                axs[n].text(steps_list[int(len(steps_list)/2)], 1.25*max_val,
                            r'Thermal diffusivity by averaging interval',
                            verticalalignment='top',
                            horizontalalignment='center',
                            fontsize=14, fontweight='bold',)
            else:
                axs[n].text(steps_list[int(len(steps_list)/2)], 1.25*max_val,
                            r'Thermal diffusivity by sampling interval',
                            verticalalignment='top',
                            horizontalalignment='center',
                            fontsize=14, fontweight='bold',)
        elif n == 1:
            axs[n].set_ylim(0, 1.15)
            axs[n].set_ylabel('Coefficient of determination R²')
            if sel_average:
                axs[n].text(steps_list[int(len(steps_list)/2)], 1.1,
                            'R² by averaging interval',
                            verticalalignment='top',
                            horizontalalignment='center',
                            fontsize=14,
                            fontweight='bold')
            else:
                axs[n].text(steps_list[int(len(steps_list)/2)], 1.1,
                            'R² by sampling interval', verticalalignment='top',
                            horizontalalignment='center', fontsize=14,
                            fontweight='bold')
    axs[1].legend(loc='lower center', shadow=True, ncol=10)
    axs[1].set_xlabel('Sampling time (minutes)')
    for ax in axs:
        ax.set_xlim(0, steps_list[-1])
        ax.label_outer()
        ax.grid(linestyle='--')
        if sel_average:
            ax.set_facecolor("#fff6ea")
        else:
            ax.set_facecolor("#eafbff")
    if not os.path.exists('saved_figures/' + site + '_figures'):
        os.makedirs('saved_figures/' + site + '_figures')
    if sel_average:
        plt.savefig('saved_figures/' + site + '_figures/' + site
                    + '_kappa_avg_dependence.pdf')
    else:
        plt.savefig('saved_figures/' + site + '_figures/' + site
                    + '_kappa_sampling_dependence.pdf')
    plt.show()


def mean_layer_temperature(df, depth_vec, idx, days, dates_list, site,
                           save_figure):
    plt.grid(linestyle='--')
    cmap = cm.viridis_r
    T_m_day = np.array(np.mean(df[df.index.date == idx]))
    plt.grid(linestyle='--')
    if len(dates_list) <= 8:
        plt.plot(T_m_day, depth_vec, color=cmap((days+1)/(len(dates_list)+1)),
                 label=idx)
    else:
        label_days = np.round(np.linspace(0, len(dates_list), num=6), 0)
        if days in label_days:
            plt.plot(T_m_day, depth_vec,
                     color=cmap((days+1)/(len(dates_list)+1)), label=idx)
            plt.plot(T_m_day, depth_vec, alpha=0,
                     color=cmap((days+1)/(len(dates_list)+1)), label='⁞')
        else:
            plt.plot(T_m_day, depth_vec,
                     color=cmap((days+1)/(len(dates_list)+1)))
    if idx == dates_list[-1]:
        plt.gca().invert_yaxis()
        plt.legend(loc='lower right', shadow=True, ncol=1)
        plt.title('Mean debris layer temperature per day and depth:')
        plt.ylabel('Depth (m)')
        plt.xlabel('Temperature (°C)')
        if not os.path.exists('saved_figures/' + site + '_figures'):
            os.makedirs('saved_figures/' + site + '_figures')
        plt.savefig('saved_figures/' + site + '_figures/' + site
                    + '_mean_layer_temperature.pdf')
        plt.show()


def mean_k_by_layer(data, depth_vec, site, save_fig):
    plt.grid(linestyle='--')
    # plt.plot(depth_vec[1:-1], k_avg[:], marker='x', label='mean value')
    k_avg = np.average(data[:, :, 0], axis=1)

    plt.plot(depth_vec[1:-1], k_avg[:], marker='x', label='mean value')
    # plt.plot(depth_vec[2:-1], k_avg[1:]*0+5e-7)
    plt.title('Thermal diffusivity by depth:')
    plt.xlabel('Depth (m)')
    plt.ylabel(r'Thermal diffusivity $\kappa$ (m²/s)')
    plt.legend(loc='upper right', shadow=True, ncol=1)
    plt.ylim(0, 1.1*np.max(k_avg))
    if not os.path.exists('saved_figures/' + site + '_figures'):
        os.makedirs('saved_figures/' + site + '_figures')
    plt.savefig('saved_figures/' + site + '_figures/' + site
                + '_mean_layer_temperature.pdf')
    plt.show()


def gradients(df_d2Tdx2, df_dTdt, depth_vec, site, step, save_figure):
    fig, ax1 = plt.subplots(figsize=(6, 4))
    color = 'tab:red'
    ax1.set_ylabel(r'$\frac{d^2T}{dx^2}$', color=color)
    df_d2Tdx2_plot = np.abs(df_d2Tdx2).mean(axis=0)
    ax1.plot(depth_vec[1:-1], df_d2Tdx2_plot[1:-1], marker='+', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    # Adding Twin Axes to plot using dataset_2
    ax2 = ax1.twinx()
    color = 'tab:green'
    ax2.set_ylabel(r'$\frac{dT}{dt}$', color=color)
    df_dTdt_plot = np.abs(df_dTdt).mean(axis=0)
    ax2.plot(depth_vec[1:-1], df_dTdt_plot[1:-1], marker='x', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax1.set_xlabel('Depth (m)')
    ax1.set_title(r'Gradients: $\frac{d^2T}{dx^2}$ vs. $\frac{dT}{dt}$')
    if save_figure:
        if not os.path.exists('saved_figures/' + site + '_figures'):
            os.makedirs('saved_figures/' + site + '_figures')
        plt.savefig('saved_figures/' + site + '_figures/' + site + '_'
                    + str(step) + '_gradients.pdf')
    plt.show()


def mean_layer_temp(Tout, nx, x, dt, days, depth, neglect_days, model_type,
                    save_figure):
    """
    Plot mean layer temperature
    """
    sec = 24*60*60  # seconds per day
    d_m = np.zeros(nx)
    plt.figure(figsize=(6, 6))
    for n in range(days):
        for m in range(nx):
            d_m[m] = np.mean(Tout[m, int((n)*sec/dt):int((n+1)*sec/dt)])
        if n > neglect_days:
            plt.plot(x, d_m, linewidth=1, label='Day ' + str(n))
    plt.title('Mean layer temperature (' + model_type + ')')
    plt.xlabel('Debris thickness (m)')
    plt.ylabel('Temperature (°C)')
    plt.legend()
    plt.grid()
    if save_figure:
        if not os.path.exists('saved_figures/' + model_type + '_figures'):
            os.makedirs('saved_figures/' + model_type + '_figures')
        plt.savefig('saved_figures/' + model_type + '_figures/' + model_type
                    + '_mean_layer_temp_' + str(depth) + 'm.pdf')
    plt.show()


def model_timeseries(Tout, nx, x_max, time, convert, depth, neglect_days,
                     model_type, save_figure):
    """
    Plot
    """
    cmap = cm.inferno_r
    for n in range(6):
        q = int(n * nx/6)
        label = str(round(n * x_max / 6, 2)) + 'm'
        plt.plot(time / (convert*24), Tout[q, :], color=cmap(n/8+2/8),
                 label=label)
    plt.plot(time / (convert*24), Tout[-1, :], color=cmap(1.),
             label=str(x_max) + 'm')
    plt.xlabel('Time (days)')
    plt.ylabel('Temperature (°C)')
    plt.tight_layout()
    plt.title('Numerical model: Daily temperature cycle in debris layer ('
              + model_type + ')')
    plt.xlim(0, time[-1] / (convert * 24))
    plt.legend(loc='lower right', shadow=True, ncol=1)
    n = np.linspace(0, neglect_days, 10)
    plt.fill_between(n, n*0+np.max(Tout)+0.5, color='#cccccc')
    plt.fill_between(n, n*0+np.min(Tout)-0.5, color='#cccccc')
    plt.ylim(np.min(Tout)-0.5, np.max(Tout)+0.5)
    plt.text(neglect_days/2, np.max(Tout)-0.5, 'Init.', fontsize=12,
             weight='bold', horizontalalignment='center',
             verticalalignment='center', color='black')
    plt.text(neglect_days/2, np.max(Tout)-2.5, 'period', fontsize=12,
             weight='bold', horizontalalignment='center',
             verticalalignment='center', color='black')
    plt.grid()
    if save_figure:
        if not os.path.exists('saved_figures/' + model_type + '_figures'):
            os.makedirs('saved_figures/' + model_type + '_figures')
        plt.savefig('saved_figures/' + model_type + '_figures/' + model_type
                    + '_temperature_timeseries_' + str(depth) + 'm.pdf')
    plt.show()
