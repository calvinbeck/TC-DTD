#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Equations:
"""

import numpy as np

def d2Tdx2(T_xm, T_x, T_xp, Delta_x):
    h_x = Delta_x[2] - Delta_x[1]
    h_xm = Delta_x[1] - Delta_x[0]
    numerator = 2 * (T_xp + h_x / h_xm * T_xm - (1 + h_x / h_xm) * T_x)
    denominator = h_x * h_xm * (1 + h_x / h_xm)
    return numerator/denominator


def dTdt(T_tm, T_tp, Delta_t):
    return (T_tp-T_tm)/(2*Delta_t)