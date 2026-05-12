"""
Translation of R process_scan functions to Python.
Provides functions to assign charge states and aggregate deconvoluted molecular weights
from an LC-MS spectrum (input as Mass, Response pairs).

Usage:
    from process_scan import process_scan
    result = process_scan(test_scan_df, polarity='Negative', baseline=1000, mz_error=0.01, ...)

Dependencies: pandas, numpy
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Optional, Dict

def cut_mmw_list(mwlist, intlist, mw_window):
    mwlist = np.asarray(mwlist)
    intlist = np.asarray(intlist)
    N = len(mwlist)
    if N == 0:
        return {'id': np.array([], dtype=int), 'mw': np.array([])}

    f = 1
    mw_feature = np.zeros(N, dtype=int)
    mw_avg = np.zeros(N, dtype=float)
    t0 = 0
    for k in range(1, N):
        ttt = np.arange(t0, k)
        min_mw = np.min(mwlist[ttt])
        idx_best = np.argmax(intlist[ttt])
        best_mw = mwlist[ttt[idx_best]]
        max_mw = np.max(mwlist[ttt])
        if (mwlist[k] - min_mw > mw_window) or (mwlist[k] - max_mw > 0.2):
            mw_feature[t0:k] = f
            mw_avg[t0:k] = np.round(best_mw, 4)
            f += 1
            t0 = k
    ttt = np.arange(t0, N)
    mw_feature[ttt] = f
    idx_best = ttt[np.argmax(intlist[ttt])] if len(ttt) > 0 else 0
    mw_avg[ttt] = np.round(mwlist[idx_best], 4) if len(ttt) > 0 else 0.0
    return {'id': mw_feature, 'mw': mw_avg}


def cut_mmw_list1(mwlist, intlist, mw_gap, mw_window):
    mwlist = np.asarray(mwlist)
    intlist = np.asarray(intlist)
    N = len(mwlist)
    if N == 0:
        return {'id': np.array([], dtype=int), 'mw': np.array([])}

    f = 1
    mw_feature = np.zeros(N, dtype=int)
    mw_avg = np.zeros(N, dtype=float)
    t0 = 0
    for k in range(1, N):
        ttt = np.arange(t0, k)
        min_mw = np.min(mwlist[ttt])
        idx_best = np.argmax(intlist[ttt])
        best_mw = mwlist[ttt[idx_best]]
        max_mw = np.max(mwlist[ttt])
        ratio_int = intlist[k] / intlist[k - 1] if k > 0 and intlist[k - 1] != 0 else 1e6
        ratio_int0 = intlist[k - 1] / intlist[k - 2] if k > 1 and intlist[k - 2] != 0 else 2

        cond1 = (mwlist[k] - best_mw > mw_window / 2 * 1.1) and (ratio_int > 1.1 and ratio_int0 > 1.1)
        cond2 = (mwlist[k] - min_mw > mw_window * 1.1) and (ratio_int > 1.1 and ratio_int0 > 1.1)
        cond3 = (mwlist[k] - max_mw > mw_gap * 1.1) and (mwlist[k] >= 1800)
        cond4 = (mwlist[k] - max_mw > mw_gap * 2.2) and (mwlist[k] < 1800)

        if cond1 or cond2 or cond3 or cond4:
            mw_feature[t0:k] = f
            mw_avg[t0:k] = np.round(best_mw, 4)
            f += 1
            t0 = k

    ttt = np.arange(t0, N)
    mw_feature[ttt] = f
    idx_best = ttt[np.argmax(intlist[ttt])] if len(ttt) > 0 else 0
    mw_avg[ttt] = np.round(mwlist[idx_best], 4) if len(ttt) > 0 else 0.0
    return {'id': mw_feature, 'mw': mw_avg}


def process_scan_high_charge_bis(scan1: pd.DataFrame, ref_charge, mz_error: float) -> pd.DataFrame:
    if scan1 is None or len(scan1) == 0:
        return pd.DataFrame(columns=['Mass', 'Response', 'z'])

    scan1 = scan1.sort_values('Mass').reset_index(drop=True)
    exp_mz = scan1['Mass'].to_numpy()
    NS = len(exp_mz)
    ref_charge = np.array(ref_charge, dtype=int)
    ref_dis = 1.0 / ref_charge

    # Validate charge state by looking at -0.5, +0.5 before/after:

    rows = []
    for j in range(NS):

        # Find candidate peaks within ±0.5 m/z and also within ±3 indices to account for potential missing peaks or noise

        lookup_range1 = np.where((exp_mz >= exp_mz[j] - 0.5) & (exp_mz <= exp_mz[j] + 0.5))[0]
        lookup_range2 = np.arange(max(0, j - 3), min(j + 3, NS - 1) + 1)
        if lookup_range1.size:
            start = min(lookup_range1[0], lookup_range2[0])
            end = max(lookup_range1[-1], lookup_range2[-1])
        else:
            start, end = lookup_range2[0], lookup_range2[-1]
        mzl = exp_mz[start:end + 1]
        
        if mzl.size < 2:
            continue
        d = np.abs(mzl.reshape(-1, 1) - mzl.reshape(1, -1))
        
        # Remove redundant pairs by taking only the lower triangle of the distance matrix
        
        tri_idx = np.tril_indices(d.shape[0], k=-1)
        dist_mzl = d[tri_idx]

        if dist_mzl.size == 0:
            continue
        
        # Find the most probable charge state by comparing the distance to the reference distances for each charge state

        dev_charge = np.abs(dist_mzl.reshape(-1, 1) - ref_dis.reshape(1, -1))
        tmp_matched = (dev_charge <= mz_error).any(axis=1).astype(int)
        best_idx = np.argmin(dev_charge, axis=1)
        charge_mz = ref_charge[best_idx] * tmp_matched
        charge_mz = charge_mz[charge_mz > 0]
        if charge_mz.size == 0:
            continue
        unique, counts = np.unique(charge_mz, return_counts=True)
        mask = counts >= 3
        if not mask.any():
            continue
        charge_labels = unique[mask].astype(int)
        charge_count = counts[mask].astype(int)
        if len(charge_labels) >= 2:
            labels = list(charge_labels)
            counts_list = list(charge_count)
            groups = []
            group_counts = []
            while len(labels) > 0:
                base = labels[0]
                to_check = [base]
                to_remove = [0]
                for idx in range(1, len(labels)):
                    other = labels[idx]
                    if (other % base == 0) or (base % other == 0):
                        to_check.append(other)
                        to_remove.append(idx)
                groups.append(to_check)
                group_counts.append(sum(counts_list[i] for i in to_remove))
                for ii in sorted(to_remove, reverse=True):
                    labels.pop(ii)
                    counts_list.pop(ii)
            best_group = groups[int(np.argmax(group_counts))]
            best_charge = max(best_group)
        else:
            best_charge = int(charge_labels[0])
        rows.append({'Mass': float(exp_mz[j]), 'Response': float(scan1.at[j, 'Response']), 'z': int(best_charge)})

    if len(rows) == 0:
        return pd.DataFrame(columns=['Mass', 'Response', 'z'])
    return pd.DataFrame(rows)[['Mass', 'Response', 'z']]


def process_scan_low_charge(scan1: pd.DataFrame, ref_charge, mz_error: float, baseline: float) -> pd.DataFrame:
    if scan1 is None or len(scan1) == 0:
        return pd.DataFrame(columns=['Mass', 'Response', 'z'])
    scan1 = scan1.sort_values('Mass').reset_index(drop=True).copy()
    scan3 = scan1.copy()
    scan4 = scan1.copy()
    scan3['z'] = 0
    scan4['z'] = 0
    exp_mz = scan1['Mass'].to_numpy()
    N = len(exp_mz)

    if 2 in ref_charge and N > 3:
        d = np.abs(exp_mz.reshape(-1, 1) - exp_mz.reshape(1, -1))
        tri_idx = np.tril_indices(d.shape[0], k=-1)
        dist_mz = d[tri_idx]
        pairs = np.array(list(zip(tri_idx[0], tri_idx[1])))
        if pairs.size:
            mask_pairs = np.where(np.abs(dist_mz - 0.5) <= mz_error)[0]
            pairs = pairs[mask_pairs]
            for a, b in pairs:
                if a >= b:
                    continue
                tmp_range = list(range(a, b + 1))
                Segment = scan1.iloc[tmp_range]
                denom = float(Segment['Response'].iloc[-1])
                SGR = float(Segment['Response'].iloc[0]) / denom if denom != 0 else 0
                if SGR > 0.5 and len(Segment) < 6:
                    scan3.at[a, 'z'] = 2
                    scan3.at[b, 'z'] = 2

    if 1 in ref_charge and N > 3:
        d = np.abs(exp_mz.reshape(-1, 1) - exp_mz.reshape(1, -1))
        tri_idx = np.tril_indices(d.shape[0], k=-1)
        dist_mz = d[tri_idx]
        pairs = np.array(list(zip(tri_idx[0], tri_idx[1])))
        if pairs.size:
            mask_pairs = np.where(np.abs(dist_mz - 1.0) <= mz_error)[0]
            pairs = pairs[mask_pairs]
            for a, b in pairs:
                if a >= b:
                    continue
                tmp_range = list(range(a, b + 1))
                Segment = scan1.iloc[tmp_range]
                denom = float(Segment['Response'].iloc[-1])
                SGR = float(Segment['Response'].iloc[0]) / denom if denom != 0 else 0
                if SGR > 0.5 and len(Segment) < 6:
                    scan4.at[a, 'z'] = 1
                    scan4.at[b, 'z'] = 1

    scan3['z1'] = scan4['z']
    for i in range(len(scan3)):
        if scan3.at[i, 'z1'] == 1:
            scan3.at[i, 'z'] = 1
        if scan3.at[i, 'z1'] == 0 and scan3.at[i, 'z'] == 0 and scan3.at[i, 'Response'] > baseline * 10:
            scan3.at[i, 'z'] = 1

    return scan3[['Mass', 'Response', 'z']]


def process_aggregation(scan_processed: pd.DataFrame) -> pd.DataFrame:
    tmp_feature = scan_processed['tmp_feature'].to_numpy()
    all_features = np.unique(tmp_feature)
    rows = []
    for f in all_features:
        valid = np.where(tmp_feature == f)[0]
        tmp_scan = scan_processed.iloc[valid]
        mass_vals = []
        for m in tmp_scan['Mass'].to_numpy():
            if isinstance(m, str) and ':' in m:
                parts = [float(x) for x in m.split(':') if x != '']
                mass_vals.extend(parts)
            else:
                mass_vals.append(float(m))
        masslist = sorted(mass_vals)
        mass_str = ':'.join([str(round(x, 4)) for x in masslist]) if masslist else '0'
        newRes = float(tmp_scan['Response'].sum())
        z_vals = []
        for z in tmp_scan['z'].to_numpy():
            if isinstance(z, str) and ':' in z:
                parts = [float(x) for x in z.split(':') if x != '']
                z_vals.extend(parts)
            else:
                z_vals.append(float(z))
        zlist = sorted(list(set(z_vals)))
        z_str = ':'.join([str(int(float(x))) for x in zlist if not np.isnan(float(x))]) if zlist else '0'
        newMW = float(np.round(np.mean(tmp_scan['MW'].astype(float)), 3)) if len(tmp_scan) > 0 else 0.0
        rows.append({'Mass': mass_str, 'Response': newRes, 'z': z_str, 'MW': newMW})
    return pd.DataFrame(rows)[['Mass', 'Response', 'z', 'MW']]


def process_scan(test_scan=None, polarity='Positive', MSMS=False, baseline=100,
                 min_charge=3, max_charge=12, min_mz=500, max_mz=1500, min_mw=4000, max_mw=12000,
                 mz_error=0.02, mw_gap=1.1, mw_window=10) -> Dict[str, pd.DataFrame]:
    if min_charge > max_charge:
        min_charge = max_charge - 1
    ref_charge_high = list(range(max(3, min_charge), max_charge + 1)) if max_charge >= 3 else []
    ref_charge_low = list(range(min_charge, 3)) if min_charge < 3 else []

    if test_scan is None:
        scan0 = pd.DataFrame(columns=['Mass', 'Response'])
    else:
        scan0 = pd.DataFrame(test_scan)
        if scan0.shape[1] >= 2:
            scan0 = scan0.iloc[:, :2]
            scan0.columns = ['Mass', 'Response']
        else:
            scan0 = pd.DataFrame(columns=['Mass', 'Response'])

    if scan0.shape[0] <= 5:
        return {'scan_processed': scan0.assign(z=0, MW=0.0), 'scan_processed_aggregated': pd.DataFrame()}

    scan0['Mass'] = pd.to_numeric(scan0['Mass'], errors='coerce')
    scan0['Response'] = pd.to_numeric(scan0['Response'], errors='coerce')
    scan0 = scan0.dropna()
    scan0 = scan0[(scan0['Response'] > baseline) & (scan0['Mass'] >= min_mz) & (scan0['Mass'] <= max_mz)].copy()

    scan_processed = pd.DataFrame()
    scan_processed_aggregated = pd.DataFrame()

    if len(scan0) > 5 and len(ref_charge_high) > 0:
        scan1 = process_scan_high_charge_bis(scan0, ref_charge_high, mz_error)
        if scan1 is not None and not scan1.empty:
            scan_processed = pd.concat([scan_processed, scan1], ignore_index=True)
            scan0 = scan0[~scan0['Mass'].isin(scan_processed['Mass'])].copy()

    if len(scan0) > 5 and len(ref_charge_low) > 0:
        scan1 = process_scan_low_charge(scan0, ref_charge_low, mz_error, baseline)
        if scan1 is not None and not scan1.empty:
            scan_processed = pd.concat([scan_processed, scan1], ignore_index=True)

    if scan_processed is None or scan_processed.empty:
        return {'scan_processed': pd.DataFrame(), 'scan_processed_aggregated': pd.DataFrame()}

    if polarity == 'Positive':
        scan_processed['MW'] = scan_processed['Mass'] * scan_processed['z'] - scan_processed['z'] * 1.00726
    else:
        scan_processed['MW'] = scan_processed['Mass'] * scan_processed['z'] + scan_processed['z'] * 1.00726

    scan_processed = scan_processed.sort_values('MW').reset_index(drop=True)
    scan_processed = scan_processed[(scan_processed['MW'] >= min_mw) & (scan_processed['MW'] <= max_mw)].copy()

    if scan_processed.shape[0] <= 3:
        annotated = pd.DataFrame(test_scan)
        if annotated.shape[1] >= 2:
            annotated = annotated.iloc[:, :2]
            annotated.columns = ['Mass', 'Response']
        annotated['z'] = np.zeros(len(annotated), dtype=int)
        annotated['MW'] = np.zeros(len(annotated), dtype=float)
        return {'scan_processed': annotated, 'scan_processed_aggregated': pd.DataFrame()}

    mmw_cutted = cut_mmw_list(scan_processed['MW'].to_numpy(), scan_processed['Response'].to_numpy(), 1)
    scan_processed['tmp_feature'] = mmw_cutted['id']
    scan_processed['MW'] = mmw_cutted['mw']

    counts = pd.Series(scan_processed['tmp_feature']).value_counts()
    if MSMS:
        frequent_feature = counts[counts >= 1].index
    else:
        frequent_feature = counts[counts > 1].index

    scan_processed = scan_processed[scan_processed['tmp_feature'].isin(frequent_feature)].copy()

    if scan_processed.shape[0] <= 3:
        annotated = pd.DataFrame(test_scan)
        if annotated.shape[1] >= 2:
            annotated = annotated.iloc[:, :2]
            annotated.columns = ['Mass', 'Response']
        annotated['z'] = 0
        annotated['MW'] = 0
        return {'scan_processed': annotated, 'scan_processed_aggregated': pd.DataFrame()}

    scan_processed_aggregated = process_aggregation(scan_processed)

    def z_positive(zval):
        try:
            parts = [float(x) for x in str(zval).split(':') if x != '']
            return any([x > 0 for x in parts])
        except Exception:
            return False
    if not scan_processed_aggregated.empty:
        scan_processed_aggregated = scan_processed_aggregated[scan_processed_aggregated['z'].apply(z_positive)].copy()
        tmp = cut_mmw_list1(scan_processed_aggregated['MW'].to_numpy(), scan_processed_aggregated['Response'].to_numpy(), mw_gap, mw_window)
        scan_processed_aggregated['Envelop'] = tmp['id']

    scan_charged = scan_processed.sort_values('Mass')[['Mass', 'Response', 'z', 'MW']].copy()
    original = pd.DataFrame(test_scan)
    if original.shape[1] >= 2:
        original = original.iloc[:, :2]
        original.columns = ['Mass', 'Response']
    original['z'] = np.zeros(len(original), dtype=int)
    original['MW'] = np.zeros(len(original), dtype=float)
    mass_to_idx = {m: i for i, m in enumerate(original['Mass'].to_numpy())}
    for _, row in scan_charged.iterrows():
        m = row['Mass']
        if m in mass_to_idx:
            idx = mass_to_idx[m]
            original.at[idx, 'z'] = row['z']
            original.at[idx, 'MW'] = row['MW']

    return {'scan_processed': original, 'scan_processed_aggregated': scan_processed_aggregated}