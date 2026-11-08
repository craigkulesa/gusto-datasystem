
#!/usr/bin/env python3
"""GUSTO UDP telemetry extractor (v1.1).

Parses the per-scanid UDP/telemetry binary files (udp_XXXXX.dat, ICD
7503-9035, 156-byte records) and distills them into small TSV files of
attitude/pointing metrics. Pure stdlib, Python 3. Run on the machine
where the udp_*.dat files live; upload only the bundle.

Record format (little-endian, 156 B, '<2I5f2I15d'):
  0  tsec   uint32  Unix seconds        1  tnsec uint32 nanoseconds
  2  alt    f32 m    3  lon f32 rad      4  lat f32 rad
  5  ra     f32 rad  6  dec f32 rad
  7  fixNum uint32 (star fix number; >0 valid)
  8  starCount uint32
  9-12  q_adj  4xf64 adjustment quaternion (scalar-LAST x,y,z,w)
  13-16 q_attT 4xf64 attitudeTaken
  17-20 q_att  4xf64 live attitude
  21-23 tel    3xf64 = (tAzD, -tElD, ~1)/norm in rad  (verified vs cmd log)

v1.1 changes (after full-mission run showed raster scan motion pollutes
the interval regression, and that C is a geometric, position-dependent
quantity):
  - stare-only regression: files with jitter < 30" (pointing stares);
    regressions use NET POINTING MOTION vs stare duration (the quantity
    the sawtooth/de-drift model needs), plus body-axis components and a
    scale-factor test vs accumulated path.
  - per-file mean pointing position (from the records' own ra/dec), so C
    can be modelled geometrically without the dataLog.
  - new outputs: udp_C_grid.tsv (0.5-deg sky grid, all files),
    udp_stare_C.tsv (clean stare C values, deduped by position+offsets).
  - upload/udp_v11.tar.gz: gzipped bundle of the five small files -
    upload ONLY this (rename to .txt if the uploader complains).

Outputs (written to -o OUTDIR):
  udp_file_summary.tsv   one row per file (keep local)
  udp_intervals.tsv      one row per within-file interval (keep local)
  udp_fix_events.tsv     fixNum/q_adj changes (keep local)
  udp_hourly.tsv         hourly aggregates              (in bundle)
  udp_regress_binned.tsv stare-only binned drift         (in bundle)
  udp_C_grid.tsv         C on a 0.5-deg sky grid         (in bundle)
  udp_stare_C.tsv        clean stare C vs position       (in bundle)
  udp_report.txt         regression results + findings   (in bundle)
  upload/udp_v11.tar.gz  the five bundle files, gzipped

Usage:
  python3 gusto_udp_extract.py /path/to/udps -o udp_out
  python3 gusto_udp_extract.py udp_03047.dat -o udp_out          (one file)
  python3 gusto_udp_extract.py /path/to/udps -o udp_out --verbose-first
"""

import argparse
import glob
import math
import os
import struct
import sys
import tarfile
from datetime import datetime, timezone

FMT = '<2I5f2I15d'
REC = struct.calcsize(FMT)          # 156
AS = 206264.80624709636            # rad -> arcsec


# ---------------- quaternion helpers (scalar-last: x,y,z,w) ----------------

def qmul(a, b):
    x1, y1, z1, w1 = a
    x2, y2, z2, w2 = b
    return (w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2)


def qconj(q):
    return (-q[0], -q[1], -q[2], q[3])


def qang(q):
    w = max(-1.0, min(1.0, abs(q[3])))
    return 2.0 * math.acos(w)


def qbetween(q1, q2):
    return qang(qmul(qconj(q1), q2))


def qrot(q, v):
    r = qmul(qmul(q, (v[0], v[1], v[2], 0.0)), qconj(q))
    return (r[0], r[1], r[2])


def qaxis(q):
    """angle (rad) and unit axis (x,y,z) of quaternion."""
    n = math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])
    th = 2.0 * math.atan2(n, abs(q[3]))
    if n < 1e-15:
        return th, (0.0, 0.0, 0.0)
    return th, (q[0] / n, q[1] / n, q[2] / n)


# ---------------- parsing ----------------

def parse_file(path):
    with open(path, 'rb') as fh:
        data = fh.read()
    if len(data) % REC:
        sys.stderr.write('WARNING %s: %d trailing bytes ignored\n'
                         % (os.path.basename(path), len(data) % REC))
    n = len(data) // REC
    return [struct.unpack_from(FMT, data, i * REC) for i in range(n)]


def utc(tsec, tnsec):
    return datetime.fromtimestamp(tsec + tnsec / 1e9, timezone.utc)


def iso(t):
    return t.strftime('%Y-%m-%dT%H:%M:%S')


# ---------------- per-file analysis ----------------

def analyze(path, recs, verbose=False):
    """Return (file_row, events, intervals)."""
    base = os.path.basename(path)
    scanid = ''
    for ch in base:
        if ch.isdigit():
            scanid += ch
        elif scanid:
            break
    row = {'file': base, 'scanid': scanid, 'nrec': len(recs), 'err': ''}
    events = []
    intervals = []
    keys = ('span_s', 'fixn0', 'fixn1', 'nfix', 'nadj', 'tazd', 'teld',
            'cmean', 'cstd', 'jit', 'netx', 'nety', 'netz',
            'rollrate', 'oscamp', 'nt_attT', 'attT_off', 'starc0', 'starc1',
            'mra', 'mde')
    if len(recs) < 2:
        row['err'] = 'nrec<2'
        for k in keys:
            row[k] = ''
        return row, events, intervals

    t0 = utc(recs[0][0], recs[0][1])
    t1 = utc(recs[-1][0], recs[-1][1])
    span = (t1 - t0).total_seconds()
    row['span_s'] = span
    row['fixn0'] = recs[0][7]
    row['fixn1'] = recs[-1][7]
    row['starc0'] = recs[0][8]
    row['starc1'] = recs[-1][8]

    q_att = [r[17:21] for r in recs]
    q_adj = [r[9:13] for r in recs]
    q_attT = [r[13:17] for r in recs]
    tel = [r[21:24] for r in recs]

    # --- fix / adj events ---
    nfix = nadj = 0
    for i in range(1, len(recs)):
        if recs[i][7] != recs[i - 1][7]:
            nfix += 1
            events.append((utc(recs[i][0], recs[i][1]), base, 'FIX',
                          '%d->%d' % (recs[i - 1][7], recs[i][7]), '',
                          recs[i][8]))
        if q_adj[i] != q_adj[i - 1]:
            a0 = qang(q_adj[i - 1]) * AS
            a1 = qang(q_adj[i]) * AS
            nadj += 1
            events.append((utc(recs[i][0], recs[i][1]), base, 'ADJ',
                          '%.2f' % a0, '%.2f' % a1, ''))
    row['nfix'] = nfix
    row['nadj'] = nadj

    # --- tel -> commanded offsets ---
    row['tazd'] = tel[0][0] * AS
    row['teld'] = -tel[0][1] * AS

    # --- C: rotation between rot(q_att,tel) and stored ra/dec ---
    offs = []
    for r, q, tv in zip(recs, q_att, tel):
        v = qrot(q, tv)
        c = (math.cos(r[5]) * math.cos(r[6]),
             math.sin(r[5]) * math.cos(r[6]),
             math.sin(r[6]))
        offs.append(math.acos(max(-1.0, min(1.0,
                v[0] * c[0] + v[1] * c[1] + v[2] * c[2]))) * AS)
    cmean = sum(offs) / len(offs)
    cstd = math.sqrt(sum((x - cmean) ** 2 for x in offs) / len(offs))
    row['cmean'] = cmean
    row['cstd'] = cstd

    # --- stored pointing: jitter about mean + mean position ---
    ra = [math.degrees(r[5]) for r in recs]
    de = [math.degrees(r[6]) for r in recs]
    mra = sum(ra) / len(ra)
    mde = sum(de) / len(de)
    cd = math.cos(math.radians(mde))
    sx = math.sqrt(sum(((r - mra) * 3600.0 * cd) ** 2 for r in ra) / len(ra))
    sy = math.sqrt(sum(((d - mde) * 3600.0) ** 2 for d in de) / len(de))
    row['jit'] = math.hypot(sx, sy)
    row['mra'] = mra
    row['mde'] = mde

    # --- net rotation & roll drift (body axes) ---
    dq = qmul(qconj(q_att[0]), q_att[-1])
    th, ax = qaxis(dq)
    row['netx'] = th * ax[0] * AS
    row['nety'] = th * ax[1] * AS
    row['netz'] = th * ax[2] * AS
    row['rollrate'] = row['netz'] / (span / 60.0)

    # --- boresight-roll oscillation: cumulative z series ---
    cum = 0.0
    lo = hi = 0.0
    for i in range(len(recs) - 1):
        d = qmul(qconj(q_att[i]), q_att[i + 1])
        t2, a2 = qaxis(d)
        cum += t2 * a2[2] * AS
        lo = min(lo, cum)
        hi = max(hi, cum)
    row['oscamp'] = (hi - lo) / 2.0

    # --- q_attT behaviour ---
    row['nt_attT'] = sum(1 for i in range(1, len(recs))
                         if q_attT[i] != q_attT[i - 1])
    row['attT_off'] = math.degrees(qbetween(q_attT[0], q_att[0]))

    # --- intervals: split at fixNum or q_adj changes ---
    bounds = [0] + [i for i in range(1, len(recs))
                    if recs[i][7] != recs[i - 1][7] or q_adj[i] != q_adj[i - 1]]
    if bounds[-1] != len(recs) - 1:
        bounds.append(len(recs) - 1)
    for b in range(len(bounds) - 1):
        i0, i1 = bounds[b], bounds[b + 1]
        dur = (utc(recs[i1][0], recs[i1][1])
               - utc(recs[i0][0], recs[i0][1])).total_seconds()
        if dur <= 0:
            continue
        dqn = qmul(qconj(q_att[i0]), q_att[i1])
        th, ax = qaxis(dqn)
        path_total = 0.0
        path_point = 0.0
        p_prev = qrot(q_att[i0], tel[i0])
        for i in range(i0, i1):
            d = qmul(qconj(q_att[i]), q_att[i + 1])
            path_total += qang(d)
            p = qrot(q_att[i + 1], tel[i + 1])
            dot = max(-1.0, min(1.0, p[0] * p_prev[0] + p[1] * p_prev[1]
                                + p[2] * p_prev[2]))
            path_point += math.acos(dot)
            p_prev = p
        # net pointing motion (angle between pointing dir at start and end)
        p0 = qrot(q_att[i0], tel[i0])
        p1 = qrot(q_att[i1], tel[i1])
        dot = max(-1.0, min(1.0, p0[0] * p1[0] + p0[1] * p1[1]
                            + p0[2] * p1[2]))
        net_point = math.acos(dot) * AS
        csub = offs[i0:i1 + 1]
        cm = sum(csub) / len(csub)
        csd = math.sqrt(sum((x - cm) ** 2 for x in csub) / len(csub))
        intervals.append((utc(recs[i0][0], recs[i0][1]), base, scanid,
                          dur, th * ax[0] * AS, th * ax[1] * AS, th * ax[2] * AS,
                          path_total * AS, path_point * AS, cm, csd,
                          recs[i0][7], recs[i1][7],
                          tel[i0][0] * AS, -tel[i0][1] * AS,
                          row['jit'], mra, mde, net_point))

    if verbose:
        print('verbose: %s  %s .. %s  C=%.2f+-%.2f"  jitter=%.2f"  '
              'rollrate=%+.0f"/min  oscamp=%.0f"'
              % (base, iso(t0), iso(t1), cmean, cstd, row['jit'],
                 row['rollrate'], row['oscamp']))
        for i, r in enumerate(recs):
            a = qbetween(q_att[i - 1], q_att[i]) * AS if i else 0.0
            print('  %s  ra %10.5f dec %10.5f  step %6.2f" fixNum %d '
                  'starCount %d'
                  % (utc(r[0], r[1]).strftime('%H:%M:%S.%f')[:12],
                     math.degrees(r[5]), math.degrees(r[6]), a, r[7], r[8]))
    return row, events, intervals


# ---------------- aggregation & regression ----------------

def median(v):
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def ols(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0:
        return 0.0, 0.0, 0.0
    slope = sxy / sxx
    icpt = my - slope * mx
    ssr = sum((y - (icpt + slope * x)) ** 2 for x, y in zip(xs, ys))
    sst = sum((y - my) ** 2 for y in ys)
    return slope, icpt, (1.0 - ssr / sst if sst > 0 else 0.0)


def hourly(rows):
    buckets = {}
    for r in rows:
        if not r.get('t0') or r.get('cmean') == '':
            continue
        buckets.setdefault(r['t0'][:13], []).append(r)
    out = []
    for key in sorted(buckets):
        rs = buckets[key]
        out.append((key, len(rs),
                    median([r['cmean'] for r in rs]),
                    max(r['cstd'] for r in rs),
                    median([r['jit'] for r in rs]),
                    median([r['rollrate'] for r in rs]),
                    rs[-1]['tazd'], rs[-1]['teld']))
    return out


def regress(intervals):
    """v1.1: STARE-ONLY bias / scale-factor regression.

    interval fields: 0 t, 1 file, 2 scanid, 3 dur_s, 4 netx, 5 nety,
    6 netz, 7 path_total, 8 path_point, 9 cmean, 10 cstd, 11 fixn0,
    12 fixn1, 13 tazd, 14 teld, 15 jit, 16 mra, 17 mde, 18 net_point
    """
    stare = [iv for iv in intervals if iv[15] < 30.0 and iv[3] >= 10.0]
    lines = []
    bins = []
    if not stare:
        return ['no stare intervals (jit<30", dur>=10s) found'], []
    ys = [iv[18] for iv in stare]
    xs = [iv[3] / 60.0 for iv in stare]
    sl, ic, r2 = ols(xs, ys)
    lines.append('POINTING DRIFT: net pointing motion vs stare duration: '
                 '%+.2f"/min (R2 %.3f, n %d)' % (sl, r2, len(stare)))
    med_durs = sorted(iv[3] for iv in stare)
    lines.append('  stare durations: median %.0f s  p90 %.0f s'
                 % (med_durs[len(med_durs) // 2],
                    med_durs[int(0.9 * len(med_durs))]))
    for ai, name in ((4, 'x'), (5, 'y'), (6, 'z')):
        sl, ic, r2 = ols(xs, [iv[ai] for iv in stare])
        lines.append('  body-%s drift: %+.2f"/min (R2 %.3f)' % (name, sl, r2))
    # scale-factor test: net pointing motion vs accumulated pointing path
    sl2, ic2, r22 = ols([iv[8] for iv in stare], ys)
    lines.append('SCALE FACTOR: net vs accumulated pointing path: '
                 '%+.4f (R2 %.3f)' % (sl2, r22))
    sl3, ic3, r23 = ols([iv[7] for iv in stare], ys)
    lines.append('  net vs total attitude path (roll-dominated): '
                 '%+.4f (R2 %.3f)' % (sl3, r23))
    # binned by duration
    nb = 20
    edges = [med_durs[int(i * (len(med_durs) - 1) / nb)] for i in range(nb + 1)]
    for b in range(nb):
        lo, hi = edges[b], edges[b + 1]
        if b == nb - 1:
            grp = [iv for iv in stare if lo <= iv[3] <= hi]
        else:
            grp = [iv for iv in stare if lo <= iv[3] < hi]
        if len(grp) < 5:
            continue
        bins.append((lo, hi, len(grp),
                     median([iv[18] for iv in grp]),
                     median([iv[4] for iv in grp]), median([iv[5] for iv in grp]),
                     median([iv[6] for iv in grp]), median([iv[9] for iv in grp])))
    return lines, bins


# ---------------- main ----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('paths', nargs='+', help='udp dir or files')
    ap.add_argument('-o', '--outdir', required=True)
    ap.add_argument('--verbose-first', action='store_true')
    ap.add_argument('--max-files', type=int, default=0)
    args = ap.parse_args()

    files = []
    for p in args.paths:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, 'udp_*.dat')))
        else:
            files.append(p)
    if args.max_files:
        files = files[:args.max_files]
    os.makedirs(args.outdir, exist_ok=True)

    all_rows, all_events, all_intervals = [], [], []
    for k, path in enumerate(files):
        try:
            recs = parse_file(path)
        except Exception as e:
            sys.stderr.write('ERROR %s: %s\n' % (path, e))
            continue
        row, events, intervals = analyze(
            path, recs, verbose=(k == 0 and args.verbose_first))
        row['t0'] = iso(utc(recs[0][0], recs[0][1])) if recs else ''
        all_rows.append(row)
        all_events.extend(events)
        all_intervals.extend(intervals)
        if (k + 1) % 500 == 0:
            sys.stderr.write('  %d/%d files\n' % (k + 1, len(files)))

    cols = ['scanid', 'file', 't0', 'nrec', 'span_s', 'err',
            'fixn0', 'fixn1', 'nfix', 'nadj', 'tazd', 'teld',
            'cmean', 'cstd', 'jit', 'netx', 'nety', 'netz',
            'rollrate', 'oscamp', 'nt_attT', 'attT_off', 'starc0', 'starc1',
            'mra', 'mde']
    with open(os.path.join(args.outdir, 'udp_file_summary.tsv'), 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in all_rows:
            fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

    with open(os.path.join(args.outdir, 'udp_intervals.tsv'), 'w') as fh:
        fh.write('t_start\tfile\tscanid\tdur_s\tnetx_as\tnety_as\tnetz_as\t'
                 'path_total_as\tpath_point_as\tcmean\tcstd\tfixn0\tfixn1\t'
                 'tazd\tteld\tjit\tmra\tmde\tnet_point_as\n')
        for iv in all_intervals:
            fh.write('%s\t%s\t%s\t%.1f\t%.2f\t%.2f\t%.2f\t%.1f\t%.2f\t'
                     '%.2f\t%.3f\t%d\t%d\t%.1f\t%.1f\t%.2f\t%.3f\t%.3f\t'
                     '%.2f\n' % ((iso(iv[0]),) + iv[1:]))

    with open(os.path.join(args.outdir, 'udp_fix_events.tsv'), 'w') as fh:
        fh.write('utc\tfile\ttype\told\tnew\tstarCount\n')
        for e in all_events:
            fh.write('%s\t%s\t%s\t%s\t%s\t%s\n'
                     % (iso(e[0]), e[1], e[2], e[3], e[4], e[5]))

    # hourly aggregate (bundle)
    with open(os.path.join(args.outdir, 'udp_hourly.tsv'), 'w') as fh:
        fh.write('hour\tnfiles\tmedC_as\tmaxCstd_as\tmedjit_as\t'
                 'medrollrate\ttAzD_as\ttElD_as\n')
        for h in hourly(all_rows):
            fh.write('%s\t%d\t%.1f\t%.3f\t%.2f\t%.0f\t%.0f\t%.0f\n' % h)

    # C on a 0.5-deg sky grid, all files (bundle)
    grid = {}
    for r in all_rows:
        if r.get('cmean') == '' or r.get('mra') == '':
            continue
        key = (round(r['mra'] * 2) / 2.0, round(r['mde'] * 2) / 2.0)
        g = grid.setdefault(key, {'C': [], 'cs': []})
        g['C'].append(r['cmean'])
        g['cs'].append(r['cstd'])
        g['tazd'], g['teld'], g['n'] = r['tazd'], r['teld'], \
            g.get('n', 0) + 1
    with open(os.path.join(args.outdir, 'udp_C_grid.tsv'), 'w') as fh:
        fh.write('ra_bin\tdec_bin\tn\tmedC\tminC\tmaxC\tmedcstd\t'
                 'tAzD_as\ttElD_as\n')
        for k in sorted(grid):
            g = grid[k]
            fh.write('%.1f\t%.1f\t%d\t%.1f\t%.1f\t%.1f\t%.3f\t%.0f\t%.0f\n'
                     % (k[0], k[1], g['n'], median(g['C']), min(g['C']),
                        max(g['C']), median(g['cs']), g['tazd'], g['teld']))

    # clean stare C vs position, deduped by 0.1-deg position + offsets (bundle)
    sg = {}
    for r in all_rows:
        if (r.get('err') or r.get('cmean') == '' or r.get('jit') == ''
                or r['jit'] >= 30.0 or r.get('cstd') == ''
                or r['cstd'] >= 0.5):
            continue
        key = (round(r['mra'] * 10) / 10.0, round(r['mde'] * 10) / 10.0,
               round(r['tazd']), round(r['teld']))
        g = sg.setdefault(key, {'C': []})
        g['C'].append(r['cmean'])
        g['cstd'], g['t0'] = r['cstd'], r['t0']
    with open(os.path.join(args.outdir, 'udp_stare_C.tsv'), 'w') as fh:
        fh.write('ra\tdec\tn\tmedC\tminC\tmaxC\tcstd\ttAzD_as\ttElD_as\t'
                 't_last\n')
        for k in sorted(sg):
            g = sg[k]
            fh.write('%.1f\t%.1f\t%d\t%.1f\t%.1f\t%.1f\t%.3f\t%.0f\t%.0f\t'
                     '%s\n' % (k[0], k[1], len(g['C']), median(g['C']),
                               min(g['C']), max(g['C']), g['cstd'],
                               k[2], k[3], g['t0']))

    # stare-only regression (bundle)
    lines, bins = regress(all_intervals)
    with open(os.path.join(args.outdir, 'udp_regress_binned.tsv'), 'w') as fh:
        fh.write('dur_lo_s\tdur_hi_s\tn\tmed_netpoint\tmed_netx\t'
                 'med_nety\tmed_netz\tmedC\n')
        for b in bins:
            fh.write('%.0f\t%.0f\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\n' % b)

    # report (bundle)
    with open(os.path.join(args.outdir, 'udp_report.txt'), 'w') as fh:
        fh.write('GUSTO UDP telemetry extraction report (v1.1)\n')
        fh.write('files processed: %d\n' % len(all_rows))
        ok = [r for r in all_rows if not r.get('err')]
        if ok:
            fh.write('intervals: %d   fix events: %d   adj events: %d\n'
                     % (len(all_intervals),
                        sum(r.get('nfix') or 0 for r in ok),
                        sum(r.get('nadj') or 0 for r in ok)))
            cs = [r['cmean'] for r in ok if r.get('cmean') != '']
            fh.write('C: median %.1f"  min %.1f  max %.1f\n'
                     % (median(cs), min(cs), max(cs)))
            stares = [r for r in ok if r.get('jit') != '' and r['jit'] < 30.0]
            clean = [r for r in stares
                     if r.get('cstd') != '' and r['cstd'] < 0.5]
            fh.write('stare files (jit<30"): %d   clean stares (cstd<0.5"): '
                     '%d\n' % (len(stares), len(clean)))
            js = [r['jit'] for r in stares]
            if js:
                fh.write('stare jitter: median %.2f"  p90 %.2f"\n'
                         % (median(js), sorted(js)[int(0.9 * len(js))]))
            fh.write('C sky grid cells populated: %d\n' % len(grid))
        fh.write('\nSTARE-ONLY bias / scale-factor regression:\n')
        for ln in lines:
            fh.write('  %s\n' % ln)
        fh.write('\nUPLOAD upload/udp_v11.tar.gz (rename .txt if needed) -'
                 ' nothing else\n')
    print('done: %d files -> %s' % (len(all_rows), args.outdir))

    # gzipped upload bundle
    updir = os.path.join(args.outdir, 'upload')
    os.makedirs(updir, exist_ok=True)
    bundle = os.path.join(updir, 'udp_v11.tar.gz')
    with tarfile.open(bundle, 'w:gz') as tar:
        for fn in ('udp_report.txt', 'udp_regress_binned.tsv',
                   'udp_C_grid.tsv', 'udp_stare_C.tsv', 'udp_hourly.tsv'):
            tar.add(os.path.join(args.outdir, fn), arcname=fn)
    print('upload bundle: %s (%d bytes)'
          % (bundle, os.path.getsize(bundle)))


if __name__ == '__main__':
    main()
