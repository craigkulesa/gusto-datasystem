#!/usr/bin/env python3
"""GUSTO UDP telemetry extractor (v1.0).

Parses the per-scanid UDP/telemetry binary files (udp_XXXXX.dat, ICD
7503-9035, 156-byte records) and distills them into small TSV files of
attitude/pointing metrics. Pure stdlib, Python 3. Run on the Mac where
the udp_*.dat files live; upload only the compact outputs.

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

Verified decode (udp_03047, 2024-01-15 22:18 UTC):
  - stored ra/dec = constant rotation C applied to rot(q_att, tel);
    C = 629.01" +- 0.01" in that file. C per file is a calibration
    constant to be tracked across the mission.
  - q_att is live at 0.4 s cadence; ra/dec jitter 2.1" 2D-rms;
    gondola rolls about the boresight (body-z) with a ~10-15 s torsional
    oscillation on top of a net drift (+342"/min in that file) -
    invisible to pointing position, visible to slit/camera.

Outputs (written to -o OUTDIR):
  udp_file_summary.tsv   one row per file (full detail; keep local)
  udp_intervals.tsv      one row per within-file fix interval (keep local)
  udp_fix_events.tsv     fixNum changes + q_adj changes (keep local)
  udp_hourly.tsv         hourly aggregates            <- UPLOAD (small)
  udp_regress_binned.tsv binned dq vs dt / path        <- UPLOAD (small)
  udp_report.txt         regression results + findings <- UPLOAD (small)

Usage:
  python3 gusto_udp_extract.py /path/to/udps -o udp_out
  python3 gusto_udp_extract.py udp_03047.dat -o udp_out          (one file)
  python3 gusto_udp_extract.py /path/to/udps -o udp_out --verbose-first
      --verbose-first prints the per-record table of the first file
"""

import argparse
import glob
import math
import os
import struct
import sys
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
            'rollrate', 'oscamp', 'nt_attT', 'attT_off', 'starc0', 'starc1')
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

    # --- stored pointing jitter about mean ---
    ra = [math.degrees(r[5]) for r in recs]
    de = [math.degrees(r[6]) for r in recs]
    mra = sum(ra) / len(ra)
    mde = sum(de) / len(de)
    cd = math.cos(math.radians(mde))
    sx = math.sqrt(sum(((r - mra) * 3600.0 * cd) ** 2 for r in ra) / len(ra))
    sy = math.sqrt(sum(((d - mde) * 3600.0) ** 2 for d in de) / len(de))
    row['jit'] = math.hypot(sx, sy)

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
        csub = offs[i0:i1 + 1]
        cm = sum(csub) / len(csub)
        csd = math.sqrt(sum((x - cm) ** 2 for x in csub) / len(csub))
        intervals.append((utc(recs[i0][0], recs[i0][1]), base, scanid,
                          dur, th * ax[0] * AS, th * ax[1] * AS, th * ax[2] * AS,
                          path_total * AS, path_point * AS, cm, csd,
                          recs[i0][7], recs[i1][7],
                          tel[i0][0] * AS, -tel[i0][1] * AS))

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
        cs = [r['cmean'] for r in rs]
        js = [r['jit'] for r in rs]
        rr = [r['rollrate'] for r in rs]
        out.append((key, len(rs), median(cs),
                    max(r['cstd'] for r in rs), median(js), median(rr),
                    rs[-1]['tazd'], rs[-1]['teld']))
    return out


def regress(intervals):
    """bias / scale-factor tests on fix intervals >= 120 s."""
    good = [iv for iv in intervals if iv[3] >= 120.0]
    lines = []
    bins = []
    if not good:
        return ['no intervals >= 120 s found'], []
    for ai, name in ((4, 'x'), (5, 'y'), (6, 'z')):
        ys = [iv[ai] for iv in good]
        sl, ic, r2 = ols([iv[3] / 60.0 for iv in good], ys)
        lines.append('%s-axis: net angle vs duration: slope %+.2f"/min '
                     '(R2 %.3f, n %d)' % (name, sl, r2, len(good)))
        sl2, ic2, r22 = ols([iv[7] for iv in good], ys)
        lines.append('%s-axis: net angle vs path:   slope %+.5f (R2 %.3f)'
                     % (name, sl2, r22))
    dC = [iv[10] - iv[9] for iv in good]
    lines.append('C change within intervals: median |dC| %.2f", p90 %.2f"'
                 % (median([abs(d) for d in dC]),
                    sorted(abs(d) for d in dC)[int(0.9 * len(dC))]))
    durs = sorted(iv[3] for iv in good)
    nb = 20
    edges = [durs[int(i * (len(durs) - 1) / nb)] for i in range(nb + 1)]
    for b in range(nb):
        lo, hi = edges[b], edges[b + 1]
        if b == nb - 1:
            grp = [iv for iv in good if lo <= iv[3] <= hi]
        else:
            grp = [iv for iv in good if lo <= iv[3] < hi]
        if len(grp) < 3:
            continue
        bins.append((lo, hi, len(grp),
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
            'rollrate', 'oscamp', 'nt_attT', 'attT_off', 'starc0', 'starc1']
    with open(os.path.join(args.outdir, 'udp_file_summary.tsv'), 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in all_rows:
            fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

    with open(os.path.join(args.outdir, 'udp_intervals.tsv'), 'w') as fh:
        fh.write('t_start\tfile\tscanid\tdur_s\tnetx_as\tnety_as\tnetz_as\t'
                 'path_total_as\tpath_point_as\tcmean\tcstd\tfixn0\tfixn1\t'
                 'tazd\tteld\n')
        for iv in all_intervals:
            fh.write('%s\t%s\t%s\t%.1f\t%.2f\t%.2f\t%.2f\t%.1f\t%.2f\t'
                     '%.2f\t%.3f\t%d\t%d\t%.1f\t%.1f\n'
                     % ((iso(iv[0]),) + iv[1:]))

    with open(os.path.join(args.outdir, 'udp_fix_events.tsv'), 'w') as fh:
        fh.write('utc\tfile\ttype\told\tnew\tstarCount\n')
        for e in all_events:
            fh.write('%s\t%s\t%s\t%s\t%s\t%s\n'
                     % (iso(e[0]), e[1], e[2], e[3], e[4], e[5]))

    with open(os.path.join(args.outdir, 'udp_hourly.tsv'), 'w') as fh:
        fh.write('hour\tnfiles\tmedC_as\tmaxCstd_as\tmedjit_as\t'
                 'medrollrate\ttAzD_as\ttElD_as\n')
        for h in hourly(all_rows):
            fh.write('%s\t%d\t%.1f\t%.3f\t%.2f\t%.0f\t%.0f\t%.0f\n' % h)

    lines, bins = regress(all_intervals)
    with open(os.path.join(args.outdir, 'udp_regress_binned.tsv'), 'w') as fh:
        fh.write('dur_lo_s\tdur_hi_s\tn\tmed_netx\tmed_nety\tmed_netz\tmedC\n')
        for b in bins:
            fh.write('%.0f\t%.0f\t%d\t%.2f\t%.2f\t%.2f\t%.1f\n' % b)

    with open(os.path.join(args.outdir, 'udp_report.txt'), 'w') as fh:
        fh.write('GUSTO UDP telemetry extraction report\n')
        fh.write('files processed: %d\n' % len(all_rows))
        ok = [r for r in all_rows if not r.get('err')]
        if ok:
            fh.write('intervals: %d\n' % len(all_intervals))
            fh.write('fix events: %d   adj events: %d\n'
                     % (sum(r.get('nfix') or 0 for r in ok),
                        sum(r.get('nadj') or 0 for r in ok)))
            cs = [r['cmean'] for r in ok if r.get('cmean') != '']
            fh.write('C (ra/dec vs attitude offset): median %.1f"  min %.1f  '
                     'max %.1f\n' % (median(cs), min(cs), max(cs)))
            unstable = [r for r in ok
                        if r.get('cstd') != '' and r['cstd'] > 0.5]
            fh.write('files with C_std > 0.5": %d\n' % len(unstable))
            js = [r['jit'] for r in ok if r.get('jit') != '']
            fh.write('pointing jitter: median %.2f"  p90 %.2f"\n'
                     % (median(js), sorted(js)[int(0.9 * len(js))]))
            rr = [r['rollrate'] for r in ok if r.get('rollrate') != '']
            fh.write('roll drift rate: median %.0f"/min  p90 |rate| %.0f"/min\n'
                     % (median(rr),
                        sorted(abs(x) for x in rr)[int(0.9 * len(rr))]))
            oa = [r['oscamp'] for r in ok if r.get('oscamp') != '']
            fh.write('roll oscillation amplitude: median %.0f"  p90 %.0f"\n'
                     % (median(oa), sorted(oa)[int(0.9 * len(oa))]))
            fh.write('q_attT changed within file (nonzero): %d files\n'
                     % sum(1 for r in ok if r.get('nt_attT')))
        fh.write('\nbias / scale-factor regression (fix intervals >= 120 s):\n')
        for ln in lines:
            fh.write('  %s\n' % ln)
        fh.write('\nuploads to chat: udp_hourly.tsv + udp_regress_binned.tsv'
                 ' + udp_report.txt (keep the rest local)\n')
    print('done: %d files -> %s' % (len(all_rows), args.outdir))


if __name__ == '__main__':
    main()