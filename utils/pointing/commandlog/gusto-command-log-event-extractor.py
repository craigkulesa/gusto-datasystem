#!/usr/bin/env python3
"""GUSTO command-log attitude-event extractor (v1.0).

Streams the ACE/MAX3 command log (~30 MB) and distills it into small TSV
files of attitude/pointing-relevant events that fit the chat upload limit:

  events.tsv      rare events: pointing-offset commands (tElD/tAzD...),
                  gyro-bias updates (estBias/updateGyroBias), interactive
                  star fits with RA/Dec diffs (fix()/estBias() output),
                  QUATERNION transitions, camera power/connection events,
                  NGSC parameter changes, DAC faults, unknown commands
  quickfix.tsv    every completed star fix (in-scan + interactive)
  bad_windows.tsv windows where the star camera was failing
                  (gyro-only drift windows)
  summary.txt     star-fix cadence stats, gap distribution, per-day counts,
                  DAng stats, interactive-fit RA/Dec error stats

Optional raw streams (big, usually NOT for upload):
  --posfix        posfix.tsv  every posFix line (DAng/DTm/DO-fix?)
  --all-commands  commands.tsv every ">>>" command instance

Usage:
  python3 gusto_cmdlog_extract.py command.txt -o cmd_events
  python3 gusto_cmdlog_extract.py command.txt -o cmd_events --posfix

Pure Python 3 stdlib only.  Timestamps in the log are UTC (verified:
"Jan 14 04:08:22" daemon lines interleave exactly with
"14-Jan-24 04:08:22 UTC" lines).  No-year "Mon DD" stamps get year 2023
for Dec, 2024 for Jan/Feb (GUSTO flew 2023-12 -> 2024-02).
"""

import argparse
import os
import re
import statistics
import sys
from collections import Counter
from datetime import datetime

MONTHS = {m: i + 1 for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])}

# 14-Jan-24 04:05:49  (explicit 2-digit year)
TS_EXPL = re.compile(r'(?<![\dA-Za-z-])(\d{1,2})-([A-Z][a-z]{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})')
# Jan 14 04:08:22  (no year)
TS_MON = re.compile(r'([A-Z][a-z]{2}) (\d{1,2}) (\d{2}):(\d{2}):(\d{2})(?!\d)')

# ---- line event patterns -------------------------------------------------
RE_CMD = re.compile(r'>>> ?([A-Za-z_][\w.]*)\s*(?:\(([^)]*)\))?')
RE_OFFSET_SET = re.compile(r'Setting (t\w+) to (-?[\d.eE+-]+) which equals (-?[\d.eE+-]+) arcsecs')
RE_SETTING_T = re.compile(r'^\s*Setting (t\w+)')
RE_DRIFT = re.compile(r'total body drift rotation \(([^)]*)\) over ([\d.]+) seconds')
RE_BIAS_P = re.compile(r'pitch bias estimate (-?[\d.eE+-]+)')
RE_BIAS_R = re.compile(r'roll bias estimate (-?[\d.eE+-]+)')
RE_BIAS_Y = re.compile(r'yaw bias estimate (-?[\d.eE+-]+)')
RE_FOG_P = re.compile(r'old pitchFogOffset = (-?[\d.eE+-]+)')
RE_FOG_R = re.compile(r'old rollFogOffset = (-?[\d.eE+-]+)')
RE_FOG_Y = re.compile(r'old yawFogOffset = (-?[\d.eE+-]+)')
RE_SUCCESS = re.compile(r'Success!\s+([\d.]+) stars recognized')
RE_ALIGNQ_VEC = re.compile(r'\[\s*(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s*\]')
RE_RMS = re.compile(r'RMS Fit Error:\s+([\d.eE+-]+)')
RE_RA_DIFF = re.compile(r'RA:.*Diff:\s+(-?[\d.eE+-]+)')
RE_DEC_DIFF = re.compile(r'Dec:.*Diff:\s+(-?[\d.eE+-]+)')
RE_PLATE = re.compile(r'Plate Scale = (-?[\d.eE+-]+)')
RE_METHOD = re.compile(r'Method:\s+(\S+)')
RE_STAR_GOOD = re.compile(r'GOOD star fix with ([\d.]+) stars recognized')
RE_STAR_BAD = re.compile(r'BAD star fix')
RE_BIAS_FAIL = re.compile(r"can't estimate bias")
RE_ATTCAL = re.compile(r'Attitude calibration with STARCAM')
RE_CAMCONN = re.compile(r'(Lost connection to MAX3|Connected to MAX3|NGSC recv Fail)')
RE_MEDPIX = re.compile(r'Med Pixel Target: ([\d.]+), Med Pixel: ([\d.]+), New Exposure: ([\d.]+)')
RE_POSFIX = re.compile(r'posFix: DAng = ([\d.]+) deg\s+DTm = ([\d.]+) min\s+(DO Star Fix!|NO Star Fix)')
RE_DAC_FAULT = re.compile(r'(DAC seems to be dead|received ERROR from DAC|'
                           r'timeout waiting for DAC reply|DAC not responding|'
                           r'DAC Command connection killed|DAC H\w* connection killed)')

ROUTINE_CMDS = {'obs2Do', 'obs2Kill', 'gotoAz', 'gotoRaDec', 'gotoEl', 'gotoGal',
                'ace', 'aceKill', 'version', 'exit', 'ngStarCamGet', 'ngPGet',
                'dacOn', 'dacOff', 'dacInit', 'scanKill', 'pointScanKill', 'gotoKill'}
ROUTINE_PREFIX = ('mike.', 'goto', 'kill', 'Kill')


def year_for(mon, args, last_expl_year):
    """Year for 'Mon DD' stamps without a year (GUSTO: Dec 2023 - Feb 2024)."""
    if mon == 12:
        return args.dec_year
    if mon in (1, 2):
        return args.dec_year + 1
    return last_expl_year or args.dec_year + 1


def stamp(dt):
    return dt.strftime('%Y-%m-%dT%H:%M:%SZ') if dt else '?'


def trim(s, n=120):
    s = ' '.join(s.split())
    return s[:n]


def pct(sorted_vals, q):
    if not sorted_vals:
        return None
    k = min(len(sorted_vals) - 1, max(0, int(round(q * (len(sorted_vals) - 1)))))
    return sorted_vals[k]


def fit_detail(f):
    d = 'n=%g' % f['n']
    for k, lab in (('rms', 'rms'), ('ra', 'ra_diff_H'), ('dec', 'dec_diff_deg'),
                   ('plate', 'plate'), ('method', 'method')):
        if k in f:
            d += ' %s=%s' % (lab, f[k])
    if 'q' in f:
        d += ' alignQ=[%s]' % f['q']
    return d


def main():
    ap = argparse.ArgumentParser(description='Extract attitude/pointing events from the GUSTO command log.')
    ap.add_argument('cmdlog', help='command log file (command.txt)')
    ap.add_argument('-o', '--outdir', default='cmd_events', help='output directory (default: cmd_events)')
    ap.add_argument('--posfix', action='store_true', help='also write raw posfix.tsv (every posFix line)')
    ap.add_argument('--all-commands', action='store_true',
                    help='also write commands.tsv (every >>> command instance)')
    ap.add_argument('--dec-year', type=int, default=2023,
                    help='December year of the mission (default 2023)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    f_ev = open(os.path.join(args.outdir, 'events.tsv'), 'w')
    f_qf = open(os.path.join(args.outdir, 'quickfix.tsv'), 'w')
    f_bw = open(os.path.join(args.outdir, 'bad_windows.tsv'), 'w')
    f_pf = open(os.path.join(args.outdir, 'posfix.tsv'), 'w') if args.posfix else None
    f_cmd = open(os.path.join(args.outdir, 'commands.tsv'), 'w') if args.all_commands else None

    def emit(dt, cls, detail):
        f_ev.write('%s\t%s\t%s\n' % (stamp(dt), cls, trim(detail, 240)))

    cur = None                 # last-seen timestamp
    last_expl_year = None
    fit = None                 # open interactive-fit record
    q_next = False             # expecting alignQ vector on next line
    bias = {}                  # pending pitch/roll/yaw bias estimate
    fog = {}                   # pending old fogOffset values
    last_mech = None           # last DAC-fault emit time (burst dedupe)
    last_fit_flush = None      # last completed interactive fit time
    line_no = 0

    cmd_counter = Counter()
    unk_seen = set()
    good_fixes = []            # (dt, n)
    bad_fixes = []             # dt
    attcals = []
    fits = []                  # interactive fit dicts
    posfix_rows = []           # (dt, dang, dtm, do_fix)
    class_counter = Counter()
    n_lines = 0
    t_first = t_last = None

    def _fit_fields(text):
        """Fill fields of an open interactive-fit record from a text segment."""
        nonlocal fit, q_next, last_fit_flush
        if q_next:
            q = RE_ALIGNQ_VEC.search(text)
            if q:
                fit['q'] = ' '.join(q.groups())
            q_next = False
        if 'alignQ' in text:
            q_next = True
        m = RE_RMS.search(text)
        if m:
            fit['rms'] = m.group(1)
        m = RE_RA_DIFF.search(text)
        if m:
            fit['ra'] = m.group(1)
        m = RE_DEC_DIFF.search(text)
        if m:
            fit['dec'] = m.group(1)
        m = RE_PLATE.search(text)
        if m:
            fit['plate'] = m.group(1)
        m = RE_METHOD.search(text)
        if m:
            fit['method'] = m.group(1)
            fits.append(fit)
            emit(fit['t'], 'FIT', fit_detail(fit))
            class_counter['FIT'] += 1
            last_fit_flush = fit['t']
            fit = None

    # ---- per-line event processing -------------------------------------
    # The log glues output text and timestamps on shared lines, e.g.
    #   "...stars recognized.Jan 14 04:40:03 (ssh)..."
    # so each line is split at its timestamps: text before a stamp belongs
    # to the previously-seen time, text after it to that stamp.
    def process_text(text, dt):
        nonlocal fit, q_next, bias, fog, last_mech, last_fit_flush

        # --- commands (>>>) ---
        for m in RE_CMD.finditer(text):
            name = m.group(1)
            argstr = m.group(2)
            if name in MONTHS:          # ">>> Jan 14 04:08:22" prompt line
                continue
            if name.isdigit():
                continue
            if m.end() < len(text) and text[m.end()] == ':':
                continue                # output like ">>> NGSC: Params Set"
            full = name + ('(' + argstr + ')' if argstr is not None else '')
            cmd_counter[name] += 1
            if f_cmd:
                f_cmd.write('%s\t%s\n' % (stamp(dt), full))
            if name.startswith('t') and re.match(r'^t[A-Z]', name):
                emit(dt, 'OFFCMD', full)
                class_counter['OFFCMD'] += 1
            elif name in ('estBias', 'updateGyroBias') or \
                    (name == 'msv' and argstr and 'fogOffset' in argstr):
                emit(dt, 'BIASCMD', full)
                class_counter['BIASCMD'] += 1
            elif name == 'fix':
                emit(dt, 'FIXCMD', full)
                class_counter['FIXCMD'] += 1
            elif name in ('mst', 'ngPSet', 'ngStarFixParams', 'msv'):
                emit(dt, 'CAMCFG', full)
                class_counter['CAMCFG'] += 1
            elif re.match(r'cams(On|Off)', name):
                emit(dt, 'CAMPWR', full)
                class_counter['CAMPWR'] += 1
            elif name in ROUTINE_CMDS or name.startswith(ROUTINE_PREFIX):
                pass
            else:
                if name not in unk_seen:
                    unk_seen.add(name)
                    emit(dt, 'UNKCMD', 'first: %s (line %d)' % (full, line_no))
                    class_counter['UNKCMD'] += 1

        # --- pointing offsets ---
        m = RE_OFFSET_SET.search(text)
        if m:
            emit(dt, 'OFFSET', '%s=%s (%s arcsec)' % (m.group(1), m.group(2), m.group(3)))
            class_counter['OFFSET'] += 1
        elif RE_SETTING_T.match(text):
            emit(dt, 'OFFRAW', text)
            class_counter['OFFRAW'] += 1

        # --- gyro bias estimates / fog offsets ---
        m = RE_BIAS_P.search(text)
        if m:
            bias['p'] = m.group(1)
        m = RE_BIAS_R.search(text)
        if m:
            bias['r'] = m.group(1)
        m = RE_BIAS_Y.search(text)
        if m:
            bias['y'] = m.group(1)
            emit(dt, 'BIASEST', 'pitch=%s roll=%s yaw=%s' %
                 (bias.get('p', '?'), bias.get('r', '?'), bias.get('y')))
            class_counter['BIASEST'] += 1
            bias = {}
        m = RE_FOG_P.search(text)
        if m:
            fog['p'] = m.group(1)
        m = RE_FOG_R.search(text)
        if m:
            fog['r'] = m.group(1)
        m = RE_FOG_Y.search(text)
        if m:
            fog['y'] = m.group(1)
            emit(dt, 'FOGOFF', 'old pitch=%s roll=%s yaw=%s' %
                 (fog.get('p', '?'), fog.get('r', '?'), fog.get('y')))
            class_counter['FOGOFF'] += 1
            fog = {}
        if RE_BIAS_FAIL.search(text):
            emit(dt, 'BIASFAIL', 'estBias could not estimate bias')
            class_counter['BIASFAIL'] += 1
        m = RE_DRIFT.search(text)
        if m:
            emit(dt, 'DRIFT', 'body drift (%s) over %s s' % (m.group(1), m.group(2)))
            class_counter['DRIFT'] += 1

        # --- interactive star-fix block (fix() / estBias()) ---
        m = RE_SUCCESS.search(text)
        if m:
            if fit and 'n' in fit:
                fits.append(fit)
                emit(fit['t'], 'FIT', fit_detail(fit))
                class_counter['FIT'] += 1
                last_fit_flush = fit['t']
            fit = {'t': dt, 'n': float(m.group(1)), 'line': line_no}
            q_next = False
            # rest of this segment may already contain fit fields
            _fit_fields(text)
        elif fit is not None:
            _fit_fields(text)

        # --- in-scan fixes ---
        m = RE_STAR_GOOD.search(text)
        if m:
            n = float(m.group(1))
            good_fixes.append((dt, n))
            src = 'interactive'
            if last_fit_flush is None or (dt and last_fit_flush and
                                         (dt - last_fit_flush).total_seconds() > 300):
                src = 'inscan'
            f_qf.write('%s\t%s\t%s\n' % (stamp(dt), ('%g' % n), src))
        if RE_STAR_BAD.search(text):
            bad_fixes.append(dt)
            f_qf.write('%s\tBAD\n' % stamp(dt))
        if RE_ATTCAL.search(text):
            attcals.append(dt)

        # --- camera connection / config / quaternion ---
        if RE_CAMCONN.search(text):
            emit(dt, 'CAMCONN', text)
            class_counter['CAMCONN'] += 1
        m = RE_MEDPIX.search(text)
        if m:
            emit(dt, 'CAMCFG', 'auto-expose: med target=%s med=%s exposure=%s us' %
                 (m.group(1), m.group(2), m.group(3)))
            class_counter['CAMCFG'] += 1
        if 'QUATERNION' in text.upper():
            emit(dt, 'QUAT', text)
            class_counter['QUAT'] += 1

        # --- DAC faults (dedupe bursts to one event / 30 s) ---
        if RE_DAC_FAULT.search(text):
            if last_mech is None or (dt and (dt - last_mech).total_seconds() > 30):
                emit(dt, 'MECH', text)
                class_counter['MECH'] += 1
            last_mech = dt

        # --- posFix lines ---
        m = RE_POSFIX.search(text)
        if m:
            posfix_rows.append((dt, float(m.group(1)), float(m.group(2)),
                                m.group(3).startswith('DO')))
            if f_pf:
                f_pf.write('%s\t%s\t%s\t%s\n' %
                           (stamp(dt), m.group(1), m.group(2),
                            'DO' if m.group(3).startswith('DO') else 'NO'))

    with open(args.cmdlog, 'r', errors='replace') as fh:
        for line in fh:
            line_no += 1
            n_lines += 1

            # --- timestamps in this line, with positions ---
            stamps = []
            for m in TS_EXPL.finditer(line):
                d, mon, yy, hh, mm, ss = m.groups()
                try:
                    dt = datetime(2000 + int(yy), MONTHS[mon],
                                  int(d), int(hh), int(mm), int(ss))
                except (KeyError, ValueError):
                    continue
                stamps.append((m.start(), dt, True))
            for m in TS_MON.finditer(line):
                mon, d, hh, mm, ss = m.groups()
                try:
                    year = year_for(MONTHS[mon], args, last_expl_year)
                    dt = datetime(year, MONTHS[mon],
                                  int(d), int(hh), int(mm), int(ss))
                except (KeyError, ValueError):
                    continue
                stamps.append((m.start(), dt, False))
            stamps.sort(key=lambda s: s[0])

            if not stamps:
                process_text(line, cur)
            else:
                # text before the first stamp belongs to the previous time
                process_text(line[:stamps[0][0]], cur)
                for i, (pos, dt, is_expl) in enumerate(stamps):
                    end = stamps[i + 1][0] if i + 1 < len(stamps) else len(line)
                    process_text(line[pos:end], dt)
                    if is_expl:
                        last_expl_year = dt.year
                cur = stamps[-1][1]
            if cur is not None:
                if t_first is None:
                    t_first = cur
                t_last = cur

    # flush trailing fit block
    if fit and 'n' in fit:
        fits.append(fit)
        emit(fit['t'], 'FIT', fit_detail(fit))
        class_counter['FIT'] += 1

    f_ev.close()
    f_qf.close()
    if f_pf:
        f_pf.close()
    if f_cmd:
        f_cmd.close()

    # ------------------------------------------------------------------
    # summary
    # ------------------------------------------------------------------
    sum_lines = []

    def S(s=''):
        sum_lines.append(s)

    S('GUSTO command-log extraction summary')
    S('input : %s' % args.cmdlog)
    S('lines : %d   time span: %s -> %s (UTC)' %
      (n_lines, stamp(t_first), stamp(t_last)))
    S()
    S('== event counts (events.tsv) ==')
    for k, v in sorted(class_counter.items()):
        S('  %-8s %d' % (k, v))
    S()
    S('== star-fix timeline ==')
    tl = sorted([t for t, _ in good_fixes] + [f['t'] for f in fits if f['t']])
    dedup = []
    for t in tl:
        if dedup and (t - dedup[-1]).total_seconds() <= 120:
            continue
        dedup.append(t)
    gaps = [(dedup[i + 1] - dedup[i]).total_seconds() / 60.0
            for i in range(len(dedup) - 1)]
    S('completed star fixes: %d in-scan GOOD + %d interactive (deduped %d)' %
      (len(good_fixes), len(fits), len(dedup)))
    S('BAD star fixes: %d ; "Attitude calibration with STARCAM": %d' %
      (len(bad_fixes), len(attcals)))
    if gaps:
        gs = sorted(gaps)
        S('minutes between consecutive fixes: median=%.1f mean=%.1f p90=%.1f max=%.1f' %
          (statistics.median(gs), sum(gs) / len(gs), pct(gs, 0.9), gs[-1]))
        S('gaps > 30 min: %d ; > 60 min: %d ; > 120 min: %d' %
          (sum(1 for g in gaps if g > 30), sum(1 for g in gaps if g > 60),
           sum(1 for g in gaps if g > 120)))
        S('longest gaps (drift-suspect windows):')
        big = sorted(range(len(gaps)), key=lambda i: -gaps[i])[:20]
        for i in sorted(big):
            S('  %s -> %s  = %6.1f min' % (stamp(dedup[i]), stamp(dedup[i + 1]), gaps[i]))
    S()
    S('== BAD-fix windows (gyro-only drift windows) ==')
    windows = []
    if bad_fixes:
        goods = sorted([t for t, _ in good_fixes] + [f['t'] for f in fits if f['t']])
        bsort = sorted(bad_fixes)
        import bisect
        i = 0
        while i < len(bsort):
            j = i
            while j + 1 < len(bsort):
                k = bisect.bisect_left(goods, bsort[j])
                nxt_good = goods[k] if k < len(goods) else None
                if nxt_good is not None and nxt_good < bsort[j + 1]:
                    break
                j += 1
            k = bisect.bisect_left(goods, bsort[j])
            end = goods[k] if k < len(goods) else None
            dur = ((end - bsort[i]).total_seconds() / 60.0) if end else None
            windows.append((bsort[i], end, dur, j - i + 1))
            i = j + 1
    for st, en, dur, nb in windows:
        f_bw.write('%s\t%s\t%s\t%d\n' % (stamp(st), stamp(en),
                   ('%.1f' % dur) if dur else 'open?', nb))
    f_bw.close()
    S('bad windows: %d (total %.1f h)' %
      (len(windows), sum(w[2] or 0.0 for w in windows) / 60.0))
    for st, en, dur, nb in sorted(windows, key=lambda w: -(w[2] or 0))[:20]:
        S('  %s -> %s   %s min  (%d bad fixes)' %
          (stamp(st), stamp(en), ('%.1f' % dur) if dur else 'open?', nb))
    S()
    S('== posFix stats (%d posFix lines) ==' % len(posfix_rows))
    if posfix_rows:
        dang = sorted(p[1] for p in posfix_rows)
        dofx = sum(1 for p in posfix_rows if p[3])
        S('DAng deg: median=%.2f p90=%.2f max=%.2f ; DO Star Fix: %d (%.1f%%)' %
          (statistics.median(dang), pct(dang, 0.9), dang[-1], dofx,
           100.0 * dofx / len(posfix_rows)))
    S()
    S('== interactive fits (fix()/estBias()): attitude error at fix time ==')
    S('n=%d ; RA Diff is in HOURS (x15->deg), Dec Diff in DEG' % len(fits))
    by_m = Counter(f.get('method', '?') for f in fits)
    S('methods: %s' % ', '.join('%s:%d' % kv for kv in sorted(by_m.items())))
    for key, lab, scale in (('ra', 'RA diff', 54000.0), ('dec', 'Dec diff', 3600.0)):
        vals = sorted(float(f[key]) * scale for f in fits if key in f)
        if vals:
            S('%s arcsec: median=%.1f p90=%.1f max_abs=%.1f' %
              (lab, statistics.median(vals), pct(vals, 0.9),
               max(abs(v) for v in vals)))
    rms_vals = sorted(float(f['rms']) for f in fits if 'rms' in f)
    if rms_vals:
        S('RMS fit error px: median=%.2f p90=%.2f' %
          (statistics.median(rms_vals), pct(rms_vals, 0.9)))
    S()
    S('== per-day counts (good fixes / bad fixes / attcal) ==')
    per_day = {}
    for t, _ in good_fixes:
        if t:
            per_day.setdefault(t.strftime('%Y-%m-%d'), [0, 0, 0])[0] += 1
    for t in bad_fixes:
        if t:
            per_day.setdefault(t.strftime('%Y-%m-%d'), [0, 0, 0])[1] += 1
    for t in attcals:
        if t:
            per_day.setdefault(t.strftime('%Y-%m-%d'), [0, 0, 0])[2] += 1
    for d in sorted(per_day):
        c = per_day[d]
        S('  %s  good=%-5d bad=%-4d attcal=%-5d' % (d, c[0], c[1], c[2]))
    S()
    S('== command frequency (all commands, incl. routine) ==')
    for name, cnt in cmd_counter.most_common(40):
        S('  %-24s %d' % (name, cnt))
    if len(cmd_counter) > 40:
        S('  ... %d distinct command names total' % len(cmd_counter))

    with open(os.path.join(args.outdir, 'summary.txt'), 'w') as f:
        f.write('\n'.join(sum_lines) + '\n')

    sizes = {}
    for fn in ('events.tsv', 'quickfix.tsv', 'bad_windows.tsv', 'summary.txt',
               'posfix.tsv', 'commands.tsv'):
        p = os.path.join(args.outdir, fn)
        if os.path.exists(p):
            sizes[fn] = os.path.getsize(p)
    print('wrote %d files to %s:' % (len(sizes), args.outdir))
    for fn, sz in sorted(sizes.items()):
        print('  %-16s %8d bytes %s' % (fn, sz, '<-- uploadable'
                                        if sz < 500000 else '(too big to upload)'))
    print('upload guidance: upload events.tsv, bad_windows.tsv, summary.txt'
          ' (and quickfix.tsv if marked uploadable)')


if __name__ == '__main__':
    main()