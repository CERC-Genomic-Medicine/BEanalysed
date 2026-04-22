/* Shared utilities for the BEanalysed web ports.
 * Assumes SheetJS (XLSX) and Plotly have been loaded on the page. */

const BE = (function () {

  // ---------- Palettes (mirror matplotlib) ----------
  const PASTEL1 = [
    '#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6',
    '#ffffcc','#e5d8bd','#fddaec','#f2f2f2'
  ];
  const TAB20 = [
    '#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c',
    '#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5',
    '#8c564b','#c49c94','#e377c2','#f7b6d2','#7f7f7f',
    '#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5'
  ];
  // Shared consequence color/marker mapping (matches the Python scripts).
  const CONSEQUENCE_COLOR = {
    'synonymous': 'green',
    'missense':   'purple',
    'non-sense':  'red',
    'splice':     'gold',
    'regulatory': 'blue',
    'non-coding': 'lightgray',
    'No predicted Mutation': 'gray',
  };
  const CONSEQUENCE_MARKER = {
    'synonymous': 'diamond',
    'missense':   'circle',
    'non-sense':  'square',
    'splice':     'triangle-up',
    'regulatory': 'x',
    'non-coding': 'triangle-down',
    'No predicted Mutation': 'cross',
  };

  function getColorPalette(n) {
    return n <= PASTEL1.length ? PASTEL1 : TAB20;
  }

  // ---------- Debounce ----------
  function debounce(fn, ms = 80) {
    let t = null;
    return function (...args) {
      clearTimeout(t);
      t = setTimeout(() => fn.apply(this, args), ms);
    };
  }

  // ---------- Status bar ----------
  function Status(el) {
    return {
      info:  msg => { el.textContent = msg; el.className = 'status info';  el.classList.remove('hidden'); },
      warn:  msg => { el.textContent = msg; el.className = 'status warn';  el.classList.remove('hidden'); },
      error: msg => { el.textContent = msg; el.className = 'status error'; el.classList.remove('hidden'); },
      clear: ()  => { el.classList.add('hidden'); },
    };
  }

  // ---------- XLSX reading ----------
  async function readWorkbook(file) {
    const buf = await file.arrayBuffer();
    const wb = XLSX.read(buf, { type: 'array' });
    const sheets = {};
    for (const name of wb.SheetNames) {
      sheets[name] = XLSX.utils.sheet_to_json(wb.Sheets[name], { defval: null });
    }
    return { workbook: wb, sheets };
  }

  // ---------- TSV / CSV reading ----------
  async function readTable(file, sep) {
    const text = await file.text();
    sep = sep || (file.name.toLowerCase().endsWith('.csv') ? ',' : '\t');
    const lines = text.split(/\r?\n/).filter(l => l.length > 0);
    if (lines.length === 0) return [];
    const header = lines[0].split(sep);
    const rows = [];
    for (let i = 1; i < lines.length; i++) {
      const parts = lines[i].split(sep);
      const obj = {};
      for (let j = 0; j < header.length; j++) {
        const v = parts[j];
        const n = Number(v);
        obj[header[j]] = (v === '' || v == null) ? null
                       : !Number.isNaN(n) && v.trim() !== '' ? n
                       : v;
      }
      rows.push(obj);
    }
    return rows;
  }

  // ---------- Checklist UI ----------
  function renderChecklist(container, items, onChange, initiallyChecked) {
    container.innerHTML = '';
    const set = new Set(
      initiallyChecked === undefined ? items : initiallyChecked
    );
    for (const item of items) {
      const lbl = document.createElement('label');
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.value = item;
      cb.checked = set.has(item);
      cb.dataset.value = item;
      cb.addEventListener('change', onChange);
      lbl.appendChild(cb);
      lbl.appendChild(document.createTextNode(' ' + item));
      container.appendChild(lbl);
    }
  }
  function getChecked(container) {
    return [...container.querySelectorAll('input[type=checkbox]')]
      .filter(cb => cb.checked)
      .map(cb => cb.dataset.value);
  }

  // ---------- Stats ----------
  function rocCurve(scores, truths) {
    const pairs = scores.map((s, i) => [s, truths[i]]);
    pairs.sort((a, b) => b[0] - a[0]);
    const P = truths.reduce((a, v) => a + (v === 1 ? 1 : 0), 0);
    const N = truths.length - P;
    if (P === 0 || N === 0) return { fpr: [0, 1], tpr: [0, 1], auc: NaN };
    const fpr = [0], tpr = [0];
    let tp = 0, fp = 0, prev = Infinity;
    for (const [s, y] of pairs) {
      if (s !== prev) { fpr.push(fp / N); tpr.push(tp / P); prev = s; }
      if (y === 1) tp++; else fp++;
    }
    fpr.push(fp / N); tpr.push(tp / P);
    let auc = 0;
    for (let i = 1; i < fpr.length; i++) {
      auc += (fpr[i] - fpr[i - 1]) * (tpr[i] + tpr[i - 1]) / 2;
    }
    return { fpr, tpr, auc };
  }

  function pearsonR(x, y) {
    const n = x.length;
    let sx = 0, sy = 0;
    for (let i = 0; i < n; i++) { sx += x[i]; sy += y[i]; }
    const mx = sx / n, my = sy / n;
    let num = 0, dx = 0, dy = 0;
    for (let i = 0; i < n; i++) {
      const a = x[i] - mx, b = y[i] - my;
      num += a * b; dx += a * a; dy += b * b;
    }
    const d = Math.sqrt(dx * dy);
    return d === 0 ? NaN : num / d;
  }

  // Quantile (linear interpolation like numpy default).
  function quantile(arr, q) {
    const a = arr.filter(v => v != null && !Number.isNaN(v)).slice().sort((x, y) => x - y);
    if (a.length === 0) return NaN;
    const pos = (a.length - 1) * q;
    const lo = Math.floor(pos), hi = Math.ceil(pos);
    if (lo === hi) return a[lo];
    return a[lo] + (a[hi] - a[lo]) * (pos - lo);
  }

  // Simple LOWESS replacement: locally weighted running mean with tricube weights.
  // Good enough for visual smoothing.
  function lowess(x, y, frac = 0.2) {
    const idx = x.map((_, i) => i).sort((a, b) => x[a] - x[b]);
    const xs = idx.map(i => x[i]);
    const ys = idx.map(i => y[i]);
    const n = xs.length;
    const span = Math.max(2, Math.floor(frac * n));
    const out = new Array(n);
    for (let i = 0; i < n; i++) {
      const lo = Math.max(0, i - Math.floor(span / 2));
      const hi = Math.min(n, lo + span);
      let wsum = 0, vsum = 0;
      const maxd = Math.max(
        Math.abs(xs[i] - xs[lo]),
        Math.abs(xs[hi - 1] - xs[i]),
        1e-12
      );
      for (let j = lo; j < hi; j++) {
        const u = Math.abs(xs[j] - xs[i]) / maxd;
        const w = u >= 1 ? 0 : Math.pow(1 - u * u * u, 3);
        wsum += w; vsum += w * ys[j];
      }
      out[i] = wsum > 0 ? vsum / wsum : ys[i];
    }
    return { x: xs, y: out };
  }

  // ---------- Navigation bar ----------
  function installNav(currentHref) {
    const nav = document.createElement('div');
    nav.className = 'toolbar-nav';
    const pages = [
      ['index.html',          'Home'],
      ['BE_ROC.html',         'ROC'],
      ['BE_ROC_PerGene.html', 'ROC per Gene'],
      ['BE_lollipop.html',    'Lollipop'],
      ['BE_scatter.html',     'Scatter'],
      ['RepeatOnRepeat.html', 'Repeat'],
      ['test_ROC_perf.html',  'ROC Quantile Sweep'],
      ['structure_viewer.html','3D Viewer'],
    ];
    for (const [href, label] of pages) {
      const a = document.createElement('a');
      a.href = href;
      a.textContent = label;
      if (currentHref && href === currentHref) a.style.fontWeight = '700';
      nav.appendChild(a);
    }
    document.body.insertBefore(nav, document.body.firstChild);
  }

  // ---------- Window-resize → Plotly.Plots.resize ----------
  // Keeps every .plot div in sync when the viewport resizes (since the
  // .plot CSS uses viewport-relative heights).
  let _autoResizeInstalled = false;
  function installAutoResize() {
    if (_autoResizeInstalled) return;
    _autoResizeInstalled = true;
    const handler = debounce(() => {
      if (typeof Plotly === 'undefined') return;
      document.querySelectorAll('.plot').forEach(div => {
        if (div && div.data) Plotly.Plots.resize(div);
      });
    }, 80);
    window.addEventListener('resize', handler);
  }

  // ---------- Misc ----------
  function uniqueSorted(arr) {
    return [...new Set(arr.filter(v => v != null && v !== ''))].sort();
  }

  return {
    PASTEL1, TAB20, CONSEQUENCE_COLOR, CONSEQUENCE_MARKER,
    getColorPalette, debounce, Status,
    readWorkbook, readTable,
    renderChecklist, getChecked,
    rocCurve, pearsonR, quantile, lowess,
    installNav, installAutoResize, uniqueSorted,
  };
})();
