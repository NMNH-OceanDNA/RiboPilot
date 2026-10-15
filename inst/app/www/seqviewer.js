// inst/app/www/seqviewer.js
// Sequence viewer for the assemble and annotate details windows, ported from
// MitoPilot's viewer. Units are always linear.
(function () {
  'use strict';

  // ---- geometry: pure functions, testable without a canvas ----
  function span(f) { return f.pos2 - f.pos1 + 1; }
  // Greedy interval packing: each feature takes the first lane it fits in.
  function lanes(feats) {
    var order = feats.slice().sort(function (a, b) { return a.pos1 - b.pos1; });
    var laneEnds = [];
    order.forEach(function (f) {
      var lane = -1;
      for (var i = 0; i < laneEnds.length; i++) if (laneEnds[i] < f.pos1) { lane = i; break; }
      if (lane < 0) { lane = laneEnds.length; laneEnds.push(0); }
      laneEnds[lane] = Math.max(laneEnds[lane], f.pos2);
      f.lane = lane;
    });
    return laneEnds.length;
  }
  // Genomic side made partial by strand + partial5/partial3: 'start' (pos1
  // end), 'end' (pos2 end), 'both', or null when both ends are real.
  function partialEdge(f) {
    var startOpen = f.dir === '-' ? f.partial3 : f.partial5;
    var endOpen = f.dir === '-' ? f.partial5 : f.partial3;
    if (startOpen && endOpen) return 'both';
    if (startOpen) return 'start';
    if (endOpen) return 'end';
    return null;
  }
  // Whether a feature shows its own pos1 / pos2 edge, vs. being cut off by the
  // visible range [vs, ve).
  function segOwns(f, vs, ve) { return { start: f.pos1 >= vs, end: f.pos2 <= ve }; }
  // Columns where the rows covering them disagree; internal gaps count,
  // terminal gaps mean the row does not cover the column.
  function diffColumns(seqs) {
    var out = [], n = seqs.length ? seqs[0].length : 0;
    var sp = seqs.map(function (s) { return [s.search(/[^-]/), s.replace(/-+$/, '').length - 1]; });
    for (var c = 0; c < n; c++) {
      var ch = null, d = false;
      for (var r = 0; r < seqs.length; r++) {
        if (c < sp[r][0] || c > sp[r][1]) continue;
        var x = seqs[r].charAt(c);
        if (ch === null) ch = x; else if (x !== ch) d = true;
      }
      if (d) out.push(c + 1);
    }
    return out;
  }
  // Native position of column `col` in a gapped string: count of non-gap
  // characters up to and including it, or null when the column is a gap.
  function nativePos(seq, col) {
    if (seq.charAt(col - 1) === '-') return null;
    var n = 0;
    for (var i = 0; i < col; i++) if (seq.charAt(i) !== '-') n++;
    return n;
  }

  window.rpseq = window.rpseq || {};
  window.rpseq.geom = { span: span, lanes: lanes, partialEdge: partialEdge, segOwns: segOwns,
                        diffColumns: diffColumns, nativePos: nativePos };
})();

(function () {
  'use strict';
  var G = window.rpseq.geom;
  var MAX_PPB = 14, NT_LETTER = 8, NT_BAR = 3;
  var RULER_H = 22, LANE_H = 22, NT_H = 20, GUTTER = 100, PAD = 4, COV_H = 60, CLIP_H = 36, IDENT_H = 36, ORIGIN_H = 18, GAP = 12;
  var DIFF_H = 14, ROW_GAP = 10, DIFF_RED = '#c0392b', GAP_INK = '#888888';
  var BASE = { A: '#4faf45', C: '#e0a53f', G: '#e0555a', T: '#4a90d9', N: '#666666' };
  var viewers = {}, last = {};

  function cssVar(name, fallback) {
    var v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return v || fallback;
  }
  function niceStep(raw) {
    var p = Math.pow(10, Math.floor(Math.log10(raw)));
    var m = raw / p;
    return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10) * p;
  }

  function Viewer(id) {
    this.id = id;
    this.canvas = document.getElementById(id);
    this.wrap = this.canvas.parentElement;
    this.tip = document.getElementById(id.replace(/-canvas$/, '-tip'));
    this.section = this.canvas.closest('details');
    this.len = 0; this.seq = ''; this.feats = []; this.selected = null; this.nLanes = 0;
    this.viewStart = 1; this.ppb = 1;
    this.showNt = true; this.showCov = true; this.showClip = true;
    this.depth = null; this.clip = null; this.depthMax = 0; this.clipFlag = 0.4;
    this.hits = [];
    this.bindControls();
    var self = this;
    if (window.ResizeObserver) new ResizeObserver(function () { self.draw(); }).observe(this.wrap);
    if (this.section) this.section.addEventListener('toggle', function () { self.draw(); });
    this.wrap.addEventListener('scroll', function () { if (self.rows) self.draw(); });
  }

  Viewer.prototype.load = function (p) {
    var seq = p.rows ? (p.rows[0] ? p.rows[0].seq : '') : (p.seq || '');
    var sameSeq = this.seq === seq && !!this.rows === !!p.rows;
    this.len = p.len; this.seq = seq; this.input = p.input; this.unit = p.unit || 'unit';
    this.rows = p.rows || null;
    if (this.rows) {
      this.rows.forEach(function (r) {
        r.feats = (r.features || []).map(function (f) { return Object.assign({}, f); });
        r.nLanes = Math.max(1, G.lanes(r.feats));
        r.first = r.seq.search(/[^-]/) + 1; r.last = r.seq.replace(/-+$/, '').length;
      });
      this.feats = [];
      this.diff = G.diffColumns(this.rows.map(function (r) { return r.seq; }));
      this.diffSet = {}; this.diff.forEach(function (c) { this.diffSet[c] = true; }, this);
      this.nLanes = 0; this.depth = null; this.clip = null; this.origin = null; this.junctions = null; this.ident = null;
    } else {
      this.feats = (p.features || []).map(function (f) { return Object.assign({}, f); });
      this.nLanes = G.lanes(this.feats);
      this.depth = p.depth || null; this.clip = p.clip || null;
      this.origin = p.origin || null; this.junctions = p.junctions || null; this.ident = p.ident || null;
      this.diff = null; this.diffSet = null;
    }
    if (p.cutoff) this.clipFlag = p.cutoff;
    this.depthMax = (this.depth || []).reduce(function (m, d) { return d !== null && d > m ? d : m; }, 1);
    this.selected = (p.selected === undefined || p.selected === null) ? null : p.selected;
    if (!sameSeq) this.whole();
    this.draw();
  };

  // ---- view state ----
  Viewer.prototype.width = function () { return Math.max(200, this.wrap.clientWidth - GUTTER); };
  Viewer.prototype.minPpb = function () { return this.width() / this.len; };
  Viewer.prototype.viewLen = function () { return this.width() / this.ppb; };
  // a whole view stays whole when the box is resized
  Viewer.prototype.clamp = function () {
    if (this.isWhole) { this.ppb = this.minPpb(); this.viewStart = 1; }
    this.ppb = Math.min(MAX_PPB, Math.max(this.minPpb(), this.ppb));
    this.viewStart = Math.max(1, Math.min(this.len - this.viewLen() + 1, this.viewStart));
  };
  Viewer.prototype.whole = function () { this.isWhole = true; this.clamp(); this.draw(); };
  Viewer.prototype.zoom = function (factor, atPos) {
    var centre = atPos || (this.viewStart + this.viewLen() / 2);
    var newPpb = Math.min(MAX_PPB, Math.max(this.minPpb(), this.ppb * factor));
    this.isWhole = false;
    this.viewStart = centre - (centre - this.viewStart) * (this.ppb / newPpb);
    this.ppb = newPpb; this.clamp(); this.draw();
  };
  Viewer.prototype.goto = function (pos) {
    this.isWhole = false;
    if (this.ppb < NT_LETTER) this.ppb = 10;
    this.viewStart = pos - this.viewLen() / 2; this.clamp(); this.draw();
  };
  Viewer.prototype.fit = function (row) {
    var f = this.feats.find(function (x) { return x.row === row; });
    if (!f) return;
    this.selected = row; this.fitRange(f.pos1, f.pos2);
  };
  Viewer.prototype.fitRange = function (pos1, pos2) {
    // open first: a collapsed section has no width to measure
    if (this.section) this.section.open = true;
    this.isWhole = false;
    var s = pos2 - pos1 + 1, margin = Math.max(3, s * 0.05);
    this.ppb = Math.min(MAX_PPB, this.width() / (s + 2 * margin));
    this.viewStart = (pos1 + pos2) / 2 - this.viewLen() / 2; this.clamp();
    this.draw();
  };
  // nearest differing column outside the view on the given side; wraps
  Viewer.prototype.stepDiff = function (dir) {
    if (!this.diff || !this.diff.length) return;
    var vs = Math.ceil(this.viewStart), ve = Math.floor(this.viewStart + this.viewLen()), col = null;
    if (dir > 0) { for (var i = 0; i < this.diff.length; i++) if (this.diff[i] > ve) { col = this.diff[i]; break; } if (col === null) col = this.diff[0]; }
    else { for (var j = this.diff.length - 1; j >= 0; j--) if (this.diff[j] < vs) { col = this.diff[j]; break; } if (col === null) col = this.diff[this.diff.length - 1]; }
    this.goto(col);
  };

  // ---- coordinate helpers ----
  Viewer.prototype.x = function (pos) { return GUTTER + (pos - this.viewStart) * this.ppb; };
  Viewer.prototype.inView = function (f) { return f.pos2 >= this.viewStart && f.pos1 <= this.viewStart + this.viewLen(); };
  Viewer.prototype.baseAt = function (pos) {
    pos = Math.floor(pos);
    if (pos < 1 || pos > this.len) return null;
    return { pos: pos, base: this.seq.charAt(pos - 1) || 'N' };
  };

  // ---- drawing ----
  Viewer.prototype.covOn = function () { return this.showCov && !!this.depth; };
  Viewer.prototype.clipOn = function () { return this.showClip && !!this.clip; };
  Viewer.prototype.identOn = function () { return this.showIdent && !!this.ident; };
  Viewer.prototype.originOn = function () { return this.showOrigin && !!this.origin; };
  Viewer.prototype.ntOn = function () { return this.showNt && this.ppb >= NT_BAR; };
  Viewer.prototype.topH = function () {
    return RULER_H + PAD + (this.covOn() ? COV_H + GAP : 0) + (this.clipOn() ? CLIP_H + GAP : 0);
  };
  Viewer.prototype.laneY = function (lane) { return this.topH() + lane * LANE_H; };
  // chimera tracks sit under the genes and bases
  Viewer.prototype.chimY = function () {
    return this.topH() + this.nLanes * LANE_H + PAD + (this.ntOn() ? NT_H + PAD : 0) + GAP;
  };
  Viewer.prototype.height = function () {
    if (this.rows) return this.alignHeight();
    return this.chimY() + (this.identOn() ? IDENT_H + GAP : 0) + (this.originOn() ? ORIGIN_H + GAP : 0);
  };
  // align mode: y of each row's first lane; the nt strip sits under its lanes
  Viewer.prototype.rowTop = function (i) {
    var y = RULER_H + PAD + DIFF_H + GAP;
    for (var k = 0; k < i; k++) y += this.rows[k].nLanes * LANE_H + (this.ntOn() ? NT_H : 0) + ROW_GAP;
    return y;
  };
  Viewer.prototype.rowNtY = function (i) { return this.rowTop(i) + this.rows[i].nLanes * LANE_H; };
  Viewer.prototype.alignHeight = function () {
    return this.rowTop(this.rows.length) - ROW_GAP + PAD;
  };
  Viewer.prototype.draw = function () {
    if (!this.len) { this.hits = []; this.canvas.width = 0; this.canvas.style.height = '0px'; return; }
    if (!this.canvas.offsetParent) return;
    this.clamp();
    var dpr = window.devicePixelRatio || 1, W = this.wrap.clientWidth, H = this.height();
    this.canvas.width = W * dpr; this.canvas.height = H * dpr;
    this.canvas.style.height = H + 'px';
    var c = this.canvas.getContext('2d'); c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.clearRect(0, 0, W, H);
    c.font = '12px ' + cssVar('--mp-font-mono', 'monospace');
    this.hits = [];
    if (this.rows) { this.drawAlign(c, W); return; }
    this.drawRuler(c);
    var y = RULER_H + PAD;
    if (this.covOn()) { this.drawTrack(c, W, y, COV_H, this.depth, this.depthMax, 'depth', false); y += COV_H + GAP; }
    if (this.clipOn()) { this.drawTrack(c, W, y, CLIP_H, this.clip, 1, 'clipped', true); y += CLIP_H + GAP; }
    this.drawLanes(c);
    if (this.ntOn()) this.drawNt(c, this.topH() + this.nLanes * LANE_H + PAD);
    y = this.chimY();
    var self = this, segColor = function (pos) { var o = self.originAt(pos); return o ? o.color : '#9e9e9e'; };
    if (this.identOn()) { this.drawTrack(c, W, y, IDENT_H, this.ident, 100, 'BLAST %id', false, { vmin: 70, colorAt: segColor, stroke: cssVar('--rp-ident-stroke', '#4d4d4d') }); y += IDENT_H + GAP; }
    if (this.originOn()) { this.drawOrigin(c, W, y); y += ORIGIN_H + GAP; }
    this.drawJunctions(c, H);
  };
  Viewer.prototype.drawRuler = function (c) {
    var step = niceStep(90 / this.ppb), vs = this.viewStart, ve = vs + this.viewLen();
    c.fillStyle = cssVar('--mp-text-muted', '#6a6a6a'); c.strokeStyle = cssVar('--mp-border', '#ccc');
    c.textAlign = 'center'; c.textBaseline = 'top';
    for (var pos = Math.ceil(vs / step) * step; pos <= ve; pos += step) {
      var x = this.x(pos);
      if (!this.baseAt(pos) || x < GUTTER) continue;
      c.beginPath(); c.moveTo(x, RULER_H - 6); c.lineTo(x, RULER_H); c.stroke();
      c.fillText(String(pos), x, 2);
    }
  };
  Viewer.prototype.drawLanes = function (c) { this.drawLanesAt(c, this.feats, this.topH()); };
  Viewer.prototype.drawLanesAt = function (c, feats, top) {
    var vs = this.viewStart, ve = this.viewStart + this.viewLen();
    c.textAlign = 'center'; c.textBaseline = 'middle';
    feats.forEach(function (f) {
      if (!this.inView(f)) return;
      var col = f.color || '#888888', y = top + f.lane * LANE_H, h = LANE_H - 6, edge = G.partialEdge(f);
      var x0 = Math.max(GUTTER, this.x(f.pos1)), x1 = Math.min(this.wrap.clientWidth, this.x(f.pos2 + 1));
      if (x1 - x0 < 1) return;
      var fwd = f.dir !== '-', head = Math.min(8, x1 - x0);
      // vertices in drawing order; pts[4]-pts[0] is the closing edge.
      var pts = fwd ? [[x0, y], [x1 - head, y], [x1, y + h / 2], [x1 - head, y + h], [x0, y + h]]
                    : [[x1, y], [x0 + head, y], [x0, y + h / 2], [x0 + head, y + h], [x1, y + h]];
      c.beginPath();
      c.moveTo(pts[0][0], pts[0][1]);
      for (var i = 1; i < pts.length; i++) c.lineTo(pts[i][0], pts[i][1]);
      c.closePath();
      c.fillStyle = col + '55'; c.fill();
      c.lineWidth = f.row === this.selected ? 2 : 1;
      c.strokeStyle = f.row === this.selected ? cssVar('--mp-primary', '#337ab7') : col;

      // the flat closing edge is the pos1 side when fwd, pos2 side otherwise;
      // partial ends are dashed
      var owns = edge ? G.segOwns(f, vs, ve) : null;
      var dashStart = !!edge && (edge === 'start' || edge === 'both') && owns.start;
      var dashEnd = !!edge && (edge === 'end' || edge === 'both') && owns.end;
      var dashFlat = fwd ? dashStart : dashEnd;
      var dashTip = fwd ? dashEnd : dashStart;

      c.setLineDash([]);
      c.beginPath(); c.moveTo(pts[0][0], pts[0][1]); c.lineTo(pts[1][0], pts[1][1]); c.stroke();
      c.beginPath(); c.moveTo(pts[3][0], pts[3][1]); c.lineTo(pts[4][0], pts[4][1]); c.stroke();
      c.setLineDash(dashTip ? [3, 2] : []);
      c.beginPath(); c.moveTo(pts[1][0], pts[1][1]); c.lineTo(pts[2][0], pts[2][1]); c.lineTo(pts[3][0], pts[3][1]); c.stroke();
      c.setLineDash(dashFlat ? [3, 2] : []);
      c.beginPath(); c.moveTo(pts[4][0], pts[4][1]); c.lineTo(pts[0][0], pts[0][1]); c.stroke();
      c.setLineDash([]);

      if (x1 - x0 > c.measureText(f.gene).width + 8) { c.fillStyle = cssVar('--mp-text', '#333'); c.fillText(f.gene, (x0 + x1) / 2, y + h / 2); }
      this.hits.push({ x0: x0, x1: x1, y0: y, y1: y + h, f: f, row: f.row });
    }, this);
  };
  Viewer.prototype.drawAlign = function (c, W) {
    var y = RULER_H + PAD;
    var self = this, H = this.alignHeight();
    // red wash over every differing column in view, behind all rows
    if (this.ppb >= NT_LETTER) {
      c.fillStyle = 'rgba(192, 57, 43, 0.12)';
      var vs = Math.floor(this.viewStart), ve = Math.ceil(this.viewStart + this.viewLen());
      for (var col = vs; col <= ve; col++) {
        if (!this.diffSet[col]) continue;
        var x = this.x(col); if (x < GUTTER) continue;
        c.fillRect(x, y + DIFF_H + GAP, Math.max(1, this.ppb), H - (y + DIFF_H + GAP));
      }
    }
    this.rows.forEach(function (r, i) {
      var top = self.rowTop(i);
      self.drawLanesAt(c, r.feats, top);
      if (self.ntOn()) self.drawRowNt(c, r, self.rowNtY(i));
      // "p1.s1 keep": the verdict takes the K/D badge colour of the contig table
      var ly = self.ntOn() ? self.rowNtY(i) + NT_H / 2 : top + LANE_H / 2, parts = r.label.split(' ');
      var verdict = parts.length > 1 ? parts[1] : '', vcol = verdict === 'keep' ? '#2d6a4f' : verdict === 'drop' ? '#c0392b' : null;
      c.textAlign = 'right'; c.textBaseline = 'middle';
      var vx = GUTTER - 6, vw = verdict ? c.measureText(verdict).width : 0;
      if (verdict) { c.fillStyle = vcol || cssVar('--mp-text', '#333'); c.fillText(verdict, vx, ly); }
      c.fillStyle = cssVar('--mp-text', '#333');
      c.fillText(parts[0], vx - vw - (verdict ? 5 : 0), ly, GUTTER - 10 - vw);
    });
    // ruler and Differences track stay pinned while the box scrolls
    var hdr = Math.max(0, Math.min(this.wrap.scrollTop, H - (y + DIFF_H + GAP)));
    c.save(); c.translate(0, hdr);
    c.fillStyle = cssVar('--mp-surface', '#ffffff'); c.fillRect(0, 0, W, y + DIFF_H + GAP / 2);
    this.drawRuler(c);
    this.drawDiffTrack(c, W, y);
    c.restore();
    this.diffTrack.y0 += hdr; this.diffTrack.y1 += hdr;
  };
  Viewer.prototype.drawDiffTrack = function (c, W, y) {
    c.fillStyle = cssVar('--mp-surface-alt', '#f5f5f5'); c.fillRect(GUTTER, y, W - GUTTER, DIFF_H);
    c.fillStyle = DIFF_RED;
    var w = Math.max(1, this.ppb), vs = this.viewStart, ve = vs + this.viewLen();
    for (var i = 0; i < this.diff.length; i++) {
      var col = this.diff[i]; if (col < vs - 1 || col > ve + 1) continue;
      var x = this.x(col); if (x < GUTTER || x > W) continue;
      c.fillRect(x, y, w, DIFF_H);
    }
    c.strokeStyle = cssVar('--mp-border', '#ccc'); c.beginPath(); c.moveTo(GUTTER, y + DIFF_H + 0.5); c.lineTo(W, y + DIFF_H + 0.5); c.stroke();
    c.fillStyle = cssVar('--mp-text-muted', '#6a6a6a'); c.textAlign = 'right'; c.textBaseline = 'middle';
    c.fillText('Differences', GUTTER - 6, y + DIFF_H / 2, GUTTER - 10);
    this.diffTrack = { y0: y, y1: y + DIFF_H };
  };
  Viewer.prototype.drawRowNt = function (c, r, y) {
    var vs = Math.floor(this.viewStart), ve = Math.ceil(this.viewStart + this.viewLen());
    var letter = this.ppb >= NT_LETTER, w = Math.max(1, this.ppb - (letter ? 1 : 0));
    c.textAlign = 'center'; c.textBaseline = 'middle';
    for (var col = vs; col <= ve; col++) {
      if (col < 1 || col > this.len) continue;
      var x = this.x(col); if (x < GUTTER) continue;
      var b = r.seq.charAt(col - 1) || 'N';
      if (b === '-') {
        if (col >= r.first && col <= r.last) { c.fillStyle = GAP_INK; c.fillRect(x + w * 0.2, y + NT_H / 2 - 1, w * 0.6, 2); }
        continue;
      }
      c.fillStyle = (BASE[b] || BASE.N) + (letter ? '99' : 'cc');
      c.fillRect(x, y + 2, w, NT_H - 4);
      if (letter) { c.fillStyle = '#000000'; c.fillText(b, x + this.ppb / 2, y + NT_H / 2); }
    }
  };
  // One bar per pixel column, the max of the positions under it. Depth is
  // blue; clipped-read fraction is grey with bars at or over the chimera
  // check's cutoff in Okabe-Ito orange. Runs of bars are filled as one shape,
  // since separate 1 px bars leave seams on scaled displays.
  Viewer.prototype.drawTrack = function (c, W, y, h, vals, vmax, label, flag, opt) {
    opt = opt || {}; var vmin = opt.vmin || 0;
    var base = flag ? '#a0a0a0' : cssVar('--mp-primary', '#337ab7'), red = '#e69f00', muted = cssVar('--mp-text-muted', '#6a6a6a');
    var cut = this.clipFlag, cols = [];
    c.fillStyle = cssVar('--mp-surface-alt', '#f5f5f5'); c.fillRect(GUTTER, y, W - GUTTER, h);
    for (var px = GUTTER; px < W; px++) {
      var a = Math.floor(this.viewStart + (px - GUTTER) / this.ppb), b = Math.floor(this.viewStart + (px + 1 - GUTTER) / this.ppb);
      var v = null, gap = false;
      for (var pos = a; pos <= b; pos++) {
        if (!this.baseAt(pos)) continue;
        var d = vals[pos - 1]; if (d === null || d === undefined) continue;
        if (v === null || d > v) v = d;
        if (d === 0) gap = true;
      }
      cols.push(v === null ? null : { v: v, gap: gap, c: opt.colorAt ? opt.colorAt(a) : base });
    }
    // columns holding any uncovered base get a full-height red wash, so gaps
    // survive the per-column max when zoomed out
    var area = function (fill, keep, full, stroke) {
      c.fillStyle = fill; c.beginPath();
      var open = false;
      for (var i = 0; i <= cols.length; i++) {
        var k = i < cols.length && cols[i] && keep(cols[i]);
        var top = k ? (full ? y : y + h - Math.min(h, Math.max(0, cols[i].v - vmin) / (vmax - vmin) * h)) : 0;
        if (k && !open) { c.moveTo(GUTTER + i, y + h); open = true; }
        if (k) { c.lineTo(GUTTER + i, top); c.lineTo(GUTTER + i + 1, top); }
        if (!k && open) { c.lineTo(GUTTER + i, y + h); c.closePath(); open = false; }
      }
      c.fill();
      if (stroke) { c.strokeStyle = stroke; c.lineWidth = 1; c.stroke(); }
    };
    if (!flag) area('rgba(217, 83, 79, 0.35)', function (o) { return o.gap; }, true);
    var fills = cols.reduce(function (acc, o) { if (o && acc.indexOf(o.c) < 0) acc.push(o.c); return acc; }, []);
    fills.forEach(function (f) { area(f, function (o) { return o.v > vmin && o.c === f; }, false, opt.stroke); });
    if (flag) area(red, function (o) { return o.v >= cut; }, false);
    c.strokeStyle = cssVar('--mp-border', '#ccc'); c.beginPath(); c.moveTo(GUTTER, y + h + 0.5); c.lineTo(W, y + h + 0.5); c.stroke();
    if (flag) { c.save(); c.strokeStyle = red; c.setLineDash([2, 3]); var fy = y + h - cut / vmax * h;
      c.beginPath(); c.moveTo(GUTTER, fy + 0.5); c.lineTo(W, fy + 0.5); c.stroke(); c.restore(); }
    c.fillStyle = muted; c.textAlign = 'right';
    c.textBaseline = 'top'; c.fillText(flag ? (vmax * 100).toFixed(0) + '%' : String(vmax), GUTTER - 6, y);
    c.textBaseline = 'bottom'; c.fillText(label, GUTTER - 6, y + h, GUTTER - 10);
  };
  // One rect per origin segment in its colour, labelled when the label fits
  // the visible part (same rule as the gene lanes).
  Viewer.prototype.drawOrigin = function (c, W, y) {
    c.textAlign = 'center'; c.textBaseline = 'middle';
    this.origin.forEach(function (s) {
      if (!this.inView(s)) return;
      var x0 = Math.max(GUTTER, this.x(s.pos1)), x1 = Math.min(W, this.x(s.pos2 + 1));
      if (x1 - x0 < 1) return;
      c.fillStyle = s.color; c.fillRect(x0, y, x1 - x0, ORIGIN_H);
      c.strokeStyle = cssVar('--mp-border', '#ccc'); c.strokeRect(x0 + 0.5, y + 0.5, x1 - x0 - 1, ORIGIN_H - 1);
      var full = originText(s), short = s.taxon === 'nohit' ? 'no hit' : s.taxon;
      var lab = x1 - x0 > c.measureText(full).width + 8 ? full : x1 - x0 > c.measureText(short).width + 8 ? short : null;
      if (lab) { c.fillStyle = '#333'; c.fillText(lab, (x0 + x1) / 2, y + ORIGIN_H / 2); }
    }, this);
    c.fillStyle = cssVar('--mp-text-muted', '#6a6a6a'); c.textAlign = 'right';
    c.fillText('BLAST origin', GUTTER - 6, y + ORIGIN_H / 2, GUTTER - 10);
  };
  Viewer.prototype.drawJunctions = function (c, H) {
    if (!this.junctions) return;
    c.save(); c.setLineDash([5, 3]); c.lineWidth = 1.5;
    this.junctions.forEach(function (j) {
      var x = this.x(j.pos); if (x < GUTTER || x > this.wrap.clientWidth) return;
      c.strokeStyle = j.severe ? '#c0392b' : '#b9770e';
      c.beginPath(); c.moveTo(x, RULER_H); c.lineTo(x, H); c.stroke();
    }, this);
    c.restore();
  };
  Viewer.prototype.originAt = function (pos) {
    return (this.origin || []).find(function (s) { return pos >= s.pos1 && pos <= s.pos2; }) || null;
  };
  function originText(s) {
    if (s.host) return s.taxon + ' (your lineage)';
    if (s.taxon === 'nohit') return 'no hit';
    return s.taxon + (s.species ? ', ' + s.species : '') + (s.pident ? ' ' + s.pident.toFixed(1) + '%' : '');
  }
  Viewer.prototype.drawNt = function (c, y) {
    var vs = Math.floor(this.viewStart), ve = Math.ceil(this.viewStart + this.viewLen());
    var letter = this.ppb >= NT_LETTER, w = Math.max(1, this.ppb - (letter ? 1 : 0));
    c.textAlign = 'center'; c.textBaseline = 'middle';
    for (var pos = vs; pos <= ve; pos++) {
      var b = this.baseAt(pos), x = this.x(pos);
      if (!b || x < GUTTER) continue;
      c.fillStyle = (BASE[b.base] || BASE.N) + (letter ? '99' : 'cc');
      c.fillRect(x, y + 2, w, NT_H - 4);
      if (letter) { c.fillStyle = '#000000'; c.fillText(b.base, x + this.ppb / 2, y + NT_H / 2); }
    }
    c.fillStyle = cssVar('--mp-text-muted', '#6a6a6a'); c.textAlign = 'right';
    c.fillText('nt', GUTTER - 6, y + NT_H / 2, GUTTER - 10);
  };

  // The bases in view as FASTA, header ">unit:start-end", 60 per line.
  Viewer.prototype.visibleFasta = function () {
    var a = Math.max(1, Math.ceil(this.viewStart)), b = Math.min(this.len, Math.floor(this.viewStart + this.viewLen()));
    var fa = function (name, s) { return '>' + name + ':' + a + '-' + b + '\n' + (s.match(/.{1,60}/g) || []).join('\n') + '\n'; };
    if (this.rows) {
      var n = 0, text = this.rows.map(function (r) { var s = r.seq.substring(a - 1, b); n += s.length; return fa(r.label.replace(/\s+/g, '_'), s); }).join('');
      return { n: n, text: text };
    }
    var s = this.seq.substring(a - 1, b);
    return { n: s.length, text: fa(this.unit, s) };
  };
  // Plain http off localhost has no Clipboard API, so fall back to a textarea.
  function copyText(t) {
    if (navigator.clipboard && window.isSecureContext) { navigator.clipboard.writeText(t); return; }
    var a = document.createElement('textarea'); a.value = t; a.style.position = 'fixed'; a.style.opacity = '0';
    document.body.appendChild(a); a.select();
    try { document.execCommand('copy'); } finally { document.body.removeChild(a); }
  }

  // ---- interaction ----
  Viewer.prototype.hitTest = function (px, py) {
    for (var i = this.hits.length - 1; i >= 0; i--) {
      var h = this.hits[i];
      if (px >= h.x0 && px <= h.x1 && py >= h.y0 && py <= h.y1) return h;
    }
    return null;
  };
  Viewer.prototype.click = function (px, py) {
    if (this.rows && this.diffTrack && py >= this.diffTrack.y0 && py <= this.diffTrack.y1 && px >= GUTTER) {
      this.goto(Math.floor(this.viewStart + (px - GUTTER) / this.ppb)); return;
    }
    var h = this.hitTest(px, py);
    if (!h) return;
    this.selected = h.f.row;
    if (this.input && window.Shiny) {
      window.Shiny.setInputValue(this.input, { row: h.f.row, nonce: Date.now() }, { priority: 'event' });
    }
    this.draw();
  };
  Viewer.prototype.bindControls = function () {
    var self = this, sec = this.section || document, prefix = this.id.replace(/-canvas$/, '');
    var cv = this.canvas, dragging = null;
    cv.addEventListener('wheel', function (e) {
      e.preventDefault();
      var rect = cv.getBoundingClientRect(), px = e.clientX - rect.left;
      var pan = e.shiftKey ? e.deltaY : Math.abs(e.deltaX) > Math.abs(e.deltaY) ? e.deltaX : 0;
      if (pan) { self.isWhole = false; self.viewStart += pan / self.ppb; self.clamp(); self.draw(); return; }
      self.zoom(Math.pow(1.15, -e.deltaY / 100), self.viewStart + (px - GUTTER) / self.ppb);
    }, { passive: false });
    cv.addEventListener('mousedown', function (e) { dragging = { x: e.clientX, start: self.viewStart, moved: false }; });
    window.addEventListener('mousemove', function (e) {
      if (!dragging) return;
      var dx = e.clientX - dragging.x; if (Math.abs(dx) > 2) { dragging.moved = true; self.isWhole = false; }
      self.viewStart = dragging.start - dx / self.ppb; self.clamp(); self.draw();
    });
    window.addEventListener('mouseup', function (e) {
      if (!dragging) return;
      var rect = cv.getBoundingClientRect();
      if (!dragging.moved) self.click(e.clientX - rect.left, e.clientY - rect.top);
      dragging = null;
    });
    cv.addEventListener('mousemove', function (e) {
      var rect = cv.getBoundingClientRect(), px = e.clientX - rect.left, py = e.clientY - rect.top;
      var h = self.hitTest(px, py), b = self.baseAt(self.viewStart + (px - GUTTER) / self.ppb);
      if (!self.tip) return;
      if (!b || px < GUTTER) { self.tip.hidden = true; return; }
      if (self.rows) {
        var t2 = 'Column ' + b.pos.toLocaleString();
        self.rows.forEach(function (r) {
          var np = G.nativePos(r.seq, b.pos), ch = r.seq.charAt(b.pos - 1);
          t2 += ' | ' + r.label.split(' ')[0] + ' ' + (np === null ? 'gap' : np.toLocaleString() + ' ' + ch);
        });
        if (h) t2 += ' | ' + h.f.gene + ' ' + h.f.pos1.toLocaleString() + '-' + h.f.pos2.toLocaleString() + ' (' + h.f.dir + ')';
        self.tip.textContent = t2; self.tip.hidden = false;
        var tw2 = self.tip.offsetWidth, th2 = self.tip.offsetHeight;
        self.tip.style.left = Math.max(0, Math.min(px + 12, cv.clientWidth - tw2)) + 'px';
        self.tip.style.top = (py + 12 + th2 > cv.clientHeight ? py - th2 - 8 : py + 12) + 'px';
        return;
      }
      var t = 'Position ' + b.pos.toLocaleString() + ', ' + b.base;
      if (self.covOn() && self.depth[b.pos - 1] !== null) t += ' | depth ' + self.depth[b.pos - 1];
      if (self.clipOn() && self.clip[b.pos - 1] !== null) t += ', clipped ' + (self.clip[b.pos - 1] * 100).toFixed(0) + '%';
      if (h) t += ' | ' + h.f.gene + ' ' + h.f.pos1.toLocaleString() + '-' + h.f.pos2.toLocaleString() + ' (' + h.f.dir + ')';
      if (self.identOn() && self.ident[b.pos - 1] !== null) t += ' | BLAST id ' + self.ident[b.pos - 1].toFixed(1) + '%';
      var o = self.originOn() && self.originAt(b.pos);
      if (o) t += ' | ' + originText(o);
      var j = (self.junctions || []).find(function (q) { return Math.abs(self.x(q.pos) - px) <= 6; });
      if (j) t += ' | junction ' + j.pos.toLocaleString() + (j.checks ? ' (' + j.checks + ')' : '');
      self.tip.textContent = t; self.tip.hidden = false;
      var tw = self.tip.offsetWidth, th = self.tip.offsetHeight;
      self.tip.style.left = Math.max(0, Math.min(px + 12, cv.clientWidth - tw)) + 'px';
      self.tip.style.top = (py + 12 + th > cv.clientHeight ? py - th - 8 : py + 12) + 'px';
    });
    cv.addEventListener('mouseleave', function () { if (self.tip) self.tip.hidden = true; });
    sec.querySelectorAll('[data-rpseq]').forEach(function (btn) {
      btn.addEventListener('click', function () {
        var a = btn.getAttribute('data-rpseq');
        if (a === 'zoom_in') self.zoom(2); else if (a === 'zoom_out') self.zoom(0.5);
        else if (a === 'whole') self.whole(); else if (a === 'fit' && self.selected !== null) self.fit(self.selected);
        else if (a === 'prev_diff') self.stepDiff(-1); else if (a === 'next_diff') self.stepDiff(1);
        else if (a === 'copy') {
          var f = self.visibleFasta(); copyText(f.text);
          var old = btn.textContent; btn.textContent = 'Copied ' + f.n.toLocaleString() + ' bp';
          setTimeout(function () { btn.textContent = old; }, 1500);
        }
      });
    });
    [['show_nt', 'showNt'], ['show_cov', 'showCov'], ['show_clip', 'showClip'], ['show_ident', 'showIdent'],
     ['show_origin', 'showOrigin']].forEach(function (pair) {
      var box = document.getElementById(prefix + '-' + pair[0]);
      if (box) { self[pair[1]] = box.checked; box.addEventListener('change', function () { self[pair[1]] = box.checked; self.draw(); }); }
    });
    var go = document.getElementById(prefix + '-goto');
    if (go) go.addEventListener('keydown', function (e) { if (e.key === 'Enter') { var v = parseInt(go.value, 10); if (v >= 1 && v <= self.len) self.goto(v); } });
  };

  function get(id) {
    var el = document.getElementById(id);
    if (!el) return null;
    if (viewers[id] && viewers[id].canvas !== el) delete viewers[id];
    return viewers[id] || (viewers[id] = new Viewer(id));
  }
  // The modal can reach the DOM after the message that fills it, so wait.
  function withViewer(id, fn, tries) {
    var v = get(id);
    if (v) { fn(v); return; }
    if ((tries || 0) >= 40) return;
    setTimeout(function () { withViewer(id, fn, (tries || 0) + 1); }, 250);
  }
  if (window.Shiny) {
    window.Shiny.addCustomMessageHandler('rpseq', function (p) {
      last[p.id] = p;
      withViewer(p.id, function (v) { v.load(p); });
    });
    window.Shiny.addCustomMessageHandler('rpseq_select', function (p) {
      withViewer(p.id, function (v) { v.fit(p.row); });
    });
    window.Shiny.addCustomMessageHandler('rpseq_goto', function (p) {
      withViewer(p.id, function (v) { v.fitRange(p.pos1, p.pos2); });
    });
  }
  // A reopened window has a fresh canvas the server may not message again.
  $(document).on('shown.bs.modal', '#shiny-modal', function () {
    Object.keys(last).forEach(function (id) {
      var el = document.getElementById(id);
      if (el && (!viewers[id] || viewers[id].canvas !== el)) get(id).load(last[id]);
    });
  });
  window.rpseq.state = function (id) {
    var v = viewers[id]; if (!v) return null;
    return { len: v.len, viewStart: v.viewStart, ppb: v.ppb, selected: v.selected, height: v.height(), nLanes: v.nLanes,
             origin: (v.origin || []).length, junctions: (v.junctions || []).length,
             rows: (v.rows || []).length, diff: (v.diff || []).length,
             features: v.feats.map(function (f) { return { row: f.row, gene: f.gene, pos1: f.pos1, pos2: f.pos2, lane: f.lane }; }) };
  };
  window.rpseq.diffAt = function (id, col) { var v = viewers[id]; return !!(v && v.diffSet && v.diffSet[col]); };
  window.rpseq.fasta = function (id) { var v = get(id); return v ? v.visibleFasta() : null; };
  window.rpseq.hitTest = function (id, px, py) { var v = get(id); return v ? v.hitTest(px, py) : null; };
  window.rpseq.click = function (id, px, py) { var v = get(id); if (v) v.click(px, py); };
  window.rpseq.laneY = function (id, lane) { var v = get(id); return v ? v.laneY(lane) + 8 : 0; };
})();
