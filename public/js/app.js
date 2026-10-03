// 浏览器端应用：录入草稿、调用联合求解器、渲染谱系图与逐段位移。
// 任一帧或限制被修改即作废当前谱系，不保留旧结果（草稿本身保留在编辑区）。
'use strict';

import { normalizeSpec, solveLineage, presentSolution } from './lineage.js';

const STORAGE_KEY = 'algal-lineage-draft-v1';

const SAMPLE = {
  frames: [
    [
      { id: 'A', x: 10, y: 50, b: 50 },
      { id: 'X0', x: 95, y: 95, b: 60 },
    ],
    [
      { id: 'B', x: 20, y: 48, b: 55 },
      { id: 'X1', x: 90, y: 90, b: 70 },
    ],
    [
      { id: 'X2a', x: 88, y: 88, b: 80 },
      { id: 'X2b', x: 60, y: 20, b: 30 },
    ],
    [
      { id: 'D', x: 36, y: 46, b: 58 },
      { id: 'X3', x: 85, y: 85, b: 90 },
    ],
    [
      { id: 'E1', x: 44, y: 44, b: 52 },
      { id: 'E2', x: 45, y: 53, b: 50 },
      { id: 'X4', x: 80, y: 80, b: 100 },
    ],
  ],
  startId: 'A',
  maxDist: 12,
  maxSkip: 1,
  target: 2,
};

const $ = (sel) => document.querySelector(sel);
const els = {
  frames: $('#frames'),
  addFrame: $('#btn-add-frame'),
  sample: $('#btn-sample'),
  startSelect: $('#start-select'),
  maxDist: $('#max-dist'),
  maxSkip: $('#max-skip'),
  target: $('#target'),
  solve: $('#btn-solve'),
  stale: $('#stale-hint'),
  errors: $('#form-errors'),
  empty: $('#result-empty'),
  infeasible: $('#result-infeasible'),
  ok: $('#result-ok'),
  stats: $('#stats'),
  chart: $('#chart'),
  tbody: $('#edges-table tbody'),
  used: $('#used-list'),
};

let draft = loadDraft();
let lastSolution = null; // 仅保存最近一次「复原」的结果；编辑即置空

function loadDraft() {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (raw) return JSON.parse(raw);
  } catch { /* 忽略损坏草稿 */ }
  return structuredClone(SAMPLE);
}

function saveDraft() {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(draft));
  } catch { /* 存储不可用时静默 */ }
}

// ---- 编辑区渲染 ----
function renderFrames() {
  els.frames.innerHTML = '';
  draft.frames.forEach((frame, t) => {
    const card = document.createElement('div');
    card.className = 'frame-card';
    card.dataset.frame = String(t);

    const head = document.createElement('div');
    head.className = 'frame-title';
    head.innerHTML = `<span><b>第 ${t + 1} 帧</b>（t=${t}）</span>`;
    if (draft.frames.length > 4) {
      const del = document.createElement('button');
      del.type = 'button';
      del.className = 'icon-btn';
      del.textContent = '删除此帧 ✕';
      del.addEventListener('click', () => {
        draft.frames.splice(t, 1);
        invalidate();
        renderAll();
      });
      head.appendChild(del);
    }
    card.appendChild(head);

    frame.forEach((spot, i) => {
      const row = document.createElement('div');
      row.className = 'spot-row';
      row.innerHTML = `
        <input class="spot-id-input" data-k="id" placeholder="编号" value="${escapeAttr(spot.id)}" />
        <input data-k="x" type="number" step="1" placeholder="x" value="${spot.x ?? ''}" />
        <input data-k="y" type="number" step="1" placeholder="y" value="${spot.y ?? ''}" />
        <input data-k="b" type="number" step="1" min="0" placeholder="亮度" value="${spot.b ?? ''}" />
        <button type="button" class="icon-btn" title="删除斑点">✕</button>`;
      row.querySelectorAll('input').forEach((input) => {
        input.addEventListener('input', () => {
          const k = input.dataset.k;
          spot[k] = k === 'id' ? input.value : (input.value === '' ? '' : Number(input.value));
          invalidate();
          if (k === 'id') renderStartOptions();
          saveDraft();
        });
      });
      row.querySelector('button').addEventListener('click', () => {
        if (frame.length <= 2) return;
        frame.splice(i, 1);
        invalidate();
        renderAll();
      });
      card.appendChild(row);
    });

    if (frame.length < 8) {
      const add = document.createElement('button');
      add.type = 'button';
      add.className = 'add-spot';
      add.textContent = '＋ 添加斑点';
      add.addEventListener('click', () => {
        frame.push({ id: nextId(t), x: '', y: '', b: '' });
        invalidate();
        renderAll();
      });
      card.appendChild(add);
    }
    els.frames.appendChild(card);
  });
  els.addFrame.disabled = draft.frames.length >= 7;
}

function nextId(frameIdx) {
  const used = new Set();
  draft.frames.forEach((fr) => fr.forEach((s) => used.add(s.id)));
  let n = frameIdx + 1;
  let cand = `F${n}-1`;
  let k = 1;
  while (used.has(cand)) { k++; cand = `F${n}-${k}`; }
  return cand;
}

function renderStartOptions() {
  const cur = draft.startId;
  els.startSelect.innerHTML = '';
  draft.frames[0].forEach((s) => {
    const opt = document.createElement('option');
    opt.value = s.id;
    opt.textContent = s.id || '（空编号）';
    els.startSelect.appendChild(opt);
  });
  if ([...els.startSelect.options].some((o) => o.value === cur)) {
    els.startSelect.value = cur;
  } else {
    draft.startId = els.startSelect.value;
  }
}

function renderParams() {
  els.maxDist.value = draft.maxDist;
  els.maxSkip.value = draft.maxSkip;
  els.target.value = draft.target;
  renderStartOptions();
}

function renderAll() {
  renderFrames();
  renderParams();
  saveDraft();
}

function escapeAttr(s) {
  return String(s ?? '').replace(/[&<>"]/g, (c) =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

// ---- 编辑即作废旧谱系 ----
function invalidate() {
  if (lastSolution) {
    lastSolution = null;
    els.ok.classList.add('hidden');
    els.infeasible.classList.add('hidden');
    els.empty.classList.remove('hidden');
  }
  els.stale.classList.remove('hidden');
  els.errors.classList.add('hidden');
  saveDraft();
}

// ---- 复原 ----
function solve() {
  els.errors.classList.add('hidden');
  els.stale.classList.add('hidden');

  const { errors, spec } = normalizeSpec(draft);
  if (errors.length) {
    showErrors(errors);
    els.ok.classList.add('hidden');
    els.infeasible.classList.add('hidden');
    els.empty.classList.remove('hidden');
    return;
  }

  // 先让浏览器完成重绘再进行可能较久的计算
  els.solve.disabled = true;
  els.solve.textContent = '求解中…';
  requestAnimationFrame(() => setTimeout(() => {
    try {
      const raw = solveLineage(spec);
      const sol = presentSolution(spec, raw);
      lastSolution = sol;

      els.empty.classList.add('hidden');
      els.ok.classList.add('hidden');
      els.infeasible.classList.add('hidden');

      if (!sol.feasible) {
        els.infeasible.innerHTML = `
          <strong>不存在可同时满足全部约束的谱系。</strong><br/>
          终帧存活数与祖先唯一性、位移或漏检限制无法同时成立；输入草稿已保留，可调整后再次复原。<br/>
          最早断开的帧间：<span class="break-frame">第 ${sol.earliestBreak.from + 1} 帧 → 第 ${sol.earliestBreak.to + 1} 帧</span>`;
        els.infeasible.classList.remove('hidden');
        return;
      }
      renderSolution(sol);
    } finally {
      els.solve.disabled = false;
      els.solve.textContent = '复　原';
    }
  }, 16));
}

function showErrors(errors) {
  els.errors.innerHTML = '<strong>输入有误：</strong><ul>' +
    errors.map((e) => `<li>${escapeAttr(e.message)}</li>`).join('') + '</ul>';
  els.errors.classList.remove('hidden');
}

// ---- 结果渲染 ----
function renderSolution(sol) {
  els.stats.innerHTML = [
    ['总亮度', sol.totalBrightness],
    ['漏检段', sol.skips],
    ['分裂次数', sol.divisions],
    ['终帧存活', sol.survivors],
  ].map(([k, v]) => `<span class="stat">${k}<b>${v}</b></span>`).join('');

  renderChart(sol);
  renderEdgesTable(sol);
  renderUsed(sol);
  els.ok.classList.remove('hidden');
}

function lookup() {
  const m = new Map();
  draft.frames.forEach((fr, t) => fr.forEach((s) => m.set(`${t}:${s.id}`, s)));
  return m;
}

function renderEdgesTable(sol) {
  const spot = lookup();
  const childCount = new Map();
  sol.edges.forEach((e) => childCount.set(`${e.fromFrame}:${e.fromId}`,
    (childCount.get(`${e.fromFrame}:${e.fromId}`) || 0) + 1));

  els.tbody.innerHTML = '';
  sol.edges.forEach((e, idx) => {
    const m = spot.get(`${e.fromFrame}:${e.fromId}`);
    const c = spot.get(`${e.toFrame}:${e.toId}`);
    const twins = childCount.get(`${e.fromFrame}:${e.fromId}`) === 2;
    let kind, badge, rowClass;
    if (e.gap === 2) {
      kind = '漏检补获'; badge = 'skip'; rowClass = 'skip-row';
    } else if (twins) {
      kind = '分裂'; badge = 'divide'; rowClass = 'divide-row';
    } else {
      kind = '保持'; badge = 'keep'; rowClass = '';
    }
    const tr = document.createElement('tr');
    tr.className = rowClass;
    tr.innerHTML = `
      <td>${idx + 1}</td>
      <td>F${e.fromFrame + 1}·${escapeAttr(e.fromId)} <em>(${m.x},${m.y})</em></td>
      <td>F${e.toFrame + 1}·${escapeAttr(e.toId)} <em>(${c.x},${c.y})</em></td>
      <td><span class="badge ${badge}">${kind}</span></td>
      <td>F${e.fromFrame + 1}→F${e.toFrame + 1}${e.gap === 2 ? '（漏 1 帧）' : ''}</td>
      <td>${e.dist}</td>`;
    els.tbody.appendChild(tr);
  });
}

function renderUsed(sol) {
  els.used.innerHTML = '';
  draft.frames.forEach((fr, t) => {
    const used = sol.used[t] || [];
    const unused = fr.map((s) => s.id).filter((id) => !used.includes(id));
    const div = document.createElement('div');
    div.className = 'used-frame';
    div.innerHTML =
      `<b>第 ${t + 1} 帧（${used.length}/${fr.length}）</b>` +
      used.map((id) => `<span>${escapeAttr(id)}</span>`).join('') +
      (unused.length ? `<em>未采用：${unused.map(escapeAttr).join('、')}</em>` : '');
    els.used.appendChild(div);
  });
}

function renderChart(sol) {
  const W = 760;
  const H = 360;
  const ML = 46, MR = 30, MT = 26, MB = 30;
  const F = draft.frames.length;

  const xs = draft.frames.flatMap((fr) => fr.map((s) => Number(s.x)));
  const ys = draft.frames.flatMap((fr) => fr.map((s) => Number(s.y)));
  const xmin = Math.min(...xs), xmax = Math.max(...xs);
  const ymin = Math.min(...ys), ymax = Math.max(...ys);
  const spanX = Math.max(1, xmax - xmin);
  const spanY = Math.max(1, ymax - ymin);

  const colX = (t) => ML + (t * (W - ML - MR)) / Math.max(1, F - 1);
  // 帧内横向散开：避免同帧斑点重叠
  const jitter = new Map();
  draft.frames.forEach((fr, t) => {
    fr.forEach((s, i) => jitter.set(`${t}:${s.id}`, (i - (fr.length - 1) / 2) * 4));
  });
  const px = (t, s) => colX(t) + jitter.get(`${t}:${s.id}`);
  const py = (s) => MT + ((ymax - Number(s.y)) / spanY) * (H - MT - MB)
    + ((Number(s.x) - xmin) / spanX) * 10 - 5;

  const usedSet = new Set();
  sol.used.forEach((ids, t) => ids.forEach((id) => usedSet.add(`${t}:${id}`)));

  let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="谱系图">`;
  for (let t = 0; t < F; t++) {
    svg += `<line x1="${colX(t)}" y1="${MT - 12}" x2="${colX(t)}" y2="${H - MB + 12}" stroke="#22304a" stroke-width="1"/>`;
    svg += `<text x="${colX(t)}" y="${H - 8}" fill="#9fb2c8" font-size="12" text-anchor="middle">第 ${t + 1} 帧</text>`;
  }

  const spot = lookup();
  sol.edges.forEach((e) => {
    const m = spot.get(`${e.fromFrame}:${e.fromId}`);
    const c = spot.get(`${e.toFrame}:${e.toId}`);
    const x1 = px(e.fromFrame, m), y1 = py(m);
    const x2 = px(e.toFrame, c), y2 = py(c);
    const mx = (x1 + x2) / 2;
    const my = (y1 + y2) / 2;
    if (e.gap === 2) {
      svg += `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}"
        stroke="#f6ad55" stroke-width="2" stroke-dasharray="6 4"/>`;
    } else {
      const cx = mx, cy = my - 6;
      svg += `<path d="M${x1},${y1} Q${cx},${cy} ${x2},${y2}"
        fill="none" stroke="#4fd1c5" stroke-width="1.8" opacity="0.85"/>`;
    }
    svg += `<text x="${mx}" y="${my - 8}" fill="#8aa0b8" font-size="9.5" text-anchor="middle">${e.dist}</text>`;
  });

  draft.frames.forEach((fr, t) => fr.forEach((s) => {
    const x = px(t, s);
    const y = py(s);
    const used = usedSet.has(`${t}:${s.id}`);
    const isStart = t === 0 && s.id === draft.startId;
    const r = 4 + Math.max(0, Number(s.b) || 0) / 40;
    svg += `<circle cx="${x}" cy="${y}" r="${Math.min(r, 9)}"
      fill="${used ? '#4fd1c5' : '#46566c'}" fill-opacity="${used ? 0.95 : 0.5}">
      <title>F${t + 1}·${escapeAttr(s.id)} (${s.x}, ${s.y}) 亮度 ${s.b}${used ? ' · 采用' : ' · 未采用'}</title></circle>`;
    if (isStart) {
      svg += `<circle cx="${x}" cy="${y}" r="${Math.min(r, 9) + 4}"
        fill="none" stroke="#f6e05e" stroke-width="2"/>`;
    }
    svg += `<text x="${x}" y="${y - Math.min(r, 9) - 4}" fill="${used ? '#cde8e5' : '#7488a0'}"
      font-size="9.5" text-anchor="middle">${escapeAttr(s.id)}</text>`;
  }));
  svg += '</svg>';
  els.chart.innerHTML = svg;
}

// ---- 事件绑定 ----
els.addFrame.addEventListener('click', () => {
  if (draft.frames.length >= 7) return;
  const t = draft.frames.length;
  draft.frames.push([
    { id: nextId(t), x: '', y: '', b: '' },
    { id: nextId(t), x: '', y: '', b: '' },
  ]);
  invalidate();
  renderAll();
});

els.sample.addEventListener('click', () => {
  draft = structuredClone(SAMPLE);
  invalidate();
  renderAll();
});

els.startSelect.addEventListener('change', () => {
  draft.startId = els.startSelect.value;
  invalidate();
});
[['maxDist', els.maxDist], ['maxSkip', els.maxSkip], ['target', els.target]]
  .forEach(([k, input]) => {
    input.addEventListener('input', () => {
      draft[k] = input.value === '' ? '' : Number(input.value);
      invalidate();
    });
  });

els.solve.addEventListener('click', solve);

renderAll();
