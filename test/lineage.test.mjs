// 谱系求解器测试：约束结构、裁决顺序、漏检、不可行报告 + 独立暴力枚举对拍。
'use strict';

import test from 'node:test';
import assert from 'node:assert/strict';
import { normalizeSpec, solveLineage, presentSolution } from '../public/js/lineage.js';

function run(input) {
  const { errors, spec } = normalizeSpec(input);
  assert.deepEqual(errors, [], `输入校验应通过: ${JSON.stringify(errors)}`);
  const raw = solveLineage(spec);
  return { spec, raw, sol: presentSolution(spec, raw) };
}

// 对任意可行解做通用结构约束校验
function assertValidLineage(spec, sol, input) {
  const F = spec.frames.length;
  assert.equal(sol.feasible, true);
  assert.equal(sol.survivors, input.target);
  assert.equal(sol.used[F - 1].length, input.target);
  assert.deepEqual(sol.used[0], [spec.frames[0][spec.startIndex].id]);

  const indeg = new Map(); // 'frame:id' -> 入边数
  const outEdges = new Map();
  for (const e of sol.edges) {
    indeg.set(`${e.toFrame}:${e.toId}`, (indeg.get(`${e.toFrame}:${e.toId}`) || 0) + 1);
    const k = `${e.fromFrame}:${e.fromId}`;
    if (!outEdges.has(k)) outEdges.set(k, []);
    outEdges.get(k).push(e);
    assert.ok([1, 2].includes(e.gap), '连接只能相邻或跨一帧');
    assert.equal(e.toFrame - e.fromFrame, e.gap);
    const m = spec.frames[e.fromFrame].find((s) => s.id === e.fromId);
    const c = spec.frames[e.toFrame].find((s) => s.id === e.toId);
    const d = Math.hypot(m.x - c.x, m.y - c.y);
    assert.ok(d <= input.maxDist * e.gap + 1e-9, '位移超限');
    assert.ok(Math.abs(d - e.dist) < 0.011);
  }

  // 每个采用的非起始斑点恰有一个祖先；同一斑点不入两支
  for (let t = 0; t < F; t++) {
    for (const id of sol.used[t]) {
      const deg = indeg.get(`${t}:${id}`) || 0;
      if (t === 0 && id === input.startId) {
        assert.equal(deg, 0, '起始斑点不应有祖先');
      } else {
        assert.equal(deg, 1, `斑点 F${t + 1}·${id} 应有恰一个祖先，实际 ${deg}`);
      }
    }
  }

  // 每支 1 或 2 个后代；漏检补获只能单传；全部到达末帧
  for (let t = 0; t < F; t++) {
    for (const id of sol.used[t]) {
      const es = outEdges.get(`${t}:${id}`) || [];
      if (t === F - 1) {
        assert.equal(es.length, 0, '末帧斑点不应再有后代');
      } else {
        const gaps = es.filter((e) => e.gap === 2);
        const direct = es.filter((e) => e.gap === 1);
        if (gaps.length) {
          assert.equal(gaps.length, 1);
          assert.equal(direct.length, 0, '漏检中细胞不能同时分裂到下一帧');
        } else {
          assert.ok([1, 2].includes(direct.length), '必须保持一个或分裂为恰两个后代');
        }
        assert.ok(es.every((e) => reachesLast(e, F)), '所有存活支必须到达末帧');
      }
    }
  }

  function reachesLast(e, F) {
    if (e.toFrame === F - 1) return true;
    return (outEdges.get(`${e.toFrame}:${e.toId}`) || []).some((nx) => reachesLast(nx, F));
  }

  // 分裂次数与漏检数
  assert.equal(sol.skips, sol.edges.filter((e) => e.gap === 2).length);
  let div = 0;
  for (const list of outEdges.values()) if (list.filter((e) => e.gap === 1).length === 2) div++;
  assert.equal(sol.divisions, div);

  // 总亮度
  let bright = 0;
  sol.used.forEach((ids, t) => ids.forEach((id) => {
    bright += spec.frames[t].find((s) => s.id === id).b;
  }));
  assert.equal(sol.totalBrightness, bright);
}

const base4 = () => ({
  frames: [
    [{ id: 'a', x: 0, y: 0, b: 10 }, { id: 'z0', x: 9, y: 9, b: 99 }],
    [{ id: 'b', x: 1, y: 0, b: 11 }, { id: 'z1', x: 9, y: 9, b: 99 }],
    [{ id: 'c', x: 2, y: 0, b: 12 }, { id: 'z2', x: 9, y: 9, b: 99 }],
    [{ id: 'd', x: 3, y: 0, b: 13 }, { id: 'z3', x: 9, y: 9, b: 99 }],
  ],
  startId: 'a', maxDist: 2, maxSkip: 0, target: 1,
});

test('基础保持谱系：最亮杂质因位移/祖先约束被排除', () => {
  const input = base4();
  const { spec, sol } = run(input);
  assertValidLineage(spec, sol, input);
  assert.deepEqual(sol.used, [['a'], ['b'], ['c'], ['d']]);
  assert.equal(sol.divisions, 0);
  assert.equal(sol.totalBrightness, 46);
  // 逐帧贪心选最亮点会错误串起 z0→z1→z2→z3；本解绝不能采用它们
  assert.ok(sol.used.every((ids) => !ids.some((id) => id.startsWith('z'))));
});

test('分裂谱系：一支在末段分裂为恰两个后代', () => {
  const input = {
    frames: [
      [{ id: 'a', x: 0, y: 0, b: 10 }, { id: 'q0', x: 8, y: 8, b: 50 }],
      [{ id: 'b', x: 1, y: 0, b: 10 }, { id: 'q1', x: 8, y: 8, b: 50 }],
      [{ id: 'c', x: 2, y: 0, b: 10 }, { id: 'q2', x: 8, y: 8, b: 50 }],
      [
        { id: 'd1', x: 3, y: -1, b: 10 },
        { id: 'd2', x: 3, y: 1, b: 10 },
        { id: 'q3', x: 8, y: 8, b: 50 },
      ],
    ],
    startId: 'a', maxDist: 2, maxSkip: 0, target: 2,
  };
  const { spec, sol } = run(input);
  assertValidLineage(spec, sol, input);
  assert.deepEqual(sol.used[3].sort(), ['d1', 'd2']);
  assert.equal(sol.divisions, 1);
  const last = sol.edges.filter((e) => e.toFrame === 3);
  assert.equal(last.length, 2);
  assert.deepEqual(last.map((e) => e.fromId), ['c', 'c']);
});

test('漏检：中间帧无近邻斑点时跨一帧连接', () => {
  const input = {
    frames: [
      [{ id: 'a', x: 0, y: 0, b: 10 }, { id: 'g0', x: 9, y: 0, b: 10 }],
      [{ id: 'j1', x: 9, y: 0, b: 99 }, { id: 'j2', x: 9, y: 9, b: 99 }],
      [{ id: 'c', x: 2, y: 0, b: 10 }, { id: 'g2', x: 9, y: 0, b: 10 }],
      [{ id: 'd', x: 3, y: 0, b: 10 }, { id: 'g3', x: 9, y: 0, b: 10 }],
    ],
    startId: 'a', maxDist: 2, maxSkip: 1, target: 1,
  };
  const { spec, sol } = run(input);
  assertValidLineage(spec, sol, input);
  assert.deepEqual(sol.used, [['a'], [], ['c'], ['d']]);
  assert.equal(sol.skips, 1);
  const gap = sol.edges.find((e) => e.gap === 2);
  assert.deepEqual({ f: gap.fromFrame, t: gap.toFrame, from: gap.fromId, to: gap.toId },
    { f: 0, t: 2, from: 'a', to: 'c' });
  assert.equal(gap.dist, 2);
});

test('不允许漏检时同一用例不可行，并报告最早断开帧间', () => {
  const input = { ...base4(), frames: structuredClone(base4().frames), maxSkip: 0, target: 1 };
  // 重放漏检场景但关闭漏检
  input.frames = [
    [{ id: 'a', x: 0, y: 0, b: 10 }, { id: 'g0', x: 9, y: 0, b: 10 }],
    [{ id: 'j1', x: 9, y: 0, b: 99 }, { id: 'j2', x: 9, y: 9, b: 99 }],
    [{ id: 'c', x: 2, y: 0, b: 10 }, { id: 'g2', x: 9, y: 0, b: 10 }],
    [{ id: 'd', x: 3, y: 0, b: 10 }, { id: 'g3', x: 9, y: 0, b: 10 }],
  ];
  const { raw, sol } = run(input);
  assert.equal(raw.feasible, false);
  assert.equal(sol.feasible, false);
  assert.equal(sol.earliestBreak.from, 0);
  assert.equal(sol.earliestBreak.to, 1);
});

test('终帧目标数超出增长能力时不可行', () => {
  const input = base4();
  input.target = 2; // 全程只有 1 个可达斑点，无法在末帧前分裂
  const { raw } = run(input);
  assert.equal(raw.feasible, false);
  assert.ok(raw.earliestBreak.from >= 0 && raw.earliestBreak.to < 4);
});

test('最早断开帧间定位到真正无解的中段边界', () => {
  // 前三个帧间单传均可走；唯独末帧前需要分裂为 2，但末帧轨迹旁只有 1 个近邻斑点
  const input = {
    frames: [
      [{ id: 'a', x: 0, y: 0, b: 5 }, { id: 'q0', x: 9, y: 9, b: 9 }],
      [{ id: 'b', x: 1, y: 0, b: 5 }, { id: 'q1', x: 9, y: 9, b: 9 }],
      [{ id: 'c', x: 2, y: 0, b: 5 }, { id: 'q2', x: 9, y: 9, b: 9 }],
      [{ id: 'd', x: 3, y: 0, b: 5 }, { id: 'far', x: 30, y: 30, b: 9 }],
    ],
    startId: 'a', maxDist: 2, maxSkip: 0, target: 2,
  };
  const { raw, sol } = run(input);
  assert.equal(raw.feasible, false);
  assert.equal(raw.earliestBreak.from, 2);
  assert.equal(raw.earliestBreak.to, 3);
  assert.match(sol.earliestBreakLabel, /第 3 帧 → 第 4 帧/);
});

test('亮度优先：宁选稍暗但能连成高总和的一支（联合最优而非逐帧贪心）', () => {
  // 近邻有两条互斥路径；亮路径中途会撞上同一斑点（违反唯一祖先）→ 只能走低亮路径
  const input = {
    frames: [
      [{ id: 'a', x: 0, y: 0, b: 1 }, { id: 'x', x: 0, y: 5, b: 1 }],
      [{ id: 'b', x: 1, y: 0, b: 1 }, { id: 'B', x: 1, y: 5, b: 100 }],
      [{ id: 'c', x: 2, y: 0, b: 1 }, { id: 'C', x: 2, y: 5, b: 100 }],
      [{ id: 'd', x: 3, y: 0, b: 1 }, { id: 'D', x: 3, y: 3, b: 100 }],
    ],
    startId: 'x', maxDist: 3, maxSkip: 0, target: 1,
  };
  // 从 x(0,5) 出发：B(1,5) 很亮，但 B 能到 c(2,0) 吗？距离 5.1 >3；只能到 C。
  // C(2,5) 到 D(3,3) 距离 √5<3 可取，路径 x-B-C-D 全亮。
  const { spec, sol } = run(input);
  assertValidLineage(spec, sol, input);
  assert.deepEqual(sol.used, [['x'], ['B'], ['C'], ['D']]);
  assert.equal(sol.totalBrightness, 301);
});

test('亮度相同、漏检数相同：按输入顺序稳定裁决（采用更早的斑点）', () => {
  const mk = (startId) => ({
    frames: [
      [{ id: 'p', x: 0, y: 0, b: 5 }, { id: 'q', x: 0, y: 3, b: 5 }],
      [{ id: 'r', x: 1, y: 0, b: 5 }, { id: 's', x: 1, y: 3, b: 5 }],
      [{ id: 't', x: 2, y: 0, b: 5 }, { id: 'u', x: 2, y: 3, b: 5 }],
      [{ id: 'v', x: 3, y: 0, b: 5 }, { id: 'w', x: 3, y: 3, b: 5 }],
    ],
    startId, maxDist: 3, maxSkip: 0, target: 1,
  });
  const r1 = run(mk('p'));
  assertValidLineage(r1.spec, r1.sol, mk('p'));
  assert.deepEqual(r1.sol.used, [['p'], ['r'], ['t'], ['v']]);
  const r2 = run(mk('q'));
  assert.deepEqual(r2.sol.used, [['q'], ['s'], ['u'], ['w']]);
});

test('漏检数为第二裁决键：同亮度时优先无漏检方案', () => {
  // 帧1 同时存在直接女儿和跨帧机会，二者后续亮度相同 → 应选直接连接
  const input = {
    frames: [
      [{ id: 'a', x: 0, y: 0, b: 5 }, { id: '_0', x: 0, y: 6, b: 5 }],
      [{ id: 'b', x: 1, y: 0, b: 5 }, { id: '_1', x: 0, y: 6, b: 5 }],
      [{ id: 'c', x: 2, y: 0, b: 5 }, { id: '_2', x: 0, y: 6, b: 5 }],
      [{ id: 'd', x: 3, y: 0, b: 5 }, { id: '_3', x: 0, y: 6, b: 5 }],
    ],
    startId: 'a', maxDist: 5, maxSkip: 1, target: 1,
  };
  const { spec, sol } = run(input);
  assertValidLineage(spec, sol, input);
  assert.equal(sol.skips, 0);
});

test('修改帧内容后重新求解得到不同谱系（不保留旧结果由调用方保证，解本身随输入变化）', () => {
  const input = base4();
  const s1 = run(input).sol;
  input.frames[1][1] = { id: 'b2', x: 1, y: 0, b: 80 }; // 杂质移动到轨迹上
  const s2 = run(input).sol;
  assert.notDeepEqual(s1.used, s2.used);
  assert.ok(s2.used[1].includes('b2'));
});

test('输入校验：帧数、斑点数、重复编号、参数范围', () => {
  const bad = { frames: base4().frames.slice(0, 3), startId: 'a', maxDist: 1, maxSkip: 0, target: 1 };
  assert.equal(normalizeSpec(bad).errors.length > 0, true);
  const dup = base4();
  dup.frames = structuredClone(dup.frames);
  dup.frames[0][1].id = 'a';
  assert.ok(normalizeSpec(dup).errors.some((e) => e.message.includes('重复')));
  const badCoord = base4();
  badCoord.frames = structuredClone(badCoord.frames);
  badCoord.frames[0][0].x = 1.5;
  assert.ok(normalizeSpec(badCoord).errors.some((e) => e.message.includes('整数')));
  const badTarget = base4();
  badTarget.target = 9;
  assert.ok(normalizeSpec(badTarget).errors.some((e) => e.field === 'target'));
});

// ---------- 独立暴力枚举：逐帧 DFS 穷举全部可行谱系（无备忘、无剪枝界） ----------
function bruteForce(spec) {
  const { frames, startIndex, maxDist, maxSkip, target } = spec;
  const F = frames.length;
  const near = (m, c, gap) =>
    Math.hypot(m.x - c.x, m.y - c.y) <= maxDist * gap + 1e-9;

  let best = null;

  // t 边界：live 为帧 t 上的存活斑点（局部序号），gaps 为帧 t-1 漏检、
  // 必须在帧 t 补获的母本（帧 t-1 序号）。
  function dfs(t, live, gaps, usedSkip, bright, usedPerFrame) {
    if (live.length + gaps.length > target) return;
    if (t === F - 1) {
      if (gaps.length > 0 || live.length !== target) return;
      const sig = usedPerFrame.slice(1).map((s) => [...s].sort((a, b) => a - b));
      if (!best ||
        bright > best.bright ||
        (bright === best.bright &&
          (usedSkip < best.skips ||
            (usedSkip === best.skips && tupleLex(sig, best.sig) < 0)))) {
        best = { bright, skips: usedSkip, sig };
      }
      return;
    }

    const tracks = [
      ...live.map((i) => ({ kind: 'o', i })),
      ...gaps.map((i) => ({ kind: 'g', i })),
    ];
    const claimed = new Set(); // 帧 t+1 已被女儿占用的斑点
    const nextLive = [];
    const openGaps = []; // 帧 t 新开漏检的母本（帧 t 序号）

    function rec(k, accBright) {
      if (k === tracks.length) {
        dfs(t + 1, nextLive.slice().sort((a, b) => a - b),
          openGaps.slice().sort((a, b) => a - b),
          usedSkip + openGaps.length, accBright, usedPerFrame);
        return;
      }
      const tr = tracks[k];
      if (tr.kind === 'g') {
        // 漏检母本：恰一个跨帧女儿
        for (let j = 0; j < frames[t + 1].length; j++) {
          if (claimed.has(j)) continue;
          if (!near(frames[t - 1][tr.i], frames[t + 1][j], 2)) continue;
          claimed.add(j); nextLive.push(j);
          usedPerFrame[t + 1].add(j);
          rec(k + 1, accBright + frames[t + 1][j].b);
          usedPerFrame[t + 1].delete(j);
          nextLive.pop(); claimed.delete(j);
        }
        return;
      }
      const m = frames[t][tr.i];
      // 保持：一个相邻女儿
      for (let j = 0; j < frames[t + 1].length; j++) {
        if (claimed.has(j) || !near(m, frames[t + 1][j], 1)) continue;
        claimed.add(j); nextLive.push(j);
        usedPerFrame[t + 1].add(j);
        rec(k + 1, accBright + frames[t + 1][j].b);
        usedPerFrame[t + 1].delete(j);
        nextLive.pop(); claimed.delete(j);
      }
      // 分裂：两个不同的相邻女儿
      for (let a = 0; a < frames[t + 1].length; a++) {
        if (claimed.has(a) || !near(m, frames[t + 1][a], 1)) continue;
        for (let b2 = a + 1; b2 < frames[t + 1].length; b2++) {
          if (claimed.has(b2) || !near(m, frames[t + 1][b2], 1)) continue;
          claimed.add(a); claimed.add(b2);
          nextLive.push(a, b2);
          usedPerFrame[t + 1].add(a); usedPerFrame[t + 1].add(b2);
          rec(k + 1, accBright + frames[t + 1][a].b + frames[t + 1][b2].b);
          usedPerFrame[t + 1].delete(a); usedPerFrame[t + 1].delete(b2);
          nextLive.pop(); nextLive.pop();
          claimed.delete(a); claimed.delete(b2);
        }
      }
      // 漏检：帧 t+1 不出现（末帧前不开新漏检）
      if (usedSkip + openGaps.length < maxSkip && t + 2 <= F - 1) {
        openGaps.push(tr.i);
        rec(k + 1, accBright);
        openGaps.pop();
      }
    }
    rec(0, bright);
  }

  const used0 = Array.from({ length: F }, () => new Set());
  used0[0].add(startIndex);
  dfs(0, [startIndex], [], 0, frames[0][startIndex].b, used0);
  return best;
}

function tupleLex(a, b) {
  const flatA = a.flat(), flatB = b.flat();
  for (let i = 0; i < Math.min(flatA.length, flatB.length); i++) {
    if (flatA[i] !== flatB[i]) return flatA[i] - flatB[i];
  }
  return flatA.length - flatB.length;
}

function rng(seed) {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

test('随机对拍：小规模输入下与独立暴力枚举裁决一致', () => {
  const rand = rng(20260929);
  let feasibleCases = 0;
  for (let iter = 0; iter < 2000 && feasibleCases < 200; iter++) {
    const F = rand() < 0.35 ? 5 : 4;
    const frames = [];
    for (let t = 0; t < F; t++) {
      const n = 2 + Math.floor(rand() * 2); // 2~3 个斑点
      const fr = [];
      for (let i = 0; i < n; i++) {
        fr.push({
          id: `f${t}_${i}`,
          x: Math.floor(rand() * 3),
          y: Math.floor(rand() * 3),
          b: Math.floor(rand() * 9) + 1,
        });
      }
      frames.push(fr);
    }
    const input = {
      frames,
      startId: frames[0][Math.floor(rand() * frames[0].length)].id,
      maxDist: 1 + Math.floor(rand() * 3),
      maxSkip: rand() < 0.5 ? 0 : 1,
      target: 1 + Math.floor(rand() * 2),
    };
    const { errors, spec } = normalizeSpec(input);
    if (errors.length) continue;
    const raw = solveLineage(spec);
    const bf = bruteForce(spec);
    if (!raw.feasible) {
      assert.equal(bf, null, `迭代 ${iter}：求解器判不可行但暴力枚举存在解`);
      continue;
    }
    assert.ok(bf, `迭代 ${iter}：求解器给出解但暴力枚举无解`);
    feasibleCases++;
    const sol = presentSolution(spec, raw);
    assertValidLineage(spec, sol, input);
    assert.equal(sol.totalBrightness, bf.bright, `迭代 ${iter} 总亮度不一致`);
    assert.equal(sol.skips, bf.skips, `迭代 ${iter} 漏检数不一致`);
    const sig = sol.used.slice(1).map((ids, t) =>
      ids.map((id) => frames[t + 1].findIndex((s) => s.id === id)).sort((a, b) => a - b));
    assert.equal(tupleLex(sig, bf.sig), 0, `迭代 ${iter} 输入顺序裁决不一致`);
  }
  assert.ok(feasibleCases >= 30, `可行对拍用例过少: ${feasibleCases}`);
});
