// 藻类细胞分裂谱系复原核心（纯函数，零依赖，浏览器 / Node 共用）
//
// 模型：
//  - 帧按时刻排列；每帧若干斑点（唯一 id、整数坐标、整数亮度）。
//  - 连接只允许相邻帧（gap=1）或跨越恰好一帧漏检（gap=2）。
//  - 一个细胞要么保持为一个后代，要么分裂为恰两个后代；不允许消亡。
//  - 每个非起始斑点恰有一个祖先（入边），同一斑点不得被两支共用。
//  - 所有存活支必须从起始斑点出发到达末帧，且末帧存活数恰为目标数。
//  - 裁决顺序：总亮度最高 → 漏检段最少 → 输入顺序（逐帧采用斑点局部序号，
//    再逐斑点母本全局序号）字典序稳定裁决。
//
// 三阶段精确求解（不缩小支持范围、不近似、不改裁决语义）：
//  阶段一 边界联合转移对 transitions(t, live, gaps)：存活母本（帧 t）与待
//    补获漏检母本（帧 t-1）在帧 t+1 上联合安排女儿，内层 DP 枚举
//    （女儿掩码, 新开漏检掩码）全部可行对。结果只依赖 (t, live, gaps)，
//    以紧凑 Uint32Array 缓存（每对一个整数，不存母本向量）。
//  阶段二 取值 DP valueDP(t, live, gaps, left)：严格按
//    (t, live, gaps, 剩余漏检额度) 备忘，每状态只存
//    （最高亮度, 最少漏检）两个小整数；配合乐观亮度上界与漏检前缀严格
//    剪枝（相等才可能并列，绝不误杀）。
//  阶段三 精确残量可行性 feas(t,…,B,S)：后缀能否恰好取得亮度 B、漏检 S，
//    备忘为布尔；随后沿最优路径在每个边界只对「该状态」的转移对按输入
//    顺序排序，取第一个残量可行的对——逐帧贪心取字典序最小即全局字典序
//    最小。被选中对的母本配对由受限于该 (女儿掩码, 漏检掩码) 的小型 DP
//    按需求出（同样只保留字典序最小配对）。

'use strict';

// 全局斑点总数 ≤ 7×8 = 56，母本序号 0..55 用 6 bit 表示；
// 未分配用 0，已分配母本 g 编码为 g+1。
const MOM_BITS = 6;
const MOM_MASK = (1 << MOM_BITS) - 1;
const momGet = (packed, j) => (((packed / 2 ** (MOM_BITS * j)) | 0) & MOM_MASK) - 1;
const momSet = (packed, j, g) => packed + (g + 1) * 2 ** (MOM_BITS * j);

/**
 * 校验并规范化输入。
 * @returns {{errors:Array<{field:string,message:string}>, spec:object|null}}
 */
export function normalizeSpec(raw) {
  const errors = [];
  const field = (name, message) => errors.push({ field: name, message });

  if (!raw || typeof raw !== 'object' || !Array.isArray(raw.frames)) {
    return { errors: [{ field: 'frames', message: '缺少帧数据' }], spec: null };
  }
  const F = raw.frames.length;
  if (F < 4 || F > 7) {
    field('frames', `帧数必须在 4 至 7 之间（当前 ${F}）`);
  }

  const frames = [];
  raw.frames.forEach((fr, t) => {
    const out = [];
    if (!Array.isArray(fr) || fr.length < 2 || fr.length > 8) {
      field(`frame${t}`, `第 ${t + 1} 帧斑点数必须在 2 至 8 之间（当前 ${Array.isArray(fr) ? fr.length : 0}）`);
      return;
    }
    const seen = new Set();
    fr.forEach((s, j) => {
      const label = `第 ${t + 1} 帧斑点 ${j + 1}`;
      if (!s || typeof s.id !== 'string' || s.id.trim() === '') {
        field(`frame${t}`, `${label} 缺少唯一编号`);
        return;
      }
      const id = s.id.trim();
      if (seen.has(id)) {
        field(`frame${t}`, `第 ${t + 1} 帧内斑点编号重复：${id}`);
        return;
      }
      seen.add(id);
      const x = Number(s.x);
      const y = Number(s.y);
      const b = Number(s.b);
      if (!Number.isInteger(x) || !Number.isInteger(y)) {
        field(`frame${t}`, `${label}（${id}）坐标必须为整数`);
        return;
      }
      if (!Number.isInteger(b) || b < 0) {
        field(`frame${t}`, `${label}（${id}）亮度必须为非负整数`);
        return;
      }
      out.push({ id, x, y, b });
    });
    frames.push(out);
  });

  if (errors.length) return { errors, spec: null };

  const startId = typeof raw.startId === 'string' ? raw.startId.trim() : '';
  const startIndex = frames[0] ? frames[0].findIndex((s) => s.id === startId) : -1;
  if (startIndex < 0) {
    field('startId', `起始斑点必须是第 1 帧中存在的编号（当前“${raw.startId}”）`);
  }

  const maxDist = Number(raw.maxDist);
  if (!Number.isFinite(maxDist) || maxDist < 0) {
    field('maxDist', '相邻帧最大位移必须为非负数');
  }

  const maxSkip = Number(raw.maxSkip);
  if (!Number.isInteger(maxSkip) || maxSkip < 0 || maxSkip > F - 2) {
    field('maxSkip', `允许漏检帧数必须为 0 至 ${Math.max(0, F - 2)} 的整数`);
  }

  const lastSize = frames[F - 1] ? frames[F - 1].length : 0;
  const target = Number(raw.target);
  if (!Number.isInteger(target) || target < 1 || target > lastSize) {
    field('target', `终帧存活细胞数必须为 1 至末帧斑点数（${lastSize}）的整数`);
  }

  if (errors.length) return { errors, spec: null };
  return {
    errors: [],
    spec: { frames, startIndex, maxDist, maxSkip, target },
  };
}

function compareList(a, b) {
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) {
    if (a[i] !== b[i]) return a[i] < b[i] ? -1 : 1;
  }
  return a.length - b.length;
}

// 8 位以内掩码的公共查表
const PC = new Uint8Array(256); // 置位数
const BITS = new Array(256); // 置位序号列表
const CTZ = new Uint8Array(256); // 最低置位序号
for (let m = 0; m < 256; m++) {
  const ls = [];
  for (let i = 0, v = m; v; i++, v >>>= 1) {
    if (v & 1) ls.push(i);
  }
  BITS[m] = ls;
  PC[m] = ls.length;
  CTZ[m] = ls.length ? ls[0] : 0;
}

/**
 * 求解谱系。
 * @returns {object} 可行时 {feasible:true, ...}；不可行时
 *   {feasible:false, earliestBreak:{from:number,to:number}}
 */
export function solveLineage(spec) {
  const { frames, startIndex, maxDist, maxSkip, target } = spec;
  const F = frames.length;
  const sizes = frames.map((fr) => fr.length);

  const offset = [0];
  for (let t = 1; t <= F; t++) offset[t] = offset[t - 1] + sizes[t - 1];
  const gi = (t, i) => offset[t] + i;
  const decode = (g) => {
    let t = 0;
    while (t + 1 < F && g >= offset[t + 1]) t++;
    return { t, i: g - offset[t] };
  };

  const d2 = (a, b) => {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    return dx * dx + dy * dy;
  };

  // 邻接位掩码：near1[t][i] 为帧 t 斑点 i 在帧 t+1 内可达的女儿掩码；
  // near2[t][i] 为跨一帧漏检后在帧 t+2 内可达的女儿掩码。
  const D2 = maxDist * maxDist;
  const G2 = 4 * D2;
  const near1 = [];
  const near2 = [];
  for (let t = 0; t < F - 1; t++) {
    near1[t] = frames[t].map((s) => {
      let mask = 0;
      frames[t + 1].forEach((q, j) => { if (d2(s, q) <= D2) mask |= 1 << j; });
      return mask;
    });
    if (t < F - 2) {
      near2[t] = frames[t].map((s) => {
        let mask = 0;
        frames[t + 2].forEach((q, j) => { if (d2(s, q) <= G2) mask |= 1 << j; });
        return mask;
      });
    }
  }

  // 各帧「掩码 → 亮度和」预计算（亮度为用户整数，用 Float64 保证大整数精确）
  const maskBright = frames.map((fr) => {
    const arr = new Float64Array(1 << fr.length);
    for (let m = 1; m < arr.length; m++) {
      const lsb = m & -m;
      arr[m] = arr[m ^ lsb] + fr[CTZ[lsb]].b;
    }
    return arr;
  });

  // 每对 (t, 母本) 的保持单女儿掩码列表与分裂双女儿掩码列表
  const keepOpts = [];
  const splitOpts = [];
  for (let t = 0; t < F - 1; t++) {
    keepOpts[t] = near1[t].map((mask) => BITS[mask].map((j) => 1 << j));
    splitOpts[t] = near1[t].map((mask) => {
      const js = BITS[mask];
      const out = [];
      for (let a = 0; a < js.length; a++) {
        for (let b = a + 1; b < js.length; b++) out.push((1 << js[a]) | (1 << js[b]));
      }
      return out;
    });
  }

  /**
   * 阶段一：边界 t 的联合转移对集合。
   * @returns {Uint32Array} 每个元素为 女儿掩码*512 + 新开漏检母本掩码。
   *   只含「女儿互不共用、漏检不在末帧开」的可行对；同一对只出现一次。
   *   结果只依赖 (t, live, gaps)，与漏检剩余额度无关，全程只算一次。
   */
  const transMemo = new Map();
  function transitions(t, live, gaps) {
    const key = (t << 20) | (live << 10) | gaps;
    const cached = transMemo.get(key);
    if (cached !== undefined) return cached;

    const liveMoms = BITS[live];
    const gapMoms = BITS[gaps];
    const totalTracks = gapMoms.length + liveMoms.length;
    let dp = new Set([0]); // 状态键低 9 位为新开漏检掩码，余为女儿掩码

    // 1) 待补获漏检母本（帧 t-1）：恰一个跨帧女儿
    let processed = 0;
    for (const mi of gapMoms) {
      const cap = near2[t - 1][mi];
      const rest = totalTracks - processed - 1; // 后续母本至少再贡献 1 支
      const ndp = new Set();
      for (const state of dp) {
        const used = state >>> 9;
        for (const b of BITS[cap & ~used]) {
          const used2 = used | (1 << b);
          if (PC[used2] + rest > target) continue;
          ndp.add(used2 * 512);
        }
      }
      dp = ndp;
      processed++;
    }

    // 2) 存活母本（帧 t）：保持一女 / 分裂两女 / 本帧漏检
    const canOpen = t + 2 <= F - 1;
    for (const mi of liveMoms) {
      const miBit = 1 << mi;
      const rest = totalTracks - processed - 1;
      const ndp = new Set();
      for (const state of dp) {
        const used = state >>> 9;
        const opened = state & 511;

        // 2a) 保持
        for (const bit of keepOpts[t][mi]) {
          if (used & bit) continue;
          const used2 = used | bit;
          if (PC[used2] + PC[opened] + rest > target) continue;
          ndp.add(used2 * 512 + opened);
        }
        // 2b) 分裂
        for (const pair of splitOpts[t][mi]) {
          if (used & pair) continue;
          const used2 = used | pair;
          if (PC[used2] + PC[opened] + rest > target) continue;
          ndp.add(used2 * 512 + opened);
        }
        // 2c) 本帧漏检（下一帧必须补获；跨帧无可达斑点则漏检无意义）
        if (canOpen && near2[t][mi] !== 0) {
          const opened2 = opened | miBit;
          if (PC[used] + PC[opened2] + rest <= target) {
            ndp.add(used * 512 + opened2);
          }
        }
      }
      dp = ndp;
      processed++;
    }

    const arr = Uint32Array.from(dp);
    transMemo.set(key, arr);
    return arr;
  }

  // 乐观亮度上界：给定边界 t 的存活支数 co 与待补获漏检支数 cg，未来各帧
  // 至多能安排多少支（每步存活母本至多翻倍、漏检补获单传；开漏检只会让
  // 包络更小，故忽略），且任何时刻支数都不可能超过终帧目标 target，取该
  // 帧最亮的 k 个斑点亮度之和作为松弛上界。
  const topSum = frames.map((fr) => {
    const sorted = fr.map((s) => s.b).sort((a, b) => b - a);
    const pref = new Float64Array(fr.length + 1);
    for (let k = 0; k < fr.length; k++) pref[k + 1] = pref[k] + sorted[k];
    return pref;
  });
  const ubCache = Array.from({ length: F }, () => new Float64Array(81).fill(-1));
  function upperBound(t, co, cg) {
    const slot = co * 9 + cg;
    const cached = ubCache[t][slot];
    if (cached >= 0) return cached;
    let sum = 0;
    let a = co;
    let g = cg;
    for (let s = t + 1; s < F; s++) {
      const k = Math.min(sizes[s], target, 2 * a + g);
      sum += topSum[s][k];
      a = k;
      g = 0;
    }
    ubCache[t][slot] = sum;
    return sum;
  }

  // 计数增长走廊（必要性）：从 (live, gaps) 起末帧存活数能否达到目标。
  function canReachTarget(t, live, gaps) {
    let co = PC[live];
    let cg = PC[gaps];
    for (let s = 1; s <= F - 1 - t; s++) {
      co = Math.min(sizes[t + s], 2 * co + cg);
      cg = 0;
    }
    return co >= target;
  }

  // ---------- 阶段二：取值 DP（最高亮度，其次最少漏检） ----------
  // 备忘值：null = 不可行；否则 {b: 最高后缀亮度, s: 最少后缀漏检}。
  const vMemo = new Map();
  const vKey = (t, live, gaps, left) => ((((t * 256 + live) * 256 + gaps) * 8) + left) * 31;

  function valueDP(t, live, gaps, left) {
    const key = vKey(t, live, gaps, left);
    const cached = vMemo.get(key);
    if (cached !== undefined) return cached;

    const count = PC[live] + PC[gaps];
    if (count > target || left < 0) { vMemo.set(key, null); return null; }
    if (t === F - 1) {
      const v = gaps === 0 && PC[live] === target ? { b: 0, s: 0 } : null;
      vMemo.set(key, v);
      return v;
    }
    if (!canReachTarget(t, live, gaps)) { vMemo.set(key, null); return null; }

    const trans = transitions(t, live, gaps);
    let bestB = -1;
    let bestS = 0;
    // 两遍枚举：先处理零新开漏检的转移，尽快以最大亮度确立现任解；
    // 随后所有「本帧必须开漏检」的分支几乎都能被漏检前缀剪枝。
    for (let pass = 0; pass < 2; pass++) {
      for (const state of trans) {
        const used = state >>> 9;
        const opened = state & 511;
        const openCount = PC[opened];
        if (pass === 0 ? openCount !== 0 : openCount === 0) continue;
        if (openCount > left) continue;
        const h = maskBright[t + 1][used];
        // 乐观上界低于现任亮度必败；持平时若本帧新开漏检已不少于现任
        // 总漏检，也不可能在「亮度持平、漏检更少」上翻盘。剪枝精确。
        if (bestB >= 0) {
          const bound = h + upperBound(t + 1, PC[used], openCount);
          if (bound < bestB) continue;
          if (bound === bestB && openCount > bestS) continue;
        }
        const sub = valueDP(t + 1, used, opened, left - openCount);
        if (sub === null) continue;
        const b = h + sub.b;
        const s = openCount + sub.s;
        if (b > bestB || (b === bestB && s < bestS)) { bestB = b; bestS = s; }
      }
    }
    const v = bestB < 0 ? null : { b: bestB, s: bestS };
    vMemo.set(key, v);
    return v;
  }

  const rootMask = 1 << startIndex;
  const rootV = valueDP(0, rootMask, 0, maxSkip);

  if (rootV === null) {
    // 最早断开帧间：逐步前向展开可达状态，以局部必要存活条件（计数走廊、
    // 漏检必须在补获帧有可达斑点、末帧计数恰为目标）筛选，找出首个
    // 所有后继都无法存活的帧间。
    const viable = (t, live, gaps, left) => {
      if (left < 0) return false;
      if (PC[live] + PC[gaps] > target) return false;
      if (t === F - 1) return gaps === 0 && PC[live] === target;
      if (!canReachTarget(t, live, gaps)) return false;
      if (t >= 1) {
        for (const mi of BITS[gaps]) {
          if (near2[t - 1][mi] === 0) return false;
        }
      }
      return true;
    };

    let reach = new Map();
    if (viable(0, rootMask, 0, maxSkip)) {
      reach.set((rootMask << 10) | 0, { live: rootMask, gaps: 0, left: maxSkip });
    }
    let earliest = 0;
    for (let t = 0; t < F - 1; t++) {
      const next = new Map();
      for (const st of reach.values()) {
        for (const state of transitions(t, st.live, st.gaps)) {
          const used = state >>> 9;
          const opened = state & 511;
          const nleft = st.left - PC[opened];
          if (!viable(t + 1, used, opened, nleft)) continue;
          const k = (used << 10) | opened;
          const prev = next.get(k);
          // 同一 (live, gaps) 保留最大剩余额度：额度越大越可能存活
          if (!prev || nleft > prev.left) next.set(k, { live: used, gaps: opened, left: nleft });
        }
      }
      if (next.size === 0) {
        earliest = t;
        break;
      }
      earliest = t + 1;
      reach = next;
    }
    earliest = Math.min(earliest, F - 2);
    return { feasible: false, earliestBreak: { from: earliest, to: earliest + 1 } };
  }

  const optB = rootV.b;
  const optS = rootV.s;

  // ---------- 阶段三a：精确残量可行性（后缀能否恰好取得亮度 B、漏检 S） ----------
  // 外层键 = (t, live, gaps, left, S)，内层 Map 为残量亮度 B → 0/1。
  const fMemo = new Map();
  const fHeadKey = (t, live, gaps, left, S) =>
    (((((t * 256 + live) * 256 + gaps) * 8 + left) * 8) + S) * 37;

  function feas(t, live, gaps, left, B, S) {
    if (B < 0 || S < 0 || left < 0 || S > left) return 0;
    const count = PC[live] + PC[gaps];
    if (count > target) return 0;
    if (t === F - 1) {
      return (gaps === 0 && PC[live] === target && B === 0 && S === 0) ? 1 : 0;
    }
    if (!canReachTarget(t, live, gaps)) return 0;
    if (B > upperBound(t, PC[live], PC[gaps])) return 0;
    // 取值上界（含额度）：残量 (B,S) 不能被该状态的最优 (亮度, 漏检) 支配。
    const v = valueDP(t, live, gaps, left);
    if (v === null) return 0;
    if (B > v.b) return 0;
    if (B === v.b && S < v.s) return 0;

    const hk = fHeadKey(t, live, gaps, left, S);
    let bucket = fMemo.get(hk);
    if (bucket !== undefined && bucket.has(B)) return bucket.get(B);

    let found = 0;
    for (const state of transitions(t, live, gaps)) {
      const used = state >>> 9;
      const opened = state & 511;
      const openCount = PC[opened];
      if (openCount > left) continue;
      const h = maskBright[t + 1][used];
      if (h > B) continue;
      if (feas(t + 1, used, opened, left - openCount, B - h, S - openCount)) {
        found = 1;
        break;
      }
    }
    if (bucket === undefined) {
      bucket = new Map();
      fMemo.set(hk, bucket);
    }
    bucket.set(B, found);
    return found;
  }

  // ---------- 阶段三b：按裁决顺序排序单个边界的转移对 ----------
  // 先比本帧采用斑点序号元组（由女儿掩码决定）；仅当女儿掩码相同（漏检
  // 掩码不同）时才需要按需计算母本配对再逐女儿比母本全局序号。
  const canonicalMemo = new Map();
  function canonicalMoms(t, live, gaps, usedTarget, openedTarget) {
    const key = ((((t * 256 + live) * 256 + gaps) * 512 + usedTarget) * 512) + openedTarget;
    const cached = canonicalMemo.get(key);
    if (cached !== undefined) return cached;

    // 受限于单一 (usedTarget, openedTarget) 的小型 DP：处理顺序与阶段一
    // 相同（漏检母本先补获，再存活母本），存活母本按 openedTarget 决定
    // 是否本帧漏检；其余必须在 usedTarget 内保持/分裂，铺满 usedTarget。
    // 中间状态只需「已占用女儿掩码」，同掩码保留字典序最小母本向量。
    let dp = new Map([[0, 0]]);

    for (const mi of BITS[gaps]) {
      const gm = gi(t - 1, mi);
      const cap = near2[t - 1][mi] & usedTarget;
      const ndp = new Map();
      for (const [assigned, packed] of dp) {
        for (const b of BITS[cap & ~assigned]) {
          const a2 = assigned | (1 << b);
          const p2 = momSet(packed, b, gm);
          const old = ndp.get(a2);
          ndp.set(a2, old === undefined ? p2 : preferPacked(p2, old, a2));
        }
      }
      dp = ndp;
    }

    for (const mi of BITS[live]) {
      const gm = gi(t, mi);
      const ndp = new Map();
      if (openedTarget & (1 << mi)) {
        // 该母本本帧漏检：不占女儿
        for (const [assigned, packed] of dp) ndp.set(assigned, packed);
      } else {
        for (const [assigned, packed] of dp) {
          // 保持一女
          for (const bit of keepOpts[t][mi]) {
            if (!(usedTarget & bit) || (assigned & bit)) continue;
            const a2 = assigned | bit;
            const p2 = momSet(packed, CTZ[bit], gm);
            const old = ndp.get(a2);
            ndp.set(a2, old === undefined ? p2 : preferPacked(p2, old, a2));
          }
          // 分裂两女
          for (const pair of splitOpts[t][mi]) {
            if ((pair & usedTarget) !== pair || (assigned & pair)) continue;
            const a2 = assigned | pair;
            let p2 = packed;
            for (const b of BITS[pair]) p2 = momSet(p2, b, gm);
            const old = ndp.get(a2);
            ndp.set(a2, old === undefined ? p2 : preferPacked(p2, old, a2));
          }
        }
      }
      dp = ndp;
    }

    const result = dp.has(usedTarget) ? dp.get(usedTarget) : null;
    canonicalMemo.set(key, result);
    return result;
  }

  // 已占用掩码相同的两条母本向量择优：按女儿序号找到首个不同位置，
  // 保留母本全局序号更小者（输入顺序稳定裁决）。
  function preferPacked(cand, old, assignedMask) {
    for (const j of BITS[assignedMask]) {
      const a = momGet(cand, j);
      const b = momGet(old, j);
      if (a !== b) return a < b ? cand : old;
    }
    return old;
  }

  // 同一 (t, live, gaps) 下，按帧裁决签名升序排列转移对（仅最优路径上
  // 的少量状态会调用）。
  function orderedTransitions(t, live, gaps) {
    const list = [...transitions(t, live, gaps)];
    list.sort((sa, sb) => {
      const ua = sa >>> 9;
      const ub = sb >>> 9;
      const c = compareList(BITS[ua], BITS[ub]);
      if (c !== 0) return c;
      // 女儿掩码相同：比较（唯一的）字典序最小母本配对
      const ma = canonicalMoms(t, live, gaps, ua, sa & 511);
      const mb = canonicalMoms(t, live, gaps, ub, sb & 511);
      if (ma === null || mb === null) return 0;
      for (const j of BITS[ua]) {
        const a = momGet(ma, j);
        const b = momGet(mb, j);
        if (a !== b) return a < b ? -1 : 1;
      }
      return 0;
    });
    return list;
  }

  // ---------- 阶段三c：沿最优 (亮度, 漏检) 重建字典序最小链 ----------
  const picks = []; // 每帧 {used, opened, moms(压缩母本向量)}
  let rt = 0;
  let rLive = rootMask;
  let rGaps = 0;
  let rLeft = maxSkip;
  let rB = optB;
  let rS = optS;
  while (rt < F - 1) {
    let chosen = null;
    for (const state of orderedTransitions(rt, rLive, rGaps)) {
      const used = state >>> 9;
      const opened = state & 511;
      const openCount = PC[opened];
      if (openCount > rLeft) continue;
      const h = maskBright[rt + 1][used];
      if (h > rB || openCount > rS) continue;
      if (feas(rt + 1, used, opened, rLeft - openCount, rB - h, rS - openCount)) {
        chosen = { used, opened, openCount, moms: canonicalMoms(rt, rLive, rGaps, used, opened) };
        break;
      }
    }
    if (!chosen) {
      // 理论不可达：(optB, optS) 已由取值 DP 与 feas 保证可行。
      throw new Error('谱系重建失败（内部错误）');
    }
    picks.push(chosen);
    rB -= maskBright[rt + 1][chosen.used];
    rS -= chosen.openCount;
    rLeft -= chosen.openCount;
    rLive = chosen.used;
    rGaps = chosen.opened;
    rt++;
  }

  // 沿链重建母女边
  const edges = [];
  for (let t = 0; t < picks.length; t++) {
    const { used, moms } = picks[t];
    for (const j of BITS[used]) {
      const g = momGet(moms, j);
      const { t: mf, i: mi } = decode(g);
      const gap = mf === t - 1 ? 2 : 1;
      edges.push({
        from: g,
        to: gi(t + 1, j),
        gap,
        dist: Math.sqrt(d2(frames[mf][mi], frames[t + 1][j])),
      });
    }
  }

  const usedPerFrame = Array.from({ length: F }, () => new Set());
  usedPerFrame[0].add(startIndex);
  for (const e of edges) {
    const { t, i } = decode(e.to);
    usedPerFrame[t].add(i);
  }

  return {
    feasible: true,
    root: gi(0, startIndex),
    totalBrightness: frames[0][startIndex].b + optB,
    skips: optS,
    survivors: target,
    edges,
    usedPerFrame: usedPerFrame.map((s) => [...s].sort((a, b) => a - b)),
    _decode: decode,
    _gi: gi,
  };
}

/**
 * 将基于序号的解翻译成带 id 的 JSON 友好结构（页面与测试共用）。
 */
export function presentSolution(spec, result) {
  if (!result.feasible) {
    return {
      feasible: false,
      earliestBreak: result.earliestBreak,
      earliestBreakLabel:
        `第 ${result.earliestBreak.from + 1} 帧 → 第 ${result.earliestBreak.to + 1} 帧`,
    };
  }
  const { frames } = spec;
  const dec = result._decode;
  const childrenOf = new Map();
  const edges = result.edges.map((e) => {
    const mf = dec(e.from);
    const cf = dec(e.to);
    if (!childrenOf.has(e.from)) childrenOf.set(e.from, []);
    childrenOf.get(e.from).push(e.to);
    return {
      fromFrame: mf.t,
      fromId: frames[mf.t][mf.i].id,
      toFrame: cf.t,
      toId: frames[cf.t][cf.i].id,
      gap: e.gap,
      dist: Math.round(e.dist * 100) / 100,
    };
  });
  edges.sort((a, b) =>
    a.fromFrame - b.fromFrame ||
    a.toFrame - b.toFrame ||
    String(a.fromId).localeCompare(String(b.fromId)) ||
    String(a.toId).localeCompare(String(b.toId)));

  let divisions = 0;
  for (const list of childrenOf.values()) if (list.length === 2) divisions++;

  return {
    feasible: true,
    totalBrightness: result.totalBrightness,
    skips: result.skips,
    survivors: result.survivors,
    divisions,
    counts: result.usedPerFrame.map((s) => s.length),
    used: result.usedPerFrame.map((list, t) => list.map((i) => frames[t][i].id)),
    edges,
  };
}
