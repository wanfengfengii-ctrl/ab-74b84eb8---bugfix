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
// 两阶段位掩码 DP：
//  1) 可行性：状态 (帧, 存活掩码, 漏检掩码, 剩余额度) 布尔备忘；后继惰性
//     递归枚举，命中任意可行后继立即返回，并带末帧可达（并集/容量匹配）与
//     逐部分的增量末帧并集剪枝。
//  2) 最优值：自顶向下备忘递归；每个状态的边界内用「逐母本安排女儿」的
//     分层合并内层 DP（复用两个定长 Float64Array 暂存，键=女儿掩码<<9|
//     新开漏检掩码，值为压入单个 48 位整数的最小母本向量）。先用可行性
//     见证与一次价值引导浅回溯取得强 incumbent，再在每层用逐支/全局走廊
//     乐观亮度上界与采用集合排名下界剪枝，终局对幸存分配做精确后缀比较。
//  裁决签名（逐帧采用斑点排名 + 母本向量）压成定宽 BigInt，整数序即
//  既有 compareTuple 的字典序。
//
// 满规模满邻接（7 帧 ×8 斑点、最大位移 1、漏检 5、末帧存活 8）在亚秒级、
// 约百 MB 内精确求解；不缩小输入规模、不删减连接、不取近似、不改裁决语义。

'use strict';

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

  const bits = (m) => {
    const out = [];
    for (let i = 0; m; i++, m >>>= 1) if (m & 1) out.push(i);
    return out;
  };
  // 斑点掩码位序号表 / 0..255 人口数表
  const bitIndex = new Int8Array(256).fill(-1);
  for (let j = 0; j < 8; j++) bitIndex[1 << j] = j;
  const PC = new Uint8Array(256);
  for (let m = 1; m < 256; m++) PC[m] = PC[m & (m - 1)] + 1;

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

  // 各帧「掩码 → 亮度和」预计算
  const maskBright = frames.map((fr) => {
    const arr = new Array(1 << fr.length).fill(0);
    for (let m = 1; m < arr.length; m++) {
      const lsb = m & -m;
      arr[m] = arr[m ^ lsb] + fr[bitIndex[lsb]].b;
    }
    return arr;
  });

  // 每对 (t, 母本) 的保持单女儿掩码列表与分裂双女儿掩码列表
  const keepOpts = [];
  const splitOpts = [];
  for (let t = 0; t < F - 1; t++) {
    keepOpts[t] = near1[t].map((mask) => bits(mask).map((j) => 1 << j));
    splitOpts[t] = near1[t].map((mask) => {
      const js = bits(mask);
      const out = [];
      for (let a = 0; a < js.length; a++) {
        for (let b = a + 1; b < js.length; b++) out.push((1 << js[a]) | (1 << js[b]));
      }
      return out;
    });
  }

  // 末帧可达掩码：reachEnd[t][i] 为帧 t 斑点 i（允许保持/分裂/恰好一帧漏检）
  // 最终能到达的末帧斑点集合。分裂只取并集即可（必要条件不要求人数翻倍细节）。
  const reachEnd = [];
  reachEnd[F - 1] = frames[F - 1].map((_, i) => 1 << i);
  for (let t = F - 2; t >= 0; t--) {
    reachEnd[t] = frames[t].map((_, i) => {
      let m = 0;
      for (let am = near1[t][i]; am; am &= am - 1) {
        m |= reachEnd[t + 1][bitIndex[am & -am]];
      }
      if (t + 2 <= F - 1) {
        for (let am = near2[t][i]; am; am &= am - 1) {
          m |= reachEnd[t + 2][bitIndex[am & -am]];
        }
      }
      return m;
    });
  }
  // 待补获漏检母本（帧 t-1 的 mi）在帧 t+1 补获后，可到末帧的斑点并集。
  const gapReachEndAt = (t, mi) => {
    let m = 0;
    for (let am = near2[t - 1][mi]; am; am &= am - 1) {
      m |= reachEnd[t + 1][bitIndex[am & -am]];
    }
    return m;
  };
  // openEnd[t][mi]：帧 t 母本开漏检（帧 t+2 补获）后可到末帧的斑点并集。
  const openEnd = [];
  for (let t = 0; t < F - 2; t++) {
    openEnd[t] = frames[t].map((_, mi) => {
      let m = 0;
      for (let am = near2[t][mi]; am; am &= am - 1) {
        m |= reachEnd[t + 2][bitIndex[am & -am]];
      }
      return m;
    });
  }

  // 母本向量压缩：nChild(≤8) 个槽、每槽 6 比特；槽 j 位于 6*(nChild-1-j)
  // 位（槽 0 在高位），向量字典序与整数大小完全一致。全局母本序号 ≤ 55。
  const POW64 = [1, 64, 4096, 262144, 16777216, 1073741824, 68719476736,
    4398046511104, 281474976710656];
  const SLOT_EMPTY = 63;

  // 采用斑点集合的「输入顺序」排名：按帧内升序序号元组字典序（长度兜底），
  // 与 compareTuple 完全一致。每帧 256 个掩码一次性排名。
  const usedRank = sizes.map((n) => {
    const rank = new Uint16Array(256);
    const masks = [];
    for (let m = 0; m < 1 << n; m++) masks.push(m);
    masks.sort((a, b) => {
      const aa = bits(a);
      const bb = bits(b);
      const k = Math.min(aa.length, bb.length);
      for (let i = 0; i < k; i++) if (aa[i] !== bb[i]) return aa[i] - bb[i];
      return aa.length - bb.length;
    });
    masks.forEach((m, r) => { rank[m] = r; });
    return rank;
  });

  // 乐观亮度潜力（忽略几何可行性与斑点唯一占用，故恒为上界）：
  //  topBright[t][k] = 帧 t 内最亮的 k 个斑点亮度和；
  //  optFuture[t][n]  = 帧 t..末帧从 n 支出发（每步至多翻倍）的亮度和上界。
  const topBright = frames.map((fr) => {
    const sorted = fr.map((s) => s.b).sort((a, b) => b - a);
    const top = new Int32Array(target + 1);
    let s = 0;
    for (let k = 1; k <= target; k++) {
      if (k - 1 < sorted.length) s += sorted[k - 1];
      top[k] = s;
    }
    return top;
  });
  const optFuture = new Array(F);
  for (let t1 = 0; t1 < F; t1++) {
    const arr = new Int32Array(target + 1);
    for (let n0 = 1; n0 <= target; n0++) {
      let n = n0;
      let sum = 0;
      for (let s = t1; s < F; s++) {
        sum += topBright[s][Math.min(n, sizes[s])];
        n = Math.min(target, 2 * n);
      }
      arr[n0] = sum;
    }
    optFuture[t1] = arr;
  }
  // topFree[t1][used*9+k]：帧 t1 未被 used 占用的斑点中最亮 k 个的亮度和。
  const topFree = sizes.map((n, t1) => {
    const arr = new Int32Array(256 * 9);
    for (let used = 0; used < 1 << n; used++) {
      const freeB = [];
      for (let j = 0; j < n; j++) if (!(used & (1 << j))) freeB.push(frames[t1][j].b);
      freeB.sort((a, b) => b - a);
      let s = 0;
      for (let k = 1; k <= 8; k++) {
        if (k - 1 < freeB.length) s += freeB[k - 1];
        arr[used * 9 + k] = s;
      }
    }
    return arr;
  });
  // lbUsedRank[t1][used*9+k]：used 并入序号最小的 k 个空闲斑点后的排名。
  const lbUsedRank = sizes.map((n, t1) => {
    const arr = new Uint16Array(256 * 9);
    const rank = usedRank[t1];
    for (let used = 0; used < 1 << n; used++) {
      const freeJ = [];
      for (let j = 0; j < n; j++) if (!(used & (1 << j))) freeJ.push(j);
      let m = used;
      for (let k = 0; k <= 8; k++) {
        if (k > 0 && k - 1 < freeJ.length) m |= 1 << freeJ[k - 1];
        arr[used * 9 + k] = rank[m];
      }
    }
    return arr;
  });

  // 逐支乐观后缀亮度（忽略斑点唯一占用、目标人数与漏检额度，故恒为上界）：
  //  K[t][j] = 单支从帧 t 斑点 j 出发（含其亮度）到末帧的最大亮度和，
  //    保持取可达女儿最大 K，分裂取可达女儿中最大两个 K 之和，漏检取
  //    跨帧可达点最大 K。
  //  S[t][j] = K[t][j] - b_t[j]：该点以后（不含自身）的乐观后缀。
  //  keepRem[t][j] = K[t][j] - b_t[j]：帧 t 母本对帧 t+1 起的乐观贡献。
  //  gapRem[t][j] ：帧 t 母本漏检一帧后（帧 t+2 起）的乐观贡献。
  const K = [];
  const S = [];
  const NEG = -1e15;
  for (let t = F - 1; t >= 0; t--) {
    const kt = new Float64Array(sizes[t]);
    for (let j = 0; j < sizes[t]; j++) {
      if (t === F - 1) {
        kt[j] = frames[t][j].b;
        continue;
      }
      let m1 = NEG;
      let m2 = NEG;
      for (let am = near1[t][j]; am; am &= am - 1) {
        const v = K[t + 1][bitIndex[am & -am]];
        if (v > m1) { m2 = m1; m1 = v; } else if (v > m2) m2 = v;
      }
      const split = m2 === NEG ? NEG : m1 + m2;
      let skip = NEG;
      if (t + 2 <= F - 1) {
        for (let am = near2[t][j]; am; am &= am - 1) {
          const v = K[t + 2][bitIndex[am & -am]];
          if (v > skip) skip = v;
        }
      }
      kt[j] = frames[t][j].b + Math.max(m1, split, skip);
    }
    K[t] = kt;
    S[t] = Float64Array.from(kt, (v, j) => v - frames[t][j].b);
  }
  // 待补获漏检母本（位于帧 t-1）在帧 t+1 的乐观贡献 = 跨帧可达最大 K。
  const gapRemAt = (t, mi) => {
    let v = NEG;
    for (let am = near2[t - 1][mi]; am; am &= am - 1) {
      const k = bitIndex[am & -am];
      if (K[t + 1][k] > v) v = K[t + 1][k];
    }
    return v;
  };

  // 静态「按乐观贡献降序」候选表（与具体状态无关，仅女儿碰撞需在搜索时过滤）：
  //  gapActs[t][mi]  边界 t 的漏检母本（帧 t-1）→ 帧 t+1 女儿位，按 K 降序；
  //  liveActs[t][mi] 存活母本（帧 t）的 保持/分裂 动作（动作位掩码），
  //                  按「涉及女儿的 K 之和」降序（并列时保持优先，再按掩码）；
  //  openGain[t][mi] 该母本开漏检后在帧 t+2 补获的乐观贡献（NEG 表示不能）。
  const gapActs = [];
  const liveActs = [];
  const openGain = [];
  for (let t = 0; t < F - 1; t++) {
    gapActs[t] = t >= 1
      ? frames[t - 1].map((_, mi) =>
          bits(near2[t - 1][mi])
            .map((j) => ({ bit: 1 << j, g: K[t + 1][j] }))
            .sort((a, b) => b.g - a.g || a.bit - b.bit))
      : null;
    liveActs[t] = frames[t].map((_, mi) => {
      const acts = [];
      for (const bit of keepOpts[t][mi]) acts.push({ mask: bit, g: K[t + 1][bitIndex[bit]] });
      for (const pair of splitOpts[t][mi]) {
        let g = 0;
        for (let bm = pair; bm; bm &= bm - 1) g += K[t + 1][bitIndex[bm & -bm]];
        acts.push({ mask: pair, g });
      }
      // 分裂两女之和通常 ≥ 保持；并列时保持（单比特）优先以贴近规范序。
      acts.sort((a, b) => b.g - a.g || (a.mask & (a.mask - 1)) - (b.mask & (b.mask - 1)) || a.mask - b.mask);
      return acts;
    });
    openGain[t] = frames[t].map((_, mi) => {
      if (t + 2 > F - 1) return NEG;
      let v = NEG;
      for (let am = near2[t][mi]; am; am &= am - 1) {
        const k = K[t + 2][bitIndex[am & -am]];
        if (k > v) v = k;
      }
      return v;
    });
  }

  // 状态索引 ((live*256+gaps)*8+left)。
  const NSTATES = 256 * 256 * 8;
  const stIdx = (live, gaps, left) => live * 2048 + gaps * 8 + left;
  const stateKey = (t, live, gaps, left) => t * NSTATES + stIdx(live, gaps, left);

  // 计数增长走廊：从 (live, gaps) 起，每步至多翻倍，漏检补获只能单传，
  // 判断末帧存活数能否达到目标。
  function canReachTarget(t, live, gaps) {
    let co = PC[live];
    let cg = PC[gaps];
    for (let s = 1; s <= F - 1 - t; s++) {
      co = Math.min(sizes[t + s], 2 * co + cg);
      cg = 0;
    }
    return co >= target;
  }

  // 局部必要条件：额度、掩码合法、计数走廊、漏检须有可达补获点、末帧恰目标。
  function stateLocallyOk(t, live, gaps, left) {
    if (left < 0) return false;
    if (live >= 1 << sizes[t]) return false;
    if (t === 0 ? gaps !== 0 : gaps >= 1 << sizes[t - 1]) return false;
    if (PC[live] + PC[gaps] > target) return false;
    if (t === F - 1) return gaps === 0 && PC[live] === target;
    if (!canReachTarget(t, live, gaps)) return false;
    if (t >= 1) {
      for (let m = gaps; m; m &= m - 1) {
        if (near2[t - 1][bitIndex[m & -m]] === 0) return false;
      }
    }
    return true;
  }

  /**
   * 惰性枚举边界 t 的全部完整分配 (used, opened)，按输入顺序生成：
   * 先处理待补获漏检母本（帧 t-1，序号升序），再处理存活母本（帧 t，
   * 序号升序）；每个存活母本依次尝试 保持（女儿升序）→ 分裂（配对升序）
   * → 本帧漏检。fn 返回 true 即提前终止。
   * useEndPrune=true 时增量维护「所有支可达末帧斑点并集」，并集不足 target
   * 个不同末帧斑点的部分分配立即剪去（可行性阶段专用；最早断帧的前向枚举
   * 须保留局部可行状态，故不传此标志）。
   */
  function forEachAssign(t, live, gaps, fn, useEndPrune = false) {
    const gapMoms = bits(gaps);
    const liveMoms = bits(live);
    const total = gapMoms.length + liveMoms.length;
    const canOpen = t + 2 <= F - 1;
    let stop = false;

    // 每位母亲（处理顺序）在尚未安排时的末帧可达潜力，并求后缀并集。
    let suffixEnd = null;
    if (useEndPrune) {
      suffixEnd = new Int32Array(total + 1);
      for (let p = total - 1; p >= 0; p--) {
        let pot;
        if (p < gapMoms.length) pot = gapReachEndAt(t, gapMoms[p]);
        else pot = reachEnd[t][liveMoms[p - gapMoms.length]];
        suffixEnd[p] = suffixEnd[p + 1] | pot;
      }
    }

    function rec(k, used, opened, endM) {
      if (stop) return;
      const rest = total - k; // 尚未处理的母本，每支至少再占 1 个存活名额
      if (PC[used] + PC[opened] + rest > target) return;
      if (useEndPrune && PC[endM | suffixEnd[k]] < target) return;
      if (k === total) {
        if (fn(used, opened)) stop = true;
        return;
      }
      if (k < gapMoms.length) {
        const mi = gapMoms[k];
        const cap = near2[t - 1][mi];
        for (let avail = cap & ~used; avail; avail &= avail - 1) {
          const bit = avail & -avail;
          const e = useEndPrune ? endM | reachEnd[t + 1][bitIndex[bit]] : 0;
          rec(k + 1, used | bit, opened, e);
          if (stop) return;
        }
        return;
      }
      const li = k - gapMoms.length;
      const mi = liveMoms[li];
      for (const bit of keepOpts[t][mi]) {
        if (used & bit) continue;
        const e = useEndPrune ? endM | reachEnd[t + 1][bitIndex[bit]] : 0;
        rec(k + 1, used | bit, opened, e);
        if (stop) return;
      }
      for (const pair of splitOpts[t][mi]) {
        if (used & pair) continue;
        let e = 0;
        if (useEndPrune) {
          const a = bitIndex[pair & -pair];
          const b = bitIndex[pair & (pair - 1)];
          e = endM | reachEnd[t + 1][a] | reachEnd[t + 1][b];
        }
        rec(k + 1, used | pair, opened, e);
        if (stop) return;
      }
      if (canOpen && near2[t][mi] !== 0) {
        const e = useEndPrune ? endM | openEnd[t][mi] : 0;
        rec(k + 1, used, opened | (1 << mi), e);
        if (stop) return;
      }
    }
    rec(0, 0, 0, 0);
  }

  // 末帧必要条件：容量二分匹配（复用静态暂存，零分配）。每条当前支未来
  // 最多翻倍 2^(剩余边界) 次；展开为至多 target 个同可达掩码的容量副本，
  // 用 Kuhn 匹配检查能否向末帧注入 target 个互异斑点。
  const matchScratch = new Int8Array(8);
  const seenScratch = new Uint8Array(8);
  const endFeasCopies = [];
  function endFeasible(t, live, gaps) {
    const rem = F - 1 - t;
    const grow = rem >= 3 ? target : Math.min(target, 1 << rem);
    // 构造容量副本掩码序列（≤8 支 × target）。
    const copies = endFeasCopies;
    copies.length = 0;
    const addTrack = (r) => {
      if (r === 0) return false;
      const cap = Math.min(grow, PC[r], target);
      for (let c = 0; c < cap; c++) copies.push(r);
      return true;
    };
    for (let m = live; m; m &= m - 1) {
      if (!addTrack(reachEnd[t][bitIndex[m & -m]])) return false;
    }
    for (let m = gaps; m; m &= m - 1) {
      if (!addTrack(gapReachEndAt(t, bitIndex[m & -m]))) return false;
    }
    const nLast = sizes[F - 1];
    matchScratch.fill(-1, 0, nLast);
    let matched = 0;
    for (let i = 0; i < copies.length && matched < target; i++) {
      seenScratch.fill(0, 0, nLast);
      const aug = (ti) => {
        for (let e = copies[ti]; e; e &= e - 1) {
          const j = bitIndex[e & -e];
          if (seenScratch[j]) continue;
          seenScratch[j] = 1;
          if (matchScratch[j] < 0 || aug(matchScratch[j])) {
            matchScratch[j] = ti;
            return true;
          }
        }
        return false;
      };
      if (aug(i)) matched++;
    }
    return matched >= target;
  }

  // ---- 第一阶段：布尔可行性备忘（惰性递归，命中即停） ----
  // 可行状态同时记录「首个找到的可行后继」(used<<9|opened)，供最优阶段取一个
  // 保证可行的初始 incumbent（剪枝界）。
  const feasMemo = new Map();
  const witness = new Map();
  const baseKey = (t, live, gaps) => (t * 256 + live) * 256 + gaps;

  // 结构（与额度无关）必要条件；false 时该状态任意额度都不可行。
  function structurallyPossible(t, live, gaps) {
    if (live >= 1 << sizes[t]) return false;
    if (t === 0 ? gaps !== 0 : gaps >= 1 << sizes[t - 1]) return false;
    if (PC[live] + PC[gaps] > target) return false;
    if (!canReachTarget(t, live, gaps)) return false;
    if (t >= 1) {
      for (let m = gaps; m; m &= m - 1) {
        if (near2[t - 1][bitIndex[m & -m]] === 0) return false;
      }
    }
    if (!endFeasible(t, live, gaps)) return false;
    return true;
  }

  function feasible(t, live, gaps, left) {
    const key = stateKey(t, live, gaps, left);
    const known = feasMemo.get(key);
    if (known !== undefined) return known;
    if (left < 0 || !structurallyPossible(t, live, gaps)) {
      feasMemo.set(key, false);
      return false;
    }
    if (t === F - 1) {
      const ok = gaps === 0 && PC[live] === target;
      feasMemo.set(key, ok);
      return ok;
    }
    let ok = false;
    forEachAssign(t, live, gaps, (used, opened) => {
      const oc = PC[opened];
      if (oc > left) return false;
      if (feasible(t + 1, used, opened, left - oc)) {
        witness.set(baseKey(t, live, gaps), (used << 9) | opened);
        ok = true;
        return true;
      }
      return false;
    }, true);
    feasMemo.set(key, ok);
    return ok;
  }

  const rootMask = 1 << startIndex;
  if (!feasible(0, rootMask, 0, maxSkip)) {
    // 最早断开帧间：逐步前向展开可达状态，以局部必要存活条件筛选，
    // 找出首个所有后继都无法存活的帧间。
    let reach = new Map();
    if (stateLocallyOk(0, rootMask, 0, maxSkip)) {
      reach.set(stIdx(rootMask, 0, maxSkip), { live: rootMask, gaps: 0, left: maxSkip });
    }
    let earliest = 0;
    for (let t = 0; t < F - 1; t++) {
      const next = new Map();
      for (const st of reach.values()) {
        forEachAssign(t, st.live, st.gaps, (used, opened) => {
          const nleft = st.left - PC[opened];
          if (!stateLocallyOk(t + 1, used, opened, nleft)) return false;
          const k = stIdx(used, opened, nleft);
          if (!next.has(k)) next.set(k, { live: used, gaps: opened, left: nleft });
          return false;
        });
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

  // ---- 第二阶段：最优值备忘递归（边界内分支定界 DFS） ----
  // 值编码 code = 亮度*16 - 漏检数（漏检 ≤5），数值序 = (亮度↑, 漏检↓)。
  // 签名：每帧 57 比特块 = 采用集合排名（9 位）| 母本向量（48 位），
  // 早帧在高位、晚帧在低位拼接；整数序即「逐帧 (采用序号, 母本序号)」
  // 字典序，末帧空块 0。
  const SIGVEC_BITS = 48n;
  const SIGFRAME_BITS = 57n;
  const frameShift = new Array(F).fill(0n);
  for (let t = F - 2; t >= 0; t--) {
    frameShift[t] = frameShift[t + 1] + SIGFRAME_BITS;
  }

  // 分层内层 DP 的稠密暂存（键 = used<<9|opened，17 位），跨所有状态复用；
  // 每层只重置 touched 到过的槽，避免逐状态 Map 分配。NaN 表示空。
  const LEV_SIZE = 1 << 17;
  const levA = new Float64Array(LEV_SIZE).fill(NaN);
  const levB = new Float64Array(LEV_SIZE).fill(NaN);
  const touchedA = [];
  const touchedB = [];

  const memo = new Map();
  function solve(t, live, gaps, left) {
    const key = stateKey(t, live, gaps, left);
    if (memo.has(key)) return memo.get(key);

    if (t === F - 1 || !stateLocallyOk(t, live, gaps, left)) {
      const leaf = (t === F - 1 && stateLocallyOk(t, live, gaps, left))
        ? { code: 0, sig: 0n, state: -1, vec: 0, usedRank: 0 }
        : null;
      memo.set(key, leaf);
      return leaf;
    }
    // 确保可行性已判定（同时记录规范序首个可行后继作为 incumbent 来源）。
    if (!feasible(t, live, gaps, left)) {
      memo.set(key, null);
      return null;
    }

    const nChild = sizes[t + 1];
    const gapMoms = bits(gaps);
    const liveMoms = bits(live);
    const totalTracks = gapMoms.length + liveMoms.length;
    const canOpen = t + 2 <= F - 1;
    const rankT = usedRank[t + 1];
    const shift = frameShift[t + 1];
    const freeTbl = topFree[t + 1];
    const rankTbl = lbUsedRank[t + 1];
    const hasFuture = t + 2 <= F - 1;

    // 既不能安排相邻女儿、又不能开漏检的存活母本 → 整状态不可行。
    for (const mi of liveMoms) {
      if (near1[t][mi] === 0 && !(canOpen && near2[t][mi] !== 0)) {
        memo.set(key, null);
        return null;
      }
    }

    let best = null;

    // 每个存活母本开漏检后，在帧 t+2 被补获的乐观贡献（最大 K[t+2]）。
    const openCap = new Float64Array(8).fill(NEG);
    if (canOpen) {
      for (const mi of liveMoms) {
        let v = NEG;
        for (let am = near2[t][mi]; am; am &= am - 1) {
          const k = bitIndex[am & -am];
          if (K[t + 2][k] > v) v = K[t + 2][k];
        }
        openCap[mi] = v;
      }
    }

    // 处理序列：先待补获漏检母本（帧 t-1），再存活母本（帧 t），均按序号。
    // 每位母亲携带其「未安排时」的独立乐观贡献，供分层 DP 的上界使用。
    const moms = [];
    for (const mi of gapMoms) moms.push({ kind: 0, mi, contrib: gapRemAt(t, mi) });
    for (const mi of liveMoms) moms.push({ kind: 1, mi, contrib: S[t][mi] });
    // remByLevel[k]：处理完前 k 位母亲后，其余母亲的独立乐观贡献之和。
    const remByLevel = new Float64Array(totalTracks + 1);
    // endSuffix[k]：其余母亲（按其最优动作）可到末帧的斑点并集。
    const endSuffix = new Int32Array(totalTracks + 1);
    for (let k = totalTracks - 1; k >= 0; k--) {
      remByLevel[k] = remByLevel[k + 1] + moms[k].contrib;
      const pot = moms[k].kind === 0
        ? gapReachEndAt(t, moms[k].mi)
        : reachEnd[t][moms[k].mi];
      endSuffix[k] = endSuffix[k + 1] | pot;
    }
    // 帧 t+1 女儿集合 → 末帧可达并集；帧 t 开漏检母本集合 → 末帧可达并集。
    const endByUsed = new Int32Array(1 << nChild);
    for (let m = 1; m < 1 << nChild; m++) {
      endByUsed[m] = endByUsed[m ^ (m & -m)] | reachEnd[t + 1][bitIndex[m & -m]];
    }
    const endByOpen = new Int32Array(256);
    if (canOpen) {
      for (let m = 1; m < 256; m++) {
        endByOpen[m] = endByOpen[m ^ (m & -m)] | openEnd[t][bitIndex[m & -m]];
      }
    }
    // 掩码 → 乐观贡献和（仅依赖掩码）。
    const sumK1 = new Float64Array(1 << nChild);
    for (let m = 1; m < 1 << nChild; m++) {
      const lsb = m & -m;
      sumK1[m] = sumK1[m ^ lsb] + K[t + 1][bitIndex[lsb]];
    }
    const sumOpen = new Float64Array(256);
    for (let m = 1; m < 256; m++) {
      const lsb = m & -m;
      sumOpen[m] = sumOpen[m ^ lsb] + openCap[bitIndex[lsb]];
    }

    // 分层部分分配键 (used<<9 | opened) 的乐观值上界（恒不小于真实完成值）。
    const boundOf = (k, used, opened) => {
      const gR = Math.max(0, gapMoms.length - k);
      const lR = Math.max(0, liveMoms.length - Math.max(0, k - gapMoms.length));
      const ubTrack = (sumK1[used] + sumOpen[opened] + remByLevel[k]) * 16 - PC[opened];
      // 全局走廊上界：已放置女儿在本帧的真实亮度 + 剩余母亲最亮空闲女儿
      // + 帧 t+2 起按人数翻倍走廊（忽略几何约束）。
      const dMaxG = Math.min(nChild - PC[used], 2 * lR + gR);
      const thisAdd = freeTbl[used * 9 + dMaxG];
      const nAtT2 = 2 * PC[used] + PC[opened] + 2 * gR + 4 * lR;
      const fut = hasFuture
        ? optFuture[t + 2][Math.max(1, Math.min(target, nAtT2))]
        : 0;
      const ubGlobal = (maskBright[t + 1][used] + thisAdd + fut) * 16 - PC[opened];
      return Math.min(ubTrack, ubGlobal);
    };

    // 给定完整分配 (used, opened)，构造字典序最小的压缩母本向量：
    // 逐母亲内层 DP（键为已占用女儿掩码），只允许落在 used 内的女儿；
    // opened 中的存活母本不安排女儿。规模 ≤8 位母亲，代价很小。
    const buildMinVec = (used, opened) => {
      const setSlot = (v, j, gm) => v + (gm - SLOT_EMPTY) * POW64[nChild - 1 - j];
      let d = new Map([[0, POW64[nChild] - 1]]);
      for (const mi of gapMoms) {
        const gm = gi(t - 1, mi);
        const nd = new Map();
        for (const [claim, v] of d) {
          for (let am = near2[t - 1][mi] & used & ~claim; am; am &= am - 1) {
            const bit = am & -am;
            const k2 = claim | bit;
            const v2 = setSlot(v, bitIndex[bit], gm);
            const old = nd.get(k2);
            if (old === undefined || v2 < old) nd.set(k2, v2);
          }
        }
        d = nd;
      }
      for (const mi of liveMoms) {
        if (opened & (1 << mi)) continue;
        const gm = gi(t, mi);
        const nd = new Map();
        const offer = (k2, v2) => {
          const old = nd.get(k2);
          if (old === undefined || v2 < old) nd.set(k2, v2);
        };
        for (const [claim, v] of d) {
          for (const bit of keepOpts[t][mi]) {
            if ((used & bit) !== bit || (claim & bit)) continue;
            offer(claim | bit, setSlot(v, bitIndex[bit], gm));
          }
          for (const pair of splitOpts[t][mi]) {
            if ((used & pair) !== pair || (claim & pair)) continue;
            let v2 = v;
            for (let bm = pair; bm; bm &= bm - 1) v2 = setSlot(v2, bitIndex[bm & -bm], gm);
            offer(claim | pair, v2);
          }
        }
        d = nd;
      }
      return d.has(used) ? d.get(used) : null;
    };

    // ---- 阶段 A：沿可行性阶段记录的「最小额度可行后继」直接取一个
    // 保证可行的完整解作为剪枝 incumbent（其后缀由 solve 正常求出）。 ----
    const wTr = witness.get(baseKey(t, live, gaps));
    if (wTr !== undefined) {
      const wUsed = wTr >> 9;
      const wOpened = wTr & 511;
      const sub = solve(t + 1, wUsed, wOpened, left - PC[wOpened]);
      if (sub) {
        const vec = buildMinVec(wUsed, wOpened);
        if (vec !== null) {
          const code = sub.code + maskBright[t + 1][wUsed] * 16 - PC[wOpened];
          const sig = (BigInt(rankT[wUsed]) << SIGVEC_BITS | BigInt(vec)) << shift | sub.sig;
          best = { code, sig, state: wTr, vec, usedRank: rankT[wUsed] };
        }
      }
    }

    // 价值引导贪心改进：沿静态「乐观贡献降序」候选表做极浅回溯，
    // 零数组分配（候选表与槽位数组复用），最多评估 budget 个完成。
    // incumbent 已达根乐观上界时跳过（价值上不可能更好）。
    const ubRoot = boundOf(0, 0, 0);
    if (!(best && best.code >= ubRoot)) {
      const slotMomG = new Int8Array(nChild).fill(-1);
      let budget = 20;
      let nodes = 64;
      const gGap = gapActs[t];
      const gLive = liveActs[t];
      const gOpen = openGain[t];
      const recG = (k, used, opened) => {
        if (budget <= 0 || nodes <= 0) return;
        nodes--;
        const rest = totalTracks - k;
        if (PC[used] + PC[opened] + rest > target) return;
        if (k === totalTracks) {
          budget--;
          if (PC[opened] > left) return;
          const sub = solve(t + 1, used, opened, left - PC[opened]);
          if (!sub) return;
          // 规范母本向量：空槽位哨兵 SLOT_EMPTY（与阶段 B 同口径）。
          let vec = POW64[nChild] - 1;
          for (let bm = used; bm; bm &= bm - 1) {
            const j = bitIndex[bm & -bm];
            vec += (slotMomG[j] - SLOT_EMPTY) * POW64[nChild - 1 - j];
          }
          const code = sub.code + maskBright[t + 1][used] * 16 - PC[opened];
          const sig = (BigInt(rankT[used]) << SIGVEC_BITS | BigInt(vec)) << shift | sub.sig;
          if (!best || code > best.code || (code === best.code && sig < best.sig)) {
            best = { code, sig, state: (used << 9) | opened, vec, usedRank: rankT[used] };
          }
          return;
        }
        const mom = moms[k];
        const gm = mom.kind === 0 ? gi(t - 1, mom.mi) : gi(t, mom.mi);
        if (mom.kind === 0) {
          for (const a of gGap[mom.mi]) {
            if (budget <= 0 || nodes <= 0) return;
            if (used & a.bit) continue;
            const j = bitIndex[a.bit];
            slotMomG[j] = gm;
            recG(k + 1, used | a.bit, opened);
            slotMomG[j] = -1;
          }
          return;
        }
        const mi = mom.mi;
        for (const a of gLive[mi]) {
          if (budget <= 0 || nodes <= 0) return;
          if (used & a.mask) continue;
          const touched = [];
          for (let bm = a.mask; bm; bm &= bm - 1) {
            const j = bitIndex[bm & -bm];
            slotMomG[j] = gm; touched.push(j);
          }
          recG(k + 1, used | a.mask, opened);
          for (const j of touched) slotMomG[j] = -1;
        }
        if (gOpen[mi] !== NEG && PC[opened] + 1 <= left) {
          recG(k + 1, used, opened | (1 << mi));
        }
      };
      recG(0, 0, 0);
    }

    // ---- 阶段 B：分层合并内层 DP（稠密暂存，键 = used<<9|opened）----
    // 同 (used,opened) 只留字典序最小（即数值最小）的压缩母本向量。
    levA[0] = POW64[nChild] - 1;
    touchedA.length = 0;
    touchedA.push(0);
    let cur = levA, nxt = levB, curT = touchedA, nxtT = touchedB;
    for (let k = 0; k < totalTracks; k++) {
      const mom = moms[k];
      const { kind, mi } = mom;
      const gm = kind === 0 ? gi(t - 1, mi) : gi(t, mi);
      const k1 = k + 1;
      const gR1 = Math.max(0, gapMoms.length - k1);
      const lR1 = Math.max(0, liveMoms.length - Math.max(0, k1 - gapMoms.length));
      const remK1 = remByLevel[k1];
      const dMaxC = 2 * lR1 + gR1;
      const cT2 = 2 * gR1 + 4 * lR1;
      const dRankMin = gR1;
      const restN = totalTracks - k1;
      nxtT.length = 0;
      const put = (key2, vec2) => {
        const old = nxt[key2];
        if (Number.isNaN(old)) { nxt[key2] = vec2; nxtT.push(key2); }
        else if (vec2 < old) nxt[key2] = vec2;
      };
      const endSuffixK1 = endSuffix[k1];
      const liveOk = (u, o) => {
        if (PC[u] + PC[o] + restN > target) return false;
        // 末帧可达斑点并集必须足以容纳 target 个互不重合的存活。
        if (PC[endByUsed[u] | endByOpen[o] | endSuffixK1] < target) return false;
        if (best) {
          const bc = best.code;
          const base = (sumK1[u] + sumOpen[o] + remK1) * 16 - PC[o];
          if (base < bc) return false;
          const dm = Math.min(nChild - PC[u], dMaxC);
          const n2 = Math.max(1, Math.min(target, 2 * PC[u] + PC[o] + cT2));
          const glob = (maskBright[t + 1][u] + freeTbl[u * 9 + dm]
            + (hasFuture ? optFuture[t + 2][n2] : 0)) * 16 - PC[o];
          if (glob < bc) return false;
          // 上界恰等当前最优：采用集合排名乐观下界必须不大于当前最优。
          if (base === bc || glob === bc) {
            const dMax = dm;
            let minRank = 0xffff;
            for (let d = dRankMin; d <= dMax; d++) {
              const r = rankTbl[u * 9 + d];
              if (r < minRank) minRank = r;
            }
            if (minRank > best.usedRank) return false;
          }
        }
        return true;
      };
      for (const state of curT) {
        const vec = cur[state];
        const used = state >> 9;
        const opened = state & 511;
        if (kind === 0) {
          for (let am = near2[t - 1][mi] & ~used; am; am &= am - 1) {
            const bit = am & -am;
            const j = bitIndex[bit];
            const used2 = used | bit;
            if (!liveOk(used2, opened)) continue;
            put((used2 << 9) | opened,
              vec + (gm - SLOT_EMPTY) * POW64[nChild - 1 - j]);
          }
          continue;
        }
        const setSlot = (v, j) => v + (gm - SLOT_EMPTY) * POW64[nChild - 1 - j];
        const mustOpen = near1[t][mi] === 0;
        if (!mustOpen) {
          for (const bit of keepOpts[t][mi]) {
            if (used & bit) continue;
            const used2 = used | bit;
            if (!liveOk(used2, opened)) continue;
            put((used2 << 9) | opened, setSlot(vec, bitIndex[bit]));
          }
          for (const pair of splitOpts[t][mi]) {
            if (used & pair) continue;
            const used2 = used | pair;
            if (!liveOk(used2, opened)) continue;
            let v2 = vec;
            for (let bm = pair; bm; bm &= bm - 1) v2 = setSlot(v2, bitIndex[bm & -bm]);
            put((used2 << 9) | opened, v2);
          }
        }
        if (canOpen && near2[t][mi] !== 0 && PC[opened] + 1 <= left) {
          const opened2 = opened | (1 << mi);
          if (liveOk(used, opened2)) put((used << 9) | opened2, vec);
        }
      }
      for (const s0 of curT) cur[s0] = NaN;
      const tmpA = cur; cur = nxt; nxt = tmpA;
      const tmpT = curT; curT = nxtT; nxtT = tmpT;
    }

    // 终局：先取出幸存的 (键, 向量) 并复位暂存，再递归求后缀值——
    // 递归的 solve 会复用同一组 levA/levB 暂存。
    const finals = [];
    for (const state of curT) finals.push(state, cur[state]);
    for (const s0 of curT) cur[s0] = NaN;
    for (let fi = 0; fi < finals.length; fi += 2) {
      const state = finals[fi];
      const vec = finals[fi + 1];
      const used = state >> 9;
      const opened = state & 511;
      if (PC[opened] > left) continue;
      const sub = solve(t + 1, used, opened, left - PC[opened]);
      if (!sub) continue;
      const code = sub.code + maskBright[t + 1][used] * 16 - PC[opened];
      const sig = (BigInt(rankT[used]) << SIGVEC_BITS | BigInt(vec)) << shift | sub.sig;
      if (!best || code > best.code || (code === best.code && sig < best.sig)) {
        best = { code, sig, state, vec, usedRank: rankT[used] };
      }
    }


    memo.set(key, best);
    return best;
  }

  const root = solve(0, rootMask, 0, maxSkip);
  // 沿最优链重建母女边
  const edges = [];
  let t = 0;
  let live = rootMask;
  let gaps = 0;
  let left = maxSkip;
  let node = root;
  while (node && node.state >= 0) {
    const state = node.state;
    const vec = node.vec;
    const used = state >> 9;
    const opened = state & 511;
    const nChild = sizes[t + 1];
    for (let bm = used; bm; bm &= bm - 1) {
      const j = bitIndex[bm & -bm];
      const g = Math.floor(vec / POW64[nChild - 1 - j]) & 63;
      const { t: mf, i: mi } = decode(g);
      const gap = mf === t - 1 ? 2 : 1;
      edges.push({
        from: g,
        to: gi(t + 1, j),
        gap,
        dist: Math.sqrt(d2(frames[mf][mi], frames[t + 1][j])),
      });
    }
    live = used;
    gaps = opened;
    left -= PC[opened];
    t++;
    node = memo.get(stateKey(t, live, gaps, left));
  }

  const usedPerFrame = Array.from({ length: F }, () => new Set());
  usedPerFrame[0].add(startIndex);
  for (const e of edges) {
    const { t: tt, i } = decode(e.to);
    usedPerFrame[tt].add(i);
  }

  const rootSkips = (16 - (Math.round(root.code) % 16)) % 16;
  const rootBright = (Math.round(root.code) + rootSkips) / 16;

  return {
    feasible: true,
    root: gi(0, startIndex),
    totalBrightness: frames[0][startIndex].b + rootBright,
    skips: rootSkips,
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
    String(b.toId).localeCompare(String(b.toId)));

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
