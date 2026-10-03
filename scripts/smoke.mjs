// 烟测：
//  1) 等待 /healthz 通过（服务由外部提供 BASE_URL，或本脚本临时拉起一个）；
//  2) 抓取页面与脚本资源，确认站点可服务；
//  3) 在同一求解器内核上跑「含一次分裂 + 一次漏检」的谱系场景并校验结果；
// 以退出码报告：0 通过，非 0 失败。
'use strict';

import { spawn } from 'node:child_process';
import { readFileSync, existsSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { dirname, resolve } from 'node:path';
import { normalizeSpec, solveLineage, presentSolution } from '../public/js/lineage.js';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const HOST = process.env.WEB_HOST || 'web';
const PORT = process.env.WEB_PORT || '8080';
const BASE_URL = process.env.BASE_URL ||
  ((HOST === 'web' || HOST === '0.0.0.0') ? `http://web:${PORT}` : `http://127.0.0.1:${PORT}`);

let ownServer = null;
let portFile = null;

function log(msg) { console.log(`[smoke] ${msg}`); }
function fail(msg) { console.error(`[smoke] 失败: ${msg}`); process.exitCode = 1; throw new Error(msg); }

async function waitHealthy(base, timeoutMs = 30000) {
  const deadline = Date.now() + timeoutMs;
  let lastErr;
  while (Date.now() < deadline) {
    try {
      const r = await fetch(`${base}/healthz`);
      if (r.ok) {
        const j = await jsonOrText(r);
        log(`健康检查通过 ${base}/healthz -> ${JSON.stringify(j)}`);
        return;
      }
    } catch (e) { lastErr = e; }
    await new Promise((r) => setTimeout(r, 300));
  }
  fail(`健康检查超时: ${lastErr?.message || 'no response'}`);
}

async function jsonOrText(r) {
  try { return await r.json(); } catch { return await r.text(); }
}

async function startOwnServer() {
  portFile = join(tmpdir(), `algal-port-${process.pid}.txt`);
  if (existsSync(portFile)) rmSync(portFile);
  const child = spawn(process.execPath, [join(ROOT, 'server.cjs')], {
    env: { ...process.env, WEB_HOST: '127.0.0.1', WEB_PORT: '0', PORT_FILE: portFile },
    stdio: ['ignore', 'pipe', 'inherit'],
  });
  ownServer = child;
  return await new Promise((res, rej) => {
    const t0 = Date.now();
    const tick = () => {
      if (existsSync(portFile)) {
        const port = Number(readFileSync(portFile, 'utf8').trim());
        if (Number.isInteger(port) && port > 0) {
          res(`http://127.0.0.1:${port}`);
          return;
        }
      }
      if (Date.now() - t0 > 10000) return rej(new Error('服务器未在 10s 内监听'));
      setTimeout(tick, 100);
    };
    tick();
  });
}
async function checkStatic(base) {
  const pages = ['/', '/index.html', '/js/lineage.js', '/js/app.js', '/css/style.css'];
  for (const p of pages) {
    const r = await fetch(`${base}${p}`);
    if (r.status !== 200) fail(`GET ${p} 状态码 ${r.status}`);
    const body = await r.text();
    if (!body.length) fail(`GET ${p} 返回空内容`);
  }
  log('静态资源全部可访问');
  const r404 = await fetch(`${base}/no-such-file`);
  if (r404.status !== 404) fail(`缺失资源应返回 404，实际 ${r404.status}`);
  log('404 行为正常');
}

// 同时含分裂与漏检的场景：
// 帧0 a → 帧1 b →（帧2 漏检）→ 帧3 c → 帧4 分裂为 e1/e2；
// 各帧还放置更亮的杂质 z*，验证不会被逐帧贪心串入。
function scenario() {
  return {
    frames: [
      [
        { id: 'a', x: 5, y: 50, b: 40 },
        { id: 'z0', x: 90, y: 90, b: 200 },
      ],
      [
        { id: 'b', x: 15, y: 50, b: 42 },
        { id: 'z1', x: 88, y: 90, b: 200 },
      ],
      [
        { id: 'z2a', x: 86, y: 90, b: 200 },
        { id: 'z2b', x: 86, y: 80, b: 190 },
      ],
      [
        { id: 'c', x: 35, y: 50, b: 44 },
        { id: 'z3', x: 84, y: 85, b: 200 },
      ],
      [
        { id: 'e1', x: 45, y: 42, b: 46 },
        { id: 'e2', x: 45, y: 58, b: 48 },
        { id: 'z4', x: 82, y: 82, b: 200 },
      ],
    ],
    startId: 'a',
    maxDist: 14,
    maxSkip: 1,
    target: 2,
  };
}

function checkScenario() {
  const input = scenario();
  const { errors, spec } = normalizeSpec(input);
  if (errors.length) fail(`场景输入校验失败: ${JSON.stringify(errors)}`);
  const raw = solveLineage(spec);
  if (!raw.feasible) fail(`含分裂与漏检的场景被误判不可行: ${JSON.stringify(raw.earliestBreak)}`);
  const sol = presentSolution(spec, raw);

  const assert = (cond, msg) => { if (!cond) fail(msg); };
  assert(sol.skips === 1, `漏检段应为 1，实际 ${sol.skips}`);
  assert(sol.divisions === 1, `分裂次数应为 1，实际 ${sol.divisions}`);
  assert(sol.survivors === 2, `终帧存活应为 2，实际 ${sol.survivors}`);
  assert(JSON.stringify(sol.used[2]) === '[]', `第 3 帧应整帧漏检，实际 ${JSON.stringify(sol.used[2])}`);
  assert(sol.used[3][0] === 'c', `第 4 帧应补获 c，实际 ${JSON.stringify(sol.used[3])}`);
  assert(JSON.stringify(sol.used[4].sort()) === JSON.stringify(['e1', 'e2']),
    `末帧应为 e1/e2，实际 ${JSON.stringify(sol.used[4])}`);
  const gap = sol.edges.find((e) => e.gap === 2);
  assert(gap && gap.fromId === 'b' && gap.toId === 'c', '漏检段应为 b→c');
  const div = sol.edges.filter((e) => e.fromFrame === 3 && e.fromId === 'c');
  assert(div.length === 2 && div.every((e) => ['e1', 'e2'].includes(e.toId)), 'c 应分裂为 e1、e2');
  assert(sol.edges.every((e) => !e.toId.startsWith('z') && !e.fromId.startsWith('z')),
    '亮杂质 z* 不得进入谱系');
  const expectedBright = 40 + 42 + 44 + 46 + 48;
  assert(sol.totalBrightness === expectedBright,
    `总亮度应为 ${expectedBright}，实际 ${sol.totalBrightness}`);
  log(`谱系烟测通过：a→b →漏检→ c →(e1,e2)，总亮度 ${sol.totalBrightness}，位移表 ${sol.edges.length} 行`);

  // 不可行场景：收紧位移使首帧间彻底断开，应报告最早断开为 帧1→帧2
  const tight = structuredClone(input);
  tight.maxDist = 2;
  const spec2 = normalizeSpec(tight).spec;
  const raw2 = solveLineage(spec2);
  assert(raw2.feasible === false, '位移收紧后应不可行');
  assert(raw2.earliestBreak.from === 0 && raw2.earliestBreak.to === 1,
    `最早断开帧间应为 1→2，实际 ${raw2.earliestBreak.from + 1}→${raw2.earliestBreak.to + 1}`);
  log('不可行报告正确：最早断开 第 1 帧 → 第 2 帧');
}

async function main() {
  let base = BASE_URL;
  if (process.env.BASE_URL) {
    log(`使用外部服务 ${base}`);
  } else {
    // Compose 的 verify 服务通过主机名 web 访问；本地直跑时自己拉起服务器
    try {
      await fetch(`${base}/healthz`, { signal: AbortSignal.timeout(800) });
    } catch {
      log(`无法连接 ${base}，改为本地临时启动服务器`);
      base = await startOwnServer();
      log(`临时服务器监听于 ${base}`);
    }
  }
  await waitHealthy(base);
  await checkStatic(base);
  checkScenario();
  log('全部烟测通过 ✔');
  if (ownServer) ownServer.kill('SIGTERM');
  if (portFile && existsSync(portFile)) rmSync(portFile);
}

main().catch((e) => {
  console.error(e.stack || e.message);
  if (ownServer) ownServer.kill('SIGTERM');
  process.exit(1);
});
