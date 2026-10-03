// 零依赖构建：语法检查全部脚本、校验页面引用的静态资源、产出构建清单。
'use strict';

import { spawnSync } from 'node:child_process';
import { readFileSync, writeFileSync, mkdirSync, existsSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
let failed = 0;

function check(file) {
  const r = spawnSync(process.execPath, ['--check', file], { encoding: 'utf8' });
  if (r.status !== 0) {
    failed++;
    console.error(`✗ 语法检查失败: ${file}\n${r.stderr}`);
  } else {
    console.log(`✓ 语法检查 ${file.replace(ROOT + '/', '')}`);
  }
}

[
  'server.cjs',
  'public/js/lineage.js',
  'public/js/app.js',
  'scripts/build.mjs',
  'scripts/smoke.mjs',
].forEach((f) => check(join(ROOT, f)));

// 校验 index.html 引用资源存在
const html = readFileSync(join(ROOT, 'public/index.html'), 'utf8');
const refs = [...html.matchAll(/(?:src|href)="\.?\/?([^"]+)"/g)].map((m) => m[1]);
for (const ref of refs) {
  const p = join(ROOT, 'public', ref);
  if (!existsSync(p)) {
    failed++;
    console.error(`✗ index.html 引用缺失: ${ref}`);
  } else {
    console.log(`✓ 资源存在 ${ref}`);
  }
}

const distDir = join(ROOT, 'dist');
mkdirSync(distDir, { recursive: true });
writeFileSync(
  join(distDir, 'build-info.json'),
  JSON.stringify({ builtAt: new Date().toISOString(), node: process.version, files: refs }, null, 2) + '\n'
);
console.log('✓ 构建清单 dist/build-info.json');

if (failed) {
  console.error(`构建失败：${failed} 项错误`);
  process.exit(1);
}
console.log('构建通过。');
