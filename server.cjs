// 零依赖静态服务器：仅提供浏览器应用；谱系计算全部在浏览器内完成。
'use strict';

const http = require('http');
const fs = require('fs');
const path = require('path');

const ROOT = path.join(__dirname, 'public');
const PORT = Number(process.env.WEB_PORT || 8080);
const HOST = process.env.WEB_HOST || '0.0.0.0';

const MIME = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.ico': 'image/x-icon',
};

function send(res, status, body, type = 'text/plain; charset=utf-8') {
  res.writeHead(status, { 'Content-Type': type, 'Cache-Control': 'no-store' });
  res.end(body);
}

const server = http.createServer((req, res) => {
  const url = new URL(req.url, `http://${req.headers.host || 'localhost'}`);
  if (url.pathname === '/healthz') {
    send(res, 200, JSON.stringify({ status: 'ok' }), 'application/json; charset=utf-8');
    return;
  }

  let rel = decodeURIComponent(url.pathname);
  if (rel === '/') rel = '/index.html';
  const file = path.normalize(path.join(ROOT, rel));
  if (file !== ROOT && !file.startsWith(ROOT + path.sep)) {
    send(res, 403, 'Forbidden');
    return;
  }

  fs.readFile(file, (err, data) => {
    if (err) {
      send(res, 404, 'Not Found');
      return;
    }
    send(res, 200, data, MIME[path.extname(file)] || 'application/octet-stream');
  });
});

server.listen(PORT, HOST, () => {
  const addr = server.address();
  const actual = typeof addr === 'object' && addr ? addr.port : PORT;
  console.log(`[web] 藻类谱系复原应用已启动: http://${HOST}:${actual}`);
  if (process.env.PORT_FILE) {
    fs.writeFileSync(process.env.PORT_FILE, String(actual));
  }
});

const shutdown = () => server.close(() => process.exit(0));
process.on('SIGTERM', shutdown);
process.on('SIGINT', shutdown);
