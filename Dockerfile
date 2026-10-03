FROM node:22-alpine

WORKDIR /app

# 零运行时依赖：仅复制清单与应用代码
COPY package.json ./
COPY server.cjs ./
COPY public ./public
COPY scripts ./scripts
COPY test ./test

ENV WEB_PORT=8080 \
    WEB_HOST=0.0.0.0 \
    NODE_ENV=production

EXPOSE 8080

HEALTHCHECK --interval=5s --timeout=3s --start-period=3s --retries=10 \
  CMD node -e "require('http').get('http://127.0.0.1:'+(process.env.WEB_PORT||8080)+'/healthz',r=>process.exit(r.statusCode===200?0:1)).on('error',()=>process.exit(1))"

CMD ["node", "server.cjs"]
