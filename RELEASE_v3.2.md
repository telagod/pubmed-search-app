# v3.2 架构优化与清理发布说明

## 亮点
- 原生 Pages 多页架构：导航更清晰、职责更单一
- 资源与数据缓存最佳实践：连接与路径绑定、数据按 mtime+令牌失效
- 配置加载更稳健：优先 Streamlit Secrets，.env 兜底
- 仓库清理：删除 archive/ 目录，保持根目录简洁

## 变更详情
- 新增 `pages/` 目录：数据管理、Dashboard、高级搜索、文献浏览、数据分析、设置、关于
- `init_db(db_path)` 以路径为键进行缓存，切换数据库自动失效
- 新增 `get_all_papers_cached` / `get_stats_cached` 与 `filter_papers_df`
- 增加缓存令牌 `db_token`：上传/清空数据库、搜索成功后自动刷新
- 移除 `archive/` 历史文件以瘦身（已在 `.gitignore` 中忽略）

## 升级指引
1. 执行 `git pull` 获取 v3.2
2. 运行 `uv run streamlit run streamlit_app.py`
3. 左侧 Pages 导航使用各功能页面
4. 首次使用先在「⚙️ 设置」配置 API，再进行「🔍 高级搜索」

## 兼容性
- 数据库结构不变；上传与导出流程保持一致
- 统计与过滤逻辑迁移至缓存化 DataFrame 层，性能更稳

---
标签：`v3.2`
标题：`v3.2 架构优化与清理`
