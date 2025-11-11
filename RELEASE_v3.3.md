# v3.3 通用化发布说明

## 亮点
- 品牌与文案全面通用化：去除 BMAL1 专属描述
- 默认查询示例改为通用主题（TP53 / cancer 等）
- About、Dashboard、标题与页脚统一为“通用 PubMed 文献检索”
- 部署文档更新为 v3.3，保持最佳实践（数据本地化、缓存）

## 代码变更
- `streamlit_app.py`
  - 标题/版本：v3.3 通用版
  - 默认查询与关键词示例调整
  - Dashboard 与 About 页面文案通用化
  - 侧栏标题与欢迎提示通用化
- `advanced_search.py`
  - 测试用默认关键词调整为 `TP53` 与 `cancer`
- `pubmed_search_v2.py`
  - 枚举策略替换为通用示例查询（肿瘤/神经/免疫/方法学）
- `STREAMLIT_DEPLOY.md`
  - 项目名与版本更新为 v3.3，检查清单与术语同步
- `README.md`
  - 改为通用项目介绍、结构、快速开始与版本记录

## 升级指引
1. `git pull` 获取最新代码
2. `uv pip install -r requirements.txt` 确保依赖一致
3. `uv run streamlit run streamlit_app.py` 启动应用
4. 使用左侧 Pages 导航进入各功能页面

## 兼容性
- 数据库结构与接口未变；BMAL1 主题仍可作为普通关键词检索

---
标签：`v3.3`
标题：`v3.3 通用化（UI/文案/示例查询）`
