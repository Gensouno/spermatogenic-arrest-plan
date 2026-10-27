# 英文论文同步系统初始化记录

日期：2026-10-27。状态：初始化已完成；正式英文论文仍是待人审的研究草稿，不是已完成科研验证或可直接投稿的成稿。本轮未新增生物实验，也未重算此前的科学分析；使用既有冻结成果建立可追溯写作与验证机制。

## 已实施的机制

项目根规则与单一 `paper-sync` 技能共同引导后续任务：先完成项目工作及必要验证，再判断论文贡献。有实质方法、设计、结果、限制或可复现性变化时，登记证据并增量编辑相关英文章节；无贡献时保留正文。日常正文唯一来源是Markdown，Word只在要求或约定里程碑导出。

已经建立同一篇英文论文的九个章节、五组证据登记（含一组明确计划中证据）、十六项主张、五十七条结构化数值断言，以及七篇经过DOI注册元数据核对的真实论文。文献身份核对不替代对被引科学内容的阅读与审查。正文保留既有负结果、描述性结果的统计边界和来源限制，未虚构创新性、署名或完成进度。

研究状态包括进度、技术决策、术语、证据/主张、文献状态、投稿清单、语义审查及版本上下文。新会话从项目根规则读取；父工作区另有不公开的路由说明，以便在工作区根目录开始任务时找到项目规则。技能位于官方支持的 `.agents/skills/` 路径，保持隐式调用可用。配置持久化并不等于可以保证每次模型判断都正确，仍以实际检查和任务报告为准。

## 本轮创建和修改的文件

### 通用基础设施

| 文件 | 实际作用 |
|---|---|
| `AGENTS.md` | 英文写作、证据范围、任务结束检查、发布边界及环境约定 |
| `.agents/skills/paper-sync/SKILL.md` | 同一篇论文的贡献判断、增量同步与语义审查流程 |
| `envs/paper_python.sh` | 加载项目本地论文工具，复用原Python启动器及代理规则 |
| `envs/requirements-paper.txt` | 本地Word导出及检查依赖 |
| `envs/requirements-paper-check.txt` | CI/轻量检查依赖 |
| `scripts/paper_support.py` | 严格YAML、路径、校验值及结构化结果读取 |
| `scripts/check_paper.py` | 完整/基础设施/投稿阶段检查，ERROR/WARNING/INFO输出 |
| `scripts/check_sync_readiness.py` | 只读Git状态、敏感路径和令牌模式检查及提交建议 |
| `scripts/export_paper.sh`、`scripts/export_paper.py` | Pandoc章节合并、citeproc引用、模板、DOCX结构检查及导出溯源 |
| `tests/test_paper_system.py` | 合成回归案例，未嵌入未发表论文或私有实验数据 |
| `.github/workflows/paper-check.yml` | 只读权限的公开基础设施CI；不要求或上传未发表正文 |
| `docs/PAPER_WORKFLOW.md` | 后续操作、数据格式、检查边界和发布方式 |
| `docs/PAPER_SYNC_INITIALIZATION.md` | 本次初始化记录 |

修改 `.gitignore` 以排除未发表学术材料、本地工具、导出和敏感配置；在README追加工作流导航，保留此前未提交的科研改动。父工作区 `AGENTS.md` 仅追加项目规则路由，保留原代理要求，未复制到公开项目。

### 本地论文文件

`paper/metadata.yaml`、`outline.md`、`references.bib`、`figures/README.md`、`templates/README.md`及七份文献身份JSON记录；正文为：

```text
paper/chapters/00-abstract.md
paper/chapters/01-introduction.md
paper/chapters/02-related-work.md
paper/chapters/03-methods.md
paper/chapters/04-study-design.md
paper/chapters/05-results.md
paper/chapters/06-discussion.md
paper/chapters/07-conclusion.md
paper/chapters/08-declarations.md
```

本地研究登记为：

```text
research/progress.md
research/decisions.md
research/evidence.yaml
research/claims.yaml
research/terminology.yaml
research/bibliography.yaml
research/paper-todo.md
research/submission-checklist.md
research/semantic-review.yaml
research/version-context.yaml
research/automation.yaml
research/sync-records/2026-10-27-initialization.yaml
```

这些正文、研究登记及文献缓存默认被Git忽略。现有GitHub远程是公开代码仓库，没有更改可见性、创建私有远程、提交、推送或上传论文。论文首次授权Git提交或私有同步目的地仍待用户决定；不能把本地文件校验记录称为已经形成论文Git历史。

## 实测检查结果

| 检查 | 实际结果与范围 |
|---|---|
| 全套代码测试 | 101项通过，含14项新增论文系统测试 |
| 完整论文检查 | 0条ERROR、15条WARNING、2条INFO |
| 基础设施模式 | 通过；明确不代表检查过私稿 |
| 投稿严格模式 | 按预期未通过：草稿、人审、署名、目标及待办仍未完成 |
| 技能规范检查 | `quick_validate.py`通过 |
| Shell语法及差异空白检查 | 通过 |
| 既有结果完整性 | 十份快照、247份文件、226项绑定代码/规格检查通过；不是原始缓存或生物有效性认证 |
| 发布准备检查 | 未发现保护路径被追踪或匹配的令牌模式；只读，不执行提交/推送，规则扫描不能替代完整隐私审核 |
| Word导出冒烟 | Pandoc 3.9成功导出、所有XML/关系文件可解析，公式/表格/图片/引用进入DOCX；八页已渲染并完成视觉审查 |

十五条正常警告包括九章仍为草稿、四项解释性主张需要人审、二十四项未完成研究/投稿事项的汇总，以及署名/目标尚未确认的提示。机器检查不能证明真实创新、科学结论、伦理/许可、动物独立性或投稿质量。

Word测试曾发现样式前缀兼容问题，已修复并补充回归测试；又核查了短表格、图题与分页。冒烟输出和修复前中间文件仅保存在忽略目录，不作为新的日常正文来源或正式Word交付。以后更换真实学校/期刊模板仍需重新渲染与审查。

论文依赖只安装到 `.paper-tools/packages`：PyYAML 6.0.3、bibtexparser 1.4.4、pypandoc_binary 1.17（Pandoc 3.9）。没有修改系统Python、共享生信虚拟环境或R库；网络下载使用现有本地代理。Zotero本地API当时未运行，本轮未改动Zotero，而以DOI注册元数据核对文献身份。

## 每周检查与GitHub边界

已由Codex应用建立真正的线程定时检查，状态为ACTIVE，计划为北京时间每周一09:00。它检查本地工作区和私稿，仅在有新的重要变更、失败或需要用户决定时提醒；不自动提交、推送、导出、上传、投稿或外部发信。该计划受本地应用及调度可用性影响，不承诺精确唤醒时间。本轮只确认计划配置，未虚构已经发生的每周执行记录。

GitHub Actions文件已在本地配置，但没有推送到远端，因而没有声称远端CI已执行。当前公开CI只检查可公开基础设施和合成测试。以后若授权私有学术仓库，可在其发布范围内增加完整私稿检查，不把未发表材料带入公开代码提交。

## 后续使用

正常下达项目任务即可；无需重复要求写新论文。每次结束时报告项目结果、论文是否增量更新、证据与剩余不确定性、自动检查及Git准备状态。明确要求Word或到约定里程碑时才导出。

常用命令及格式说明见 [工作流指南](PAPER_WORKFLOW.md)。接下来优先补充来源核实、相关研究比较和人审，设计真正适用的统计后续，而不是为了扩写论文制造结果。作者、目标期刊、私有Git目的地和正式披露均保留为待确认事项。
