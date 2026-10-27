# 项目与英文论文同步维护

本项目维护一篇按章节增量更新的英文科研论文。先完成研究任务及必要验证，再判断论文贡献；常规格式化、小重构不强行增加正文。规则入口为项目根目录 [AGENTS.md](../AGENTS.md)，技能为 [paper-sync](../.agents/skills/paper-sync/SKILL.md)。项目规则和技能按官方支持的仓库位置组织，见 [OpenAI 技能文档](https://learn.chatgpt.com/docs/build-skills)与 [AGENTS 文档](https://learn.chatgpt.com/docs/agent-configuration/agents-md)。

## 文件与隐私

- `paper/metadata.yaml`：唯一论文身份、章节顺序、英文语言、未确认作者及投稿目标。
- `paper/chapters/*.md`：日常正文唯一来源，包含章节目的、状态和待补证据。
- `paper/references.bib`：真实文献；`research/bibliography.yaml`记录核实状态、发表状态和元数据校验值。
- `research/evidence.yaml`：证据、校验值、配置、版本、复现命令、验证范围和限制。
- `research/claims.yaml`：主张、证据关联、章节位置和关键数值断言。
- `research/progress.md`、`decisions.md`、`terminology.yaml`、`paper-todo.md`：持续研究状态。
- `research/semantic-review.yaml`：Codex 语义审查的范围、发现及被审查文件校验值。
- `research/submission-checklist.md`：署名、伦理/许可、利益冲突、语言、匿名、AI披露和投稿材料待办。

现有代码远程仓库是公开仓库。未发表正文、研究登记、文献缓存、导出及本地工具默认被 Git 忽略，仍由本地论文检查读取；不自动上传、投稿或改变仓库设置。通用规则和检查基础设施可以单独准备代码提交。首次论文 Git 提交或私有同步目的地尚需用户授权；不能将本地文件校验记录说成已经形成论文 Git 历史。

## 本地依赖与检查

论文工具安装在忽略的 `.paper-tools/packages`，不修改共享生信环境：

```bash
# 使用既有代理配置；只安装到项目目录
bash envs/python.sh -m pip install --target .paper-tools/packages -r envs/requirements-paper.txt

# 完整检查本地未发表论文及证据
bash envs/paper_python.sh scripts/check_paper.py

# 仅检查可公开的基础设施；不代表检查过论文
bash envs/paper_python.sh scripts/check_paper.py --infrastructure-only

# 投稿阶段严格模式；草稿待办或缺少人审会阻止通过
bash envs/paper_python.sh scripts/check_paper.py --submission-ready

# 只读准备提交建议，不执行提交或推送
bash envs/paper_python.sh scripts/check_sync_readiness.py --report logs/paper-sync/git-readiness.json
```

公开克隆默认不含论文正文和研究登记，完整检查会明确报缺失，不能自动生成替代稿；获得授权的本地/私有副本后才运行完整检查。CI使用基础设施模式和独立的合成测试，不上传私稿，不创建PR、提交或其他外部消息。工作流须经授权推送后才会在GitHub执行，不能把本地配置说成远端CI已运行。

## 证据与主张格式

证据状态为 `planned/draft/verified/superseded`；`verified`必须带有具体范围，如原始文件身份、数值复现或冻结结果完整性，而不是默认生物验证。记录来源路径/SHA-256、基准Git提交与未提交状态、配置、依赖、适用随机种子/不适用原因、日期或不可变结果版本、命令和限制。

主张通过 `evidence_ids`关联证据，`chapters`记录正文位置。正文段落前用 `<!-- claim: C-ID -->` 标记；重要结果的 `assertions` 选择冻结CSV单行/计数/求和或JSON字段，记录预期值及该主张段落中的实际文本。检查器不能穷尽自然语言里所有未登记数字或保证段落语义，所以仍需审查完整正文。

来源代码、规格或结果校验值变化会使相关证据失效。应复跑或降级主张；不得单纯更新哈希来接受变化。无法确认的历史运行日期保留为空，以实际不可变版本标识，不编造日期或声称未提交代码已经存在于基准commit。

## 增量同步与语义审查

完成主任务后，判断它是否改变方法、实验设计、结果、解释、重要限制或可复现性。有贡献才登记证据并编辑相关英文章节；无贡献则保留正文，只维护必要状态。重要更新后，逐项比较方法与实现、配置与实际运行、结果与表格、摘要与结论、术语和否定结论；在语义审查记录中绑定被审查文件哈希。新的解释仍需人类科研审核。

检查输出分为ERROR、WARNING、INFO。错误必须修复；草稿、未完成投稿项目及人工审查提示是正常的可见警告。结构通过不等于内容真实、创新成立或可以投稿。

## Word 导出

```bash
# 默认创建新的本地导出；不覆盖已有版本
bash scripts/export_paper.sh --output paper/build/manuscript-milestone-01.docx

# 指定真实学校/期刊样式模板
bash scripts/export_paper.sh --output paper/build/manuscript-milestone-02.docx --reference-doc paper/templates/reference.docx
```

Pandoc按元数据章节顺序合并Markdown，处理表格、数学公式、图片、中文字符及citeproc参考文献；本项目正式论文检查另要求英文内容。没有模板时使用通用A4黑色标题样式，不编造学校或期刊规范。导出会核对DOCX包、图片和未解析引用，并生成来源校验记录；交付前还必须渲染并逐页审查。初始化可在忽略目录做一次冒烟测试，正常任务不重复导出Word。

## GitHub 与定时检查

每个重要里程碑先准备有意义的英文 Conventional Commit 摘要，运行适当代码测试、论文检查和敏感路径扫描，再核对用户授权与实际发布范围。提交和推送是两个动作；不强推、不覆盖他人变更、不把旧结果写成新验证、不默认公开论文。

每周本地检查由真正的调度工具承担，检查本地Git状态和未发表论文，仅在有值得处理的变化或失败时提醒。调度是否已创建及其ID保存在本地研究状态；AGENTS规则和远端CI不等于本地定时任务。离线或应用未运行会影响调度执行，不能保证精确到时。自动任务不执行提交、推送、投稿、外部发信或材料上传。
