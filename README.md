# 2026peroxisome

IJMS 投稿稿件配套的分析代码与数据（脓毒症 / 过氧化物酶体相关基因）。

私有仓库，用于稿件修订期间的版本管理与审稿人材料留存。

## 目录结构

```
.
├── code.R                  # 全部分析主脚本
├── README.md
├── .gitignore
└── data/                   # 输入数据与中间结果
```

## 已纳入版本管理的文件

| 文件 | 大小 | 说明 |
|---|---|---|
| `code.R` | 192 KB | 分析主脚本 |
| `data/interactions.tsv` | 11.6 MB | 互作关系表 |
| `data/expr_all_main.txt` | 3.5 MB | 全基因表达矩阵 |
| `data/patient_pseudobulk_pDC.tsv` | 1.7 MB | 病人层面 pseudobulk 表达 |
| `data/GSE290679-test-expr.txt` | 0.9 MB | 验证集表达矩阵 |
| `data/GSE290679-test-metadata.txt` | 1 KB | 验证集样本信息 |
| `data/cibersort.ref.txt` | 147 KB | CIBERSORT 参考矩阵 |
| `data/patient_metadata.tsv` | 5 KB | 病人临床信息 |
| `data/GenesetInfo.txt` | 2 KB | 基因集注释 |
| `data/113_ML_methods.txt` | 2 KB | 机器学习方法清单 |
| `data/Immunomodulator_and_chemokines.txt` | 4 KB | 免疫调节因子清单 |
| `data/hallmark.gs.rdata` | 12 KB | Hallmark 基因集 |
| `data/signature.txt` | 459 B | 过氧化物酶体相关基因签名 |
| `data/metadata_all_main.txt` | 4 KB | 主队列样本信息 |
| `data/final_pdf_manifest.txt` | 2 KB | 图表产出清单 |
| `data/input_manifest.tsv` | 2 KB | 输入文件清单 |

## 大数据文件（未纳入 git）

以下文件超过 GitHub 单文件 100 MB 硬性上限，或体积过大不宜入库，请从 **Releases** 页面（tag `v1.0`）下载获取完整数据：

| 文件 | 大小 | Release 附件名 |
|---|---|---|
| 全部文件打包 | 1281 MB | `data_code.zip` |
| `data/seurat.annotated.rds` | 3006 MB | `seurat.annotated.rds.part-00` + `part-01` |
| `data/mircode_highconsfamilies.txt` | 121 MB | `mircode_highconsfamilies.txt` |
| `data/keygene_fullcell_lognormalized.rds` | 105 MB | `keygene_fullcell_lognormalized.rds` |
| `data/pDC_fullgene_counts.rds` | 97.5 MB | `pDC_fullgene_counts.rds` |
| `data/fullgene_patient_log2cpm.rds` | 47.9 MB | `fullgene_patient_log2cpm.rds` |
| `data/pDC_fullcell_no_sampling.rds` | 19.0 MB | `pDC_fullcell_no_sampling.rds` |

`seurat.annotated.rds` 约 3 GB，超过 Release 单文件 2 GB 上限，已拆分为两段。下载后在同目录执行：

```bash
cat seurat.annotated.rds.part-00 seurat.annotated.rds.part-01 > seurat.annotated.rds
```

校验（拼接后 sha256 应为 `9c0497ae825ba428c571356709fc6018ff8adb61c1355874dbc2000a4d1cba9a`）：

```bash
sha256sum seurat.annotated.rds
```

若只需完整数据、不关心分文件，直接下载 `data_code.zip` 即可，其中已包含全部文件。

## 运行环境

- R（版本见 `code.R` 头部注释）
- 主要依赖：Seurat、CIBERSORT、clusterProfiler 等
- 工作目录需设为仓库根目录，`data/` 下文件按上表就位后可直接运行 `code.R`

## 数据可获得性

原始单细胞数据见 GEO 编号 GSE290679。
