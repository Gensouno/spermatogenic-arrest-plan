#!/usr/bin/env python3
"""Present frozen results and documentary readiness without inventing completion scores."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from run_independent_contexts import ROOT,COHORTS
from check_all_snapshots import check


def plot(out,language):
    out.mkdir(parents=True,exist_ok=True);zh=language=='zh'
    if zh:
        matches=[p for p in font_manager.findSystemFonts() if 'Arial Unicode' in p]
        if not matches:raise RuntimeError('A CJK font is required for the Chinese presentation')
        font_manager.fontManager.addfont(matches[0]);plt.rcParams['font.family']=font_manager.FontProperties(fname=matches[0]).get_name()
    plt.rcParams.update({'axes.unicode_minus':False,'font.size':11,'savefig.facecolor':'white'})
    design=json.loads((ROOT/'replicability_snapshot/design.json').read_text())
    priorities=pd.read_csv(ROOT/'replicability_snapshot/stable_descriptive_priorities.csv')
    assert len(priorities)==36 and design['primary_BY_discovery_hypotheses']==0
    integrity=check(ROOT)
    fig,ax=plt.subplots(figsize=(13.8,7.3));ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)
    title='项目进度：计算成果已形成，生物学结论仍待验证' if zh else 'Research progress: computational results established; biological validation pending'
    ax.text(.01,.98,title,va='top',fontsize=17,weight='bold',color='#17334a')
    stats=[('0','主要检验正式发现' if zh else 'Primary adjusted discoveries'),('36','描述性下调候选' if zh else 'Descriptive downward priorities'),(str(integrity['snapshots']),'冻结分析分支' if zh else 'Frozen analysis branches')]
    for x,(value,label) in zip([.01,.35,.69],stats):
        ax.text(x,.86,value,fontsize=28,color='#245a72',va='top',weight='bold');ax.text(x,.785,label,fontsize=11,color='#45566a')
    rows=[
        ('数据与代码','已运行','主要计数分析：5个队列、40份样本、4个遗传单元','#267967'),
        ('证据与复现','计算核验通过',f"{integrity['files']}个冻结产物；{integrity['source_specification_checks']}项来源/规格校验；不等于生物验证",'#267967'),
        ('英文论文','初稿形成','单一英文稿；署名、科学人审与目标期刊仍待确认','#a16b16'),
        ('共享机制','尚未建立','候选方向稳定性 ≠ 正式复制发现 ≠ 细胞内共同机制','#a16b16'),
        ('新资源 GSE72833','暂缓定量使用','表达文件与补充资料已下载；幼年数据单位和样本对应关系待核实','#a16b16'),
        ('投稿准备','待完成','独立科学审查、伦理/许可、作者确认和投稿材料','#6b7480')
    ] if zh else [
        ('Data and code','Executed','Primary counts: 5 cohorts, 40 samples, 4 genetic units','#267967'),
        ('Evidence and reproduction','Computation checked',f"{integrity['files']} frozen artifacts; {integrity['source_specification_checks']} source/spec checks; not biological validation",'#267967'),
        ('English manuscript','Draft available','Single English draft; authorship, human review and venue pending','#a16b16'),
        ('Shared mechanism','Not established','Effect stability is not an adjusted discovery or a cell-intrinsic mechanism','#a16b16'),
        ('New GSE72833 resource','Quantitative use gated','Files acquired; juvenile units and cytology/sample mapping unresolved','#a16b16'),
        ('Submission preparation','Pending','Independent review, ethics/licenses, authorship and submission materials','#6b7480')
    ]
    for i,(label,status,detail,color) in enumerate(rows):
        y=.67-i*.102
        ax.scatter(.02,y,s=95,color=color);ax.text(.045,y,label,va='center',weight='bold',fontsize=12)
        ax.text(.31,y,status,va='center',color=color,fontsize=11)
        ax.text(.31,y-.032,detail,va='center',fontsize=10,color='#45566a')
        if i<5:ax.plot([.045,.98],[y-.058,y-.058],color='#e0e5ea',lw=.8)
    ax.text(.01,.015,'状态截至 2026-10-31；阶段状态不是完成百分比。' if zh else 'Status at 2026-10-31. Stage labels are not completion percentages.',color='#45566a',fontsize=10)
    fig.tight_layout();fig.savefig(out/'research_progress.png',dpi=170);plt.close(fig)

    retention=pd.read_csv(ROOT/'independent_contexts_snapshot/priority_summary.csv').set_index('cohort').loc[COHORTS]
    camera=pd.read_csv(ROOT/'competitive_set_snapshot/camera_results.csv')
    fig,axes=plt.subplots(1,2,figsize=(14,6.0),gridspec_kw={'width_ratios':[1,1.12]})
    labels=['Rif1 P14','Hormad SFB P19','Mael P16','Morc1 P14.5'];y=np.arange(4)
    stable=retention.all_down_1p5fold.to_numpy();measured=retention.full_measured.to_numpy()
    axes[0].barh(y,stable,color='#267967',height=.5)
    axes[0].barh(y,measured-stable,left=stable,color='#b7c9d7',height=.5)
    axes[0].barh(y,36-measured,left=measured,color='#eeeeee',edgecolor='#b8bfc5',hatch='///',height=.5,lw=.3)
    for i,(a,b) in enumerate(zip(stable,measured)):
        axes[0].text(37,i,f'{a} / {b} / 36',va='center',fontsize=10)
    axes[0].set_yticks(y,labels);axes[0].invert_yaxis();axes[0].set_xlim(0,49);axes[0].set_xticks([0,12,24,36]);axes[0].grid(axis='x',alpha=.15)
    axes[0].set_xlabel('基因数；标注：稳定 / 测得 / 固定候选' if zh else 'Genes; label = retained / measured / fixed')
    axes[0].set_title('候选保留随来源不同' if zh else 'Descriptive retention varies across sources',fontsize=13,pad=14)
    handles=[Patch(facecolor='#267967'),Patch(facecolor='#b7c9d7'),Patch(facecolor='#eeeeee',hatch='///',edgecolor='#b8bfc5')]
    legend=['全模型及6次留一均≥1.5倍下调','测得但未满足上述条件','缺失/过滤/非唯一映射'] if zh else ['Down >=1.5-fold in full + all 6 omissions','Measured; criterion not retained','Unavailable / filtered / nonunique']
    axes[0].legend(handles,legend,loc='upper left',bbox_to_anchor=(-.05,-.16),frameon=False,fontsize=9)
    for policy,offset,color,marker in [('estimated',-.09,'#245a72','o'),('fixed_0p01',.09,'#ba641b','s')]:
        d=camera[camera.policy==policy].set_index('cohort').loc[COHORTS]
        label=('估计相关性（主解释）' if zh else 'Estimated correlation (main)') if policy=='estimated' else ('固定相关性0.01（敏感性）' if zh else 'Assumed correlation 0.01 (sensitivity)')
        axes[1].scatter(-np.log10(d.BY),y+offset,c=color,marker=marker,s=55,label=label,zorder=3)
    axes[1].set_yticks(y,labels);axes[1].invert_yaxis();axes[1].set_xlim(0,16);axes[1].grid(axis='x',alpha=.15)
    axes[1].axvline(-np.log10(.05),color='#666666',lw=1,ls='--');axes[1].text(1.5,.98,'BY = 0.05',transform=axes[1].get_xaxis_transform(),va='top',fontsize=9)
    axes[1].set_xlabel('-log10（联合BY校正双侧P值）' if zh else '-log10(joint BY-adjusted two-sided P value)')
    axes[1].set_title('显著性依赖模型设置' if zh else 'Adjusted support depends on model policy',fontsize=13,pad=14)
    axes[1].legend(loc='upper left',bbox_to_anchor=(0,-.16),frameon=False,fontsize=9)
    fig.suptitle('目前得到了什么：描述性线索与明确的证据边界' if zh else 'Current findings: descriptive leads with explicit evidence limits',fontsize=16,y=.99)
    note='左：固定36个候选的描述性统计，不能作为候选集FDR发现。右：8项检验统一校正；两种设置同时改变方差膨胀和检验自由度。\n所有基因集方向均偏下；估计相关性时均未通过校正。上述结果不改变原主要检验“零个正式发现”的结论。' if zh else 'Left: descriptive follow-up of a fixed 36-gene list; no candidate-subset FDR guarantee. Right: eight jointly adjusted tests; policies change variance inflation and df.\nAll set directions are Down; none passes adjustment with estimated correlation. The original primary zero-discovery result is unchanged.'
    fig.text(.025,.025,note,fontsize=9,color='#45566a')
    fig.subplots_adjust(left=.125,right=.985,top=.80,bottom=.29,wspace=.32)
    fig.savefig(out/'research_findings.png',dpi=170);plt.close(fig)
    summary=dict(date='2026-10-31',primary_BY_discoveries=design['primary_BY_discovery_hypotheses'],fixed_descriptive_priorities=len(priorities),
        snapshot_integrity={k:integrity[k] for k in ['snapshots','files','source_specification_checks']},
        additional_contexts=retention.reset_index()[['cohort','full_measured','all_down_1p5fold']].to_dict('records'),
        estimated_policy_BY_below_0p05=int(camera.query("policy == 'estimated'").BY.lt(.05).sum()),
        assumed_policy_BY_below_0p05=int(camera.query("policy == 'fixed_0p01'").BY.lt(.05).sum()),
        scope='Presentation of frozen results and documentary stage status; no new experiment or biological validation.')
    (out/'presentation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--language',choices=['en','zh'],default='en')
    parser.add_argument('--output',type=Path,default=ROOT/'docs/figures/progress-20261031');args=parser.parse_args();plot(args.output,args.language)
