#!/usr/bin/env python3
"""Integrate two-condition CTAT fusion isoforms and explore isoform changes.

Input is a tab-separated table made by concatenating per-sample
``01_clustered_all.tsv`` files and adding a ``group`` column.  The program
keeps the analysis traceable to the source cluster IDs and writes core
isoforms, exploratory expression candidates, within-fusion usage changes,
dominant switches, and an aggregation-sensitivity audit.

This is an effect-size workflow.  It does not calculate p-values or FDR and
must not be presented as statistical differential expression when a condition
has no biological replicates.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

BLOCK_RE = re.compile(r"^(.*):(\d+)-(\d+):([+.-])$")
BP_RE = re.compile(r"^(.*):(\d+):([+.-])$")


def integer(value):
    try: return int(float(str(value).strip()))
    except (ValueError, TypeError): return 0


def parse_blocks(text):
    blocks=[]
    for token in filter(None,map(str.strip,str(text).split(";"))):
        m=BLOCK_RE.match(token)
        if not m: raise ValueError(f"Invalid block token: {token}")
        blocks.append((m.group(1),int(m.group(2)),int(m.group(3))))
    if not blocks: raise ValueError("Empty block chain")
    return tuple(blocks)


def junctions(blocks):
    return tuple((a[0],a[2],b[0],b[1]) for a,b in zip(blocks,blocks[1:]))


def bp(value):
    m=BP_RE.match(str(value).strip())
    return (m.group(1),int(m.group(2)),m.group(3)) if m else None


def bp_ok(a,b,tol):
    x,y=bp(a),bp(b)
    return bool(x and y and x[0]==y[0] and x[2]==y[2] and abs(x[1]-y[1])<=tol)


def j_ok(a,b,tol):
    return a[0]==b[0] and a[2]==b[2] and abs(a[1]-b[1])<=tol and abs(a[3]-b[3])<=tol


def subset(rep,cand,tol):
    if len(cand)>len(rep): return False
    return any(all(j_ok(a,b,tol) for a,b in zip(rep[o:o+len(cand)],cand))
               for o in range(len(rep)-len(cand)+1))


def chain_text(chain):
    return ";".join(f"{a}:{b}->{c}:{d}" for a,b,c,d in chain)


def write_tsv(path,rows):
    if not rows: return
    with path.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,required=True,help="Combined clustered TSV with a group column")
    p.add_argument("--group-a",required=True)
    p.add_argument("--group-b",required=True)
    p.add_argument("--mapped-reads-a",type=float,required=True)
    p.add_argument("--mapped-reads-b",type=float,required=True)
    p.add_argument("--output-prefix",type=Path,required=True)
    p.add_argument("--junction-tolerance",type=int,default=5)
    p.add_argument("--breakpoint-tolerance",type=int,default=5)
    p.add_argument("--min-expression-reads",type=int,default=5)
    p.add_argument("--min-abs-log2fc",type=float,default=1.0)
    p.add_argument("--min-fusion-reads-per-group",type=int,default=10)
    p.add_argument("--min-isoform-reads",type=int,default=3)
    p.add_argument("--min-abs-delta-usage",type=float,default=0.20)
    a=p.parse_args()
    csv.field_size_limit(sys.maxsize)
    with a.input.open(encoding="utf-8",newline="") as h: source=list(csv.DictReader(h,delimiter="\t"))
    required={"group","cluster_id","fusion_name","left_breakpoint","right_breakpoint","splice_type","representative_block_chain","total_support_reads"}
    missing=required-set(source[0] if source else [])
    if missing: p.error("Missing columns: "+", ".join(sorted(missing)))
    groups={a.group_a,a.group_b}
    valid=[]; excluded=[]
    for row in source:
        if row["group"] not in groups: continue
        blocks=parse_blocks(row["representative_block_chain"]); full=junctions(blocks)
        if any(j[1]>=j[3] for j in full):
            bad=dict(row); bad["_support"]=integer(row["total_support_reads"])
            excluded.append(bad); continue
        x=dict(row); x["_full"]=full; x["_core"]=full[1:-1] if len(full)>2 else full
        x["_length"]=sum(e-s+1 for _,s,e in blocks); x["_support"]=integer(row["total_support_reads"])
        valid.append(x)

    def compatible(seed,cand):
        return (seed["fusion_name"]==cand["fusion_name"] and seed["splice_type"]==cand["splice_type"]
                and bp_ok(seed["left_breakpoint"],cand["left_breakpoint"],a.breakpoint_tolerance)
                and bp_ok(seed["right_breakpoint"],cand["right_breakpoint"],a.breakpoint_tolerance)
                and subset(seed["_core"],cand["_core"],a.junction_tolerance))

    bins=defaultdict(list)
    for x in valid: bins[(x["fusion_name"],x["splice_type"])].append(x)
    clusters=[]
    for _,g in sorted(bins.items()):
        remaining=set(range(len(g)))
        order=sorted(remaining,key=lambda i:(-len(g[i]["_core"]),-g[i]["_length"],-g[i]["_support"],g[i]["cluster_id"]))
        while remaining:
            si=next(i for i in order if i in remaining)
            members=[i for i in order if i in remaining and compatible(g[si],g[i])]
            remaining-=set(members); clusters.append((g[si],[g[i] for i in members]))

    core=[]; member_map={}
    for seed,members in clusters:
        digest=hashlib.sha256("\t".join(sorted(x["cluster_id"] for x in members)).encode()).hexdigest()[:20]
        cid="COREISO|"+digest; member_map[cid]=members
        counts={g:sum(x["_support"] for x in members if x["group"]==g) for g in groups}
        present=[g for g in (a.group_a,a.group_b) if counts[g]]
        core.append({"core_isoform_id":cid,"presence":"SHARED" if len(present)==2 else present[0]+"_ONLY",
          "fusion_name":seed["fusion_name"],"lncRNA_involved":seed.get("lncRNA_involved",""),
          "left_breakpoint":seed["left_breakpoint"],"right_breakpoint":seed["right_breakpoint"],
          "splice_type":seed["splice_type"],"core_junction_chain":chain_text(seed["_core"]),
          "representative_full_chain":chain_text(seed["_full"]),f"{a.group_a}_reads":counts[a.group_a],
          f"{a.group_b}_reads":counts[a.group_b],"total_reads":sum(counts.values())})

    import math
    expression=[]
    for x in core:
        ca=x[f"{a.group_a}_reads"]; cb=x[f"{a.group_b}_reads"]
        log2fc=math.log2(((cb+0.5)/a.mapped_reads_b)/((ca+0.5)/a.mapped_reads_a))
        if ca+cb>=a.min_expression_reads and abs(log2fc)>=a.min_abs_log2fc:
            expression.append({**x,f"{a.group_a}_FFPM":ca/a.mapped_reads_a*1e6,f"{a.group_b}_FFPM":cb/a.mapped_reads_b*1e6,
              f"log2FC_{a.group_b}_vs_{a.group_a}":log2fc,"direction":f"{a.group_b}_UP" if log2fc>0 else f"{a.group_b}_DOWN"})

    fusion_totals=defaultdict(lambda:defaultdict(int))
    for x in core:
        for g in (a.group_a,a.group_b): fusion_totals[x["fusion_name"]][g]+=x[f"{g}_reads"]
    usage=[]
    for x in core:
        totals=fusion_totals[x["fusion_name"]]
        if min(totals[a.group_a],totals[a.group_b])<a.min_fusion_reads_per_group: continue
        ca=x[f"{a.group_a}_reads"]; cb=x[f"{a.group_b}_reads"]
        delta=cb/totals[a.group_b]-ca/totals[a.group_a]
        if ca+cb>=a.min_isoform_reads and abs(delta)>=a.min_abs_delta_usage:
            usage.append({"fusion_name":x["fusion_name"],"core_isoform_id":x["core_isoform_id"],"lncRNA_involved":x["lncRNA_involved"],
              f"{a.group_a}_isoform_reads":ca,f"{a.group_b}_isoform_reads":cb,f"{a.group_a}_fusion_reads":totals[a.group_a],
              f"{a.group_b}_fusion_reads":totals[a.group_b],f"{a.group_a}_usage":ca/totals[a.group_a],
              f"{a.group_b}_usage":cb/totals[a.group_b],"delta_usage":delta})

    switches=[]
    for fusion,totals in fusion_totals.items():
        xs=[x for x in core if x["fusion_name"]==fusion]
        if min(totals[a.group_a],totals[a.group_b])<a.min_fusion_reads_per_group: continue
        da=max(xs,key=lambda x:x[f"{a.group_a}_reads"]); db=max(xs,key=lambda x:x[f"{a.group_b}_reads"])
        if da["core_isoform_id"]==db["core_isoform_id"]: continue
        ua=da[f"{a.group_a}_reads"]/totals[a.group_a]; ub=db[f"{a.group_b}_reads"]/totals[a.group_b]
        if ua<0.30 or ub<0.30: continue
        switches.append({"fusion_name":fusion,"lncRNA_involved":da["lncRNA_involved"],f"{a.group_a}_fusion_reads":totals[a.group_a],f"{a.group_b}_fusion_reads":totals[a.group_b],
          f"{a.group_a}_dominant_isoform":da["core_isoform_id"],f"{a.group_b}_dominant_isoform":db["core_isoform_id"],
          f"{a.group_a}_dominant_usage":ua,f"{a.group_b}_dominant_usage":ub})

    prefix=a.output_prefix
    paths={"core":Path(str(prefix)+".core_isoforms.tsv"),"expression":Path(str(prefix)+".expression_candidates.tsv"),
      "usage":Path(str(prefix)+".isoform_usage_changes.tsv"),"switches":Path(str(prefix)+".dominant_switches.tsv")}
    for path in paths.values(): path.parent.mkdir(parents=True,exist_ok=True)
    write_tsv(paths["core"],core); write_tsv(paths["expression"],expression); write_tsv(paths["usage"],usage); write_tsv(paths["switches"],switches)
    summary={"source_rows":len(source),"excluded_nonmonotonic":len(excluded),"core_isoforms":len(core),
      "expression_candidates":len(expression),"usage_changes":len(usage),"dominant_switches":len(switches),
      "support_reconciles":sum(integer(x["total_support_reads"]) for x in source if x["group"] in groups)==sum(x["_support"] for x in excluded)+sum(x["total_reads"] for x in core),
      "outputs":{k:str(v) for k,v in paths.items()}}
    Path(str(prefix)+".summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__=="__main__": main()
