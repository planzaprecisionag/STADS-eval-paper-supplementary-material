"""Reproduce S3 tables from de-identified inputs. Run: python code/reproduce_s3.py

The linkage audit documents delivered associations; it does not re-link private
raw scouting records. Historical reversal tests and supplied CIs are retained
separately; regenerated tables use the explicit seeds written in each output.
"""
from pathlib import Path
import json, math, itertools
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu
from analysis_functions import krippendorff_alpha, fleiss_kappa_from_matrix, gwet_ac1_from_matrix, pairwise_kappas
ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'data'; OUT=ROOT/'results'
SCORES={'Confirmed':1.,'Possibly Confirmed':.5,'Not Confirmed':0.}
WEIGHTS={'High':1.,'Medium':2/3,'Low':1/3}
RULES=['Strict','Moderate','Liberal','Confidence-weighted']
ratings=pd.read_csv(DATA/'blind_reviewer_ratings.csv',keep_default_na=False)
inv=pd.read_csv(DATA/'blind_pair_inventory.csv')
points=pd.read_csv(DATA/'blind_point_linkage_audit.csv')
base=ratings[ratings.presentation=='target'].merge(inv,on='pair_id',validate='many_to_one')
def mode(v):
 c=pd.Series(list(v)).value_counts()
 return '' if c.empty else ('Tie' if (c==c.iloc[0]).sum()>1 else c.index[0])
def consensus(d):
 rows=[]
 for pid,g in d.groupby('pair_id'):
  v=g[g.presence.isin(SCORES)];n=len(v);w=v.presence_confidence.map(WEIGHTS);s=v.presence.map(SCORES);ok=w.notna()
  cw=np.average(s[ok].astype(float),weights=w[ok].astype(float)) if ok.any() else np.nan
  direct=g[g.direction.isin(['Positive','Negative','Same'])].direction;dm=mode(direct)
  row=dict(pair_id=pid,field_id=g.field_id.iloc[0],crop=g.crop.iloc[0],stads_direction=g.stads_direction.iloc[0],n_presence_valid=n,n_direction_valid=len(direct),presence_evaluable=n>=2,confidence_weighted_support=cw,consensus_direction=dm,direction_evaluable=dm in ['Positive','Negative'],direction_match=dm==g.stads_direction.iloc[0],presence_mode=mode(v.presence))
  row.update({'Strict':n>=2 and (v.presence=='Confirmed').sum()>=math.ceil(.75*n),'Moderate':n>=2 and v.presence.isin(['Confirmed','Possibly Confirmed']).sum()>=math.ceil(.75*n),'Liberal':n>=2 and v.presence.isin(['Confirmed','Possibly Confirmed']).sum()>=math.ceil(.5*n),'Confidence-weighted':n>=2 and pd.notna(cw) and cw>=.5})
  row['composite']='' if n<2 else ('Full correspondence' if row['Moderate'] and row['direction_match'] else ('Presence-only correspondence' if row['Moderate'] else 'No correspondence'))
  rows.append(row)
 return pd.DataFrame(rows)
def stats(d,col,seed):
 n=len(d);k=int(d[col].sum()) if n else 0
 fs=d.groupby('field_id')[col].agg(['sum','size']);nf=len(fs)
 lo=hi=np.nan
 if nf>=2:
  rng=np.random.default_rng(seed);idx=rng.integers(0,nf,size=(10000,nf));boot=fs['sum'].to_numpy()[idx].sum(1)/fs['size'].to_numpy()[idx].sum(1);lo,hi=np.quantile(boot,[.025,.975])*100
 return dict(n=n,successes=k,percent=100*k/n if n else np.nan,n_fields=nf,ci_lower_percent=lo,ci_upper_percent=hi,bootstrap_replicates=10000 if nf>=2 else 0,seed=seed)
allc=consensus(base);allc.to_csv(OUT/'blind_pair_consensus.csv',index=False)
rows=[];seed=20260812
# Original seed progression: each rule overall, then alphabetically sorted crops.
for rule in RULES:
 e=allc[allc.presence_evaluable]
 for group,d in [('Overall',e)]+[(c,e[e.crop==c]) for c in sorted(e.crop.unique())]:
  rows.append(dict(subset='All five reviewers',endpoint=rule,group=group,**stats(d,rule,seed)));seed+=1
# Secondary analyses are newly regenerated, retaining original definitions.
for subset,d in [('No identified competing interests',base[base.reviewer_id.isin(['R1','R3','R5'])]),('All-five valid presence complete cases',base[base.pair_id.isin(allc.loc[allc.n_presence_valid==5,'pair_id'])])]:
 c=consensus(d);c.insert(0,'subset',subset)
 c.to_csv(OUT/('blind_no_competing_interests_consensus.csv' if subset.startswith('No') else 'blind_complete_case_consensus.csv'),index=False)
 e=c[c.presence_evaluable]
 for rule in RULES:
  rows.append(dict(subset=subset,endpoint=rule,group='Overall',**stats(e,rule,seed)));seed+=1
for subset,c in [('All five reviewers',allc),('No identified competing interests',consensus(base[base.reviewer_id.isin(['R1','R3','R5'])]))]:
 e=c[c.presence_evaluable]
 for category in ['Full correspondence','Presence-only correspondence','No correspondence']:
  d=e.assign(success=e.composite==category);rows.append(dict(subset=subset,endpoint=category,group='Overall',**stats(d,'success',seed)));seed+=1
 e=c[c.direction_evaluable]
 rows.append(dict(subset=subset,endpoint='Pair consensus direction agreement',group='Overall',**stats(e,'direction_match',seed)));seed+=1
# Individual binary direction ratings; same/indeterminate/excluded omitted.
d=base[base.direction.isin(['Positive','Negative'])].copy();d['success']=d.direction==d.stads_direction
for group,g in [('Overall',d)]+[(c,d[d.crop==c]) for c in sorted(d.crop.unique())]:
 rows.append(dict(subset='Reviewer-level ratings',endpoint='Direction agreement',group=group,**stats(g,'success',seed)));seed+=1
summary=pd.DataFrame(rows);summary.to_csv(OUT/'blind_endpoint_summary.csv',index=False)
# Reliability, using the supplied v5 functions.
rel=[];pairwise=[]
for endpoint,col,cats in [('Presence, all categories','presence',['Confirmed','Possibly Confirmed','Not Confirmed','Indeterminate','Excluded']),('Presence, evaluable categories','presence',list(SCORES)),('Direction, evaluable categories','direction',['Positive','Negative','Same'])]:
 d=base[base[col].isin(cats)];m=d.pivot(index='pair_id',columns='reviewer_id',values=col);counts=m.notna().sum(1);modal=int(counts.mode().iloc[0]);q=pairwise_kappas(m,cats);q.insert(0,'endpoint',endpoint);pairwise.append(q)
 rel.append(dict(endpoint=endpoint,items_any_rating=int((counts>0).sum()),items_two_or_more=int((counts>=2).sum()),fleiss_modal_rating_count=modal,fleiss_items=int((counts==modal).sum()),alpha_nominal=krippendorff_alpha(m,cats),alpha_interval=krippendorff_alpha(m,cats,'interval',SCORES) if col=='presence' and len(cats)==3 else np.nan,fleiss_kappa=fleiss_kappa_from_matrix(m,cats),gwet_ac1=gwet_ac1_from_matrix(m,cats)))
pd.DataFrame(rel).to_csv(OUT/'reviewer_reliability.csv',index=False);pd.concat(pairwise).to_csv(OUT/'reviewer_pairwise_agreement.csv',index=False)
# Confidence relationships are descriptive; majority agreement includes each reviewer.
conf=[]
for level in ['Low','Medium','High']:
 d=base[(base.presence_confidence==level)&base.presence.isin(SCORES)].copy();d['success']=d.presence.isin(['Confirmed','Possibly Confirmed']);conf.append(dict(endpoint='Presence support',confidence=level,**stats(d,'success',seed)));seed+=1
 d=d.merge(allc[['pair_id','presence_mode']],on='pair_id');d=d[~d.presence_mode.isin(['','Tie'])];d['success']=d.presence==d.presence_mode;conf.append(dict(endpoint='Agreement with panel presence mode (self-included)',confidence=level,**stats(d,'success',seed)));seed+=1
 d=base[(base.direction_confidence==level)&base.direction.isin(['Positive','Negative'])].copy();d['success']=d.direction==d.stads_direction;conf.append(dict(endpoint='Binary direction agreement with STADS',confidence=level,**stats(d,'success',seed)));seed+=1
pd.DataFrame(conf).to_csv(OUT/'reviewer_confidence.csv',index=False)
# Retention comparison: prescout attributes, never outcome-based selection.
ret=[]
for col in ['crop_source','window_id','farm_id','scene_date','stads_direction']:
 for val,g in inv.groupby(col,dropna=False):ret.append(dict(attribute=col,category=str(val),planned=len(g),retained=int(g.retained.sum()),excluded=int((~g.retained).sum()),retained_percent=100*g.retained.mean()))
pd.DataFrame(ret).to_csv(OUT/'blind_retention_comparison.csv',index=False)
sep=[]
for label,g in inv.groupby('retained'):
 x=g.planned_separation_m;sep.append(dict(retained=label,n=len(g),minimum_m=x.min(),q1_m=x.quantile(.25),median_m=x.median(),q3_m=x.quantile(.75),maximum_m=x.max()))
pd.DataFrame(sep).to_csv(OUT/'blind_retention_separation.csv',index=False)
tests=[]
for col in ['crop_source','stads_direction']:
 tab=pd.crosstab(inv[col],inv.retained);odds,p=fisher_exact(tab);tests.append(dict(attribute=col,test='Two-sided Fisher exact',statistic=odds,p_value=p))
u,p=mannwhitneyu(inv.loc[inv.retained,'planned_separation_m'],inv.loc[~inv.retained,'planned_separation_m'],alternative='two-sided');tests.append(dict(attribute='planned_separation_m',test='Two-sided Mann-Whitney U',statistic=u,p_value=p))
pd.DataFrame(tests).to_csv(OUT/'blind_retention_tests.csv',index=False)
# Regenerate reversal agreement directly from both presentations, no narrative text.
a=ratings[ratings.presentation=='target'];b=ratings[ratings.presentation=='reference'];rev=a.merge(b,on=['pair_id','reviewer_id'],suffixes=('_target','_reference'),validate='one_to_one')
flip={'Positive':'Negative','Negative':'Positive','Same':'Same','Indeterminate':'Indeterminate','Excluded':'Excluded','':''}
rev['reference_direction_reoriented']=rev.direction_reference.map(flip);rev['presence_exact_agreement']=rev.presence_target==rev.presence_reference;rev['direction_reoriented_agreement']=rev.direction_target==rev.reference_direction_reoriented
rev.to_csv(OUT/'reversal_comparisons_regenerated.csv',index=False)
rr=[]
for label,g in list(rev.groupby('reviewer_id'))+[('All reviewers',rev),('Historical masking subset R1 R2 R3 R5',rev[rev.reviewer_id!='R4']),('No competing interests R1 R3 R5',rev[rev.reviewer_id.isin(['R1','R3','R5'])])]:
 rr.append(dict(reviewer_set=label,n=len(g),presence_exact_agreement=int(g.presence_exact_agreement.sum()),presence_percent=100*g.presence_exact_agreement.mean(),direction_reoriented_agreement=int(g.direction_reoriented_agreement.sum()),direction_percent=100*g.direction_reoriented_agreement.mean()))
pd.DataFrame(rr).to_csv(OUT/'reversal_summary_regenerated.csv',index=False)
# RTK authoritative summary validation and exact bootstrap regeneration.
rtk=pd.read_csv(DATA/'rtk_analysis_records.csv');rs=pd.read_csv(OUT/'rtk_endpoint_summary.csv');rv=[]
for _,row in rs.iterrows():
 if row.analysis=='anomaly_confirmation':d=rtk[rtk.anomaly_confirmation_evaluable].copy();col='anomaly_confirmed'
 elif row.analysis=='direction_agreement':d=rtk[rtk.direction_evaluable].copy();col='direction_agreement'
 else:d=rtk[rtk.composite_evaluable].copy();d['success']=d.composite_classification==row['group'];col='success'
 if row['group'] in ['Corn silage','Soybean']:d=d[d.crop==row['group']]
 if row['group'] in ['Negative deviance','Positive deviance']:d=d[d.stads_direction==row['group'].split()[0]]
 d=d.rename(columns={'public_field_cluster_id':'field_id'});v=stats(d,col,int(row.random_seed));ok=v['n']==row.n_evaluable and v['successes']==row.n_success and np.allclose([v['ci_lower_percent'],v['ci_upper_percent']],[row.field_cluster_bootstrap_ci_lower_percent,row.field_cluster_bootstrap_ci_upper_percent])
 rv.append(dict(analysis=row.analysis,group=row['group'],verified=bool(ok),**v))
pd.DataFrame(rv).to_csv(OUT/'rtk_summary_verification.csv',index=False)
checks={'planned_blind_pairs':len(inv),'retained_blind_pairs':int(inv.retained.sum()),'linked_points':int(points.observation_id.notna().sum()),'direct_links':int((points.match_method=='direct').sum()),'inferred_links':int(points.match_method.str.startswith('inferred').sum()),'target_ratings':len(base),'reversed_pair_comparisons':len(rev),'presence_evaluable':int(allc.presence_evaluable.sum()),'moderate_confirmed':int(allc.Moderate.sum()),'composite_counts':allc.composite.value_counts().to_dict(),'rtk_summary_rows_verified':int(sum(x['verified'] for x in rv))}
assert checks['planned_blind_pairs']==99 and checks['retained_blind_pairs']==55
assert checks['direct_links']==96 and checks['inferred_links']==14
assert checks['presence_evaluable']==53 and checks['moderate_confirmed']==25
assert checks['composite_counts']=={'No correspondence':28,'Full correspondence':18,'Presence-only correspondence':7,'':2}
assert checks['rtk_summary_rows_verified']==11,rv
(OUT/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2));print(summary[(summary.group=='Overall')].to_string(index=False))
