"""Regenerate Figure 4 (self-efficacy trajectories + HLM fits) and Figure 5 (MI indicators) from Psyclient_data.sav.
Run from the Psyclient培训 folder:  python figures_v30/make_figures.py
"""
import os,sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from style import *
import pyreadstat, numpy as np, pandas as pd, statsmodels.formula.api as smf
setup()
df,_=pyreadstat.read_sav(os.path.join(HERE,'..','Psyclient_data.sav'))
df['grp']=df['Group'].map({1:'Feedback',2:'No-feedback',3:'Lecture-only'})
t=np.arange(1,7); tt=np.linspace(1,6,60)
# ---------- Figure 4 ----------
fig,axes=plt.subplots(1,3,figsize=(180/25.4,58/25.4),sharey=True)
for ax,(code,title),letter in zip(axes,[('HS','Helping Skills'),('SM','Session Management'),('CC','Counseling Challenges')],'abc'):
    long=df.melt(id_vars=['ID','grp'],value_vars=[f'{code}T{i}' for i in t],var_name='wave',value_name='y').dropna()
    long['time']=long['wave'].str.extract(r'T(\d)').astype(float)-1
    long['grp']=pd.Categorical(long['grp'],categories=ORDER)
    m=smf.mixedlm('y ~ (time + I(time**2)) * C(grp)',long,groups=long['ID']).fit(reml=True)
    for k,g in enumerate(ORDER):
        sub=df[df.grp==g]; cols=[f'{code}T{i}' for i in t]
        mean=sub[cols].mean().values; se=sub[cols].std().values/np.sqrt(sub[cols].count().values)
        ax.errorbar(t+(k-1)*0.07,mean,yerr=se,color=COL[g],lw=0,marker='o',ms=3.2,mfc='white',mew=1.1,capsize=1.6,elinewidth=0.7,label=g,zorder=3)
        X=pd.DataFrame({'time':tt-1,'grp':pd.Categorical([g]*len(tt),categories=ORDER)})
        ax.plot(tt,m.predict(X),color=COL[g],lw=1.0,ls=(0,(4,2)),zorder=2)
    ax.set_title(title,pad=4); ax.set_xticks(t); ax.set_xticklabels([f'T{i}' for i in t]); ax.set_xlim(0.6,6.4)
    ax.set_ylim(4,8); ax.set_yticks([4,5,6,7,8]); ax.set_xlabel('Assessment week'); panel_label(ax,letter,x=-0.16 if letter=='a' else -0.08)
axes[0].set_ylabel('Self-efficacy (item mean, 0–9)'); axes[0].legend(loc='lower left',handlelength=1.6,borderaxespad=0.2,labelspacing=0.3)
fig.subplots_adjust(left=0.075,right=0.99,top=0.88,bottom=0.2,wspace=0.12)
fig.savefig(os.path.join(HERE,'figure4_fit.png'),dpi=400); fig.savefig(os.path.join(HERE,'figure4_fit.pdf'))
# ---------- Figure 5 ----------
rng=np.random.default_rng(7)
fig,axes=plt.subplots(1,2,figsize=(130/25.4,62/25.4))
for ax,(col,ylab,refs,ylim),letter in zip(axes,[('RQR_new_post','Reflection-to-question ratio',[(1.0,'fair'),(2.0,'good')],(0,8.6)),('questioning_post','Open questions (%)',[(50,'fair'),(70,'good')],(0,103))],'ab'):
    for k,g in enumerate(ORDER):
        v=df.loc[df.grp==g,col].dropna().values; x=k+rng.uniform(-0.17,0.17,len(v))
        ax.scatter(x,v,s=9,color=COL[g],alpha=0.55,linewidths=0,zorder=2,clip_on=False)
        mu=v.mean(); se=v.std(ddof=1)/np.sqrt(len(v))
        ax.errorbar(k,mu,yerr=se,fmt='D',ms=3.6,color='black',capsize=2.2,elinewidth=0.8,lw=0.8,zorder=4); ax.hlines(mu,k-0.27,k+0.27,color='black',lw=0.9,zorder=4)
    for yv,lab in refs: ax.axhline(yv,ls=(0,(3,2)),lw=0.6,color='#999999',zorder=1); ax.text(2.55,yv,lab,fontsize=6,color='#777777',va='bottom',ha='right')
    ax.set_xticks(range(3)); ax.set_xticklabels(ORDER,fontsize=6.8); ax.set_xlim(-0.55,2.55); ax.set_ylabel(ylab); ax.set_ylim(*ylim); panel_label(ax,letter,x=-0.2 if letter=='a' else -0.24)
fig.subplots_adjust(left=0.1,right=0.99,top=0.9,bottom=0.12,wspace=0.4)
fig.savefig(os.path.join(HERE,'figure5_new.png'),dpi=400); fig.savefig(os.path.join(HERE,'figure5_new.pdf'))
# ---------- Figure 6 ----------
# needs output/analysis/rct_scale_totals_long.csv and output/r_analysis/feedback_group_prediction_dataset.csv
rct=pd.read_csv(os.path.join(HERE,'..','output','analysis','rct_scale_totals_long.csv'),dtype={'full_id':str})
fb=pd.read_csv(os.path.join(HERE,'..','output','r_analysis','feedback_group_prediction_dataset.csv'),dtype={'full_id':str})
t1=rct[rct.timepoint=='T1'][['full_id','self_efficacy_total']].rename(columns={'self_efficacy_total':'se_t1'})
t6=rct[rct.timepoint=='T6'][['full_id','self_efficacy_total']].rename(columns={'self_efficacy_total':'se_t6'})
m6=fb[['full_id','overall_mean_score_weighted']].merge(t1,on='full_id').merge(t6,on='full_id').dropna()
fit=smf.ols('se_t6 ~ se_t1 + overall_mean_score_weighted',m6).fit()
xs=np.linspace(m6.overall_mean_score_weighted.min(),m6.overall_mean_score_weighted.max(),50)
pred=fit.get_prediction(pd.DataFrame({'se_t1':m6.se_t1.mean(),'overall_mean_score_weighted':xs})).summary_frame(alpha=0.05)
fig,ax=plt.subplots(figsize=(80/25.4,62/25.4))
ax.fill_between(xs,pred['mean_ci_lower'],pred['mean_ci_upper'],color=COL['Feedback'],alpha=0.12,lw=0); ax.plot(xs,pred['mean'],color=COL['Feedback'],lw=1.2)
ax.scatter(m6.overall_mean_score_weighted,m6.se_t6,s=16,color=COL['Feedback'],alpha=0.85,linewidths=0,zorder=3)
ax.set_xlabel('Weighted copying score'); ax.set_ylabel('T6 total self-efficacy'); ax.set_xlim(0.1,1.0); ax.set_xticks([0.2,0.4,0.6,0.8,1.0]); ax.set_ylim(150,400)
ax.text(0.98,0.97,'b = %.2f, p = %s\nadj. R² = %.3f'%(fit.params['overall_mean_score_weighted'],('%.3f'%fit.pvalues['overall_mean_score_weighted']).lstrip('0'),fit.rsquared_adj),transform=ax.transAxes,ha='right',va='top',fontsize=6.5)
fig.subplots_adjust(left=0.17,right=0.95,top=0.97,bottom=0.17)
fig.savefig(os.path.join(HERE,'figure6_new.png'),dpi=400); fig.savefig(os.path.join(HERE,'figure6_new.pdf'))
print('done')
