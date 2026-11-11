import pickle, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyBboxPatch
plt.rcParams.update({'font.size':10,'font.family':'serif','mathtext.fontset':'cm'})
d=pickle.load(open('../../code/figdata.pkl','rb'))
cls=['2A','3A','4B','5A','6A','7A','8A']
cols=plt.cm.viridis(np.linspace(0.05,0.9,7))
# Fig 1: ratio sum|d_g^b|/d* vs D (3D ribbons)
fig=plt.figure(figsize=(7.2,4.6)); ax=fig.add_subplot(111,projection='3d')
for i,c in enumerate(cls):
    D=np.array([x for x,_ in d['ratio'][c]]); r=np.log10([y for _,y in d['ratio'][c]])
    ax.plot(D, np.full_like(D,i,dtype=float), r, color=cols[i], lw=1.4)
    zb=np.full_like(r,r.min()-0.5); ax.fill_between(D,np.full_like(D,i,dtype=float),r,D,np.full_like(D,i,dtype=float),zb,color=cols[i],alpha=0.25)
ax.set_yticks(range(7)); ax.set_yticklabels(cls); ax.set_xlabel('discriminant $D$'); ax.set_zlabel(r'$\log_{10}\,\sum_b|d_{g^b}|/d$')
ax.view_init(24,-58); ax.set_box_aspect((2.2,1.4,1))
plt.tight_layout(); plt.savefig('fig_ratio.png',dpi=220); plt.close()
# Fig 2: S_a/d at D=599/600 as 3D bars
fig=plt.figure(figsize=(7.2,4.4)); ax=fig.add_subplot(111,projection='3d')
for i,c in enumerate(cls):
    D,S=[x for x in d['S'][c] if x[0]==3][0]
    N=len(S)
    for a,s in enumerate(S):
        v=100*(s*N-1); ax.bar3d(a-0.35,i-0.3,min(v,0),0.7,0.6,abs(v),color=cols[i],alpha=0.85,shade=True)
ax.set_yticks(range(7)); ax.set_yticklabels(cls); ax.set_xticks(range(8)); ax.set_xlabel('eigenvalue index $a$'); ax.set_zlabel(r'$100\,(N S_a/d-1)$ at $D=3$')
ax.view_init(28,-50); ax.set_box_aspect((1.6,1.4,0.9))
plt.tight_layout(); plt.savefig('fig_sectors.png',dpi=220); plt.close()
# Fig 3: margins
fig,axs=plt.subplots(1,2,figsize=(7.4,3.3),sharey=True)
for j in (0,1):
    ax=axs[j]
    for i,c in enumerate(cls):
        if c=='2A' and j==1: continue
        m=np.array([x for x,_ in d['marg'][(c,j)]][:60]); r=np.array([y for _,y in d['marg'][(c,j)]][:60])
        ax.loglog(m,r,color=cols[i],lw=1.1,label=c)
    ax.axhline(1,color='k',ls='--',lw=0.8); ax.set_xlabel('$m$'); ax.set_title(f'$j={j}$')
axs[0].set_ylabel(r'$a_m\,\max(m,2)\,/\,(24|a_0|\sigma(m))$'); axs[0].legend(fontsize=7,ncol=2)
plt.tight_layout(); plt.savefig('fig_margins.png',dpi=220); plt.close()
# Fig 4: proof flow
fig,ax=plt.subplots(figsize=(7.2,3.2)); ax.axis('off'); ax.set_xlim(0,10); ax.set_ylim(-0.25,4.4)
def box(x,y,w,h,t,fc):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.08',fc=fc,ec='0.25',lw=0.9)); ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=8.6)
box(0.1,3.0,2.6,1.1,'fixed points of $g$\n$\\Rightarrow$ twining genus $\\phi_g$','#e8f0fa')
box(3.6,3.0,2.8,1.1,'$p^2$-coefficient $X_{g,2}$\nof $1/\\Phi_g$','#e8f0fa')
box(7.2,3.0,2.7,1.1,'polar part $d_g(2)A_{2,1}/\\eta_g$\n(only $D<0$ terms)','#fbeee6')
box(0.6,1.4,4.4,1.1,'finite part in closed form\n$\\eta_g\\psi^F_{g,1}=f_g\\phi_{-2,1}+h_g\\phi_{0,1}-2d_g(2)H$','#eaf6ea')
box(5.6,1.4,3.6,1.1,'majorant: $\\sum_b|d_{g^b}|\\leq$ coefficients\nof explicit series over $(q)_\\infty^{24}$','#eaf6ea')
box(2.2,0.0,5.6,0.95,'$d(D)-\\sum_{b\\ne 0}|d_{g^b}(D)|>0$ for every $D>0$ $\\Rightarrow$ $S_a>0$','#f6f0d8')
for (x1,y1,x2,y2) in [(2.7,3.55,3.6,3.55),(6.4,3.55,7.2,3.55),(5.0,3.0,3.4,2.5),(8.5,3.0,4.9,2.5),(5.0,1.95,5.6,1.95),(3.4,1.4,4.2,0.95),(7.4,1.4,6.0,0.95)]:
    ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='->',lw=0.9,color='0.3'))
plt.savefig('fig_flow.png',dpi=220,bbox_inches='tight'); plt.close()
# Fig 5: 3D surface of fixed-point local genus |f_alpha| coefficient growth: log10 |d_g(D)| by class
C=pickle.load(open('../../code/C_600.pkl','rb'))
fig=plt.figure(figsize=(7.2,4.4)); ax=fig.add_subplot(111,projection='3d')
allc=['1A']+cls
Ds=[D for D in range(3,600) if D%4 in (0,3)]
X,Y=np.meshgrid(Ds,range(len(allc)))
Z=np.array([[np.log10(abs(C[c][D][0]/C[c][D][1])) for D in Ds] for c in allc])
ax.plot_surface(X,Y,Z,cmap='viridis',linewidth=0,alpha=0.9,rstride=1,cstride=2)
ax.set_yticks(range(8)); ax.set_yticklabels(allc); ax.set_xlabel('discriminant $D$'); ax.set_zlabel(r'$\log_{10}|d_g(D)|$')
ax.view_init(22,-62); ax.set_box_aspect((2.2,1.4,1))
plt.tight_layout(); plt.savefig('fig_growth.png',dpi=220); plt.close()
