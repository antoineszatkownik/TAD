import matplotlib as mtl
from matplotlib.colors import LogNorm

"""
score contient les score des positions 16 à n
label contient les labels des scores 16 à n
"""

#x_norm = (x - np.min(x)) / (np.max(x) - np.min(x))
#x_rot = ndimage.rotate(x_norm,45)[0:(len(contact_map)//2)+1]
#x_rot = rotate45(np.triu(x[3500:4001,3500:4001]))
s1,s2 = 500,700
ax = sns.heatmap(np.tril(x[s1:s2,s1:s2]),norm=LogNorm())#cmap = sns.color_palette("flare", as_cmap=True))
#ax.invert_yaxis()
for i in range(len(arrec)):
    if ml[i]>=0.:
        rect = mtl.patches.Rectangle((arrec[i,0],arrec[i,0]),
                                     arrec[i,1]-arrec[i,0], arrec[i,1]-arrec[i,0],
                                     linewidth=1, edgecolor='b', facecolor='none')
        #ax.add_patch(rect)
rec = mtl.patches.Rectangle((20,20),33,33, color = 'gray' )
cent = mtl.patches.Rectangle((36,36),1,1, color = 'black' )
#ax.add_patch(rec)
#ax.add_patch(cent)

for i,s in enumerate(score[s1:s2]):
    if label[i]==1:
        ax.plot((i+16,
                 i+16+np.cos(-np.pi/4)*(50*score[i])),
                (i+16,
                 i+16+np.sin(-np.pi/4)*(50*score[i]) +1)  , c="g")
        
        #rec = mtl.patches.Rectangle((i+16,i+16),1,1, color = "g" )
        #ax.add_patch(rec)
        
    elif label[i]==0  :
        ax.plot((i+16,
         i+16+np.cos(-np.pi/4)*(50*score[i])),
        (i+16,
         i+16+np.sin(-np.pi/4)*(50*score[i]) +1)  , c="r")

        #rec = mtl.patches.Rectangle((i+16,i+16),1,1, color = "r" )
        #ax.add_patch(rec)

ax.plot((0,s2-s1),(np.sin(-np.pi/4)*(80)+1,(s2-s1)+np.sin(-np.pi/4)*(80) +1), c="black", linestyle = '--')
ax.annotate("Score threshold: 0.8", (140,120), rotation=-45)

ax.set_xticklabels([x*25/1000 for x in range(s1,s2,8)])
ax.set_yticklabels([x*25/1000 for x in range(s1,s2,8)])
ax.set_title("Position scoring on chromosome X 12.5 - 17.5 Mb")
ax.set_xlabel("Position [Mb]")
ax.set_ylabel("Position [Mb]")

mtl.pyplot.show()