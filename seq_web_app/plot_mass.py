import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
import glob
import errno

#path = '../Result/*.txt'

x = []
y = []
n = []
x1 = []
y1 = []
n1 = []
x2 = []
y2 = []
n2 = []
allx=[]
ally=[]
color=[]
color1=[]
color2=[]
allcolor=[]
count = 1

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


fig, ax = plt.subplots(figsize=(16,8))




f2 = open("uploadedFiles/LCMS_data.txt")
for line2 in f2:
    splitArray2 = line2.strip("\n").split("\t")
    #print splitArray2
    x2.append(float(splitArray2[0]))
    #allx.append(float(splitArray2[0]))
    y2.append(float(splitArray2[1]))
    #ally.append(float(splitArray2[1]))
    color2.append(float(splitArray2[2]))
    #allcolor.append(float(splitArray2[2]))
f2.close()

#s2 = ax.scatter(x2, y2, c=color2, cmap='cool')
s2 = ax.scatter(x2, y2, facecolors='b', edgecolors='b')

f1 = open("Result/finalRead.txt")
for line1 in f1:
	splitArray = line1.strip("\n").split("\t")
	x.append(float(splitArray[0]))
	allx.append(float(splitArray[0]))
	y.append(float(splitArray[1]))
	ally.append(float(splitArray[1]))
	n.append(splitArray[2])
	color.append(float(splitArray[3]))
	allcolor.append(float(splitArray[3]))
f1.close()

s = ax.scatter(x, y, marker="s", facecolors='r', edgecolor = 'k')
ax.plot(x, y, label="Final Draft Sequence ")

for i, txt in enumerate(n):
    if (i % 2 == 0):
        ax.annotate(txt, (x[i],y[i]), (x[i]-50,y[i]-1), fontsize=12, weight="bold")
    else:
        ax.annotate(txt, (x[i],y[i]), (x[i]-150,y[i]+0.5), fontsize=12, weight="bold")

#files = glob.glob(path)

#for name in files: # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
#	try:
#		with open(name) as f: # No need to specify 'r': this is the default.
#			#sys.stdout.write(f.read())
#			for line in f:
#				splitArray = line.strip("\n").split("\t")
#				#print splitArray
#				x.append(float(splitArray[0]))
#				allx.append(float(splitArray[0]))
#				y.append(float(splitArray[1]))
#				ally.append(float(splitArray[1]))
#				n.append(splitArray[2])
#				color.append(float(splitArray[3]))
#				allcolor.append(float(splitArray[3]))
#			#s = ax.scatter(x, y, marker="s", edgecolor = 'k', c=color, cmap='cool')
#			s = ax.scatter(x, y, marker="s", edgecolor = 'k')
#			ax.plot(x, y, label="Draft Sequence " + str(count))
#			
#			for i, txt in enumerate(n):
#				ax.annotate(txt, (x[i],y[i]), (x[i]-150,y[i]+0.4), fontsize=12, weight="bold")
#			count+=1
#			x=[]
#			y=[]
#			n=[]
#			color=[]
#	except IOError as exc:
#		if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
#			raise # Propagate other kinds of IOError.



#print len(x)
#print len(y)
#print len(n)



#ax.plot(x, y, label="3' Biotin labeled 19nt RNA Sequence")
#ax.plot(x1, y1, label = 'Unlabeled 19nt RNA Sequence')


#plt.title("b) 3' biotinylated ladder; 5' unlabeled ladder", loc='left', fontsize=14)

#alls = ax.scatter(allx, ally, c=allcolor, cmap='cool')

#s1 = ax.scatter(x, y, marker="o", edgecolor = 'k', c=color, cmap='cool')

#s2 = ax.scatter(x1, y1, marker="s", edgecolor = 'k', c=color1, cmap='cool')


#colorbar = plt.colorbar(s2, format=ticker.FuncFormatter(fmt))
#colorbar.set_label('Intensity (vol)', fontsize=16)

#colorbar = plt.colorbar(s2, format=ticker.FuncFormatter(fmt))
#colorbar.set_label('Intensity (vol)', fontsize=16)

#for font_objects in colorbar.ax.yaxis.get_ticklabels():
#    font_objects.set_size(12)


ax.legend(loc='upper left')

plt.xlabel("Mass (Da)", fontsize=16)
plt.ylabel("Retention Time (Min)", fontsize=16)



#ax.set_ylim(ymin=0)
#ax.set_ylim(ymax=40)

#for i, txt in enumerate(n):
#    ax.annotate(txt, (x[i],y[i]), (x[i]-150,y[i]+0.4), fontsize=20, weight="bold")

#for i, txt in enumerate(n1):
#    ax.annotate(txt, (x1[i],y1[i]), (x1[i]-150,y1[i]+0.4), fontsize=20, weight="bold")


#ax1.set_ylim(ymin=0)
#ax1.set_ylim(ymax=11)
#ax2.set_ylim(ymin=0)
#ax2.set_ylim(ymax=11)

#plt.show()
plt.savefig('Result/finalRead.png', dpi = 100)




