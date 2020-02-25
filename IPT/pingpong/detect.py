import matplotlib.pyplot as plt
import numpy as np
import cv2

vidcap = cv2.VideoCapture('bon.mp4')
maxt = int(vidcap.get(cv2.CAP_PROP_FRAME_COUNT))

position = np.empty((maxt,3))

lower_orange = np.array([200,180,0])
upper_orange = np.array([255,255,80])

for count in np.arange(maxt):
    success,image = vidcap.read()
    image = cv2.cvtColor(image,cv2.COLOR_BGR2RGB)
    if success != True :
        position[count] = [count,-1,-1]
        maxt = count
        break
    mask = cv2.inRange(image, lower_orange, upper_orange)
    truemask = np.where(mask!=0)
    xbar = 0 
    ybar = 0 
    nbar = len(truemask[0])
    if nbar > 0 :
        ybar = np.sum(truemask[0])/nbar
        xbar = np.sum(truemask[1])/nbar
        position[count] = [count,xbar,ybar]
    else :
        position[count] = [count,-1,-1]

position = position[np.where(position[:,1] > 0)]

#plt.plot(position[:,0],position[:,1],label="Position x")
#plt.plot(position[:,0],position[:,2],label="Position y")
#plt.savefig('graphe.pdf')
#plt.show()
#plt.close()
