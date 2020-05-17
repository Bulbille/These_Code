#!/usr/bin/python
# -*- coding: utf8
import numpy as np
from PIL import Image, ImageDraw,ImageFont

im = Image.new('RGB',(1000,1000),(255,255,255))
d0 = ImageDraw.Draw(im)
blue = (0,0,255)

L = 8
dl = 15

def s(a) :
    return (a+1)*100

## Hamiltonien normal
for i in np.arange(L):
    for j in np.arange(L) :
        i,j
        d0.ellipse((s(i)-dl,s(j)-dl,s(i)+dl,s(j)+dl), fill=blue, outline=(0, 0, 0))
        if i < L-1 :
            d0.line((s(i),s(j),s(i+1),s(j)),fill=blue,width=int(dl/3))
        if j < L-1 :
            d0.line((s(i),s(j),s(i),s(j+1)),fill=blue,width=int(dl/3))

d0.text((500,500), 'Test',textsize=140, fill=(255,0,0) )
print(d0.textsize("Test"))


im.show()
im.save('test.pdf')
exit()

dH0.arrow(-0.5,-0.5,L,0,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH0.arrow(L,-0.5,-L,0,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH0.text(L/2-0.5,-1.5,'$L\'$',color="blue",fontsize=20)

dH0.arrow(L,-0.5,0,L-0.5,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH0.arrow(L,L,0,-L,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH0.text(L+0.5,L/2-0.5,'$L$',color="blue",fontsize=20)

dH0.axis('off')
#dH0.axis('square')
#dH0.set_aspect('equal')
plt.tight_layout(pad=-1.0) 
dH0.set_xlim([-1,L+1])
dH0.set_ylim([-1,L+1])
plt.savefig('cross-h0.pdf',bbox_inches='tight')
plt.show()
