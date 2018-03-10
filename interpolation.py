from bspline import *
from numpy.linalg import inv

#ctr = np.array([(1,2), (2,3), (3, -3), (4,4), (5,5), (6,-5), (7,-6), (8,-5), (8, -3), (9, -4), (10 , 1)]).transpose()

ctr =np.array( [(3 , 1), (2.5, 4), (0, 1), (-2.5, 4), (-3, 0), (-2.5, -4), (0, -1), (2.5, -4), (3 , -1), (1, 1)]).transpose()


x=ctr[0,:]
y=ctr[1,:]

m = len(x) #points to interpolate
p = 3; #polygrade

n = m

s = np.linspace(0, n, n)

#knots = np.zeros(m+2*p)
knots = np.zeros(n+2*p)
knots[:(p)]=s[0]
knots[-(p):]=s[-1]
knots[(p+1):-(p)]= s[1:]

t=np.zeros((ctr.shape[0], ctr.shape[1])) #or m

t_0 = np.zeros((1,ctr.shape[0]))
t_0[:,0] = (x[1]-x[0])/(s[1]-s[0])
t_0[:,1] = (y[1]-y[0])/(s[1]-s[0])

t_n = np.zeros((1,ctr.shape[0]))
t_n[:,0] = (x[-1]-x[-2])/(s[-1]-s[-2])
t_n[:,1] = (y[-1]-y[-2])/(s[-1]-s[-2])


B_ = np.zeros((n,n+2))

for k in range(1, n):
	s_k = s[k]
	for i in range(3):
		B_[k,k+i] = BasisFunc(s_k, k+i, knots, p)  

B = B_[1:-1, 2:-2]

R = np.zeros((2, n-2))
R[:, 0] = ctr[:,1] - B_[1,1]*ctr[:,0]# BasisFunc(s[1], 1, knots, p)*(ctr[:,0])
R[:, -1] = ctr[:,-2] - B_[-2,-2]*ctr[:,-1]# BasisFunc(s[-2], n+1, knots, p)*(ctr[:,-1])
#R[:,1:-1] = ctr[:,2:-3]
R[:,1:-1] = ctr[:,2:n-2]

#P = np.zeros((m+2, 2)).transpose()
P = np.zeros((n+2, 2)).transpose()

P[:,0] = ctr[:,0]

P[:,1] = ctr[:,0]# + s[4]/3*t_0

P[:,2:-2] = np.dot(inv(B),R.transpose()).transpose()
#np.linalg.pinv(B)
P[:, -2] = ctr[:,-1]# - (1 - s[-1])/3*t_n

P[:,-1] = ctr[:,-1]

P_x = P[0, :]
P_y = P[1, :]

u3=np.linspace(0,n,(max(n*2,70)))

out = np.zeros((2, len(u3)))

for i in range(len(u3)):
	for j in range(len(P_x)):
		out[0,i] += P_x[j]*BasisFunc(u3[i], j, knots, p)
		out[1,i] += P_y[j]*BasisFunc(u3[i], j, knots, p)


plt.plot(out[0],out[1],'b',linewidth=2.0,label='B-spline curve')
plt.plot(P_x, P_y,'--',label='Control polygon',marker='x',markerfacecolor='green')
plt.plot(x,y,'k--',label='Control polygon',marker='o',markerfacecolor='red')

plt.savefig('bspline_interp.pdf')
plt.show()



